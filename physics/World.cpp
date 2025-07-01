#include "World.h"

#include <algorithm>

static _forceinline float Clamp(float x, float lo, float hi)
{
    return std::max(lo, std::min(x, hi));
}

static void ComputeConstraint(World &World, ConstraintIndex Ci, BodyIndex ForBi, float &C, Vector2 &Jacobian, Vector2 &Gd, float &NewLambda)
{
    const Constraint &con = *World.GetConstraint(Ci);
    switch (con.Type)
    {
        case Normal:
        {
            const Body &a = *World.GetBody(con.A);

            if (con.B == INVALID)
            {
                // negative = penetrating
                C = std::min(a.Pos.Dot(con.Normal) - (con.Rest + a.Radius), 0.f);
                Jacobian = C < 0 ? con.Normal : Vector2{0, 0};
            }
            else
            {
                const Body &b = *World.GetBody(con.B);

                const float multiplier = ForBi == con.A ? 1 : -1;

                const Vector2 Vn = (b.Pos - a.Pos).GetNormal();

                C = multiplier * std::min((b.Pos - a.Pos).Dot(Vn) - a.Radius - b.Radius, 0.f);
                Jacobian = C < 0 ? Vn * multiplier : Vector2{0, 0};
            }

            NewLambda = con.Lambda + con.K * C;
            Gd = {0.f, 0.f};
        }
        break;

        default:
            throw std::runtime_error("Unknown Constraint Type");
    }
}

// Solve A x = b for a symmetric 2×2 matrix
static _forceinline Vector2 SolveSymmetric2x2(const Matrix2x2& A, const Vector2& b)
{
    const float det = A.Values[0] * A.Values[3] - A.Values[1] * A.Values[1];      // A is [m00 m01; m01 m11]
    if (fabsf(det) < 1e-20f)
    {
        throw std::runtime_error("Bad determinant");
    }

    const float invDet = 1.0f / det;
    return {
        invDet * (  A.Values[3] * b.X - A.Values[1] * b.Y ),
        invDet * ( -A.Values[1] * b.X + A.Values[0] * b.Y )
    };
}


void UpdateBody(World *World, WorldSettings &Settings, BodyIndex Bi, float InvDt2)
{
    Body &B = *World->GetBody(Bi);
    Vector2 F = (B.Pos - B.Scratch.Y) * -1 * B.Mass * InvDt2;
    Matrix2x2 H(B.Mass * InvDt2);

    ConstraintIndex Ci = World->GetStartConstraint(Bi);
    while (Ci != INVALID)
    {
        const Constraint &Con = *World->GetConstraint(Ci);
        const bool IsA = (Con.A == Bi);

        float C;
        Vector2 CJacobian;
        Vector2 Gd;
        float UnclampedNewLambda;
        ComputeConstraint(*World, Ci, Bi, C, CJacobian, Gd, UnclampedNewLambda);

        F -= CJacobian * Clamp(UnclampedNewLambda, Con.LambdaMin, Con.LambdaMax);
        H += CJacobian.OuterProduct(CJacobian) * Con.K;
        H.Values[0] += Gd.X * UnclampedNewLambda;
        H.Values[3] += Gd.Y * UnclampedNewLambda;

        Ci = IsA ? Con.NextA : Con.NextB;
    }

    B.Scratch.NewPos = B.Pos + SolveSymmetric2x2(H, F);
}

void World::Step(float Dt)
{
    DoCollisionDetection();
    ComputeInertials(Dt);

    // warm start
    ConstraintIndex Ci = NextActiveConstraintIndexPlusOne;
    while (Ci != INVALID)
    {
        Constraint &Con = *GetConstraint(Ci);
        Con.Lambda *= Settings.Alpha * Settings.Gamma;
        //Con.K = std::max(Settings.KStart, Con.K * Settings.Gamma);
        Con.K = Con.K * Settings.Gamma;
        Ci = Con.NextA;
    }

    const float InvDt = 1.0f / Dt;
    const float InvDt2 = InvDt * InvDt;

    // set OldPos
    for (BodyIndex bi = 0; bi < NumBodies; ++bi)
    {
        Bodies[bi].Scratch.OldPos.X = Bodies[bi].Pos.X;
        Bodies[bi].Scratch.OldPos.Y = Bodies[bi].Pos.Y;
    }

    // iterate
    for (int iter = 0; iter < Settings.Iterations; ++iter)
    {
        for (BodyIndex bi = 0; bi < NumBodies; ++bi)
        {
            UpdateBody(this, Settings, bi + 1, InvDt2);
        }

        // update after to avoid races
        for (BodyIndex bi = 0; bi < NumBodies; ++bi)
        {
            Bodies[bi].Pos = Bodies[bi].Scratch.NewPos;
        }

        // update the constraints
        ConstraintIndex Ci = NextActiveConstraintIndexPlusOne;
        while (Ci != INVALID)
        {
            Constraint &Con = *GetConstraint(Ci);

            float C;
            Vector2 _CJacobian;
            Vector2 _Gd;
            float UnclampedNewLambda;
            ComputeConstraint(*this, Ci, Con.A, C, _CJacobian, _Gd, UnclampedNewLambda);

            Con.Lambda = Clamp(UnclampedNewLambda, Con.LambdaMin, Con.LambdaMax);
            if ((Con.Lambda > Con.LambdaMin) & (Con.Lambda < Con.LambdaMax))
            {
                Con.K += Settings.Beta * std::abs(C);
            }

            Ci = Con.Next;
        }
    }

    // update velocities
    for (BodyIndex bi = 0; bi < NumBodies; ++bi)
    {
        Bodies[bi].Vel = (Bodies[bi].Pos - Bodies[bi].Scratch.OldPos) * InvDt;
    }
}

void World::AddConstraint(BodyIndex A, BodyIndex B, ConstraintType Type, Vector2 Normal, float Rest, float K,
                          float LambdaMin, float LambdaMax)
{
    if (NextFreeConstraintIndexPlusOne == INVALID)
    {
        throw std::runtime_error("Too many bodies");
    }

    Constraint& New = Constraints[NextFreeConstraintIndexPlusOne - 1];
    New.A = A;
    New.B = B;
    New.Type = Type;
    New.Normal = Normal;
    New.Rest = Rest;
    New.K = K;
    New.Lambda = 0.f;
    New.LambdaMin = LambdaMin;
    New.LambdaMax = LambdaMax;

    // splice into linked lists
    const ConstraintIndex NextFree = New.Next;
    New.Next = NextActiveConstraintIndexPlusOne;
    NextActiveConstraintIndexPlusOne = NextFreeConstraintIndexPlusOne;
    NextFreeConstraintIndexPlusOne = NextFree;

    if (A != INVALID)
    {
        New.NextA = ConstraintHeads[A - 1];
        ConstraintHeads[A - 1] = NextActiveConstraintIndexPlusOne;
    }
    else
    {
        New.NextA = INVALID;
    }
    if (B != INVALID)
    {
        New.NextB = ConstraintHeads[B - 1];
        ConstraintHeads[B - 1] = NextActiveConstraintIndexPlusOne;
    }
    else
    {
        New.NextB = INVALID;
    }
}

void World::ComputeInertials(float Dt)
{
    const Vector2 ATSq = Vector2(0, Settings.Gravity * Dt * Dt);
    for (int i = 0; i < NumBodies; i++)
    {
        Body& Body = Bodies[i];
        Body.Pos += Body.Vel * Dt + ATSq;
        Body.Scratch.Y = Body.Pos;
    }
}

void World::DoCollisionDetection()
{
}
