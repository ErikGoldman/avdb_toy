#include "World.h"

static _forceinline float Clamp(float x, float lo, float hi)
{
    return std::max(lo, std::min(x, hi));
}

void ComputeConstraint(World *World, ConstraintIndex Ci, BodyIndex ForBi, float &C, Vector2 &Jacobian, Vector2 &Gd)
{
    const Constraint &Con = *World->GetConstraint(Ci);
    switch (Con.Type)
    {
        case Normal:
        {

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
        ComputeConstraint(World, Ci, Bi, C, CJacobian, Gd);

        const float NewLambda = Con.K * C + Con.Lambda;

        //F -= CJacobian * Clamp(NewLambda, Settings.LambdaMin, Settings.LambdaMax);
        F -= CJacobian * NewLambda;
        H += CJacobian.OuterProduct(CJacobian) * Con.K;
        H.Values[0] += Gd.X * NewLambda;
        H.Values[3] += Gd.Y * NewLambda;

        Ci = IsA ? Con.NextA : Con.NextB;
    }

    B.Scratch.NewPos = B.Pos + SolveSymmetric2x2(H, F);
}

void World::Step(float Dt)
{
    DoCollisionDetection();
    ComputeInertials(Dt);

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
            Vector2 CJacobian;
            Vector2 Gd;
            ComputeConstraint(this, Ci, true, C, CJacobian, Gd);

            Con.Lambda = Clamp(Con.Lambda + Con.K * C, Settings.LambdaMin, Settings.LambdaMax);
            if (Con.Lambda < Settings.LambdaMax && Con.Lambda > Settings.LambdaMin)
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

void World::ComputeInertials(float Dt)
{
    const Vector2 FTSq = Vector2(0, Settings.Gravity * Dt * Dt);
    for (int i = 0; i < NumBodies; i++)
    {
        Body& Body = Bodies[i];
        Body.Scratch.Y = Body.Pos + Body.Vel * Dt + FTSq * Body.InvMass;
    }
}
