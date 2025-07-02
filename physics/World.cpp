#include "World.h"

#include <algorithm>
#include <SDL_log.h>

static _forceinline float Clamp(float x, float lo, float hi)
{
    return std::max(lo, std::min(x, hi));
}

static void ComputeConstraint(const World &World, ConstraintIndex Ci, BodyIndex ForBi, float &C, Vector2 &Jacobian)
{
    const Constraint &con = *World.GetConstraint(Ci);
    switch (con.Type)
    {
        case Normal:
        {
            const Body &a = *World.GetBody(con.A);
            const Vector2 dPosA = a.Pos - a.Scratch.OldPos;

            // first two terms of the taylor series
            C = con.C0 * (1 - World.Settings.Alpha) + con.JA.Dot(dPosA);

            if (con.B != INVALID)
            {
                const Body &b = *World.GetBody(con.B);
                const Vector2 dPosB = b.Pos - b.Scratch.OldPos;
                C += con.JB.Dot(dPosB);
            }

            Jacobian = ForBi == con.A ? con.JA : con.JB;
        }
        break;

        default:
            throw std::runtime_error("Unknown Constraint Type");
    }
}

// Solve A x = b for a symmetric 2×2 matrix
static _forceinline Vector2 SolveSymmetric2x2(const Matrix2x2& A, const Vector2& b)
{
    /* --- LDLᵀ factorisation --------------------------------------- */
    const float D1  = A.Values[0];                  // pivot 1
    const float L21 = A.Values[2] / D1;          // A21 / D1
    const float D2  = A.Values[3] - L21*L21*D1;     // Schur complement

    /* --- Forward substitution  (Ly = b) --------------------------- */
    const float y1 = b.X;                       // y1 = b1
    const float y2 = b.Y - L21*y1;              // y2 = b2 – L21·y1

    /* --- Diagonal solve        (Dz = y) --------------------------- */
    const float z1 = y1 / D1;
    const float z2 = y2 / D2;

    /* --- Back substitution     (Lᵀx = z) -------------------------- */
    Vector2 x;
    x.Y = z2;                             // x2 = z2
    x.X = z1 - L21*x.Y;                   // x1 = z1 – L21·x2

    return x;                             // solution of A x = b
}


void UpdateBody(World *World, WorldSettings &Settings, BodyIndex Bi, float InvDt2)
{
    Body &B = *World->GetBody(Bi);
    Vector2 F = (B.Pos - B.Scratch.Inertial) * B.Mass * InvDt2;
    Matrix2x2 H(B.Mass * InvDt2);

    ConstraintIndex Ci = World->GetStartConstraint(Bi);
    while (Ci != INVALID)
    {
        const Constraint &Con = *World->GetConstraint(Ci);
        const bool IsA = (Con.A == Bi);

        float C;
        Vector2 CJacobian;
        ComputeConstraint(*World, Ci, Bi, C, CJacobian);

        //SDL_Log("C = %.6f   CJ = %.6f, %.6f", C, CJacobian.X, CJacobian.Y);

        const float ConstraintFMag = Clamp(Con.K * C + Con.Lambda, Con.FMin, Con.FMax);
        F += CJacobian * ConstraintFMag;

        // TODO: this formula ignores G (sec 3.5)
        H += CJacobian.OuterProduct(CJacobian * Con.K);

        /* G only needed for soft springs
        Matrix2x2 G(length(force->H[i].col(0)), length(force->H[i].col(1)), length(force->H[i].col(2)) * std::abs(ConstraintFMag));
        H.Values[0] += Gd.X * UnclampedNewLambda;
        H.Values[3] += Gd.Y * UnclampedNewLambda;
        */

        Ci = IsA ? Con.NextA : Con.NextB;
    }

    B.Pos -= SolveSymmetric2x2(H, F);
}

void World::Step(float Dt)
{
    const float InvDt = 1.0f / Dt;
    const float InvDt2 = InvDt * InvDt;

    DoCollisionDetection();

    PrepareBodiesForStep(Dt, InvDt);
    PrepareConstraintsForStep();

    // iterate
    for (int iter = 0; iter < Settings.Iterations; ++iter)
    {
        // primal (position) update
        for (BodyIndex bi = 0; bi < NumBodies; ++bi)
        {
            UpdateBody(this, Settings, bi, InvDt2);
        }

        // dual (constraint) update
        ConstraintIndex Ci = NextActiveConstraintIndex;
        while (Ci != INVALID)
        {
            Constraint &Con = *GetConstraint(Ci);

            float C;
            Vector2 _CJacobian;
            ComputeConstraint(*this, Ci, Con.A, C, _CJacobian);

            Con.Lambda = Clamp(Con.Lambda + Con.K * C, Con.FMin, Con.FMax);
            if ((Con.Lambda > Con.FMin) & (Con.Lambda < Con.FMax))
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

ConstraintIndex World::AddConstraint(BodyIndex A, BodyIndex B, ConstraintType Type, Vector2 Normal, float Rest, float K,
                                     float FMin, float FMax)
{
    if (NextFreeConstraintIndex == INVALID)
    {
        throw std::runtime_error("Too many bodies");
    }

    ConstraintIndex CreatedIndex = NextFreeConstraintIndex;
    Constraint& New = Constraints[NextFreeConstraintIndex];
    New.A = A;
    New.B = B;
    New.Plane = INVALID;
    New.Type = Type;
    New.Normal = Normal;
    New.Rest = Rest;
    New.K = K;
    New.Lambda = 0.f;
    New.FMin = FMin;
    New.FMax = FMax;

    // splice into linked lists
    const ConstraintIndex NextFree = New.Next;
    New.Next = NextActiveConstraintIndex;
    NextActiveConstraintIndex = NextFreeConstraintIndex;
    NextFreeConstraintIndex = NextFree;

    if (A != INVALID)
    {
        New.NextA = ConstraintHeads[A];
        ConstraintHeads[A] = CreatedIndex;
    }
    else
    {
        New.NextA = INVALID;
    }
    if (B != INVALID)
    {
        New.NextB = ConstraintHeads[B];
        ConstraintHeads[B] = CreatedIndex;
    }
    else
    {
        New.NextB = INVALID;
    }

    return CreatedIndex;
}

void World::RemoveConstraint(ConstraintIndex Index)
{
    for (ConstraintIndex i = NextActiveConstraintIndex; i != INVALID; i = Constraints[i].Next)
    {
        if (Constraints[i].Next == Index)
        {
            RemoveConstraint(Index, i);
            return;
        }
    }

    throw std::runtime_error("Constraint doesn't exist");
}

ConstraintIndex World::RemoveConstraint(ConstraintIndex Index, ConstraintIndex PrevCi)
{
    const ConstraintIndex NextActiveCi = Constraints[Index].Next;

    // remove from body linked lists
    BodyIndex ai = Constraints[Index].A;
    if (ai != INVALID)
    {
        RemoveFromLinkedList(Index, ai);
    }
    BodyIndex bi = Constraints[Index].B;
    if (bi != INVALID)
    {
        RemoveFromLinkedList(Index, bi);
    }

    // remove collision registration
    if (Constraints[Index].Type == Normal)
    {
        if (Constraints[Index].B != INVALID)
        {
            for (int i = 0; i < MAX_BODY_COLLISIONS; ++i)
            {
                BodyIndex &bodyIndex = BodyCollisions[ai][i];
                if (bodyIndex == bi)
                {
                    bodyIndex = INVALID;
                    break;
                }
            }
        }
        else if (Constraints[Index].Plane != INVALID)
        {
            for (int i = 0; i < MAX_PLANE_COLLISIONS; i++)
            {
                BodyIndex &bodyIndex = PlaneCollisions[Constraints[Index].Plane][i];
                if (bodyIndex == ai)
                {
                    bodyIndex = INVALID;
                    break;
                }
            }
        }
    }

    if (PrevCi == INVALID)
    {
        NextActiveConstraintIndex = NextActiveCi;
    }
    else
    {
        Constraints[PrevCi].Next = NextActiveCi;
    }

    Constraints[Index].Next = NextFreeConstraintIndex;
    NextFreeConstraintIndex = Index;

    Constraints[Index].NextA = INVALID;
    Constraints[Index].NextB = INVALID;
    Constraints[Index].Plane = INVALID;

    return NextActiveCi;
}

void World::RemoveFromLinkedList(ConstraintIndex Index, BodyIndex Body)
{
    ConstraintIndex *iterIndex = &ConstraintHeads[Body];
    while (*iterIndex != INVALID)
    {
        ConstraintIndex *nextIndex = Constraints[*iterIndex].A == Body ? &Constraints[*iterIndex].NextA : &Constraints[*iterIndex].NextB;
        if (*iterIndex == Index)
        {
            *iterIndex = *nextIndex;
            break;
        }
        iterIndex = nextIndex;
    }
}

void World::DoCollisionDetection()
{
    // planes
    for (BodyIndex p = 0; p < NumPlanes; p++)
    {
        const float constraintDistance = Planes[p].Distance + Settings.CollisionMargin;
        for (BodyIndex b = 0; b < NumBodies; b++)
        {
            int CollideIndex = -1;
            for (int i = 0; i < MAX_PLANE_COLLISIONS; i++)
            {
                if (PlaneCollisions[p][i] == b)
                {
                    CollideIndex = i;
                    break;
                }
            }

            const float distance = Bodies[b].Pos.Dot(Planes[p].Normal) - Bodies[b].Radius;
            if (distance < constraintDistance)
            {
                if (CollideIndex == -1)
                {
                    for (int i = 0; i < MAX_PLANE_COLLISIONS; i++)
                    {
                        if (PlaneCollisions[p][i] == INVALID)
                        {
                            PlaneCollisions[p][i] = b;

                            const ConstraintIndex ci = AddConstraint(b, INVALID, Normal, Planes[p].Normal,
                                Planes[p].Distance, Settings.KStart, -1e8, 0);
                            // hack: hide the plane index in NextB (B is invalid so it's never used)
                            Constraints[ci].Plane = p;
                            break;
                        }
                    }
                }
            }
        }
    }

    // bodies. start with no broad phase for now
    for (BodyIndex b1 = 0; b1 < NumBodies; b1++)
    {
        for (BodyIndex b2 = b1 + 1; b2 < NumBodies; b2++)
        {
            const float sqDistance = (Bodies[b1].Pos - Bodies[b2].Pos).SqLength();
            const float collisionDistance = Bodies[b1].Radius + Bodies[b2].Radius + Settings.CollisionMargin;
            if (sqDistance < collisionDistance * collisionDistance)
            {
                bool AlreadyCollided = false;
                for (int i = 0; i < MAX_BODY_COLLISIONS; i++)
                {
                    if (BodyCollisions[b1][i] == b2)
                    {
                        AlreadyCollided = true;
                        break;
                    }
                }

                if (!AlreadyCollided)
                {
                    for (int i = 0; i < MAX_BODY_COLLISIONS; i++)
                    {
                        if (BodyCollisions[b1][i] == INVALID)
                        {
                            AddConstraint(b1, b2, Normal, Vector2(0, 0),
                                Bodies[b1].Radius + Bodies[b2].Radius, Settings.KStart, -1e8, 0);
                            BodyCollisions[b1][i] = b2;
                            break;
                        }
                    }
                }
            }
        }
    }
}

void World::PrepareBodiesForStep(const float Dt, const float InvDt)
{
    // arbitrary restriction from our acceleration estimation code. can disable that if needed
    if (abs(Settings.Gravity) < 0.01)
    {
        throw std::runtime_error("Gravity too small");
    }

    const float DtSq = Dt * Dt;
    const Vector2 ATSq = Vector2(0, Settings.Gravity * DtSq);
    for (int i = 0; i < NumBodies; i++)
    {
        Body& Body = Bodies[i];

        // compute inertial. this is where the body would go if unconstrained
        Body.Scratch.Inertial = Body.Pos + Body.Vel * Dt + ATSq;

        // warm start the body positions by extrapolating from the last acceleration. this saves us some iteration time
        // spent just following the inertial since it's probably going to lead here anyway
        const Vector2 LastAcceleration = (Body.Vel - Body.Scratch.OldVel) * InvDt;

        // LastAcceleration is the sum of the external forces (i.e. gravity) and any constraint forces that pushed
        // us. measure it just in the direction of external forces to estimate how much those forces will probably
        // accelerate us again in practice (e.g. if we're on a surface, gravity will accelerate us 0)
        // TODO: obviously I don't need to recompute this every time/this is too generalized. keeping for readability
        const Vector2 UnconstrainedAccelExt(0, Settings.Gravity);
        const float ActualAccelExt = LastAcceleration.Dot(UnconstrainedAccelExt.GetNormal());
        const float AccelWeight = Clamp(ActualAccelExt / Settings.Gravity, 0.0f, 1.0f);

        // save old pos and compute warm started position
        Body.Scratch.OldPos = Body.Pos;
        Body.Pos += Body.Vel * Dt + UnconstrainedAccelExt * (AccelWeight * DtSq);
    }
}

void World::PrepareConstraintsForStep()
{
    ConstraintIndex PrevCi = INVALID;
    ConstraintIndex Ci = NextActiveConstraintIndex;

    while (Ci != INVALID)
    {
        Constraint& Con = *GetConstraint(Ci);
        if (!PrepareConstraintForStep(Con))
        {
            Ci = RemoveConstraint(Ci, PrevCi);
        }
        else
        {
            // warm start parameters
            Con.Lambda *= Settings.Alpha * Settings.Gamma;
            Con.K = std::max(Con.K * Settings.Gamma, Settings.KStart);
            // TODO: soft constraint clamping logic goes here (don't exceed soft constraint stiffness)

            PrevCi = Ci;
            Ci = Con.Next;
        }
    }
}

bool World::PrepareConstraintForStep(Constraint& Con) const
{
    switch (Con.Type)
    {
    case Normal:
        {
            Vector2 PosB;
            float BRadius;
            if (Con.A != INVALID && Con.B != INVALID)
            {
                Con.Normal = (Bodies[Con.A].Pos - Bodies[Con.B].Pos).GetNormal();
                PosB = Bodies[Con.B].Pos;
                BRadius = Bodies[Con.B].Radius;
            }
            else
            {
                PosB = Vector2(Con.Normal * Con.Rest);
                BRadius = 0;
            }

            // TODO: if I'm not using shapes that are rotationally symmetric then I need to transform things here
            Con.C0 = (Bodies[Con.A].Pos - PosB).Dot(Con.Normal) - Bodies[Con.A].Radius - BRadius - Settings.CollisionMargin;
            Con.JA = Con.Normal;
            Con.JB = Con.Normal * -1;

            // the constraint should remove itself if there's no penetration
            return Con.C0 < 0.0f;
        }
    default:
        throw std::runtime_error("Invalid constraint type");
    }
}
