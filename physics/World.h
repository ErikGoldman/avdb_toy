#pragma once

#include <iostream>
#include "BasicObjects.h"

struct WorldSettings
{
    float Gravity = -9.81f;
    float KStart = 1e4;
    int Iterations = 4;

    // how much to increase K each step
    float Beta = 10;

    float LambdaMin = 10, LambdaMax = 2e8;
};

class World final
{
private:
    static constexpr int MAX_BODIES = 3000;
    static constexpr int MAX_CONSTRAINTS = 5 * MAX_BODIES;

public:
    World(): Constraints{}, ConstraintHeads{}, Bodies{}
    {
        memset(ConstraintHeads, INVALID, sizeof(ConstraintHeads));
    }

    void Step(float Dt);

    void AddBody(const Vector2 &Pos, const float Radius, const float Mass)
    {
        if (NumBodies == MAX_BODIES)
        {
            throw std::runtime_error("Too many bodies");
        }

        Body &New = Bodies[NumBodies++];
        New.Pos = Pos;
        New.Radius = Radius;
        New.Mass = Mass;
        New.InvMass = 1.0f / Mass;
    }

    Body *GetBody(const BodyIndex Bi)
    {
        return &Bodies[Bi - 1];
    }
    const Body *GetBody(const BodyIndex Bi) const
    {
        return &Bodies[Bi - 1];
    }

    Constraint *GetConstraint(const ConstraintIndex Ci)
    {
        return &Constraints[Ci - 1];
    }
    const Constraint *GetConstraint(const ConstraintIndex Ci) const
    {
        return &Constraints[Ci - 1];
    }

    ConstraintIndex GetStartConstraint(const BodyIndex Bi) const
    {
        return ConstraintHeads[Bi];
    }

    void AddConstraint(BodyIndex A, BodyIndex B, ConstraintType Type, Vector2 Normal, float Rest, float K,
        float Lambda = 0.0f)
    {
        if (NextFreeConstraintIndexPlusOne == INVALID)
        {
            throw std::runtime_error("Too many bodies");
        }

        Constraint &New = Constraints[NextFreeConstraintIndexPlusOne - 1];
        New.A = A;
        New.B = B;
        New.Type = Type;
        New.Normal = Normal;
        New.Rest = Rest;
        New.K = K;
        New.Lambda = Lambda;

        // splice into linked lists
        const ConstraintIndex NextFree = New.Next;
        New.Next = NextActiveConstraintIndexPlusOne;
        NextActiveConstraintIndexPlusOne = NextFreeConstraintIndexPlusOne;
        NextFreeConstraintIndexPlusOne = NextFree;

        if (A != INVALID)
        {
            New.NextA = ConstraintHeads[A];
            ConstraintHeads[A] = NextActiveConstraintIndexPlusOne;
        }
        else
        {
            New.NextA = INVALID;
        }
        if (B != INVALID)
        {
            New.NextB = ConstraintHeads[B];
            ConstraintHeads[B] = NextActiveConstraintIndexPlusOne;
        }
        else
        {
            New.NextB = INVALID;
        }
    }

    int GetNumBodies() const
    {
        return NumBodies;
    }

    WorldSettings Settings;

private:
    void ComputeInertials(float Dt);
    void DoCollisionDetection();

    // constraint info
    Constraint Constraints[MAX_CONSTRAINTS];
    ConstraintIndex NextFreeConstraintIndexPlusOne = 1;
    ConstraintIndex NextActiveConstraintIndexPlusOne = 0;
    ConstraintIndex ConstraintHeads[MAX_BODIES];

    // can't remove bodies for now so just use NumBodies to keep track
    Body Bodies[MAX_BODIES];
    BodyIndex NumBodies = 0;
};
