#pragma once

#include <iostream>
#include <random>

#include "BasicObjects.h"

struct WorldSettings
{
    float Gravity = -9.81f;
    float KStart = 1e4;
    int Iterations = 4;

    float CollisionMargin = 0.01f;

    // how much to increase K each step
    float Beta = 10.f;

    // how much to decrease k during a warm start
    float Gamma = 0.99f;

    // how much extra to decrease lambda during a warm start
    float Alpha = 0.95f;

    float LambdaMin = 1000.f, LambdaMax = 1000000000.f;
};

class World final
{
private:
    static constexpr int MAX_BODIES = 3000;
    static constexpr int MAX_PLANES = 5;
    static constexpr int MAX_PLANE_COLLISIONS = 30;
    static constexpr int MAX_BODY_COLLISIONS = 10;
    static constexpr int MAX_CONSTRAINTS = 5 * MAX_BODIES;

public:
    World(): Constraints{}, ConstraintHeads{}, Bodies{}, Planes{}
    {
        memset(ConstraintHeads, INVALID, sizeof(ConstraintHeads));
        memset(PlaneCollisions, INVALID, sizeof(PlaneCollisions));
        memset(BodyCollisions, INVALID, sizeof(BodyCollisions));

        for (ConstraintIndex i = 0; i < MAX_CONSTRAINTS - 1; ++i)
        {
            Constraints[i].Next = i + 1;
        }
        Constraints[MAX_CONSTRAINTS - 1].Next = INVALID;
    }

    void Step(float Dt);

    void SetRandomColor(Body &b)
    {
        static std::mt19937                    rng{std::random_device{}()};
        static std::uniform_real_distribution<>  pickH(0.0, 360.0);
        static std::uniform_real_distribution<>  pickS(0.55, 0.90);
        static std::uniform_real_distribution<>  pickV(0.75, 1.00);

        double H = pickH(rng);    // degrees
        double S = pickS(rng);
        double V = pickV(rng);

        double C = V * S;
        double Hp = H / 60.0;
        double X = C * (1.0 - std::fabs(std::fmod(Hp, 2.0) - 1.0));

        double r1, g1, b1;
        switch (static_cast<int>(std::floor(Hp)) % 6)
        {
        case 0: r1=C; g1=X; b1=0; break;
        case 1: r1=X; g1=C; b1=0; break;
        case 2: r1=0; g1=C; b1=X; break;
        case 3: r1=0; g1=X; b1=C; break;
        case 4: r1=X; g1=0; b1=C; break;
        default:r1=C; g1=0; b1=X; break;   // case 5
        }

        double m = V - C;
        b.Color[0] = static_cast<uint8_t>(std::round((r1 + m) * 255.0));
        b.Color[1] = static_cast<uint8_t>(std::round((g1 + m) * 255.0));
        b.Color[2] = static_cast<uint8_t>(std::round((b1 + m) * 255.0));
    }

    void AddBody(const Vector2 &Pos, const Vector2 &Vel, const float Radius, const float Mass)
    {
        if (NumBodies == MAX_BODIES)
        {
            throw std::runtime_error("Too many bodies");
        }

        Body &New = Bodies[NumBodies++];
        New.Pos = Pos;
        New.Vel = Vel;
        New.Radius = Radius;
        New.Mass = Mass;
        New.InvMass = 1.0f / Mass;
        SetRandomColor(New);
    }

    void AddPlane(const Vector2 &Normal, const float Distance)
    {
        if (NumPlanes == MAX_PLANES)
        {
            throw std::runtime_error("Too many planes");
        }

        Plane &New = Planes[NumPlanes++];
        New.Normal = Normal;
        New.Distance = Distance;
    }

    Body *GetBody(const BodyIndex Bi)
    {
        return &Bodies[Bi];
    }
    const Body *GetBody(const BodyIndex Bi) const
    {
        return &Bodies[Bi];
    }

    Constraint *GetConstraint(const ConstraintIndex Ci)
    {
        return &Constraints[Ci];
    }
    const Constraint *GetConstraint(const ConstraintIndex Ci) const
    {
        return &Constraints[Ci];
    }

    ConstraintIndex GetStartConstraint(const BodyIndex Bi) const
    {
        return ConstraintHeads[Bi];
    }

    ConstraintIndex AddConstraint(BodyIndex A, BodyIndex B, ConstraintType Type, Vector2 Normal, float Rest, float K,
        float FMin, float FMax);
    void RemoveConstraint(ConstraintIndex Index);

    BodyIndex GetNumBodies() const
    {
        return NumBodies;
    }

    BodyIndex GetNumPlanes() const
    {
        return NumPlanes;
    }
    Plane &GetPlane(const BodyIndex Index)
    {
        return Planes[Index];
    }

    WorldSettings Settings;

private:
    void DoCollisionDetection();
    void PrepareBodiesForStep(float Dt, float InvDt);
    void PrepareConstraintsForStep();
    bool PrepareConstraintForStep(Constraint& Con) const;
    ConstraintIndex RemoveConstraint(ConstraintIndex Index, ConstraintIndex PrevCi);

    void RemoveFromLinkedList(ConstraintIndex Index, BodyIndex Body);

    // constraint info
    Constraint Constraints[MAX_CONSTRAINTS];
    ConstraintIndex NextFreeConstraintIndex = 0;
    ConstraintIndex NextActiveConstraintIndex = INVALID;
    ConstraintIndex ConstraintHeads[MAX_BODIES];

    // can't remove bodies for now so just use NumBodies to keep track
    Body Bodies[MAX_BODIES];
    BodyIndex NumBodies = 0;

    Plane Planes[MAX_PLANES];
    BodyIndex NumPlanes = 0;

    // collision
    BodyIndex PlaneCollisions[MAX_PLANES][MAX_PLANE_COLLISIONS];
    BodyIndex BodyCollisions[MAX_BODIES][MAX_BODY_COLLISIONS];
};
