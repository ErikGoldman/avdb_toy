#pragma once
#include <cstdint>

typedef uint16_t BodyIndex;
typedef uint16_t ConstraintIndex;
static constexpr uint16_t INVALID = ~0;

enum ConstraintType
{
    Normal
};
struct Matrix2x2
{
    float Values[4];

    Matrix2x2()
    {
        memset(Values, 0, sizeof(Values));
    }

    Matrix2x2(float V1, float V2, float V3, float V4)
    {
        Values[0] = V1;
        Values[1] = V2;
        Values[2] = V3;
        Values[3] = V4;
    }

    explicit Matrix2x2(float Diagonal)
    {
        Values[0] = Diagonal;
        Values[1] = 0.f;
        Values[2] = 0.f;
        Values[3] = Diagonal;
    }

    Matrix2x2 &operator+=(const Matrix2x2 &Other)
    {
        Values[0] += Other.Values[0];
        Values[1] += Other.Values[1];
        Values[2] += Other.Values[2];
        Values[3] += Other.Values[3];
        return *this;
    }

    Matrix2x2 operator+(const Matrix2x2 &Other) const
    {
        Matrix2x2 Out;
        Out.Values[0] = Values[0] + Other.Values[0];
        Out.Values[1] = Values[1] + Other.Values[1];
        Out.Values[2] = Values[2] + Other.Values[2];
        Out.Values[3] = Values[3] + Other.Values[3];
        return Out;
    }

    Matrix2x2 operator*(float Other) const
    {
        Matrix2x2 Out;
        Out.Values[0] = Values[0] * Other;
        Out.Values[1] = Values[1] * Other;
        Out.Values[2] = Values[2] * Other;
        Out.Values[3] = Values[3] * Other;
        return Out;
    }
};

struct Vector2
{
    float X, Y;

    Vector2()
    : X(0), Y(0)
    {}

    Vector2(float InX, float InY)
    : X(InX), Y(InY)
    {}

    Vector2 operator+(const Vector2 &Other) const
    {
        return {X + Other.X, Y + Other.Y};
    }

    Vector2 &operator+=(const Vector2 &Other)
    {
        X += Other.X;
        Y += Other.Y;
        return *this;
    }
    Vector2 &operator-=(const Vector2 &Other)
    {
        X -= Other.X;
        Y -= Other.Y;
        return *this;
    }

    Vector2 operator-(const Vector2 &Other) const
    {
        return {X - Other.X, Y - Other.Y};
    }

    Vector2 operator*(const float F) const
    {
        return {X * F, Y * F};
    }

    Vector2 GetNormal() const
    {
        const float length = sqrtf(Dot(*this));
        return { X / length, Y / length };
    }

    float Dot(const Vector2 &Other) const
    {
        return X * Other.X + Y * Other.Y;
    }

    float Length() const
    {
        return sqrtf(X * X + Y * Y);
    }

    float SqLength() const
    {
        return X * X + Y * Y;
    }

    Matrix2x2 OuterProduct(const Vector2 &Other) const
    {
        return {
            X * Other.X, X * Other.Y,
            Y * Other.X, Y * Other.Y
        };
    }
};



struct BodyScratch
{
    // inertial
    Vector2 Inertial;

    // old pos (used for velocity)
    Vector2 OldPos;

    // old velocity (used to warm start the body with an acceleration guess)
    Vector2 OldVel;
};

struct Body
{
    Vector2 Pos, Vel;
    float Radius;
    float Mass;
    float InvMass;
    uint8_t Color[3];

    BodyScratch Scratch;
};

struct Plane
{
    Vector2 Normal;
    float Distance;
};

struct Constraint
{
    BodyIndex A, B, Plane;
    ConstraintIndex NextA, NextB, Next;

    ConstraintType Type;
    float Rest;
    float K, Lambda;
    float FMin, FMax;

    // "scratch" data precomputed during the initialize step
    Vector2 Normal, JA, JB;
    float C0;
};
