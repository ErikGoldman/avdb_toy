#pragma once

#include "BasicObjects.h"

class World final
{
public:

private:
    static constexpr int MAX_BODIES = 3000;

    // can't remove bodies for now so just use NumBodies to keep track
    Body Bodies[MAX_BODIES];
    int NumBodies;
};
