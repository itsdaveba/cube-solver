#include "defs.h"
#include "thistlethwaite.h"

int thistle_get_phase_coords_size(int phase)
{
    if (phase == 0)
        return 1;
    if (phase == 1)
        return 2;
    if (phase == 2)
        return 3;
    if (phase == 3)
        return 5;
    return NONE;
}

void thistle_get_phase_coords(int *phase_coords, int *coords, int phase)
{
    if (phase == 0)
    {
        phase_coords[0] = coords[1];
    }
    else if (phase == 1)
    {
        phase_coords[0] = coords[0];
        phase_coords[1] = coords[5] / OP_SIZE;
    }
    else if (phase == 2)
    {
        phase_coords[0] = coords[2] / OP_SIZE;
        phase_coords[1] = CORNER_THREAD[coords[2] % OP_SIZE][coords[3] % OP_SIZE];
        phase_coords[2] = coords[4] / OP_SIZE;
    }
    else if (phase == 3)
    {
        phase_coords[0] = coords[2];
        phase_coords[1] = (coords[3] % OP_SIZE) / NUM_THREADS;
        phase_coords[2] = coords[4];
        phase_coords[3] = coords[5] % OP_SIZE;
        phase_coords[4] = (coords[6] % OP_SIZE) / 2;
    }
}

void thistle_set_coords(int *coords, int *phase_coords, int phase)
{
    if (phase == 0)
    {
        coords[1] = phase_coords[0];
    }
    else if (phase == 1)
    {
        coords[0] = phase_coords[0];
        coords[5] = phase_coords[1] * OP_SIZE;
    }
    else if (phase == 2)
    {
        coords[2] = phase_coords[0] * OP_SIZE;
        coords[3] = (COMB_8C4_MINUS_ONE - phase_coords[0]) * OP_SIZE + phase_coords[1];
        coords[4] = phase_coords[2] * OP_SIZE;
    }
    else if (phase == 3)
    {
        coords[2] = phase_coords[0];
        coords[3] = COMB_8C4_MINUS_ONE * OP_SIZE + TETRAD_M11_PERM[phase_coords[0]][phase_coords[1]];
        coords[4] = phase_coords[2];
        coords[5] = COMB_12C4_MINUS_ONE * OP_SIZE + phase_coords[3];
        coords[6] = COMB_8C4_MINUS_ONE * OP_SIZE + SLICE_STANDING_PERM[SLICE_PARITY[SLICE_PARITY[0][phase_coords[2]]][phase_coords[3]]][phase_coords[4]];
    }
}
