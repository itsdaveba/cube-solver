#include "defs.h"
#include "kociemba.h"

int kociemba_get_phase_coords_size(int phase)
{
    if (phase == 0)
        return 3;
    if (phase == 1)
        return 3;
    return NONE;
}

void kociemba_get_phase_coords(int *phase_coords, int *coords, int phase)
{
    if (phase == 0)
    {
        phase_coords[0] = coords[0];
        phase_coords[1] = coords[1];
        phase_coords[2] = coords[4] / OP_SIZE;
    }
    else if (phase == 1)
    {
        phase_coords[0] = coords[2];
        phase_coords[1] = coords[5] + (coords[3] - COMB_8C4_MINUS_ONE + coords[3] / OP_SIZE) * OP_SIZE;
        phase_coords[2] = coords[4] % OP_SIZE;
    }
}

void kociemba_set_coords(int *coords, int *phase_coords, int phase)
{
    if (phase == 0)
    {
        coords[0] = phase_coords[0];
        coords[1] = phase_coords[1];
        coords[4] = phase_coords[2] * OP_SIZE;
    }
    else if (phase == 1)
    {
        coords[2] = phase_coords[0];
        coords[3] = phase_coords[1] / OP_SIZE;
        coords[5] = phase_coords[1] + (COMB_8C4_MINUS_ONE - coords[3] - coords[3] / OP_SIZE) * OP_SIZE;
        coords[4] = COMB_12C4_MINUS_ONE * OP_SIZE + phase_coords[2];
    }
}
