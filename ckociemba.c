#include <stdlib.h>

void get_phase_coords(int *phase_coords, int phase, int *coords)
{
    if (phase == 0)
    {
        phase_coords[0] = coords[0];
        phase_coords[1] = coords[1];
        phase_coords[2] = coords[3] / 24;
    }
    else if (phase == 1)
    {
        phase_coords[0] = coords[2];
        phase_coords[1] = coords[4] + (coords[5] + coords[5] / 24 - 69) * 24;
        phase_coords[2] = coords[3] % 24;
    }
}

int get_phase_coords_len(int phase)
{
    return 3;
}

void get_coords_from_phase_coords(int *coords, int phase, int *phase_coords)
{
    if (phase == 0)
    {
        coords[0] = phase_coords[0];
        coords[1] = phase_coords[1];
        coords[3] = phase_coords[2] * 24;
    }
    else if (phase == 1)
    {
        coords[2] = phase_coords[0];
        coords[3] = 11856 + phase_coords[2];
        coords[4] = (24 * (69 - phase_coords[1] / 576)) + phase_coords[1] % 24;
        coords[5] = phase_coords[1] / 24;
    }
}
