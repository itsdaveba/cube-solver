#include "korf.h"

int korf_get_phase_coords_size(int phase)
{
    return 6;
}

void korf_get_phase_coords(int *phase_coords, int *coords, int phase)
{
    for (int i = 0; i < 6; i++)
        phase_coords[i] = coords[i];
}

void korf_set_coords(int *coords, int *phase_coords, int phase)
{
    for (int i = 0; i < 6; i++)
        coords[i] = phase_coords[i];
}
