#ifndef THISTLE_H
#define THISTLE_H

int thistle_get_phase_coords_size(int phase);
void thistle_get_phase_coords(int *phase_coords, int *coords, int phase);
void thistle_set_coords(int *coords, int *phase_coords, int phase);

#endif
