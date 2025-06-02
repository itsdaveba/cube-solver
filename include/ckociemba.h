#ifndef CKOCIEMBA_H
#define CKOCIEMBA_H

void get_phase_coords(int *phase_coords, int phase, int *coords);

int get_phase_coords_len(int phase);

void get_coords_from_phase_coords(int *coords, int phase, int *phase_coords);

#endif