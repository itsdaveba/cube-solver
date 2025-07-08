#ifndef KOCIEMBA_H
#define KOCIEMBA_H

int kociemba_get_phase_coords_size(int phase);
void kociemba_get_phase_coords(int *phase_coords, int *coords, int phase);
void kociemba_set_coords(int *coords, int *phase_coords, int phase);

#endif
