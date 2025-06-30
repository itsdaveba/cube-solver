#ifndef KORF_H
#define KORF_H

int korf_get_phase_coords_size(int phase);
void korf_get_phase_coords(int *phase_coords, int *coords, int phase);
void korf_set_coords(int *coords, int *phase_coords, int phase);

#endif
