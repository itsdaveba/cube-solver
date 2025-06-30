#ifndef CUTILS_H
#define CUTILS_H

int get_orientation_coord(int *array, int v, int n);
int get_permutation_coord(int *array, int n);
int get_combination_coord(int *array, int n);
int get_partial_permutation_coord(int *array, const int *orbits, int orbit, int n);

void set_orientation_coord(int *array, int coord, int v, int n);
void set_permutation_coord(int *array, int coord, int n);
void set_combination_coord(int *array, int coord, int n);
void set_partial_permutation_coord(int *array, int coord, int n);

npy_intp get_table_index(int *phase_coords, int *indexes, npy_intp *dims, int nd);
void set_phase_coords(int *phase_coords, npy_intp index, int *indexes, npy_intp *dims, int nd);

#endif
