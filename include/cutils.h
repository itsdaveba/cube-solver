#ifndef CUTILS_H
#define CUTILS_H

#define EMPTY -1
#define NUM_FACES 6
#define NUM_MOVES 18
#define NUM_CORNERS 8
#define NUM_EDGES 12
#define NUM_AXIS_ELEMS 4
#define NUM_AXIS_ELEMS_PERM 24

extern const int CORNER_AXES[NUM_CORNERS];
extern const int EDGE_AXES[NUM_EDGES];
extern const int CORNER_AXIS_OFFSET[2];
extern const int EDGE_AXIS_OFFSET[3];
extern const int MOVES[NUM_FACES][2][4];
extern const int ORIENTATION[NUM_FACES][2][4];

int get_orientation_coord(int *array, int v, int n);
int get_permutation_coord(int *array, int n);
int get_partial_permutation_coord(int *array, int n, const int *axes, int axis);

void set_orientation_coord(int *array, int coord, int v, int n);
void set_permutation_coord(int *array, int coord, int n);
void set_partial_permutation_coord(int *array, int coord, int n, int axis_offset);

#endif
