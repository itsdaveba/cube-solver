#ifndef DEFS_H
#define DEFS_H

#define NONE -1

#define NUM_CORNERS 8
#define NUM_EDGES 12
#define NUM_FACES 6
#define NUM_MOVES 18
#define NUM_THREADS 6
#define NUM_ORBIT_ELEMS 4

#define OP_SIZE 24
#define COMB_8C4_MINUS_ONE 69
#define COMB_12C4_MINUS_ONE 494

extern const int PERM[NUM_FACES][2][4];
extern const int ORIENT_DIFF[NUM_FACES][2][4];
extern const int FACTORIAL[NUM_EDGES + 1];
extern const int COMBINATION[NUM_EDGES + 1][NUM_ORBIT_ELEMS + 1];
extern const int CORNER_ORBITS[NUM_CORNERS];
extern const int EDGE_ORBITS[NUM_EDGES];
extern const int CORNER_THREAD[OP_SIZE][OP_SIZE];
extern const int TETRAD_M11_PERM[OP_SIZE][OP_SIZE / NUM_THREADS];
extern const int SLICE_PARITY[2][OP_SIZE];
extern const int SLICE_STANDING_PERM[2][OP_SIZE / 2];

#endif