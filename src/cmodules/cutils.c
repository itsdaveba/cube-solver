// #include "cutils.h"

// const int CORNER_AXES[NUM_CORNERS] = {0, 0, 0, 0, 1, 1, 1, 1};
// const int EDGE_AXES[NUM_EDGES] = {2, 2, 2, 2, 1, 1, 1, 1, 0, 0, 0, 0};
// const int CORNER_AXIS_OFFSET[2] = {0, 4};
// const int EDGE_AXIS_OFFSET[3] = {8, 4, 0};

// const int MOVES[NUM_FACES][2][4] = {
//     {{0, 5, 1, 4}, {0, 4, 1, 5}},
//     {{1, 5, 3, 7}, {1, 10, 3, 11}},
//     {{1, 7, 2, 4}, {5, 11, 7, 9}},
//     {{2, 7, 3, 6}, {2, 7, 3, 6}},
//     {{0, 4, 2, 6}, {0, 9, 2, 8}},
//     {{0, 6, 3, 5}, {4, 8, 6, 10}}};

// const int ORIENTATION[NUM_FACES][2][4] = {
//     {{0, 0, 0, 0}, {0, 0, 0, 0}},
//     {{1, 2, 1, 2}, {1, 1, 1, 1}},
//     {{2, 1, 2, 1}, {0, 0, 0, 0}},
//     {{0, 0, 0, 0}, {0, 0, 0, 0}},
//     {{1, 2, 1, 2}, {1, 1, 1, 1}},
//     {{2, 1, 2, 1}, {0, 0, 0, 0}}};

// static const int COMBINATION[NUM_AXIS_ELEMS][NUM_EDGES] = {
//     {0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11},
//     {0, 0, 1, 3, 6, 10, 15, 21, 28, 36, 45, 55},
//     {0, 0, 0, 1, 4, 10, 20, 35, 56, 84, 120, 165},
//     {0, 0, 0, 0, 1, 5, 15, 35, 70, 126, 210, 330}};

// int get_orientation_coord(int *array, int v, int n)
// {
//     int coord = 0;
//     for (int i = 0; i < n - 1; i++)
//     {
//         coord *= v;
//         coord += array[i];
//     }
//     return coord;
// }

// int get_permutation_coord(int *array, int n)
// {
//     int coord = 0;
//     for (int i = 0; i < n - 1; i++)
//     {
//         coord *= n - i;
//         for (int j = i + 1; j < n; j++)
//         {
//             if (array[i] > array[j])
//             {
//                 coord += 1;
//             }
//         }
//     }
//     return coord;
// }

// int get_partial_permutation_coord(int *array, int n, const int *axes, int axis)
// {
//     int i = 0;
//     int comb[NUM_AXIS_ELEMS];
//     int perm[NUM_AXIS_ELEMS];
//     for (int c = 0; c < n; c++)
//     {
//         if (array[c] != -1 && axes[array[c]] == axis)
//         {
//             comb[i] = c;
//             perm[i++] = array[c];
//             if (i == NUM_AXIS_ELEMS)
//             {
//                 break;
//             }
//         }
//     }

//     int comb_coord = 0;
//     for (int i = 0; i < NUM_AXIS_ELEMS; i++)
//     {
//         comb_coord += COMBINATION[i][comb[i]];
//     }
//     int perm_coord = get_permutation_coord(perm, NUM_AXIS_ELEMS);

//     return comb_coord * NUM_AXIS_ELEMS_PERM + perm_coord;
// }

// void set_orientation_coord(int *array, int coord, int v, int n)
// {
//     int o = 0;
//     for (int i = n - 2; i >= 0; i--)
//     {
//         array[i] = coord % v;
//         o -= array[i];
//         coord = coord / v;
//     }
//     o = o % v;
//     array[n - 1] = o < 0 ? o + v : o;
// }

// void set_permutation_coord(int *array, int coord, int n)
// {
//     array[n - 1] = 0;
//     for (int i = n - 2; i >= 0; i--)
//     {
//         array[i] = coord % (n - i);
//         for (int j = i + 1; j < n; j++)
//         {
//             if (array[j] >= array[i])
//             {
//                 array[j] += 1;
//             }
//         }
//         coord = coord / (n - i);
//     }
// }

// void set_partial_permutation_coord(int *array, int coord, int n, int axis_offset)
// {
//     for (int i = 0; i < n; i++)
//     {
//         array[i] = -1;
//     }

//     int comb_coord = coord / NUM_AXIS_ELEMS_PERM;
//     int perm_coord = coord % NUM_AXIS_ELEMS_PERM;

//     int perm[NUM_AXIS_ELEMS];
//     set_permutation_coord(perm, perm_coord, NUM_AXIS_ELEMS);

//     int i = NUM_AXIS_ELEMS - 1;
//     for (int c = n - 1; c >= 0; c--)
//     {
//         if (comb_coord >= COMBINATION[i][c])
//         {
//             comb_coord -= COMBINATION[i][c];
//             array[c] = perm[i--] + axis_offset;
//             if (i < 0)
//             {
//                 break;
//             }
//         }
//     }
// }
