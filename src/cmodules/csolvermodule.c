// #define PY_SSIZE_T_CLEAN
// #define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
// #include <Python.h>
// #include <numpy/ndarrayobject.h>
// #include "queue.h"
// #include "cutils.h"
// #include "ckociemba.h"

// static int corner_orientation[NUM_CORNERS];
// static int edge_orientation[NUM_EDGES];
// static int corner_permutation[NUM_CORNERS];
// static int edge_permutation[NUM_EDGES];

// static void reset_position()
// {
//     set_orientation_coord(corner_orientation, 0, 3, NUM_CORNERS);
//     set_orientation_coord(edge_orientation, 0, 2, NUM_EDGES);
//     set_permutation_coord(corner_permutation, 0, NUM_CORNERS);
//     set_permutation_coord(edge_permutation, 0, NUM_EDGES);
// }

// static void apply_move(int *co, int *eo, int *cp, int *ep, int move)
// {
//     int move_div = move / 3;
//     int move_mod = move % 3;
//     for (int i = 0; i < 4; i++)
//     {
//         co[MOVES[move_div][0][i]] = (corner_orientation[MOVES[move_div][0][(i + 1 + move_mod) % 4]] + (move_mod == 1 ? 0 : ORIENTATION[move_div][0][i])) % 3;
//         eo[MOVES[move_div][1][i]] = (edge_orientation[MOVES[move_div][1][(i + 1 + move_mod) % 4]] + (move_mod == 1 ? 0 : ORIENTATION[move_div][1][i])) % 2;
//         cp[MOVES[move_div][0][i]] = corner_permutation[MOVES[move_div][0][(i + 1 + move_mod) % 4]];
//         ep[MOVES[move_div][1][i]] = edge_permutation[MOVES[move_div][1][(i + 1 + move_mod) % 4]];
//     }
// }

// static void get_coords(int *coords, int partial_corner_perm, int partial_edge_perm)
// {
//     int i = 0;
//     coords[i++] = get_orientation_coord(corner_orientation, 3, NUM_CORNERS);
//     coords[i++] = get_orientation_coord(edge_orientation, 2, NUM_EDGES);
//     if (partial_corner_perm)
//     {
//         for (int axis = 0; axis < 2; axis++)
//         {
//             coords[i++] = get_partial_permutation_coord(corner_permutation, NUM_CORNERS, CORNER_AXES, axis);
//         }
//     }
//     else
//     {
//         coords[i++] = get_permutation_coord(corner_permutation, NUM_CORNERS);
//     }
//     if (partial_edge_perm)
//     {
//         for (int axis = 0; axis < 3; axis++)
//         {
//             coords[i++] = get_partial_permutation_coord(edge_permutation, NUM_EDGES, EDGE_AXES, axis);
//         }
//     }
//     else
//     {
//         coords[i++] = get_permutation_coord(edge_permutation, NUM_EDGES);
//     }
// }

// static void next_position(int *next_coords, int partial_corner_perm, int partial_edge_perm, npy_uint16 *transition_arrays[], int *coords, int move)
// {
//     int i = 0;
//     next_coords[i] = transition_arrays[0][coords[i] * NUM_MOVES + move];
//     i++;
//     next_coords[i] = transition_arrays[1][coords[i] * NUM_MOVES + move];
//     i++;
//     next_coords[i] = transition_arrays[2][coords[i] * NUM_MOVES + move];
//     i++;
//     if (partial_corner_perm)
//     {
//         next_coords[i] = transition_arrays[2][coords[i] * NUM_MOVES + move];
//         i++;
//     }
//     next_coords[i] = transition_arrays[3][coords[i] * NUM_MOVES + move];
//     i++;
//     if (partial_edge_perm)
//     {
//         next_coords[i] = transition_arrays[3][coords[i] * NUM_MOVES + move];
//         i++;
//         next_coords[i] = transition_arrays[3][coords[i] * NUM_MOVES + move];
//         i++;
//     }
// }

// static PyObject *generate_transition_table(PyObject *self, PyObject *args)
// {
//     int co[NUM_CORNERS];
//     int eo[NUM_EDGES];
//     int cp[NUM_CORNERS];
//     int ep[NUM_EDGES];

//     int size;
//     char *coord_type;
//     if (!PyArg_ParseTuple(args, "si", &coord_type, &size))
//     {
//         return NULL;
//     }

//     npy_intp const shape[2] = {size, NUM_MOVES};
//     PyObject *transition_table = PyArray_SimpleNew(2, shape, NPY_UINT16);
//     npy_uint16 *transition_array = PyArray_DATA((PyArrayObject *)transition_table);

//     if (!strcmp(coord_type, "co"))
//     {
//         for (int coord = 0; coord < size; coord++)
//         {
//             set_orientation_coord(corner_orientation, coord, 3, NUM_CORNERS);
//             for (int move = 0; move < NUM_MOVES; move++)
//             {
//                 memcpy(co, corner_orientation, NUM_CORNERS * sizeof(int));
//                 apply_move(co, eo, cp, ep, move);
//                 transition_array[coord * NUM_MOVES + move] = (npy_uint16)get_orientation_coord(co, 3, NUM_CORNERS);
//             }
//         }
//     }
//     else if (!strcmp(coord_type, "eo"))
//     {
//         for (int coord = 0; coord < size; coord++)
//         {
//             set_orientation_coord(edge_orientation, coord, 2, NUM_EDGES);
//             for (int move = 0; move < NUM_MOVES; move++)
//             {
//                 memcpy(eo, edge_orientation, NUM_EDGES * sizeof(int));
//                 apply_move(co, eo, cp, ep, move);
//                 transition_array[coord * NUM_MOVES + move] = (npy_uint16)get_orientation_coord(eo, 2, NUM_EDGES);
//             }
//         }
//     }
//     else if (!strcmp(coord_type, "cp"))
//     {
//         for (int coord = 0; coord < size; coord++)
//         {
//             set_permutation_coord(corner_permutation, coord, NUM_CORNERS);
//             for (int move = 0; move < NUM_MOVES; move++)
//             {
//                 memcpy(cp, corner_permutation, NUM_CORNERS * sizeof(int));
//                 apply_move(co, eo, cp, ep, move);
//                 transition_array[coord * NUM_MOVES + move] = (npy_uint16)get_permutation_coord(cp, NUM_CORNERS);
//             }
//         }
//     }
//     else if (!strcmp(coord_type, "ep"))
//     {
//         for (int coord = 0; coord < size; coord++)
//         {
//             set_permutation_coord(edge_permutation, coord, NUM_EDGES);
//             for (int move = 0; move < NUM_MOVES; move++)
//             {
//                 memcpy(ep, edge_permutation, NUM_EDGES * sizeof(int));
//                 apply_move(co, eo, cp, ep, move);
//                 transition_array[coord * NUM_MOVES + move] = (npy_uint16)get_permutation_coord(ep, NUM_EDGES);
//             }
//         }
//     }
//     else if (!strcmp(coord_type, "pcp"))
//     {
//         for (int coord = 0; coord < size; coord++)
//         {
//             set_partial_permutation_coord(corner_permutation, coord, NUM_CORNERS, CORNER_AXIS_OFFSET[0]);
//             for (int move = 0; move < NUM_MOVES; move++)
//             {
//                 memcpy(cp, corner_permutation, NUM_CORNERS * sizeof(int));
//                 apply_move(co, eo, cp, ep, move);
//                 transition_array[coord * NUM_MOVES + move] = (npy_uint16)get_partial_permutation_coord(cp, NUM_CORNERS, CORNER_AXES, 0);
//             }
//         }
//     }
//     else if (!strcmp(coord_type, "pep"))
//     {
//         for (int coord = 0; coord < size; coord++)
//         {
//             set_partial_permutation_coord(edge_permutation, coord, NUM_EDGES, EDGE_AXIS_OFFSET[0]);
//             for (int move = 0; move < NUM_MOVES; move++)
//             {
//                 memcpy(ep, edge_permutation, NUM_EDGES * sizeof(int));
//                 apply_move(co, eo, cp, ep, move);
//                 transition_array[coord * NUM_MOVES + move] = (npy_uint16)get_partial_permutation_coord(ep, NUM_EDGES, EDGE_AXES, 0);
//             }
//         }
//     }
//     else
//     {
//         return NULL;
//     }

//     return transition_table;
// }

// static PyObject *generate_pruning_table(PyObject *self, PyObject *args)
// {
//     int phase;
//     PyObject *shape_tuple;
//     PyObject *index_list;
//     if (!PyArg_ParseTuple(args, "OiOO", &self, &phase, &shape_tuple, &index_list))
//     {
//         return NULL;
//     }

//     PyObject *partial_corner_perm_attr = PyObject_GetAttrString(self, "partial_corner_perm");
//     PyObject *partial_edge_perm_attr = PyObject_GetAttrString(self, "partial_edge_perm");
//     int partial_corner_perm = PyObject_IsTrue(partial_corner_perm_attr);
//     int partial_edge_perm = PyObject_IsTrue(partial_edge_perm_attr);

//     reset_position();

//     int coords_size = 4 + partial_corner_perm + 2 * partial_edge_perm;
//     int coords[coords_size];
//     get_coords(coords, partial_corner_perm, partial_edge_perm);
//     int phase_coords_len = get_phase_coords_len(phase);
//     int phase_coords[phase_coords_len];
//     get_phase_coords(phase_coords, phase, coords);
//     Py_ssize_t index_len = PyList_Size(index_list);
//     long indexes[index_len];
//     int prune_coords[index_len];
//     for (Py_ssize_t i = 0; i < index_len; i++)
//     {
//         PyObject *item = PyList_GetItem(index_list, i);
//         indexes[i] = PyLong_AsLong(item);
//         prune_coords[i] = phase_coords[indexes[i]];
//     }

//     long size = 1;
//     Py_ssize_t num_dim = PyTuple_Size(shape_tuple);
//     long shape[num_dim];
//     for (Py_ssize_t i = 0; i < num_dim; i++)
//     {
//         PyObject *item = PyTuple_GetItem(shape_tuple, i);
//         shape[i] = PyLong_AsLong(item);
//         size *= shape[i];
//     }
//     PyObject *pruning_table = PyArray_SimpleNew(num_dim, shape, NPY_INT8);
//     npy_int8 *pruning_array = PyArray_DATA((PyArrayObject *)pruning_table);

//     for (long i = 0; i < size; i++)
//     {
//         pruning_array[i] = EMPTY;
//     }

//     long index = 0;
//     for (int i = 0; i < num_dim; i++)
//     {
//         index *= shape[i];
//         index += prune_coords[i];
//     }
//     pruning_array[index] = 0;

//     PyObject *next_moves = PyObject_GetAttrString(self, "next_moves");
//     PyObject *phase_next_moves = PyList_GetItem(next_moves, phase);
//     PyObject *keys = PyDict_Keys(phase_next_moves);
//     Py_ssize_t keys_len = PyList_Size(keys);
//     PyObject *key = PyList_GetItem(keys, keys_len - 1);  // TODO rename
//     PyObject *phase_moves = PyDict_GetItem(phase_next_moves, key);
//     Py_ssize_t phase_moves_len = PyList_Size(phase_moves);
//     int moves[phase_moves_len];
//     for (Py_ssize_t i = 0; i < phase_moves_len; i++)
//     {
//         PyObject *item = PyList_GetItem(phase_moves, i);
//         moves[i] = PyLong_AsLong(item);
//     }

//     PyObject *transition_tables = PyObject_GetAttrString(self, "transition_tables");  // TODO check if use_transition_table?
//     Py_ssize_t transition_tables_len = PyList_Size(transition_tables);
//     npy_uint16 *transition_arrays[transition_tables_len];
//     for (Py_ssize_t i = 0; i < transition_tables_len; i++)
//     {
//         PyObject *table = PyList_GetItem(transition_tables, i);
//         transition_arrays[i] = PyArray_DATA((PyArrayObject *)table);
//     }

//     int depth = 0;
//     int i, m;
//     long idx;
//     int next_coords[coords_size];
//     do
//     {
//         printf("depth: %d ", depth);
//         clock_t start = clock();
//         idx = -1;
//         for (index = 0; index < size; index++)
//         {
//             if (pruning_array[index] == depth)
//             {
//                 idx = index;
//                 for (i = num_dim - 1; i >= 0; i--)
//                 {
//                     prune_coords[i] = idx % shape[i];
//                     idx /= shape[i];
//                 }
//                 for (i = 0; i < index_len; i++)
//                 {
//                     phase_coords[indexes[i]] = prune_coords[i];
//                 }
//                 get_coords_from_phase_coords(coords, phase, phase_coords);  // TODO transition table per phase (test fisrt thist)
//                 for (m = 0; m < phase_moves_len; m++)
//                 {
//                     next_position(next_coords, partial_corner_perm, partial_edge_perm, transition_arrays, coords, moves[m]);
//                     get_phase_coords(phase_coords, phase, next_coords);
//                     for (i = 0; i < index_len; i++)
//                     {
//                         prune_coords[i] = phase_coords[indexes[i]];
//                     }
//                     idx = 0;
//                     for (i = 0; i < num_dim; i++)
//                     {
//                         idx *= shape[i];
//                         idx += prune_coords[i];
//                     }
//                     if (pruning_array[idx] == EMPTY)
//                     {
//                         pruning_array[idx] = depth + 1;
//                     }
//                 }
//             }
//         }
//         depth++;
//         printf("time: %ld seconds\n", (clock() - start) / CLOCKS_PER_SEC);
//     } while (idx >= 0);

//     return pruning_table;
// }

// static PyMethodDef SolverMethods[] = {
//     {"generate_transition_table", generate_transition_table, METH_VARARGS, NULL}, // TODO doc
//     {"generate_pruning_table", generate_pruning_table, METH_VARARGS, NULL},
//     {NULL, NULL, 0, NULL}};

// static struct PyModuleDef csolvermodule = {
//     PyModuleDef_HEAD_INIT, "csolver", NULL, -1, SolverMethods, NULL, NULL, NULL, NULL};

// PyObject *PyInit_csolver()
// {
//     if (PyArray_ImportNumPyAPI() < 0)
//     {
//         return NULL;
//     }
//     return PyModule_Create(&csolvermodule);
// }