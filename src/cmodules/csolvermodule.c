#define PY_SSIZE_T_CLEAN
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/ndarrayobject.h>

#include "defs.h"
#include "utils.h"

#include "korf.h"
#include "kociemba.h"
#include "thistlethwaite.h"

static int corner_orientation[NUM_CORNERS];
static int edge_orientation[NUM_EDGES];
static int corner_permutation[NUM_CORNERS];
static int edge_permutation[NUM_EDGES];

static void reset()
{
    set_orientation_coord(corner_orientation, 0, 3, NUM_CORNERS);
    set_orientation_coord(edge_orientation, 0, 2, NUM_EDGES);
    set_permutation_coord(corner_permutation, 0, NUM_CORNERS);
    set_permutation_coord(edge_permutation, 0, NUM_EDGES);
}

static void apply_move(int *co, int *eo, int *cp, int *ep, int move)
{
    int face = move / 3;
    int shift = move % 3 + 1;
    for (int i = 0; i < 4; i++)
    {
        co[PERM[face][0][i]] = (corner_orientation[PERM[face][0][(i + shift) % 4]] + (shift == 2 ? 0 : ORIENT_DIFF[face][0][i])) % 3;
        eo[PERM[face][1][i]] = (edge_orientation[PERM[face][1][(i + shift) % 4]] + (shift == 2 ? 0 : ORIENT_DIFF[face][1][i])) % 2;
        cp[PERM[face][0][i]] = corner_permutation[PERM[face][0][(i + shift) % 4]];
        ep[PERM[face][1][i]] = edge_permutation[PERM[face][1][(i + shift) % 4]];
    }
}

static void get_coords(int *coords, int partial_corner_perm, int partial_edge_perm)
{
    coords[0] = get_orientation_coord(corner_orientation, 3, NUM_CORNERS);
    coords[1] = get_orientation_coord(edge_orientation, 2, NUM_EDGES);
    int i = 2;
    if (!partial_corner_perm)
        coords[i++] = get_permutation_coord(corner_permutation, NUM_CORNERS);
    else
        for (int orbit = 0; orbit < 2; orbit++)
            coords[i++] = get_partial_permutation_coord(corner_permutation, CORNER_ORBITS, orbit, NUM_CORNERS);
    if (!partial_edge_perm)
        coords[i++] = get_permutation_coord(edge_permutation, NUM_EDGES);
    else
        for (int orbit = 0; orbit < 3; orbit++)
            coords[i++] = get_partial_permutation_coord(edge_permutation, EDGE_ORBITS, orbit, NUM_EDGES);
}

static void next_position(int *next_coords, int *coords, int move, npy_uint16 *transition_tables[], int partial_corner_perm, int partial_edge_perm)
{
    next_coords[0] = transition_tables[0][NUM_MOVES * coords[0] + move];
    next_coords[1] = transition_tables[1][NUM_MOVES * coords[1] + move];
    int i = 2;
    int j = 2;
    if (!partial_corner_perm)
        next_coords[i++] = transition_tables[2][NUM_MOVES * coords[j++] + move];
    else
        for (int orbit = 0; orbit < 2; orbit++)
            next_coords[i++] = transition_tables[2][NUM_MOVES * coords[j++] + move];
    if (!partial_edge_perm)
        next_coords[i++] = transition_tables[3][NUM_MOVES * coords[j++] + move];
    else
        for (int orbit = 0; orbit < 3; orbit++)
            next_coords[i++] = transition_tables[3][NUM_MOVES * coords[j++] + move];
}

static PyObject *generate_transition_table(PyObject *self, PyObject *args, PyObject *kwargs)
{
    // parse args
    int coord_size;
    const char *coord_name;
    char *keywords[] = {"coord_name", "coord_size", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "si", keywords, &coord_name, &coord_size))
        return NULL;
    if (coord_size - 1 > NPY_MAX_UINT16)
        return NULL;

    // create array
    const npy_intp dims[] = {coord_size, NUM_MOVES};
    PyObject *table = PyArray_SimpleNew(2, dims, NPY_UINT16);
    if (!PyArray_IS_C_CONTIGUOUS((PyArrayObject *)table))
        return NULL;
    npy_uint16 *data = PyArray_DATA((PyArrayObject *)table);

    int co[NUM_CORNERS];
    int eo[NUM_EDGES];
    int cp[NUM_CORNERS];
    int ep[NUM_EDGES];

    // fill table
    if (!strcmp(coord_name, "co"))
        for (npy_uint16 coord = 0; coord < coord_size; coord++)
        {
            set_orientation_coord(corner_orientation, coord, 3, NUM_CORNERS);
            for (int move = 0; move < NUM_MOVES; move++)
            {
                memcpy(co, corner_orientation, NUM_CORNERS * sizeof(int));
                apply_move(co, eo, cp, ep, move);
                *data++ = (npy_uint16)get_orientation_coord(co, 3, NUM_CORNERS);
            }
        }
    else if (!strcmp(coord_name, "eo"))
        for (npy_uint16 coord = 0; coord < coord_size; coord++)
        {
            set_orientation_coord(edge_orientation, coord, 2, NUM_EDGES);
            for (int move = 0; move < NUM_MOVES; move++)
            {
                memcpy(eo, edge_orientation, NUM_EDGES * sizeof(int));
                apply_move(co, eo, cp, ep, move);
                *data++ = (npy_uint16)get_orientation_coord(eo, 2, NUM_EDGES);
            }
        }
    else if (!strcmp(coord_name, "cp"))
        for (npy_uint16 coord = 0; coord < coord_size; coord++)
        {
            set_permutation_coord(corner_permutation, coord, NUM_CORNERS);
            for (int move = 0; move < NUM_MOVES; move++)
            {
                memcpy(cp, corner_permutation, NUM_CORNERS * sizeof(int));
                apply_move(co, eo, cp, ep, move);
                *data++ = (npy_uint16)get_permutation_coord(cp, NUM_CORNERS);
            }
        }
    else if (!strcmp(coord_name, "ep"))
        for (npy_uint16 coord = 0; coord < coord_size; coord++)
        {
            set_permutation_coord(edge_permutation, coord, NUM_EDGES);
            for (int move = 0; move < NUM_MOVES; move++)
            {
                memcpy(ep, edge_permutation, NUM_EDGES * sizeof(int));
                apply_move(co, eo, cp, ep, move);
                *data++ = (npy_uint16)get_permutation_coord(ep, NUM_EDGES);
            }
        }
    else if (!strcmp(coord_name, "pcp"))
        for (npy_uint16 coord = 0; coord < coord_size; coord++)
        {
            set_partial_permutation_coord(corner_permutation, coord, NUM_CORNERS);
            for (int move = 0; move < NUM_MOVES; move++)
            {
                memcpy(cp, corner_permutation, NUM_CORNERS * sizeof(int));
                apply_move(co, eo, cp, ep, move);
                *data++ = (npy_uint16)get_partial_permutation_coord(cp, CORNER_ORBITS, 0, NUM_CORNERS);
            }
        }
    else if (!strcmp(coord_name, "pep"))
        for (npy_uint16 coord = 0; coord < coord_size; coord++)
        {
            set_partial_permutation_coord(edge_permutation, coord, NUM_EDGES);
            for (int move = 0; move < NUM_MOVES; move++)
            {
                memcpy(ep, edge_permutation, NUM_EDGES * sizeof(int));
                apply_move(co, eo, cp, ep, move);
                *data++ = (npy_uint16)get_partial_permutation_coord(ep, EDGE_ORBITS, 0, NUM_EDGES);
            }
        }
    else
        return NULL;

    return table;
}

static PyObject *generate_pruning_table(PyObject *self, PyObject *args, PyObject *kwargs)
{
    // parse args
    int phase;
    const char *name = "";
    PyObject *solver;
    PyObject *shape;
    PyObject *idxs;
    char *keywords[] = {"solver", "phase", "shape", "indexes", "name", NULL};
    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "OiOO|s", keywords, &solver, &phase, &shape, &idxs, &name))
        return NULL;
    // parse dims
    PyObject *tmp;
    if (!PyTuple_Check(shape))
    {
        tmp = shape;
        shape = PyTuple_Pack(1, shape);
        Py_DECREF(tmp);
    }
    int nd = PyTuple_Size(shape);
    npy_intp dims[nd];
    npy_intp table_size = 1;
    for (int i = 0; i < nd; i++)
    {
        tmp = PyTuple_GetItem(shape, i);
        dims[i] = PyLong_AsLong(tmp);
        table_size *= dims[i];
    }
    // parse indexes
    if (!PyTuple_Check(idxs))
    {
        tmp = idxs;
        idxs = PyTuple_Pack(1, idxs);
        Py_DECREF(tmp);
    }
    int ni = PyTuple_Size(idxs);
    int indexes[nd];
    for (int i = 0; i < nd; i++)
    {
        if (i < ni)
        {
            tmp = PyTuple_GetItem(idxs, i);
            indexes[i] = PyLong_Check(tmp) ? PyLong_AsLong(tmp) : NONE;
        }
        else
            indexes[i] = NONE;
    }

    // solver functions
    int (*get_phase_coords_size)(int);
    void (*get_phase_coords)(int *, int *, int);
    void (*set_coords)(int *, int *, int);
    name = solver->ob_type->tp_name;
    if (!strcmp(name, "Korf"))
    {
        get_phase_coords_size = &korf_get_phase_coords_size;
        get_phase_coords = &korf_get_phase_coords;
        set_coords = &korf_set_coords;
    }
    else if (!strcmp(name, "Thistlethwaite"))
    {
        get_phase_coords_size = &thistle_get_phase_coords_size;
        get_phase_coords = &thistle_get_phase_coords;
        set_coords = &thistle_set_coords;
    }
    else if (!strcmp(name, "Kociemba"))
    {
        get_phase_coords_size = &kociemba_get_phase_coords_size;
        get_phase_coords = &kociemba_get_phase_coords;
        set_coords = &kociemba_set_coords;
    }
    else
        return NULL;

    // class attributes
    PyObject *pcp = PyObject_GetAttrString(solver, "partial_corner_perm");
    PyObject *pep = PyObject_GetAttrString(solver, "partial_edge_perm");
    int partial_corner_perm = PyLong_AsLong(pcp);
    int partial_edge_perm = PyLong_AsLong(pep);
    PyObject *phase_moves = PyObject_GetAttrString(solver, "phase_moves");
    phase_moves = PyList_GetItem(phase_moves, phase);
    int n_moves = PyList_Size(phase_moves);
    int moves[n_moves];
    for (int i = 0; i < n_moves; i++)
    {
        tmp = PyList_GetItem(phase_moves, i);
        moves[i] = PyLong_AsLong(tmp);
    }

    // transition tables
    PyObject *trans_table;
    npy_uint16 *transition_tables[4];
    PyObject *trans_tables = PyObject_GetAttrString(solver, "transition_tables");
    trans_table = PyDict_GetItemString(trans_tables, "co");
    transition_tables[0] = PyArray_DATA((PyArrayObject *)trans_table);
    trans_table = PyDict_GetItemString(trans_tables, "eo");
    transition_tables[1] = PyArray_DATA((PyArrayObject *)trans_table);
    trans_table = PyDict_GetItemString(trans_tables, partial_corner_perm ? "pcp" : "cp");
    transition_tables[2] = PyArray_DATA((PyArrayObject *)trans_table);
    trans_table = PyDict_GetItemString(trans_tables, partial_edge_perm ? "pep" : "ep");
    transition_tables[3] = PyArray_DATA((PyArrayObject *)trans_table);

    // create array
    PyObject *table = PyArray_SimpleNew(nd, dims, NPY_INT8);
    if (!PyArray_IS_C_CONTIGUOUS((PyArrayObject *)table))
        return NULL;
    PyArray_FILLWBYTE((PyArrayObject *)table, NONE);
    npy_int8 *data = PyArray_DATA((PyArrayObject *)table);

    // fill table
    reset();
    npy_intp counter;
    npy_int8 depth = 0;
    int coords_size = 4 + partial_corner_perm + 2 * partial_edge_perm;
    int coords[coords_size];
    int next_coords[coords_size];
    get_coords(coords, partial_corner_perm, partial_edge_perm);
    int phase_coords_size = get_phase_coords_size(phase);
    if (phase_coords_size == NONE)
        return NULL;
    int phase_coords[phase_coords_size];
    get_phase_coords(phase_coords, coords, phase);
    npy_intp table_index = get_table_index(phase_coords, indexes, dims, nd);
    data[table_index] = depth;
    do
    {
        counter = 0;
        for (npy_intp i = 0; i < table_size; i++)
            if (data[i] == depth)
            {
                set_phase_coords(phase_coords, i, indexes, dims, nd);
                set_coords(coords, phase_coords, phase);
                for (int m = 0; m < n_moves; m++)
                {
                    next_position(next_coords, coords, moves[m], transition_tables, partial_corner_perm, partial_edge_perm);
                    get_phase_coords(phase_coords, next_coords, phase);
                    table_index = get_table_index(phase_coords, indexes, dims, nd);
                    if (data[table_index] == NONE)
                    {
                        data[table_index] = depth + 1;
                        counter++;
                    }
                }
            }
        depth++;
    } while (counter);
    return table;
}

static struct PyMethodDef csolver_methods[] = {
    {"generate_transition_table", (PyCFunction)(void (*)(void))generate_transition_table, METH_VARARGS | METH_KEYWORDS,
     "Generate the cube coordinate transition table."},
    {"generate_pruning_table", (PyCFunction)(void (*)(void))generate_pruning_table, METH_VARARGS | METH_KEYWORDS,
     "Generate the phase coordinates pruning table."},
    {NULL, NULL, 0, NULL} /* Sentinel */
};

static struct PyModuleDef csolver_module = {
    .m_base = PyModuleDef_HEAD_INIT,
    .m_name = "csolver",
    .m_size = 0, // non-negative
    .m_methods = csolver_methods};

PyMODINIT_FUNC PyInit_csolver(void)
{
    if (PyArray_ImportNumPyAPI() < 0)
        return NULL;
    return PyModuleDef_Init(&csolver_module);
}
