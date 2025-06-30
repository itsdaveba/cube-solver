#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <numpy/ndarrayobject.h>

#include "defs.h"
#include "utils.h"


int get_orientation_coord(int *array, int v, int n)
{
    int coord = 0;
    for (int i = 0; i < n - 1; i++)
    {
        coord *= v;
        coord += array[i];
    }
    return coord;
}

int get_permutation_coord(int *array, int n)
{
    int coord = 0;
    for (int i = 0; i < n - 1; i++)
    {
        coord *= n - i;
        for (int j = i + 1; j < n; j++)
            if (array[i] > array[j])
                coord++;
    }
    return coord;
}

int get_combination_coord(int *array, int n)
{
    int coord = 0;
    for (int i = 0; i < n; i++)
        coord += COMBINATION[array[i]][i + 1];
    return coord;
}

int get_partial_permutation_coord(int *array, const int *orbits, int orbit, int n)
{
    int i = 0;
    int perm[NUM_ORBIT_ELEMS];
    int comb[NUM_ORBIT_ELEMS];
    for (int j = 0; j < n && i < NUM_ORBIT_ELEMS; j++)
        if (array[j] != NONE && orbits[array[j]] == orbit)
        {
            perm[i] = array[j];
            comb[i++] = j;
        }
    int perm_coord = get_permutation_coord(perm, NUM_ORBIT_ELEMS);
    int comb_coord = get_combination_coord(comb, NUM_ORBIT_ELEMS);
    return perm_coord + comb_coord * OP_SIZE;
}

void set_orientation_coord(int *array, int coord, int v, int n)
{
    int o = 0;
    for (int i = n - 2; i >= 0; i--)
    {
        array[i] = coord % v;
        o -= array[i];
        coord /= v;
    }
    o %= v;
    array[n - 1] = o < 0 ? o + v : o;
}

void set_permutation_coord(int *array, int coord, int n)
{
    array[n - 1] = 0;
    for (int i = n - 2; i >= 0; i--)
    {
        array[i] = coord % (n - i);
        coord /= n - i;
        for (int j = i + 1; j < n; j++)
        {
            if (array[j] >= array[i])
                array[j]++;
        }
    }
}

void set_combination_coord(int *array, int coord, int n)
{
    int i = n - 1;
    int m = 1;
    while (coord >= COMBINATION[m][n])
        m++;
    for (int c = m - 1; c >= 0 && i >= 0; c--)
        if (coord >= COMBINATION[c][i + 1])
        {
            coord -= COMBINATION[c][i + 1];
            array[i--] = c;
        }
}

void set_partial_permutation_coord(int *array, int coord, int n)
{
    int perm[NUM_ORBIT_ELEMS];
    int comb[NUM_ORBIT_ELEMS];
    int comb_coord = coord / FACTORIAL[NUM_ORBIT_ELEMS];
    int perm_coord = coord % FACTORIAL[NUM_ORBIT_ELEMS];
    set_permutation_coord(perm, perm_coord, NUM_ORBIT_ELEMS);
    set_combination_coord(comb, comb_coord, NUM_ORBIT_ELEMS);
    for (int i = 0; i < n; i++)
        array[i] = NONE;
    for (int i = 0; i < NUM_ORBIT_ELEMS; i++)
        array[comb[i]] = perm[i];
}

npy_intp get_table_index(int *phase_coords, int *indexes, npy_intp *dims, int nd)
{
    npy_intp index = 0;
    for (int i = 0; i < nd; i++)
    {
        index *= dims[i];
        index += indexes[i] == NONE ? phase_coords[i] : phase_coords[indexes[i]];
    }
    return index;
}

void set_phase_coords(int *phase_coords, npy_intp index, int *indexes, npy_intp *dims, int nd)
{
    int idx;
    for (int i = nd - 1; i >= 0; i--)
    {
        idx = indexes[i] == NONE ? i : indexes[i];
        phase_coords[idx] = index % dims[i];
        index /= dims[i];
    }
}
