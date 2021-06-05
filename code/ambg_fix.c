#include "ambg_fix.h"

int8_t LDL_decomposition(matrix_t *A_in, matrix_t *L, matrix_t *D)
{
    uint32_t dem = A_in->col;
    if (dem <= 0)
    {
        return 0;
    }

    matrix_t B;
    matrix_init(&B, dem, dem);
    matrix_copy(A_in, &B);
    matrix_t *A = &B;

    // symmetric check
    for (uint32_t i = 0; i < dem; ++i)
    {
        for (uint32_t j = 0; j < dem; ++j)
        {
            if (A->element[i][j] != A->element[j][i])
            {
                matrix_free(A);
                return 0;
            }
        }
    }

    for (uint32_t k = 0; k < dem; ++k)
    {
        if (k == 0)
        {
            D->element[k][k] = A->element[k][k];

            for (uint32_t m = 0; m < dem; ++m)
            {
                L->element[m][k] = A->element[m][k] / D->element[k][k];
            }
        }
        else
        {
            fp64 sum1 = 0;
            for (uint32_t m = 0; m < k; ++m)
            {
                sum1 += A->element[k][m] * L->element[k][m];
            }
            D->element[k][k] = A->element[k][k] - sum1;

            for (uint32_t j = k + 1; j < dem; ++j)
            {
                fp64 sum2 = 0;
                for (uint32_t m = 0; m < k; ++m)
                {
                    sum2 += A->element[j][m] * L->element[k][m];
                }
                A->element[j][k] -= sum2;
                L->element[j][k] = A->element[j][k] / D->element[k][k];
            }
        }
    }
    return 1;
}
static int8_t lower_matrix_check(const matrix_t *L)
{
    if (!L->is_valid || L->row <= 0 || L->col <= 0 || L->col != L->row)
    {
        return -1;
    }
    uint32_t i = 0;
    uint32_t j = 0;
    for (; i < L->row; ++i)
    {
        for (; j < L->col; ++j)
        {
            if (i == j)
            {
                fp64 tmp = L->element[i][j] - 1;
                if (!IS_ZEROS(tmp))
                {
                    return 0;
                }
            }
            else if (j > i)
            {
                if (!IS_ZEROS(L->element[i][j]))
                {
                    return 0;
                }
            }
        }
    }
}
static fp64 get_integer(fp64 a)
{
    fp64 tmp = floor(a);
    fp64 frac = a - tmp;
    
    if (frac >= 0.5)
    {
        tmp += 1;
    }

    return tmp;
}
int8_t gauss_transform(matrix_t *L, matrix_t *Z)
{
    if (!lower_matrix_check(L))
    {
        return -1;
    }

    if (L->col <= 2)
    {
        return -1;
    }
    uint32_t col = 0;
    uint32_t row = 0;
    for (row = 1; row < L->row; ++row)
    {
        for (col = 0; col < row; ++col)
        {
            fp64 u = get_integer(L->element[row][col]);
            if (!IS_ZEROS(u))
            {
                Z->element[row][col] = -u;
                uint32_t i = 0;
                for (i = row; i < L->row; ++i)
                {
                    L->element[i][col] -= u * L->element[i][row];
                }
            }
        }
    }

    return 1;

}