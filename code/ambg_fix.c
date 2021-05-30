#include "ambg_fix.h"

int32_t LDL_decomposition(matrix_t *A_in, matrix_t *L, matrix_t *D)
{
    uint32_t dem = A_in->col;
    if (dem <= 0)
    {
        return -1;
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
                return -1;
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
    return 0;
}