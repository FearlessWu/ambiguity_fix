#include "ambg_fix.h"

void main(void)
{
    fp64 q[36] = {
         0.6098,  0.8399, -0.8271, -0.4150,  0.5950,  0.2838,
         0.8399,  1.1692, -1.1241, -0.5666,  0.8282,  0.3987,
        -0.8271, -1.1241,  2.0589,  0.9191, -0.6979, -0.1592,
        -0.4150, -0.5666,  0.9191,  0.4232, -0.3634, -0.1075,
         0.5950,  0.8282, -0.6979, -0.3634,  0.6035,  0.3061,
         0.2838,  0.3987, -0.1592, -0.1075,  0.3061,  0.1972
    };

    matrix_t A;
    matrix_t L;
    matrix_t D;
    matrix_init(&A, 6, 6);
    matrix_eye(&L, 6);
    matrix_init(&D, 6, 6);

    
    for (uint8_t i = 0; i < 6; ++i)
    {
        for (uint8_t j = 0; j < 6; ++j)
        {
            A.element[i][j] = q[i * 6 + j];
        }
    }
    matrix_print(A);
    //A.element[0][0] = 5;
    //A.element[0][1] = -4;
    //A.element[0][2] = 1;
    //A.element[1][0] = -4;
    //A.element[1][1] = 6;
    //A.element[1][2] = -4;
    //A.element[2][0] = 1;
    //A.element[2][1] = -4;
    //A.element[2][2] = -6;
    LDL_decomposition(&A, &L, &D);

    // matrix_print(D);
    // matrix_print(L);

    matrix_t Z;
    matrix_eye(&Z, 6);
    if (gauss_transform(&L, &Z))
    {
        matrix_print(Z);
        matrix_print(L);
    }
    matrix_t ZT;
    matrix_t tmp;
    matrix_t tmp1;
    matrix_init(&ZT, 6, 6);
    matrix_init(&tmp, 6, 6);
    matrix_init(&tmp1, 6, 6);
    matrix_trs(&Z, &ZT);
    matrix_mlt(&ZT, &A, &tmp);
    matrix_mlt(&tmp, &Z, &tmp1);
    matrix_print(tmp1);
    
    system("pause");
}