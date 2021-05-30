#include "ambg_fix.h"

void main(void)
{
    matrix_t A;
    matrix_t L;
    matrix_t D;
    matrix_init(&A, 3, 3);
    matrix_eye(&L, 3);
    matrix_init(&D, 3, 3);
    A.element[0][0] = 5;
    A.element[0][1] = -4;
    A.element[0][2] = 1;
    A.element[1][0] = -4;
    A.element[1][1] = 6;
    A.element[1][2] = -4;
    A.element[2][0] = 1;
    A.element[2][1] = -4;
    A.element[2][2] = -6;
    LDL_decomposition(&A, &L, &D);
    matrix_print(A);
    matrix_print(L);
    matrix_print(D);

    matrix_t LT;
    matrix_t tmp1;
    matrix_init(&tmp1, L.col, L.row);
    matrix_init(&LT, L.row, L.col);
    matrix_trs(&L, &LT);
    matrix_mlt(&L, &D, &tmp1);
    matrix_mlt(&tmp1, &LT, &L);
    matrix_print(L);
    system("pause");
}