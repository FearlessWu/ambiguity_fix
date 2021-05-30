#ifndef _AMBG_FIX_H_
#define _AMBG_FIX_H_
#include "Matrix.h"
/**
 *@brief            LDL decompodition
 *@param[in]        A: symmetric matrix
 *@param[out]       L: lower triangular matrix
 *@param[out]       D: dialog matrix
 *@retval           0: success  -1: fail
 *@*/

extern int32_t LDL_decomposition(matrix_t *A, matrix_t *L, matrix_t *D);

#endif