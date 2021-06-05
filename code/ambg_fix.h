#ifndef _AMBG_FIX_H_
#define _AMBG_FIX_H_
#include "Matrix.h"
/**
 *@brief            LDL decompodition
 *@param[in]        A: symmetric matrix
 *@param[out]       L: lower triangular matrix
 *@param[out]       D: dialog matrix
 *@retval           1: success  0: fail
 *@*/

extern int8_t LDL_decomposition(matrix_t *A, matrix_t *L, matrix_t *D);

/*
 *@brief            guass transformation
 *@param[in]        L: lower triangluar matrix
 *@param[out]       Z: Z transformation matrix
 *@retval           1: success  0: fial
 */
extern int8_t gauss_transform(matrix_t *L, matrix_t *Z);

#endif