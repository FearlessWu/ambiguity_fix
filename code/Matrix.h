/**
  **********************************(C) COPYRIGHT 2020 Wyatt Wu***********************************
  * @file        Matrix.c/h
  * @brief       This is a Matrix header file, all about matrix operations are written here. matrix 
  *              multiply and inverse opetation actually rewrite from RTKLIB. Using the form of 
  *              structure to manage the matrix, try to avoid the wrong use of the matrix. For convenience,
  *              we used a two-dimensional pointer performs matrix operation. 
                 
  * @note
  * @history
  * Version      Date               Author               Modification
  * V1.0.0       Nov-29-2020        Wyatt Wu             1. establih this file
  ************************************************************************************************
*/
#ifndef _MATRIX_H_
#define _MATRIX_H_

#include "common.h"

/**
*@brief Matrix struct
*@note  is_valid means whether Matrix is ready to use or not
*/
typedef struct
{
    uint32_t row;
    uint32_t col;
    uint8_t  is_valid;
    fp64     **element;
} matrix_t;

/**
 *@brief            multiply matrix
 *@author           quote RTKLIB
 *@param[in]        tr   : transpose flag. the first char means the first input matrix transpose state,
 *                         the second char means the second input matrix transpose state
 *@param[in]        n    : row size of matrix A
 *@param[in]        k    : col size of matrix B
 *@param[in]        m    : col/row size of matrix A/B
 *@param[in]        alpha: mlutiply factor
 *@param[in]        A    : front matrix
 *@param[in]        B    : back matrix
 *@param[in]        beta : mlutiply factor
 *@param[inout]     C    : multiply matrix result
 *@retval           none
 *@note             if tr = "NN", C = alpha * (A  * B)  + beta * C
 *                  if tr = "TN", C = alpha * (AT * B)  + beta * C
 *                  if tr = "NT", C = alpha * (A  * BT) + beta * C
 *                  if tr = "TT", C = alpha * (AT * BT) + beta * C
 **/
extern void matmul(const char* tr, int32_t n, int32_t k, int32_t m, fp64 alpha, const fp64 *A, 
                   const fp64 *B, fp64 beta, fp64 *C);

/**
 *@brief            inverse operation of matrix
 *@author           quote RTKLIB
 *@param[inout]     A: input matrix and output matrix
 *@param[in]        n: row and col of A
 *@retval           0: success  -1: fail
 *@*/
extern int32_t matinv(fp64 *A, int32_t n);

/**
 *@brief            Initialize matrix, its elements are set to 0. Matrix must be initialized before using.
 *@author           wyatt.wu
 *@param[inout]     row   : matrix row
 *@param[in]        col   : matrix col
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_init(matrix_t *matrix, const uint32_t row, const uint32_t col);

/**
 *@brief            resize matix size
 *@author           wyatt.wu
 *@param[inout]     mat   : matrix need to be resized
 *@param[in]        row   : matrix row
 *@param[in]        col   : matrix col
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_resize(matrix_t *mat, const uint32_t row, const uint32_t col);

/**
 *@brief            free matrix memory. if matrix has been initialized, it must be free.
 *@author           wyatt.wu
 *@param[in]        matrix: matrix need to be free
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_free(matrix_t *matrix);

/**
 *@brief            matrix multiply operation
 *@author           wyatt.wu
 *@note:            mat_out = mat_in_1 * mat_in_2
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_mlt(matrix_t *mat_in_1, matrix_t *mat_in_2, matrix_t *mat_out);

/**
 *@brief            matrix tranpose  operation
 *@author           wyatt.wu
 *@note:            mat_out = mat_in ^ T
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_trs(matrix_t *mat_in, matrix_t *mat_out);

/**
 *@brief            matrix inverse  operation
 *@author           wyatt.wu
 *@note:            mat_out = mat_in ^ -1
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_inv(matrix_t *mat_in, matrix_t *mat_out);

/**
 *@brief            matrix add  operation
 *@author           wyatt.wu
 *@note:            mat_out = mat_in_1 + mat_in_2
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_add(matrix_t *mat_in_1, matrix_t *mat_in_2, matrix_t *mat_out);

/**
 *@brief            matrix add  operation
 *@author           wyatt.wu
 *@note:            mat_out = mat_in_1 - mat_in_2
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_miu(matrix_t *mat_in_1, matrix_t *mat_in_2, matrix_t *mat_out);

/**
 *@brief            extend n colums of matrix, and original elements still store in matrix. the new colums elements are set to 0.
 *@author           wyatt.wu
 *@note:            mat_out(row, col) --->   mat_out(row, col + n)
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_extend_col(matrix_t *mat, uint32_t n);

/**
 *@brief            extend n row of matrix, and original elements still store in matrix. the new rows elements are set to 0.
 *@author           wyatt.wu
 *@note:            mat_out(row, col) --->   mat_out(row + n, col)
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_extend_row(matrix_t *mat, uint32_t n);

/**
 *@brief            matrix copy  operation
 *@author           wyatt.wu
 *@note:            mat_out = mat_in;
 *@retval           RET_FAIL: fail, RET_SUCCESS: sucesss
 **/
extern RETURN_STATUS matrix_copy(const matrix_t *mat_in, matrix_t *mat_out);

/**
 *@brief            print matrix to screen according to its row and col,
 *@author           wyatt.wu
 *@note:
 *@retval           none
 **/
extern void matrix_print(matrix_t matrix);

/**
 *@brief            print matrix to log file according to its row and col,
 *@author           wyatt.wu
 *@note:
 *@retval           none
 **/
extern void matrix_log(matrix_t matrix, file_t *logger, char *message);

/**
 *@brief            create a unit matrix
 *@author           wyatt.wu
 *@param[in]        mat: input and output matrix
 *@param[in]        n: demenstion of mat
 *@note:
 *@retval           none
 **/
extern void matrix_eye(matrix_t *mat, uint32_t n);
#endif
