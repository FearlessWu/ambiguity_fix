
#include "Matrix.h"

static void matcpy(fp64* A, const fp64* B, int32_t n, int32_t m)
{
    memcpy(A, B, sizeof(fp64) * n * m);
}

static fp64* mat(int32_t n, int32_t m)
{
    fp64* p;

    if (n <= 0 || m <= 0) return NULL;
    if (!(p = (fp64*)malloc(sizeof(fp64) * n * m))) {
        printf("matrix memory allocation error: n=%d,m=%d\n", n, m);
    }
    return p;
}

static  int32_t* imat(int32_t n, int32_t m)
{
    int32_t* p;

    if (n <= 0 || m <= 0) return NULL;
    if (!(p = (int32_t*)malloc(sizeof(int32_t) * n * m)))
    {
        // TODO: print32_t log
        printf("int32_teger matrix memory allocation error: n=%d,m=%d\n", n, m);
    }
    return p;
}

static void lubksb(const fp64* A, int32_t n, const int32_t* indx, fp64* b)
{
    fp64 s;
    int32_t i, ii = -1, ip, j;

    for (i = 0; i < n; i++)
    {
        ip = indx[i];
        s = b[ip];
        b[ip] = b[i];
        if (ii >= 0)
        {
            for (j = ii; j < i; j++)
            {
                s -= A[i + j * n] * b[j];
            }
        }
        else if (s)
        {
            ii = i;
        }
        b[i] = s;
    }
    for (i = n - 1; i >= 0; i--)
    {
        s = b[i];
        for (j = i + 1; j < n; j++) s -= A[i + j * n] * b[j]; b[i] = s / A[i + i * n];
    }
}

static int32_t ludcmp(fp64 *A, int32_t n, int32_t *indx, fp64 *d)
{
    fp64 big, s, tmp, * vv = mat(n, 1);
    int32_t i, imax = 0, j, k;

    *d = 1.0;
    for (i = 0; i < n; i++)
    {
        big = 0.0; for (j = 0; j < n; j++) if ((tmp = fabs(A[i + j * n])) > big) big = tmp;
        if (big > 0.0) vv[i] = 1.0 / big; else { free(vv); return -1; }
    }
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < j; i++)
        {
            s = A[i + j * n]; for (k = 0; k < i; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
        }
        big = 0.0;
        for (i = j; i < n; i++)
        {
            s = A[i + j * n];
            for (k = 0; k < j; k++) s -= A[i + k * n] * A[k + j * n]; A[i + j * n] = s;
            if ((tmp = vv[i] * fabs(s)) >= big)
            {
                big = tmp;
                imax = i;
            }
        }
        if (j != imax)
        {
            for (k = 0; k < n; k++)
            {
                tmp = A[imax + k * n];
                A[imax + k * n] = A[j + k * n];
                A[j + k * n] = tmp;
            }
            *d = -(*d);
            vv[imax] = vv[j];
        }
        indx[j] = imax;
        if (A[j + j * n] == 0.0)
        {
            free(vv);
            return -1;
        }
        if (j != n - 1)
        {
            tmp = 1.0 / A[j + j * n];
            for (i = j + 1; i < n; i++)
            {
                A[i + j * n] *= tmp;
            }
        }
    }
    free(vv);

    return 0;
}

extern int32_t matinv(fp64* A, int32_t n)
{
    fp64 d, * B;
    int64_t i, j;
    int32_t* indx;

    indx = imat(n, 1);
    B = mat(n, n);

    matcpy(B, A, n, n);
    if (ludcmp(B, n, indx, &d))
    {
        free(indx);
        free(B);

        return -1;
    }
    for (j = 0; j < n; j++)
    {
        for (i = 0; i < n; i++)
        {
            A[i + j * n] = 0.0;
            A[j + j * n] = 1.0;
        }

        lubksb(B, n, indx, A + j * n);
    }
    free(indx); free(B);

    return 0;
}

void matmul(const char *tr, int32_t n, int32_t k, int32_t m, fp64 alpha, const fp64 *A, const fp64 *B, fp64 beta, fp64 *C)
{
    fp64 d;
    int32_t i, j, x, f = tr[0] == 'N' ? (tr[1] == 'N' ? 1 : 2) : (tr[1] == 'N' ? 3 : 4);

    for (i = 0; i < n; i++) for (j = 0; j < k; j++) {
        d = 0.0;
        switch (f)
        {
        case 1: for (x = 0; x < m; x++) d += A[i + x * n] * B[x + j * m]; break;
        case 2: for (x = 0; x < m; x++) d += A[i + x * n] * B[j + x * k]; break;
        case 3: for (x = 0; x < m; x++) d += A[x + i * m] * B[x + j * m]; break;
        case 4: for (x = 0; x < m; x++) d += A[x + i * m] * B[j + x * k]; break;
        }
        if (beta == 0.0) C[i + j * n] = alpha * d; else C[i + j * n] = alpha * d + beta * C[i + j * n];
    }
}

static RETURN_STATUS trans_matrix_to_rtklib_mat(matrix_t *mat, fp64 *m)
{
    if (!mat->is_valid || mat->col <= 0 || mat->row <= 0 || m == NULL)
    {
        // TODO: print32_t log
        return RET_FAIL;
    }

    // in rtklib, matrix elements are saved according colum
    uint32_t i, j;

    for (j = 0; j < mat->col; ++j)
    {
        for (i = 0; i < mat->row; ++i)
        {
            m[j * mat->row + i] = mat->element[i][j];
        }
    }

    return RET_SUCCESS;
}

static RETURN_STATUS trans_rtklib_mat_to_matrix(fp64 *m, matrix_t *mat)
{
    if (!mat->is_valid || mat->col <= 0 || mat->row <= 0 || m == NULL)
    {
        // TODO: print32_t log
        return RET_FAIL;
    }

    uint32_t i, j;
    for (j = 0; j < mat->col; ++j)
    {
        for (i = 0; i < mat->row; ++i)
        {
            mat->element[i][j] = m[j * mat->row + i];
        }
    }
    return RET_SUCCESS;
}

RETURN_STATUS matrix_mlt(matrix_t *mat_in_1, matrix_t *mat_in_2, matrix_t *mat_out)
{
    if ((mat_in_1->row != mat_out->row) || (mat_in_1->col != mat_in_2->row) 
     || (mat_in_2->col != mat_out->col) || (!mat_in_1->is_valid) || (!mat_in_2->is_valid) 
     || (!mat_out->is_valid) || (mat_in_1->col <= 0) || (mat_in_1->row <= 0)
     || (mat_in_1->col <= 0) || (mat_in_2->col <=0))
     {
         // TODO: print32_t log
         return RET_FAIL;
     }

    fp64 *rtklib_mat_in_1;
    fp64 *rtklib_mat_in_2;
    fp64 *rtklib_mat_out;
    if (!(rtklib_mat_in_1 = (fp64*)malloc(sizeof(fp64) * mat_in_1->col * mat_in_1->row)))
    {
        return RET_FAIL;
    }
    if (!(rtklib_mat_in_2 = (fp64*)malloc(sizeof(fp64) * mat_in_2->col * mat_in_2->row)))
    {
        free(rtklib_mat_in_1);

        return RET_FAIL;
    }
    if (!(rtklib_mat_out = (fp64*)malloc(sizeof(fp64) * mat_out->col * mat_out->row)))
    {
        free(rtklib_mat_in_1);
        free(rtklib_mat_in_2);

        return RET_FAIL;
    }

    trans_matrix_to_rtklib_mat(mat_in_1, rtklib_mat_in_1);
    trans_matrix_to_rtklib_mat(mat_in_2, rtklib_mat_in_2);
    trans_matrix_to_rtklib_mat(mat_out, rtklib_mat_out);

    matmul("NN", mat_in_1->row, mat_in_2->col, mat_in_1->col, 1.0, rtklib_mat_in_1, rtklib_mat_in_2, 0.0, rtklib_mat_out);

    trans_rtklib_mat_to_matrix(rtklib_mat_out,  mat_out);

    free(rtklib_mat_in_1);
    free(rtklib_mat_in_2);
    free(rtklib_mat_out);

    return RET_SUCCESS;
}

RETURN_STATUS matrix_init(matrix_t *matrix, const uint32_t row, const uint32_t col)
{
    if (row <= 0 || col <= 0)
    {
        // TODO: send warning to log
        return RET_FAIL;
    }

    uint32_t i;
    uint32_t j;
    if (matrix->element = (fp64**)malloc(sizeof(fp64*) * row))
    {
        for (i = 0; i < row; ++i)
        {
            if (!(matrix->element[i] = (fp64*)malloc(sizeof(fp64) * col)))
            {
                for (j = 0; j < i; ++j)
                {
                    free(matrix->element[j]);
                }
                free(matrix->element);
                matrix->is_valid = 0;
                matrix->col = 0;
                matrix->row = 0;

                // TODO: print32_t log
                return RET_FAIL;
            }
        }
    }
    else
    {
        matrix->col      = 0;
        matrix->row      = 0;
        matrix->is_valid = 0;
        // TODO: print32_t log
        return RET_FAIL;
    }
    
    matrix->col      = col;
    matrix->row      = row;
    matrix->is_valid = 1;

    for (i = 0; i < row; ++i)
    {
        for (j = 0; j < col; ++j)
        {
            matrix->element[i][j] = 0.0;
        }
    }

    return RET_SUCCESS;
}

RETURN_STATUS matrix_free(matrix_t *matrix)
{
    if (matrix->row <= 0 || matrix->col <= 0 || matrix->is_valid == 0)
    {
        return RET_FAIL;
    }

    uint32_t i;

    for (i = 0; i < matrix->row; ++i)
    {
        free(matrix->element[i]);
    }
    free(matrix->element);
    matrix->row      = 0;
    matrix->col      = 0;
    matrix->element  = NULL;
    matrix->is_valid = 0;
    
    return RET_SUCCESS;
}

RETURN_STATUS matrix_resize(matrix_t* mat, const uint32_t row, const uint32_t col)
{
    if (mat->is_valid || mat->col <= 0 || mat->row <= 0)
    {
        if (!matrix_init(mat, row, col))
        {
            return RET_FAIL;
        }

        return RET_SUCCESS;
    }
    else
    {
        matrix_free(mat);
        if (!matrix_init(mat, row, col))
        {
            return RET_FAIL;
        }

        return RET_SUCCESS;
    }
}

RETURN_STATUS matrix_add(matrix_t *mat_in_1, matrix_t *mat_in_2, matrix_t *mat_out)
{
    if ((mat_in_1->row != mat_in_2->row) || (mat_in_1->col != mat_in_2->col) 
     || (mat_in_1->row != mat_out->row)  || (mat_in_1->col != mat_out->col) 
     || (!mat_in_1->is_valid) || (!mat_in_2->is_valid) || (!mat_out->is_valid)
     || (mat_in_1->col <= 0)  || (mat_in_1->row <= 0))
     {
         // TODO: print32_t log
         return RET_FAIL;
     }

     uint32_t i;
     uint32_t j;
     for (i = 0; i < mat_out->row; ++i)
     {
         for (j = 0; j <mat_out->col; ++j)
         {
             mat_out->element[i][j] = mat_in_1->element[i][j] + mat_in_2->element[i][j];
         }
     }

     return RET_SUCCESS;
}

RETURN_STATUS matrix_miu(matrix_t *mat_in_1, matrix_t *mat_in_2, matrix_t *mat_out)
{
    if ((mat_in_1->row != mat_in_2->row) || (mat_in_1->col != mat_in_2->col) 
     || (mat_in_1->row != mat_out->row)  || (mat_in_1->col != mat_out->col) 
     || (!mat_in_1->is_valid) || (!mat_in_2->is_valid) || (!mat_out->is_valid)
     || (mat_in_1->col <= 0)  || (mat_in_1->row <= 0))
     {
         // TODO: print32_t log
         return RET_FAIL;
     }

     uint32_t i;
     uint32_t j;
     for (i = 0; i < mat_out->row; ++i)
     {
         for (j = 0; j <mat_out->col; ++j)
         {
             mat_out->element[i][j] = mat_in_1->element[i][j] - mat_in_2->element[i][j];
         }
     }

     return RET_SUCCESS;
}

RETURN_STATUS matrix_inv(matrix_t *mat_in, matrix_t *mat_out)
{
    if (mat_in->row  <= 0 || mat_in->col  <= 0 || !mat_in->is_valid  || mat_in->col  !=mat_in->row
     || mat_out->row <= 0 || mat_out->col <= 0 || !mat_out->is_valid || mat_out->col != mat_out->row
     || mat_in->col != mat_out->col)
    {
        // TODO: print32_t log
        return RET_FAIL;
    }
    
    fp64* temp = (fp64*)malloc(sizeof(fp64) * mat_in->col * mat_in->row);
    if (temp == NULL)
    {
        // TODO: print32_t log
        return RET_FAIL;
    }

    trans_matrix_to_rtklib_mat(mat_in, temp);

    if (matinv(temp, mat_in->row))
    {
        // TODO: print32_t log
        return RET_FAIL;
    }

    trans_rtklib_mat_to_matrix(temp, mat_out);

    free(temp);

    return RET_SUCCESS;
}

RETURN_STATUS matrix_trs(matrix_t *mat_in, matrix_t *mat_out)
{
    if (mat_in->row  <= 0 || mat_in->col  <= 0 || !mat_in->is_valid  || mat_in->col !=mat_out->row
     || mat_out->row <= 0 || mat_out->col <= 0 || !mat_out->is_valid || mat_in->row != mat_out->col)
    {
        // TODO: print log
        return RET_FAIL;
    }

    uint32_t i, j;

    for (i = 0; i < mat_in->row; ++i)
    {
        for (j = 0; j < mat_in->col; ++j)
        {
            mat_out->element[j][i] = mat_in->element[i][j];
        }
    }

    return RET_SUCCESS;
}

RETURN_STATUS matrix_copy(const matrix_t *mat_in, matrix_t *mat_out)
{
    if (mat_in->row <= 0 || mat_in->col <= 0 || !mat_in->is_valid || mat_out->row <= 0 || mat_out->col <= 0 || !mat_out->is_valid)
    {
        // TODO: report error
        return RET_FAIL;
    }

    if (mat_in->col != mat_out->col || mat_in->row != mat_out->row)
    {
        // TODO: report error
        return RET_FAIL;
    }

    uint32_t i, j;
    for (i = 0; i < mat_in->row; ++i)
    {
        for (j = 0; j < mat_in->col; ++j)
        {
            mat_out->element[i][j] = mat_in->element[i][j];
        }
    }
  
    return RET_SUCCESS;
}

RETURN_STATUS matrix_extend_col(matrix_t *mat, uint32_t n)
{
    if (mat->col <= 0 || mat->row <= 0 || !mat->is_valid || n <= 0)
    {
        // TODO: report error
        return RET_FAIL;
    }

    matrix_t tmp;
    if (!matrix_init(&tmp, mat->row, mat->col))
    {
        // TODO: report error
        return RET_FAIL;
    }
    matrix_copy(mat, &tmp);
    if (!matrix_resize(mat, mat->row, (mat->col + n)))
    {
        // TODO: report error
        return RET_FAIL;
    }

    uint32_t i, j;
    for (i = 0; i < tmp.row; ++i)
    {
        for (j = 0; j < tmp.col; ++j)
        {
            mat->element[i][j] = tmp.element[i][j];
        }
    }
    matrix_free(&tmp);

    return RET_SUCCESS;
}

RETURN_STATUS matrix_extend_row(matrix_t *mat, uint32_t n)
{
    if (mat->col <= 0 || mat->row <= 0 || !mat->is_valid || n <= 0)
    {
        // TODO: report error
        return RET_FAIL;
    }

    matrix_t tmp;
    if (!matrix_init(&tmp, mat->row, mat->col))
    {
        // TODO: report error
        return RET_FAIL;
    }
    matrix_copy(mat, &tmp);
    if (!matrix_resize(mat, (mat->row + n), mat->col))
    {
        // TODO: report error
        return RET_FAIL;
    }

    uint32_t i, j;
    for (i = 0; i < tmp.row; ++i)
    {
        for (j = 0; j < tmp.col; ++j)
        {
            mat->element[i][j] = tmp.element[i][j];
        }
    }
    matrix_free(&tmp);

    return RET_SUCCESS;
}

void matrix_print(matrix_t matrix)
{
    printf("\n");
    uint32_t i;
    uint32_t j;
    for (i = 0; i < matrix.row; ++i)
    {
        for (j = 0; j < matrix.col; ++j)
        {
            if (j == matrix.col - 1)
            {
                printf("%8.3f\n", matrix.element[i][j]);
            }
            else
            {
                printf("%8.3f", matrix.element[i][j]);
            }
        }
    }
}

void matrix_log(matrix_t matrix, file_t *logger, char *message)
{
    uint32_t i;
    uint32_t j;

    if (!logger->is_open)
    {
        return;
    }
    fprintf(logger->fp, "%s\n", message);
    for (i = 0; i < matrix.row; ++i)
    {
        for (j = 0; j < matrix.col; ++j)
        {
            if (j == matrix.col - 1)
            {
                fprintf(logger->fp, "%8.3f\n", matrix.element[i][j]);
            }
            else
            {
                fprintf(logger->fp, "%8.3f", matrix.element[i][j]);
            }
        }
    }

    fflush(logger->fp);
}

void matrix_eye(matrix_t *mat, uint32_t n)
{
    if (matrix_init(mat, n, n))
    {
        for (uint32_t i = 0; i < n; ++i)
        {
            mat->element[i][i] = 1;
        }
    }

}