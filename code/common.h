/**
  ***************************************(C) COPYRIGHT 2020 Wyatt Wu***********************************
  * @file        common.c/h
  * @brief       This head/source file is used to redefine data type to adapt other platform
  * @note
  * @history
  * Version      Date            Author          Modification
  * V1.0.0       Nov-07-2020     Wyatt Wu        1. build this file and some foundational content.
  * V1.0.1       Nov-18-2020     Wyatt Wu        1. move content to lib.c/h
  ***************************************(C) COPYRIGHT 2020 Wyatt Wu***********************************
*/
#ifndef _COMMMON_H_
#define _COMMMON_H_

#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>

#define SQR(x) ((x) * (x))
/* redefine data type to adapt other platform */
typedef signed char         int8_t;
typedef unsigned char       uint8_t;
typedef signed short        int16_t;
typedef unsigned short      uint16_t;
typedef int                 int32_t;
typedef unsigned int        uint32_t;
typedef long long           int64_t;
typedef unsigned long long  uint64_t;
typedef float               fp32;
typedef double              fp64;
typedef unsigned char       bool_t;
typedef unsigned char       RETURN_STATUS;

typedef enum
{
    RET_FAIL    = 0,
    RET_SUCCESS = 1,
} ret_status_t;

typedef struct
{
    bool_t  is_open;
    FILE    *fp;
} file_t;

typedef struct
{
  file_t logger;
  file_t out_pos;
} files_manager_t;

#endif