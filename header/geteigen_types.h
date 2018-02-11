/*
 * geteigen_types.h
 *
 * Code generation for function 'geteigen'
 *
 * C source code generated on: Sat Feb 11 19:17:44 2017
 *
 */

#ifndef __GETEIGEN_TYPES_H__
#define __GETEIGEN_TYPES_H__

/* Type Definitions */
#ifndef struct_emxArray__common
#define struct_emxArray__common
typedef struct emxArray__common
{
    void *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray__common;
#endif
#ifndef struct_emxArray_creal32_T
#define struct_emxArray_creal32_T
typedef struct emxArray_creal32_T
{
    creal32_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_creal32_T;
#endif
#ifndef struct_emxArray_int32_T
#define struct_emxArray_int32_T
typedef struct emxArray_int32_T
{
    int32_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_int32_T;
#endif
#ifndef struct_emxArray_real32_T
#define struct_emxArray_real32_T
typedef struct emxArray_real32_T
{
    real32_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_real32_T;
#endif

#endif
/* End of code generation (geteigen_types.h) */
