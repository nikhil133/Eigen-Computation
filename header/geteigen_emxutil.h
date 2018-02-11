/*
 * geteigen_emxutil.h
 *
 * Code generation for function 'geteigen_emxutil'
 *
 * C source code generated on: Sat Feb 11 19:17:44 2017
 *
 */

#ifndef __GETEIGEN_EMXUTIL_H__
#define __GETEIGEN_EMXUTIL_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "geteigen_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void b_emxInit_creal32_T(emxArray_creal32_T **pEmxArray, int32_T numDimensions);
extern void b_emxInit_real32_T(emxArray_real32_T **pEmxArray, int32_T numDimensions);
extern void c_emxInit_real32_T(emxArray_real32_T **pEmxArray, int32_T numDimensions);
extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
extern void emxFree_creal32_T(emxArray_creal32_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real32_T(emxArray_real32_T **pEmxArray);
extern void emxInit_creal32_T(emxArray_creal32_T **pEmxArray, int32_T numDimensions);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions);
extern void emxInit_real32_T(emxArray_real32_T **pEmxArray, int32_T numDimensions);
#endif
/* End of code generation (geteigen_emxutil.h) */
