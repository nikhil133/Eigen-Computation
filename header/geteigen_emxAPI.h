/*
 * geteigen_emxAPI.h
 *
 * Code generation for function 'geteigen_emxAPI'
 *
 * C source code generated on: Sat Feb 11 19:17:45 2017
 *
 */

#ifndef __GETEIGEN_EMXAPI_H__
#define __GETEIGEN_EMXAPI_H__
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
extern emxArray_real32_T *emxCreateND_real32_T(int32_T numDimensions, int32_T *size);
extern emxArray_real32_T *emxCreateWrapperND_real32_T(real32_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_real32_T *emxCreateWrapper_real32_T(real32_T *data, int32_T rows, int32_T cols);
extern emxArray_real32_T *emxCreate_real32_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_real32_T(emxArray_real32_T *emxArray);
#endif
/* End of code generation (geteigen_emxAPI.h) */
