/*
 * geteigen_emxAPI.c
 *
 * Code generation for function 'geteigen_emxAPI'
 *
 * C source code generated on: Sat Feb 11 19:17:44 2017
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "geteigen.h"
#include "geteigen_emxAPI.h"
#include "geteigen_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
emxArray_real32_T *emxCreateND_real32_T(int32_T numDimensions, int32_T *size)
{
  emxArray_real32_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_real32_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real32_T *)calloc((uint32_T)numEl, sizeof(real32_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}

emxArray_real32_T *emxCreateWrapperND_real32_T(real32_T *data, int32_T
  numDimensions, int32_T *size)
{
  emxArray_real32_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  c_emxInit_real32_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;
  for (i = 0; i <= loop_ub; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real32_T *emxCreateWrapper_real32_T(real32_T *data, int32_T rows,
  int32_T cols)
{
  emxArray_real32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_real32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}

emxArray_real32_T *emxCreate_real32_T(int32_T rows, int32_T cols)
{
  emxArray_real32_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  c_emxInit_real32_T(&emx, 2);
  numEl = 1;
  for (i = 0; i < 2; i++) {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real32_T *)calloc((uint32_T)numEl, sizeof(real32_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}

void emxDestroyArray_real32_T(emxArray_real32_T *emxArray)
{
  emxFree_real32_T(&emxArray);
}

/* End of code generation (geteigen_emxAPI.c) */
