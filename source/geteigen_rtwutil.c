/*
 * geteigen_rtwutil.c
 *
 * Code generation for function 'geteigen_rtwutil'
 *
 * C source code generated on: Sat Feb 11 19:17:44 2017
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "geteigen.h"
#include "geteigen_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
real32_T rt_hypotf_snf(real32_T u0, real32_T u1)
{
  real32_T y;
  real32_T a;
  a = (real32_T)fabs(u0);
  y = (real32_T)fabs(u1);
  if (a < y) {
    a /= y;
    y *= (real32_T)sqrt(a * a + 1.0F);
  } else if (a > y) {
    y /= a;
    y = a * (real32_T)sqrt(y * y + 1.0F);
  } else if (rtIsNaNF(y)) {
  } else {
    y = a * 1.41421354F;
  }

  return y;
}

/* End of code generation (geteigen_rtwutil.c) */
