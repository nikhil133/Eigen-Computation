/*
 * sqrt.c
 *
 * Code generation for function 'sqrt'
 *
 * C source code generated on: Sat Feb 11 19:17:44 2017
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "geteigen.h"
#include "sqrt.h"
#include "geteigen_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_sqrt(creal32_T *x)
{
  real32_T absxi;
  real32_T absxr;
  if (x->im == 0.0F) {
    if (x->re < 0.0F) {
      absxi = 0.0F;
      absxr = (real32_T)sqrt((real32_T)fabs(x->re));
    } else {
      absxi = (real32_T)sqrt(x->re);
      absxr = 0.0F;
    }
  } else if (x->re == 0.0F) {
    if (x->im < 0.0F) {
      absxi = (real32_T)sqrt(-x->im / 2.0F);
      absxr = -absxi;
    } else {
      absxi = (real32_T)sqrt(x->im / 2.0F);
      absxr = absxi;
    }
  } else if (rtIsNaNF(x->re) || rtIsNaNF(x->im)) {
    absxi = ((real32_T)rtNaN);
    absxr = ((real32_T)rtNaN);
  } else if (rtIsInfF(x->im)) {
    absxi = ((real32_T)rtInf);
    absxr = x->im;
  } else if (rtIsInfF(x->re)) {
    if (x->re < 0.0F) {
      absxi = 0.0F;
      absxr = ((real32_T)rtInf);
    } else {
      absxi = ((real32_T)rtInf);
      absxr = 0.0F;
    }
  } else {
    absxr = (real32_T)fabs(x->re);
    absxi = (real32_T)fabs(x->im);
    if ((absxr > 8.50705867E+37F) || (absxi > 8.50705867E+37F)) {
      absxr *= 0.5F;
      absxi *= 0.5F;
      absxi = rt_hypotf_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = (real32_T)sqrt(absxi) * (real32_T)sqrt(1.0F + absxr / absxi);
      } else {
        absxi = (real32_T)sqrt(absxi) * 1.41421354F;
      }
    } else {
      absxi = (real32_T)sqrt((rt_hypotf_snf(absxr, absxi) + absxr) * 0.5F);
    }

    if (x->re > 0.0F) {
      absxr = 0.5F * (x->im / absxi);
    } else {
      if (x->im < 0.0F) {
        absxr = -absxi;
      } else {
        absxr = absxi;
      }

      absxi = 0.5F * (x->im / absxr);
    }
  }

  x->re = absxi;
  x->im = absxr;
}

/* End of code generation (sqrt.c) */
