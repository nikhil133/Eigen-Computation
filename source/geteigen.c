/*
 * geteigen.c
 *
 * Code generation for function 'geteigen'
 *
 * C source code generated on: Sat Feb 11 19:17:44 2017
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "geteigen.h"
#include "geteigen_emxutil.h"
#include "sqrt.h"
#include "geteigen_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_matlab_zlartg(const creal32_T f, const creal32_T g, real32_T
  *cs, creal32_T *sn);
static int32_T div_s32_floor(int32_T numerator, int32_T denominator);
static creal32_T eml_div(const creal32_T x, real_T y);
static void eml_matlab_zggbal(emxArray_creal32_T *A, int32_T *ilo, int32_T *ihi,
  emxArray_int32_T *rscale);
static void eml_matlab_zggev(emxArray_creal32_T *A, real_T *info,
  emxArray_creal32_T *alpha1, emxArray_creal32_T *beta1, emxArray_creal32_T *V);
static void eml_matlab_zhgeqz(emxArray_creal32_T *A, int32_T ilo, int32_T ihi,
  emxArray_creal32_T *Z, real_T *info, emxArray_creal32_T *alpha1,
  emxArray_creal32_T *beta1);
static real32_T eml_matlab_zlanhs(const emxArray_creal32_T *A, int32_T ilo,
  int32_T ihi);
static void eml_matlab_zlartg(const creal32_T f, const creal32_T g, real32_T *cs,
  creal32_T *sn, creal32_T *r);
static void eml_matlab_ztgevc(const emxArray_creal32_T *A, emxArray_creal32_T *V);

/* Function Definitions */
static void b_eml_matlab_zlartg(const creal32_T f, const creal32_T g, real32_T
  *cs, creal32_T *sn)
{
  real32_T scale;
  real32_T f2s;
  real32_T g2;
  real32_T fs_re;
  real32_T fs_im;
  real32_T gs_re;
  real32_T gs_im;
  boolean_T guard1 = FALSE;
  real32_T g2s;
  scale = (real32_T)fabs(f.re);
  f2s = (real32_T)fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = (real32_T)fabs(g.re);
  g2 = (real32_T)fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = FALSE;
  if (scale >= 5.49755814E+11F) {
    do {
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    guard1 = TRUE;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
    } else {
      do {
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0F > g2) {
      f2s = 1.0F;
    }

    if (scale <= f2s * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        scale = rt_hypotf_snf((real32_T)fabs(gs_re), (real32_T)fabs(gs_im));
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = (real32_T)sqrt(g2);
        *cs = rt_hypotf_snf((real32_T)fabs(fs_re), (real32_T)fabs(fs_im)) / g2s;
        f2s = (real32_T)fabs(f.re);
        g2 = (real32_T)fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0F) {
          scale = rt_hypotf_snf((real32_T)fabs(f.re), (real32_T)fabs(f.im));
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2s = 5.49755814E+11F * f.re;
          g2 = 5.49755814E+11F * f.im;
          scale = rt_hypotf_snf((real32_T)fabs(f2s), (real32_T)fabs(g2));
          fs_re = f2s / scale;
          fs_im = g2 / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = (real32_T)sqrt(1.0F + g2 / scale);
      *cs = 1.0F / f2s;
      scale += g2;
      fs_re = f2s * fs_re / scale;
      fs_im = f2s * fs_im / scale;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

static int32_T div_s32_floor(int32_T numerator, int32_T denominator)
{
  int32_T quotient;
  uint32_T absNumerator;
  uint32_T absDenominator;
  int32_T quotientNeedsNegation;
  uint32_T tempAbsQuotient;
  if (denominator == 0) {
    quotient = numerator >= 0 ? MAX_int32_T : MIN_int32_T;
  } else {
    absNumerator = (uint32_T)(numerator >= 0 ? numerator : -numerator);
    absDenominator = (uint32_T)(denominator >= 0 ? denominator : -denominator);
    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if ((uint32_T)quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > (uint32_T)0) {
        tempAbsQuotient++;
      }
    }

    quotient = (uint32_T)quotientNeedsNegation ? -(int32_T)tempAbsQuotient :
      (int32_T)tempAbsQuotient;
  }

  return quotient;
}

static creal32_T eml_div(const creal32_T x, real_T y)
{
  creal32_T z;
  if (x.im == 0.0F) {
    z.re = x.re / (real32_T)y;
    z.im = 0.0F;
  } else if (x.re == 0.0F) {
    z.re = 0.0F;
    z.im = x.im / (real32_T)y;
  } else {
    z.re = x.re / (real32_T)y;
    z.im = x.im / (real32_T)y;
  }

  return z;
}

static void eml_matlab_zggbal(emxArray_creal32_T *A, int32_T *ilo, int32_T *ihi,
  emxArray_int32_T *rscale)
{
  emxArray_creal32_T *b_A;
  int32_T n;
  int32_T nzcount;
  int32_T ii;
  int32_T exitg2;
  int32_T i;
  int32_T j;
  boolean_T found;
  boolean_T exitg5;
  boolean_T exitg6;
  real_T A_re;
  real_T A_im;
  boolean_T guard2 = FALSE;
  real32_T atmp_re;
  real32_T atmp_im;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = FALSE;
  emxInit_creal32_T(&b_A, 2);
  n = A->size[0];
  nzcount = rscale->size[0];
  rscale->size[0] = n;
  emxEnsureCapacity((emxArray__common *)rscale, nzcount, (int32_T)sizeof(int32_T));
  ii = n - 1;
  for (nzcount = 0; nzcount <= ii; nzcount++) {
    rscale->data[nzcount] = 0;
  }

  *ilo = 1;
  *ihi = n;
  if (n <= 1) {
    *ihi = 1;
    rscale->data[0] = 1;
  } else {
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      found = FALSE;
      ii = *ihi;
      exitg5 = FALSE;
      while ((exitg5 == 0U) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = *ihi;
        n = 1;
        exitg6 = FALSE;
        while ((exitg6 == 0U) && (n <= *ihi)) {
          A_re = A->data[(ii + A->size[0] * (n - 1)) - 1].re;
          A_im = A->data[(ii + A->size[0] * (n - 1)) - 1].im;
          guard2 = FALSE;
          if ((A_re != 0.0) || (A_im != 0.0) || (ii == n)) {
            if (nzcount == 0) {
              j = n;
              nzcount = 1;
              guard2 = TRUE;
            } else {
              nzcount = 2;
              exitg6 = TRUE;
            }
          } else {
            guard2 = TRUE;
          }

          if (guard2 == TRUE) {
            n++;
          }
        }

        if (nzcount < 2) {
          found = TRUE;
          exitg5 = TRUE;
        } else {
          ii--;
        }
      }

      if (!found) {
        exitg2 = 2;
      } else {
        n = A->size[0];
        if (i != *ihi) {
          for (nzcount = 0; nzcount + 1 <= n; nzcount++) {
            atmp_re = A->data[(i + A->size[0] * nzcount) - 1].re;
            atmp_im = A->data[(i + A->size[0] * nzcount) - 1].im;
            A->data[(i + A->size[0] * nzcount) - 1] = A->data[(*ihi + A->size[0]
              * nzcount) - 1];
            A->data[(*ihi + A->size[0] * nzcount) - 1].re = atmp_re;
            A->data[(*ihi + A->size[0] * nzcount) - 1].im = atmp_im;
          }
        }

        if (j != *ihi) {
          for (nzcount = 0; nzcount + 1 <= *ihi; nzcount++) {
            atmp_re = A->data[nzcount + A->size[0] * (j - 1)].re;
            atmp_im = A->data[nzcount + A->size[0] * (j - 1)].im;
            A->data[nzcount + A->size[0] * (j - 1)] = A->data[nzcount + A->size
              [0] * (*ihi - 1)];
            A->data[nzcount + A->size[0] * (*ihi - 1)].re = atmp_re;
            A->data[nzcount + A->size[0] * (*ihi - 1)].im = atmp_im;
          }
        }

        rscale->data[*ihi - 1] = j;
        (*ihi)--;
        if (*ihi == 1) {
          rscale->data[0] = 1;
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0U);

    if (exitg2 == 1U) {
    } else {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        found = FALSE;
        n = *ilo;
        exitg3 = FALSE;
        while ((exitg3 == 0U) && (n <= *ihi)) {
          nzcount = 0;
          i = *ihi;
          j = n;
          ii = *ilo;
          exitg4 = FALSE;
          while ((exitg4 == 0U) && (ii <= *ihi)) {
            A_re = A->data[(ii + A->size[0] * (n - 1)) - 1].re;
            A_im = A->data[(ii + A->size[0] * (n - 1)) - 1].im;
            guard1 = FALSE;
            if ((A_re != 0.0) || (A_im != 0.0) || (ii == n)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                guard1 = TRUE;
              } else {
                nzcount = 2;
                exitg4 = TRUE;
              }
            } else {
              guard1 = TRUE;
            }

            if (guard1 == TRUE) {
              ii++;
            }
          }

          if (nzcount < 2) {
            found = TRUE;
            exitg3 = TRUE;
          } else {
            n++;
          }
        }

        if (!found) {
          exitg1 = 1;
        } else {
          n = A->size[0];
          if (i != *ilo) {
            for (nzcount = *ilo - 1; nzcount + 1 <= n; nzcount++) {
              atmp_re = A->data[(i + A->size[0] * nzcount) - 1].re;
              atmp_im = A->data[(i + A->size[0] * nzcount) - 1].im;
              A->data[(i + A->size[0] * nzcount) - 1] = A->data[(*ilo + A->size
                [0] * nzcount) - 1];
              A->data[(*ilo + A->size[0] * nzcount) - 1].re = atmp_re;
              A->data[(*ilo + A->size[0] * nzcount) - 1].im = atmp_im;
            }
          }

          if (j != *ilo) {
            for (nzcount = 0; nzcount + 1 <= *ihi; nzcount++) {
              atmp_re = A->data[nzcount + A->size[0] * (j - 1)].re;
              atmp_im = A->data[nzcount + A->size[0] * (j - 1)].im;
              A->data[nzcount + A->size[0] * (j - 1)] = A->data[nzcount +
                A->size[0] * (*ilo - 1)];
              A->data[nzcount + A->size[0] * (*ilo - 1)].re = atmp_re;
              A->data[nzcount + A->size[0] * (*ilo - 1)].im = atmp_im;
            }
          }

          rscale->data[*ilo - 1] = j;
          (*ilo)++;
          if (*ilo == *ihi) {
            rscale->data[*ilo - 1] = *ilo;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0U);
    }
  }

  emxFree_creal32_T(&b_A);
}

static void eml_matlab_zggev(emxArray_creal32_T *A, real_T *info,
  emxArray_creal32_T *alpha1, emxArray_creal32_T *beta1, emxArray_creal32_T *V)
{
  int32_T n;
  int32_T jrow;
  int32_T jcol;
  int32_T jcolp1;
  real32_T anrm;
  boolean_T exitg1;
  real32_T absxk;
  boolean_T ilascl;
  real32_T anrmto;
  real32_T ctoc;
  boolean_T notdone;
  real32_T cfrom1;
  real32_T cto1;
  real32_T mul;
  emxArray_creal32_T *b_A;
  emxArray_int32_T *rscale;
  emxArray_real32_T *I;
  int32_T ihi;
  int32_T ilo;
  int32_T b_n;
  int32_T i;
  int32_T jrowm1;
  creal32_T c_A;
  creal32_T d_A;
  creal32_T tmp;
  creal32_T s;
  int32_T j;
  real32_T y_re;
  real32_T y_im;
  real_T b_info;
  *info = 0.0;
  n = A->size[0];
  jrow = alpha1->size[0];
  alpha1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, jrow, (int32_T)sizeof(creal32_T));
  jcol = A->size[0] - 1;
  for (jrow = 0; jrow <= jcol; jrow++) {
    alpha1->data[jrow].re = 0.0F;
    alpha1->data[jrow].im = 0.0F;
  }

  jrow = beta1->size[0];
  beta1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)beta1, jrow, (int32_T)sizeof(creal32_T));
  jcol = A->size[0] - 1;
  for (jrow = 0; jrow <= jcol; jrow++) {
    beta1->data[jrow].re = 0.0F;
    beta1->data[jrow].im = 0.0F;
  }

  jcol = A->size[0];
  jrow = V->size[0] * V->size[1];
  V->size[0] = jcol;
  emxEnsureCapacity((emxArray__common *)V, jrow, (int32_T)sizeof(creal32_T));
  jcolp1 = A->size[0];
  jrow = V->size[0] * V->size[1];
  V->size[1] = jcolp1;
  emxEnsureCapacity((emxArray__common *)V, jrow, (int32_T)sizeof(creal32_T));
  jcol = A->size[0] * A->size[0] - 1;
  for (jrow = 0; jrow <= jcol; jrow++) {
    V->data[jrow].re = 0.0F;
    V->data[jrow].im = 0.0F;
  }

  if ((A->size[0] == 0) || (A->size[1] == 0)) {
  } else {
    anrm = 0.0F;
    jrow = A->size[0] * A->size[1];
    jcol = 0;
    exitg1 = FALSE;
    while ((exitg1 == 0U) && (jcol <= jrow - 1)) {
      absxk = rt_hypotf_snf((real32_T)fabs(A->data[jcol].re), (real32_T)fabs
                            (A->data[jcol].im));
      if (rtIsNaNF(absxk)) {
        anrm = ((real32_T)rtNaN);
        exitg1 = TRUE;
      } else {
        if (absxk > anrm) {
          anrm = absxk;
        }

        jcol++;
      }
    }

    if (!((!rtIsInfF(anrm)) && (!rtIsNaNF(anrm)))) {
      jrow = alpha1->size[0];
      alpha1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)alpha1, jrow, (int32_T)sizeof
                        (creal32_T));
      jcol = A->size[0] - 1;
      for (jrow = 0; jrow <= jcol; jrow++) {
        alpha1->data[jrow].re = ((real32_T)rtNaN);
        alpha1->data[jrow].im = 0.0F;
      }

      jrow = beta1->size[0];
      beta1->size[0] = A->size[0];
      emxEnsureCapacity((emxArray__common *)beta1, jrow, (int32_T)sizeof
                        (creal32_T));
      jcol = A->size[0] - 1;
      for (jrow = 0; jrow <= jcol; jrow++) {
        beta1->data[jrow].re = ((real32_T)rtNaN);
        beta1->data[jrow].im = 0.0F;
      }

      jcol = A->size[0];
      jrow = V->size[0] * V->size[1];
      V->size[0] = jcol;
      emxEnsureCapacity((emxArray__common *)V, jrow, (int32_T)sizeof(creal32_T));
      jcolp1 = A->size[0];
      jrow = V->size[0] * V->size[1];
      V->size[1] = jcolp1;
      emxEnsureCapacity((emxArray__common *)V, jrow, (int32_T)sizeof(creal32_T));
      jcol = A->size[0] * A->size[0] - 1;
      for (jrow = 0; jrow <= jcol; jrow++) {
        V->data[jrow].re = ((real32_T)rtNaN);
        V->data[jrow].im = 0.0F;
      }
    } else {
      ilascl = FALSE;
      anrmto = anrm;
      if ((anrm > 0.0F) && (anrm < 9.09494702E-13F)) {
        anrmto = 9.09494702E-13F;
        ilascl = TRUE;
      } else {
        if (anrm > 1.09951163E+12F) {
          anrmto = 1.09951163E+12F;
          ilascl = TRUE;
        }
      }

      if (ilascl) {
        absxk = anrm;
        ctoc = anrmto;
        notdone = TRUE;
        while (notdone) {
          cfrom1 = absxk * 1.97215226E-31F;
          cto1 = ctoc / 5.0706024E+30F;
          if ((cfrom1 > ctoc) && (ctoc != 0.0F)) {
            mul = 1.97215226E-31F;
            absxk = cfrom1;
          } else if (cto1 > absxk) {
            mul = 5.0706024E+30F;
            ctoc = cto1;
          } else {
            mul = ctoc / absxk;
            notdone = FALSE;
          }

          jrow = A->size[0] * A->size[1];
          A->size[0] = A->size[0];
          A->size[1] = A->size[1];
          emxEnsureCapacity((emxArray__common *)A, jrow, (int32_T)sizeof
                            (creal32_T));
          jcol = A->size[0];
          jcolp1 = A->size[1];
          jcol = jcol * jcolp1 - 1;
          for (jrow = 0; jrow <= jcol; jrow++) {
            A->data[jrow].re *= mul;
            A->data[jrow].im *= mul;
          }
        }
      }

      emxInit_creal32_T(&b_A, 2);
      jrow = b_A->size[0] * b_A->size[1];
      b_A->size[0] = A->size[0];
      b_A->size[1] = A->size[1];
      emxEnsureCapacity((emxArray__common *)b_A, jrow, (int32_T)sizeof(creal32_T));
      jcol = A->size[0] * A->size[1] - 1;
      for (jrow = 0; jrow <= jcol; jrow++) {
        b_A->data[jrow] = A->data[jrow];
      }

      emxInit_int32_T(&rscale, 1);
      emxInit_real32_T(&I, 2);
      eml_matlab_zggbal(b_A, &ilo, &ihi, rscale);
      b_n = b_A->size[0];
      jrow = I->size[0] * I->size[1];
      I->size[0] = b_n;
      I->size[1] = b_n;
      emxEnsureCapacity((emxArray__common *)I, jrow, (int32_T)sizeof(real32_T));
      jcol = b_n * b_n - 1;
      for (jrow = 0; jrow <= jcol; jrow++) {
        I->data[jrow] = 0.0F;
      }

      if (b_n > 0) {
        for (i = 0; i + 1 <= b_n; i++) {
          I->data[i + I->size[0] * i] = 1.0F;
        }
      }

      jrow = V->size[0] * V->size[1];
      V->size[0] = I->size[0];
      V->size[1] = I->size[1];
      emxEnsureCapacity((emxArray__common *)V, jrow, (int32_T)sizeof(creal32_T));
      jcol = I->size[0] * I->size[1] - 1;
      for (jrow = 0; jrow <= jcol; jrow++) {
        V->data[jrow].re = I->data[jrow];
        V->data[jrow].im = 0.0F;
      }

      emxFree_real32_T(&I);
      if ((b_n <= 1) || (ihi < ilo + 2)) {
      } else {
        jcol = ilo - 1;
        while (jcol + 1 < ihi - 1) {
          jcolp1 = jcol + 1;
          jrow = ihi - 1;
          while (jrow + 1 > jcolp1 + 1) {
            jrowm1 = jrow - 1;
            c_A = b_A->data[jrowm1 + b_A->size[0] * jcol];
            d_A = b_A->data[jrow + b_A->size[0] * jcol];
            eml_matlab_zlartg(c_A, d_A, &mul, &s, &tmp);
            b_A->data[jrowm1 + b_A->size[0] * jcol] = tmp;
            b_A->data[jrow + b_A->size[0] * jcol].re = 0.0F;
            b_A->data[jrow + b_A->size[0] * jcol].im = 0.0F;
            for (j = jcolp1; j + 1 <= ihi; j++) {
              tmp.re = mul * b_A->data[jrowm1 + b_A->size[0] * j].re;
              tmp.im = mul * b_A->data[jrowm1 + b_A->size[0] * j].im;
              y_re = s.re * b_A->data[jrow + b_A->size[0] * j].re - s.im *
                b_A->data[jrow + b_A->size[0] * j].im;
              y_im = s.re * b_A->data[jrow + b_A->size[0] * j].im + s.im *
                b_A->data[jrow + b_A->size[0] * j].re;
              absxk = b_A->data[jrowm1 + b_A->size[0] * j].re;
              ctoc = b_A->data[jrowm1 + b_A->size[0] * j].im;
              cfrom1 = b_A->data[jrowm1 + b_A->size[0] * j].im;
              cto1 = b_A->data[jrowm1 + b_A->size[0] * j].re;
              b_A->data[jrow + b_A->size[0] * j].re = mul * b_A->data[jrow +
                b_A->size[0] * j].re - (s.re * absxk + s.im * ctoc);
              b_A->data[jrow + b_A->size[0] * j].im = mul * b_A->data[jrow +
                b_A->size[0] * j].im - (s.re * cfrom1 - s.im * cto1);
              b_A->data[jrowm1 + b_A->size[0] * j].re = tmp.re + y_re;
              b_A->data[jrowm1 + b_A->size[0] * j].im = tmp.im + y_im;
            }

            s.re = -s.re;
            s.im = -s.im;
            for (i = ilo - 1; i + 1 <= ihi; i++) {
              tmp.re = mul * b_A->data[i + b_A->size[0] * jrow].re;
              tmp.im = mul * b_A->data[i + b_A->size[0] * jrow].im;
              y_re = s.re * b_A->data[i + b_A->size[0] * jrowm1].re - s.im *
                b_A->data[i + b_A->size[0] * jrowm1].im;
              y_im = s.re * b_A->data[i + b_A->size[0] * jrowm1].im + s.im *
                b_A->data[i + b_A->size[0] * jrowm1].re;
              absxk = b_A->data[i + b_A->size[0] * jrow].re;
              ctoc = b_A->data[i + b_A->size[0] * jrow].im;
              cfrom1 = b_A->data[i + b_A->size[0] * jrow].im;
              cto1 = b_A->data[i + b_A->size[0] * jrow].re;
              b_A->data[i + b_A->size[0] * jrowm1].re = mul * b_A->data[i +
                b_A->size[0] * jrowm1].re - (s.re * absxk + s.im * ctoc);
              b_A->data[i + b_A->size[0] * jrowm1].im = mul * b_A->data[i +
                b_A->size[0] * jrowm1].im - (s.re * cfrom1 - s.im * cto1);
              b_A->data[i + b_A->size[0] * jrow].re = tmp.re + y_re;
              b_A->data[i + b_A->size[0] * jrow].im = tmp.im + y_im;
            }

            for (i = 0; i + 1 <= b_n; i++) {
              tmp.re = mul * V->data[i + V->size[0] * jrow].re;
              tmp.im = mul * V->data[i + V->size[0] * jrow].im;
              y_re = s.re * V->data[i + V->size[0] * jrowm1].re - s.im * V->
                data[i + V->size[0] * jrowm1].im;
              y_im = s.re * V->data[i + V->size[0] * jrowm1].im + s.im * V->
                data[i + V->size[0] * jrowm1].re;
              absxk = V->data[i + V->size[0] * jrow].re;
              ctoc = V->data[i + V->size[0] * jrow].im;
              cfrom1 = V->data[i + V->size[0] * jrow].im;
              cto1 = V->data[i + V->size[0] * jrow].re;
              V->data[i + V->size[0] * jrowm1].re = mul * V->data[i + V->size[0]
                * jrowm1].re - (s.re * absxk + s.im * ctoc);
              V->data[i + V->size[0] * jrowm1].im = mul * V->data[i + V->size[0]
                * jrowm1].im - (s.re * cfrom1 - s.im * cto1);
              V->data[i + V->size[0] * jrow].re = tmp.re + y_re;
              V->data[i + V->size[0] * jrow].im = tmp.im + y_im;
            }

            jrow = jrowm1;
          }

          jcol = jcolp1;
        }
      }

      eml_matlab_zhgeqz(b_A, ilo, ihi, V, &b_info, alpha1, beta1);
      *info = b_info;
      if (b_info != 0.0) {
      } else {
        eml_matlab_ztgevc(b_A, V);
        b_n = V->size[0];
        jcolp1 = V->size[1];
        if (ilo > 1) {
          for (i = ilo - 2; i + 1 >= 1; i--) {
            if (rscale->data[i] != i + 1) {
              for (j = 0; j + 1 <= jcolp1; j++) {
                tmp = V->data[i + V->size[0] * j];
                V->data[i + V->size[0] * j] = V->data[(rscale->data[i] + V->
                  size[0] * j) - 1];
                V->data[(rscale->data[i] + V->size[0] * j) - 1] = tmp;
              }
            }
          }
        }

        if (ihi < b_n) {
          while (ihi + 1 <= b_n) {
            if (rscale->data[ihi] != ihi + 1) {
              for (j = 0; j + 1 <= jcolp1; j++) {
                tmp = V->data[ihi + V->size[0] * j];
                V->data[ihi + V->size[0] * j] = V->data[(rscale->data[ihi] +
                  V->size[0] * j) - 1];
                V->data[(rscale->data[ihi] + V->size[0] * j) - 1] = tmp;
              }
            }

            ihi++;
          }
        }

        for (jcol = 0; jcol <= n - 1; jcol++) {
          absxk = (real32_T)fabs(V->data[V->size[0] * jcol].re) + (real32_T)fabs
            (V->data[V->size[0] * jcol].im);
          if (n > 1) {
            for (jcolp1 = 1; jcolp1 - 1 <= n - 2; jcolp1++) {
              ctoc = (real32_T)fabs(V->data[jcolp1 + V->size[0] * jcol].re) +
                (real32_T)fabs(V->data[jcolp1 + V->size[0] * jcol].im);
              if (ctoc > absxk) {
                absxk = ctoc;
              }
            }
          }

          if (absxk >= 9.09494702E-13F) {
            absxk = 1.0F / absxk;
            for (jcolp1 = 0; jcolp1 <= n - 1; jcolp1++) {
              V->data[jcolp1 + V->size[0] * jcol].re *= absxk;
              V->data[jcolp1 + V->size[0] * jcol].im *= absxk;
            }
          }
        }

        if (ilascl) {
          notdone = TRUE;
          while (notdone) {
            cfrom1 = anrmto * 1.97215226E-31F;
            cto1 = anrm / 5.0706024E+30F;
            if ((cfrom1 > anrm) && (anrm != 0.0F)) {
              mul = 1.97215226E-31F;
              anrmto = cfrom1;
            } else if (cto1 > anrmto) {
              mul = 5.0706024E+30F;
              anrm = cto1;
            } else {
              mul = anrm / anrmto;
              notdone = FALSE;
            }

            jrow = alpha1->size[0];
            alpha1->size[0] = alpha1->size[0];
            emxEnsureCapacity((emxArray__common *)alpha1, jrow, (int32_T)sizeof
                              (creal32_T));
            jcol = alpha1->size[0] - 1;
            for (jrow = 0; jrow <= jcol; jrow++) {
              alpha1->data[jrow].re *= mul;
              alpha1->data[jrow].im *= mul;
            }
          }
        }
      }

      emxFree_int32_T(&rscale);
      emxFree_creal32_T(&b_A);
    }
  }
}

static void eml_matlab_zhgeqz(emxArray_creal32_T *A, int32_T ilo, int32_T ihi,
  emxArray_creal32_T *Z, real_T *info, emxArray_creal32_T *alpha1,
  emxArray_creal32_T *beta1)
{
  boolean_T b0;
  int32_T n;
  int32_T i;
  int32_T loop_ub;
  real32_T eshift_re;
  real32_T eshift_im;
  creal32_T ctemp;
  real32_T rho_re;
  real32_T rho_im;
  real32_T anorm;
  real32_T temp;
  real32_T b_atol;
  real32_T ascale;
  boolean_T failed;
  int32_T j;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  int32_T ifirst;
  int32_T istart;
  int32_T ilast;
  int32_T ilastm1;
  int32_T ifrstm;
  int32_T ilastm;
  int32_T iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int32_T jiter;
  int32_T exitg1;
  boolean_T exitg3;
  int32_T jm1;
  boolean_T ilazro;
  int32_T b_j;
  boolean_T b_guard1 = FALSE;
  creal32_T t1;
  creal32_T d;
  creal32_T sigma1;
  real32_T sigma2_re;
  real32_T sigma2_im;
  int32_T jp1;
  boolean_T exitg2;
  real32_T tempr;
  creal32_T fc0;
  if (!((Z->size[0] == 0) || (Z->size[1] == 0))) {
    b0 = TRUE;
  } else {
    b0 = FALSE;
  }

  n = A->size[0];
  i = alpha1->size[0];
  alpha1->size[0] = n;
  emxEnsureCapacity((emxArray__common *)alpha1, i, (int32_T)sizeof(creal32_T));
  loop_ub = n - 1;
  for (i = 0; i <= loop_ub; i++) {
    alpha1->data[i].re = 0.0F;
    alpha1->data[i].im = 0.0F;
  }

  i = beta1->size[0];
  beta1->size[0] = n;
  emxEnsureCapacity((emxArray__common *)beta1, i, (int32_T)sizeof(creal32_T));
  loop_ub = n - 1;
  for (i = 0; i <= loop_ub; i++) {
    beta1->data[i].re = 1.0F;
    beta1->data[i].im = 0.0F;
  }

  eshift_re = 0.0F;
  eshift_im = 0.0F;
  ctemp.re = 0.0F;
  ctemp.im = 0.0F;
  rho_re = 0.0F;
  rho_im = 0.0F;
  anorm = eml_matlab_zlanhs(A, ilo, ihi);
  temp = 1.1920929E-7F * anorm;
  b_atol = 1.17549435E-38F;
  if (temp > 1.17549435E-38F) {
    b_atol = temp;
  }

  temp = 1.17549435E-38F;
  if (anorm > 1.17549435E-38F) {
    temp = anorm;
  }

  ascale = 1.0F / temp;
  failed = TRUE;
  for (j = ihi; j + 1 <= n; j++) {
    alpha1->data[j] = A->data[j + A->size[0] * j];
  }

  guard1 = FALSE;
  guard2 = FALSE;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    if (b0) {
      ifrstm = 1;
      ilastm = n;
    } else {
      ifrstm = ilo;
      ilastm = ihi;
    }

    iiter = 0;
    goto60 = FALSE;
    goto70 = FALSE;
    goto90 = FALSE;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = TRUE;
        } else if ((real32_T)fabs(A->data[ilast + A->size[0] * ilastm1].re) +
                   (real32_T)fabs(A->data[ilast + A->size[0] * ilastm1].im) <=
                   b_atol) {
          A->data[ilast + A->size[0] * ilastm1].re = 0.0F;
          A->data[ilast + A->size[0] * ilastm1].im = 0.0F;
          goto60 = TRUE;
        } else {
          j = ilastm1;
          exitg3 = FALSE;
          while ((exitg3 == 0U) && (j + 1 >= ilo)) {
            jm1 = j - 1;
            if (j + 1 == ilo) {
              ilazro = TRUE;
            } else if ((real32_T)fabs(A->data[j + A->size[0] * jm1].re) +
                       (real32_T)fabs(A->data[j + A->size[0] * jm1].im) <=
                       b_atol) {
              A->data[j + A->size[0] * jm1].re = 0.0F;
              A->data[j + A->size[0] * jm1].im = 0.0F;
              ilazro = TRUE;
            } else {
              ilazro = FALSE;
            }

            if (ilazro) {
              ifirst = j + 1;
              goto70 = TRUE;
              exitg3 = TRUE;
            } else {
              j = jm1;
            }
          }
        }

        if (goto60 || goto70) {
          ilazro = TRUE;
        } else {
          ilazro = FALSE;
        }

        if (!ilazro) {
          jm1 = alpha1->size[0];
          i = alpha1->size[0];
          alpha1->size[0] = jm1;
          emxEnsureCapacity((emxArray__common *)alpha1, i, (int32_T)sizeof
                            (creal32_T));
          loop_ub = jm1 - 1;
          for (i = 0; i <= loop_ub; i++) {
            alpha1->data[i].re = ((real32_T)rtNaN);
            alpha1->data[i].im = 0.0F;
          }

          jm1 = beta1->size[0];
          i = beta1->size[0];
          beta1->size[0] = jm1;
          emxEnsureCapacity((emxArray__common *)beta1, i, (int32_T)sizeof
                            (creal32_T));
          loop_ub = jm1 - 1;
          for (i = 0; i <= loop_ub; i++) {
            beta1->data[i].re = ((real32_T)rtNaN);
            beta1->data[i].im = 0.0F;
          }

          if (b0) {
            i = Z->size[0] * Z->size[1];
            Z->size[0] = Z->size[0];
            Z->size[1] = Z->size[1];
            emxEnsureCapacity((emxArray__common *)Z, i, (int32_T)sizeof
                              (creal32_T));
            loop_ub = Z->size[1] - 1;
            for (i = 0; i <= loop_ub; i++) {
              jm1 = Z->size[0] - 1;
              for (b_j = 0; b_j <= jm1; b_j++) {
                Z->data[b_j + Z->size[0] * i].re = ((real32_T)rtNaN);
                Z->data[b_j + Z->size[0] * i].im = 0.0F;
              }
            }
          }

          *info = -1.0;
          exitg1 = 1;
        } else {
          b_guard1 = FALSE;
          if (goto60) {
            goto60 = FALSE;
            alpha1->data[ilast] = A->data[ilast + A->size[0] * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = FALSE;
              guard2 = TRUE;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0F;
              eshift_im = 0.0F;
              if (!b0) {
                ilastm = ilast + 1;
                if (ifrstm > ilast + 1) {
                  ifrstm = ilo;
                }
              }

              b_guard1 = TRUE;
            }
          } else {
            if (goto70) {
              goto70 = FALSE;
              iiter++;
              if (!b0) {
                ifrstm = ifirst;
              }

              if (iiter - div_s32_floor(iiter, 10) * 10 != 0) {
                t1.re = -(A->data[ilast + A->size[0] * ilast].re - A->
                          data[ilastm1 + A->size[0] * ilastm1].re);
                t1.im = -(A->data[ilast + A->size[0] * ilast].im - A->
                          data[ilastm1 + A->size[0] * ilastm1].im);
                t1 = eml_div(t1, 2.0);
                temp = A->data[ilastm1 + A->size[0] * ilast].re * A->data[ilast
                  + A->size[0] * ilastm1].re - A->data[ilastm1 + A->size[0] *
                  ilast].im * A->data[ilast + A->size[0] * ilastm1].im;
                anorm = A->data[ilastm1 + A->size[0] * ilast].re * A->data[ilast
                  + A->size[0] * ilastm1].im + A->data[ilastm1 + A->size[0] *
                  ilast].im * A->data[ilast + A->size[0] * ilastm1].re;
                d.re = (t1.re * t1.re - t1.im * t1.im) + temp;
                d.im = (t1.re * t1.im + t1.im * t1.re) + anorm;
                b_sqrt(&d);
                sigma1.re = A->data[ilastm1 + A->size[0] * ilastm1].re - (t1.re
                  - d.re);
                sigma1.im = A->data[ilastm1 + A->size[0] * ilastm1].im - (t1.im
                  - d.im);
                sigma2_re = A->data[ilastm1 + A->size[0] * ilastm1].re - (t1.re
                  + d.re);
                sigma2_im = A->data[ilastm1 + A->size[0] * ilastm1].im - (t1.im
                  + d.im);
                rho_re = sigma1.re - A->data[ilast + A->size[0] * ilast].re;
                rho_im = sigma1.im - A->data[ilast + A->size[0] * ilast].im;
                temp = sigma2_re - A->data[ilast + A->size[0] * ilast].re;
                anorm = sigma2_im - A->data[ilast + A->size[0] * ilast].im;
                if (rt_hypotf_snf((real32_T)fabs(rho_re), (real32_T)fabs(rho_im))
                    <= rt_hypotf_snf((real32_T)fabs(temp), (real32_T)fabs(anorm)))
                {
                  sigma2_re = sigma1.re;
                  sigma2_im = sigma1.im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += A->data[ilast + A->size[0] * ilastm1].re;
                eshift_im += A->data[ilast + A->size[0] * ilastm1].im;
                sigma2_re = eshift_re;
                sigma2_im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = FALSE;
              while ((exitg2 == 0U) && (j + 1 > ifirst)) {
                jm1 = j - 1;
                istart = j + 1;
                ctemp.re = A->data[j + A->size[0] * j].re - sigma2_re;
                ctemp.im = A->data[j + A->size[0] * j].im - sigma2_im;
                temp = ascale * ((real32_T)fabs(ctemp.re) + (real32_T)fabs
                                 (ctemp.im));
                anorm = ascale * ((real32_T)fabs(A->data[jp1 + A->size[0] * j].
                  re) + (real32_T)fabs(A->data[jp1 + A->size[0] * j].im));
                tempr = temp;
                if (anorm > temp) {
                  tempr = anorm;
                }

                if ((tempr < 1.0F) && (tempr != 0.0F)) {
                  temp /= tempr;
                  anorm /= tempr;
                }

                if (((real32_T)fabs(A->data[j + A->size[0] * jm1].re) +
                     (real32_T)fabs(A->data[j + A->size[0] * jm1].im)) * anorm <=
                    temp * b_atol) {
                  goto90 = TRUE;
                  exitg2 = TRUE;
                } else {
                  jp1 = j;
                  j = jm1;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = A->data[(ifirst + A->size[0] * (ifirst - 1)) - 1].
                    re - sigma2_re;
                  ctemp.im = A->data[(ifirst + A->size[0] * (ifirst - 1)) - 1].
                    im - sigma2_im;
                }

                goto90 = TRUE;
              }
            }

            if (goto90) {
              goto90 = FALSE;
              t1 = A->data[istart + A->size[0] * (istart - 1)];
              b_eml_matlab_zlartg(ctemp, t1, &sigma2_im, &sigma1);
              j = istart - 1;
              jm1 = istart - 2;
              while (j + 1 < ilast + 1) {
                jp1 = j + 1;
                if (j + 1 > istart) {
                  t1 = A->data[j + A->size[0] * jm1];
                  d = A->data[jp1 + A->size[0] * jm1];
                  eml_matlab_zlartg(t1, d, &sigma2_im, &sigma1, &fc0);
                  A->data[j + A->size[0] * jm1] = fc0;
                  A->data[jp1 + A->size[0] * jm1].re = 0.0F;
                  A->data[jp1 + A->size[0] * jm1].im = 0.0F;
                }

                for (b_j = j; b_j + 1 <= ilastm; b_j++) {
                  t1.re = sigma2_im * A->data[j + A->size[0] * b_j].re;
                  t1.im = sigma2_im * A->data[j + A->size[0] * b_j].im;
                  d.re = sigma1.re * A->data[jp1 + A->size[0] * b_j].re -
                    sigma1.im * A->data[jp1 + A->size[0] * b_j].im;
                  d.im = sigma1.re * A->data[jp1 + A->size[0] * b_j].im +
                    sigma1.im * A->data[jp1 + A->size[0] * b_j].re;
                  temp = A->data[j + A->size[0] * b_j].re;
                  anorm = A->data[j + A->size[0] * b_j].im;
                  tempr = A->data[j + A->size[0] * b_j].im;
                  sigma2_re = A->data[j + A->size[0] * b_j].re;
                  A->data[jp1 + A->size[0] * b_j].re = sigma2_im * A->data[jp1 +
                    A->size[0] * b_j].re - (sigma1.re * temp + sigma1.im * anorm);
                  A->data[jp1 + A->size[0] * b_j].im = sigma2_im * A->data[jp1 +
                    A->size[0] * b_j].im - (sigma1.re * tempr - sigma1.im *
                    sigma2_re);
                  A->data[j + A->size[0] * b_j].re = t1.re + d.re;
                  A->data[j + A->size[0] * b_j].im = t1.im + d.im;
                }

                sigma1.re = -sigma1.re;
                sigma1.im = -sigma1.im;
                loop_ub = jp1 + 2;
                if (ilast + 1 < loop_ub) {
                  loop_ub = ilast + 1;
                }

                for (i = ifrstm - 1; i + 1 <= loop_ub; i++) {
                  t1.re = sigma2_im * A->data[i + A->size[0] * jp1].re;
                  t1.im = sigma2_im * A->data[i + A->size[0] * jp1].im;
                  d.re = sigma1.re * A->data[i + A->size[0] * j].re - sigma1.im *
                    A->data[i + A->size[0] * j].im;
                  d.im = sigma1.re * A->data[i + A->size[0] * j].im + sigma1.im *
                    A->data[i + A->size[0] * j].re;
                  temp = A->data[i + A->size[0] * jp1].re;
                  anorm = A->data[i + A->size[0] * jp1].im;
                  tempr = A->data[i + A->size[0] * jp1].im;
                  sigma2_re = A->data[i + A->size[0] * jp1].re;
                  A->data[i + A->size[0] * j].re = sigma2_im * A->data[i +
                    A->size[0] * j].re - (sigma1.re * temp + sigma1.im * anorm);
                  A->data[i + A->size[0] * j].im = sigma2_im * A->data[i +
                    A->size[0] * j].im - (sigma1.re * tempr - sigma1.im *
                    sigma2_re);
                  A->data[i + A->size[0] * jp1].re = t1.re + d.re;
                  A->data[i + A->size[0] * jp1].im = t1.im + d.im;
                }

                if (b0) {
                  for (i = 0; i + 1 <= n; i++) {
                    t1.re = sigma2_im * Z->data[i + Z->size[0] * jp1].re;
                    t1.im = sigma2_im * Z->data[i + Z->size[0] * jp1].im;
                    d.re = sigma1.re * Z->data[i + Z->size[0] * j].re -
                      sigma1.im * Z->data[i + Z->size[0] * j].im;
                    d.im = sigma1.re * Z->data[i + Z->size[0] * j].im +
                      sigma1.im * Z->data[i + Z->size[0] * j].re;
                    anorm = Z->data[i + Z->size[0] * jp1].re;
                    temp = Z->data[i + Z->size[0] * jp1].im;
                    tempr = Z->data[i + Z->size[0] * jp1].im;
                    sigma2_re = Z->data[i + Z->size[0] * jp1].re;
                    Z->data[i + Z->size[0] * j].re = sigma2_im * Z->data[i +
                      Z->size[0] * j].re - (sigma1.re * anorm + sigma1.im * temp);
                    Z->data[i + Z->size[0] * j].im = sigma2_im * Z->data[i +
                      Z->size[0] * j].im - (sigma1.re * tempr - sigma1.im *
                      sigma2_re);
                    Z->data[i + Z->size[0] * jp1].re = t1.re + d.re;
                    Z->data[i + Z->size[0] * jp1].im = t1.im + d.im;
                  }
                }

                jm1 = j;
                j = jp1;
              }
            }

            b_guard1 = TRUE;
          }

          if (b_guard1 == TRUE) {
            jiter++;
          }
        }
      } else {
        guard2 = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0U);
  } else {
    guard1 = TRUE;
  }

  if (guard2 == TRUE) {
    if (failed) {
      *info = (real_T)(ilast + 1);
      for (jm1 = 0; jm1 + 1 <= ilast + 1; jm1++) {
        alpha1->data[jm1].re = ((real32_T)rtNaN);
        alpha1->data[jm1].im = 0.0F;
        beta1->data[jm1].re = ((real32_T)rtNaN);
        beta1->data[jm1].im = 0.0F;
      }

      if (b0) {
        i = Z->size[0] * Z->size[1];
        Z->size[0] = Z->size[0];
        Z->size[1] = Z->size[1];
        emxEnsureCapacity((emxArray__common *)Z, i, (int32_T)sizeof(creal32_T));
        loop_ub = Z->size[1] - 1;
        for (i = 0; i <= loop_ub; i++) {
          jm1 = Z->size[0] - 1;
          for (b_j = 0; b_j <= jm1; b_j++) {
            Z->data[b_j + Z->size[0] * i].re = ((real32_T)rtNaN);
            Z->data[b_j + Z->size[0] * i].im = 0.0F;
          }
        }
      }
    } else {
      guard1 = TRUE;
    }
  }

  if (guard1 == TRUE) {
    for (j = 0; j + 1 <= ilo - 1; j++) {
      alpha1->data[j] = A->data[j + A->size[0] * j];
    }

    *info = 0.0;
  }
}

static real32_T eml_matlab_zlanhs(const emxArray_creal32_T *A, int32_T ilo,
  int32_T ihi)
{
  real32_T f;
  real32_T scale;
  real32_T sumsq;
  boolean_T firstNonZero;
  int32_T j;
  int32_T c;
  int32_T i;
  real32_T reAij;
  real32_T imAij;
  real32_T temp2;
  f = 0.0F;
  if (ilo > ihi) {
  } else {
    scale = 0.0F;
    sumsq = 0.0F;
    firstNonZero = TRUE;
    for (j = ilo; j <= ihi; j++) {
      c = j + 1;
      if (ihi < c) {
        c = ihi;
      }

      for (i = ilo; i <= c; i++) {
        reAij = A->data[(i + A->size[0] * (j - 1)) - 1].re;
        imAij = A->data[(i + A->size[0] * (j - 1)) - 1].im;
        if (reAij != 0.0F) {
          reAij = (real32_T)fabs(reAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            scale = reAij;
            firstNonZero = FALSE;
          } else if (scale < reAij) {
            temp2 = scale / reAij;
            sumsq = 1.0F + sumsq * temp2 * temp2;
            scale = reAij;
          } else {
            temp2 = reAij / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (imAij != 0.0F) {
          reAij = (real32_T)fabs(imAij);
          if (firstNonZero) {
            sumsq = 1.0F;
            scale = reAij;
            firstNonZero = FALSE;
          } else if (scale < reAij) {
            temp2 = scale / reAij;
            sumsq = 1.0F + sumsq * temp2 * temp2;
            scale = reAij;
          } else {
            temp2 = reAij / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    f = scale * (real32_T)sqrt(sumsq);
  }

  return f;
}

static void eml_matlab_zlartg(const creal32_T f, const creal32_T g, real32_T *cs,
  creal32_T *sn, creal32_T *r)
{
  real32_T scale;
  real32_T f2s;
  real32_T g2;
  real32_T fs_re;
  real32_T fs_im;
  real32_T gs_re;
  real32_T gs_im;
  int32_T count;
  int32_T rescaledir;
  boolean_T guard1 = FALSE;
  real32_T g2s;
  scale = (real32_T)fabs(f.re);
  f2s = (real32_T)fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = (real32_T)fabs(g.re);
  g2 = (real32_T)fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = 0;
  rescaledir = 0;
  guard1 = FALSE;
  if (scale >= 5.49755814E+11F) {
    do {
      count++;
      fs_re *= 1.8189894E-12F;
      fs_im *= 1.8189894E-12F;
      gs_re *= 1.8189894E-12F;
      gs_im *= 1.8189894E-12F;
      scale *= 1.8189894E-12F;
    } while (!(scale < 5.49755814E+11F));

    rescaledir = 1;
    guard1 = TRUE;
  } else if (scale <= 1.8189894E-12F) {
    if ((g.re == 0.0F) && (g.im == 0.0F)) {
      *cs = 1.0F;
      sn->re = 0.0F;
      sn->im = 0.0F;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 5.49755814E+11F;
        fs_im *= 5.49755814E+11F;
        gs_re *= 5.49755814E+11F;
        gs_im *= 5.49755814E+11F;
        scale *= 5.49755814E+11F;
      } while (!(scale > 1.8189894E-12F));

      rescaledir = -1;
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0F > g2) {
      f2s = 1.0F;
    }

    if (scale <= f2s * 1.97215226E-31F) {
      if ((f.re == 0.0F) && (f.im == 0.0F)) {
        *cs = 0.0F;
        r->re = rt_hypotf_snf((real32_T)fabs(g.re), (real32_T)fabs(g.im));
        r->im = 0.0F;
        f2s = rt_hypotf_snf((real32_T)fabs(gs_re), (real32_T)fabs(gs_im));
        sn->re = gs_re / f2s;
        sn->im = -gs_im / f2s;
      } else {
        g2s = (real32_T)sqrt(g2);
        *cs = rt_hypotf_snf((real32_T)fabs(fs_re), (real32_T)fabs(fs_im)) / g2s;
        f2s = (real32_T)fabs(f.re);
        g2 = (real32_T)fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0F) {
          f2s = rt_hypotf_snf((real32_T)fabs(f.re), (real32_T)fabs(f.im));
          fs_re = f.re / f2s;
          fs_im = f.im / f2s;
        } else {
          g2 = 5.49755814E+11F * f.re;
          scale = 5.49755814E+11F * f.im;
          f2s = rt_hypotf_snf((real32_T)fabs(g2), (real32_T)fabs(scale));
          fs_re = g2 / f2s;
          fs_im = scale / f2s;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = (real32_T)sqrt(1.0F + g2 / scale);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0F / f2s;
      f2s = scale + g2;
      g2 = r->re / f2s;
      f2s = r->im / f2s;
      sn->re = g2 * gs_re - f2s * -gs_im;
      sn->im = g2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 1; rescaledir <= count; rescaledir++) {
          r->re *= 5.49755814E+11F;
          r->im *= 5.49755814E+11F;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 1; rescaledir <= count; rescaledir++) {
            r->re *= 1.8189894E-12F;
            r->im *= 1.8189894E-12F;
          }
        }
      }
    }
  }
}

static void eml_matlab_ztgevc(const emxArray_creal32_T *A, emxArray_creal32_T *V)
{
  emxArray_creal32_T *b_V;
  emxArray_creal32_T *work1;
  int32_T je;
  int32_T loop_ub;
  emxArray_creal32_T *work2;
  emxArray_real32_T *rworka;
  real32_T SMALL;
  real32_T BIG;
  real32_T BIGNUM;
  real32_T anorm;
  int32_T j;
  real32_T y;
  real32_T scale;
  real32_T ascale;
  int32_T b_je;
  real32_T temp;
  real32_T temp_im;
  real32_T salpha_re;
  real32_T salpha_im;
  real32_T acoeff;
  boolean_T b1;
  boolean_T b2;
  real_T x;
  real32_T acoefa;
  int32_T jr;
  real32_T dmin;
  real32_T d_re;
  real32_T d_im;
  real32_T A_im;
  emxInit_creal32_T(&b_V, 2);
  b_emxInit_creal32_T(&work1, 1);
  je = work1->size[0];
  work1->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)work1, je, (int32_T)sizeof(creal32_T));
  loop_ub = A->size[0] - 1;
  for (je = 0; je <= loop_ub; je++) {
    work1->data[je].re = 0.0F;
    work1->data[je].im = 0.0F;
  }

  b_emxInit_creal32_T(&work2, 1);
  je = work2->size[0];
  work2->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)work2, je, (int32_T)sizeof(creal32_T));
  loop_ub = A->size[0] - 1;
  for (je = 0; je <= loop_ub; je++) {
    work2->data[je].re = 0.0F;
    work2->data[je].im = 0.0F;
  }

  b_emxInit_real32_T(&rworka, 1);
  SMALL = 1.17549435E-38F * (real32_T)A->size[0] / 1.1920929E-7F;
  BIG = 1.0F / SMALL;
  BIGNUM = 1.0F / (1.17549435E-38F * (real32_T)A->size[0]);
  je = rworka->size[0];
  rworka->size[0] = A->size[0];
  emxEnsureCapacity((emxArray__common *)rworka, je, (int32_T)sizeof(real32_T));
  loop_ub = A->size[0] - 1;
  for (je = 0; je <= loop_ub; je++) {
    rworka->data[je] = 0.0F;
  }

  anorm = (real32_T)fabs(A->data[0].re) + (real32_T)fabs(A->data[0].im);
  for (j = 1; j - 1 <= A->size[0] - 2; j++) {
    for (je = 0; je <= j - 1; je++) {
      rworka->data[j] += (real32_T)fabs(A->data[je + A->size[0] * j].re) +
        (real32_T)fabs(A->data[je + A->size[0] * j].im);
    }

    y = rworka->data[j] + ((real32_T)fabs(A->data[j + A->size[0] * j].re) +
      (real32_T)fabs(A->data[j + A->size[0] * j].im));
    if (y > anorm) {
      anorm = y;
    }
  }

  scale = anorm;
  if (1.17549435E-38F > anorm) {
    scale = 1.17549435E-38F;
  }

  ascale = 1.0F / scale;
  for (je = 1; je - 1 <= A->size[0] - 1; je++) {
    b_je = A->size[0] - je;
    y = ((real32_T)fabs(A->data[b_je + A->size[0] * b_je].re) + (real32_T)fabs
         (A->data[b_je + A->size[0] * b_je].im)) * ascale;
    if (1.0F > y) {
      y = 1.0F;
    }

    temp = 1.0F / y;
    scale = temp * A->data[b_je + A->size[0] * b_je].re;
    temp_im = temp * A->data[b_je + A->size[0] * b_je].im;
    salpha_re = ascale * scale;
    salpha_im = ascale * temp_im;
    acoeff = temp * ascale;
    if (((real32_T)fabs(temp) >= 1.17549435E-38F) && ((real32_T)fabs(acoeff) <
         SMALL)) {
      b1 = TRUE;
    } else {
      b1 = FALSE;
    }

    if (((real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im) >=
         1.17549435E-38F) && ((real32_T)fabs(salpha_re) + (real32_T)fabs
         (salpha_im) < SMALL)) {
      b2 = TRUE;
    } else {
      b2 = FALSE;
    }

    scale = 1.0F;
    if (b1) {
      scale = anorm;
      if (BIG < anorm) {
        scale = BIG;
      }

      scale *= SMALL / (real32_T)fabs(temp);
    }

    if (b2) {
      x = 1.0;
      if (BIG < 1.0F) {
        x = BIG;
      }

      y = SMALL / ((real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im)) *
        (real32_T)x;
      if (y > scale) {
        scale = y;
      }
    }

    if (b1 || b2) {
      y = (real32_T)fabs(acoeff);
      temp_im = (real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im);
      if (1.0F > y) {
        y = 1.0F;
      }

      if (temp_im > y) {
        y = temp_im;
      }

      y = 1.0F / (1.17549435E-38F * y);
      if (y < scale) {
        scale = y;
      }

      if (b1) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      if (b2) {
        salpha_re *= scale;
        salpha_im *= scale;
      } else {
        salpha_re *= scale;
        salpha_im *= scale;
      }
    }

    acoefa = (real32_T)fabs(acoeff);
    for (jr = 0; jr <= A->size[0] - 1; jr++) {
      work1->data[jr].re = 0.0F;
      work1->data[jr].im = 0.0F;
    }

    work1->data[b_je].re = 1.0F;
    work1->data[b_je].im = 0.0F;
    dmin = 1.1920929E-7F * acoefa * anorm;
    y = 1.1920929E-7F * ((real32_T)fabs(salpha_re) + (real32_T)fabs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }

    if (1.17549435E-38F > dmin) {
      dmin = 1.17549435E-38F;
    }

    for (jr = 0; jr <= b_je - 1; jr++) {
      work1->data[jr].re = acoeff * A->data[jr + A->size[0] * b_je].re;
      work1->data[jr].im = acoeff * A->data[jr + A->size[0] * b_je].im;
    }

    work1->data[b_je].re = 1.0F;
    work1->data[b_je].im = 0.0F;
    for (j = 1; j - 1 <= b_je - 1; j++) {
      loop_ub = b_je - j;
      d_re = acoeff * A->data[loop_ub + A->size[0] * loop_ub].re - salpha_re;
      d_im = acoeff * A->data[loop_ub + A->size[0] * loop_ub].im - salpha_im;
      if ((real32_T)fabs(d_re) + (real32_T)fabs(d_im) <= dmin) {
        d_re = dmin;
        d_im = 0.0F;
      }

      if (((real32_T)fabs(d_re) + (real32_T)fabs(d_im) < 1.0F) && ((real32_T)
           fabs(work1->data[loop_ub].re) + (real32_T)fabs(work1->data[loop_ub].
            im) >= BIGNUM * ((real32_T)fabs(d_re) + (real32_T)fabs(d_im)))) {
        temp = 1.0F / ((real32_T)fabs(work1->data[loop_ub].re) + (real32_T)fabs
                       (work1->data[loop_ub].im));
        for (jr = 0; jr <= b_je; jr++) {
          work1->data[jr].re *= temp;
          work1->data[jr].im *= temp;
        }
      }

      temp = -work1->data[loop_ub].re;
      A_im = -work1->data[loop_ub].im;
      if (d_im == 0.0F) {
        if (A_im == 0.0F) {
          work1->data[loop_ub].re = temp / d_re;
          work1->data[loop_ub].im = 0.0F;
        } else if (temp == 0.0F) {
          work1->data[loop_ub].re = 0.0F;
          work1->data[loop_ub].im = A_im / d_re;
        } else {
          work1->data[loop_ub].re = temp / d_re;
          work1->data[loop_ub].im = A_im / d_re;
        }
      } else if (d_re == 0.0F) {
        if (temp == 0.0F) {
          work1->data[loop_ub].re = A_im / d_im;
          work1->data[loop_ub].im = 0.0F;
        } else if (A_im == 0.0F) {
          work1->data[loop_ub].re = 0.0F;
          work1->data[loop_ub].im = -(temp / d_im);
        } else {
          work1->data[loop_ub].re = A_im / d_im;
          work1->data[loop_ub].im = -(temp / d_im);
        }
      } else {
        y = (real32_T)fabs(d_re);
        temp_im = (real32_T)fabs(d_im);
        if (y > temp_im) {
          scale = d_im / d_re;
          temp_im = d_re + scale * d_im;
          work1->data[loop_ub].re = (temp + scale * A_im) / temp_im;
          work1->data[loop_ub].im = (A_im - scale * temp) / temp_im;
        } else if (temp_im == y) {
          scale = d_re > 0.0F ? 0.5F : -0.5F;
          temp_im = d_im > 0.0F ? 0.5F : -0.5F;
          work1->data[loop_ub].re = (temp * scale + A_im * temp_im) / y;
          work1->data[loop_ub].im = (A_im * scale - temp * temp_im) / y;
        } else {
          scale = d_re / d_im;
          temp_im = d_im + scale * d_re;
          work1->data[loop_ub].re = (scale * temp + A_im) / temp_im;
          work1->data[loop_ub].im = (scale * A_im - temp) / temp_im;
        }
      }

      if (loop_ub + 1 > 1) {
        if ((real32_T)fabs(work1->data[loop_ub].re) + (real32_T)fabs(work1->
             data[loop_ub].im) > 1.0F) {
          temp = 1.0F / ((real32_T)fabs(work1->data[loop_ub].re) + (real32_T)
                         fabs(work1->data[loop_ub].im));
          if (acoefa * rworka->data[loop_ub] >= BIGNUM * temp) {
            for (jr = 0; jr <= b_je; jr++) {
              work1->data[jr].re *= temp;
              work1->data[jr].im *= temp;
            }
          }
        }

        d_re = acoeff * work1->data[loop_ub].re;
        d_im = acoeff * work1->data[loop_ub].im;
        for (jr = 0; jr <= loop_ub - 1; jr++) {
          temp_im = d_re * A->data[jr + A->size[0] * loop_ub].re - d_im *
            A->data[jr + A->size[0] * loop_ub].im;
          scale = d_re * A->data[jr + A->size[0] * loop_ub].im + d_im * A->
            data[jr + A->size[0] * loop_ub].re;
          work1->data[jr].re += temp_im;
          work1->data[jr].im += scale;
        }
      }
    }

    for (jr = 0; jr <= A->size[0] - 1; jr++) {
      work2->data[jr].re = 0.0F;
      work2->data[jr].im = 0.0F;
    }

    for (loop_ub = 0; loop_ub <= b_je; loop_ub++) {
      for (jr = 0; jr <= A->size[0] - 1; jr++) {
        scale = V->data[jr + V->size[0] * loop_ub].re * work1->data[loop_ub].re
          - V->data[jr + V->size[0] * loop_ub].im * work1->data[loop_ub].im;
        temp_im = V->data[jr + V->size[0] * loop_ub].re * work1->data[loop_ub].
          im + V->data[jr + V->size[0] * loop_ub].im * work1->data[loop_ub].re;
        work2->data[jr].re += scale;
        work2->data[jr].im += temp_im;
      }
    }

    scale = (real32_T)fabs(work2->data[0].re) + (real32_T)fabs(work2->data[0].im);
    if (A->size[0] > 1) {
      for (jr = 1; jr - 1 <= A->size[0] - 2; jr++) {
        y = (real32_T)fabs(work2->data[jr].re) + (real32_T)fabs(work2->data[jr].
          im);
        if (y > scale) {
          scale = y;
        }
      }
    }

    if (scale > 1.17549435E-38F) {
      temp = 1.0F / scale;
      for (jr = 0; jr <= A->size[0] - 1; jr++) {
        V->data[jr + V->size[0] * b_je].re = temp * work2->data[jr].re;
        V->data[jr + V->size[0] * b_je].im = temp * work2->data[jr].im;
      }
    } else {
      for (jr = 0; jr <= A->size[0] - 1; jr++) {
        V->data[jr + V->size[0] * b_je].re = 0.0F;
        V->data[jr + V->size[0] * b_je].im = 0.0F;
      }
    }
  }

  emxFree_real32_T(&rworka);
  emxFree_creal32_T(&work2);
  emxFree_creal32_T(&work1);
  emxFree_creal32_T(&b_V);
}

void geteigen(const emxArray_real32_T *A, emxArray_real32_T *eigenvector,
              emxArray_real32_T *eigenvalue)
{
  emxArray_creal32_T *b_A;
  int32_T i0;
  int32_T kend;
  emxArray_creal32_T *c_A;
  emxArray_creal32_T *alpha1;
  emxArray_creal32_T *beta1;
  emxArray_creal32_T *V;
  real_T info;
  int32_T n;
  int32_T nm1;
  int32_T lastcol;
  int32_T coltop;
  real32_T colnorm;
  real32_T scale;
  int32_T nv;
  real32_T t;
  real32_T absxk;
  real32_T alpha1_re;
  real32_T alpha1_im;
  emxInit_creal32_T(&b_A, 2);
  i0 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity((emxArray__common *)b_A, i0, (int32_T)sizeof(creal32_T));
  kend = A->size[0] * A->size[1] - 1;
  for (i0 = 0; i0 <= kend; i0++) {
    b_A->data[i0].re = A->data[i0];
    b_A->data[i0].im = 0.0F;
  }

  emxInit_creal32_T(&c_A, 2);
  i0 = c_A->size[0] * c_A->size[1];
  c_A->size[0] = b_A->size[0];
  c_A->size[1] = b_A->size[1];
  emxEnsureCapacity((emxArray__common *)c_A, i0, (int32_T)sizeof(creal32_T));
  kend = b_A->size[0] * b_A->size[1] - 1;
  for (i0 = 0; i0 <= kend; i0++) {
    c_A->data[i0] = b_A->data[i0];
  }

  b_emxInit_creal32_T(&alpha1, 1);
  b_emxInit_creal32_T(&beta1, 1);
  emxInit_creal32_T(&V, 2);
  eml_matlab_zggev(c_A, &info, alpha1, beta1, V);
  n = b_A->size[0];
  emxFree_creal32_T(&c_A);
  if (n > 0) {
    nm1 = n - 1;
    lastcol = nm1 * n + 1;
    for (coltop = 0; coltop + 1 <= lastcol; coltop += n) {
      colnorm = 0.0F;
      if (n == 1) {
        colnorm = rt_hypotf_snf((real32_T)fabs(V->data[coltop].re), (real32_T)
          fabs(V->data[coltop].im));
      } else {
        scale = 1.17549435E-38F;
        kend = coltop + n;
        for (nv = coltop; nv + 1 <= kend; nv++) {
          t = V->data[nv].re;
          absxk = (real32_T)fabs(t);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = 1.0F + colnorm * t * t;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }

          t = V->data[nv].im;
          absxk = (real32_T)fabs(t);
          if (absxk > scale) {
            t = scale / absxk;
            colnorm = 1.0F + colnorm * t * t;
            scale = absxk;
          } else {
            t = absxk / scale;
            colnorm += t * t;
          }
        }

        colnorm = scale * (real32_T)sqrt(colnorm);
      }

      i0 = (coltop + nm1) + 1;
      for (kend = coltop; kend + 1 <= i0; kend++) {
        t = V->data[kend].re;
        absxk = V->data[kend].im;
        if (absxk == 0.0F) {
          V->data[kend].re = t / colnorm;
          V->data[kend].im = 0.0F;
        } else if (t == 0.0F) {
          V->data[kend].re = 0.0F;
          V->data[kend].im = absxk / colnorm;
        } else {
          V->data[kend].re = t / colnorm;
          V->data[kend].im = absxk / colnorm;
        }
      }
    }
  }

  i0 = alpha1->size[0];
  alpha1->size[0] = alpha1->size[0];
  emxEnsureCapacity((emxArray__common *)alpha1, i0, (int32_T)sizeof(creal32_T));
  kend = alpha1->size[0] - 1;
  for (i0 = 0; i0 <= kend; i0++) {
    alpha1_re = alpha1->data[i0].re;
    alpha1_im = alpha1->data[i0].im;
    t = beta1->data[i0].re;
    scale = beta1->data[i0].im;
    if (scale == 0.0F) {
      if (alpha1_im == 0.0F) {
        alpha1->data[i0].re = alpha1_re / t;
        alpha1->data[i0].im = 0.0F;
      } else if (alpha1_re == 0.0F) {
        alpha1->data[i0].re = 0.0F;
        alpha1->data[i0].im = alpha1_im / t;
      } else {
        alpha1->data[i0].re = alpha1_re / t;
        alpha1->data[i0].im = alpha1_im / t;
      }
    } else if (t == 0.0F) {
      if (alpha1_re == 0.0F) {
        alpha1->data[i0].re = alpha1_im / scale;
        alpha1->data[i0].im = 0.0F;
      } else if (alpha1_im == 0.0F) {
        alpha1->data[i0].re = 0.0F;
        alpha1->data[i0].im = -(alpha1_re / scale);
      } else {
        alpha1->data[i0].re = alpha1_im / scale;
        alpha1->data[i0].im = -(alpha1_re / scale);
      }
    } else {
      colnorm = (real32_T)fabs(t);
      absxk = (real32_T)fabs(scale);
      if (colnorm > absxk) {
        absxk = scale / t;
        t += absxk * scale;
        alpha1->data[i0].re = (alpha1_re + absxk * alpha1_im) / t;
        alpha1->data[i0].im = (alpha1_im - absxk * alpha1_re) / t;
      } else if (absxk == colnorm) {
        t = t > 0.0F ? 0.5F : -0.5F;
        absxk = scale > 0.0F ? 0.5F : -0.5F;
        alpha1->data[i0].re = (alpha1_re * t + alpha1_im * absxk) / colnorm;
        alpha1->data[i0].im = (alpha1_im * t - alpha1_re * absxk) / colnorm;
      } else {
        absxk = t / scale;
        t = scale + absxk * t;
        alpha1->data[i0].re = (absxk * alpha1_re + alpha1_im) / t;
        alpha1->data[i0].im = (absxk * alpha1_im - alpha1_re) / t;
      }
    }
  }

  emxFree_creal32_T(&beta1);
  nv = alpha1->size[0];
  i0 = b_A->size[0] * b_A->size[1];
  b_A->size[0] = nv;
  emxEnsureCapacity((emxArray__common *)b_A, i0, (int32_T)sizeof(creal32_T));
  i0 = b_A->size[0] * b_A->size[1];
  b_A->size[1] = nv;
  emxEnsureCapacity((emxArray__common *)b_A, i0, (int32_T)sizeof(creal32_T));
  kend = nv * nv - 1;
  for (i0 = 0; i0 <= kend; i0++) {
    b_A->data[i0].re = 0.0F;
    b_A->data[i0].im = 0.0F;
  }

  for (kend = 0; kend + 1 <= nv; kend++) {
    b_A->data[kend + b_A->size[0] * kend] = alpha1->data[kend];
  }

  emxFree_creal32_T(&alpha1);
  i0 = eigenvector->size[0] * eigenvector->size[1];
  eigenvector->size[0] = V->size[0];
  eigenvector->size[1] = V->size[1];
  emxEnsureCapacity((emxArray__common *)eigenvector, i0, (int32_T)sizeof
                    (real32_T));
  kend = V->size[0] * V->size[1] - 1;
  for (i0 = 0; i0 <= kend; i0++) {
    eigenvector->data[i0] = V->data[i0].re;
  }

  emxFree_creal32_T(&V);
  i0 = eigenvalue->size[0] * eigenvalue->size[1];
  eigenvalue->size[0] = b_A->size[0];
  eigenvalue->size[1] = b_A->size[1];
  emxEnsureCapacity((emxArray__common *)eigenvalue, i0, (int32_T)sizeof(real32_T));
  kend = b_A->size[0] * b_A->size[1] - 1;
  for (i0 = 0; i0 <= kend; i0++) {
    eigenvalue->data[i0] = b_A->data[i0].re;
  }

  emxFree_creal32_T(&b_A);
}

/* End of code generation (geteigen.c) */
