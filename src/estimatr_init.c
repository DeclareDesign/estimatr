#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _estimatr_lm_robust_helper(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _estimatr_mat_sqrt_inv(SEXP);
extern SEXP _estimatr_mult_diag(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_estimatr_lm_robust_helper", (DL_FUNC) &_estimatr_lm_robust_helper, 6},
  {"_estimatr_mat_sqrt_inv",     (DL_FUNC) &_estimatr_mat_sqrt_inv,     1},
  {"_estimatr_mult_diag",        (DL_FUNC) &_estimatr_mult_diag,        2},
  {NULL, NULL, 0}
};

void R_init_estimatr(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, TRUE);
}
