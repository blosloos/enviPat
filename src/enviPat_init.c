#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP iso_centroid_Call(SEXP, SEXP, SEXP);
extern SEXP iso_pattern(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP iso_pattern_4(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP iso_ppm_Call(SEXP, SEXP, SEXP);
extern SEXP iso_profile_with_trace_Call(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"iso_centroid_Call",           (DL_FUNC) &iso_centroid_Call,           3},
    {"iso_pattern",                 (DL_FUNC) &iso_pattern,                 9},
    {"iso_pattern_4",               (DL_FUNC) &iso_pattern_4,               9},
    {"iso_ppm_Call",                (DL_FUNC) &iso_ppm_Call,                3},
    {"iso_profile_with_trace_Call", (DL_FUNC) &iso_profile_with_trace_Call, 7},
    {NULL, NULL, 0}
};

void R_init_enviPat(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
