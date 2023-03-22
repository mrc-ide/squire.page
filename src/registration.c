#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void nimue_booster_diff_rhs_dde(void *);
extern void nimue_booster_initmod_desolve(void *);
extern void nimue_booster_output_dde(void *);
extern void nimue_booster_rhs_dde(void *);
extern void nimue_booster_rhs_desolve(void *);

/* .Call calls */
extern SEXP nimue_booster_contents(SEXP);
extern SEXP nimue_booster_create(SEXP);
extern SEXP nimue_booster_diff_contents(SEXP);
extern SEXP nimue_booster_diff_create(SEXP);
extern SEXP nimue_booster_diff_initial_conditions(SEXP, SEXP);
extern SEXP nimue_booster_diff_metadata(SEXP);
extern SEXP nimue_booster_diff_rhs_r(SEXP, SEXP, SEXP);
extern SEXP nimue_booster_diff_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP nimue_booster_diff_set_user(SEXP, SEXP);
extern SEXP nimue_booster_initial_conditions(SEXP, SEXP);
extern SEXP nimue_booster_metadata(SEXP);
extern SEXP nimue_booster_rhs_r(SEXP, SEXP, SEXP);
extern SEXP nimue_booster_set_initial(SEXP, SEXP, SEXP, SEXP);
extern SEXP nimue_booster_set_user(SEXP, SEXP);

static const R_CMethodDef CEntries[] = {
    {"nimue_booster_diff_rhs_dde",    (DL_FUNC) &nimue_booster_diff_rhs_dde,    1},
    {"nimue_booster_initmod_desolve", (DL_FUNC) &nimue_booster_initmod_desolve, 1},
    {"nimue_booster_output_dde",      (DL_FUNC) &nimue_booster_output_dde,      1},
    {"nimue_booster_rhs_dde",         (DL_FUNC) &nimue_booster_rhs_dde,         1},
    {"nimue_booster_rhs_desolve",     (DL_FUNC) &nimue_booster_rhs_desolve,     1},
    {NULL, NULL, 0}
};

static const R_CallMethodDef CallEntries[] = {
    {"nimue_booster_contents",                (DL_FUNC) &nimue_booster_contents,                1},
    {"nimue_booster_create",                  (DL_FUNC) &nimue_booster_create,                  1},
    {"nimue_booster_diff_contents",           (DL_FUNC) &nimue_booster_diff_contents,           1},
    {"nimue_booster_diff_create",             (DL_FUNC) &nimue_booster_diff_create,             1},
    {"nimue_booster_diff_initial_conditions", (DL_FUNC) &nimue_booster_diff_initial_conditions, 2},
    {"nimue_booster_diff_metadata",           (DL_FUNC) &nimue_booster_diff_metadata,           1},
    {"nimue_booster_diff_rhs_r",              (DL_FUNC) &nimue_booster_diff_rhs_r,              3},
    {"nimue_booster_diff_set_initial",        (DL_FUNC) &nimue_booster_diff_set_initial,        4},
    {"nimue_booster_diff_set_user",           (DL_FUNC) &nimue_booster_diff_set_user,           2},
    {"nimue_booster_initial_conditions",      (DL_FUNC) &nimue_booster_initial_conditions,      2},
    {"nimue_booster_metadata",                (DL_FUNC) &nimue_booster_metadata,                1},
    {"nimue_booster_rhs_r",                   (DL_FUNC) &nimue_booster_rhs_r,                   3},
    {"nimue_booster_set_initial",             (DL_FUNC) &nimue_booster_set_initial,             4},
    {"nimue_booster_set_user",                (DL_FUNC) &nimue_booster_set_user,                2},
    {NULL, NULL, 0}
};

void R_init_squire_page(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
