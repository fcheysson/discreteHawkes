#include <Rinternals.h>

/* Error messages */
#define R_MSG_NA        _("NaNs produced")

/* Interfaces to routines from package expint */
double(*pkg_expint_E1)(double,int);
