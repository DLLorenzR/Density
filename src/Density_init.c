/* Save as Density_init.c */
#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
 Check these declarations against the C/Fortran source code.
 */

/* .Fortran calls */
extern void F77_NAME(dsden)(double *, long *, double *, double *, double *, long *);


static const R_FortranMethodDef FortranEntries[] = {
  {"dsden", (DL_FUNC) &F77_NAME(dsden), 6},
  {NULL, NULL, 0}
};

void R_init_Density(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
}