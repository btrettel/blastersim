# `-check uninit` has false positives with namelists for the moment.
# <https://community.intel.com/t5/Intel-Fortran-Compiler/Known-bug-with-check-all-or-check-uninit-in-ifx-2024-0-0-for/m-p/1545825>

FFLAGS   = -warn errors -warn all -diag-error=remark,warn,error -fltconsistency -stand f18 -init=snan,arrays
DFLAGS   = -diag-disable=10440 -O0 -g -traceback -debug full -check all,nouninit -fpe0 -standard-semantics
RFLAGS   = -O2 -qopt-report -flto
AFLAGS   = 
NFLAGS   = -xHost
SFLAGS   = -static
OMPFLAGS = -fiopenmp

# TODO: Ask question about `do concurrent` on Intel Fortran forum. What's the right compiler flag? `-parallel` isn't for ifx contrary to what was said there. Update this afterward.

# Removed: `-diag-error-limit=1`

# `-xHost`: <https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/xhost-qxhost.html>

# static linking: <https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/static-002.html>
# `-static`

# `-standard-semantics` is only used for BUILD=debug due to a potential performance reduction.
# <https://github.com/fortran-lang/fpm/issues/868#issuecomment-1720960031>

# `-diag-disable=10440`: This is to avoid a warning about debug mode disabling optimizations. This warning might only appear on older versions of ifx.

# `-nostandard-realloc-lhs`: This could be useful to catch undesired reallocation on assignment. Memory allocation could be slow.
# <https://www.intel.com/content/www/us/en/docs/fortran-compiler/developer-guide-reference/2024-2/standard-realloc-lhs.html>
