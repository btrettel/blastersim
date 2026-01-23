# `-Werror=uninitialized` seems to require `-fcheck=all` (not used in the release build), or else there's a false-positive compile-time error, at least for gfortran 13. Check if this is true for later versions. `-Wno-uninitialized` eliminates this false-positive and should be removed in the future if gfortran eliminates this problem.

FFLAGS   = -Wall -Wextra -Werror -pedantic-errors -Wno-maybe-uninitialized -Wno-do-subscript -std=f2018 -Wconversion -Wconversion-extra -fimplicit-none -fno-unsafe-math-optimizations -finit-real=snan -finit-integer=-2147483647 -finit-logical=true -finit-derived -Wimplicit-interface -Wunused -ffree-line-length-132
DFLAGS   = -Og -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal -fsanitize=leak --coverage
RFLAGS   = -O2 -Wno-uninitialized -fopt-info-missed=$(MISSED) -flto
AFLAGS   = 
NFLAGS   = -march=native
SFLAGS   = -static
OMPFLAGS = -fopenmp

# <https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>

# TODO: `-funroll`
# Removed: `-fmax-errors=1`

# static linking:
# `-static`
# <https://gcc.gnu.org/onlinedocs/gfortran/Link-Options.html>
# <https://fortran-lang.discourse.group/t/problems-with-build-of-static-exe/8228/5>
# <https://fortran-lang.discourse.group/t/makefile-errors-using-gfortran-static-option/3491>

# `-Wno-do-subscript`
# <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=97320>
