# `-Werror=uninitialized` seems to require `-fcheck=all` (not used in the release build), or else there's a false-positive compile-time error, at least for gfortran 13. Check if this is true for later versions. `-Wno-uninitialized` eliminates this false-positive and should be removed in the future if gfortran eliminates this problem.
# `-fsanitize=leak` seems to continue to have false positives as of gfortran 13.
# `-Wdo-subscript` is misleadingly named and warns for things that are not problems. <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=90237>
# I think I had some false positives before that made me add `-Wno-maybe-uninitialized`, but it seems unnecessary now.

FFLAGS   = -Wall -Wextra -Werror -pedantic-errors -Wno-do-subscript -std=f2018 -Wconversion -Wconversion-extra -fimplicit-none -fno-unsafe-math-optimizations -finit-real=snan -finit-integer=2147483647 -finit-logical=true -finit-derived -Wimplicit-interface -Wunused -Wcharacter-truncation -Winteger-division -ffree-line-length-132
DFLAGS   = -Og -g -fbacktrace -fcheck=all -ffpe-trap=invalid,zero,overflow,underflow,denormal --coverage
RFLAGS   = -O2 -Wno-uninitialized -fopt-info-missed=$(MISSED) -flto
AFLAGS   = 
NFLAGS   = -march=native
SFLAGS   = -static
OMPFLAGS = -fopenmp -Wdeprecated-openmp

# Making `-finit-integer=2147483647` is better than making it a large negative number.
# Say you have a loop: `do i = 1, n` and you forget to initialize `n`.
# If `n` is arbitrarily large then this will eventually return an out-of-bounds array error.
# gfortran's uninitialized variable detection seems to not detect this.

# <https://gcc.gnu.org/onlinedocs/gcc/Optimize-Options.html>
# <https://gcc.gnu.org/onlinedocs/gfortran/Error-and-Warning-Options.html>

# `-march=native` may make the code not run on different machines, but that's not an issue for me.
# <https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html>
# <https://stackoverflow.com/questions/52653025/why-is-march-native-used-so-rarely>

# `-Warray-temporaries` might be useful when optimizing code.

# TODO: `-funroll`
# Removed: `-fmax-errors=1`

# static linking:
# `-static`
# <https://gcc.gnu.org/onlinedocs/gfortran/Link-Options.html>
# <https://fortran-lang.discourse.group/t/problems-with-build-of-static-exe/8228/5>
# <https://fortran-lang.discourse.group/t/makefile-errors-using-gfortran-static-option/3491>

# `-Wno-do-subscript`
# <https://gcc.gnu.org/bugzilla/show_bug.cgi?id=97320>
