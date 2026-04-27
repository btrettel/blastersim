# <https://flang.llvm.org/docs/FlangCommandLineReference.html>
# <https://gist.github.com/scivision/43841c17cc47dfaef2817a331a31c90a>: > note that Flang doesn't have a `-Wall` like flag yet.
# <https://discourse.llvm.org/t/rfc-support-for-std-f2023/89608>: > The only choice today for an explicit standard version is `-std=f2018`, but as far as I can tell it doesn’t really do anything.

FFLAGS   = -pedantic -std=f2018 -Werror
DFLAGS   = -g
RFLAGS   = -O2 -flto
AFLAGS   = 
NFLAGS   = -march=native
SFLAGS   = -static
OMPFLAGS = -fopenmp
