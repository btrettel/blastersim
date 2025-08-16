FFLAGS   = -Error=Any_Warning -f2018 -info -quiet -I src -w=ques
DFLAGS   = -g
RFLAGS   = -O2
AFLAGS   = 
NFLAGS   = -target=native
SFLAGS   = -Bstatic
OMPFLAGS = -openmp

# `-I src`: Required because unlike other compilers, nagfor won't look for `include`d files in the same directory as the source code calling it.

# `-w=ques`: Added because nagfor doesn't like mixing `real`s and `integer`s in `system_clock`:
# `Questionable: src/timer.f90, line 63: Argument COUNT (no. 1) to intrinsic subroutine SYSTEM_CLOCK is of type INTEGER(int64), but argument COUNT_RATE (no. 2) is of type DOUBLE PRECISION`
# I don't agree that this is questionable and might ask them to remove this message.
