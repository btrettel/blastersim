# ifort is here due to possible performance benefits on x86.
# `-init=snan,arrays` seems to lead to false positives with ifort.

FFLAGS   = -diag-disable=10448 -warn errors -warn all -diag-error=remark,warn,error -fltconsistency -stand f18
DFLAGS   = -O0 -g -traceback -debug full -check all,noudio_iostat -fpe0 -standard-semantics
RFLAGS   = -O2 -parallel -ipo
AFLAGS   = 
NFLAGS   = -xHost
SFLAGS   = -static
OMPFLAGS = -qopenmp

# Removed: `-diag-error-limit=1`

# `-diag-disable=10448`: This is to avoid a warning about ifort being deprecated.
