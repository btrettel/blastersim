# Note: The Intel compilers will not compile BlasterSim as of this writing. This file is present for testing and does not imply that the Intel compilers work with BlasterSim.

FFLAGS    = /Qdiag-disable:10448 /warn:errors /warn:all /Qdiag-error:remark,warn,error /fltconsistency /stand:f18
DFLAGS    = /O0 /traceback /debug:full /check:all,noudio_iostat /fpe:0 /standard-semantics
RFLAGS    = /O2 /Qparallel /Qipo
FUZZFLAGS = /O2 /Qparallel /Qipo /traceback /debug:full /check:all,noudio_iostat /fpe:0
AFLAGS    = 
NFLAGS    = /QxHost
SFLAGS    = /static
OMPFLAGS  = /Qopenmp
OFLAG     = /o
OBJFLAGS  = /c /o

# Removed: `/Qdiag-error-limit:1`
