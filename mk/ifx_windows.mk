FFLAGS    = /warn:errors /warn:all /Qdiag-error:remark,warn,error /fltconsistency /stand:f18
DFLAGS    = /Qdiag-disable:10440 /O0 /traceback /debug:full /check:all /fpe:0 /standard-semantics
RFLAGS    = /O2 /Qopt-report /Qflto
FUZZFLAGS = /O2 /Qflto /traceback /debug:full /check:all /fpe:0
AFLAGS    = 
NFLAGS    = /QxHost
SFLAGS    = /static
OMPFLAGS  = /Qiopenmp
OFLAG     = /o
OBJFLAGS  = /c /o

# Removed: `/Qdiag-error-limit:1`
