FFLAGS   = /Qdiag-disable:10448 /warn:errors /warn:all /Qdiag-error:remark,warn,error /fltconsistency /stand:f18
DFLAGS   = /O0 /traceback /debug:full /check:all,noudio_iostat /fpe:0 /standard-semantics
RFLAGS   = /O2 /Qparallel /Qipo
AFLAGS   = 
NFLAGS   = /QxHost
SFLAGS   = /static
OMPFLAGS = /Qopenmp
OFLAG    = /o
OBJFLAGS = /c /o

# Removed: `/Qdiag-error-limit:1`
