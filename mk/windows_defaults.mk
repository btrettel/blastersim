BINEXT   = .exe
RUN      = 
RM       = del /q
OFLAG    = -o
OBJEXT   = obj
OBJFLAGS = -c -o
MISSED   = missed.txt
CLEAN    = *.nml *.mod *.$(OBJEXT) src$(DIR_SEP)*.$(OBJEXT) $(MISSED) src$(DIR_SEP)rev.f90 *$(BINEXT) *.pdb *.opt.yaml *.optrpt
CMP      = fc
PYTHON   = python
ARCH     = 

# The `^` escapes the backslash. Otherwise, `\` would be a new line in NMAKE.
DIR_SEP = ^\
