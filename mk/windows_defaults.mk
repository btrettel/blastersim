BINEXT   = .exe
RUN      = 
RM       = del /q
OFLAG    = -o
OBJEXT   = obj
OBJFLAGS = -c -o
MISSED   = missed.txt
CLEAN    = *.nml *.mod *.$(OBJEXT) src$(DIR_SEP)*.$(OBJEXT) $(MISSED) src$(DIR_SEP)rev.f90 rev.tex *$(BINEXT) *.pdb *.opt.yaml *.optrpt *.csv
CMP      = fc
PYTHON   = python
ARCH     = 
# The `^` escapes the backslash. Otherwise, `\` would be a new line in NMAKE.
DIR_SEP = ^\
GREP    = findstr /n
CP      = copy
