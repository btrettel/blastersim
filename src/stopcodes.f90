! Module to set program stop (exit) codes.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [flt](https://github.com/btrettel/flt)
! License: [GPLv3](https://www.gnu.org/licenses/gpl-3.0.en.html)

module stopcodes

implicit none

! Following the [sysexits standard](https://man.freebsd.org/cgi/man.cgi?query=sysexits).
! I didn't copy all over as not all are relevant.

! > The successful exit is always indicated by a status of 0, or `EX_OK`.
integer, parameter :: EX_OK = 0

! > The command was used incorrectly, e.g., with the wrong number of arguments, a bad flag, a bad syntax in a parameter, or
! > whatever.
integer, parameter :: EX_USAGE = 64

! > The input data was incorrect in some way. This should only be used for user's data and not system files.
integer, parameter :: EX_DATAERR = 65

! > An input file (not a system file) did not exist or was not readable. This could also include errors like "No message" to a
! > mailer (if it cared to catch it).
integer, parameter :: EX_NOINPUT = 66

! > A service is unavailable. This can occur if  a support program or file does not exist. This can also be used as a catchall
! > message when something you wanted to do does not work, but you do not know why.
integer, parameter :: EX_UNAVAILABLE = 70

! > An internal software error has been detected. This should be limited to non-operating system related errors as possible.
integer, parameter :: EX_SOFTWARE = 70

! > A (user specified) output file cannot be created.
integer, parameter :: EX_CANTCREAT = 73

! > An error occurred while doing I/O on some file.
integer, parameter :: EX_IOERR = 74

! > Something was found in an unconfigured or misconfigured state.
integer, parameter :: EX_CONFIG = 78

! Later I may add exit codes not from the sysexits standard below. These will be chosen to not conflict.

end module stopcodes
