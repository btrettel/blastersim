program test_units_fail

use units
implicit none

type(si_length) :: x, v
type(si_time)   :: t

x%v = 1.0
t%v = 1.0

v = x / t

print *, v

end program test_units_fail
