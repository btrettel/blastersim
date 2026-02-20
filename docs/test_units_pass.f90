program test_units_pass

use units
implicit none

type(si_length)   :: x
type(si_time)     :: t
type(si_velocity) :: v

x%v = 1.0
t%v = 1.0

v = x / t

print *, v

end program test_units_pass
