! tests for unit conversions
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_convert

use prec, only: WP
use units
use convert
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

type(si_temperature) :: temp
type(si_pressure)    :: p
type(si_length)      :: x
type(si_volume)      :: vol

call tests%start_tests("convert.nml")

temp = celsius_const(20.0_WP, 0)
p    = psi_const(14.696_WP, 0)
x    = inch_const(39.37_WP, 0)
vol  = cubic_inches_const(231.0_WP, 0) ! 1 gallon

call tests%real_eq(temp%v%v, 293.15_WP, "celsius_const")
call tests%real_eq(p%v%v, 101325.0_WP, "psi_const", abs_tol=1.0_WP)
call tests%real_eq(x%v%v, 1.0_WP, "inch_const", abs_tol = 1.0e-5_WP)
call tests%real_eq(vol%v%v, 0.003785411784_WP, "cubic_inches_const", abs_tol=1.0e-7_WP)

call tests%end_tests()

end program test_convert
