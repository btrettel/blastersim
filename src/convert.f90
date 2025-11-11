! Module to for unit conversions.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module convert

use prec, only: WP
use units
use checks, only: assert
implicit none
private

public :: set_celsius_const

real(WP), public, parameter :: TEMP_C_TO_K = 273.15_WP ! K, temperature to add to convert from C to K

contains

pure function set_celsius_const(temp, n_d)
    real(WP), intent(in) :: temp
    integer, intent(in)  :: n_d
    
    type(si_temperature) :: set_celsius_const
    
    call assert(temp > (-TEMP_C_TO_K), "convert (set_celsius_const): temperature below absolute zero")
    
    call set_celsius_const%v%init_const(TEMP_C_TO_K + temp, n_d)
end function set_celsius_const

end module convert
