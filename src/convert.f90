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

public :: celsius_const, psi_const, inch_const, cubic_inches_const, fps_const, lbf_per_in_const

real(WP), public, parameter :: CONVERT_C_TO_K = 273.15_WP ! K, temperature to add to convert from C to K

real(WP), parameter :: CONVERT_PSI_TO_PA = 6894.8_WP    ! Pa/psi
real(WP), parameter :: CONVERT_IN_TO_M   = 2.54e-2_WP   ! m/in
real(WP), parameter :: CONVERT_IN3_TO_M3 = 16.387e-6_WP ! m3/in3
real(WP), parameter :: CONVERT_FT_TO_M   = 0.3048_WP    ! m/ft
real(WP), parameter :: CONVERT_LB_TO_N   = 0.22481_WP   ! N/lbf

contains

pure function celsius_const(temp, n_d)
    real(WP), intent(in) :: temp
    integer, intent(in)  :: n_d
    
    type(si_temperature) :: celsius_const
    
    call celsius_const%v%init_const(CONVERT_C_TO_K + temp, n_d)
    
    call assert(celsius_const%v%v > 0.0_WP, "convert (celsius_const): temperature below or equal to absolute zero")
end function celsius_const

pure function psi_const(p, n_d)
    real(WP), intent(in) :: p
    integer, intent(in)  :: n_d
    
    type(si_pressure) :: psi_const
    
    call psi_const%v%init_const(CONVERT_PSI_TO_PA*p, n_d)
end function psi_const

pure function inch_const(x, n_d)
    real(WP), intent(in) :: x
    integer, intent(in)  :: n_d
    
    type(si_length) :: inch_const
    
    call inch_const%v%init_const(CONVERT_IN_TO_M*x, n_d)
end function inch_const

pure function cubic_inches_const(vol, n_d)
    real(WP), intent(in) :: vol
    integer, intent(in)  :: n_d
    
    type(si_volume) :: cubic_inches_const
    
    call cubic_inches_const%v%init_const(CONVERT_IN3_TO_M3*vol, n_d)
end function cubic_inches_const

pure function fps_const(v, n_d)
    real(WP), intent(in) :: v
    integer, intent(in)  :: n_d
    
    type(si_velocity) :: fps_const
    
    call fps_const%v%init_const(CONVERT_FT_TO_M*v, n_d)
end function fps_const

pure function lbf_per_in_const(k, n_d)
    real(WP), intent(in) :: k
    integer, intent(in)  :: n_d
    
    type(si_stiffness) :: lbf_per_in_const
    
    call lbf_per_in_const%v%init_const(CONVERT_LB_TO_N*k/CONVERT_IN_TO_M, n_d)
end function lbf_per_in_const

end module convert
