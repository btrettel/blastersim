! tests for the port module
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_thermo

use fmad, only: ad
use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                 si_temperature  => unit_p00_p00_p00_p10, &
                 si_pressure     => unit_m10_p10_m20_p00
use prec, only: WP
use thermo, only: p_eos
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("thermo.nml")

call test_p_eos(tests)

call tests%end_tests()

contains

subroutine test_p_eos(tests)
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    
    ! <https://en.wikipedia.org/wiki/Density_of_air>
    call rho%v%init_const(1.2250_WP, 0)
    call temp%v%init_const(273.15_WP + 15.0_WP, 0)
    
    p = p_eos(rho, temp)
    
    call tests%real_eq(p%v%v, 101.325e3_WP, "p_eos, atmospheric", abs_tol=1.0_WP)
end subroutine test_p_eos

end program test_thermo
