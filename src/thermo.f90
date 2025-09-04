! Module for thermodynamic calculations.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module thermo

use prec, only: WP
implicit none
private

! <https://en.wikipedia.org/wiki/Gas_constant>
real(WP), public, parameter :: R_BAR = 8.31446261815324_WP ! J/(mol*K)

! <https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html>
real(WP), public, parameter :: M_AIR = 28.9647e-3_WP ! kg/mol

public :: p_eos

contains

pure function p_eos(rho, temp)
    use fmad, only: ad
    use units, only: si_mass_density  => unit_m30_p10_p00_p00, &
                     si_temperature   => unit_p00_p00_p00_p10, &
                     si_pressure      => unit_m10_p10_m20_p00, &
                     si_specific_heat => unit_p20_p00_m20_m10
    
    type(si_mass_density), intent(in) :: rho
    type(si_temperature), intent(in)  :: temp
    type(si_pressure)                 :: p_eos
    
    integer :: n_dv
    
    type(si_specific_heat) :: R
    
    n_dv = size(rho%v%d)
    
    call R%v%init_const(R_BAR/M_AIR, n_dv)
    
    p_eos = rho * R * temp
end function p_eos

end module thermo
