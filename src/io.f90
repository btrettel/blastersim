! Module for input and output subroutines.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module io

use prec, only: WP
use units
use checks, only: assert
implicit none

public :: create_barrel
!public :: read_springer_namelist

contains

subroutine create_barrel(vol_dead, d_barrel, p_atm, temp_atm, m_p, p_fs, p_fd, l_travel, n_d, cv)
    use prec, only: PI
    use gasdata, only: DRY_AIR
    use cva, only: cv_type
    
    type(si_volume), intent(in)      :: vol_dead   ! dead volume
    type(si_length), intent(in)      :: d_barrel   ! barrel diameter
    type(si_pressure), intent(in)    :: p_atm      ! atmospheric pressure
    type(si_temperature), intent(in) :: temp_atm   ! atmospheric temperature
    type(si_mass), intent(in)        :: m_p        ! projectile mass
    type(si_pressure)                :: p_fs, p_fd ! pressures of static and dynamic friction, respectively
    type(si_length), intent(in)      :: l_travel   ! projectile travel
    integer, intent(in)              :: n_d
    type(cv_type), intent(out)       :: cv
    
    type(si_length)    :: x_d, x_stop, x_z
    type(si_area)      :: csa_barrel
    type(si_velocity)  :: x_dot
    type(unitless)     :: y(1)
    type(si_stiffness) :: k
    
    call x_dot%v%init_const(0.0_WP, n_d)
    call y(1)%v%init_const(1.0_WP, n_d)
    csa_barrel = (PI/4.0_WP)*square(d_barrel)
    x_d = vol_dead/csa_barrel
    call k%v%init_const(0.0_WP, n_d)
    call x_z%v%init_const(0.0_WP, n_d)
    x_stop = x_d + l_travel
    
    call cv%set(x_d, x_dot, y, p_atm, temp_atm, csa_barrel, 1.0_WP/m_p, p_fs, p_fd, p_atm, k, x_z, [DRY_AIR], x_stop)
end subroutine create_barrel

!subroutine read_springer_namelist(input_file, sys, rc)
!    use gasdata, only: P_ATM_ => P_ATM, T_ATM_ => T_ATM
    
!    include "geninput_springer_subroutine.f90"
!end subroutine read_springer_namelist

end module io
