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

!public :: read_springer_namelist

contains

subroutine create_barrel(vol_d, d, p_atm, temp_atm, m_p, p_fs, p_fd, l, n_d, cv)
    use prec, only: PI
    use gasdata, only: DRY_AIR
    use cva, only: cv_type
    
    real(WP), intent(in)       :: vol_d      ! dead volume (m3)
    real(WP), intent(in)       :: d          ! barrel diameter (m)
    real(WP), intent(in)       :: p_atm      ! atmospheric pressure (Pa), temp_atm, m_p, p_fs, p_fd, l
    real(WP), intent(in)       :: temp_atm   ! atmospheric temperature (K)
    real(WP), intent(in)       :: m_p        ! projectile mass (kg)p_fs, p_fd, l
    real(WP), intent(in)       :: p_fs, p_fd ! pressures of static and dynamic friction, respectively (Pa)
    real(WP), intent(in)       :: l          ! barrel length (projectile travel) (m)
    integer, intent(in)        :: n_d
    type(cv_type), intent(out) :: cv
    
    type(si_volume)      :: vol_d_
    type(si_length)      :: d_, x_d, l_, x_stop, x_z
    type(si_area)        :: csa
    type(si_velocity)    :: x_dot
    type(unitless)       :: y(1)
    type(si_pressure)    :: p_atm_, p_fs_, p_fd_
    type(si_temperature) :: temp_atm_
    type(si_mass)        :: m_p_
    type(si_stiffness)   :: k
    
    call vol_d_%v%init_const(vol_d, n_d)
    call d_%v%init_const(d, n_d)
    csa = (PI/4.0_WP)*square(d_)
    x_d = vol_d_/csa
    call x_dot%v%init_const(0.0_WP, n_d)
    call y(1)%v%init_const(1.0_WP, n_d)
    call p_atm_%v%init_const(p_atm, n_d)
    call temp_atm_%v%init_const(temp_atm, n_d)
    call m_p_%v%init_const(m_p, n_d)
    call p_fs_%v%init_const(p_fs, n_d)
    call p_fd_%v%init_const(p_fd, n_d)
    call l_%v%init_const(l, n_d)
    x_stop = x_d + l_
    call k%v%init_const(0.0_WP, n_d)
    call x_z%v%init_const(0.0_WP, n_d)
    
    call cv%set(x_d, x_dot, y, p_atm_, temp_atm_, csa, 1.0_WP/m_p_, p_fs_, p_fd_, p_atm_, k, x_z, [DRY_AIR], x_stop)
end subroutine create_barrel

!subroutine read_springer_namelist(input_file, sys, rc)
!    use gasdata, only: P_ATM_ => P_ATM, T_ATM_ => T_ATM
    
!    include "geninput_springer_subroutine.f90"
!end subroutine read_springer_namelist

end module io
