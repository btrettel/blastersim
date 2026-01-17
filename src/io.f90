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

public :: write_latex_engineering
public :: create_barrel
public :: read_springer_namelist

contains

subroutine write_latex_engineering(tex_unit, x, macro_name)
    integer, intent(in)          :: tex_unit
    real(WP), intent(in)         :: x
    character(len=*), intent(in) :: macro_name
    
    logical  :: tex_unit_is_open
    integer  :: n
    real(WP) :: m
    
    inquire(unit=tex_unit, opened=tex_unit_is_open)
    call assert(tex_unit_is_open, "io (write_latex_engineering): tex_unit must be open")
    
    n = 3*nint(log(x) / (3.0_WP*log(10.0_WP)))
    m = x / (10.0_WP**n)
    
    write(unit=tex_unit, fmt="(3a, f8.3, a, i0, a)") "\newcommand*{\", macro_name, "}{", m, " \cdot 10^{", n, "}}"
end subroutine write_latex_engineering

subroutine create_barrel(vol_dead, csa_barrel, p_atm, temp_atm, m_p, p_fs, p_fd, l_travel, gas, i_cv_mirror, cv)
    use gasdata, only: gas_type
    use cva, only: cv_type
    
    type(si_volume), intent(in)      :: vol_dead    ! dead volume
    type(si_area), intent(in)        :: csa_barrel  ! cross-sectional area
    type(si_pressure), intent(in)    :: p_atm       ! atmospheric pressure
    type(si_temperature), intent(in) :: temp_atm    ! atmospheric temperature
    type(si_mass), intent(in)        :: m_p         ! projectile mass
    type(si_pressure)                :: p_fs, p_fd  ! pressures of static and dynamic friction, respectively
    type(si_length), intent(in)      :: l_travel    ! projectile travel
    type(gas_type), intent(in)       :: gas(:)      ! gas data
    integer, intent(in)              :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    type(cv_type), intent(out)       :: cv
    
    type(si_length)    :: x_d, x_stop, x_pre
    type(si_velocity)  :: x_dot
    type(unitless)     :: y(1)
    type(si_stiffness) :: k
    integer            :: n_d
    
    n_d = size(vol_dead%v%d)
    
    call x_dot%v%init_const(0.0_WP, n_d)
    call y(1)%v%init_const(1.0_WP, n_d)
    x_d = vol_dead/csa_barrel
    call k%v%init_const(0.0_WP, n_d)
    call x_pre%v%init_const(0.0_WP, n_d)
    x_stop = x_d + l_travel
    
    call cv%set(x_d, x_dot, y, p_atm, temp_atm, "barrel", csa_barrel, 1.0_WP/m_p, p_fs, p_fd, k, x_pre, gas, &
                    i_cv_mirror, x_stop=x_stop)
end subroutine create_barrel

subroutine read_springer_namelist(input_file, sys, rc)
    use, intrinsic :: iso_fortran_env, only: IOSTAT_END, ERROR_UNIT
    use cva, only: cv_system_type
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM_ => TEMP_ATM
    use checks, only: is_close, check
    use prec, only: CL
    
    character(len=*), intent(in)      :: input_file
    type(cv_system_type), intent(out) :: sys
    integer, intent(out)              :: rc
    
    include "geninput_springer_subroutine.f90"
    
    ! TODO: construct `sys`
end subroutine read_springer_namelist

end module io
