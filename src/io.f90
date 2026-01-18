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
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM_ => TEMP_ATM, DRY_AIR
    use checks, only: is_close, check
    use prec, only: CL, PI
    
    character(len=*), intent(in)                   :: input_file
    type(cv_system_type), allocatable, intent(out) :: sys
    integer, intent(out)                           :: rc
    
    integer :: n_d
    
    type(si_length)      :: d_plunger_u
    type(si_area)        :: csa_plunger
    type(si_pressure)    :: p_atm_u
    type(si_temperature) :: temp_atm_u
    
    include "geninput_springer_subroutine.f90"
    
    n_d = 0
    
    ! construct `sys`
    
    allocate(sys)
    allocate(sys%cv(4))
    allocate(sys%con(4, 4))
    
    ! `sys%con`
    
    sys%con(1, 1)%active = .false.
    sys%con(1, 2)%active = .false.
    sys%con(1, 3)%active = .false.
    sys%con(1, 4)%active = .false.
    
    sys%con(2, 1)%active = .false.
    sys%con(2, 2)%active = .false.
    sys%con(2, 3)%active = .false.
    sys%con(2, 4)%active = .false.
    
    sys%con(3, 1)%active = .false.
    sys%con(3, 2)%active = .false.
    sys%con(3, 3)%active = .false.
    sys%con(3, 4)%active = .true.
    call sys%con(3, 4)%a_e%v%init_const(a_e, n_d)
    call sys%con(3, 4)%b%v%init_const(b, n_d)
    
    sys%con(4, 1)%active = .false.
    sys%con(4, 2)%active = .false.
    sys%con(4, 3) = sys%con(3, 4)
    sys%con(4, 4)%active = .false.
    
    ! TODO: `sys%cv(1)`: atmosphere for plunger tube
    call d_plunger_u%v%init_const(d_plunger, n_d)
    csa_plunger = (PI/4.0_WP)*square(d_plunger_u)
    call p_atm_u%v%init_const(p_atm, n_d)
    call temp_atm_u%v%init_const(temp_atm, n_d)
    call sys%cv(1)%set_const("atmosphere for chamber", csa_plunger, p_atm_u, temp_atm_u, [DRY_AIR], 3)
    
    ! TODO: `sys%cv(2)`: atmosphere for barrel
    
    ! TODO: `sys%cv(3)`: plunger tube
    
    ! TODO: `sys%cv(4)`: barrel
end subroutine read_springer_namelist

end module io
