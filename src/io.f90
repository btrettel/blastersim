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
public :: read_springer_namelist, read_pneumatic_namelist

contains

subroutine write_latex_engineering(tex_unit, x, macro_name, m_fmt)
    integer, intent(in)          :: tex_unit
    real(WP), intent(in)         :: x
    character(len=*), intent(in) :: macro_name, m_fmt
    
    logical  :: tex_unit_is_open
    integer  :: n
    real(WP) :: m
    
    inquire(unit=tex_unit, opened=tex_unit_is_open)
    call assert(tex_unit_is_open, "io (write_latex_engineering): tex_unit must be open")
    
    n = 3*nint(log(x) / (3.0_WP*log(10.0_WP)))
    m = x / (10.0_WP**n)
    
    write(unit=tex_unit, fmt="(3a, " // trim(m_fmt) // ", a, i0, a)") "\newcommand*{\", macro_name, "}{", &
                                                                        m, " \cdot 10^{", n, "}}"
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

subroutine read_springer_namelist(input_file, sys, config, rc)
    use, intrinsic :: iso_fortran_env, only: IOSTAT_END, ERROR_UNIT
    use cva, only: cv_system_type, run_config_type, DT_DEFAULT
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM_ => TEMP_ATM, DRY_AIR
    use checks, only: is_close, check
    use prec, only: CL, PI
    
    character(len=*), intent(in)                   :: input_file
    type(cv_system_type), allocatable, intent(out) :: sys
    type(run_config_type), intent(out)             :: config
    integer, intent(out)                           :: rc
    
    type(si_velocity) :: x_dot
    type(unitless)    :: y(1)
    type(si_area)     :: csa_plunger, csa_barrel
    
    integer, parameter :: I_PLUNGER_ATM = 1, I_BARREL_ATM  = 2, I_PLUNGER = 3, I_BARREL = 4
    
    include "geninput_springer.f90"
    
    ! construct `sys`
    
    allocate(sys)
    allocate(sys%cv(4))
    allocate(sys%con(4, 4))
    
    ! `sys%con`
    
    sys%con(I_PLUNGER_ATM, I_PLUNGER_ATM)%active = .false.
    sys%con(I_PLUNGER_ATM, I_BARREL_ATM)%active  = .false.
    sys%con(I_PLUNGER_ATM, I_PLUNGER)%active     = .false.
    sys%con(I_PLUNGER_ATM, I_BARREL)%active      = .false.
    
    sys%con(I_BARREL_ATM, I_PLUNGER_ATM)%active = .false.
    sys%con(I_BARREL_ATM, I_BARREL_ATM)%active  = .false.
    sys%con(I_BARREL_ATM, I_PLUNGER)%active     = .false.
    sys%con(I_BARREL_ATM, I_BARREL)%active      = .false.
    
    sys%con(I_PLUNGER, I_PLUNGER_ATM)%active = .false.
    sys%con(I_PLUNGER, I_BARREL_ATM)%active  = .false.
    sys%con(I_PLUNGER, I_PLUNGER)%active     = .false.
    sys%con(I_PLUNGER, I_BARREL)%active      = .true.
    sys%con(I_PLUNGER, I_BARREL)%a_e         = (PI/4.0_WP)*square(d_e_u)
    sys%con(I_PLUNGER, I_BARREL)%b           = b_u
    call sys%con(I_PLUNGER, I_BARREL)%t_opening%v%init_const(0.0_WP, 0)
    call sys%con(I_PLUNGER, I_BARREL)%alpha_0%v%init_const(1.0_WP, 0)
    call sys%con(I_PLUNGER, I_BARREL)%alpha_dot_0%v%init_const(0.0_WP, 0)
    call sys%con(I_PLUNGER, I_BARREL)%m_dot_0%v%init_const(0.0_WP, 0)
    
    sys%con(I_BARREL, I_PLUNGER_ATM)%active = .false.
    sys%con(I_BARREL, I_BARREL_ATM)%active  = .false.
    sys%con(I_BARREL, I_PLUNGER)            = sys%con(I_PLUNGER, I_BARREL)
    sys%con(I_BARREL, I_BARREL)%active      = .false.
    
    ! The same for every control volume.
    call x_dot%v%init_const(0.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    
    ! `sys%cv(I_PLUNGER_ATM)`: atmosphere for plunger tube
    csa_plunger = (PI/4.0_WP)*square(d_plunger_u)
    call sys%cv(I_PLUNGER_ATM)%set_const("atmosphere for plunger", csa_plunger, p_atm_u, temp_atm_u, [DRY_AIR], y, I_PLUNGER)
    
    ! `sys%cv(I_BARREL_ATM)`: atmosphere for barrel
    csa_barrel = (PI/4.0_WP)*square(d_barrel_u)
    call sys%cv(I_BARREL_ATM)%set_const("atmosphere for barrel", csa_barrel, p_atm_u, temp_atm_u, [DRY_AIR], y, I_BARREL)
    
    ! `sys%cv(I_PLUNGER)`: plunger tube
    call sys%cv(I_PLUNGER)%set(l_draw_u, x_dot, y, p_atm_u, temp_atm_u, "plunger tube", csa_plunger, &
                        1.0_WP/m_plunger_u, p_fs_plunger_u, p_fd_plunger_u, k_u, l_pre_u, [DRY_AIR], I_PLUNGER_ATM, &
                        m_spring=m_spring_u)
    
    ! `sys%cv(I_BARREL)`: barrel
    call create_barrel(vol_dead_u, csa_barrel, p_atm_u, temp_atm_u, m_proj_u, p_fs_proj_u, p_fd_proj_u, l_travel_u, &
                        [DRY_AIR], I_BARREL_ATM, sys%cv(I_BARREL))
    
    call config%set(id, csv_output=.true., dt=dt_u, n_d=0)
end subroutine read_springer_namelist

subroutine read_pneumatic_namelist(input_file, sys, config, rc)
    use, intrinsic :: iso_fortran_env, only: IOSTAT_END, ERROR_UNIT
    use cva, only: cv_system_type, run_config_type, DT_DEFAULT
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM_ => TEMP_ATM, DRY_AIR
    use checks, only: is_close, check
    use prec, only: CL, PI
    
    character(len=*), intent(in)                   :: input_file
    type(cv_system_type), allocatable, intent(out) :: sys
    type(run_config_type), intent(out)             :: config
    integer, intent(out)                           :: rc
    
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_area)         :: csa_chamber, csa_barrel
    type(si_inverse_mass) :: rm_p
    type(si_stiffness)    :: k
    type(si_length)       :: l_pre, x_chamber
    type(si_pressure)     :: p_f_chamber
    
    integer, parameter :: I_BARREL_ATM  = 1, I_CHAMBER = 2, I_BARREL = 3
    
    include "geninput_pneumatic.f90"
    
    ! construct `sys`
    
    allocate(sys)
    allocate(sys%cv(3))
    allocate(sys%con(3, 3))
    
    ! `sys%con`
    
    sys%con(I_BARREL_ATM, I_BARREL_ATM)%active  = .false.
    sys%con(I_BARREL_ATM, I_CHAMBER)%active     = .false.
    sys%con(I_BARREL_ATM, I_BARREL)%active      = .false.
    
    sys%con(I_CHAMBER, I_BARREL_ATM)%active  = .false.
    sys%con(I_CHAMBER, I_CHAMBER)%active     = .false.
    sys%con(I_CHAMBER, I_BARREL)%active      = .true.
    sys%con(I_CHAMBER, I_BARREL)%a_e         = (PI/4.0_WP)*square(d_e_u)
    sys%con(I_CHAMBER, I_BARREL)%b           = b_u
    sys%con(I_CHAMBER, I_BARREL)%t_opening = t_opening_u
    call sys%con(I_CHAMBER, I_BARREL)%alpha_0%v%init_const(0.0_WP, 0)
    call sys%con(I_CHAMBER, I_BARREL)%alpha_dot_0%v%init_const(1.0_WP, 0)
    call sys%con(I_CHAMBER, I_BARREL)%m_dot_0%v%init_const(0.0_WP, 0)
    
    sys%con(I_BARREL, I_BARREL_ATM)%active  = .false.
    sys%con(I_BARREL, I_CHAMBER)            = sys%con(I_CHAMBER, I_BARREL)
    sys%con(I_BARREL, I_BARREL)%active      = .false.
    
    call x_dot%v%init_const(0.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call rm_p%v%init_const(0.0_WP, 0) ! immobile
    call k%v%init_const(0.0_WP, 0)
    call l_pre%v%init_const(0.0_WP, 0)
    call p_f_chamber%v%init_const(0.0_WP, 0)
    
    ! `sys%cv(I_BARREL_ATM)`: atmosphere
    csa_barrel = (PI/4.0_WP)*square(d_barrel_u)
    call sys%cv(I_BARREL_ATM)%set_const("atmosphere", csa_barrel, p_atm_u, temp_atm_u, [DRY_AIR], y, I_BARREL)
    
    ! `sys%cv(I_CHAMBER)`: chamber
    csa_chamber = (PI/4.0_WP)*square(d_chamber_u)
    x_chamber   = vol_chamber_u/csa_chamber
    call sys%cv(I_CHAMBER)%set(x_chamber, x_dot, y, p_chamber_u, temp_atm_u, "chamber", csa_chamber, &
                        rm_p, p_f_chamber, p_f_chamber, k, l_pre, [DRY_AIR], 0, &
                        isentropic_filling=.true., p_atm=p_atm_u)
    
    ! `sys%cv(I_BARREL)`: barrel
    call create_barrel(vol_dead_u, csa_barrel, p_atm_u, temp_atm_u, m_proj_u, p_fs_proj_u, p_fd_proj_u, l_travel_u, &
                        [DRY_AIR], I_BARREL_ATM, sys%cv(I_BARREL))
    
    call config%set(id, csv_output=.true., dt=dt_u, n_d=0)
end subroutine read_pneumatic_namelist

end module io
