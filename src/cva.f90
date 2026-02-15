! Module for control volume analysis.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module cva

use prec, only: WP
use units
use checks, only: assert, assert_dimension, is_close
use gasdata, only: gas_type
implicit none
private

public :: smooth_min
public :: d_x_d_t, d_x_dot_d_t, d_m_k_d_t, d_e_d_t, d_e_f_d_t
public :: f_m_dot, g_m_dot
public :: time_step, run, check_sys, write_csv_row

! pressure ratio laminar flow nominally starts at
! based on first part of beater_pneumatic_2007 eq. 5.4
real(WP), public, parameter :: P_RL = 0.999_WP ! unitless

! fraction of spring mass to add to plunger mass to get effective mass
! ruby_equivalent_2000
real(WP), public, parameter :: C_MS = 1.0_WP/3.0_WP ! unitless

real(WP), public, parameter :: X_STOP_DEFAULT         = 1.0e3_WP   ! m (If you have a barrel that's a km long, that's probably wrong.)
real(WP), public, parameter :: DT_DEFAULT             = 1.0e-5_WP  ! s
real(WP), public, parameter :: T_STOP_DEFAULT         = 0.5_WP     ! s
real(WP), public, parameter :: MASS_TOLERANCE         = 1.0e-5_WP  ! unitless
real(WP), public, parameter :: ENERGY_TOLERANCE       = 1.0e-4_WP  ! unitless
real(WP), public, parameter :: MASS_DERIV_TOLERANCE   = 1.0e-8_WP  ! unitless
real(WP), public, parameter :: ENERGY_DERIV_TOLERANCE = 1.0e-5_WP  ! unitless (TODO: decrease later and see what breaks)
real(WP), public, parameter :: MIRROR_X_TOLERANCE     = 1.0e-12_WP ! unitless

integer, public, parameter :: IDEAL_EOS = 1 ! ideal gas equation of state
integer, public, parameter :: CONST_EOS = 2 ! constant pressure, temperature, density
integer, public, parameter :: RK_EOS    = 3 ! Redlich–Kwong equation of state (not yet implemented)
integer, public, parameter :: MAX_EOS   = 2

integer, public, parameter :: NORMAL_CV_TYPE = 1
integer, public, parameter :: MIRROR_CV_TYPE = 2
integer, public, parameter :: MAX_CV_TYPE    = 2

integer, public, parameter :: CONTINUE_RUN_RC               = -1
integer, public, parameter :: SUCCESS_RUN_RC                = 0
integer, public, parameter :: TIMEOUT_RUN_RC                = 1
integer, public, parameter :: NEGATIVE_CV_M_TOTAL_RUN_RC    = 2
integer, public, parameter :: NEGATIVE_CV_TEMP_RUN_RC       = 3
integer, public, parameter :: MASS_TOLERANCE_RUN_RC         = 4
integer, public, parameter :: ENERGY_TOLERANCE_RUN_RC       = 5
integer, public, parameter :: MASS_DERIV_TOLERANCE_RUN_RC   = 6
integer, public, parameter :: ENERGY_DERIV_TOLERANCE_RUN_RC = 7
integer, public, parameter :: IDEAL_EOS_RUN_RC              = 8
integer, public, parameter :: MIRROR_X_TOLERANCE_RUN_RC     = 9
integer, public, parameter :: X_BLOW_UP_RUN_RC              = 10
integer, public, parameter :: X_DOT_BLOW_UP_RUN_RC          = 11
integer, public, parameter :: M_BLOW_UP_RUN_RC              = 12
integer, public, parameter :: E_BLOW_UP_RUN_RC              = 13
integer, public, parameter :: E_F_BLOW_UP_RUN_RC            = 14

integer, public, parameter :: HEADER_ROW_TYPE = 1
integer, public, parameter :: NUMBER_ROW_TYPE = 2

type, public :: cv_type ! control volume
    ! time varying
    type(si_length)            :: x     ! location of projectile/plunger
    type(si_velocity)          :: x_dot ! velocity of projectile/plunger
    type(si_mass), allocatable :: m(:)  ! mass(es) of gas(es) in control volume
    type(si_energy)            :: e     ! energy of gas in control volume
    type(si_energy)            :: e_f   ! energy lost to projectile/plunger friction in control volume
    
    ! constants
    character(len=32)           :: label       ! human-readable label for control volume
    integer                     :: eos         ! equation of state to use for control volume
    integer                     :: type        ! type of control volume
    type(si_area)               :: csa         ! cross-sectional area
    type(si_inverse_mass)       :: rm_p        ! reciprocal mass of projectile/plunger
    type(si_pressure)           :: p_fs, p_fd  ! static and dynamic friction pressure
    type(si_stiffness)          :: k           ! stiffness of spring attached to plunger
    type(si_length)             :: l_pre       ! spring precompression length
    type(gas_type), allocatable :: gas(:)      ! gas data
    integer                     :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    ! `i_cv_mirror = 0` disables mirror CVs. Use that for constant volume chambers.
    type(si_pressure)           :: p_const     ! if `cv%eos = CONST_EOS`, then `cv%p() = p_const`
    type(si_temperature)        :: temp_const  ! if `cv%eos = CONST_EOS`, then `cv%temp() = temp_const`
    type(unitless), allocatable :: y_const(:)  ! if `cv%eos = CONST_EOS`, then this mass fraction will be used
    type(si_length)             :: x_stop      ! `x` location where simulation will stop
    type(si_mass)               :: m_spring    ! mass of spring
    logical                     :: constant_friction ! whether `p_f` will be constant or not
contains
    procedure :: m_total
    procedure :: e_s
    procedure :: e_p
    procedure :: e_total
    procedure :: p_eos
    procedure :: rho_eos
    procedure :: p_c
    procedure :: y
    procedure :: chi
    procedure :: r     => r_cv
    procedure :: temp  => temp_cv
    procedure :: vol   => vol_cv
    procedure :: rho   => rho_cv
    procedure :: p     => p_cv
    procedure :: u     => u_cv
    procedure :: h     => h_cv
    procedure :: gamma => gamma_cv
    procedure :: set
    procedure :: set_const
    procedure :: p_f
    procedure :: p_f0
end type cv_type

type :: cv_delta_type
    type(si_length)            :: x     ! delta of location of projectile/plunger
    type(si_velocity)          :: x_dot ! delta of velocity of projectile/plunger
    type(si_mass), allocatable :: m(:)  ! delta of mass(es) of gas(es) in control volume
    type(si_energy)            :: e     ! delta of energy of gas in control volume
    type(si_energy)            :: e_f   ! delta of energy lost to projectile/plunger friction in control volume
end type cv_delta_type

!tripwire$ begin 75CD9D0A Update sections of docs listed in source when adding valve opening model using poppet motion.
! usage.tex `\secref{pneumatic}`
! theory.tex `\secref{valve-opening-model}`
type, public :: con_type ! connection between control volumes
    logical        :: active
    type(si_area)  :: a_e         ! effective area
    type(unitless) :: b           ! critical pressure ratio
    type(si_time)  :: t_opening   ! valve opening time
    type(unitless) :: alpha_0     ! valve opening fraction at time zero
    type(unitless) :: alpha_0_dot ! valve opening rate at time zero
contains
    procedure :: m_dot
end type con_type
!tripwire$ end

type, public :: cv_system_type
    type(cv_type), allocatable  :: cv(:)
    type(con_type), allocatable :: con(:, :)
contains
    procedure :: calculate_flows
    procedure :: m_total => m_total_sys
    procedure :: e_total => e_total_sys
end type cv_system_type

type, public :: run_config_type
    character(len=128) :: id ! CSV file name
    logical            :: csv_output
    integer            :: csv_frequency
    type(si_time)      :: t_stop, dt
    logical            :: tolerance_checks
contains
    procedure :: set => set_run_config
end type run_config_type

type, public :: run_status_type
    integer       :: rc
    type(si_time) :: t
    integer, allocatable  :: i_cv(:) ! control volume(s) with associated error
    real(WP), allocatable :: data(:) ! additional error data
end type run_status_type

contains

!!!!!!!!!!!!!!!!!!!!!!!!!
! methods for `cv_type` !
!!!!!!!!!!!!!!!!!!!!!!!!!

pure function m_total(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_mass) :: m_total
    
    integer :: k_gas
    
    call m_total%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do k_gas = 1, size(cv%m)
        m_total = m_total + cv%m(k_gas)
    end do
end function m_total

pure function e_s(cv)
    ! Spring potential energy
    
    class(cv_type), intent(in) :: cv
    
    type(si_energy) :: e_s
    
    if (cv%type == MIRROR_CV_TYPE) then
        call e_s%v%init_const(0.0_WP, size(cv%x_dot%v%d))
    else
        e_s = 0.5_WP*cv%k*square(cv%x + cv%l_pre)
    end if
end function e_s

pure function e_p(cv)
    ! Projectile/plunger kinetic energy including added mass of spring
    
    class(cv_type), intent(in) :: cv
    
    type(si_energy)       :: e_p
    type(si_inverse_mass) :: r_mp_eff ! effective mass of projectile/plunger
    
    if (is_close(cv%rm_p%v%v, 0.0_WP) .or. (cv%type == MIRROR_CV_TYPE)) then
        if (is_close(cv%rm_p%v%v, 0.0_WP)) then
            ! `MIRROR_CV_TYPE` might simply copy what the other CV has, so it might not have no inverse mass.
            call assert(is_close(cv%k%v%v, 0.0_WP), "cva (e_p): if the plunger is immobile, k should be zero", &
                            print_real=[cv%k%v%v])
        end if
        
        call e_p%v%init_const(0.0_WP, size(cv%x_dot%v%d))
    else
        r_mp_eff = cv%rm_p / (1.0_WP + C_MS * cv%m_spring * cv%rm_p)
        e_p = 0.5_WP*square(cv%x_dot)/r_mp_eff
    end if
end function e_p

pure function e_total(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_energy) :: e_total
    
    e_total = cv%e + cv%e_f + cv%e_s() + cv%e_p()
end function e_total

pure function p_eos(cv, rho, temp)
    ! Calculate pressure using the equation of state.
    
    class(cv_type), intent(in)        :: cv
    type(si_mass_density), intent(in) :: rho
    type(si_temperature), intent(in)  :: temp
    type(si_pressure)                 :: p_eos
    
    integer :: n_d ! number of derivatives
    
    select case (cv%eos)
        case (IDEAL_EOS)
            call assert(rho%v%v  > 0.0_WP, "cva (p_eos, IDEAL_EOS): rho%v > 0 violated", print_real=[rho%v%v])
            call assert(temp%v%v > 0.0_WP, "cva (p_eos, IDEAL_EOS): temp%v > 0 violated", print_real=[temp%v%v])
            call assert_dimension(rho%v%d, temp%v%d)
            
            n_d   = size(rho%v%d)
            p_eos = rho * cv%r() * temp
        case (CONST_EOS)
            call assert(is_close(temp%v%v, cv%temp_const%v%v), "cva (p_eos, CONST_EOS): temp /= temp_const", &
                            print_real=[temp%v%v, cv%temp_const%v%v])
            p_eos = cv%p_const
        case default
            error stop "cva (p_eos): invalid cv%eos"
    end select
    
    call assert(p_eos%v%v > 0.0_WP, "cva (p_eos): p_eos%v > 0 violated", print_real=[p_eos%v%v])
end function p_eos

pure function rho_eos(cv, p, temp, y)
    ! Calculate density using the equation of state.
    ! In the future, if an EOS more complex than the ideal gas law is used, it might make sense to calculate `rho_eos` from `p_eos`
    ! with the Newton method.
    
    use gasdata, only: R_BAR
    
    class(cv_type), intent(in)       :: cv
    type(si_pressure), intent(in)    :: p
    type(si_temperature), intent(in) :: temp
    type(unitless), intent(in)       :: y(:)
    
    type(si_mass_density) :: rho_eos
    
    integer                     :: n_d, i, j
    type(si_molar_mass)         :: mm, gas_mm ! molar mass
    type(si_molar_gas_constant) :: r_bar_ ! universal gas constant
    type(unitless)              :: denominator, y_sum
    type(si_specific_heat)      :: r_cv
    
    ! Don't check `p_c` here as that requires the masses.
    
    call assert(cv%eos == IDEAL_EOS, "cva (rho_eos): ideal equation of state required", print_integer=[cv%eos])
    call assert(p%v%v    >  0.0_WP, "cva (rho_eos): p%v > 0 violated", print_real=[p%v%v])
    call assert(temp%v%v >  0.0_WP, "cva (rho_eos): temp%v > 0 violated", print_real=[temp%v%v])
    call assert(size(y)  >= 1,      "cva (rho_eos): size(y) >= 1 violated", print_integer=[size(y)])
    call assert_dimension(p%v%d, temp%v%d)
    call assert_dimension(p%v%d, y(1)%v%d)
    call assert_dimension(y, cv%gas)
    
    n_d = size(p%v%d)
    
    ! <https://en.wikipedia.org/wiki/Molar_mass#Average_molar_mass_of_mixtures>
    ! Why not use the `r_cv` function? That is based around masses, and this is based around mass fractions.
    ! In the `set` method, the masses are not yet available when this method is called.
    call mm%v%init_const(0.0_WP, n_d)
    call y_sum%v%init_const(0.0_WP, n_d)
    do i = 1, size(y)
        call denominator%v%init_const(0.0_WP, n_d)
        do j = 1, size(y)
            denominator = denominator + y(j)*(cv%gas(i)%mm/cv%gas(j)%mm) ! MAYBE: change so that the molar masses have units?
        end do
        call assert(denominator%v%v > 0.0_WP, "cva (rho_eos): denominator is zero", print_real=[denominator%v%v])
        call gas_mm%v%init_const(cv%gas(i)%mm, n_d)
        mm    = mm + y(i)*gas_mm/denominator
        y_sum = y_sum + y(i)
    end do
    
    call assert(is_close(y_sum%v%v, 1.0_WP), "cva (rho_eos): y does not sum to 1", print_real=[y_sum%v%v])
    
    call r_bar_%v%init_const(R_BAR, n_d)
    r_cv = r_bar_ / mm
    
    rho_eos = p / (r_cv * temp)
    
    call assert(rho_eos%v%v > 0.0_WP, "cva (rho_eos): rho_eos%v > 0 violated", print_real=[rho_eos%v%v])
end function rho_eos

pure function p_c(cv)
    ! Estimate critical pressure of a mixture using Kay's rule.
    ! moran_fundamentals_2008 p. 613, eq. 11.97
    
    class(cv_type), intent(in) :: cv
    
    type(si_pressure) :: p_c
    
    integer           :: i, j
    type(unitless)    :: chi ! mole fraction for gas in loop
    type(si_mass)     :: denominator
    type(si_pressure) :: p_ci
    type(unitless)    :: chi_sum
    
    call p_c%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    call chi_sum%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        call denominator%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        do j = 1, size(cv%m)
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm) ! MAYBE: change so that the molar masses have units?
        end do
        call assert(denominator%v%v > 0.0_WP, "cva (p_c): denominator is zero", print_real=[denominator%v%v])
        chi = cv%m(i) / denominator
        call p_ci%v%init_const(cv%gas(i)%p_c, size(cv%m(1)%v%d))
        p_c = p_c + chi*p_ci
        chi_sum = chi_sum + chi
    end do
    
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (p_c): chi does not sum to 1", print_real=[chi_sum%v%v])
end function p_c

pure function y(cv)
    ! mass fractions
    ! Why not use `y` and `chi` to get mass and mole fractions in other procedures here?
    ! Doing so frequently leads to weird run-time errors in gfortran.
    
    class(cv_type), intent(in) :: cv
    
    type(unitless) :: y(size(cv%m))
    
    integer        :: i
    type(si_mass)  :: m_total
    type(unitless) :: y_sum
    
    call y_sum%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    m_total = cv%m_total()
    do i = 1, size(cv%m)
        y(i)  = cv%m(i) / m_total
        y_sum = y_sum + y(i)
        call assert(y(i)%v%v >= 0.0_WP, "cva (y): y >= 0 violated", print_real=[y(i)%v%v])
        call assert(y(i)%v%v <= 1.0_WP, "cva (y): y <= 1 violated", print_real=[y(i)%v%v])
    end do
    
    call assert(is_close(y_sum%v%v, 1.0_WP), "cva (y): y does not sum to 1")
end function y

pure function chi(cv)
    ! mole fractions
    ! Why not use `y` and `chi` to get mass and mole fractions in other procedures here?
    ! Doing so frequently leads to weird run-time errors in gfortran.
    
    class(cv_type), intent(in) :: cv
    
    type(unitless) :: chi(size(cv%m))
    
    integer        :: i, j
    type(si_mass)  :: denominator
    type(unitless) :: chi_sum
    
    call assert(size(cv%m) >= 1, "cva (chi): size(cv%m) >= 1 violated", print_integer=[size(cv%m)])
    
    call chi_sum%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        call denominator%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        do j = 1, size(cv%m)
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm) ! MAYBE: change so that the molar masses have units?
            call assert(denominator%v%v >= 0.0_WP, "cva (chi): denominator >= 0 violated", print_real=[denominator%v%v])
        end do
        call assert(denominator%v%v > 0.0_WP, "cva (chi): denominator is zero", print_real=[denominator%v%v])
        chi(i)  = cv%m(i) / denominator
        chi_sum = chi_sum + chi(i)
        call assert(chi(i)%v%v >= 0.0_WP, "cva (chi): chi >= 0 violated", print_real=[chi(i)%v%v])
        call assert(chi(i)%v%v <= 1.0_WP, "cva (chi): chi <= 1 violated", print_real=[chi(i)%v%v])
    end do
    
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (chi): chi does not sum to 1", print_real=[chi_sum%v%v])
end function chi

pure function r_cv(cv)
    ! Gas constant for a gas *mixture* in a control volume.
    ! Some of the units are intentionally wrong here.
    ! This is done to avoid adding mol to the unit system, which would make compilation much slower.
    
    use gasdata, only: R_BAR
    
    class(cv_type), intent(in) :: cv
    
    type(si_specific_heat) :: r_cv
    
    integer                     :: n_d, i, j
    type(unitless)              :: chi ! mole fraction for gas in loop
    type(si_molar_mass)         :: mm, gas_mm ! molar mass
    type(si_mass)               :: denominator
    type(si_molar_gas_constant) :: r_bar_ ! universal gas constant
    type(unitless)              :: chi_sum
    
    n_d = size(cv%m(1)%v%d)
    
    call assert_mass(cv, "r_cv")
    
    ! <https://en.wikipedia.org/wiki/Molar_mass#Average_molar_mass_of_mixtures>
    call mm%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    call chi_sum%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        call denominator%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        do j = 1, size(cv%m)
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm) ! MAYBE: change so that the molar masses have units?
        end do
        call assert(denominator%v%v > 0.0_WP, "cva (r_cv): denominator is zero", print_real=[denominator%v%v])
        chi = cv%m(i) / denominator
        chi_sum = chi_sum + chi
        call gas_mm%v%init_const(cv%gas(i)%mm, n_d)
        mm  = mm + chi*gas_mm
    end do
    
    call r_bar_%v%init_const(R_BAR, n_d)
    r_cv = r_bar_ / mm
    
    call assert(r_cv%v%v > 0.0_WP, "cva (r_cv): r_cv > 0 violated", print_real=[r_cv%v%v])
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (r_cv): chi does not sum to 1", print_real=[chi_sum%v%v])
end function r_cv

pure function temp_cv(cv)
    use gasdata, only: TEMP_0
    
    class(cv_type), intent(in) :: cv
    
    type(si_temperature) :: temp_cv
    
    integer :: n_d, i
    type(si_energy)          :: e_0
    type(si_specific_energy) :: u_0
    type(si_heat_capacity)   :: heat_capacity
    type(si_temperature)     :: temp_0_
    
    ! Constant specific heats assumed for now. Will improve later.
    
    call assert_mass(cv, "temp_cv")
    
    select case (cv%eos)
        case (IDEAL_EOS)
            n_d = size(cv%m(1)%v%d)
            call e_0%v%init_const(0.0_WP, n_d)
            call heat_capacity%v%init_const(0.0_WP, n_d)
            do i = 1, size(cv%m)
                call u_0%v%init_const(cv%gas(i)%u_0, n_d)
                e_0           = e_0           + cv%m(i)*u_0
                heat_capacity = heat_capacity + cv%m(i)*cv%gas(i)%c_v(n_d)
            end do
            
            call temp_0_%v%init_const(TEMP_0, n_d)
            temp_cv = temp_0_ + (cv%e - e_0) / heat_capacity
        case (CONST_EOS)
            temp_cv = cv%temp_const
        case default
            error stop "cva (temp_cv): invalid cv%eos"
    end select
    
    ! Disabled in favor of `check_sys` temperature check.
    ! If this is active, `check_sys`'s temperature check won't work in debug compilation.
    !call assert(temp_cv%v%v > 0.0_WP, "cva (temp_cv): temp_cv > 0 violated", print_real=[temp_cv%v%v])
end function temp_cv

pure function vol_cv(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_volume) :: vol_cv
    
    call assert(cv%x%v%v > 0.0_WP, "cva (vol_cv): cv%x > 0 violated", print_real=[cv%x%v%v])
    call assert(cv%csa%v%v > 0.0_WP, "cva (vol_cv): cv%csa > 0 violated", print_real=[cv%csa%v%v])
    
    vol_cv = cv%x * cv%csa
    
    call assert(vol_cv%v%v > 0.0_WP, "cva (vol_cv): vol_cv > 0 violated", print_real=[vol_cv%v%v])
end function vol_cv

pure function rho_cv(cv)
    ! Returns simply the mass divided by the volume.
    ! For `CONST_EOS`, this may not be the desired mass density.
    
    class(cv_type), intent(in) :: cv
    
    type(si_mass_density) :: rho_cv
    
    type(si_volume) :: vol
    
    call assert_mass(cv, "rho_cv")
    
    vol = cv%vol()
    call assert(vol%v%v > 0.0_WP, "cva (rho_cv): volume is zero", print_real=[vol%v%v])
    rho_cv = cv%m_total() / vol
    
    select case (cv%eos)
        case (IDEAL_EOS)
            call assert(rho_cv%v%v > 0.0_WP, "cva (rho_cv, IDEAL_EOS): rho_cv > 0 violated", print_real=[rho_cv%v%v])
        case (CONST_EOS)
            call assert(rho_cv%v%v >= 0.0_WP, "cva (rho_cv, CONST_EOS): rho_cv >= 0 violated", print_real=[rho_cv%v%v])
        case default
            error stop "cva (rho_cv): invalid cv%eos"
    end select
end function rho_cv

pure function p_cv(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_pressure) :: p_cv
    
    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    
    rho  = cv%rho()
    temp = cv%temp()
    
    p_cv = cv%p_eos(rho, temp)
    
    call assert(p_cv%v%v > 0.0_WP, "cva (p_cv): p_cv > 0 violated", print_real=[p_cv%v%v])
end function p_cv

pure function u_cv(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_specific_energy) :: u_cv
    
    integer              :: i
    type(si_mass)        :: m_total
    type(si_temperature) :: temp
    
    call assert_mass(cv, "u_cv")
    
    m_total = cv%m_total()
    temp    = cv%temp()
    
    call u_cv%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        u_cv = u_cv + cv%m(i)*cv%gas(i)%u(temp)/m_total
    end do
    
    call assert(u_cv%v%v > 0.0_WP, "cva (u_cv): u_cv > 0 violated", print_real=[u_cv%v%v])
end function u_cv

pure function h_cv(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_specific_energy) :: h_cv
    
    integer              :: i
    type(si_mass)        :: m_total
    type(si_temperature) :: temp
    type(unitless)       :: y
    
    call assert_mass(cv, "h_cv")
    
    m_total = cv%m_total()
    temp    = cv%temp()
    
    call h_cv%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        if (.not. is_close(m_total%v%v, 0.0_WP)) then
            y = cv%m(i)/m_total
        else
            ! TODO: Not sure the derivatives of this should be zero.
            call y%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        end if
        h_cv = h_cv + y*cv%gas(i)%h(temp)
    end do
    
    call assert(h_cv%v%v >= 0.0_WP, "cva (h_cv): h_cv > 0 violated", print_real=[h_cv%v%v])
end function h_cv

pure function gamma_cv(cv, y)
    class(cv_type), intent(in) :: cv
    type(unitless), intent(in) :: y(:)
    
    type(unitless) :: gamma_cv
    
    integer                :: n_d, i
    type(si_specific_heat) :: c_p_cv, c_v_cv
    
    call assert_dimension(y, cv%gas)
    
    n_d = size(y(1)%v%d)
    
    call c_p_cv%v%init_const(0.0_WP, n_d)
    call c_v_cv%v%init_const(0.0_WP, n_d)
    do i = 1, size(y)
        c_p_cv = c_p_cv + y(i)*cv%gas(i)%c_p(n_d)
        c_v_cv = c_v_cv + y(i)*cv%gas(i)%c_v(n_d)
    end do
    
    gamma_cv = c_p_cv / c_v_cv
    
    call assert(gamma_cv%v%v > 1.0_WP, "cva (gamma_cv): gamma_cv > 1 violated", print_real=[gamma_cv%v%v])
end function gamma_cv

pure subroutine set(cv, x, x_dot, y, p, temp_atm, label, csa, rm_p, p_fs, p_fd, k, l_pre, gas, &
                            i_cv_mirror, x_stop, isentropic_filling, p_atm, eos, type, m_spring, constant_friction)
    class(cv_type), intent(in out) :: cv
    
    ! time varying
    type(si_length), intent(in)      :: x        ! location of projectile/plunger
    type(si_velocity), intent(in)    :: x_dot    ! velocity of projectile/plunger
    type(unitless), intent(in)       :: y(:)     ! mass fractions of each gas
    type(si_pressure), intent(in)    :: p        ! pressure
    
    ! constant
    type(si_temperature), intent(in)  :: temp_atm    ! atmospheric temperature
    character(len=*), intent(in)      :: label       ! human-readable label for control volume
    type(si_area), intent(in)         :: csa         ! cross-sectional area
    type(si_inverse_mass), intent(in) :: rm_p        ! reciprocal mass of projectile/plunger
    type(si_pressure), intent(in)     :: p_fs, p_fd  ! static and dynamic friction pressure
    type(si_stiffness), intent(in)    :: k           ! stiffness of spring attached to plunger
    type(si_length), intent(in)       :: l_pre       ! spring precompression length
    type(gas_type), intent(in)        :: gas(:)      ! gas data
    integer, intent(in)               :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    
    type(si_length), intent(in), optional   :: x_stop   ! `x` location where simulation will stop
    logical, intent(in), optional           :: isentropic_filling
    type(si_pressure), intent(in), optional :: p_atm    ! atmospheric pressure (only requried if `isentropic_filling = .true.`
    integer, intent(in), optional           :: eos      ! equation of state to use
    integer, intent(in), optional           :: type     ! type of CV to use
    type(si_mass), intent(in), optional     :: m_spring ! mass of spring
    logical, intent(in), optional           :: constant_friction
    
    integer              :: i, n_d
    type(si_temperature) :: temp
    type(si_mass)        :: m_total
    type(unitless)       :: gamma_cv, y_sum
    logical              :: isentropic_filling_
    
    n_d = size(x%v%d)
    
    cv%x     = x
    cv%x_dot = x_dot
    call cv%e_f%v%init_const(0.0_WP, n_d)
    ! `p` and `temp` will be handled below
    
    cv%label       = label
    cv%csa         = csa
    cv%rm_p        = rm_p
    cv%p_fs        = p_fs
    cv%p_fd        = p_fd
    cv%k           = k
    cv%l_pre       = l_pre
    cv%gas         = gas
    cv%i_cv_mirror = i_cv_mirror
    
    if (present(x_stop)) then
        cv%x_stop = x_stop
        
        call assert(cv%x%v%v < cv%x_stop%v%v, "cva (set): x >= x_stop will cause immediate termination of run", &
                    print_real=[cv%x%v%v, cv%x_stop%v%v])
    else
        call cv%x_stop%v%init_const(X_STOP_DEFAULT, n_d)
    end if
    
    if (present(isentropic_filling)) then
        isentropic_filling_ = isentropic_filling
    else
        isentropic_filling_ = .false.
    end if
    
    if (present(constant_friction)) then
        cv%constant_friction = constant_friction
        
        if (cv%constant_friction) then
            call assert(is_close(p_fs%v%v, p_fd%v%v), &
                    "cva (set): constant_friction = .true. requires that p_fs = p_fd as otherwise p_f would not be constant", &
                    print_real=[p_fs%v%v, p_fd%v%v])
        end if
    else
        cv%constant_friction = .false.
    end if
    
    if (present(eos)) then
        cv%eos = eos
    else
        cv%eos = IDEAL_EOS
    end if
    
    if (present(type)) then
        cv%type = type
    else
        cv%type = NORMAL_CV_TYPE
    end if
    
    if (present(m_spring)) then
        cv%m_spring = m_spring
        
        if (cv%k%v%v > 0.0_WP) then
            call assert(cv%m_spring%v%v > 0.0_WP, "cva (set): m_spring must be set to > 0 if k > 0", &
                            print_real=[cv%m_spring%v%v])
        else
            call assert(is_close(cv%m_spring%v%v, 0.0_WP), "cva (set): m_spring must be set to 0 if k == 0", &
                            print_real=[cv%m_spring%v%v])
        end if
    else
        call cv%m_spring%v%init_const(0.0_WP, n_d)
    end if
    
    call assert(cv%x%v%v            >  0.0_WP, "cva (set): x > 0 violated", print_real=[cv%x%v%v])
    call assert(p%v%v               >  0.0_WP, "cva (set): p > 0 violated", print_real=[p%v%v])
    call assert(temp_atm%v%v        >  0.0_WP, "cva (set): temp_atm > 0 violated", print_real=[temp_atm%v%v])
    call assert(len(trim(cv%label)) >       0, "cva (set): len(label) > 0 violated", print_integer=[len(trim(cv%label))])
    call assert(cv%csa%v%v          >  0.0_WP, "cva (set): csa > 0 violated", print_real=[cv%csa%v%v])
    call assert(cv%p_fs%v%v         >= 0.0_WP, "cva (set): p_fs >= 0 violated", print_real=[cv%p_fs%v%v])
    call assert(cv%p_fd%v%v         >= 0.0_WP, "cva (set): p_fd >= 0 violated", print_real=[cv%p_fd%v%v])
    call assert(cv%k%v%v            >= 0.0_WP, "cva (set): k >= 0 violated", print_real=[cv%k%v%v])
    call assert(cv%i_cv_mirror      >= 0,      "cva (set): i_cv_mirror >= 0 violated", print_integer=[cv%i_cv_mirror])
    call assert(cv%m_spring%v%v     >= 0.0_WP, "cva (set): m_spring >= 0 violated", print_real=[cv%m_spring%v%v])
    
    call assert((cv%eos  >= 1) .and. (cv%eos <= MAX_EOS),      "cva (set): invalid EOS", print_integer=[cv%eos])
    call assert((cv%type >= 1) .and. (cv%type <= MAX_CV_TYPE), "cva (set): invalid control volume type", print_integer=[cv%type])
    
    call assert_dimension(y, cv%gas)
    allocate(cv%m(size(y)))
    call y_sum%v%init_const(0.0_WP, n_d)
    call m_total%v%init_const(1.0_WP, n_d)
    do i = 1, size(y)
        call assert(gas(i)%gamma > 1.0_WP, "cva (set): gas%gamma > 1 violated", print_real=[gas(i)%gamma])
        call assert(gas(i)%mm    > 0.0_WP, "cva (set): gas%mm > 0 violated", print_real=[gas(i)%mm])
        
        ! to catch using g/mol by mistake
        call assert(gas(i)%mm    < 0.1_WP, "cva (set): gas%mm < 0.1 violated", print_real=[gas(i)%mm])
        
        call assert(gas(i)%p_c   > 0.0_WP, "cva (set): gas%p_c > 0 violated", print_real=[gas(i)%p_c])
        
        y_sum = y_sum + y(i)
    end do
    
    call assert(is_close(y_sum%v%v, 1.0_WP), "cva (set): mass fractions do not sum to 1", print_real=[y_sum%v%v])
    
    ! Get correct temperature depending on how the chamber is filled.
    if (isentropic_filling_) then
        call assert(present(p_atm), "cva (set): isentropic_filling = .true. requires p_atm")
        call assert(p_atm%v%v > 0.0_WP, "cva (set): p_atm > 0 required for isentropic_filling", print_real=[p_atm%v%v])
        gamma_cv = cv%gamma(y)
        temp = temp_atm * ((p / p_atm)**((gamma_cv - 1.0_WP)/gamma_cv))
    else
        call assert(.not. present(p_atm), &
                        "cva (set): p_atm only affects isentropic_filling so it should not be set otherwise")
        
        ! isothermal
        temp = temp_atm
    end if
    
    m_total = cv%vol() * cv%rho_eos(p, temp, y)
    
    ! Now set `m` and `e`
    call cv%e%v%init_const(0.0_WP, n_d)
    do i = 1, size(y)
        cv%m(i) = y(i)*m_total
        cv%e    = cv%e + cv%m(i)*cv%gas(i)%u(temp)
    end do
    
    call assert_mass(cv, "set")
    call assert_dimension(cv%x%v%d, cv%x_dot%v%d)
    call assert_dimension(cv%x%v%d, cv%e%v%d)
    
    if (cv%eos == IDEAL_EOS) then
        ! This is checked here and not in `rho_eos` as the masses are not defined when `rho_eos` is called.
        call assert(p < cv%p_c(), "cva (set): ideal gas law validity is questionable", print_real=[p%v%v])
    end if
end subroutine set

pure subroutine set_const(cv, label, csa, p_const, temp_const, gas, y_const, i_cv_mirror, type, x_dot)
    class(cv_type), intent(in out) :: cv
    
    character(len=*), intent(in)     :: label       ! human-readable label for control volume
    type(si_area), intent(in)        :: csa         ! cross-sectional area
    type(si_pressure), intent(in)    :: p_const     ! pressure
    type(si_temperature), intent(in) :: temp_const  ! temperature
    type(gas_type), intent(in)       :: gas(:)      ! gas data
    type(unitless), intent(in)       :: y_const(:)  ! mass fractions of each gas
    integer, intent(in)              :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    
    integer, intent(in), optional           :: type  ! type of CV to use
    type(si_velocity), intent(in), optional :: x_dot ! initial velocity of projectile/plunger
    
    integer :: n_d, n_gas, k_gas
    
    call assert(len(trim(cv%label)) > 0, "cva (set_const): len(label) > 0 violated", print_integer=[len(trim(cv%label))])
    
    n_d   = size(p_const%v%d)
    n_gas = size(gas)
    
    allocate(cv%m(n_gas))
    do k_gas = 1, n_gas
        call cv%m(k_gas)%v%init_const(0.0_WP, n_d)
        call assert(y_const(k_gas)%v%v >= 0.0_WP, "cva (set_const): y_const must be >= 0")
        call assert(y_const(k_gas)%v%v <= 1.0_WP, "cva (set_const): y_const must be <= 1")
    end do
    
    call cv%x%v%init_const(1.0_WP, n_d)
    call cv%e%v%init_const(0.0_WP, n_d)
    call cv%e_f%v%init_const(0.0_WP, n_d)
    call cv%rm_p%v%init_const(0.0_WP, n_d)
    call cv%p_fs%v%init_const(0.0_WP, n_d)
    call cv%p_fd%v%init_const(0.0_WP, n_d)
    call cv%k%v%init_const(0.0_WP, n_d)
    call cv%l_pre%v%init_const(0.0_WP, n_d)
    call cv%x_stop%v%init_const(X_STOP_DEFAULT, n_d)
    
    call assert(cv%x%v%v < cv%x_stop%v%v, "cva (set_const): x >= x_stop will cause immediate termination of run", &
                    print_real=[cv%x%v%v, cv%x_stop%v%v])
    
    cv%label       = label
    cv%csa         = csa
    cv%eos         = CONST_EOS
    cv%gas         = gas
    cv%i_cv_mirror = i_cv_mirror
    cv%p_const     = p_const
    cv%temp_const  = temp_const
    cv%y_const     = y_const
    
    call assert(cv%p_const%v%v    > 0.0_WP, "cva (set_const): p_const > 0 violated", print_real=[cv%p_const%v%v])
    call assert(cv%temp_const%v%v > 0.0_WP, "cva (set_const): temp_const > 0 violated", print_real=[cv%temp_const%v%v])
    
    if (present(type)) then
        cv%type = type
    else
        cv%type = MIRROR_CV_TYPE
    end if
    
    if (present(x_dot)) then
        cv%x_dot = x_dot
    else
        call cv%x_dot%v%init_const(0.0_WP, n_d)
    end if
    
    select case (cv%type)
        case (NORMAL_CV_TYPE)
            call assert(cv%i_cv_mirror >= 0, "cva (set_const): i_cv_mirror >= 0 violated", print_integer=[cv%i_cv_mirror])
        case (MIRROR_CV_TYPE)
            call assert(cv%i_cv_mirror >= 1, "cva (set_const, MIRROR_CV_TYPE): i_cv_mirror >= 1 violated", &
                            print_integer=[cv%i_cv_mirror])
        case default
            error stop "cva (set_const): invalid cv%type"
    end select
end subroutine set_const

!tripwire$ begin 28F1DB4A Update `\secref{friction}` of theory.tex when changing `d_x_dot_d_t_normal` if necessary.
pure function p_f(cv, p_fe)
    ! Returns pressure of friction.
    
    class(cv_type), intent(in)    :: cv
    type(si_pressure), intent(in) :: p_fe ! equilibrium pressure
    
    type(si_pressure) :: p_f
    
    type(si_velocity) :: v_scale
    
    call v_scale%v%init_const(0.1_WP, size(cv%x%v%d))
    
    call assert(cv%p_fs%v%v >= 0.0_WP, "cva (p_f): cv%p_fs%v > 0 violated", print_real=[cv%p_fs%v%v])
    call assert(cv%p_fd%v%v >= 0.0_WP, "cva (p_f): cv%p_fd%v > 0 violated", print_real=[cv%p_fd%v%v])
    
    if (cv%constant_friction) then
        p_f = cv%p_fs
    else
        p_f = p_f0(cv, p_fe) + (cv%p_fd - tanh(cv%x_dot/v_scale)*p_f0(cv, p_fe))*tanh(cv%x_dot/v_scale)
    end if
    
    ! The 1.1 factor was added as I guess the inequality with a factor of 1.0 isn't guaranteed for numerical reasons?
    ! But even that is sometimes violated?
    call assert(abs(p_f%v%v) <= 1.1_WP*max(cv%p_fs%v%v, cv%p_fd%v%v), "cva (p_f): abs(p_f) <= 1.1*max(p_fs, p_fd) violated", &
                    print_real=[p_f%v%v, cv%p_fs%v%v, cv%p_fd%v%v])
end function p_f

pure function p_f0(cv, p_fe)
    ! Returns actual static pressure of friction.
    ! `p_fs` is the *maximum* static pressure of friction.
    
    class(cv_type), intent(in)    :: cv
    type(si_pressure), intent(in) :: p_fe ! equilibrium pressure
    
    type(si_pressure) :: p_f0, p_s
    
    call assert(cv%p_fs%v%v >= 0.0_WP, "cva (p_f0): cv%p_fs%v > 0 violated", print_real=[cv%p_fs%v%v])
    
    p_s = 0.1_WP*cv%p_fs ! TODO: make a function of `dt`
    call assert(p_s <= cv%p_fs, "cva (p_f0): p_s <= p_fs violated", print_real=[p_s%v%v, cv%p_fs%v%v])
    
    if (p_fe <= -p_s) then
        p_f0 = -p_f0_high(p_fe, cv%p_fs, p_s)
        
        call assert(p_f0 <= -p_s, "cva (p_f0), first branch: p_f0 <= -p_s violated", print_real=[p_f0%v%v, p_s%v%v])
    else if (p_fe <= p_s) then
        p_f0 = p_fe
    else
        p_f0 = p_f0_high(p_fe, cv%p_fs, p_s)
        
        call assert(p_f0 >= p_s, "cva (p_f0), third branch: p_f0 >= p_s violated", print_real=[p_f0%v%v, p_s%v%v])
    end if
    
    call assert(p_f0 >= -cv%p_fs, "cva (p_f0): p_f0 >= -p_fs violated", print_real=[p_f0%v%v, cv%p_fs%v%v])
    call assert(p_f0 <= cv%p_fs, "cva (p_f0): p_f0 <= p_fs violated", print_real=[p_f0%v%v, cv%p_fs%v%v])
    
    contains
    
    pure function p_f0_high(p_fe, p_fs, p_s)
        type(si_pressure), intent(in) :: p_fe, p_fs, p_s
        type(si_pressure) :: p_f0_high
        
        if (is_close(p_fs%v%v, 0.0_WP)) then
            p_f0_high = p_fs
        else
            p_f0_high = p_fs * tanh((abs(p_fe) - p_s)/(p_fs - p_s) + atanh(p_s/p_fs))
        end if
    end function p_f0_high
end function p_f0
!tripwire$ end

pure function d_x_d_t(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_velocity) :: d_x_d_t
    
    call assert(sys%cv(i_cv)%i_cv_mirror >= 0, "cva (d_x_d_t): i_cv_mirror must be a positive integer or zero", &
                print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
    call assert(sys%cv(i_cv)%i_cv_mirror <= size(sys%cv), "cva (d_x_d_t): i_cv_mirror can't be larger than the cv array", &
                print_integer=[sys%cv(i_cv)%i_cv_mirror, size(sys%cv), i_cv])
    call assert(sys%cv(i_cv)%i_cv_mirror /= i_cv, "cva (d_x_d_t): i_cv_mirror can not equal i_cv", &
                    print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
    
    select case (sys%cv(i_cv)%type)
        case (NORMAL_CV_TYPE)
            d_x_d_t = sys%cv(i_cv)%x_dot
        case (MIRROR_CV_TYPE)
            call assert(sys%cv(i_cv)%i_cv_mirror >= 1, "cva (d_x_d_t): i_cv_mirror must be defined for a mirror CV", &
                            print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
            call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == NORMAL_CV_TYPE, &
                            "cva (d_x_d_t): mirror CV is not a NORMAL_CV_TYPE", &
                            print_integer=[sys%cv(sys%cv(i_cv)%i_cv_mirror)%type, NORMAL_CV_TYPE, i_cv])
            d_x_d_t = -sys%cv(sys%cv(i_cv)%i_cv_mirror)%x_dot
        case default
            error stop "cva (d_x_d_t): invalid cv%type"
    end select
end function d_x_d_t

pure function d_x_dot_d_t(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_acceleration) :: d_x_dot_d_t
    
    type(si_area)     :: csa_i_mirror
    character(len=3)  :: i_cv_string, i_cv_mirror_string
    
    call assert(sys%cv(i_cv)%i_cv_mirror >= 0, "cva (d_x_dot_d_t): i_cv_mirror must be a positive integer or zero", &
                print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
    call assert(sys%cv(i_cv)%i_cv_mirror <= size(sys%cv), "cva (d_x_dot_d_t): i_cv_mirror can't be larger than the cv array", &
                print_integer=[sys%cv(i_cv)%i_cv_mirror, size(sys%cv), i_cv])
    call assert(sys%cv(i_cv)%i_cv_mirror /= i_cv, "cva (d_x_dot_d_t): i_cv_mirror can not equal i_cv", &
                print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
    
    if (sys%cv(i_cv)%i_cv_mirror >= 1) then
        csa_i_mirror = sys%cv(sys%cv(i_cv)%i_cv_mirror)%csa
        write(unit=i_cv_string, fmt="(i0)") i_cv
        write(unit=i_cv_mirror_string, fmt="(i0)") sys%cv(i_cv)%i_cv_mirror
        call assert(is_close(sys%cv(i_cv)%csa%v%v, csa_i_mirror%v%v), &
                        "cva (d_x_dot_d_t): CSA for mirror CV (" // trim(i_cv_mirror_string) &
                        // ") assumed equal to CSA for current CV (" // trim(i_cv_string) // ")", &
                        print_real=[sys%cv(i_cv)%csa%v%v, csa_i_mirror%v%v], print_integer=[i_cv])
    end if
    
    select case (sys%cv(i_cv)%type)
        case (NORMAL_CV_TYPE)
            d_x_dot_d_t = d_x_dot_d_t_normal(sys, i_cv)
        case (MIRROR_CV_TYPE)
            call assert(sys%cv(i_cv)%i_cv_mirror >= 1, "cva (d_x_dot_d_t): i_cv_mirror must be defined for a mirror CV", &
                            print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
            call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == NORMAL_CV_TYPE, &
                            "cva (d_x_dot_d_t): mirror CV is not a NORMAL_CV_TYPE", &
                            print_integer=[sys%cv(sys%cv(i_cv)%i_cv_mirror)%type, NORMAL_CV_TYPE, i_cv])
            call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%i_cv_mirror == i_cv, &
                            "cva (d_x_dot_d_t): i_cv_mirror for the mirror CV is not i_cv (not matching)", &
                            print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
            d_x_dot_d_t = -d_x_dot_d_t_normal(sys, sys%cv(i_cv)%i_cv_mirror)
        case default
            error stop "cva (d_x_dot_d_t): invalid cv%type"
    end select
end function d_x_dot_d_t

!tripwire$ begin 0D2F5ED8 Update `\secref{plunger-impact}` of theory.tex when changing `d_x_dot_d_t_normal` if necessary.
pure function d_x_dot_d_t_normal(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_acceleration) :: d_x_dot_d_t_normal
    
    type(si_pressure)     :: p_fe ! friction pressure at equilibrium ($\partial \dot{x}/\partial t = 0$)
    type(si_pressure)     :: p_mirror
    type(si_inverse_mass) :: r_mp_eff ! effective mass of projectile/plunger
    
    call assert(sys%cv(i_cv)%csa%v%v > 0.0_WP, "cva (d_x_dot_d_t_normal): cv%csa > 0 violated", &
                    print_real=[sys%cv(i_cv)%csa%v%v], print_integer=[i_cv])
    call assert(sys%cv(i_cv)%type == NORMAL_CV_TYPE, "cva (d_x_dot_d_t_normal): CV needs to be NORMAL_CV_TYPE", &
                        print_integer=[sys%cv(i_cv)%type, NORMAL_CV_TYPE, i_cv])
    
    if (sys%cv(i_cv)%i_cv_mirror >= 1) then
        p_mirror = sys%cv(sys%cv(i_cv)%i_cv_mirror)%p()
        
        call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == MIRROR_CV_TYPE, &
                        "cva (d_x_dot_d_t_normal): mirror CV not MIRROR_CV_TYPE", &
                        print_integer=[sys%cv(sys%cv(i_cv)%i_cv_mirror)%type, MIRROR_CV_TYPE, i_cv])
    else
        ! If `i_cv_mirror == 0` then there is no mirror CV.
        call p_mirror%v%init_const(0.0_WP, size(sys%cv(i_cv)%csa%v%d))
    end if
    
    p_fe = sys%cv(i_cv)%p() - p_mirror - (sys%cv(i_cv)%k/sys%cv(i_cv)%csa)*(sys%cv(i_cv)%x + sys%cv(i_cv)%l_pre)
    
    ! This calculates the effective (inverse) mass of the projectile/plunger factoring in the spring mass.
    ! ruby_equivalent_2000
    r_mp_eff = sys%cv(i_cv)%rm_p / (1.0_WP + C_MS * sys%cv(i_cv)%m_spring * sys%cv(i_cv)%rm_p)
    
    d_x_dot_d_t_normal = sys%cv(i_cv)%csa*r_mp_eff*(sys%cv(i_cv)%p() - p_mirror - sys%cv(i_cv)%p_f(p_fe)) &
                            - sys%cv(i_cv)%k*r_mp_eff*(sys%cv(i_cv)%x + sys%cv(i_cv)%l_pre)
end function d_x_dot_d_t_normal
!tripwire$ end

pure function d_m_k_d_t(sys, m_dot, k_gas, i_cv)
    type(cv_system_type), intent(in)    :: sys
    type(si_mass_flow_rate), intent(in) :: m_dot(:, :)
    integer, intent(in)                 :: k_gas, i_cv
    
    type(si_mass_flow_rate) :: d_m_k_d_t
    
    integer :: n_d, n_cv, j_cv
    type(si_mass)  :: m_total
    type(unitless) :: y_k_i, y_k_j
    
    call assert(size(m_dot, 1) == size(m_dot, 2), "cva (d_m_k_d_t): m_dots must be square", &
                    print_integer=[size(m_dot, 1), size(m_dot, 2)])
    call assert_dimension(m_dot(1, 1)%v%d, sys%cv(i_cv)%x%v%d)
    
    n_d = size(sys%cv(i_cv)%x%v%d)
    call d_m_k_d_t%v%init_const(0.0_WP, n_d)
    
    n_cv = size(m_dot, 1)
    if (sys%cv(i_cv)%eos == CONST_EOS) then
        y_k_i = sys%cv(i_cv)%y_const(k_gas)
    else
        m_total = sys%cv(i_cv)%m_total()
        if (.not. is_close(m_total%v%v, 0.0_WP)) then
            y_k_i = sys%cv(i_cv)%m(k_gas) / m_total
        else
            call y_k_i%v%init_const(0.0_WP, n_d) ! TODO: Not sure the derivatives of this should be zero.
        end if
    end if
    
    do j_cv = 1, n_cv
        call assert(is_close(m_dot(j_cv, j_cv)%v%v, 0.0_WP), "cva (d_m_k_d_t): mass can not flow from self to self", &
                        print_real=[m_dot(j_cv, j_cv)%v%v])
        
        if (sys%cv(j_cv)%eos == CONST_EOS) then
            y_k_j = sys%cv(j_cv)%y_const(k_gas)
        else
            m_total = sys%cv(j_cv)%m_total()
            if (.not. is_close(m_total%v%v, 0.0_WP)) then
                y_k_j = sys%cv(j_cv)%m(k_gas) / m_total
            else
                call y_k_j%v%init_const(0.0_WP, n_d) ! TODO: Not sure the derivatives of this should be zero.
            end if
        end if
        
        d_m_k_d_t = d_m_k_d_t + y_k_j*m_dot(j_cv, i_cv) - y_k_i*m_dot(i_cv, j_cv)
    end do
end function d_m_k_d_t

pure function d_e_d_t(cv, h_dot, i_cv)
    type(cv_type), intent(in)             :: cv
    type(si_energy_flow_rate), intent(in) :: h_dot(:, :)
    integer, intent(in)                   :: i_cv
    
    type(si_energy_flow_rate) :: d_e_d_t
    
    integer :: n_cv, j_cv
    
    call assert(size(h_dot, 1) == size(h_dot, 2), "cva (d_e_d_t): h_dots must be square", &
                    print_integer=[size(h_dot, 1), size(h_dot, 2)])
    call assert_dimension(h_dot(1, 1)%v%d, cv%x%v%d)
    call assert(cv%csa%v%v > 0.0_WP, "cva (d_e_d_t): cv%csa > 0 violated", print_real=[cv%csa%v%v])
    
    d_e_d_t = -cv%p() * cv%csa * cv%x_dot
    
    n_cv = size(h_dot, 1)
    do j_cv = 1, n_cv
        call assert(is_close(h_dot(j_cv, j_cv)%v%v, 0.0_WP), "cva (d_e_d_t): energy can not flow from self to self", &
                        print_real=[h_dot(j_cv, j_cv)%v%v])
        d_e_d_t = d_e_d_t + h_dot(j_cv, i_cv) - h_dot(i_cv, j_cv)
    end do
end function d_e_d_t

pure function d_e_f_d_t(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_energy_flow_rate) :: d_e_f_d_t
    
    type(si_pressure) :: p_fe ! friction pressure at equilibrium ($\partial \dot{x}/\partial t = 0$)
    type(si_pressure) :: p_mirror
    
    call assert(sys%cv(i_cv)%i_cv_mirror >= 0, "cva (d_e_f_d_t): i_cv_mirror must be a positive integer or zero", &
                    print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
    call assert(sys%cv(i_cv)%i_cv_mirror <= size(sys%cv), "cva (d_e_f_d_t): i_cv_mirror can't be larger than the cv array", &
                    print_integer=[sys%cv(i_cv)%i_cv_mirror, size(sys%cv), i_cv])
    call assert(sys%cv(i_cv)%i_cv_mirror /= i_cv, "cva (d_e_f_d_t): i_cv_mirror can not equal i_cv", &
                    print_integer=[sys%cv(i_cv)%i_cv_mirror, i_cv])
    
    call assert(sys%cv(i_cv)%csa%v%v > 0.0_WP, "cva (d_x_dot_d_t_normal): cv%csa > 0 violated", &
                    print_real=[sys%cv(i_cv)%csa%v%v])
    
    select case (sys%cv(i_cv)%type)
        case (NORMAL_CV_TYPE)
            if (sys%cv(i_cv)%i_cv_mirror >= 1) then
                p_mirror = sys%cv(sys%cv(i_cv)%i_cv_mirror)%p()
                
                call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == MIRROR_CV_TYPE, &
                                "cva (d_x_dot_d_t_normal): mirror CV not MIRROR_CV_TYPE", &
                                print_integer=[sys%cv(sys%cv(i_cv)%i_cv_mirror)%type, MIRROR_CV_TYPE])
            else
                ! If `i_cv_mirror == 0` then there is no mirror CV.
                call p_mirror%v%init_const(0.0_WP, size(sys%cv(i_cv)%csa%v%d))
            end if
            
            p_fe = sys%cv(i_cv)%p() - p_mirror - (sys%cv(i_cv)%k/sys%cv(i_cv)%csa)*(sys%cv(i_cv)%x + sys%cv(i_cv)%l_pre)
            
            d_e_f_d_t = sys%cv(i_cv)%p_f(p_fe) * sys%cv(i_cv)%csa * sys%cv(i_cv)%x_dot
        case (MIRROR_CV_TYPE)
            call d_e_f_d_t%v%init_const(0.0_WP, size(sys%cv(i_cv)%csa%v%d))
            call assert(is_close(sys%cv(i_cv)%e_f%v%v, 0.0_WP), "cva (d_e_f_d_t): e_f must be zero for MIRROR_CV_TYPE", &
                            print_real=[sys%cv(i_cv)%e_f%v%v])
        case default
            error stop "cva (d_e_f_d_t): invalid cv%type"
    end select
end function d_e_f_d_t

pure subroutine assert_mass(cv, procedure_name)
    ! Why not make this a type-bound operator?
    ! That would make my assertion counting Python program not count these.
    
    type(cv_type), intent(in)    :: cv
    character(len=*), intent(in) :: procedure_name
    
    integer       :: i
    type(si_mass) :: m_total
    
    call assert(len(trim(procedure_name)) > 0, "cva (assert_mass): procedure name should not be empty")
    
    do i = 1, size(cv%m)
        call assert(cv%m(i)%v%v >= 0.0_WP, "cva (" // trim(procedure_name) // "): cv%m >= 0 violated", &
                        print_real=[cv%m(i)%v%v], print_integer=[i])
    end do
    
    m_total = cv%m_total()
    select case (cv%eos)
        case (IDEAL_EOS)
            call assert(m_total%v%v > 0.0_WP, "cva (" // trim(procedure_name) // ", IDEAL_EOS): cv%m_total > 0 violated", &
                            print_real=[m_total%v%v])
        case (CONST_EOS)
            call assert(m_total%v%v >= 0.0_WP, "cva (" // trim(procedure_name) // ", CONST_EOS): cv%m_total >= 0 violated", &
                            print_real=[m_total%v%v])
        case default
            error stop "cva (assert_mass): invalid cv%eos"
    end select
    
    ! I would also check that `cv%e` is positive here, but...
    ! Strictly speaking, the people making the thermodynamic tables might not have made internal energy always positive.
    ! Given that, violating `cv%e > 0` might be okay.
    ! So I'm making the assertion based on temperature as that should always be positive.
    ! I can't put this assertion here because it would make this a recursive subroutine.
    ! The assertion in `temp_cv` is sufficient.
end subroutine assert_mass

!!!!!!!!!!!!!!!!!!!!!!!!!!
! methods for `con_type` !
!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function smooth_min(x, y)
    ! Exponential from <https://iquilezles.org/articles/smin/>.
    ! Also see: <https://en.wikipedia.org/wiki/Smooth_maximum>
    
    type(unitless), intent(in) :: x, y
    
    type(unitless) :: smooth_min
    
    type(unitless) :: k
    
    call assert_dimension(x%v%d, y%v%d)
    
    call k%v%init_const(0.01_WP, size(x%v%d))
    
    smooth_min = -k*log(exp(-x/k) + exp(-y/k))
end function smooth_min

pure function f_m_dot(p_r, b)
    ! See beater_pneumatic_2007 eq. 5.4
    ! This is a replacement for the ((p_2/p_1 - b)/(1-b))**2 term, smoothly going between the various cases.
    
    type(unitless), intent(in) :: p_r, b
    
    type(unitless) :: f_m_dot
    
    type(unitless) :: p_rs, p_rl_ ! scales used to make function differentiable
    
    call assert(p_r%v%v >= 0.0_WP, "cva (f_m_dot): p_r >= 0 violated", print_real=[p_r%v%v])
    call assert(p_r%v%v <= 1.0_WP, "cva (f_m_dot): p_r <= 1 violated", print_real=[p_r%v%v])
    
    call assert_dimension(p_r%v%d, b%v%d)
    
    call p_rs%v%init_const(0.01_WP, size(p_r%v%d)) ! TODO: make `p_rs` a function of `dt`
    call p_rl_%v%init_const(P_RL, size(p_r%v%d)) ! based on first part of beater_pneumatic_2007 eq. 5.4
    
    f_m_dot = square((smooth_min(p_r, p_rl_) - b) / (1.0_WP - b)) &
                * 0.5_WP * (1.0_WP + tanh((p_r - b) / p_rs))
end function f_m_dot

pure function g_m_dot(p_r)
    ! This is a multiplier on the `p_2` term, smoothly going between the various cases.
    
    type(unitless), intent(in) :: p_r
    
    type(unitless) :: g_m_dot
    
    call assert(p_r%v%v >= 0.0_WP, "cva (g_m_dot): p_r >= 0 violated", print_real=[p_r%v%v])
    call assert(p_r%v%v <= 1.0_WP, "cva (g_m_dot): p_r <= 1 violated", print_real=[p_r%v%v])
    
    if (p_r%v%v < P_RL) then
        call g_m_dot%v%init_const(0.0_WP, size(p_r%v%d))
    else
        g_m_dot = 2.0_WP*((1.0_WP - p_r) / (1.0_WP - P_RL))*((1.0_WP - p_r) / (1.0_WP - P_RL))*((1.0_WP - p_r) / (1.0_WP - P_RL)) &
                    - 3.0_WP*((1.0_WP - p_r) / (1.0_WP - P_RL))*((1.0_WP - p_r) / (1.0_WP - P_RL)) + 1.0_WP
    end if
    
    call assert(g_m_dot%v%v >= 0.0_WP, "cva (g_m_dot): g_m_dot >= 0 violated", print_real=[g_m_dot%v%v])
    call assert(g_m_dot%v%v <= 1.0_WP, "cva (g_m_dot): g_m_dot <= 1 violated", print_real=[g_m_dot%v%v])
end function g_m_dot

pure function m_dot(con, cv_from, cv_to)
    ! Modified con flow rate model from beater_pneumatic_2007 ch. 5.
    ! Modified to be differentiable.
    
    class(con_type), intent(in) :: con
    type(cv_type), intent(in)   :: cv_from, cv_to
    
    type(si_mass_flow_rate) :: m_dot
    
    integer        :: n_d
    type(unitless) :: p_r
    
    call assert_dimension(cv_from%x%v%d, cv_to%x%v%d)
    
    if (con%active .and. (cv_from%p() >= cv_to%p())) then
        call assert_dimension(con%a_e%v%d, con%b%v%d)
        call assert_dimension(cv_from%x%v%d, con%a_e%v%d)
        
        p_r = cv_to%p() / cv_from%p()
        call assert(p_r%v%v >= 0.0_WP, "cva (m_dot): p_r >= 0 violated", print_real=[p_r%v%v])
        call assert(p_r%v%v <= 1.0_WP, "cva (m_dot): p_r <= 1 violated", print_real=[p_r%v%v])
        
        m_dot = con%a_e * (cv_from%p() - g_m_dot(p_r) * cv_to%p()) &
                    * sqrt((1.0_WP - con%b) / (cv_from%r() * cv_from%temp())) &
                    * sqrt(1.0_WP - f_m_dot(p_r, con%b))
    else
        n_d = size(cv_from%x%v%d)
        call m_dot%v%init_const(0.0_WP, n_d)
    end if
end function m_dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! methods for control volume systems !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure subroutine calculate_flows(sys, m_dot, h_dot)
    class(cv_system_type), intent(in)                   :: sys
    type(si_mass_flow_rate), allocatable, intent(out)   :: m_dot(:, :)
    type(si_energy_flow_rate), allocatable, intent(out) :: h_dot(:, :)
    
    integer :: n_cv, i_from_cv, i_to_cv
    
    call assert(size(sys%cv) == size(sys%con, 1), "cva (calculate_flows): inconsistent sys%cv and sys%con sizes", &
                    print_integer=[size(sys%cv), size(sys%con, 1)])
    
    n_cv = size(sys%cv)
    call assert(n_cv > 1, "cva (calculate_flows): there needs to be at least 2 control volumes to have flows", &
                    print_integer=[n_cv])
    
    allocate(m_dot(n_cv, n_cv))
    allocate(h_dot(n_cv, n_cv))
    do i_from_cv = 1, n_cv
        do i_to_cv = 1, n_cv
            if (i_from_cv == i_to_cv) call assert(.not. sys%con(i_from_cv, i_to_cv)%active, &
                                                    "cva (calculate_flows): can't flow from self to self", &
                                                    print_integer=[i_from_cv, i_to_cv])
            
            m_dot(i_from_cv, i_to_cv) = sys%con(i_from_cv, i_to_cv)%m_dot(sys%cv(i_from_cv), sys%cv(i_to_cv))
            h_dot(i_from_cv, i_to_cv) = sys%cv(i_from_cv)%h() * m_dot(i_from_cv, i_to_cv)
        end do
    end do
end subroutine calculate_flows

pure function m_total_sys(sys)
    class(cv_system_type), intent(in) :: sys
    
    type(si_mass) :: m_total_sys
    
    integer :: i_cv
    
    call m_total_sys%v%init_const(0.0_WP, size(sys%cv(1)%m(1)%v%d))
    do i_cv = 1, size(sys%cv)
        m_total_sys = m_total_sys + sys%cv(i_cv)%m_total()
    end do
end function m_total_sys

pure function e_total_sys(sys)
    class(cv_system_type), intent(in) :: sys
    
    type(si_energy) :: e_total_sys
    
    integer :: i_cv
    
    call e_total_sys%v%init_const(0.0_WP, size(sys%cv(1)%e%v%d))
    do i_cv = 1, size(sys%cv)
        e_total_sys = e_total_sys + sys%cv(i_cv)%e_total()
    end do
end function e_total_sys

pure subroutine time_step(sys_old, dt, sys_new)
    ! Advances by one time step.
    
    type(cv_system_type), allocatable, intent(in)  :: sys_old
    type(si_time), intent(in)                      :: dt
    type(cv_system_type), allocatable, intent(out) :: sys_new
    
    type(cv_delta_type), allocatable :: cv_delta_0(:), cv_delta_1(:), cv_delta_2(:), cv_delta_3(:), cv_delta_4(:)
    
    integer :: i_cv, n_cv, k_gas, n_gas, n_d
    
    n_cv  = size(sys_old%cv)
    n_gas = size(sys_old%cv(1)%m)
    n_d   = size(sys_old%cv(1)%m(1)%v%d)
    
    ! stage 1
    allocate(cv_delta_0(n_cv))
    do i_cv = 1, n_cv
        allocate(cv_delta_0(i_cv)%m(n_gas))
    end do
    do i_cv = 1, n_cv
        call cv_delta_0(i_cv)%x%v%init_const(0.0_WP, n_d)
        call cv_delta_0(i_cv)%x_dot%v%init_const(0.0_WP, n_d)
        call cv_delta_0(i_cv)%e%v%init_const(0.0_WP, n_d)
        call cv_delta_0(i_cv)%e_f%v%init_const(0.0_WP, n_d)
        do k_gas = 1, n_gas
            call cv_delta_0(i_cv)%m(k_gas)%v%init_const(0.0_WP, n_d)
        end do
    end do
    call rk_stage(dt, 1.0_WP, sys_old, cv_delta_0, cv_delta_1)
    
    ! stage 2
    call rk_stage(dt, 0.5_WP, sys_old, cv_delta_1, cv_delta_2)
    
    ! stage 3
    call rk_stage(dt, 0.5_WP, sys_old, cv_delta_2, cv_delta_3)
    
    ! stage 4
    call rk_stage(dt, 1.0_WP, sys_old, cv_delta_3, cv_delta_4)
    
    ! Put it all together.
    sys_new = sys_old
    do i_cv = 1, n_cv
        sys_new%cv(i_cv)%x     = sys_old%cv(i_cv)%x + (cv_delta_1(i_cv)%x &
                                                        + 2.0_WP*cv_delta_2(i_cv)%x &
                                                        + 2.0_WP*cv_delta_3(i_cv)%x &
                                                        + cv_delta_4(i_cv)%x &
                                                        )/6.0_WP
        sys_new%cv(i_cv)%x_dot = sys_old%cv(i_cv)%x_dot + (cv_delta_1(i_cv)%x_dot &
                                                        + 2.0_WP*cv_delta_2(i_cv)%x_dot &
                                                        + 2.0_WP*cv_delta_3(i_cv)%x_dot &
                                                        + cv_delta_4(i_cv)%x_dot &
                                                        )/6.0_WP
        sys_new%cv(i_cv)%e     = sys_old%cv(i_cv)%e + (cv_delta_1(i_cv)%e &
                                                        + 2.0_WP*cv_delta_2(i_cv)%e &
                                                        + 2.0_WP*cv_delta_3(i_cv)%e &
                                                        + cv_delta_4(i_cv)%e &
                                                        )/6.0_WP
        sys_new%cv(i_cv)%e_f   = sys_old%cv(i_cv)%e_f + (cv_delta_1(i_cv)%e_f &
                                                        + 2.0_WP*cv_delta_2(i_cv)%e_f &
                                                        + 2.0_WP*cv_delta_3(i_cv)%e_f &
                                                        + cv_delta_4(i_cv)%e_f &
                                                        )/6.0_WP
        do k_gas = 1, n_gas
            sys_new%cv(i_cv)%m(k_gas) = sys_old%cv(i_cv)%m(k_gas) + (cv_delta_1(i_cv)%m(k_gas) &
                                                                        + 2.0_WP*cv_delta_2(i_cv)%m(k_gas) &
                                                                        + 2.0_WP*cv_delta_3(i_cv)%m(k_gas) &
                                                                        + cv_delta_4(i_cv)%m(k_gas) &
                                                                        )/6.0_WP
        end do
        
        call assert_mass(sys_new%cv(i_cv), "time_step")
    end do
end subroutine time_step

pure subroutine rk_stage(dt, a, sys_old, cv_delta_in, cv_delta_out)
    type(si_time), intent(in)                      :: dt
    real(WP), intent(in)                           :: a
    type(cv_system_type), allocatable, intent(in)  :: sys_old
    type(cv_delta_type), allocatable, intent(in)   :: cv_delta_in(:)
    type(cv_delta_type), allocatable, intent(out)  :: cv_delta_out(:)
    
    type(cv_system_type), allocatable      :: sys
    type(si_mass_flow_rate), allocatable   :: m_dot(:, :)
    type(si_energy_flow_rate), allocatable :: h_dot(:, :)
    integer :: i_cv, n_cv, k_gas, n_gas
    
    n_cv  = size(sys_old%cv)
    n_gas = size(sys_old%cv(1)%m)
    
    allocate(cv_delta_out(n_cv))
    do i_cv = 1, n_cv
        allocate(cv_delta_out(i_cv)%m(n_gas))
    end do
    
    sys = sys_old
    do i_cv = 1, n_cv
        sys%cv(i_cv)%x     = sys_old%cv(i_cv)%x     + a*cv_delta_in(i_cv)%x
        sys%cv(i_cv)%x_dot = sys_old%cv(i_cv)%x_dot + a*cv_delta_in(i_cv)%x_dot
        sys%cv(i_cv)%e     = sys_old%cv(i_cv)%e     + a*cv_delta_in(i_cv)%e
        sys%cv(i_cv)%e_f   = sys_old%cv(i_cv)%e_f   + a*cv_delta_in(i_cv)%e_f
        do k_gas = 1, n_gas
            sys%cv(i_cv)%m(k_gas) = sys_old%cv(i_cv)%m(k_gas) + a*cv_delta_in(i_cv)%m(k_gas)
        end do
    end do
    call sys%calculate_flows(m_dot, h_dot)
    do i_cv = 1, n_cv
        cv_delta_out(i_cv)%x     = dt*d_x_d_t(sys, i_cv)
        cv_delta_out(i_cv)%x_dot = dt*d_x_dot_d_t(sys, i_cv)
        cv_delta_out(i_cv)%e     = dt*d_e_d_t(sys%cv(i_cv), h_dot, i_cv)
        cv_delta_out(i_cv)%e_f   = dt*d_e_f_d_t(sys, i_cv)
        do k_gas = 1, n_gas
            cv_delta_out(i_cv)%m(k_gas) = dt*d_m_k_d_t(sys, m_dot, k_gas, i_cv)
        end do
    end do
end subroutine rk_stage

subroutine set_run_config(config, id, n_d, csv_output, csv_frequency, t_stop, dt, tolerance_checks)
    class(run_config_type), intent(out) :: config
    character(len=*), intent(in)        :: id ! CSV file name
    integer, intent(in)                 :: n_d
    
    logical, intent(in), optional       :: csv_output
    integer, intent(in), optional       :: csv_frequency
    type(si_time), intent(in), optional :: t_stop, dt
    logical, intent(in), optional       :: tolerance_checks
    
    config%id = id
    
    if (present(csv_output)) then
        config%csv_output = csv_output
    else
        config%csv_output = .false.
    end if
    
    if (present(csv_frequency)) then
        call assert(csv_output, "cva (set_run_config): Why set csv_frequency if csv_output=.false.?", &
                        print_logical=[csv_output])
        config%csv_frequency = csv_frequency
    else
        config%csv_frequency = 10
    end if
    
    if (present(t_stop)) then
        config%t_stop = t_stop
    else
        call config%t_stop%v%init_const(T_STOP_DEFAULT, n_d)
    end if
    
    if (present(dt)) then
        config%dt = dt
    else
        call config%dt%v%init_const(DT_DEFAULT, n_d)
    end if
    
    if (present(tolerance_checks)) then
        config%tolerance_checks = tolerance_checks
    else
        config%tolerance_checks = .true.
    end if
end subroutine set_run_config

subroutine run(config, sys_start, sys_end, status)
    use, intrinsic :: iso_fortran_env, only: ERROR_UNIT
    use prec, only: CL
    
    type(run_config_type), intent(in)              :: config
    type(cv_system_type), allocatable, intent(in)  :: sys_start
    type(cv_system_type), allocatable, intent(out) :: sys_end
    type(run_status_type), intent(out)             :: status
    
    type(cv_system_type), allocatable :: sys_old, sys_new, sys_temp
    
    character(len=CL)     :: error_message
    integer               :: n_d, i, csv_unit
    type(si_time)         :: t, t_old
    type(run_status_type) :: status_final
    
    n_d = size(sys_start%cv(1)%x%v%d)
    
    sys_old = sys_start
    call t%v%init_const(0.0_WP, n_d)
    i = 0
    
    if (config%csv_output) then
        open(newunit=csv_unit, action="write", status="replace", position="rewind", &
                file=trim(config%id) // ".csv", iostat=status%rc, iomsg=error_message)
        if (status%rc /= 0) then
            write(unit=ERROR_UNIT, fmt="(a)") trim(error_message)
            status%t = t
            return
        end if
        
        status%rc = CONTINUE_RUN_RC
        call write_csv_row(csv_unit, sys_old, t, status, HEADER_ROW_TYPE)
        call write_csv_row(csv_unit, sys_old, t, status, NUMBER_ROW_TYPE)
    end if
    
    time_loop: do
        call time_step(sys_old, config%dt, sys_new)
        t_old = t
        t     = t + config%dt
        i     = i + 1
        
        !print *, t%v%v
        
        call check_sys(config, sys_new, sys_start, t, status)
        if (status%rc >= 0) exit time_loop
        if ((config%csv_output) .and. (mod(i, config%csv_frequency) == 0)) then
            call write_csv_row(csv_unit, sys_new, t, status, NUMBER_ROW_TYPE)
        end if
        
        call move_alloc(from=sys_old,  to=sys_temp)
        call move_alloc(from=sys_new,  to=sys_old)
        call move_alloc(from=sys_temp, to=sys_new)
    end do time_loop
    
    if (status%rc == SUCCESS_RUN_RC) then
        ! If successful, use Hénon's trick to ensure that `x == x_stop`.
        call sys_interp(t_old, config%dt, status%i_cv(1), sys_old, sys_new, t, sys_end)
        
        call check_sys(config, sys_end, sys_start, t, status_final)
        ! This will make sure that the status rc stays the same if everything is okay but `x` is a bit short.
        if (status_final%rc /= CONTINUE_RUN_RC) then
            status = status_final
        end if
    else
        sys_end = sys_new
    end if
    
    if (config%csv_output) then
        call write_csv_row(csv_unit, sys_end, t, status, NUMBER_ROW_TYPE)
        close(unit=csv_unit)
    end if
    
    status%t = t
end subroutine run

pure subroutine check_sys(config, sys, sys_start, t, status)
    type(run_config_type), intent(in)             :: config
    type(cv_system_type), allocatable, intent(in) :: sys, sys_start
    type(si_time), intent(in)                     :: t
    type(run_status_type), intent(out)            :: status
    
    integer              :: n_cv, i_cv, j_cv, n_bad_cv, n_d, i_d, i_d_max
    type(si_mass)        :: m_start, m_total_i, m_total_j, rel_m
    type(si_energy)      :: e_start, rel_e
    type(si_temperature) :: temp_i, temp_j
    type(si_pressure)    :: p_j
    type(unitless)       :: rel_delta
    real(WP)             :: max_abs_m_deriv, max_abs_e_deriv
    type(si_length)      :: x_sum_start, x_sum
    
    n_cv      = size(sys%cv)
    n_d       = size(sys%cv(1)%x%v%d)
    status%rc = CONTINUE_RUN_RC
    
    m_start = sys_start%m_total()
    e_start = sys_start%e_total()
    
    do i_cv = 1, n_cv
        ! Check whether the projectile left the barrel.
        if (sys%cv(i_cv)%x >= sys%cv(i_cv)%x_stop) then
            status%rc = SUCCESS_RUN_RC
            allocate(status%i_cv(1))
            status%i_cv(1) = i_cv
            return
        end if
        
        ! Check that masses of each control volume are positive.
        m_total_i = sys%cv(i_cv)%m_total()
        if (m_total_i%v%v < 0.0_WP) then
            status%rc = NEGATIVE_CV_M_TOTAL_RUN_RC
            
            allocate(status%data(n_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                m_total_j = sys%cv(j_cv)%m_total()
                if (m_total_j%v%v < 0.0_WP) n_bad_cv = n_bad_cv + 1
                status%data(j_cv) = m_total_j%v%v
            end do
            call assert(n_bad_cv >= 1, "cva (check_sys): number of control volumes with negative mass should be >= 1", &
                            print_integer=[n_bad_cv])
            
            allocate(status%i_cv(n_bad_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                m_total_j = sys%cv(j_cv)%m_total()
                if (m_total_j%v%v < 0.0_WP) then
                    n_bad_cv = n_bad_cv + 1
                    status%i_cv(n_bad_cv) = j_cv
                end if
            end do
            
            return
        end if
        
        ! Check that temperatures of each control volume are positive.
        temp_i = sys%cv(i_cv)%temp()
        if (temp_i%v%v <= 0.0_WP) then
            status%rc = NEGATIVE_CV_TEMP_RUN_RC
            
            allocate(status%data(n_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                temp_j = sys%cv(j_cv)%temp()
                if (temp_j%v%v < 0.0_WP) n_bad_cv = n_bad_cv + 1
                status%data(j_cv) = temp_j%v%v
            end do
            call assert(n_bad_cv >= 1, &
                "cva (check_sys): number of control volumes with negative or zero temperature should be >= 1", &
                            print_integer=[n_bad_cv])
            
            allocate(status%i_cv(n_bad_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                temp_j = sys%cv(j_cv)%temp()
                if (temp_j%v%v < 0.0_WP) then
                    n_bad_cv = n_bad_cv + 1
                    status%i_cv(n_bad_cv) = j_cv
                end if
            end do
            
            return
        end if
        
        if (sys%cv(i_cv)%eos == IDEAL_EOS) then
            ! Check that the ideal gas law is valid if it is used.
            if (sys%cv(i_cv)%p() >= sys%cv(i_cv)%p_c()) then
                status%rc = IDEAL_EOS_RUN_RC
                
                allocate(status%data(n_cv))
                n_bad_cv = 0
                do j_cv = 1, n_cv
                    p_j = sys%cv(j_cv)%p()
                    if (p_j >= sys%cv(j_cv)%p_c()) n_bad_cv = n_bad_cv + 1
                    status%data(j_cv) = p_j%v%v
                end do
                call assert(n_bad_cv >= 1, "cva (check_sys): number of control volumes with p >= p_c should be >= 1", &
                            print_integer=[n_bad_cv])
                
                allocate(status%i_cv(n_bad_cv))
                n_bad_cv = 0
                do j_cv = 1, n_cv
                    p_j = sys%cv(j_cv)%p()
                    if (p_j >= sys%cv(j_cv)%p_c()) then
                        n_bad_cv = n_bad_cv + 1
                        status%i_cv(n_bad_cv) = j_cv
                    end if
                end do
                
                return
            end if
        end if
        
        ! Check that the mirror CV position invariant is satisfied.
        ! If one CV's plunger moves by a certain amount, the mirror CV's plunger must move the same amount in the opposite direction.
        ! I could also add a tolerance for `x_dot`, but given that it's set to exactly the same, an assertion would make more sense.
        if ((sys%cv(i_cv)%type == NORMAL_CV_TYPE) .and. (sys%cv(i_cv)%i_cv_mirror >= 1)) then
            call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == MIRROR_CV_TYPE, &
                            "cva (check_sys): mirror CV is not a mirror CV?", &
                            print_integer=[sys%cv(sys%cv(i_cv)%i_cv_mirror)%type, MIRROR_CV_TYPE, i_cv, sys%cv(i_cv)%i_cv_mirror])
            
            x_sum_start = sys_start%cv(i_cv)%x + sys_start%cv(sys%cv(i_cv)%i_cv_mirror)%x
            x_sum       = sys%cv(i_cv)%x       + sys%cv(sys%cv(i_cv)%i_cv_mirror)%x
            rel_delta   = abs(x_sum - x_sum_start) / x_sum_start
            
            if (rel_delta%v%v > MIRROR_X_TOLERANCE) then
                status%rc = MIRROR_X_TOLERANCE_RUN_RC
                allocate(status%data(1))
                status%data(1) = rel_delta%v%v
                allocate(status%i_cv(1))
                status%i_cv(1) = i_cv
                return
            end if
        end if
    end do
    
    if (config%tolerance_checks) then
        ! Check that total system mass is staying constant.
        call assert(m_start%v%v > 0.0_WP, "cva (check_sys): m_start must be greater than zero", print_real=[m_start%v%v])
        rel_delta = abs(sys%m_total() - m_start) / m_start
        if (rel_delta%v%v > MASS_TOLERANCE) then
            status%rc = MASS_TOLERANCE_RUN_RC
            allocate(status%data(1))
            status%data(1) = rel_delta%v%v
            return
        end if
        
        ! Why not check that the total system mass *of each species* is staying constant?
        ! I might want to add combustion to this later, and mass of each species won't remain constant in that case.
        
        ! Check that derivatives of mass are staying constant.
        ! It appears that dividing by `m_start` like with `rel_delta` makes the derivatives too small.
        rel_m = sys%m_total() - m_start
        max_abs_m_deriv = 0.0_WP
        do i_d = 1, n_d
            if (abs(rel_m%v%d(i_d)) > max_abs_m_deriv) then
                max_abs_m_deriv = abs(rel_m%v%d(i_d))
                i_d_max = i_d
            end if
        end do
        if (max_abs_m_deriv > MASS_DERIV_TOLERANCE) then
            status%rc = MASS_DERIV_TOLERANCE_RUN_RC
            allocate(status%data(2))
            status%data(1) = max_abs_m_deriv
            status%data(2) = real(i_d_max, WP)
            return
        end if
        
        ! Check that total system energy is staying constant.
        call assert(e_start%v%v > 0.0_WP, "cva (check_sys): e_start must be greater than zero", print_real=[e_start%v%v])
        rel_delta = abs(sys%e_total() - e_start) / e_start
        if (rel_delta%v%v > ENERGY_TOLERANCE) then
            status%rc = ENERGY_TOLERANCE_RUN_RC
            allocate(status%data(1))
            status%data(1) = rel_delta%v%v
            return
        end if
        
        ! Check that derivatives of energy are staying constant.
        rel_e = sys%e_total() - e_start
        max_abs_e_deriv = 0.0_WP
        do i_d = 1, n_d
            if (abs(rel_e%v%d(i_d)) > max_abs_e_deriv) then
                max_abs_e_deriv = abs(rel_e%v%d(i_d))
                i_d_max = i_d
            end if
        end do
        if (max_abs_e_deriv > ENERGY_DERIV_TOLERANCE) then
            status%rc = ENERGY_DERIV_TOLERANCE_RUN_RC
            allocate(status%data(2))
            status%data(1) = max_abs_e_deriv
            status%data(2) = real(i_d_max, WP)
            return
        end if
    end if
    
    ! Check if the simulation has timed out.
    if (t >= config%t_stop) then
        status%rc = TIMEOUT_RUN_RC
        return
    end if
end subroutine check_sys

pure subroutine sys_interp(t_old, dt, i_cv_interp, sys_old, sys_new, t, sys_end)
    ! Use Hénon's trick to ensure that `x == x_stop`.
    ! Uses the secant method to find `dt` where `x == x_stop` for control volume number `i_cv_interp`.
    
    type(si_time), intent(in)                      :: t_old, dt
    integer, intent(in)                            :: i_cv_interp ! control volume number to interpolate based on
    type(cv_system_type), allocatable, intent(in)  :: sys_old, sys_new
    type(si_time), intent(out)                     :: t
    type(cv_system_type), allocatable, intent(out) :: sys_end
    
    integer       :: i
    type(si_time) :: dt_i, dt_im1, dt_im2
    real(WP)      :: x_tol
    type(cv_system_type), allocatable :: sys_i, sys_im1, sys_im2, sys_temp
    integer, parameter :: MAX_ITERS = 10
    
    x_tol = 100.0_WP*spacing(sys_old%cv(i_cv_interp)%x_stop%v%v)
    
    call assert(is_close(sys_old%cv(i_cv_interp)%x_stop%v%v, sys_new%cv(i_cv_interp)%x_stop%v%v), &
                    "cva (sys_interp): x_stop is inconsistent", &
                    print_real=[sys_old%cv(i_cv_interp)%x_stop%v%v, sys_new%cv(i_cv_interp)%x_stop%v%v])
    
    call dt_im2%v%init_const(0.0_WP, size(sys_new%cv(i_cv_interp)%x_stop%v%d))
    dt_im1 = dt
    
    sys_im2 = sys_old
    sys_im1 = sys_new
    
    do i = 1, MAX_ITERS
        ! The stopping criteria is based on the difference between iterates due to risk of catastrophic cancellation.
        if (abs(sys_im1%cv(i_cv_interp)%x%v%v - sys_im2%cv(i_cv_interp)%x%v%v) < x_tol) exit 
        
        ! <https://en.wikipedia.org/wiki/Secant_method#The_method>
        dt_i = dt_im1 - (dt_im1 - dt_im2) * (sys_im1%cv(i_cv_interp)%x - sys_old%cv(i_cv_interp)%x_stop) &
                    / (sys_im1%cv(i_cv_interp)%x - sys_im2%cv(i_cv_interp)%x)
        
        !print *, dt_i%v%v, dt_im1%v%v, dt_im2%v%v
        
        call assert(dt_i%v%v >= 0.0_WP, "cva (sys_interp): dt_i can not be negative", &
                        print_real=[dt_i%v%v], print_integer=[i])
        call assert(dt_i <= dt, "cva (sys_interp): dt_i can not be greater than dt", &
                        print_real=[dt_i%v%v], print_integer=[i])
        
        call time_step(sys_old, dt_i, sys_i)
        
        call assert(sys_i%cv(i_cv_interp)%x >= sys_old%cv(i_cv_interp)%x, &
                        "cva (sys_interp): x_i >= x_old violated", &
                        print_real=[sys_old%cv(i_cv_interp)%x%v%v, sys_i%cv(i_cv_interp)%x%v%v, &
                        sys_new%cv(i_cv_interp)%x%v%v, dt_i%v%v], print_integer=[i])
        call assert(sys_i%cv(i_cv_interp)%x <= sys_new%cv(i_cv_interp)%x, &
                        "cva (sys_interp): x_i <= x_new violated", &
                        print_real=[sys_old%cv(i_cv_interp)%x%v%v, sys_i%cv(i_cv_interp)%x%v%v, &
                        sys_new%cv(i_cv_interp)%x%v%v, dt_i%v%v], print_integer=[i])
        
        dt_im2 = dt_im1
        dt_im1 = dt_i
        
        call move_alloc(from=sys_im2,  to=sys_temp)
        call move_alloc(from=sys_im1,  to=sys_im2)
        call move_alloc(from=sys_i,    to=sys_im1)
        call move_alloc(from=sys_temp, to=sys_i)
    end do
    
    sys_end = sys_i
    call assert(is_close(sys_end%cv(i_cv_interp)%x%v%v, sys_end%cv(i_cv_interp)%x_stop%v%v, &
                            abs_tol=max(100.0_WP*x_tol, 1.0e-6_WP)), &
                    "cva (sys_interp): x_end is not close to x_stop", &
                    print_real=[sys_end%cv(i_cv_interp)%x%v%v, sys_end%cv(i_cv_interp)%x_stop%v%v])
    
    t = t_old + dt_i
end subroutine sys_interp

!tripwire$ begin C4615087 Update `\secref{csv}` of usage.tex when changing `write_csv_row`.
subroutine write_csv_row(csv_unit, sys, t, status, row_type)
    integer, intent(in)                           :: csv_unit
    type(cv_system_type), allocatable, intent(in) :: sys
    type(si_time), intent(in)                     :: t
    type(run_status_type), intent(in)             :: status
    integer, intent(in)                           :: row_type
    
    integer :: i_cv, n_cv, k_gas, n_gas
    logical :: csv_unit_opened
    type(si_pressure)     :: p
    type(si_temperature)  :: temp
    type(si_mass_density) :: rho
    type(si_mass)         :: m_total
    type(si_energy)       :: e_total, e_s, e_p
    
    inquire(unit=csv_unit, opened=csv_unit_opened)
    call assert(csv_unit_opened, "cva (write_csv_row): csv_unit needs to be open", print_logical=[csv_unit_opened])
    
    n_cv  = size(sys%cv)
    n_gas = size(sys%cv(1)%m)
    
    ! for all CVs
    
    ! `t`, time
    select case (row_type)
        case (HEADER_ROW_TYPE)
            write(unit=csv_unit, fmt="(a)", advance="no") '"t (s)",'
        case (NUMBER_ROW_TYPE)
            write(unit=csv_unit, fmt="(g0, a)", advance="no") t%v%v, ","
        case default
            error stop "cva (write_csv_row, t): invalid row_type"
    end select
    
    ! `m_total`
    m_total = sys%m_total()
    select case (row_type)
        case (HEADER_ROW_TYPE)
            write(unit=csv_unit, fmt="(a)", advance="no") '"m_total (kg)",'
        case (NUMBER_ROW_TYPE)
            write(unit=csv_unit, fmt="(g0, a)", advance="no") m_total%v%v, ","
        case default
            error stop "cva (write_csv_row, m_total): invalid row_type"
    end select
    
    ! `e_total`
    e_total = sys%e_total()
    select case (row_type)
        case (HEADER_ROW_TYPE)
            write(unit=csv_unit, fmt="(a)", advance="no") '"e_total (J)",'
        case (NUMBER_ROW_TYPE)
            write(unit=csv_unit, fmt="(g0, a)", advance="no") e_total%v%v, ","
        case default
            error stop "cva (write_csv_row, e_total): invalid row_type"
    end select
    
    ! per CV
    do i_cv = 1, n_cv
        if (sys%cv(i_cv)%type == NORMAL_CV_TYPE) then
            ! `x`, location of projectile/plunger
            select case (row_type)
                case (HEADER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(3a)", advance="no") '"x (m, ', trim(sys%cv(i_cv)%label), ')",'
                case (NUMBER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(g0, a)", advance="no") sys%cv(i_cv)%x%v%v, ","
                case default
                    error stop "cva (write_csv_row, x): invalid row_type"
            end select
            
            ! `x_dot`, velocity of projectile/plunger
            select case (row_type)
                case (HEADER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(3a)", advance="no") '"x_dot (m/s, ', trim(sys%cv(i_cv)%label), ')",'
                case (NUMBER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(g0, a)", advance="no") sys%cv(i_cv)%x_dot%v%v, ","
                case default
                    error stop "cva (write_csv_row, x_dot): invalid row_type"
            end select
        end if
        
        ! `m(:)`, mass(es) of gas(es) in control volume
        do k_gas = 1, n_gas
            select case (row_type)
                case (HEADER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(4a)", advance="no") '"m (kg, ', &
                            trim(sys%cv(i_cv)%gas(k_gas)%label), ', ', &
                            trim(sys%cv(i_cv)%label), ')",'
                case (NUMBER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(g0, a)", advance="no") sys%cv(i_cv)%m%v%v, ","
                case default
                    error stop "cva (write_csv_row, m): invalid row_type"
            end select
        end do
        
        ! `e`, energy of gas in control volume
        select case (row_type)
            case (HEADER_ROW_TYPE)
                write(unit=csv_unit, fmt="(3a)", advance="no") '"e (J, ', trim(sys%cv(i_cv)%label), ')",'
            case (NUMBER_ROW_TYPE)
                write(unit=csv_unit, fmt="(g0, a)", advance="no") sys%cv(i_cv)%e%v%v, ","
            case default
                error stop "cva (write_csv_row, e): invalid row_type"
        end select
        
        ! `p`
        select case (row_type)
            case (HEADER_ROW_TYPE)
                write(unit=csv_unit, fmt="(3a)", advance="no") '"p (Pa, absolute, ', trim(sys%cv(i_cv)%label), ')",'
            case (NUMBER_ROW_TYPE)
                p = sys%cv(i_cv)%p()
                write(unit=csv_unit, fmt="(g0, a)", advance="no") p%v%v, ","
            case default
                error stop "cva (write_csv_row, p): invalid row_type"
        end select
        
        ! `temp`
        select case (row_type)
            case (HEADER_ROW_TYPE)
                write(unit=csv_unit, fmt="(3a)", advance="no") '"temp (K, ', trim(sys%cv(i_cv)%label), ')",'
            case (NUMBER_ROW_TYPE)
                temp = sys%cv(i_cv)%temp()
                write(unit=csv_unit, fmt="(g0, a)", advance="no") temp%v%v, ","
            case default
                error stop "cva (write_csv_row, temp): invalid row_type"
        end select
        
        ! `rho`
        select case (row_type)
            case (HEADER_ROW_TYPE)
                write(unit=csv_unit, fmt="(3a)", advance="no") '"rho (kg/m3, ', trim(sys%cv(i_cv)%label), ')",'
            case (NUMBER_ROW_TYPE)
                rho = sys%cv(i_cv)%rho()
                write(unit=csv_unit, fmt="(g0, a)", advance="no") rho%v%v, ","
            case default
                error stop "cva (write_csv_row, rho): invalid row_type"
        end select
        
        if (sys%cv(i_cv)%type == NORMAL_CV_TYPE) then
            ! `e_f`, energy lost to projectile/plunger friction in control volume
            select case (row_type)
                case (HEADER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(3a)", advance="no") '"e_f (J, ', trim(sys%cv(i_cv)%label), ')",'
                case (NUMBER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(g0, a)", advance="no") sys%cv(i_cv)%e_f%v%v, ","
                case default
                    error stop "cva (write_csv_row, e_f): invalid row_type"
            end select
            
            ! `e_s`
            select case (row_type)
                case (HEADER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(3a)", advance="no") '"e_s (J, ', trim(sys%cv(i_cv)%label), ')",'
                case (NUMBER_ROW_TYPE)
                    e_s = sys%cv(i_cv)%e_s()
                    write(unit=csv_unit, fmt="(g0, a)", advance="no") e_s%v%v, ","
                case default
                    error stop "cva (write_csv_row, e_s): invalid row_type"
            end select
            
            ! `e_p`
            select case (row_type)
                case (HEADER_ROW_TYPE)
                    write(unit=csv_unit, fmt="(3a)", advance="no") '"e_p (J, ', trim(sys%cv(i_cv)%label), ')",'
                case (NUMBER_ROW_TYPE)
                    e_p = sys%cv(i_cv)%e_p()
                    write(unit=csv_unit, fmt="(g0, a)", advance="no") e_p%v%v, ","
                case default
                    error stop "cva (write_csv_row, e_p): invalid row_type"
            end select
        end if
    end do
    
    ! between CVs
    
    ! TODO: m_dot
    ! TODO: h_dot

    select case (row_type)
        case (HEADER_ROW_TYPE)
            write(unit=csv_unit, fmt="(a)") '"return code"'
        case (NUMBER_ROW_TYPE)
            write(unit=csv_unit, fmt="(i0)") status%rc
        case default
            error stop "cva (write_csv_row, rc): invalid row_type"
    end select
end subroutine write_csv_row
!tripwire$ end

end module cva
