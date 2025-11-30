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
public :: d_x_d_t, d_xdot_d_t, d_m_k_d_t, d_e_d_t
public :: f_m_dot, g_m_dot
public :: time_step, run

! pressure ratio laminar flow nominally starts at
! based on first part of beater_pneumatic_2007 eq. 5.4
real(WP), public, parameter :: P_RL = 0.999_WP ! unitless

real(WP), public, parameter :: X_STOP_DEFAULT = 1.0e3_WP  ! m (If you have a barrel that's a km long, that's probably wrong.)
real(WP), public, parameter :: DT_DEFAULT     = 1.0e-5_WP ! s
real(WP), public, parameter :: T_STOP_DEFAULT = 0.5_WP    ! s

integer, public, parameter :: IDEAL_EOS = 1 ! ideal gas equation of state
integer, public, parameter :: CONST_EOS = 2 ! constant pressure, temperature, density
integer, public, parameter :: RK_EOS    = 3 ! Redlich–Kwong equation of state
integer, public, parameter :: MAX_EOS   = 2

integer, public, parameter :: NORMAL_CV_TYPE = 1
integer, public, parameter :: MIRROR_CV_TYPE = 2
integer, public, parameter :: MAX_CV_TYPE    = 2

integer, public, parameter :: NORMAL_RUN_RC  = 0
integer, public, parameter :: TIMEOUT_RUN_RC = 1
integer, public, parameter :: MASS_RUN_RC    = 2
integer, public, parameter :: TEMP_RUN_RC    = 3

type, public :: cv_type ! control volume
    ! time varying
    type(si_length)            :: x     ! location of piston/projectile
    type(si_velocity)          :: x_dot ! velocity of piston/projectile
    type(si_mass), allocatable :: m(:)  ! mass(es) of gas(es) in control volume
    type(si_energy)            :: e     ! energy of gas in control volume
    
    ! constants
    character(len=32)           :: label       ! human-readable label for control volume
    integer                     :: eos         ! equation of state to use for control volume
    integer                     :: type        ! type of control volume
    type(si_area)               :: csa         ! cross-sectional area
    type(si_inverse_mass)       :: rm_p        ! reciprocal mass of piston/projectile
    type(si_pressure)           :: p_fs, p_fd  ! static and dynamic friction pressure
    type(si_stiffness)          :: k           ! stiffness of spring attached to piston
    type(si_length)             :: x_z         ! zero force location for spring
    type(gas_type), allocatable :: gas(:)      ! gas data
    integer                     :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    type(si_pressure)           :: p_const    ! if `cv%eos = CONST_EOS`, then `cv%p() = p_const`
    type(si_temperature)        :: temp_const ! if `cv%eos = CONST_EOS`, then `cv%temp() = temp_const`
    type(si_length)             :: x_stop     ! `x` location where simulation will stop
contains
    procedure :: m_total
    procedure :: spring_pe
    procedure :: m_p_ke
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
    type(si_length)            :: x     ! delta of location of piston/projectile
    type(si_velocity)          :: x_dot ! delta of velocity of piston/projectile
    type(si_mass), allocatable :: m(:)  ! delta of mass(es) of gas(es) in control volume
    type(si_energy)            :: e     ! delta of energy of gas in control volume
end type cv_delta_type

type, public :: con_type ! connection between control volumes
    logical        :: active
    type(si_area)  :: a_e ! effective area
    type(unitless) :: b   ! critical pressure ratio
contains
    procedure :: m_dot
end type con_type

type, public :: cv_system_type
    type(cv_type), allocatable  :: cv(:)
    type(con_type), allocatable :: con(:, :)
contains
    procedure :: calculate_flows
    procedure :: m_total => m_total_sys
    procedure :: e_total => e_total_sys
end type cv_system_type

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
    
    integer :: i_gas
    
    call m_total%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i_gas = 1, size(cv%m)
        m_total = m_total + cv%m(i_gas)
    end do
end function m_total

pure function spring_pe(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_energy) :: spring_pe
    
    if (cv%type == MIRROR_CV_TYPE) then
        call spring_pe%v%init_const(0.0_WP, size(cv%x_dot%v%d))
    else
        spring_pe = 0.5_WP*cv%k*square(cv%x - cv%x_z)
    end if
end function spring_pe

pure function m_p_ke(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_energy) :: m_p_ke
    
    if (is_close(cv%rm_p%v%v, 0.0_WP) .or. (cv%type == MIRROR_CV_TYPE)) then
        if (is_close(cv%rm_p%v%v, 0.0_WP)) then
            ! `MIRROR_CV_TYPE` might simply copy what the other CV has, so it might not have no inverse mass.
            call assert(is_close(cv%k%v%v, 0.0_WP), "cva (m_p_ke): if the piston is immobile, k should be zero")
        end if
        
        call m_p_ke%v%init_const(0.0_WP, size(cv%x_dot%v%d))
    else
        m_p_ke = 0.5_WP*square(cv%x_dot)/cv%rm_p
    end if
end function m_p_ke

pure function e_total(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_energy) :: e_total
    
    e_total = cv%e + cv%spring_pe() + cv%m_p_ke()
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
            call assert(rho%v%v  > 0.0_WP, "cva (p_eos, IDEAL_EOS): rho%v > 0 violated")
            call assert(temp%v%v > 0.0_WP, "cva (p_eos, IDEAL_EOS): temp%v > 0 violated")
            call assert_dimension(rho%v%d, temp%v%d)
            
            n_d   = size(rho%v%d)
            p_eos = rho * cv%r() * temp
            
            call assert(p_eos < cv%p_c(), "cva (p_eos, IDEAL_EOS): ideal gas law validity is questionable")
            call assert(p_eos%v%v > 0.0_WP, "cva (p_eos, IDEAL_EOS): p_eos%v > 0 violated")
        case (CONST_EOS)
            call assert(is_close(temp%v%v, cv%temp_const%v%v), "cva (p_eos, CONST_EOS): temp /= temp_const")
            p_eos = cv%p_const
            call assert(p_eos%v%v >= 0.0_WP, "cva (p_eos, CONST_EOS): p_eos%v >= 0 violated")
        case default
            error stop "cva (p_eos): invalid cv%eos"
    end select
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
    
    call assert(cv%eos == IDEAL_EOS, "cva (rho_eos): ideal equation of state required")
    call assert(p%v%v    >  0.0_WP, "cva (rho_eos): p%v > 0 violated")
    call assert(temp%v%v >  0.0_WP, "cva (rho_eos): temp%v > 0 violated")
    call assert(size(y)  >= 1,      "cva (rho_eos): size(y) >= 1 violated")
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
        call assert(denominator%v%v > 0.0_WP, "cva (rho_eos): denominator is zero")
        call gas_mm%v%init_const(cv%gas(i)%mm, n_d)
        mm    = mm + y(i)*gas_mm/denominator
        y_sum = y_sum + y(i)
    end do
    
    call assert(is_close(y_sum%v%v, 1.0_WP), "cva (rho_eos): y does not sum to 1")
    
    call r_bar_%v%init_const(R_BAR, n_d)
    r_cv = r_bar_ / mm
    
    rho_eos = p / (r_cv * temp)
    
    call assert(rho_eos%v%v > 0.0_WP, "cva (rho_eos): rho_eos%v > 0 violated")
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
        call assert(denominator%v%v > 0.0_WP, "cva (p_c): denominator is zero")
        chi = cv%m(i) / denominator
        call p_ci%v%init_const(cv%gas(i)%p_c, size(cv%m(1)%v%d))
        p_c = p_c + chi*p_ci
        chi_sum = chi_sum + chi
    end do
    
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (p_c): chi does not sum to 1")
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
        call assert(y(i)%v%v >= 0.0_WP, "cva (y): y >= 0 violated")
        call assert(y(i)%v%v <= 1.0_WP, "cva (y): y <= 1 violated")
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
    
    call assert(size(cv%m) >= 1, "cva (chi): size(cv%m) >= 1 violated")
    
    call chi_sum%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        call denominator%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        do j = 1, size(cv%m)
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm) ! MAYBE: change so that the molar masses have units?
            call assert(denominator%v%v >= 0.0_WP, "cva (chi): denominator >= 0 violated")
        end do
        call assert(denominator%v%v > 0.0_WP, "cva (chi): denominator is zero")
        chi(i)  = cv%m(i) / denominator
        chi_sum = chi_sum + chi(i)
        call assert(chi(i)%v%v >= 0.0_WP, "cva (chi): chi >= 0 violated")
        call assert(chi(i)%v%v <= 1.0_WP, "cva (chi): chi <= 1 violated")
    end do
    
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (chi): chi does not sum to 1")
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
        call assert(denominator%v%v > 0.0_WP, "cva (r_cv): denominator is zero")
        chi = cv%m(i) / denominator
        chi_sum = chi_sum + chi
        call gas_mm%v%init_const(cv%gas(i)%mm, n_d)
        mm  = mm + chi*gas_mm
    end do
    
    call r_bar_%v%init_const(R_BAR, n_d)
    r_cv = r_bar_ / mm
    
    call assert(r_cv%v%v > 0.0_WP, "cva (r_cv): r_cv > 0 violated")
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (r_cv): chi does not sum to 1")
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
    
    call assert(temp_cv%v%v > 0.0_WP, "cva (temp_cv): temp_cv > 0 violated")
end function temp_cv

pure function vol_cv(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_volume) :: vol_cv
    
    call assert(cv%x%v%v > 0.0_WP, "cva (vol_cv): cv%x > 0 violated")
    call assert(cv%csa%v%v > 0.0_WP, "cva (vol_cv): cv%csa > 0 violated")
    
    vol_cv = cv%x * cv%csa
    
    call assert(vol_cv%v%v > 0.0_WP, "cva (vol_cv): vol_cv > 0 violated")
end function vol_cv

pure function rho_cv(cv)
    ! Returns simply the mass divided by the volume.
    ! For `CONST_EOS`, this may not be the desired mass density.
    
    class(cv_type), intent(in) :: cv
    
    type(si_mass_density) :: rho_cv
    
    type(si_volume) :: vol
    
    call assert_mass(cv, "rho_cv")
    
    vol = cv%vol()
    call assert(vol%v%v > 0.0_WP, "cva (rho_cv): volume is zero")
    rho_cv = cv%m_total() / vol
    
    call assert(rho_cv%v%v >= 0.0_WP, "cva (rho_cv): rho_cv > 0 violated")
end function rho_cv

pure function p_cv(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_pressure) :: p_cv
    
    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    
    rho  = cv%rho()
    temp = cv%temp()
    
    p_cv = cv%p_eos(rho, temp)
    
    call assert(p_cv%v%v >= 0.0_WP, "cva (p_cv): p_cv >= 0 violated")
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
    
    call assert(u_cv%v%v > 0.0_WP, "cva (u_cv): u_cv > 0 violated")
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
            call y%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        end if
        h_cv = h_cv + y*cv%gas(i)%h(temp)
    end do
    
    call assert(h_cv%v%v >= 0.0_WP, "cva (h_cv): h_cv > 0 violated")
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
    
    call assert(gamma_cv%v%v > 1.0_WP, "cva (gamma_cv): gamma_cv > 1 violated")
end function gamma_cv

pure subroutine set(cv, x, x_dot, y, p, temp_atm, label, csa, rm_p, p_fs, p_fd, k, x_z, gas, &
                            i_cv_mirror, x_stop, isentropic_filling, p_atm, eos, type)
    class(cv_type), intent(in out) :: cv
    
    ! time varying
    type(si_length), intent(in)      :: x        ! location of piston/projectile
    type(si_velocity), intent(in)    :: x_dot    ! velocity of piston/projectile
    type(unitless), intent(in)       :: y(:)     ! mass fractions of each gas
    type(si_pressure), intent(in)    :: p        ! pressure
    
    ! constant
    type(si_temperature), intent(in)  :: temp_atm    ! atmospheric temperature
    character(len=*), intent(in)      :: label       ! human-readable label for control volume
    type(si_area), intent(in)         :: csa         ! cross-sectional area
    type(si_inverse_mass), intent(in) :: rm_p        ! reciprocal mass of piston/projectile
    type(si_pressure), intent(in)     :: p_fs, p_fd  ! static and dynamic friction pressure
    type(si_stiffness), intent(in)    :: k           ! stiffness of spring attached to piston
    type(si_length), intent(in)       :: x_z         ! zero force location for spring
    type(gas_type), intent(in)        :: gas(:)      ! gas data
    integer, intent(in)               :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    
    type(si_length), intent(in), optional   :: x_stop ! `x` location where simulation will stop
    logical, intent(in), optional           :: isentropic_filling
    type(si_pressure), intent(in), optional :: p_atm  ! atmospheric pressure (only requried if `isentropic_filling = .true.`
    integer, intent(in), optional           :: eos    ! equation of state to use
    integer, intent(in), optional           :: type   ! type of CV to use
    
    integer              :: i, n_d
    type(si_temperature) :: temp
    type(si_mass)        :: m_total
    type(unitless)       :: gamma_cv, y_sum
    logical              :: isentropic_filling_
    
    n_d = size(x%v%d)
    
    cv%x     = x
    cv%x_dot = x_dot
    ! `p` and `temp` will be handled below
    
    cv%label       = label
    cv%csa         = csa
    cv%rm_p        = rm_p
    cv%p_fs        = p_fs
    cv%p_fd        = p_fd
    cv%k           = k
    cv%x_z         = x_z
    cv%gas         = gas
    cv%i_cv_mirror = i_cv_mirror
    
    if (present(x_stop)) then
        cv%x_stop = x_stop
    else
        call cv%x_stop%v%init_const(X_STOP_DEFAULT, n_d)
    end if
    
    if (present(isentropic_filling)) then
        isentropic_filling_ = isentropic_filling
    else
        isentropic_filling_ = .false.
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
    
    call assert(cv%x%v%v            >  0.0_WP, "cva (set): x > 0 violated")
    call assert(p%v%v               >  0.0_WP, "cva (set): p > 0 violated")
    call assert(temp_atm%v%v        >  0.0_WP, "cva (set): temp_atm > 0 violated")
    call assert(len(trim(cv%label)) >       0, "cva (set): len(label) > 0 violated")
    call assert(cv%csa%v%v          >  0.0_WP, "cva (set): csa > 0 violated")
    call assert(cv%p_fs%v%v         >= 0.0_WP, "cva (set): p_fs >= 0 violated")
    call assert(cv%p_fd%v%v         >= 0.0_WP, "cva (set): p_fd >= 0 violated")
    call assert(cv%k%v%v            >= 0.0_WP, "cva (set): k >= 0 violated")
    call assert(cv%i_cv_mirror      >= 0,      "cva (set): i_cv_mirror >= 0 violated")
    
    call assert((cv%eos  >= 1) .and. (cv%eos <= MAX_EOS),      "cva (set): invalid EOS")
    call assert((cv%type >= 1) .and. (cv%type <= MAX_CV_TYPE), "cva (set): invalid control volume type")
    
    call assert_dimension(y, cv%gas)
    allocate(cv%m(size(y)))
    call y_sum%v%init_const(0.0_WP, n_d)
    call m_total%v%init_const(1.0_WP, n_d)
    do i = 1, size(y)
        call assert(gas(i)%gamma > 1.0_WP, "cva (set): gas%gamma > 1 violated")
        call assert(gas(i)%mm    > 0.0_WP, "cva (set): gas%mm > 0 violated")
        call assert(gas(i)%mm    < 0.1_WP, "cva (set): gas%mm < 0.1 violated") ! to catch using g/mol by mistake
        call assert(gas(i)%p_c   > 0.0_WP, "cva (set): gas%p_c > 0 violated")
        
        y_sum = y_sum + y(i)
    end do
    
    call assert(is_close(y_sum%v%v, 1.0_WP), "cva (set): mass fractions do not sum to 1")
    
    ! Get correct temperature depending on how the chamber is filled.
    if (isentropic_filling_) then
        call assert(present(p_atm), "cva (set): isentropic_filling = .true. requires p_atm")
        call assert(p_atm%v%v > 0.0_WP, "cva (set): p_atm > 0 required for isentropic_filling")
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
        call assert(p < cv%p_c(), "cva (set): ideal gas law validity is questionable")
    end if
end subroutine set

pure subroutine set_const(cv, label, csa, p_const, temp_const, gas, i_cv_mirror, type)
    class(cv_type), intent(in out) :: cv
    
    character(len=*), intent(in)     :: label       ! human-readable label for control volume
    type(si_area), intent(in)        :: csa         ! cross-sectional area
    type(si_pressure), intent(in)    :: p_const     ! pressure
    type(si_temperature), intent(in) :: temp_const  ! temperature
    type(gas_type), intent(in)       :: gas(:)      ! gas data
    integer, intent(in)              :: i_cv_mirror ! index of control volume to use in pressure difference calculation
    
    integer, intent(in), optional :: type ! type of CV to use
    
    integer :: n_d, n_gas, i_gas
    
    call assert(len(trim(cv%label)) > 0, "cva (set_const): len(label) > 0 violated")
    
    n_d   = size(p_const%v%d)
    n_gas = size(gas)
    
    allocate(cv%m(n_gas))
    do i_gas = 1, n_gas
        call cv%m(i_gas)%v%init_const(0.0_WP, n_d)
    end do
    
    call cv%x%v%init_const(1.0_WP, n_d)
    call cv%x_dot%v%init_const(0.0_WP, n_d)
    call cv%e%v%init_const(0.0_WP, n_d)
    call cv%rm_p%v%init_const(0.0_WP, n_d)
    call cv%p_fs%v%init_const(0.0_WP, n_d)
    call cv%p_fd%v%init_const(0.0_WP, n_d)
    call cv%k%v%init_const(0.0_WP, n_d)
    call cv%x_z%v%init_const(0.0_WP, n_d)
    call cv%x_stop%v%init_const(X_STOP_DEFAULT, n_d)
    
    call assert(.not. is_close(cv%x%v%v, cv%x_stop%v%v), "cva (set_const): x = x_stop will cause immediate termination of run")
    
    cv%label       = label
    cv%csa         = csa
    cv%eos         = CONST_EOS
    cv%gas         = gas
    cv%i_cv_mirror = i_cv_mirror
    cv%p_const     = p_const
    cv%temp_const  = temp_const
    
    call assert(cv%p_const%v%v    > 0.0_WP, "cva (set_const): p_const > 0 violated")
    call assert(cv%temp_const%v%v > 0.0_WP, "cva (set_const): temp_const > 0 violated")
    
    if (present(type)) then
        cv%type = type
    else
        cv%type = MIRROR_CV_TYPE
    end if
    
    select case (cv%type)
        case (NORMAL_CV_TYPE)
            call assert(cv%i_cv_mirror >= 0, "cva (set_const): i_cv_mirror >= 0 violated")
        case (MIRROR_CV_TYPE)
            call assert(cv%i_cv_mirror >= 1, "cva (set_const, MIRROR_CV_TYPE): i_cv_mirror >= 1 violated")
        case default
            error stop "cva (set_const): invalid cv%type"
    end select
end subroutine set_const

pure function p_f(cv, p_fe)
    ! Returns pressure of friction.
    
    class(cv_type), intent(in)    :: cv
    type(si_pressure), intent(in) :: p_fe ! equilibrium pressure
    
    type(si_pressure) :: p_f
    
    type(si_velocity) :: v_scale
    
    call v_scale%v%init_const(0.1_WP, size(cv%x%v%d))
    
    call assert(cv%p_fs%v%v >= 0.0_WP, "cva (p_f): cv%p_fs%v > 0 violated")
    call assert(cv%p_fd%v%v >= 0.0_WP, "cva (p_f): cv%p_fd%v > 0 violated")
    
    p_f = p_f0(cv, p_fe) + (cv%p_fd - tanh(cv%x_dot/v_scale)*p_f0(cv, p_fe))*tanh(cv%x_dot/v_scale)
    
    call assert(p_f%v%v <= max(cv%p_fs%v%v, cv%p_fd%v%v), "cva (p_f): p_f <= max(p_fs, p_fd) violated")
end function p_f

pure function p_f0(cv, p_fe)
    ! Returns actual static pressure of friction.
    ! `p_fs` is the *maximum* static pressure of friction.
    
    class(cv_type), intent(in)    :: cv
    type(si_pressure), intent(in) :: p_fe ! equilibrium pressure
    
    type(si_pressure) :: p_f0, p_s
    
    call assert(cv%p_fs%v%v >= 0.0_WP, "cva (p_f0): cv%p_fs%v > 0 violated")
    
    p_s = 0.1_WP*cv%p_fs ! TODO: make a function of `dt`
    call assert(p_s <= cv%p_fs, "cva (p_f0): p_s <= p_fs violated")
    
    if (p_fe <= -p_s) then
        p_f0 = -p_f0_high(p_fe, cv%p_fs, p_s)
        
        call assert(p_f0 <= -p_s, "cva (p_f0), first branch: p_f0 <= -p_s violated")
    else if (p_fe <= p_s) then
        p_f0 = p_fe
    else
        p_f0 = p_f0_high(p_fe, cv%p_fs, p_s)
        
        call assert(p_f0 >= p_s, "cva (p_f0), third branch: p_f0 >= p_s violated")
    end if
    
    call assert(p_f0 >= -cv%p_fs, "cva (p_f0): p_f0 >= -p_fs violated")
    call assert(p_f0 <= cv%p_fs, "cva (p_f0): p_f0 <= p_fs violated")
    
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

pure function d_x_d_t(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_velocity) :: d_x_d_t
    
    call assert(sys%cv(i_cv)%i_cv_mirror >= 0, "cva (d_x_d_t): i_cv_mirror must be a positive integer or zero")
    call assert(sys%cv(i_cv)%i_cv_mirror /= i_cv, "cva (d_x_d_t): i_cv_mirror can not equal i_cv")
    
    select case (sys%cv(i_cv)%type)
        case (NORMAL_CV_TYPE)
            d_x_d_t = sys%cv(i_cv)%x_dot
        case (MIRROR_CV_TYPE)
            call assert(sys%cv(i_cv)%i_cv_mirror >= 1, "cva (d_xdot_d_t): i_cv_mirror must be defined for a mirror CV")
            call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == NORMAL_CV_TYPE, &
                            "cva (d_x_d_t): mirror CV is not a NORMAL_CV_TYPE")
            d_x_d_t = -sys%cv(sys%cv(i_cv)%i_cv_mirror)%x_dot
        case default
            error stop "cva (d_x_d_t): invalid cv%type"
    end select
end function d_x_d_t

pure function d_xdot_d_t(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_acceleration) :: d_xdot_d_t
    
    type(si_area)     :: csa_i_mirror
    character(len=3)  :: i_cv_string, i_cv_mirror_string
    
    call assert(sys%cv(i_cv)%i_cv_mirror >= 0, "cva (d_xdot_d_t): i_cv_mirror must be a positive integer or zero")
    call assert(sys%cv(i_cv)%i_cv_mirror /= i_cv, "cva (d_xdot_d_t): i_cv_mirror can not equal i_cv")
    
    if (sys%cv(i_cv)%i_cv_mirror >= 1) then
        csa_i_mirror = sys%cv(sys%cv(i_cv)%i_cv_mirror)%csa
        write(unit=i_cv_string, fmt="(i0)") i_cv
        write(unit=i_cv_mirror_string, fmt="(i0)") sys%cv(i_cv)%i_cv_mirror
        call assert(is_close(sys%cv(i_cv)%csa%v%v, csa_i_mirror%v%v), &
                        "cva (d_xdot_d_t): CSA for mirror CV (" // trim(i_cv_mirror_string) &
                        // ") assumed equal to CSA for current CV (" // trim(i_cv_string) // ")")
    end if
    
    select case (sys%cv(i_cv)%type)
        case (NORMAL_CV_TYPE)
            d_xdot_d_t = d_xdot_d_t_normal(sys, i_cv)
        case (MIRROR_CV_TYPE)
            call assert(sys%cv(i_cv)%i_cv_mirror >= 1, "cva (d_xdot_d_t): i_cv_mirror must be defined for a mirror CV")
            call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == NORMAL_CV_TYPE, &
                            "cva (d_xdot_d_t): mirror CV is not a NORMAL_CV_TYPE")
            d_xdot_d_t = -d_xdot_d_t_normal(sys, sys%cv(i_cv)%i_cv_mirror)
        case default
            error stop "cva (d_xdot_d_t): invalid cv%type"
    end select
end function d_xdot_d_t

pure function d_xdot_d_t_normal(sys, i_cv)
    type(cv_system_type), intent(in) :: sys
    integer, intent(in)              :: i_cv
    
    type(si_acceleration) :: d_xdot_d_t_normal
    
    type(si_pressure) :: p_fe ! friction pressure at equilibrium ($\partial \dot{x}/\partial t = 0$)
    type(si_pressure) :: p_other
    
    call assert(sys%cv(i_cv)%csa%v%v > 0.0_WP, "cva (d_xdot_d_t_normal): cv%csa > 0 violated")
    call assert(sys%cv(i_cv)%type == NORMAL_CV_TYPE, "cva (d_xdot_d_t_normal): CV needs to be NORMAL_CV_TYPE")
    
    if (sys%cv(i_cv)%i_cv_mirror >= 1) then
        p_other = sys%cv(sys%cv(i_cv)%i_cv_mirror)%p()
        
        call assert(sys%cv(sys%cv(i_cv)%i_cv_mirror)%type == MIRROR_CV_TYPE, &
                        "cva (d_xdot_d_t_normal): mirror CV not MIRROR_CV_TYPE")
    else
        ! If `i_cv_mirror == 0` then there is no mirror CV.
        call p_other%v%init_const(0.0_WP, size(sys%cv(i_cv)%csa%v%d))
    end if
    
    p_fe = sys%cv(i_cv)%p() - p_other - (sys%cv(i_cv)%k/sys%cv(i_cv)%csa)*(sys%cv(i_cv)%x - sys%cv(i_cv)%x_z)
    
    d_xdot_d_t_normal = sys%cv(i_cv)%csa*sys%cv(i_cv)%rm_p &
                            * (sys%cv(i_cv)%p() - p_other - sys%cv(i_cv)%p_f(p_fe)) &
                            - sys%cv(i_cv)%k*sys%cv(i_cv)%rm_p*(sys%cv(i_cv)%x - sys%cv(i_cv)%x_z)
end function d_xdot_d_t_normal

pure function d_m_k_d_t(cv, m_dot, k_gas, i_cv)
    type(cv_type), intent(in)           :: cv
    type(si_mass_flow_rate), intent(in) :: m_dot(:, :)
    integer, intent(in)                 :: k_gas, i_cv
    
    type(si_mass_flow_rate) :: d_m_k_d_t
    
    integer :: n_d, n_cv, j_cv
    type(si_mass)  :: m_total
    type(unitless) :: y_k
    
    call assert(size(m_dot, 1) == size(m_dot, 2), "cva (d_m_k_d_t): m_dots must be square")
    call assert_dimension(m_dot(1, 1)%v%d, cv%x%v%d)
    
    n_d = size(cv%x%v%d)
    call d_m_k_d_t%v%init_const(0.0_WP, n_d)
    
    n_cv = size(m_dot, 1)
    m_total = cv%m_total()
    if (.not. is_close(m_total%v%v, 0.0_WP)) then
        y_k = cv%m(k_gas) / m_total
    else
        call y_k%v%init_const(0.0_WP, n_d)
    end if
    
    do j_cv = 1, n_cv
        call assert(is_close(m_dot(j_cv, j_cv)%v%v, 0.0_WP), "cva (d_m_k_d_t): mass can not flow from self to self")
        d_m_k_d_t = d_m_k_d_t + y_k*m_dot(j_cv, i_cv) - y_k*m_dot(i_cv, j_cv)
    end do
end function d_m_k_d_t

pure function d_e_d_t(cv, h_dot, i_cv)
    type(cv_type), intent(in)             :: cv
    type(si_energy_flow_rate), intent(in) :: h_dot(:, :)
    integer, intent(in)                   :: i_cv
    
    type(si_energy_flow_rate) :: d_e_d_t
    
    integer :: n_cv, j_cv
    
    call assert(size(h_dot, 1) == size(h_dot, 2), "cva (d_e_d_t): h_dots must be square")
    call assert_dimension(h_dot(1, 1)%v%d, cv%x%v%d)
    call assert(cv%csa%v%v > 0.0_WP, "cva (d_e_d_t): cv%csa > 0 violated")
    
    d_e_d_t = -cv%p() * cv%csa * cv%x_dot
    
    n_cv = size(h_dot, 1)
    do j_cv = 1, n_cv
        call assert(is_close(h_dot(j_cv, j_cv)%v%v, 0.0_WP), "cva (d_e_d_t): energy can not flow from self to self")
        d_e_d_t = d_e_d_t + h_dot(j_cv, i_cv) - h_dot(i_cv, j_cv)
    end do
end function d_e_d_t

pure subroutine assert_mass(cv, procedure_name)
    ! Why not make this a type-bound operator?
    ! That would make my assertion counting Python program not count these.
    
    type(cv_type), intent(in)    :: cv
    character(len=*), intent(in) :: procedure_name
    
    integer       :: i
    type(si_mass) :: m_total
    
    call assert(len(trim(procedure_name)) > 0, "cva (assert_mass): procedure name should not be empty")
    
    do i = 1, size(cv%m)
        call assert(cv%m(i)%v%v >= 0.0_WP, "cva (" // trim(procedure_name) // "): cv%m >= 0 violated")
    end do
    
    m_total = cv%m_total()
    select case (cv%eos)
        case (IDEAL_EOS)
            call assert(m_total%v%v > 0.0_WP, "cva (" // trim(procedure_name) // ", IDEAL_EOS): cv%m_total > 0 violated")
        case (CONST_EOS)
            call assert(m_total%v%v >= 0.0_WP, "cva (" // trim(procedure_name) // ", CONST_EOS): cv%m_total >= 0 violated")
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
    
    call assert(p_r%v%v >= 0.0_WP, "cva (f_m_dot): p_r >= 0 violated")
    call assert(p_r%v%v <= 1.0_WP, "cva (f_m_dot): p_r <= 1 violated")
    
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
    
    call assert(p_r%v%v >= 0.0_WP, "cva (g_m_dot): p_r >= 0 violated")
    call assert(p_r%v%v <= 1.0_WP, "cva (g_m_dot): p_r <= 1 violated")
    
    if (p_r%v%v < P_RL) then
        call g_m_dot%v%init_const(0.0_WP, size(p_r%v%d))
    else
        g_m_dot = 2.0_WP*((1.0_WP - p_r) / (1.0_WP - P_RL))*((1.0_WP - p_r) / (1.0_WP - P_RL))*((1.0_WP - p_r) / (1.0_WP - P_RL)) &
                    - 3.0_WP*((1.0_WP - p_r) / (1.0_WP - P_RL))*((1.0_WP - p_r) / (1.0_WP - P_RL)) + 1.0_WP
    end if
    
    call assert(g_m_dot%v%v >= 0.0_WP, "cva (g_m_dot): g_m_dot >= 0 violated")
    call assert(g_m_dot%v%v <= 1.0_WP, "cva (g_m_dot): g_m_dot <= 1 violated")
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
        call assert(p_r%v%v >= 0.0_WP, "cva (m_dot): p_r >= 0 violated")
        call assert(p_r%v%v <= 1.0_WP, "cva (m_dot): p_r <= 1 violated")
        
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
    
    call assert(size(sys%cv) == size(sys%con, 1), "cva (calculate_flows): inconsistent sys%cv and sys%con sizes")
    
    n_cv = size(sys%cv)
    call assert(n_cv > 1, "cva (calculate_flows): there needs to be at least 2 control volumes to have flows")
    
    allocate(m_dot(n_cv, n_cv))
    allocate(h_dot(n_cv, n_cv))
    do i_from_cv = 1, n_cv
        do i_to_cv = 1, n_cv
            if (i_from_cv == i_to_cv) call assert(.not. sys%con(i_from_cv, i_to_cv)%active, &
                                                    "cva (calculate_flows): can't flow from self to self")
            
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
        do k_gas = 1, n_gas
            sys%cv(i_cv)%m(k_gas) = sys_old%cv(i_cv)%m(k_gas) + a*cv_delta_in(i_cv)%m(k_gas)
        end do
    end do
    call sys%calculate_flows(m_dot, h_dot)
    do i_cv = 1, n_cv
        cv_delta_out(i_cv)%x     = dt*d_x_d_t(sys, i_cv)
        cv_delta_out(i_cv)%x_dot = dt*d_xdot_d_t(sys, i_cv)
        cv_delta_out(i_cv)%e     = dt*d_e_d_t(sys%cv(i_cv), h_dot, i_cv)
        do k_gas = 1, n_gas
            cv_delta_out(i_cv)%m(k_gas) = dt*d_m_k_d_t(sys%cv(i_cv), m_dot, k_gas, i_cv)
        end do
    end do
end subroutine rk_stage

subroutine run(sys_start, sys_end, status, t_stop)
    type(cv_system_type), allocatable, intent(in)  :: sys_start
    type(cv_system_type), allocatable, intent(out) :: sys_end
    type(run_status_type), intent(out)             :: status
    type(si_time), intent(in), optional            :: t_stop
    
    type(cv_system_type), allocatable :: sys_old, sys_new, sys_temp
    
    type(si_time)        :: t_stop_ ! time where simulation will stop
    integer              :: n_d
    type(si_time)        :: t, dt, t_old
    logical              :: exit_time_loop
    
    n_d = size(sys_start%cv(1)%x%v%d)
    
    if (present(t_stop)) then
        t_stop_ = t_stop
    else
        call t_stop_%v%init_const(T_STOP_DEFAULT, n_d)
    end if
    
    sys_old = sys_start
    call t%v%init_const(0.0_WP, n_d)
    call dt%v%init_const(DT_DEFAULT, n_d)
    
    time_loop: do
        call time_step(sys_old, dt, sys_new)
        t_old = t
        t     = t + dt
        
        !print *, t%v%v
        
        call check_sys(sys_new, t_stop_, t, status, exit_time_loop)
        if (exit_time_loop) exit time_loop
        
        call move_alloc(from=sys_old,  to=sys_temp)
        call move_alloc(from=sys_new,  to=sys_old)
        call move_alloc(from=sys_temp, to=sys_new)
    end do time_loop
    
    if (status%rc == 0) then
        ! If successful, interpolate to correct `x` value.
        call sys_interp(t_old, dt, status%i_cv(1), sys_old, sys_new, t, sys_end)
    else
        sys_end = sys_new
    end if
    
    status%t = t
end subroutine run

pure subroutine check_sys(sys, t_stop, t, status, exit_time_loop)
    type(cv_system_type), allocatable, intent(in) :: sys
    type(si_time), intent(in)                     :: t_stop, t
    type(run_status_type), intent(out)            :: status
    logical, intent(out)                          :: exit_time_loop
    
    integer              :: n_cv, i_cv, j_cv, n_bad_cv
    type(si_mass)        :: m_total_i, m_total_j
    type(si_temperature) :: temp_i, temp_j
    
    n_cv = size(sys%cv)
    exit_time_loop = .false.
    
    do i_cv = 1, n_cv
        if (sys%cv(i_cv)%x >= sys%cv(i_cv)%x_stop) then
            status%rc = 0
            allocate(status%i_cv(1))
            status%i_cv(1) = i_cv
            
            exit_time_loop = .true.
            return
        end if
        
        m_total_i = sys%cv(i_cv)%m_total()
        if (m_total_i%v%v < 0.0_WP) then
            status%rc = MASS_RUN_RC
            
            allocate(status%data(n_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                m_total_j = sys%cv(j_cv)%m_total()
                if (m_total_j%v%v < 0.0_WP) n_bad_cv = n_bad_cv + 1
                status%data(j_cv) = m_total_j%v%v
            end do
            call assert(n_bad_cv >= 1, "cva (check_sys): number of control volumes with negative mass should be >= 1")
            
            allocate(status%i_cv(n_bad_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                n_bad_cv = n_bad_cv + 1
                status%i_cv(n_bad_cv) = j_cv
            end do
            
            exit_time_loop = .true.
            return
        end if
        
        temp_i = sys%cv(i_cv)%temp()
        if (temp_i%v%v <= 0.0_WP) then
            status%rc = TEMP_RUN_RC
            
            allocate(status%data(n_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                temp_j = sys%cv(j_cv)%temp()
                if (temp_j%v%v < 0.0_WP) n_bad_cv = n_bad_cv + 1
                status%data(j_cv) = temp_j%v%v
            end do
            call assert(n_bad_cv >= 1, &
                "cva (check_sys): number of control volumes with negative or zero temperature should be >= 1")
            
            allocate(status%i_cv(n_bad_cv))
            n_bad_cv = 0
            do j_cv = 1, n_cv
                n_bad_cv = n_bad_cv + 1
                status%i_cv(n_bad_cv) = j_cv
            end do
            
            exit_time_loop = .true.
            return
        end if
    end do
    
    if (t >= t_stop) then
        status%rc = TIMEOUT_RUN_RC
        exit_time_loop = .true.
        return
    end if
end subroutine check_sys

pure subroutine sys_interp(t_old, dt, i_cv_interp, sys_old, sys_new, t, sys_end)
    ! Linearly interpolate system to state between two systems based on `x_stop` for control volume number `i_cv_interp`.
    
    type(si_time), intent(in)                      :: t_old, dt
    integer, intent(in)                            :: i_cv_interp ! control volume number to interpolate based on
    type(cv_system_type), allocatable, intent(in)  :: sys_old, sys_new
    type(si_time), intent(out)                     :: t
    type(cv_system_type), allocatable, intent(out) :: sys_end
    
    type(unitless) :: frac
    integer        :: i_cv, n_cv, k_gas, n_gas
    
    call assert(is_close(sys_old%cv(i_cv_interp)%x_stop%v%v, sys_new%cv(i_cv_interp)%x_stop%v%v), &
                    "cva (sys_interp): x_stop is inconsistent")
    
    n_cv  = size(sys_old%cv)
    n_gas = size(sys_old%cv(1)%m)
    
    frac = (sys_new%cv(i_cv_interp)%x_stop - sys_old%cv(i_cv_interp)%x) &
                / (sys_new%cv(i_cv_interp)%x - sys_old%cv(i_cv_interp)%x)
    
    call assert(frac%v%v >= 0.0_WP, "cva (sys_interp): fraction >= 0 violated")
    call assert(frac%v%v <= 1.0_WP, "cva (sys_interp): fraction <= 1 violated")
    
    sys_end = sys_old
    do i_cv = 1, n_cv
        sys_end%cv(i_cv)%x     = sys_old%cv(i_cv)%x     + frac*(sys_new%cv(i_cv)%x     - sys_old%cv(i_cv)%x)
        sys_end%cv(i_cv)%x_dot = sys_old%cv(i_cv)%x_dot + frac*(sys_new%cv(i_cv)%x_dot - sys_old%cv(i_cv)%x_dot)
        sys_end%cv(i_cv)%e     = sys_old%cv(i_cv)%e     + frac*(sys_new%cv(i_cv)%e     - sys_old%cv(i_cv)%e)
        do k_gas = 1, n_gas
            sys_end%cv(i_cv)%m(k_gas) = sys_old%cv(i_cv)%m(k_gas) + frac*(sys_new%cv(i_cv)%m(k_gas) - sys_old%cv(i_cv)%m(k_gas))
        end do
    end do
    t = t_old + frac*dt
end subroutine sys_interp

end module cva
