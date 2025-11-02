! Module for control volume analysis.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module cva

use prec, only: WP
use units, only: si_length       => unit_p10_p00_p00_p00_p00, &
                 si_velocity     => unit_p10_p00_m10_p00_p00, &
                 si_mass         => unit_p00_p10_p00_p00_p00, &
                 si_energy       => unit_p20_p10_m20_p00_p00, &
                 si_area         => unit_p20_p00_p00_p00_p00, &
                 si_inverse_mass => unit_p00_m10_p00_p00_p00, &
                 si_pressure     => unit_m10_p10_m20_p00_p00, &
                 unitless        => unit_p00_p00_p00_p00_p00, &
                 si_stiffness    => unit_p00_p10_m20_p00_p00
implicit none
private

public :: smooth_min
public :: f_m_dot, g_m_dot, m_dot
public :: p_v_h2o

! <https://en.wikipedia.org/wiki/Gas_constant>
real(WP), public, parameter :: R_BAR = 8.31446261815324_WP ! J/(mol*K)

real(WP), public, parameter :: TEMP_C_TO_K = 273.15_WP ! K, temperature to add to convert from C to K

! <https://en.wikipedia.org/wiki/Density_of_air>
real(WP), public, parameter :: P_ATM   = 101325.0_WP           ! Pa
real(WP), public, parameter :: T_ATM   = TEMP_C_TO_K + 15.0_WP ! K
real(WP), public, parameter :: RHO_ATM = 1.2250_WP             ! kg/m3

! pressure ratio laminar flow nominally starts at
! based on first part of beater_pneumatic_2007 eq. 5.4
real(WP), public, parameter :: P_RL = 0.999_WP ! unitless

real(WP), public, parameter :: TEMP_0 = 300.0_WP ! K, temperature that `gamma`, `u_0`, and `h_0` are taken at in `gas_type`

type, public :: gas_type
    real(WP) :: gamma ! ratio of specific heats, unitless
    real(WP) :: u_0   ! internal energy at `TEMP_0`, J/kg
    real(WP) :: h_0   ! enthalphy at `TEMP_0`, J/kg
    real(WP) :: mm    ! molar mass, kg/mol
    real(WP) :: p_c   ! critical pressure, Pa
contains
    procedure :: u => u_gas
    procedure :: h => h_gas
    procedure, private :: r   => r_gas
    procedure, private :: c_v => c_v_gas
    procedure, private :: c_p => c_p_gas
end type gas_type

! Molecular mass and critical pressure of air (consistent with dry air): moran_fundamentals_2008 table A-1
! Specific heat ratio of air (consistent with dry air): moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of air: moran_fundamentals_2008 table A-22
type(gas_type), public, parameter :: DRY_AIR = gas_type(gamma = 1.400_WP, &
                                                        u_0   = 214.07e3_WP, &
                                                        h_0   = 300.19e3_WP, &
                                                        mm    = 28.97e-3_WP, &
                                                        p_c   = 37.7e5_WP)

! Molecular mass and critical pressure of N2: moran_fundamentals_2008 table A-1
! Specific heat ratio of N2: moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of N2: moran_fundamentals_2008 table A-23
type(gas_type), public, parameter :: N2      = gas_type(gamma = 1.400_WP, &
                                                        u_0   = 28.01e-3_WP*6229.0e3_WP, &
                                                        h_0   = 28.01e-3_WP*8723.0e3_WP, &
                                                        mm    = 28.01e-3_WP, &
                                                        p_c   = 33.9e5_WP)

! Molecular mass and critical pressure of O2: moran_fundamentals_2008 table A-1
! Specific heat ratio of O2: moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of O2: moran_fundamentals_2008 table A-23
type(gas_type), public, parameter :: O2      = gas_type(gamma = 1.395_WP, &
                                                        u_0   = 6242.0e3_WP, &
                                                        h_0   = 8736.0e3_WP, &
                                                        mm    = 32.0e-3_WP, &
                                                        p_c   = 50.5e5_WP)

! Molecular mass and critical pressure of Ar: moran_fundamentals_2008 table A-1
! `gamma`, `u_0`, `h_0` from <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/Ar> ("Fluid Properties")
type(gas_type), public, parameter :: AR      = gas_type(gamma = 0.52154_WP/0.31239_WP, &
                                                        u_0   = 155.90e3_WP, &
                                                        h_0   = 93.497e3_WP, &
                                                        mm    = 39.94e-3_WP, &
                                                        p_c   = 48.6e5_WP)

! Molecular mass and critical pressure of CO2: moran_fundamentals_2008 table A-1
! Specific heat ratio of CO2: moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of CO2: moran_fundamentals_2008 table A-23
type(gas_type), public, parameter :: CO2     = gas_type(gamma = 1.288_WP, &
                                                        u_0   = 44.01e-3_WP*6939.0e3_WP, &
                                                        h_0   = 44.01e-3_WP*9431.0e3_WP, &
                                                        mm    = 44.01e-3_WP, &
                                                        p_c   = 73.9e5_WP)

! Molecular mass and critical pressure of H2O: moran_fundamentals_2008 table A-1
! Specific heat ratio of H2O: <https://en.wikipedia.org/wiki/Heat_capacity_ratio>
! Internal energy and enthalpy of H2O: moran_fundamentals_2008 table A-23
! TODO: Switch to stream tables for `u` and `h` if adding liquid water; see moran_thermodynamics_2008 pp. 666--667
type(gas_type), public, parameter :: H2O     = gas_type(gamma = 1.330_WP, & ! at 20 C, not 300 K, but close enough
                                                        u_0   = 18.02e-3_WP*7472.0e3_WP, &
                                                        h_0   = 18.02e-3_WP*9966.0e3_WP, &
                                                        mm    = 18.02e-3_WP, &
                                                        p_c   = 220.9e5_WP)

type, public :: cv_type ! control volume
    ! time varying
    type(si_length)            :: x     ! location of piston/projectile
    type(si_velocity)          :: x_dot ! velocity of piston/projectile
    type(si_mass), allocatable :: m(:)  ! mass(es) of gas(es) in control volume
    type(si_energy)            :: e     ! energy of gas in control volume
    
    ! constants
    type(si_area)               :: csa        ! cross-sectional area
    type(si_inverse_mass)       :: rm_p       ! reciprocal mass of piston/projectile
    type(si_pressure)           :: p_fs, p_fd ! static and dynamic friction pressure
    type(si_pressure)           :: p_atm      ! atmospheric pressure
    type(si_stiffness)          :: k          ! stiffness of spring attached to piston
    type(si_length)             :: x_z        ! zero force location for spring
    type(gas_type), allocatable :: gas(:)     ! gas data
contains
    procedure :: m_total
    procedure :: p_eos
    procedure :: rho_eos
    procedure :: p_c
    procedure :: y
    procedure :: chi
    procedure :: r    => r_cv
    procedure :: temp => temp_cv
    procedure :: vol  => vol_cv
    procedure :: rho  => rho_cv
    procedure :: p    => p_cv
    procedure :: set
    procedure :: p_f
    procedure :: p_f0
end type cv_type

type, public :: con_type ! connection between control volumes
    logical        :: active
    type(si_area)  :: a_e ! effective area
    type(unitless) :: b   ! critical pressure ratio
contains
    procedure :: m_dot
end type con_type

type, public :: cv_system_type
    type(cv_type), allocatable  :: cvs(:)
    type(con_type), allocatable :: cons(:, :)
end type cv_system_type

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!
! methods for `gas_type` !
!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function u_gas(gas, temp)
    ! Constant specific heats assumed for now. Will improve later.
    
    use units, only: si_temperature     => unit_p00_p00_p00_p10_p00, &
                     si_specific_energy => unit_p20_p00_m20_p00_p00
    
    class(gas_type), intent(in)      :: gas
    type(si_temperature), intent(in) :: temp
    
    type(si_specific_energy) :: u_gas, u_0
    type(si_temperature)     :: temp_0_
    
    integer :: n_d ! number of derivatives
    
    n_d = size(temp%v%d)
    call temp_0_%v%init_const(TEMP_0, n_d)
    call u_0%v%init_const(gas%u_0, n_d)
    u_gas = u_0 + gas%c_v(n_d) * (temp - temp_0_)
end function u_gas

pure function h_gas(gas, temp)
    ! Constant specific heats assumed for now. Will improve later.
    
    use units, only: si_temperature     => unit_p00_p00_p00_p10_p00, &
                     si_specific_energy => unit_p20_p00_m20_p00_p00
    
    class(gas_type), intent(in)      :: gas
    type(si_temperature), intent(in) :: temp
    
    type(si_specific_energy) :: h_gas, h_0
    type(si_temperature)     :: temp_0_
    
    integer :: n_d ! number of derivatives
    
    n_d = size(temp%v%d)
    call temp_0_%v%init_const(TEMP_0, n_d)
    call h_0%v%init_const(gas%h_0, n_d)
    h_gas = h_0 + gas%c_p(n_d) * (temp - temp_0_)
end function h_gas

pure function r_gas(gas, n_d)
    ! Gas constant for a *pure* gas.
    
    use units, only: si_specific_heat => unit_p20_p00_m20_m10_p00
    
    class(gas_type), intent(in) :: gas
    integer, intent(in)         :: n_d ! number of derivatives
    
    type(si_specific_heat) :: r_gas
    
    call r_gas%v%init_const(R_BAR/gas%mm, n_d)
end function r_gas

pure function c_v_gas(gas, n_d)
    use units, only: si_specific_heat => unit_p20_p00_m20_m10_p00
    use checks, only: assert
    
    class(gas_type), intent(in) :: gas
    integer, intent(in)         :: n_d ! number of derivatives
    
    type(si_specific_heat) :: c_v_gas
    
    call assert(gas%gamma > 1.0_WP, "cva (temp_cv): gamma > 1 violated")
    
    ! moran_fundamentals_2008 eq. 3.47b, p. 119
    c_v_gas = gas%r(n_d) / (gas%gamma - 1.0_WP)
end function c_v_gas

pure function c_p_gas(gas, n_d)
    use units, only: si_specific_heat => unit_p20_p00_m20_m10_p00
    use checks, only: assert
    
    class(gas_type), intent(in) :: gas
    integer, intent(in)         :: n_d ! number of derivatives
    
    type(si_specific_heat) :: c_p_gas
    
    call assert(gas%gamma > 1.0_WP, "cva (temp_cv): gamma > 1 violated")
    
    ! moran_fundamentals_2008 eq. 3.47a, p. 119
    c_p_gas = gas%gamma * gas%r(n_d) / (gas%gamma - 1.0_WP)
end function c_p_gas

!!!!!!!!!!!!!!!!!!!!!!!!!
! methods for `cv_type` !
!!!!!!!!!!!!!!!!!!!!!!!!!

pure function m_total(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_mass) :: m_total
    
    integer :: i
    
    call m_total%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        m_total = m_total + cv%m(i)
    end do
end function m_total

pure function p_eos(cv, rho, temp)
    ! Calculate pressure using the equation of state.
    
    use units, only: si_mass_density  => unit_m30_p10_p00_p00_p00, &
                     si_temperature   => unit_p00_p00_p00_p10_p00, &
                     si_specific_heat => unit_p20_p00_m20_m10_p00
    use checks, only: assert, assert_dimension
    
    class(cv_type), intent(in)        :: cv
    type(si_mass_density), intent(in) :: rho
    type(si_temperature), intent(in)  :: temp
    type(si_pressure)                 :: p_eos
    
    integer :: n_d ! number of derivatives
    
    call assert(rho%v%v  > 0.0_WP, "cva (p_eos): rho%v > 0 violated")
    call assert(temp%v%v > 0.0_WP, "cva (p_eos): temp%v > 0 violated")
    call assert_dimension(rho%v%d, temp%v%d)
    
    n_d   = size(rho%v%d)
    p_eos = rho * cv%r() * temp
    
    call assert(p_eos%v%v > 0.0_WP, "cva (p_eos): p_eos%v > 0 violated")
    
    call assert(p_eos < cv%p_c(), "cva (p_eos): ideal gas law validity is questionable")
end function p_eos

pure function rho_eos(cv, p, temp, y)
    ! Calculate density using the equation of state.
    ! In the future, if an EOS more complex than the ideal gas law is used, it might make sense to calculate `rho_eos` from `p_eos`
    ! with the Newton method.
    
    use units, only: si_mass_density  => unit_m30_p10_p00_p00_p00, &
                     si_temperature   => unit_p00_p00_p00_p10_p00, &
                     si_specific_heat => unit_p20_p00_m20_m10_p00
    use checks, only: assert, assert_dimension, is_close
    
    class(cv_type), intent(in)       :: cv
    type(si_pressure), intent(in)    :: p
    type(si_temperature), intent(in) :: temp
    type(unitless), intent(in)       :: y(:)
    
    type(si_mass_density) :: rho_eos
    
    integer                :: n_d, i, j
    type(unitless)         :: mm ! molar mass (mol/kg), not actually unitless!
    type(si_specific_heat) :: r_bar_ ! universal gas constant, which also doesn't have the same units as specific heat!
    type(unitless)         :: denominator, y_sum
    type(si_specific_heat) :: r_cv
    
    ! Don't check `p_c` here as that requires the masses.
    
    call assert(p%v%v    > 0.0_WP, "cva (rho_eos): p%v > 0 violated")
    call assert(temp%v%v > 0.0_WP, "cva (rho_eos): temp%v > 0 violated")
    call assert(size(y)  >= 1,     "cva (rho_eos): size(y) >= 1 violated")
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
            denominator = denominator + y(j)*(cv%gas(i)%mm/cv%gas(j)%mm)
        end do
        mm    = mm + y(i)*cv%gas(i)%mm/denominator
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
    
    use checks, only: assert, is_close
    
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
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm)
        end do
        chi = cv%m(i) / denominator
        call p_ci%v%init_const(cv%gas(i)%p_c, size(cv%m(1)%v%d))
        p_c = p_c + chi*p_ci
        chi_sum = chi_sum + chi
    end do
    
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (p_c): chi does not sum to 1")
end function p_c

pure function y(cv)
    ! mass fractions
    
    use checks, only: assert, is_close
    
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
    
    use checks, only: assert, is_close
    
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
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm)
            call assert(denominator%v%v >= 0.0_WP, "cva (chi): denominator >= 0 violated")
        end do
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
    
    use units, only: si_specific_heat => unit_p20_p00_m20_m10_p00
    use checks, only: assert, is_close
    
    class(cv_type), intent(in) :: cv
    
    type(si_specific_heat) :: r_cv
    
    integer        :: n_d, i, j
    type(unitless) :: chi ! mole fraction for gas in loop
    type(unitless) :: mm  ! molar mass (mol/kg), not actually unitless!
    type(si_mass)  :: denominator
    type(si_specific_heat) :: r_bar_ ! universal gas constant, which also doesn't have the same units as specific heat!
    type(unitless) :: chi_sum
    
    n_d = size(cv%m(1)%v%d)
    
    call assert_mass(cv, "r_cv")
    
    ! <https://en.wikipedia.org/wiki/Molar_mass#Average_molar_mass_of_mixtures>
    call mm%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    call chi_sum%v%init_const(0.0_WP, size(cv%m(1)%v%d))
    do i = 1, size(cv%m)
        call denominator%v%init_const(0.0_WP, size(cv%m(1)%v%d))
        do j = 1, size(cv%m)
            denominator = denominator + cv%m(j)*(cv%gas(i)%mm/cv%gas(j)%mm)
        end do
        chi = cv%m(i) / denominator
        chi_sum = chi_sum + chi
        mm  = mm + chi*cv%gas(i)%mm
    end do
    
    call r_bar_%v%init_const(R_BAR, n_d)
    r_cv = r_bar_ / mm
    
    call assert(r_cv%v%v > 0.0_WP, "cva (r_cv): r_cv > 0 violated")
    call assert(is_close(chi_sum%v%v, 1.0_WP), "cva (r_cv): chi does not sum to 1")
end function r_cv

pure function temp_cv(cv)
    use units, only: si_temperature     => unit_p00_p00_p00_p10_p00, &
                     si_heat_capacity   => unit_p20_p10_m20_m10_p00, &
                     si_specific_energy => unit_p20_p00_m20_p00_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_temperature) :: temp_cv
    
    integer :: n_d, i
    type(si_energy)          :: e_0
    type(si_specific_energy) :: u_0
    type(si_heat_capacity)   :: heat_capacity
    type(si_temperature)     :: temp_0_
    
    ! Constant specific heats assumed for now. Will improve later.
    
    call assert_mass(cv, "temp_cv")
    
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
    
    call assert(temp_cv%v%v > 0.0_WP, "cva (temp_cv): temp_cv > 0 violated")
end function temp_cv

pure function vol_cv(cv)
    use units, only: si_volume => unit_p30_p00_p00_p00_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_volume) :: vol_cv
    
    call assert(cv%x%v%v > 0.0_WP, "cva (vol_cv): cv%x > 0 violated")
    call assert(cv%csa%v%v > 0.0_WP, "cva (vol_cv): cv%csa > 0 violated")
    
    vol_cv = cv%x * cv%csa
    
    call assert(vol_cv%v%v > 0.0_WP, "cva (vol_cv): vol_cv > 0 violated")
end function vol_cv

pure function rho_cv(cv)
    use units, only: si_mass_density => unit_m30_p10_p00_p00_p00, &
                     si_volume       => unit_p30_p00_p00_p00_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_mass_density) :: rho_cv
    
    call assert_mass(cv, "rho_cv")
    
    rho_cv = cv%m_total() / cv%vol()
    
    call assert(rho_cv%v%v > 0.0_WP, "cva (rho_cv): rho_cv > 0 violated")
end function rho_cv

pure function p_cv(cv)
    use units, only: si_mass_density => unit_m30_p10_p00_p00_p00, &
                     si_temperature  => unit_p00_p00_p00_p10_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_pressure) :: p_cv
    
    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    
    rho  = cv%rho()
    temp = cv%temp()
    
    p_cv = cv%p_eos(rho, temp)
    
    call assert(p_cv%v%v > 0.0_WP, "cva (p_cv): p_cv > 0 violated")
end function p_cv

pure subroutine set(cv, x, x_dot, y, p, temp, csa, rm_p, p_fs, p_fd, p_atm, k, x_z, gas)
    use units, only: si_temperature => unit_p00_p00_p00_p10_p00
    use checks, only: assert, assert_dimension, is_close
    
    class(cv_type), intent(in out) :: cv
    
    ! time varying
    type(si_length), intent(in)      :: x     ! location of piston/projectile
    type(si_velocity), intent(in)    :: x_dot ! velocity of piston/projectile
    type(unitless), intent(in)       :: y(:)  ! mass fractions of each gas
    type(si_pressure), intent(in)    :: p     ! pressure
    type(si_temperature), intent(in) :: temp  ! temperature
    
    ! constant
    type(si_area), intent(in)         :: csa        ! cross-sectional area
    type(si_inverse_mass), intent(in) :: rm_p       ! reciprocal mass of piston/projectile
    type(si_pressure), intent(in)     :: p_fs, p_fd ! static and dynamic friction pressure
    type(si_pressure), intent(in)     :: p_atm      ! atmospheric pressure
    type(si_stiffness), intent(in)    :: k          ! stiffness of spring attached to piston
    type(si_length), intent(in)       :: x_z        ! zero force location for spring
    type(gas_type), intent(in)        :: gas(:)     ! gas data
    
    integer        :: i, n_d
    type(si_mass)  :: m_total
    type(unitless) :: y_sum
    
    n_d = size(x%v%d)
    
    cv%x     = x
    cv%x_dot = x_dot
    ! `p` and `temp` will be handled below
    
    cv%csa   = csa
    cv%rm_p  = rm_p
    cv%p_fs  = p_fs
    cv%p_fd  = p_fd
    cv%p_atm = p_atm
    cv%k     = k
    cv%x_z   = x_z
    cv%gas   = gas
    
    call assert(cv%x%v%v     >  0.0_WP, "cva (set): x > 0 violated")
    call assert(p%v%v        >  0.0_WP, "cva (set): p > 0 violated")
    call assert(temp%v%v     >  0.0_WP, "cva (set): temp > 0 violated")
    call assert(csa%v%v      >  0.0_WP, "cva (set): csa > 0 violated")
    call assert(cv%p_fs%v%v  >= 0.0_WP, "cva (set): p_fs >= 0 violated")
    call assert(cv%p_fd%v%v  >= 0.0_WP, "cva (set): p_fd >= 0 violated")
    call assert(cv%p_atm%v%v >  0.0_WP, "cva (set): p_atm > 0 violated")
    call assert(cv%k%v%v     >= 0.0_WP, "cva (set): k >= 0 violated")
    call assert(cv%x_z%v%v   >= 0.0_WP, "cva (set): x_z >= 0 violated")
    
    call assert_dimension(y, cv%gas)
    allocate(cv%m(size(y)))
    call y_sum%v%init_const(0.0_WP, n_d)
    call m_total%v%init_const(1.0_WP, n_d)
    do i = 1, size(y)
        call assert(gas(i)%gamma > 1.0_WP, "cva (set): gas%gamma > 1 violated")
        call assert(gas(i)%mm    > 0.0_WP, "cva (set): gas%mm > 0 violated")
        call assert(gas(i)%mm    < 0.1_WP, "cva (set): gas%mm < 0.1 violated") ! to catch using g/mol by mistake
        call assert(gas(i)%p_c   > 0.0_WP, "cva (set): gas%p_c > 0 violated") ! to catch using g/mol by mistake
        
        y_sum = y_sum + y(i)
    end do
    
    call assert(is_close(y_sum%v%v, 1.0_WP), "cva (set): mass fractions do not sum to 1")
    
    m_total = cv%vol() * cv%rho_eos(p, temp, y)
    
    ! Now set `m` to the correct values and calculate `e`
    call cv%e%v%init_const(0.0_WP, n_d)
    do i = 1, size(y)
        cv%m(i) = y(i)*m_total
        cv%e    = cv%e + cv%m(i)*cv%gas(i)%u(temp)
    end do
    
    call assert_mass(cv, "set")
    call assert_dimension(cv%x%v%d, cv%x_dot%v%d)
    call assert_dimension(cv%x%v%d, cv%e%v%d)
    
    ! This is checked here and not in `rho_eos` as the masses are not defined when `rho_eos` is called.
    call assert(p < cv%p_c(), "cva (set): ideal gas law validity is questionable")
end subroutine set

pure function p_f(cv, p_fe)
    ! Returns pressure of friction.
    
    use units, only: tanh
    use checks, only: assert
    
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
    
    use checks, only: assert
    
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
        use units, only: tanh, abs, atanh
        use checks, only: is_close
        
        type(si_pressure), intent(in) :: p_fe, p_fs, p_s
        type(si_pressure) :: p_f0_high
        
        if (is_close(p_fs%v%v, 0.0_WP)) then
            p_f0_high = p_fs
        else
            p_f0_high = p_fs * tanh((abs(p_fe) - p_s)/(p_fs - p_s) + atanh(p_s/p_fs))
        end if
    end function p_f0_high
end function p_f0

pure function d_x_d_t(cv)
    class(cv_type), intent(in) :: cv
    
    type(si_velocity) :: d_x_d_t
    
    d_x_d_t = cv%x_dot
end function d_x_d_t

pure function d_xdot_d_t(cv)
    use units, only: si_acceleration => unit_p10_p00_m20_p00_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_acceleration) :: d_xdot_d_t
    
    type(si_pressure) :: p_fe ! friction pressure at equilibrium ($\partial \dot{x}/\partial t = 0$)
    
    call assert(cv%csa%v%v > 0.0_WP, "cva (d_xdot_d_t): cv%csa > 0 violated")
    
    p_fe = cv%p() - cv%p_atm - (cv%k/cv%csa)*(cv%x - cv%x_z)
    
    d_xdot_d_t = cv%csa*cv%rm_p*(cv%p() - cv%p_atm - cv%p_f(p_fe)) - cv%k*cv%rm_p*(cv%x - cv%x_z)
end function d_xdot_d_t

pure subroutine assert_mass(cv, procedure_name)
    ! Why not make this a type-bound operator?
    ! That would make my assertion counting Python program not count these.
    
    use units, only: si_temperature => unit_p00_p00_p00_p10_p00
    use checks, only: assert, assert_dimension
    
    type(cv_type), intent(in)    :: cv
    character(len=*), intent(in) :: procedure_name
    
    integer       :: i
    type(si_mass) :: m_total
    
    call assert(len(trim(procedure_name)) > 0, "cva (assert_mass): procedure name should not be empty")
    
    do i = 1, size(cv%m)
        call assert(cv%m(i)%v%v >= 0.0_WP, "cva (" // trim(procedure_name) // "): cv%m >= 0 violated")
        call assert_dimension(cv%e%v%d, cv%m(i)%v%d)
    end do
    
    m_total = cv%m_total()
    call assert(m_total%v%v > 0.0_WP, "cva (" // trim(procedure_name) // "): cv%m_total > 0 violated")
    
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
    
    use checks, only: assert, assert_dimension
    use units, only: log, exp
    
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
    
    use units, only: tanh, square
    use checks, only: assert, assert_dimension
    
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
    
    use units, only: tanh
    use checks, only: assert
    
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
    
    use units, only: si_mass_flow_rate => unit_p00_p10_m10_p00_p00, &
                     si_specific_heat  => unit_p20_p00_m20_m10_p00, &
                     sqrt
    use checks, only: assert, assert_dimension
    
    class(con_type), intent(in) :: con
    type(cv_type), intent(in)   :: cv_from, cv_to
    
    type(si_mass_flow_rate) :: m_dot
    
    integer        :: n_d
    type(unitless) :: p_r
    
    call assert_dimension(con%a_e%v%d, con%b%v%d)
    
    call assert(con%active, "cva (m_dot): connection must be active")
    call assert(cv_from%p() >= cv_to%p(), "cva (m_dot): cv_from%p >= cv_to%p violated")
    
    n_d = size(con%a_e%v%d)
    
    p_r = cv_to%p() / cv_from%p()
    call assert(p_r%v%v >= 0.0_WP, "cva (m_dot): p_r >= 0 violated")
    call assert(p_r%v%v <= 1.0_WP, "cva (m_dot): p_r <= 1 violated")
    
    m_dot = con%a_e * (cv_from%p() - g_m_dot(p_r) * cv_to%p()) &
                * sqrt((1.0_WP - con%b) / (cv_from%r() * cv_from%temp())) &
                * sqrt(1.0_WP - f_m_dot(p_r, con%b))
end function m_dot

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions for air humidity !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function p_v_h2o(temp)
    ! Using the Antoine equation with NIST's parameters for the experiments of Stull, 1947.
    ! Data ranges from 255.9 to 373.0 K. Probably okay to go a bit outside of those ranges.
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2>
    
    use units, only: si_temperature => unit_p00_p00_p00_p10_p00, &
                     exp
    
    type(si_temperature), intent(in) :: temp
    
    type(si_pressure) :: p_v_h2o
    
    integer              :: n_d
    type(si_pressure)    :: a
    type(si_temperature) :: b, c
    
    n_d = size(temp%v%d)
    
    call a%v%init_const((10.0_WP**4.6543_WP)*1.0e5_WP, n_d)
    call b%v%init_const(1435.264_WP*log(10.0_WP), n_d)
    call c%v%init_const(-64.848_WP, n_d)
    
    p_v_h2o = a*exp(-b/(temp + c))
end function p_v_h2o

end module cva
