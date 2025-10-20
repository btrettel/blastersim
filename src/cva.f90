! Module for control volume analysis.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module cva

use prec, only: WP
use units, only: si_length       => unit_p10_p00_p00_p00, &
                 si_velocity     => unit_p10_p00_m10_p00, &
                 si_mass         => unit_p00_p10_p00_p00, &
                 si_energy       => unit_p20_p10_m20_p00, &
                 si_area         => unit_p20_p00_p00_p00, &
                 si_inverse_mass => unit_p00_m10_p00_p00, &
                 si_pressure     => unit_m10_p10_m20_p00, &
                 unitless        => unit_p00_p00_p00_p00, &
                 si_stiffness    => unit_p00_p10_m20_p00, &
                 si_area         => unit_p20_p00_p00_p00
implicit none
private

public :: p_eos, rho_eos
public :: smooth_min
public :: f_m_dot, g_m_dot, m_dot

! <https://en.wikipedia.org/wiki/Gas_constant>
real(WP), public, parameter :: R_BAR = 8.31446261815324_WP ! J/(mol*K)

! <https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html>
! Molecular mass of air
real(WP), public, parameter :: M_AIR = 28.9647e-3_WP ! kg/mol

! moran_fundamentals_2008 table A-1
! Critical pressure of air
real(WP), public, parameter :: P_C_AIR = 37.7e5_WP ! Pa

! moran_fundamentals_2008 table A-20
! Specific heat ratio of air at 300 K
real(WP), public, parameter :: K_AIR = 1.400_WP ! unitless

! <https://en.wikipedia.org/wiki/Density_of_air>
real(WP), public, parameter :: P_ATM   = 101325.0_WP         ! Pa
real(WP), public, parameter :: T_ATM   = 273.15_WP + 15.0_WP ! K
real(WP), public, parameter :: RHO_ATM = 1.2250_WP           ! kg/m3

! pressure ratio laminar flow nominally starts at
! based on first part of beater_pneumatic_2007 eq. 5.4
real(WP), public, parameter :: P_RL = 0.999_WP ! unitless

type, public :: cv_type
    ! time varying
    type(si_length)   :: x     ! location of piston/projectile
    type(si_velocity) :: x_dot ! velocity of piston/projectile
    type(si_mass)     :: m     ! mass of gas in control volume
    type(si_energy)   :: e     ! energy of gas in control volume
    
    ! constants
    type(si_area)         :: csa        ! cross-sectional area
    type(si_inverse_mass) :: rm_p       ! reciprocal mass of piston/projectile
    type(si_pressure)     :: p_fs, p_fd ! static and dynamic friction pressure
    type(si_stiffness)    :: k          ! stiffness of spring attached to piston
    type(si_length)       :: x_z        ! zero force location for spring
contains
    procedure :: p_f  => p_f
    procedure :: p_f0 => p_f0
    procedure :: temp => temp_cv
    procedure :: vol  => vol_cv
    procedure :: rho  => rho_cv
    procedure :: p    => p_cv
    procedure :: set  => set
end type cv_type

type, public :: valve_type
    type(si_area)  :: a_e ! effective area
    type(unitless) :: b   ! critical pressure ratio
contains
    procedure :: m_dot => m_dot
end type valve_type

contains

pure function p_eos(rho, temp)
    ! Calculate pressure using the equation of state.
    
    use units, only: si_mass_density  => unit_m30_p10_p00_p00, &
                     si_temperature   => unit_p00_p00_p00_p10, &
                     si_specific_heat => unit_p20_p00_m20_m10
    use checks, only: assert, assert_dimension
    
    type(si_mass_density), intent(in) :: rho
    type(si_temperature), intent(in)  :: temp
    type(si_pressure)                 :: p_eos
    
    integer :: n_d ! number of derivatives
    
    type(si_specific_heat) :: r_air
    
    call assert(rho%v%v  > 0.0_WP, "cva (p_eos): rho%v > 0 violated")
    call assert(temp%v%v > 0.0_WP, "cva (p_eos): temp%v > 0 violated")
    call assert_dimension(rho%v%d, temp%v%d)
    
    n_d = size(rho%v%d)
    
    call r_air%v%init_const(R_BAR/M_AIR, n_d)
    
    p_eos = rho * r_air * temp
    
    call assert(p_eos%v%v > 0.0_WP, "cva (p_eos): p_eos%v > 0 violated")
    call assert(p_eos%v%v < P_C_AIR, "cva (p_eos): ideal gas law validity is questionable")
end function p_eos

pure function rho_eos(p, temp)
    ! Calculate density using the equation of state.
    ! In the future, if an EOS more complex than the ideal gas law is used, it might make sense to calculate `rho_eos` from `p_eos`
    ! with the Newton method.
    
    use units, only: si_mass_density  => unit_m30_p10_p00_p00, &
                     si_temperature   => unit_p00_p00_p00_p10, &
                     si_specific_heat => unit_p20_p00_m20_m10
    use checks, only: assert, assert_dimension
    
    type(si_pressure), intent(in)    :: p
    type(si_temperature), intent(in) :: temp
    
    type(si_mass_density) :: rho_eos
    
    integer :: n_d ! number of derivatives
    
    type(si_specific_heat) :: r_air
    
    call assert(p%v%v    > 0.0_WP, "cva (rho_eos): p%v > 0 violated")
    call assert(p%v%v    < P_C_AIR, "cva (rho_eos): ideal gas law validity is questionable")
    call assert(temp%v%v > 0.0_WP, "cva (rho_eos): temp%v > 0 violated")
    call assert_dimension(p%v%d, temp%v%d)
    
    n_d = size(p%v%d)
    
    call r_air%v%init_const(R_BAR/M_AIR, n_d)
    
    rho_eos = p / (r_air * temp)
    
    call assert(rho_eos%v%v > 0.0_WP, "cva (rho_eos): rho_eos%v > 0 violated")
end function rho_eos

pure function temp_cv(cv)
    use units, only: si_temperature     => unit_p00_p00_p00_p10, &
                     si_specific_energy => unit_p20_p00_m20_p00, &
                     si_specific_heat   => unit_p20_p00_m20_m10
    use checks, only: assert, assert_dimension
    
    class(cv_type), intent(in) :: cv
    
    type(si_temperature) :: temp_cv
    
    integer                  :: n_d  ! number of derivatives
    type(si_specific_heat)   :: r_air, c_v
    type(si_specific_energy) :: u
    
    ! Constant specific heats assumed for now. Will improve later.
    
    call assert(K_AIR > 1.0_WP, "cva (temp_cv): K_AIR > 1 violated")
    call assert(cv%m%v%v > 0.0_WP, "cva (temp_cv): cv%m > 0 violated")
    call assert(cv%e%v%v > 0.0_WP, "cva (temp_cv): cv%e > 0 violated")
    call assert_dimension(cv%e%v%d, cv%m%v%d)
    
    n_d = size(cv%m%v%d)
    
    call r_air%v%init_const(R_BAR/M_AIR, n_d)
    
    ! moran_fundamentals_2008 eq. 3.47b, p. 119
    c_v = r_air / (K_AIR - 1.0_WP)
    
    u = cv%e / cv%m
    
    temp_cv = u / c_v
    
    call assert(temp_cv%v%v > 0.0_WP, "cva (temp_cv): temp_cv > 0 violated")
end function temp_cv

pure function vol_cv(cv)
    use units, only: si_volume => unit_p30_p00_p00_p00, &
                     si_length => unit_p10_p00_p00_p00, &
                     si_area   => unit_p20_p00_p00_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_volume) :: vol_cv
    
    call assert(cv%x%v%v > 0.0_WP, "cva (vol_cv): cv%x > 0 violated")
    call assert(cv%csa%v%v > 0.0_WP, "cva (vol_cv): cv%csa > 0 violated")
    
    vol_cv = cv%x * cv%csa
    
    call assert(vol_cv%v%v > 0.0_WP, "cva (vol_cv): vol_cv > 0 violated")
end function vol_cv

pure function rho_cv(cv)
    use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                     si_volume       => unit_p30_p00_p00_p00
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_mass_density) :: rho_cv
    type(si_volume)       :: vol
    
    call assert(cv%m%v%v > 0.0_WP, "cva (rho_cv): cv%m > 0 violated")
    
    vol    = cv%vol()
    rho_cv = cv%m / vol
    
    call assert(rho_cv%v%v > 0.0_WP, "cva (rho_cv): rho_cv > 0 violated")
end function rho_cv

pure function p_cv(cv)
    use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                     si_temperature  => unit_p00_p00_p00_p10
    use checks, only: assert
    
    class(cv_type), intent(in) :: cv
    
    type(si_pressure) :: p_cv
    
    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    
    rho  = cv%rho()
    temp = cv%temp()
    
    p_cv = p_eos(rho, temp)
    
    call assert(p_cv%v%v > 0.0_WP, "cva (p_cv): p_cv > 0 violated")
end function p_cv

pure subroutine set(cv, x, x_dot, p, temp, csa, m_p, p_fs, p_fd, k, x_z)
    use units, only: si_temperature     => unit_p00_p00_p00_p10, &
                     si_mass            => unit_p00_p10_p00_p00, &
                     si_specific_energy => unit_p20_p00_m20_p00, &
                     si_specific_heat   => unit_p20_p00_m20_m10
    use checks, only: assert, assert_dimension
    
    class(cv_type), intent(in out) :: cv
    
    ! time varying
    type(si_length), intent(in)      :: x     ! location of piston/projectile
    type(si_velocity), intent(in)    :: x_dot ! velocity of piston/projectile
    type(si_pressure), intent(in)    :: p     ! pressure
    type(si_temperature), intent(in) :: temp  ! temperature
    
    ! constant
    type(si_area), intent(in)      :: csa        ! cross-sectional area
    type(si_mass), intent(in)      :: m_p        ! mass of piston/projectile
    type(si_pressure), intent(in)  :: p_fs, p_fd ! static and dynamic friction pressure
    type(si_stiffness), intent(in) :: k          ! stiffness of spring attached to piston
    type(si_length), intent(in)    :: x_z        ! zero force location for spring
    
    integer                  :: n_d
    type(si_specific_heat)   :: r_air, c_v
    type(si_specific_energy) :: u
    
    ! Constant specific heats assumed for now. Will improve later.
    
    call assert(K_AIR > 1.0_WP, "cva (temp_cv): K_AIR > 1 violated")
    
    n_d = size(x%v%d)
    
    cv%x     = x
    cv%x_dot = x_dot
    ! `p` and `temp` will be handled below
    
    cv%csa  = csa
    cv%rm_p = 1.0_WP/m_p ! reciprocal mass of piston/projectile
    cv%p_fs = p_fs
    cv%p_fd = p_fd
    cv%k    = k
    cv%x_z  = x_z
    
    call assert(cv%x%v%v    >  0.0_WP, "cva (set): x > 0 violated")
    call assert(p%v%v       >  0.0_WP, "cva (set): p > 0 violated")
    call assert(temp%v%v    >  0.0_WP, "cva (set): temp > 0 violated")
    call assert(csa%v%v     >  0.0_WP, "cva (set): csa > 0 violated")
    call assert(cv%p_fs%v%v >= 0.0_WP, "cva (set): p_fs > 0 violated")
    call assert(cv%p_fd%v%v >= 0.0_WP, "cva (set): p_fd > 0 violated")
    call assert(cv%k%v%v    >= 0.0_WP, "cva (set): k >= 0 violated")
    call assert(cv%x_z%v%v  >= 0.0_WP, "cva (set): x_z >= 0 violated")
    
    cv%m = cv%vol() * rho_eos(p, temp)
    
    call r_air%v%init_const(R_BAR/M_AIR, n_d)
    
    ! moran_fundamentals_2008 eq. 3.47b, p. 119
    c_v = r_air / (K_AIR - 1.0_WP)
    
    u = c_v * temp
    
    cv%e = u * cv%m
    
    call assert(cv%m%v%v > 0.0_WP, "cva (set): cv%m > 0 violated")
    call assert(cv%e%v%v > 0.0_WP, "cva (set): cv%e > 0 violated")
    
    call assert_dimension(cv%x%v%d, cv%x_dot%v%d)
    call assert_dimension(cv%x%v%d, cv%m%v%d)
    call assert_dimension(cv%x%v%d, cv%e%v%d)
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
    
    call assert(smooth_min <= x, "cva (smooth_min): smooth_min <= x violated")
    call assert(smooth_min <= y, "cva (smooth_min): smooth_min <= y violated")
end function smooth_min

pure function f_m_dot(p_r, b)
    ! See beater_pneumatic_2007 eq. 5.4
    ! This is a replacement for the ((p_2/p_1 - b)/(1-b))**2 term, smoothly going between the various cases.
    
    use units, only: tanh, square
    use checks, only: assert_dimension
    
    type(unitless), intent(in) :: p_r, b
    
    type(unitless) :: f_m_dot
    
    type(unitless) :: p_rs, p_rl_ ! scales used to make function differentiable
    
    call assert_dimension(p_r%v%d, b%v%d)
    
    call p_rs%v%init_const(0.01_WP, size(p_r%v%d)) ! TODO: make `p_rs` a function of `dt`
    call p_rl_%v%init_const(P_RL, size(p_r%v%d)) ! based on first part of beater_pneumatic_2007 eq. 5.4
    
    f_m_dot = square((smooth_min(p_r, p_rl_) - b) / (1.0_WP - b)) &
                * 0.5_WP * (1.0_WP + tanh((p_r - b) / p_rs))
end function f_m_dot

pure function g_m_dot(p_r)
    ! See beater_pneumatic_2007 eq. 5.4
    ! This is a replacement for the p_1 term, smoothly going between the various cases.
    
    use units, only: tanh
    use checks, only: assert
    
    type(unitless), intent(in) :: p_r
    
    type(unitless) :: g_m_dot
    
    type(unitless) :: p_rs, p_rl_ ! scales used to make function differentiable
    
    call p_rs%v%init_const((1.0_WP - P_RL) / 10.0_WP, size(p_r%v%d)) ! TODO: make `p_rs` a function of `dt`
    call p_rl_%v%init_const(P_RL, size(p_r%v%d)) ! based on first part of beater_pneumatic_2007 eq. 5.4
    
    g_m_dot = 0.5_WP * (1.0_WP + tanh((p_r - p_rl_) / p_rs))
    
    call assert(g_m_dot%v%v >= 0.0_WP, "cva (g_m_dot): g_m_dot >= 0 violated")
    call assert(g_m_dot%v%v <= 1.0_WP, "cva (g_m_dot): g_m_dot <= 1 violated")
end function g_m_dot

pure function m_dot(valve, cv_from, cv_to)
    ! Modified valve flow rate model from beater_pneumatic_2007 ch. 5.
    ! Modified to be differentiable.
    
    use units, only: si_mass_flow_rate => unit_p00_p10_m10_p00, &
                     si_specific_heat  => unit_p20_p00_m20_m10, &
                     sqrt
    use checks, only: assert, assert_dimension
    
    class(valve_type), intent(in) :: valve
    type(cv_type), intent(in)     :: cv_from, cv_to
    
    type(si_mass_flow_rate) :: m_dot
    
    integer                :: n_d
    type(si_specific_heat) :: r_air
    type(unitless)         :: p_r
    
    call assert_dimension(valve%a_e%v%d, valve%b%v%d)
    
    call assert(cv_from%p() >= cv_to%p(), "cva (m_dot): cv_from%p >= cv_to%p violated")
    
    n_d = size(valve%a_e%v%d)
    call r_air%v%init_const(R_BAR/M_AIR, n_d)
    
    p_r = cv_to%p() / cv_from%p()
    call assert(p_r%v%v >  0.0_WP, "cva (m_dot): p_r > 0 violated")
    call assert(p_r%v%v <= 0.0_WP, "cva (m_dot): p_r <= 1 violated")
    
    m_dot = valve%a_e * (cv_from%p() - g_m_dot(p_r) * cv_to%p()) &
                * sqrt((1.0_WP - valve%b) / (r_air * cv_from%temp())) &
                * sqrt(1.0_WP - f_m_dot(p_r, valve%b))
end function m_dot

end module cva
