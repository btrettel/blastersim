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
                 si_stiffness    => unit_p00_p10_m20_p00
implicit none
private

! <https://en.wikipedia.org/wiki/Gas_constant>
real(WP), public, parameter :: R_BAR = 8.31446261815324_WP ! J/(mol*K)

! <https://www.engineeringtoolbox.com/molecular-mass-air-d_679.html>
! Molecular mass of air
real(WP), public, parameter :: M_AIR = 28.9647e-3_WP ! kg/mol

! moran_fundamentals_2008 table A-1
! Critical pressure of air
real(WP), public, parameter :: P_C_AIR = 37.7e5_WP ! Pa

! <https://en.wikipedia.org/wiki/Density_of_air>
real(WP), public, parameter :: P_ATM   = 101325.0_WP         ! Pa
real(WP), public, parameter :: T_ATM   = 273.15_WP + 15.0_WP ! K
real(WP), public, parameter :: RHO_ATM = 1.2250_WP           ! kg/m3

public :: p_eos

type, public :: cv_type
    ! time varying
    type(si_length)   :: x     ! location of piston/projectile
    type(si_velocity) :: x_dot ! velocity of piston/projectile
    type(si_mass)     :: m     ! mass of gas in control volume
    type(si_energy)   :: e     ! energy of gas in control volume
    
    ! constants
    type(si_area)         :: csa        ! cross-sectional area
    type(si_inverse_mass) :: rmp        ! reciprocol mass of piston/projectile
    type(si_pressure)     :: p_fs, p_fd ! static and dynamic friction pressure
    type(si_stiffness)    :: k          ! stiffness of spring attached to piston
    type(si_length)       :: x_z        ! zero force location for spring
contains
    !procedure :: p => p_cv
    procedure :: p_f  => p_f
    procedure :: p_f0 => p_f0
end type cv_type

contains

pure function p_eos(rho, temp)
    use fmad, only: ad
    use units, only: si_mass_density  => unit_m30_p10_p00_p00, &
                     si_temperature   => unit_p00_p00_p00_p10, &
                     si_specific_heat => unit_p20_p00_m20_m10
    use checks, only: assert
    
    type(si_mass_density), intent(in) :: rho
    type(si_temperature), intent(in)  :: temp
    type(si_pressure)                 :: p_eos
    
    integer :: n_d
    
    type(si_specific_heat) :: R
    
    call assert(rho%v%v > 0.0_WP, "cva (p_eos): rho%v > 0 violated")
    call assert(temp%v%v > 0.0_WP, "cva (p_eos): temp%v > 0 violated")
    
    n_d = size(rho%v%d)
    
    call R%v%init_const(R_BAR/M_AIR, n_d)
    
    p_eos = rho * R * temp
    
    call assert(p_eos%v%v > 0.0_WP, "cva (p_eos): p_eos%v > 0 violated")
    call assert(p_eos%v%v < P_C_AIR, "cva (p_eos): ideal gas law validity is questionable")
end function p_eos

!pure function p_cv(cv)
!    use units, only: si_mass_density  => unit_m30_p10_p00_p00
    
!    class(cv_type), intent(in) :: cv
    
    
!end function p_cv(cv)

pure function p_f(cv, p_fe)
    ! Returns pressure of friction.
    
    use units, only: tanh
    use checks, only: assert
    
    class(cv_type), intent(in)    :: cv
    type(si_pressure), intent(in) :: p_fe ! equilibrium pressure
    
    type(si_pressure) :: p_f
    
    type(si_velocity) :: v_scale
    
    call v_scale%v%init_const(0.1_WP, size(cv%x%v%d))
    
    call assert(cv%p_fs%v%v > 0.0_WP, "cva (p_f): cv%p_fs%v > 0 violated")
    call assert(cv%p_fd%v%v > 0.0_WP, "cva (p_f): cv%p_fd%v > 0 violated")
    
    p_f = p_f0(cv, p_fe) + (cv%p_fd - tanh(cv%x_dot/v_scale)*p_f0(cv, p_fe))*tanh(cv%x_dot/v_scale)
end function p_f

pure function p_f0(cv, p_fe)
    ! Returns actual static pressure of friction.
    ! `p_fs` is the *maximum* static pressure of friction.
    
    use units, only: tanh
    use checks, only: assert
    
    class(cv_type), intent(in)    :: cv
    type(si_pressure), intent(in) :: p_fe ! equilibrium pressure
    
    type(si_pressure) :: p_f0
    
    call assert(cv%p_fs%v%v > 0.0_WP, "cva (p_f0): cv%p_fs%v > 0 violated")
    call assert(cv%p_fd%v%v > 0.0_WP, "cva (p_f0): cv%p_fd%v > 0 violated")
    
    p_f0 = cv%p_fs * tanh(p_fe/cv%p_fs)
end function p_f0

end module cva
