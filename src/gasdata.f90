! Module for gas thermochemical data.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module gasdata

use prec, only: WP
use convert, only: CONVERT_C_TO_K, CONVERT_PA_TO_KPA
use units
use checks, only: assert, assert_dimension, is_close
implicit none
private

public :: p_v_h2o

! <https://en.wikipedia.org/wiki/Gas_constant>
real(WP), public, parameter :: R_BAR = 8.31446261815324_WP ! J/(mol*K)

! <https://en.wikipedia.org/wiki/Density_of_air>
real(WP), public, parameter :: P_ATM     = 101325.0_WP              ! Pa
real(WP), public, parameter :: P_ATM_KPA = P_ATM*CONVERT_PA_TO_KPA  ! kPa
real(WP), public, parameter :: TEMP_ATM  = CONVERT_C_TO_K + 15.0_WP ! K
real(WP), public, parameter :: RHO_ATM   = 1.2250_WP                ! kg/m3

real(WP), public, parameter :: TEMP_0 = 300.0_WP ! K, temperature that `gamma`, `u_0`, and `h_0` are taken at in `gas_type`

type, public :: gas_type
    character(len=32) :: label ! human-readable label for gas
    real(WP)          :: gamma ! ratio of specific heats, unitless
    real(WP)          :: u_0   ! internal energy at `TEMP_0`, J/kg
    real(WP)          :: h_0   ! enthalphy at `TEMP_0`, J/kg
    real(WP)          :: mm    ! molar mass, kg/mol
    real(WP)          :: p_c   ! critical pressure, Pa
contains
    procedure :: u   => u_gas
    procedure :: h   => h_gas
    procedure :: r   => r_gas
    procedure :: c_v => c_v_gas
    procedure :: c_p => c_p_gas
end type gas_type

! Molecular mass and critical pressure of air (consistent with dry air): moran_fundamentals_2008 table A-1
! Specific heat ratio of air (consistent with dry air): moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of air: moran_fundamentals_2008 table A-22
type(gas_type), public, parameter :: DRY_AIR = gas_type(label = "dry air", &
                                                        gamma = 1.400_WP, &
                                                        u_0   = 214.07e3_WP, &
                                                        h_0   = 300.19e3_WP, &
                                                        mm    = 28.97e-3_WP, &
                                                        p_c   = 37.7e5_WP)

! Molecular mass and critical pressure of N2: moran_fundamentals_2008 table A-1
! Specific heat ratio of N2: moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of N2: moran_fundamentals_2008 table A-23
type(gas_type), public, parameter :: N2      = gas_type(label = "N2", &
                                                        gamma = 1.400_WP, &
                                                        u_0   = 28.01e-3_WP*6229.0e3_WP, &
                                                        h_0   = 28.01e-3_WP*8723.0e3_WP, &
                                                        mm    = 28.01e-3_WP, &
                                                        p_c   = 33.9e5_WP)

! Molecular mass and critical pressure of O2: moran_fundamentals_2008 table A-1
! Specific heat ratio of O2: moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of O2: moran_fundamentals_2008 table A-23
type(gas_type), public, parameter :: O2      = gas_type(label = "O2", &
                                                        gamma = 1.395_WP, &
                                                        u_0   = 6242.0e3_WP, &
                                                        h_0   = 8736.0e3_WP, &
                                                        mm    = 32.0e-3_WP, &
                                                        p_c   = 50.5e5_WP)

! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/Ar> ("Fluid Properties")
type(gas_type), public, parameter :: AR      = gas_type(label = "AR", &
                                                        gamma = 0.52154_WP/0.31239_WP, &
                                                        u_0   = 155.90e3_WP, &
                                                        h_0   = 93.497e3_WP, &
                                                        mm    = 39.948e-3_WP, &
                                                        p_c   = 4.8630e6_WP)

! Molecular mass and critical pressure of CO2: moran_fundamentals_2008 table A-1
! Specific heat ratio of CO2: moran_fundamentals_2008 table A-20
! Internal energy and enthalpy of CO2: moran_fundamentals_2008 table A-23
type(gas_type), public, parameter :: CO2     = gas_type(label = "CO2", &
                                                        gamma = 1.288_WP, &
                                                        u_0   = 44.01e-3_WP*6939.0e3_WP, &
                                                        h_0   = 44.01e-3_WP*9431.0e3_WP, &
                                                        mm    = 44.01e-3_WP, &
                                                        p_c   = 73.9e5_WP)

! Molecular mass and critical pressure of H2O: moran_fundamentals_2008 table A-1
! Specific heat ratio of H2O: <https://en.wikipedia.org/wiki/Heat_capacity_ratio>
! Internal energy and enthalpy of H2O: moran_fundamentals_2008 table A-23
! TODO: Switch to stream tables for `u` and `h` if adding liquid water; see moran_thermodynamics_2008 pp. 666--667
type(gas_type), public, parameter :: H2O     = gas_type(label = "H2O", &
                                                        gamma = 1.330_WP, & ! at 20 C, not 300 K, but close enough
                                                        u_0   = 18.02e-3_WP*7472.0e3_WP, &
                                                        h_0   = 18.02e-3_WP*9966.0e3_WP, &
                                                        mm    = 18.02e-3_WP, &
                                                        p_c   = 220.9e5_WP)

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!
! methods for `gas_type` !
!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function u_gas(gas, temp)
    ! Constant specific heats assumed for now. Will improve later.
    
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
    
    class(gas_type), intent(in) :: gas
    integer, intent(in)         :: n_d ! number of derivatives
    
    type(si_specific_heat) :: r_gas
    
    call r_gas%v%init_const(R_BAR/gas%mm, n_d)
end function r_gas

pure function c_v_gas(gas, n_d)
    class(gas_type), intent(in) :: gas
    integer, intent(in)         :: n_d ! number of derivatives
    
    type(si_specific_heat) :: c_v_gas
    
    call assert(gas%gamma > 1.0_WP, "gasdata (temp_cv): gamma > 1 violated", print_real=[gas%gamma])
    
    ! moran_fundamentals_2008 eq. 3.47b, p. 119
    c_v_gas = gas%r(n_d) / (gas%gamma - 1.0_WP)
end function c_v_gas

pure function c_p_gas(gas, n_d)
    class(gas_type), intent(in) :: gas
    integer, intent(in)         :: n_d ! number of derivatives
    
    type(si_specific_heat) :: c_p_gas
    
    call assert(gas%gamma > 1.0_WP, "gasdata (temp_cv): gamma > 1 violated")
    
    ! moran_fundamentals_2008 eq. 3.47a, p. 119
    c_p_gas = gas%gamma * gas%r(n_d) / (gas%gamma - 1.0_WP)
end function c_p_gas

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! functions for air humidity !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

pure function p_v_h2o(temp)
    ! Using the Antoine equation with NIST's parameters for the experiments of Stull, 1947.
    ! Data ranges from 255.9 to 373.0 K. Probably okay to go a bit outside of those ranges.
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/H2O/h1H2>
    
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

end module gasdata
