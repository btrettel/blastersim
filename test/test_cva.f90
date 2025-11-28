! tests for control volume analysis
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_cva

use prec, only: WP
use units
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("cva.nml")

call test_m_total(tests)
call test_p_eos(tests)
call test_rho_eos(tests)
call test_r_cv(tests)
call test_p_c(tests)
call test_p_f_1(tests)
call test_p_f_2(tests)
call test_p_f0_1(tests)
call test_p_f0_2(tests)
call test_temp_cv(tests)
call test_set_1(tests)
call test_set_2(tests)
call test_set_3(tests)
call test_rates(tests)
call test_u_h_cv(tests)
call test_gamma_cv(tests)

call test_smooth_min(tests)
call test_f_m_dot(tests)
call test_g_m_dot(tests)
call test_m_dot_1(tests)
call test_m_dot_2(tests)
call test_m_dot_3(tests)
call test_m_dot_4(tests)

call test_calculate_flows(tests)
call test_conservation(tests)

call tests%end_tests()

contains

subroutine test_m_total(tests)
    use gasdata, only: DRY_AIR
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests
    
    type(cv_type) :: cv
    
    type(si_mass)  :: m_total
    type(unitless) :: y(2), chi(2)
    
    allocate(cv%m(2))
    call cv%m(1)%v%init_const(1.0_WP, 0)
    call cv%m(2)%v%init_const(5.0_WP, 0)
    
    allocate(cv%gas(2))
    cv%gas(1) = DRY_AIR
    cv%gas(2) = DRY_AIR
    
    m_total = cv%m_total()
    
    call tests%real_eq(m_total%v%v, 6.0_WP, "m_total")
    
    y = cv%y()
    call tests%real_eq(y(1)%v%v, 1.0_WP/6.0_WP, "y(1)")
    call tests%real_eq(y(2)%v%v, 5.0_WP/6.0_WP, "y(2)")
    
    chi = cv%chi()
    call tests%real_eq(chi(1)%v%v, 1.0_WP/6.0_WP, "chi(1)")
    call tests%real_eq(chi(2)%v%v, 5.0_WP/6.0_WP, "chi(2)")
end subroutine test_m_total

subroutine test_p_eos(tests)
    use gasdata, only: P_ATM, TEMP_ATM, RHO_ATM, DRY_AIR
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    type(cv_type)         :: cv
    
    call rho%v%init_const(RHO_ATM, 0)
    call temp%v%init_const(TEMP_ATM, 0)
    allocate(cv%m(1))
    call cv%m(1)%v%init_const(1.0_WP, 0)
    call cv%e%v%init_const(1.0_WP, 0)
    cv%i_cv_other = 2
    
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    p = cv%p_eos(rho, temp)
    
    call tests%real_eq(p%v%v, P_ATM, "p_eos, atmospheric", abs_tol=20.0_WP)
end subroutine test_p_eos

subroutine test_rho_eos(tests)
    use gasdata, only: P_ATM, TEMP_ATM, RHO_ATM, DRY_AIR
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    type(unitless)        :: y(1)
    type(cv_type)         :: cv
    
    call p%v%init_const(P_ATM, 0)
    call temp%v%init_const(TEMP_ATM, 0)
    call y%v%init_const(1.0_WP, 0)
    allocate(cv%m(1))
    call cv%m(1)%v%init_const(1.0_WP, 0)
    call cv%e%v%init_const(1.0_WP, 0)
    
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    rho = cv%rho_eos(p, temp, y)
    call tests%real_eq(rho%v%v, RHO_ATM, "rho_eos, atmospheric (1)", abs_tol=1.0e-3_WP)
end subroutine test_rho_eos

subroutine test_r_cv(tests)
    ! Test using ambient air.
    ! Tests both the `r_cv` and `rho_eos`
    
    use gasdata, only: P_ATM, TEMP_ATM, RHO_ATM, R_BAR, gas_type, DRY_AIR, N2, O2, AR, CO2
    use cva, only: cv_type
    use checks, only: assert, is_close
    
    type(test_results_type), intent(in out) :: tests
    
    real(WP)               :: chi(4)
    type(unitless)         :: y(4)
    type(cv_type)          :: cv
    type(si_specific_heat) :: r_cv
    type(si_mass_density)  :: rho
    type(si_temperature)   :: temp
    type(si_pressure)      :: p
    type(unitless)         :: gamma_cv
    
    ! Dry air
    ! <https://en.wikipedia.org/wiki/Atmosphere_of_Earth>
    chi(1) = 0.7808_WP ! N2
    chi(2) = 0.2095_WP ! O2
    chi(3) = 0.0093_WP ! Ar
    chi(4) = 0.0004_WP ! CO2
    
    call assert(is_close(sum(chi), 1.0_WP), "test_cva (test_r_cv): chi array doesn't sum to 1")
    
    cv%gas = [N2, O2, AR, CO2]
    
    allocate(cv%m(4))
    call cv%m(1)%v%init_const(chi(1)*cv%gas(1)%mm, 0)
    call cv%m(2)%v%init_const(chi(2)*cv%gas(2)%mm, 0)
    call cv%m(3)%v%init_const(chi(3)*cv%gas(3)%mm, 0)
    call cv%m(4)%v%init_const(chi(4)*cv%gas(4)%mm, 0)
    
    y(1) = cv%m(1) / (cv%m(1) + cv%m(2) + cv%m(3) + cv%m(4))
    y(2) = cv%m(2) / (cv%m(1) + cv%m(2) + cv%m(3) + cv%m(4))
    y(3) = cv%m(3) / (cv%m(1) + cv%m(2) + cv%m(3) + cv%m(4))
    y(4) = cv%m(4) / (cv%m(1) + cv%m(2) + cv%m(3) + cv%m(4))
    
    ! Needed to avoid an assertion failing.
    call cv%e%v%init_const(1.0_WP, 0)
    
    r_cv = cv%r()
    call tests%real_eq(r_cv%v%v, 0.2870e3_WP, "cv%r for mixture (moran_fundamentals_2008 table 3.1)", abs_tol=0.1_WP)
    call tests%real_eq(r_cv%v%v, R_BAR/DRY_AIR%mm, "cv%r for mixture (DRY_AIR)", abs_tol=0.1_WP)
    
    ! Test for `rho_eos` with multiple gas species.
    call p%v%init_const(P_ATM, 0)
    call temp%v%init_const(TEMP_ATM, 0)
    rho = cv%rho_eos(p, temp, y)
    call tests%real_eq(rho%v%v, RHO_ATM, "rho_eos, atmospheric (2)", abs_tol=1.0e-3_WP)
    
    gamma_cv = cv%gamma(y)
    call tests%real_eq(gamma_cv%v%v, DRY_AIR%gamma, "cv%gamma for mixture (1)", abs_tol=0.5e-3_WP)
end subroutine test_r_cv

subroutine test_gamma_cv(tests)
    use gasdata, only: AR, CO2
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)        :: cv
    type(si_temperature) :: temp
    type(unitless)       :: y(2), gamma_cv
    real(WP)             :: gamma_expected
    
    allocate(cv%gas(2))
    allocate(cv%m(2))
    
    cv%gas(1) = AR
    cv%gas(2) = CO2
    
    call y(1)%v%init_const(0.5_WP, 0)
    call y(2)%v%init_const(0.5_WP, 0)
    call temp%v%init_const(300.0_WP, 0)
    
    gamma_cv = cv%gamma(y)
    ! `AR` specific heat values from NIST; see where `AR` is defined.
    ! CO2 specific heat values from moran_fundamentals_2008 table A-20 at 300 K.
    gamma_expected = (0.52154_WP + 0.846_WP) / (0.31239_WP + 0.657_WP)
    call tests%real_eq(gamma_cv%v%v, gamma_expected, "cv%gamma for mixture (2)", abs_tol=1.0e-4_WP)
end subroutine test_gamma_cv

subroutine test_p_c(tests)
    ! Based on moran_fundamentals_2008 example 11.10, pp. 615--617.
    
    use gasdata, only: gas_type
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)     :: cv
    real(WP)          :: n(2), m(2)
    type(si_pressure) :: p_c_cv
    
    allocate(cv%gas(2))
    
    ! > 0.18 kmol of methane (CH4) and 0.274 of butane (C4H10)
    ! methane (CH4)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/CH4/h1H4> ("Fluid Properties" at 0.101325 MPa)
    cv%gas(1) = gas_type(gamma = 2.2359_WP/1.7129_WP, &
                         u_0   = 758.87e3_WP, &
                         h_0   = 914.09e3_WP, &
                         mm    = 16.04e-3_WP, & ! moran_fundamentals_2008 table A-1
                         p_c   = 46.4e5_WP) ! moran_fundamentals_2008 table A-1
    
    ! butane (C4H10)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/C4H10/c1-3-4-2/h3-4H2%2C1-2H3> ("Fluid Properties" at 0.101325 MPa)
    cv%gas(2) = gas_type(gamma = 1.7388_WP/1.5744_WP, &
                         u_0   = 589.09e3_WP, &
                         h_0   = 630.74e3_WP, &
                         mm    = 58.12e-3_WP, & ! moran_fundamentals_2008 table A-1
                         p_c   = 38.0e5_WP) ! moran_fundamentals_2008 table A-1
    
    n(1) = 0.18e3_WP
    n(2) = 0.274e3_WP
    m(1) = n(1)*cv%gas(1)%mm
    m(2) = n(2)*cv%gas(2)%mm
    
    allocate(cv%m(2))
    call cv%m(1)%v%init_const(m(1), 0)
    call cv%m(2)%v%init_const(m(2), 0)
    
    ! p. 616: critical pressure
    p_c_cv = cv%p_c()
    call tests%real_eq(p_c_cv%v%v, 41.33e5_WP, "cv%p_c", abs_tol=50.0_WP)
    
    ! I was going to test `cv%gamma()` for this case. However, something is off about butane.
    ! The specific heats calculated from moran_fundamentals_2008 eqs. 3.47a-3.47b are about 13% off what NIST says.
end subroutine test_p_c

subroutine test_p_f_1(tests)
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)     :: cv
    type(si_pressure) :: p_fe, p_f
    
    call cv%x%v%init_const(0.0_WP, 0)
    call cv%p_fs%v%init_const(0.1e5_WP, 0)
    call cv%p_fd%v%init_const(0.05e5_WP, 0)
    
    call cv%x_dot%v%init_const(-10.0_WP, 0)
    call p_fe%v%init_const(-1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, -cv%p_fd%v%v, "p_f, x_dot < 0 (dynamic friction), p_fe < 0")
    
    call cv%x_dot%v%init_const(-10.0_WP, 0)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, -cv%p_fd%v%v, "p_f, x_dot < 0 (dynamic friction), p_fe == 0")
    
    call cv%x_dot%v%init_const(-10.0_WP, 0)
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, -cv%p_fd%v%v, "p_f, x_dot < 0 (dynamic friction), p_fe > 0")
    
    call cv%x_dot%v%init_const(0.0_WP, 0)
    p_fe = -cv%p_fs/10.0_WP
    p_f  = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, p_fe%v%v, "p_f, x_dot == 0 (static friction), p_fe < 0")
    
    call cv%x_dot%v%init_const(0.0_WP, 0)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, 0.0_WP, "p_f, x_dot == 0 (static friction), p_fe == 0")
    
    call cv%x_dot%v%init_const(0.0_WP, 0)
    p_fe = cv%p_fs/10.0_WP
    p_f  = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, p_fe%v%v, "p_f, x_dot == 0 (static friction), p_fe > 0")
    
    call cv%x_dot%v%init_const(10.0_WP, 0)
    call p_fe%v%init_const(-1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, cv%p_fd%v%v, "p_f, x_dot > 0 (dynamic friction), p_fe < 0")
    
    call cv%x_dot%v%init_const(10.0_WP, 0)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, cv%p_fd%v%v, "p_f, x_dot > 0 (dynamic friction), p_fe == 0")
    
    call cv%x_dot%v%init_const(10.0_WP, 0)
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, cv%p_fd%v%v, "p_f, x_dot > 0 (dynamic friction), p_fe > 0")
end subroutine test_p_f_1

subroutine test_p_f_2(tests)
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)     :: cv
    type(si_pressure) :: p_fe, p_f
    
    call cv%x%v%init_const(0.0_WP, 0)
    call cv%p_fs%v%init_const(0.0_WP, 0)
    call cv%p_fd%v%init_const(0.0_WP, 0)
    
    call cv%x_dot%v%init_const(10.0_WP, 0)
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, 0.0_WP, "p_f, p_fs = p_fd = 0, x_dot > 0")
    
    call cv%x_dot%v%init_const(0.0_WP, 0)
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, 0.0_WP, "p_f, p_fs = p_fd = 0, x_dot == 0")
end subroutine test_p_f_2

subroutine test_p_f0_1(tests)
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)     :: cv
    type(si_pressure) :: p_fe, p_f0
    
    call cv%x%v%init_const(0.0_WP, 0)
    call cv%p_fs%v%init_const(0.1e5_WP, 0)
    call cv%p_fd%v%init_const(0.05e5_WP, 0)
    
    call p_fe%v%init_const(-1.0e5_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, -cv%p_fs%v%v, "p_f0 (1)", abs_tol=1.0e-4_WP)
    
    p_fe = -cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, p_fe%v%v, "p_f0 (2)")
    
    call p_fe%v%init_const(0.0_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (3)")
    
    p_fe = cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, p_fe%v%v, "p_f0 (4)")
    
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, cv%p_fs%v%v, "p_f0 (5)", abs_tol=1.0e-4_WP)
end subroutine test_p_f0_1

subroutine test_p_f0_2(tests)
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)     :: cv
    type(si_pressure) :: p_fe, p_f0
    
    call cv%x%v%init_const(0.0_WP, 0)
    call cv%p_fs%v%init_const(0.0_WP, 0)
    call cv%p_fd%v%init_const(0.0_WP, 0)
    
    call p_fe%v%init_const(-1.0e5_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (1)")
    
    p_fe = -cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (2)")
    
    call p_fe%v%init_const(0.0_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (3)")
    
    p_fe = cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (4)")
    
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (5)")
end subroutine test_p_f0_2

! TODO: Plot `p_f0` to test it.

subroutine test_temp_cv(tests)
    use gasdata, only: DRY_AIR
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)        :: cv
    type(si_temperature) :: temp
    
    ! 1 kg of mass
    allocate(cv%m(1))
    call cv%m(1)%v%init_const(1.0_WP, 0)
    cv%i_cv_other = 1
    
    ! air
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    ! moran_fundamentals_2008 table A-22, 300 K with 1 kg of mass
    call cv%e%v%init_const(214.07e3_WP, 0)
    
    temp = cv%temp()
    call tests%real_eq(temp%v%v, 300.0_WP, "temp_cv (qualitative)", abs_tol=5.0_WP)
end subroutine test_temp_cv

subroutine test_set_1(tests)
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM, RHO_ATM, DRY_AIR
    use cva, only: X_STOP_DEFAULT, cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type) :: cv
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1)
    type(si_pressure)        :: p, p_cv
    type(si_temperature)     :: temp, temp_cv
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_pressure)        :: p_fs, p_fd, p_atm
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    type(si_specific_energy) :: u
    
    type(si_volume)       :: vol_cv
    type(si_mass_density) :: rho_cv
    
    call x%v%init_const(0.5_WP, 0)
    call x_dot%v%init_const(0.5_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p%v%init_const(2.0_WP*P_ATM_, 0)
    call temp%v%init_const(sqrt(2.0_WP)*TEMP_ATM, 0)
    call csa%v%init_const(0.1_WP, 0)
    call rm_p%v%init_const(1.0_WP/0.5_WP, 0)
    call p_fs%v%init_const(0.2e5_WP, 0)
    call p_fd%v%init_const(0.1e5_WP, 0)
    call p_atm%v%init_const(P_ATM_, 0)
    call k%v%init_const(10.0_WP, 0)
    call x_z%v%init_const(3.0_WP, 0)
    
    call cv%set(x, x_dot, y, p, temp, "test_set_1", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    
    call tests%real_eq(cv%x%v%v, x%v%v, "set 1, x")
    call tests%real_eq(cv%x_dot%v%v, x_dot%v%v, "set 1, x_dot")
    ! no `p` or `temp` member variables
    call tests%real_eq(cv%csa%v%v, csa%v%v, "set 1, csa")
    call tests%real_eq(cv%rm_p%v%v, rm_p%v%v, "set 1, rm_p")
    call tests%real_eq(cv%p_fs%v%v, p_fs%v%v, "set 1, p_fs")
    call tests%real_eq(cv%p_fd%v%v, p_fd%v%v, "set 1, p_fd")
    call tests%integer_eq(cv%i_cv_other, 2, "set 1, i_cv_other")
    
    call tests%real_eq(cv%k%v%v, k%v%v, "set 1, k")
    call tests%real_eq(cv%x_z%v%v, x_z%v%v, "set 1, x_z")
    call tests%real_eq(cv%x_stop%v%v, X_STOP_DEFAULT, "set 1, x_stop")
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, temp%v%v, "set 1, temp")
    
    vol_cv = cv%vol()
    call tests%real_eq(vol_cv%v%v, 0.05_WP, "set 1, vol")
    
    ! Not exact, based on multiple of ambient density for simplicity.
    rho_cv = cv%rho()
    call tests%real_eq(rho_cv%v%v, sqrt(2.0_WP)*RHO_ATM, "set 1, rho", abs_tol=1.0e-3_WP)
    
    p_cv = cv%p()
    call tests%real_eq(p_cv%v%v, p%v%v, "set 1, p")
    
    ! Not exact, based on multiple of ambient density for simplicity.
    u = DRY_AIR%u(temp)
    call tests%real_eq(cv%m(1)%v%v, sqrt(2.0_WP)*RHO_ATM*0.05_WP, "set 1, m(1)", abs_tol=1.0e-4_WP)
    call tests%real_eq(cv%e%v%v, u%v%v*sqrt(2.0_WP)*RHO_ATM*0.05_WP, "set 1, e", abs_tol=10.0_WP)
end subroutine test_set_1

subroutine test_set_2(tests)
    ! Intended to test gas mixtures.
    ! Based on moran_fundamentals_2008 example 11.10, pp. 615--617.
    
    use convert, only: celsius_const
    use gasdata, only: P_ATM_ => P_ATM, gas_type
    use cva, only: X_STOP_DEFAULT, cv_type
    use checks, only: assert, is_close
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type) :: cv
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(2), chi_cv(2), y_cv(2)
    type(si_pressure)     :: p, p_cv
    type(si_temperature)  :: temp, temp_cv
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_mass)         :: m_total
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: x_z, x_stop
    type(gas_type)        :: gas(2)
    
    real(WP)        :: n(2), m(2)
    type(si_volume) :: vol_cv
    
    ! > 0.18 kmol of methane (CH4) and 0.274 of butane (C4H10)
    ! methane (CH4)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/CH4/h1H4> ("Fluid Properties" at 0.101325 MPa)
    gas(1) = gas_type(gamma = 2.2359_WP/1.7129_WP, &
                      u_0   = 758.87e3_WP, &
                      h_0   = 914.09e3_WP, &
                      mm    = 16.04e-3_WP, & ! moran_fundamentals_2008 table A-1
                      p_c   = 10.0_WP*46.4e5_WP) ! moran_fundamentals_2008 table A-1, set larger to avoid ideal gas validity check
    
    ! butane (C4H10)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/C4H10/c1-3-4-2/h3-4H2%2C1-2H3> ("Fluid Properties" at 0.101325 MPa)
    gas(2) = gas_type(gamma = 1.7388_WP/1.5744_WP, &
                      u_0   = 589.09e3_WP, &
                      h_0   = 630.74e3_WP, &
                      mm    = 58.12e-3_WP, & ! moran_fundamentals_2008 table A-1
                      p_c   = 10.0_WP*38.0e5_WP) ! moran_fundamentals_2008 table A-1, set larger to avoid ideal gas validity check
    
    n(1) = 0.18e3_WP
    n(2) = 0.274e3_WP
    m(1) = n(1)*gas(1)%mm
    m(2) = n(2)*gas(2)%mm
    call y(1)%v%init_const(m(1) / (m(1) + m(2)), 0)
    call y(2)%v%init_const(m(2) / (m(1) + m(2)), 0)
    
    ! > volume of 0.241 m3
    call x%v%init_const(0.241_WP, 0)
    call csa%v%init_const(1.0_WP, 0)
    
    ! > temperature of 238 C
    temp = celsius_const(238.0_WP, 0)
    
    ! p. 615: > The experimental value for the pressure is 68.9 bar.
    ! But the ideal gas law is off here. I'm testing the ideal gas law right now so I'll use the ideal gas value.
    ! p. 616: > 80.01 bar
    call p%v%init_const(80.01e5_WP, 0)
    
    ! These are set to basically random values as they are irrelevant to this test.
    call x_dot%v%init_const(0.0_WP, 0)
    call rm_p%v%init_const(0.0_WP, 0)
    call p_fs%v%init_const(0.0_WP, 0)
    call p_fd%v%init_const(0.0_WP, 0)
    call p_atm%v%init_const(P_ATM_, 0)
    call k%v%init_const(0.0_WP, 0)
    call x_z%v%init_const(0.0_WP, 0)
    call x_stop%v%init_const(1.0_WP, 0) ! to test non-default `x_stop`
    call assert(.not. is_close(x_stop%v%v, X_STOP_DEFAULT), "x_stop must not equal the default for this test")
    
    call cv%set(x, x_dot, y, p, temp, "test_set_2", csa, rm_p, p_fs, p_fd, k, x_z, gas, 2, x_stop)
    
    ! All the masses are slightly off. This is expected, as total mass is not an input here.
    ! Pressure is the input setting the total mass, and it's only approximate.
    call tests%real_eq(cv%m(1)%v%v, m(1), "set 2, cv%m(1)", abs_tol=1.0e-1_WP)
    call tests%real_eq(cv%m(2)%v%v, m(2), "set 2, cv%m(2)", abs_tol=1.0e-1_WP)
    m_total = cv%m_total()
    call tests%real_eq(m_total%v%v, m(1) + m(2), "set 2, cv%m_total", abs_tol=1.0e-1_WP)
    m_total = cv%vol() * cv%rho_eos(p, temp, y)
    call tests%real_eq(m_total%v%v, m(1) + m(2), "set 2, alt m_total", abs_tol=1.0e-1_WP)
    
    ! p. 616: mole fractions, book is probably only accurate to 3 decimal points
    chi_cv = cv%chi()
    call tests%real_eq(chi_cv(1)%v%v, 0.396_WP, "set 2, chi_cv(1)", abs_tol=1.0e-3_WP)
    call tests%real_eq(chi_cv(2)%v%v, 0.604_WP, "set 2, chi_cv(2)", abs_tol=1.0e-3_WP)
    
    ! mass fractions
    y_cv = cv%y()
    call tests%real_eq(y_cv(1)%v%v, y(1)%v%v, "set 2, y_cv(1)")
    call tests%real_eq(y_cv(2)%v%v, y(2)%v%v, "set 2, y_cv(2)")
    
    vol_cv = cv%vol()
    call tests%real_eq(vol_cv%v%v, 0.241_WP, "set 2, temp")
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, temp%v%v, "set 2, temp")
    
    p_cv = cv%p()
    call tests%real_eq(p_cv%v%v, p%v%v, "set 2, p")
    
    call tests%real_eq(cv%x_stop%v%v, 1.0_WP, "set 2, x_stop")
end subroutine test_set_2

subroutine test_set_3(tests)
    ! Testing isentropic filling
    
    use gasdata, only: DRY_AIR
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type) :: cv
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_pressure)     :: p
    type(si_temperature)  :: temp_atm, temp_cv
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: x_z
    
    call x%v%init_const(0.5_WP, 0)
    call x_dot%v%init_const(0.5_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p%v%init_const(2.0e5_WP, 0)
    call temp_atm%v%init_const(300.0_WP, 0)
    call csa%v%init_const(0.1_WP, 0)
    call rm_p%v%init_const(1.0_WP/0.5_WP, 0)
    call p_fs%v%init_const(0.2e5_WP, 0)
    call p_fd%v%init_const(0.1e5_WP, 0)
    call p_atm%v%init_const(1.0e5_WP, 0)
    call k%v%init_const(10.0_WP, 0)
    call x_z%v%init_const(3.0_WP, 0)
    
    call cv%set(x, x_dot, y, p, temp_atm, "test_set_3", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2, &
                    isentropic_filling=.true., p_atm=p_atm)
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, 300.0_WP*(2.0_WP**(0.4_WP/1.4_WP)), "set 3, isentropic_filling=.true., cv%temp")
end subroutine test_set_3

subroutine test_rates(tests)
    use gasdata, only: DRY_AIR
    use cva, only: cv_system_type, d_x_d_t, d_xdot_d_t, d_m_k_d_t, d_e_d_t
    
    type(test_results_type), intent(in out) :: tests

    type(cv_system_type) :: sys
    
    type(si_length)           :: x
    type(si_velocity)         :: x_dot, d_x_d_t_
    type(unitless)            :: y(1)
    type(si_pressure)         :: p, p_fs, p_fd, p_atm
    type(si_temperature)      :: temp
    type(si_area)             :: csa
    type(si_inverse_mass)     :: rm_p
    type(si_stiffness)        :: k
    type(si_length)           :: x_z
    type(si_acceleration)     :: d_xdot_d_t_
    type(si_mass_flow_rate)   :: m_dots(2, 2), d_m_1_d_t
    type(si_energy_flow_rate) :: h_dots(2, 2), d_e_d_t_
    
    call x%v%init_const(1.5_WP, 0)
    call x_dot%v%init_const(10.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p%v%init_const(12.0e5_WP, 0)
    call temp%v%init_const(300.0_WP, 0)
    call csa%v%init_const(4.0_WP, 0)
    call rm_p%v%init_const(1.0_WP/2.0_WP, 0)
    call p_fs%v%init_const(2.0e5_WP, 0)
    call p_fd%v%init_const(1.0e5_WP, 0)
    call p_atm%v%init_const(1.0e5_WP, 0)
    call k%v%init_const(1.0e6_WP, 0)
    call x_z%v%init_const(0.5_WP, 0)
    
    allocate(sys%cv(2))
    
    call sys%cv(1)%set_const("test_rates 1 (constant pressure)", p_atm, temp, [DRY_AIR])
    call sys%cv(2)%set(x, x_dot, y, p, temp, "test_rates 2", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 1)
    
    d_x_d_t_ = d_x_d_t(sys%cv(2))
    call tests%real_eq(d_x_d_t_%v%v, x_dot%v%v, "d_x_d_t")
    
    d_xdot_d_t_ = d_xdot_d_t(sys, 2)
    call tests%real_eq(d_xdot_d_t_%v%v, 1.5e6_WP, "d_xdot_d_t")
    
    call m_dots(1, 1)%v%init_const(0.0_WP, 0)
    call m_dots(1, 2)%v%init_const(2.0_WP, 0)
    call m_dots(2, 1)%v%init_const(0.0_WP, 0)
    call m_dots(2, 2)%v%init_const(0.0_WP, 0)
    
    d_m_1_d_t = d_m_k_d_t(sys%cv(2), m_dots, 1, 1)
    call tests%real_eq(d_m_1_d_t%v%v, -2.0_WP, "d_m_d_t(1)")
    
    d_m_1_d_t = d_m_k_d_t(sys%cv(2), m_dots, 1, 2)
    call tests%real_eq(d_m_1_d_t%v%v, 2.0_WP, "d_m_d_t(2)")
    
    call h_dots(1, 1)%v%init_const(0.0_WP, 0)
    call h_dots(1, 2)%v%init_const(0.0_WP, 0)
    call h_dots(2, 1)%v%init_const(2.0e5_WP, 0)
    call h_dots(2, 2)%v%init_const(0.0_WP, 0)
    
    d_e_d_t_ = d_e_d_t(sys%cv(2), h_dots, 1)
    call tests%real_eq(d_e_d_t_%v%v, -p%v%v*csa%v%v*x_dot%v%v + 2.0e5_WP, "d_e_d_t(1)")
    
    d_e_d_t_ = d_e_d_t(sys%cv(2), h_dots, 2)
    call tests%real_eq(d_e_d_t_%v%v, -p%v%v*csa%v%v*x_dot%v%v - 2.0e5_WP, "d_e_d_t(2)")
end subroutine test_rates

subroutine test_u_h_cv(tests)
    use gasdata, only: N2, H2O
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type) :: cv
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(2)
    type(si_pressure)        :: p, p_fs, p_fd, p_atm
    type(si_temperature)     :: temp
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    type(si_specific_energy) :: u, h
    
    real(WP) :: u_exact, h_exact
    
    call x%v%init_const(1.5_WP, 0)
    call x_dot%v%init_const(10.0_WP, 0)
    call y(1)%v%init_const(0.25_WP, 0)
    call y(2)%v%init_const(0.75_WP, 0)
    call p%v%init_const(12.0e5_WP, 0)
    call temp%v%init_const(300.0_WP, 0)
    call csa%v%init_const(4.0_WP, 0)
    call rm_p%v%init_const(1.0_WP/2.0_WP, 0)
    call p_fs%v%init_const(2.0e5_WP, 0)
    call p_fd%v%init_const(1.0e5_WP, 0)
    call p_atm%v%init_const(1.0e5_WP, 0)
    call k%v%init_const(1.0e6_WP, 0)
    call x_z%v%init_const(0.5_WP, 0)
    
    call cv%set(x, x_dot, y, p, temp, "test_u_h_cv", csa, rm_p, p_fs, p_fd, k, x_z, [N2, H2O], 2)
    
    u_exact = y(1)%v%v*N2%u_0 + y(2)%v%v*H2O%u_0
    h_exact = y(1)%v%v*N2%h_0 + y(2)%v%v*H2O%h_0
    
    u = cv%u()
    h = cv%h()
    
    ! Data from moran_fundamentals_2008 table A-22.
    call tests%real_eq(u%v%v, u_exact, "u (cv), 300 K")
    call tests%real_eq(h%v%v, h_exact, "h (cv), 300 K")
end subroutine test_u_h_cv

subroutine test_smooth_min(tests)
    use cva, only: smooth_min
    
    type(test_results_type), intent(in out) :: tests
    
    type(unitless) :: x, y, z
    
    call x%v%init_const(0.0_WP, 0)
    call y%v%init_const(1.0_WP, 0)
    z = smooth_min(x, y)
    call tests%real_eq(z%v%v, 0.0_WP, "smooth_min (1)")
    
    call x%v%init_const(1.1_WP, 0)
    call y%v%init_const(0.1_WP, 0)
    z = smooth_min(x, y)
    call tests%real_eq(z%v%v, 0.1_WP, "smooth_min (1)")
end subroutine test_smooth_min

subroutine test_f_m_dot(tests)
    use cva, only: f_m_dot
    
    type(test_results_type), intent(in out) :: tests
    
    type(unitless) :: p_r, b, f
    
    call b%v%init_const(0.5_WP, 0)
    
    call p_r%v%init_const(0.0_WP, 0)
    f = f_m_dot(p_r, b)
    call tests%real_eq(f%v%v, 0.0_WP, "f_m_dot (1)")
    
    call p_r%v%init_const(1.0_WP, 0)
    f = f_m_dot(p_r, b)
    call tests%real_eq(f%v%v, 1.0_WP, "f_m_dot (2)", abs_tol=0.04_WP)
    call tests%real_lt(f%v%v, 1.0_WP, "f_m_dot (3)")
end subroutine test_f_m_dot

! TODO: Plot `f_m_dot` to test it.

subroutine test_g_m_dot(tests)
    use cva, only: P_RL, g_m_dot
    
    type(test_results_type), intent(in out) :: tests
    
    type(unitless) :: p_r, g
    
    call p_r%v%init(0.0_WP, 1, 1)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 0.0_WP, "g_m_dot (value, 1)")
    call tests%real_eq(g%v%d(1), 0.0_WP, "g_m_dot (derivative, 1)")
    
    call p_r%v%init(0.9_WP, 1, 1)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 0.0_WP, "g_m_dot (value, 2)")
    call tests%real_eq(g%v%d(1), 0.0_WP, "g_m_dot (derivative, 2)")
    
    call p_r%v%init(P_RL, 1, 1)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 0.0_WP, "g_m_dot (value, 3)")
    call tests%real_eq(g%v%d(1), 0.0_WP, "g_m_dot (derivative, 3)")
    
    call p_r%v%init(0.5_WP*(P_RL + 1.0_WP), 1, 1)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 0.5_WP, "g_m_dot (value, 4)", abs_tol=1.0e-12_WP)
    call tests%real_eq(g%v%d(1), 3.0_WP/(2.0_WP*(1.0_WP - P_RL)), "g_m_dot (derivative, 4)")
    
    call p_r%v%init(1.0_WP, 1, 1)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 1.0_WP, "g_m_dot (value, 5)")
    call tests%real_eq(g%v%d(1), 0.0_WP, "g_m_dot (derivative, 5)")
end subroutine test_g_m_dot

! TODO: Plot `g_m_dot` to test it.

subroutine test_m_dot_1(tests)
    ! Test to check that $\dot{m} \varpropto \Delta p$ at small $\Delta p$.
    ! Tests the laminar branch (`cv_from%p()/cv_to%p() >= P_RL`).
    
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM, DRY_AIR, R_BAR
    use cva, only: P_RL, cv_type, con_type
    
    type(test_results_type), intent(in out) :: tests
    
    type(cv_type)     :: cv_from, cv_to
    type(con_type)    :: con
    type(si_pressure) :: delta_p
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1)
    type(si_pressure)        :: p
    type(si_temperature)     :: temp
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_pressure)        :: p_fs, p_fd, p_atm
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    
    type(si_mass_flow_rate) :: m_dot_con
    real(WP) :: d_m_dot_d_delta_p
    
    call delta_p%v%init(0.0_WP, 1, 1)
    
    con%active = .true.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    
    call x%v%init_const(0.5_WP, 1)
    call x_dot%v%init_const(0.5_WP, 1)
    call y(1)%v%init_const(1.0_WP, 1)
    call p%v%init_const(2.0_WP*P_ATM_, 1)
    call temp%v%init_const(sqrt(2.0_WP)*TEMP_ATM, 1)
    call csa%v%init_const(1.0_WP, 1)
    call rm_p%v%init_const(0.0_WP, 1)
    call p_fs%v%init_const(0.0_WP, 1)
    call p_fd%v%init_const(0.0_WP, 1)
    call p_atm%v%init_const(P_ATM_, 1)
    call k%v%init_const(0.0_WP, 1)
    call x_z%v%init_const(0.0_WP, 1)
    
    call cv_from%set(x, x_dot, y, p, temp, "from", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p - delta_p, temp, "to", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(cv_from, cv_to)
    
    ! It is important that `m_dot` goes to zero at zero pressure difference, so this is exact.
    call tests%real_eq(m_dot_con%v%v, 0.0_WP, "m_dot, small delta_p, v")
    
    ! This won't be exact due to the smoothing I applied.
    ! The important part is that it is finite.
    d_m_dot_d_delta_p = con%a_e%v%v * sqrt((1.0_WP - con%b%v%v) / ((R_BAR/DRY_AIR%mm) * sqrt(2.0_WP)*TEMP_ATM)) &
                            * sqrt(1.0_WP - ((P_RL - con%b%v%v) / (1.0_WP - con%b%v%v))**2)
    call tests%real_eq(m_dot_con%v%d(1), d_m_dot_d_delta_p, "m_dot, small delta_p, d", abs_tol=1.0e-3_WP)
end subroutine test_m_dot_1

subroutine test_m_dot_2(tests)
    ! Tests the subsonic branch (`P_RL > cv_from%p()/cv_to%p() > b`).
    
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM, DRY_AIR, R_BAR
    use cva, only: P_RL, cv_type, con_type
    use checks, only: assert
    
    type(test_results_type), intent(in out) :: tests
    
    type(cv_type)  :: cv_from, cv_to
    type(con_type) :: con
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1), p_r
    type(si_pressure)        :: p_in, p_out
    type(si_temperature)     :: temp
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_pressure)        :: p_fs, p_fd, p_atm
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    
    type(si_mass_flow_rate) :: m_dot_con
    real(WP) :: m_dot_con_, d_m_dot_d_p_r
    
    con%active = .true.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    
    call p_r%v%init(0.5_WP*(P_RL + con%b%v%v), 1, 1)
    call assert(p_r%v%v < P_RL,      "test_cva (test_m_dot_2): p_r < P_RL violated")
    call assert(p_r%v%v > con%b%v%v, "test_cva (test_m_dot_2): p_r > b violated")
    
    call x%v%init_const(0.5_WP, 1)
    call x_dot%v%init_const(0.5_WP, 1)
    call y(1)%v%init_const(1.0_WP, 1)
    call p_in%v%init_const(2.0_WP*P_ATM_, 1)
    call temp%v%init_const(sqrt(2.0_WP)*TEMP_ATM, 1)
    call csa%v%init_const(1.0_WP, 1)
    call rm_p%v%init_const(0.0_WP, 1)
    call p_fs%v%init_const(0.0_WP, 1)
    call p_fd%v%init_const(0.0_WP, 1)
    call p_atm%v%init_const(P_ATM_, 1)
    call k%v%init_const(0.0_WP, 1)
    call x_z%v%init_const(0.0_WP, 1)
    
    p_out = p_in*p_r
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(cv_from, cv_to)
    
    ! These should not be exact due to the smoothing applied, but they are still very close!
    
    m_dot_con_ = con%a_e%v%v * p_in%v%v * sqrt((1.0_WP - con%b%v%v)/((R_BAR/DRY_AIR%mm) * temp%v%v)) &
                        * sqrt(1.0_WP-((p_r%v%v-con%b%v%v)/(1.0_WP-con%b%v%v))**2)
    call tests%real_eq(m_dot_con%v%v, m_dot_con_, "m_dot, subsonic, v", abs_tol=1.0e-10_WP)
    
    d_m_dot_d_p_r = -con%a_e%v%v * p_in%v%v * sqrt((1.0_WP - con%b%v%v)/((R_BAR/DRY_AIR%mm) * temp%v%v)) &
                        * ((p_r%v%v-con%b%v%v) / (((1.0_WP-con%b%v%v)**2) &
                                * sqrt(1.0_WP-((p_r%v%v-con%b%v%v)/(1.0_WP-con%b%v%v))**2)))
    call tests%real_eq(m_dot_con%v%d(1), d_m_dot_d_p_r, "m_dot, subsonic, d", abs_tol=1.0e-8_WP)
end subroutine test_m_dot_2

subroutine test_m_dot_3(tests)
    ! Tests the choked branch (`cv_from%p()/cv_to%p() <= b`).
    
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM, DRY_AIR, R_BAR
    use cva, only: cv_type, con_type
    use checks, only: assert
    
    type(test_results_type), intent(in out) :: tests
    
    type(cv_type)  :: cv_from, cv_to
    type(con_type) :: con
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1), p_r
    type(si_pressure)        :: p_in, p_out
    type(si_temperature)     :: temp
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_pressure)        :: p_fs, p_fd, p_atm
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    
    type(si_mass_flow_rate) :: m_dot_con
    real(WP) :: m_dot_con_
    
    con%active = .true.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    
    call p_r%v%init(1.0e-6_WP, 1, 1)
    call assert(p_r%v%v < con%b%v%v, "test_cva (test_m_dot_3): p_r < b violated")
    
    call x%v%init_const(0.5_WP, 1)
    call x_dot%v%init_const(0.5_WP, 1)
    call y(1)%v%init_const(1.0_WP, 1)
    call p_in%v%init_const(2.0_WP*P_ATM_, 1)
    call temp%v%init_const(sqrt(2.0_WP)*TEMP_ATM, 1)
    call csa%v%init_const(1.0_WP, 1)
    call rm_p%v%init_const(0.0_WP, 1)
    call p_fs%v%init_const(0.0_WP, 1)
    call p_fd%v%init_const(0.0_WP, 1)
    call p_atm%v%init_const(P_ATM_, 1)
    call k%v%init_const(0.0_WP, 1)
    call x_z%v%init_const(0.0_WP, 1)
    
    p_out = p_in*p_r
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(cv_from, cv_to)
    
    m_dot_con_ = con%a_e%v%v * p_in%v%v * sqrt((1.0_WP - con%b%v%v)/((R_BAR/DRY_AIR%mm) * temp%v%v))
    call tests%real_eq(m_dot_con%v%v, m_dot_con_, "m_dot, choked, v") !, abs_tol=1.0e-10_WP)
    
    call tests%real_eq(m_dot_con%v%d(1), 0.0_WP, "m_dot, choked, d") !, abs_tol=1.0e-8_WP)
end subroutine test_m_dot_3

subroutine test_m_dot_4(tests)
    ! Test inactive connection.
    
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM, DRY_AIR
    use cva, only: cv_type, con_type
    use checks, only: assert
    
    type(test_results_type), intent(in out) :: tests
    
    type(cv_type)  :: cv_from, cv_to
    type(con_type) :: con
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1)
    type(si_pressure)        :: p_in, p_out
    type(si_temperature)     :: temp
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_pressure)        :: p_fs, p_fd, p_atm
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    
    type(si_mass_flow_rate) :: m_dot_con
    
    con%active = .false.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    
    call x%v%init_const(0.5_WP, 1)
    call x_dot%v%init_const(0.5_WP, 1)
    call y(1)%v%init_const(1.0_WP, 1)
    call p_in%v%init_const(2.0_WP*P_ATM_, 1)
    call p_out%v%init_const(2.0_WP*P_ATM_, 1)
    call temp%v%init_const(sqrt(2.0_WP)*TEMP_ATM, 1)
    call csa%v%init_const(1.0_WP, 1)
    call rm_p%v%init_const(0.0_WP, 1)
    call p_fs%v%init_const(0.0_WP, 1)
    call p_fd%v%init_const(0.0_WP, 1)
    call p_atm%v%init_const(P_ATM_, 1)
    call k%v%init_const(0.0_WP, 1)
    call x_z%v%init_const(0.0_WP, 1)
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(cv_from, cv_to)
    
    call tests%real_eq(m_dot_con%v%v, 0.0_WP, "m_dot, inactive, value")
    call tests%real_eq(m_dot_con%v%d(1), 0.0_WP, "m_dot, inactive, derivative")
end subroutine test_m_dot_4

subroutine test_calculate_flows(tests)
    use gasdata, only: P_ATM_ => P_ATM, R_BAR, DRY_AIR
    use cva, only: cv_system_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_system_type) :: sys
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y(1), p_r
    type(si_pressure)        :: p_in, p_out
    type(si_temperature)     :: temp
    type(si_area)            :: csa
    type(si_inverse_mass)    :: rm_p
    type(si_pressure)        :: p_fs, p_fd, p_atm
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    real(WP)                 :: m_dot_12
    type(si_specific_energy) :: h_1
    
    type(si_mass_flow_rate), allocatable   :: m_dot(:, :)
    type(si_energy_flow_rate), allocatable :: h_dot(:, :)
    
    allocate(sys%cv(2))
    allocate(sys%con(2, 2))
    
    sys%con(1, 1)%active = .false.
    sys%con(2, 1)%active = .false.
    sys%con(2, 2)%active = .false.
    
    sys%con(1, 2)%active = .true.
    call sys%con(1, 2)%a_e%v%init_const(0.25_WP, 0)
    call sys%con(1, 2)%b%v%init_const(0.5_WP, 0)
    
    call p_r%v%init_const(1.0e-6_WP, 0)
    
    call x%v%init_const(0.5_WP, 0)
    call x_dot%v%init_const(0.5_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p_in%v%init_const(2.0_WP*P_ATM_, 0)
    call temp%v%init_const(300.0_WP, 0)
    call csa%v%init_const(1.0_WP, 0)
    call rm_p%v%init_const(0.0_WP, 0)
    call p_fs%v%init_const(0.0_WP, 0)
    call p_fd%v%init_const(0.0_WP, 0)
    call p_atm%v%init_const(P_ATM_, 0)
    call k%v%init_const(0.0_WP, 0)
    call x_z%v%init_const(0.0_WP, 0)
    
    p_out = p_in*p_r
    
    call sys%cv(1)%set(x, x_dot, y, p_in, temp, "1", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 3)
    call sys%cv(2)%set(x, x_dot, y, p_out, temp, "2", csa, rm_p, p_fs, p_fd, k, x_z, [DRY_AIR], 3)
    
    m_dot_12 = sys%con(1, 2)%a_e%v%v * p_in%v%v * sqrt((1.0_WP - sys%con(1, 2)%b%v%v)/((R_BAR/DRY_AIR%mm) * temp%v%v))
    
    call sys%calculate_flows(m_dot, h_dot)
    
    call tests%integer_eq(size(m_dot, 1), 2, "size(m_dot, 1)")
    call tests%integer_eq(size(m_dot, 2), 2, "size(m_dot, 2)")
    call tests%integer_eq(size(h_dot, 1), 2, "size(h_dot, 1)")
    call tests%integer_eq(size(h_dot, 2), 2, "size(h_dot, 2)")
    
    call tests%real_eq(m_dot(1, 1)%v%v, 0.0_WP, "m_dots(1, 1)")
    call tests%real_eq(m_dot(1, 2)%v%v, m_dot_12, "m_dots(1, 2)")
    call tests%real_eq(m_dot(2, 1)%v%v, 0.0_WP, "m_dots(2, 1)")
    call tests%real_eq(m_dot(2, 2)%v%v, 0.0_WP, "m_dots(2, 2)")
    
    h_1 = DRY_AIR%h(temp)
    call tests%real_eq(h_dot(1, 1)%v%v, 0.0_WP, "h_dots(1, 1)")
    call tests%real_eq(h_dot(1, 2)%v%v, h_1%v%v*m_dot_12, "h_dots(1, 2)")
    call tests%real_eq(h_dot(2, 1)%v%v, 0.0_WP, "m_dots(2, 1)")
    call tests%real_eq(h_dot(2, 2)%v%v, 0.0_WP, "m_dots(2, 2)")
end subroutine test_calculate_flows

subroutine test_conservation(tests)
    use convert
    use gasdata, only: DRY_AIR
    use prec, only: PI
    use cva, only: cv_system_type, run_status_type, run
    
    type(test_results_type), intent(in out) :: tests
    
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    type(si_length)       :: d_e, x_1, x_2, d_1, d_2, x_stop_2
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_pressure)     :: p_atm, p_1, p_2, p_fs_1, p_fd_1, p_fs_2, p_fd_2
    type(si_temperature)  :: temp_atm
    type(si_area)         :: csa_1, csa_2
    type(si_mass)         :: m_p_1, m_p_2, m_start, m_end
    type(si_stiffness)    :: k
    type(si_length)       :: x_z
    type(si_energy)       :: e_start, e_end, spring_pe_2_start, m_p_ke_2_start, e_start_2
    
    allocate(sys_start)
    allocate(sys_start%cv(3))
    allocate(sys_start%con(3, 3))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(1, 2)%active = .false.
    sys_start%con(1, 3)%active = .false.
    
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    sys_start%con(2, 3)%active = .true.
    d_e = inch_const(0.1_WP, 0)
    sys_start%con(2, 3)%a_e = (PI/4.0_WP)*square(d_e)
    call sys_start%con(2, 3)%b%v%init_const(0.5_WP, 0)
    
    sys_start%con(3, 1)%active = .false.
    sys_start%con(3, 2) = sys_start%con(2, 3)
    sys_start%con(3, 3)%active = .false.
    
    ! The same for every control volume.
    call y(1)%v%init_const(1.0_WP, 0)
    call temp_atm%v%init_const(300.0_WP, 0)
    
    ! 1: atmosphere
    call p_atm%v%init_const(0.0_WP, 0)
    
    call sys_start%cv(1)%set_const("atmosphere", p_atm, temp_atm, [DRY_AIR])
    
    ! 2: chamber
    call x_dot%v%init_const(-2.0_WP, 0)
    call d_1%v%init_const(2.0e-2_WP, 0)
    csa_1 = (PI/4.0_WP)*square(d_1)
    call x_1%v%init_const(10.0e-2_WP, 0)
    call p_1%v%init_const(5.0e5_WP, 0)
    call m_p_1%v%init_const(30.0e-3_WP, 0)
    call p_fs_1%v%init_const(0.0_WP, 0)
    call p_fd_1%v%init_const(0.0_WP, 0)
    call k%v%init_const(700.0_WP, 0)
    call x_z%v%init_const(1.0e-2_WP, 0)
    
    call sys_start%cv(2)%set(x_1, x_dot, y, p_1, temp_atm, "pressure chamber", csa_1, 1.0_WP/m_p_1, p_fs_1, p_fd_1, k, &
                                x_z, [DRY_AIR], 1)
    ! `isentropic_filling=.true.` requires that `p_atm > 0`, so it's not used here.
    
    ! 3: barrel
    
    call x_dot%v%init_const(0.0_WP, 0)
    call d_2%v%init_const(2.0e-2_WP, 0)
    csa_2 = (PI/4.0_WP)*square(d_2)
    call x_2%v%init_const(10.0e-2_WP, 0)
    call p_2%v%init_const(1.0e5_WP, 0)
    call m_p_2%v%init_const(1.0e-3_WP, 0)
    call p_fs_2%v%init_const(0.0_WP, 0)
    call p_fd_2%v%init_const(0.0_WP, 0)
    x_stop_2 = x_2 + inch_const(12.0_WP, 0)
    call k%v%init_const(0.0_WP, 0)
    call x_z%v%init_const(0.0_WP, 0)
    
    call sys_start%cv(3)%set(x_2, x_dot, y, p_2, temp_atm, "barrel", csa_2, 1.0_WP/m_p_2, p_fs_2, p_fd_2, k, &
                                x_z, [DRY_AIR], 1, x_stop=x_stop_2)
    
    call run(sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, 0, "test_conservation, status%rc")
    
    spring_pe_2_start = sys_start%cv(2)%spring_pe()
    call tests%real_eq(spring_pe_2_start%v%v, 0.5_WP*(700.0_WP)*(9.0e-2_WP)**2, "test_conservation, spring_pe_1_start")
    
    m_p_ke_2_start = sys_start%cv(2)%m_p_ke()
    call tests%real_eq(m_p_ke_2_start%v%v, 0.5_WP*(30.0e-3_WP)*(2.0_WP)**2, "test_conservation, m_p_ke_1_start")
    
    e_start_2 = sys_start%cv(2)%e_total()
    call tests%real_eq(e_start_2%v%v, sys_start%cv(2)%e%v%v + spring_pe_2_start%v%v + m_p_ke_2_start%v%v, &
                            "test_conservation, sys_start%cv(2)%e_total()")
    
    m_start = sys_start%m_total()
    m_end   = sys_end%m_total()
    call tests%real_eq(m_start%v%v, m_end%v%v, "test_conservation, m_start == m_end")
    
    e_start = sys_start%e_total()
    e_end   = sys_end%e_total()
    call tests%real_eq(e_start%v%v, e_end%v%v, "test_conservation, e_start == e_end", abs_tol=1.0e-5_WP)
end subroutine test_conservation

end program test_cva
