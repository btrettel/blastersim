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

real(WP), parameter :: TEST_EXACT_T_STOP = 0.0273_WP

call tests%start_tests("cva.nml")

call write_defaults()

call test_m_total(tests)
call test_p_eos_ideal(tests)
call test_constant_eos(tests)
call test_rho_eos(tests)
call test_r_cv(tests)
call test_p_c(tests)
call test_p_f_1(tests)
call test_p_f_2(tests)
call test_p_f0_1(tests)
call test_p_f0_2(tests)
call test_temp_cv_ideal(tests)
call test_set_normal_1(tests)
call test_set_normal_2(tests)
call test_set_normal_3(tests)
call test_set_const(tests)
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
call test_alpha_m_dot(tests)

call test_calculate_flows(tests)
call test_conservation_1(tests)
call test_conservation_2(tests)
call test_mirror_1(tests)
call test_mirror_2(tests)
call test_exact(tests)

call test_check_sys(tests)

call tests%end_tests()

contains

subroutine write_defaults()
    use cva, only: DT_DEFAULT, MASS_TOLERANCE, ENERGY_TOLERANCE, MASS_DERIV_TOLERANCE, ENERGY_DERIV_TOLERANCE, &
                    MIRROR_X_TOLERANCE
    use gasdata, only: P_ATM, TEMP_ATM
    use convert, only: CONVERT_C_TO_K
    use io, only: write_latex_engineering
    
    integer :: tex_unit
    
    open(newunit=tex_unit, action="write", status="replace", position="rewind", file="defaults.tex", delim="quote")
    write(unit=tex_unit, fmt="(a)") "% auto-generated"
    call write_latex_engineering(tex_unit, DT_DEFAULT, "dtdefault", "f4.1")
    write(unit=tex_unit, fmt="(a, f8.1, a)") "\newcommand*{\patmdefault}{", P_ATM, "}"
    write(unit=tex_unit, fmt="(a, f6.2, a)") "\newcommand*{\tempatmdefaultk}{", TEMP_ATM, "}"
    write(unit=tex_unit, fmt="(a, f5.2, a)") "\newcommand*{\tempatmdefaultc}{", TEMP_ATM - CONVERT_C_TO_K, "}"
    write(unit=tex_unit, fmt="(a, f5.3, a)") "\newcommand*{\masstolerance}{", 100.0_WP*MASS_TOLERANCE, "\%}"
    write(unit=tex_unit, fmt="(a, f4.2, a)") "\newcommand*{\energytolerance}{", 100.0_WP*ENERGY_TOLERANCE, "\%}"
    call write_latex_engineering(tex_unit, MASS_DERIV_TOLERANCE, "massderivtolerance", "f4.1")
    call write_latex_engineering(tex_unit, ENERGY_DERIV_TOLERANCE, "energyderivtolerance", "f4.1")
    call write_latex_engineering(tex_unit, MIRROR_X_TOLERANCE, "mirrorxtolerance", "f3.1")
    close(tex_unit)
end subroutine write_defaults

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

subroutine test_p_eos_ideal(tests)
    use gasdata, only: P_ATM, TEMP_ATM, RHO_ATM, DRY_AIR
    use cva, only: cv_type, NORMAL_CV_TYPE, IDEAL_EOS
    
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
    cv%i_cv_mirror = 2
    cv%type = NORMAL_CV_TYPE
    cv%eos  = IDEAL_EOS
    
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    p = cv%p_eos(rho, temp)
    
    call tests%real_eq(p%v%v, P_ATM, "p_eos, normal CV, atmospheric", abs_tol=20.0_WP)
end subroutine test_p_eos_ideal

subroutine test_constant_eos(tests)
    use gasdata, only: P_ATM, TEMP_ATM, DRY_AIR
    use cva, only: cv_type, CONST_EOS
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)         :: cv
    type(si_pressure)     :: p
    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    
    call cv%p_const%v%init_const(P_ATM, 0)
    call cv%temp_const%v%init_const(TEMP_ATM, 0)
    call rho%v%init_const(1.0_WP, 0) ! Doesn't matter what this is set to initially.
    allocate(cv%m(1))
    call cv%m(1)%v%init_const(4.0_WP, 0)
    call cv%e%v%init_const(1.0_WP, 0)
    call cv%x%v%init_const(2.0_WP, 0)
    call cv%csa%v%init_const(1.0_WP, 0)
    cv%eos = CONST_EOS
    
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    p    = cv%p_eos(rho, cv%temp_const)
    temp = cv%temp()
    rho  = cv%rho()
    
    call tests%real_eq(p%v%v, P_ATM, "p_eos, constant CV")
    call tests%real_eq(temp%v%v, TEMP_ATM, "temp_cv, constant CV")
    call tests%real_eq(rho%v%v, 2.0_WP, "rho_cv, constant CV") ! calculated from mass, not meaningful
end subroutine test_constant_eos

subroutine test_rho_eos(tests)
    use gasdata, only: P_ATM, TEMP_ATM, RHO_ATM, DRY_AIR
    use cva, only: cv_type, IDEAL_EOS
    
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
    cv%eos = IDEAL_EOS
    
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    rho = cv%rho_eos(p, temp, y)
    call tests%real_eq(rho%v%v, RHO_ATM, "rho_eos, atmospheric (1)", abs_tol=1.0e-3_WP)
end subroutine test_rho_eos

subroutine test_r_cv(tests)
    ! Test using ambient air.
    ! Tests both the `r_cv` and `rho_eos`
    
    use gasdata, only: P_ATM, TEMP_ATM, RHO_ATM, R_BAR, gas_type, DRY_AIR, N2, O2, AR, CO2
    use cva, only: cv_type, IDEAL_EOS
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
    
    cv%eos = IDEAL_EOS
    
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
    cv%gas(1) = gas_type(label = "CH4", &
                         gamma = 2.2359_WP/1.7129_WP, &
                         u_0   = 758.87e3_WP, &
                         h_0   = 914.09e3_WP, &
                         mm    = 16.04e-3_WP, & ! moran_fundamentals_2008 table A-1
                         p_c   = 46.4e5_WP) ! moran_fundamentals_2008 table A-1
    
    ! butane (C4H10)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/C4H10/c1-3-4-2/h3-4H2%2C1-2H3> ("Fluid Properties" at 0.101325 MPa)
    cv%gas(2) = gas_type(label = "C4H10", &
                         gamma = 1.7388_WP/1.5744_WP, &
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
    cv%constant_friction = .false.
    
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
    cv%constant_friction = .false.
    
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
    cv%constant_friction = .false.
    
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
    cv%constant_friction = .false.
    
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

subroutine test_temp_cv_ideal(tests)
    use gasdata, only: DRY_AIR
    use cva, only: cv_type, NORMAL_CV_TYPE, IDEAL_EOS
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)        :: cv
    type(si_temperature) :: temp
    
    ! 1 kg of mass
    allocate(cv%m(1))
    call cv%m(1)%v%init_const(1.0_WP, 0)
    cv%i_cv_mirror = 2
    cv%type = NORMAL_CV_TYPE
    cv%eos  = IDEAL_EOS
    
    ! air
    allocate(cv%gas(1))
    cv%gas(1) = DRY_AIR
    
    ! moran_fundamentals_2008 table A-22, 300 K with 1 kg of mass
    call cv%e%v%init_const(214.07e3_WP, 0)
    
    temp = cv%temp()
    call tests%real_eq(temp%v%v, 300.0_WP, "temp_cv (qualitative)", abs_tol=5.0_WP)
end subroutine test_temp_cv_ideal

subroutine test_set_normal_1(tests)
    use gasdata, only: P_ATM_ => P_ATM, TEMP_ATM, RHO_ATM, DRY_AIR
    use cva, only: X_STOP_DEFAULT, cv_type, NORMAL_CV_TYPE
    
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
    type(si_length)          :: delta_pre
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
    call delta_pre%v%init_const(3.0_WP, 0)
    
    call cv%set(x, x_dot, y, p, temp, "test_set_normal_1", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    
    call tests%real_eq(cv%x%v%v, x%v%v, "set 1, x")
    call tests%real_eq(cv%x_dot%v%v, x_dot%v%v, "set 1, x_dot")
    ! no `p` or `temp` member variables
    call tests%real_eq(cv%csa%v%v, csa%v%v, "set 1, csa")
    call tests%real_eq(cv%rm_p%v%v, rm_p%v%v, "set 1, rm_p")
    call tests%real_eq(cv%p_fs%v%v, p_fs%v%v, "set 1, p_fs")
    call tests%real_eq(cv%p_fd%v%v, p_fd%v%v, "set 1, p_fd")
    call tests%integer_eq(cv%i_cv_mirror, 2, "set 1, i_cv_mirror")
    call tests%integer_eq(cv%type, NORMAL_CV_TYPE, "set 1, type")
    call tests%character_eq(cv%label, "test_set_normal_1", "test_set_normal_1, label")
    
    call tests%integer_eq(size(cv%gas), 1, "test_set_normal_1, size(gas)")
    call tests%character_eq(cv%gas(1)%label, "dry air", "test_set_normal_1, gas label")
    
    call tests%real_eq(cv%k%v%v, k%v%v, "test_set_normal_1, k")
    call tests%real_eq(cv%delta_pre%v%v, delta_pre%v%v, "test_set_normal_1, delta_pre")
    call tests%real_eq(cv%x_stop%v%v, X_STOP_DEFAULT, "test_set_normal_1, x_stop")
    call tests%real_eq(cv%m_spring%v%v, 0.0_WP, "test_set_normal_1, m_spring")
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, temp%v%v, "test_set_normal_1, temp")
    
    vol_cv = cv%vol()
    call tests%real_eq(vol_cv%v%v, 0.05_WP, "test_set_normal_1, vol")
    
    ! Not exact, based on multiple of ambient density for simplicity.
    rho_cv = cv%rho()
    call tests%real_eq(rho_cv%v%v, sqrt(2.0_WP)*RHO_ATM, "test_set_normal_1, rho", abs_tol=1.0e-3_WP)
    
    p_cv = cv%p()
    call tests%real_eq(p_cv%v%v, p%v%v, "test_set_normal_1, p")
    
    ! Not exact, based on multiple of ambient density for simplicity.
    u = DRY_AIR%u(temp)
    call tests%real_eq(cv%m(1)%v%v, sqrt(2.0_WP)*RHO_ATM*0.05_WP, "test_set_normal_1, m(1)", abs_tol=1.0e-4_WP)
    call tests%real_eq(cv%e%v%v, u%v%v*sqrt(2.0_WP)*RHO_ATM*0.05_WP, "test_set_normal_1, e", abs_tol=10.0_WP)
end subroutine test_set_normal_1

subroutine test_set_normal_2(tests)
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
    type(si_length)       :: delta_pre, x_stop
    type(gas_type)        :: gas(2)
    
    real(WP)        :: n(2), m(2)
    type(si_volume) :: vol_cv
    
    ! > 0.18 kmol of methane (CH4) and 0.274 of butane (C4H10)
    ! methane (CH4)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/CH4/h1H4> ("Fluid Properties" at 0.101325 MPa)
    gas(1) = gas_type(label = "CH4", &
                      gamma = 2.2359_WP/1.7129_WP, &
                      u_0   = 758.87e3_WP, &
                      h_0   = 914.09e3_WP, &
                      mm    = 16.04e-3_WP, & ! moran_fundamentals_2008 table A-1
                      p_c   = 10.0_WP*46.4e5_WP) ! moran_fundamentals_2008 table A-1, set larger to avoid ideal gas validity check
    
    ! butane (C4H10)
    ! <https://webbook.nist.gov/cgi/inchi/InChI%3D1S/C4H10/c1-3-4-2/h3-4H2%2C1-2H3> ("Fluid Properties" at 0.101325 MPa)
    gas(2) = gas_type(label = "C4H10", &
                      gamma = 1.7388_WP/1.5744_WP, &
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
    call delta_pre%v%init_const(0.0_WP, 0)
    call x_stop%v%init_const(1.0_WP, 0) ! to test non-default `x_stop`
    call assert(.not. is_close(x_stop%v%v, X_STOP_DEFAULT), "x_stop must not equal the default for this test")
    
    call cv%set(x, x_dot, y, p, temp, "test_set_normal_2", csa, rm_p, p_fs, p_fd, k, delta_pre, gas, 2, x_stop=x_stop)
    
    ! All the masses are slightly off. This is expected, as total mass is not an input here.
    ! Pressure is the input setting the total mass, and it's only approximate.
    call tests%real_eq(cv%m(1)%v%v, m(1), "test_set_normal_2, cv%m(1)", abs_tol=1.0e-1_WP)
    call tests%real_eq(cv%m(2)%v%v, m(2), "test_set_normal_2, cv%m(2)", abs_tol=1.0e-1_WP)
    m_total = cv%m_total()
    call tests%real_eq(m_total%v%v, m(1) + m(2), "test_set_normal_2, cv%m_total", abs_tol=1.0e-1_WP)
    m_total = cv%vol() * cv%rho_eos(p, temp, y)
    call tests%real_eq(m_total%v%v, m(1) + m(2), "test_set_normal_2, alt m_total", abs_tol=1.0e-1_WP)
    
    ! p. 616: mole fractions, book is probably only accurate to 3 decimal points
    chi_cv = cv%chi()
    call tests%real_eq(chi_cv(1)%v%v, 0.396_WP, "test_set_normal_2, chi_cv(1)", abs_tol=1.0e-3_WP)
    call tests%real_eq(chi_cv(2)%v%v, 0.604_WP, "test_set_normal_2, chi_cv(2)", abs_tol=1.0e-3_WP)
    
    ! mass fractions
    y_cv = cv%y()
    call tests%real_eq(y_cv(1)%v%v, y(1)%v%v, "test_set_normal_2, y_cv(1)")
    call tests%real_eq(y_cv(2)%v%v, y(2)%v%v, "test_set_normal_2, y_cv(2)")
    
    vol_cv = cv%vol()
    call tests%real_eq(vol_cv%v%v, 0.241_WP, "test_set_normal_2, vol")
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, temp%v%v, "test_set_normal_2, temp")
    
    p_cv = cv%p()
    call tests%real_eq(p_cv%v%v, p%v%v, "test_set_normal_2, p")
    
    call tests%real_eq(cv%x_stop%v%v, 1.0_WP, "test_set_normal_2, x_stop")
end subroutine test_set_normal_2

subroutine test_set_normal_3(tests)
    ! Testing isentropic filling and `m_spring`.
    
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
    type(si_length)       :: delta_pre
    type(si_mass)         :: m_spring
    
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
    call delta_pre%v%init_const(3.0_WP, 0)
    call m_spring%v%init_const(2.5_WP, 0)
    
    call cv%set(x, x_dot, y, p, temp_atm, "test_set_normal_3", csa, rm_p, p_fs, p_fd, k, delta_pre, &
                        [DRY_AIR], 2, isentropic_filling=.true., p_atm=p_atm, m_spring=m_spring)
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, 300.0_WP*(2.0_WP**(0.4_WP/1.4_WP)), &
                        "test_set_normal_3, isentropic_filling=.true., cv%temp")
    call tests%real_eq(cv%m_spring%v%v, 2.5_WP, "test_set_normal_3, cv%m_spring")
end subroutine test_set_normal_3

subroutine test_set_const(tests)
    use gasdata, only: DRY_AIR
    use cva, only: cv_type, MIRROR_CV_TYPE, CONST_EOS, X_STOP_DEFAULT
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)         :: cv
    type(si_area)         :: csa
    type(si_pressure)     :: p
    type(si_temperature)  :: temp
    type(unitless)        :: y(1)
    
    call csa%v%init_const(2.0_WP, 0)
    call p%v%init_const(2.0e5_WP, 0)
    call temp%v%init_const(350.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    
    call cv%set_const("test_set_const", csa, p, temp, [DRY_AIR], y, 3)
    
    call tests%real_eq(cv%x%v%v, 1.0_WP, "test_set_const, x")
    call tests%real_eq(cv%x_dot%v%v, 0.0_WP, "test_set_const, x_dot")
    call tests%real_eq(cv%e%v%v, 0.0_WP, "test_set_const, e")
    call tests%real_eq(cv%rm_p%v%v, 0.0_WP, "test_set_const, rm_p")
    call tests%real_eq(cv%p_fs%v%v, 0.0_WP, "test_set_const, p_fs")
    call tests%real_eq(cv%p_fd%v%v, 0.0_WP, "test_set_const, p_fd")
    call tests%real_eq(cv%k%v%v, 0.0_WP, "test_set_const, k")
    call tests%real_eq(cv%delta_pre%v%v, 0.0_WP, "test_set_const, delta_pre")
    call tests%real_eq(cv%x_stop%v%v, X_STOP_DEFAULT, "test_set_const, x_stop")
    
    call tests%character_eq(cv%label, "test_set_const", "test_set_const, label")
    call tests%integer_eq(cv%type, MIRROR_CV_TYPE, "test_set_const, type")
    call tests%integer_eq(cv%eos, CONST_EOS, "test_set_const, eos")
    call tests%integer_eq(size(cv%gas), 1, "test_set_const, size(gas)")
    call tests%character_eq(cv%gas(1)%label, "dry air", "test_set_const, gas label")
    call tests%real_eq(cv%csa%v%v, csa%v%v, "test_set_const, csa")
    call tests%real_eq(cv%p_const%v%v, p%v%v, "test_set_const, p_const")
    call tests%real_eq(cv%temp_const%v%v, temp%v%v, "test_set_const, temp_const")
    call tests%integer_eq(cv%i_cv_mirror, 3, "test_set_const, i_cv_mirror")
end subroutine test_set_const

subroutine test_rates(tests)
    use gasdata, only: DRY_AIR
    use cva, only: cv_system_type, d_x_d_t, d_x_dot_d_t, d_m_k_d_t, d_e_d_t, d_e_f_d_t
    
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
    type(si_length)           :: delta_pre
    type(si_acceleration)     :: d_x_dot_d_t_
    type(si_mass_flow_rate)   :: m_dots(2, 2), d_m_1_d_t
    type(si_energy_flow_rate) :: h_dots(2, 2), d_e_d_t_, d_e_f_d_t_
    
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
    call delta_pre%v%init_const(-0.5_WP, 0)
    
    allocate(sys%cv(2))
    
    call sys%cv(1)%set_const("test_rates CV 1 (constant pressure)", csa, p_atm, temp, [DRY_AIR], y, 2)
    call sys%cv(2)%set(x, x_dot, y, p, temp, "test_rates CV 2", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 1)
    
    d_x_d_t_ = d_x_d_t(sys, 2)
    call tests%real_eq(d_x_d_t_%v%v, x_dot%v%v, "d_x_d_t")
    
    d_x_dot_d_t_ = d_x_dot_d_t(sys, 2)
    call tests%real_eq(d_x_dot_d_t_%v%v, 1.5e6_WP, "d_x_dot_d_t")
    
    call m_dots(1, 1)%v%init_const(0.0_WP, 0)
    call m_dots(1, 2)%v%init_const(2.0_WP, 0)
    call m_dots(2, 1)%v%init_const(0.0_WP, 0)
    call m_dots(2, 2)%v%init_const(0.0_WP, 0)
    
    d_m_1_d_t = d_m_k_d_t(sys, m_dots, 1, 1)
    call tests%real_eq(d_m_1_d_t%v%v, -2.0_WP, "d_m_d_t, CV 1")
    
    d_m_1_d_t = d_m_k_d_t(sys, m_dots, 1, 2)
    call tests%real_eq(d_m_1_d_t%v%v, 2.0_WP, "d_m_d_t, CV 2")
    
    call h_dots(1, 1)%v%init_const(0.0_WP, 0)
    call h_dots(1, 2)%v%init_const(0.0_WP, 0)
    call h_dots(2, 1)%v%init_const(2.0e5_WP, 0)
    call h_dots(2, 2)%v%init_const(0.0_WP, 0)
    
    d_e_d_t_ = d_e_d_t(sys%cv(2), h_dots, 1)
    call tests%real_eq(d_e_d_t_%v%v, -p%v%v*csa%v%v*x_dot%v%v + 2.0e5_WP, "d_e_d_t(1)")
    
    d_e_d_t_ = d_e_d_t(sys%cv(2), h_dots, 2)
    call tests%real_eq(d_e_d_t_%v%v, -p%v%v*csa%v%v*x_dot%v%v - 2.0e5_WP, "d_e_d_t(2)")
    
    d_e_f_d_t_ = d_e_f_d_t(sys, 2)
    call tests%real_eq(d_e_f_d_t_%v%v, 4.0_WP*10.0_WP*1.0e5_WP, "d_e_f_d_t")
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
    type(si_length)          :: delta_pre
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
    call delta_pre%v%init_const(0.5_WP, 0)
    
    call cv%set(x, x_dot, y, p, temp, "test_u_h_cv", csa, rm_p, p_fs, p_fd, k, delta_pre, [N2, H2O], 2)
    
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
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_pressure)     :: p
    type(si_temperature)  :: temp
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    type(si_time)         :: t
    
    type(si_mass_flow_rate) :: m_dot_con
    real(WP) :: d_m_dot_d_delta_p
    
    call t%v%init_const(0.0_WP, 1)
    
    call delta_p%v%init(0.0_WP, 1, 1)
    
    con%active = .true.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    call con%t_opening%v%init_const(0.0_WP, 1)
    call con%alpha_0%v%init_const(1.0_WP, 1)
    call con%alpha_dot_0%v%init_const(0.0_WP, 1)
    call con%m_dot_0%v%init_const(0.0_WP, 1)
    
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
    call delta_pre%v%init_const(0.0_WP, 1)
    
    call cv_from%set(x, x_dot, y, p, temp, "from", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p - delta_p, temp, "to", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(t, cv_from, cv_to)
    
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
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1), p_r
    type(si_pressure)     :: p_in, p_out
    type(si_temperature)  :: temp
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    type(si_time)         :: t
    
    type(si_mass_flow_rate) :: m_dot_con
    real(WP) :: m_dot_con_, d_m_dot_d_p_r
    
    call t%v%init_const(0.0_WP, 1)
    
    con%active = .true.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    call con%t_opening%v%init_const(0.0_WP, 1)
    call con%alpha_0%v%init_const(1.0_WP, 1)
    call con%alpha_dot_0%v%init_const(0.0_WP, 1)
    call con%m_dot_0%v%init_const(0.0_WP, 1)
    
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
    call delta_pre%v%init_const(0.0_WP, 1)
    
    p_out = p_in*p_r
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(t, cv_from, cv_to)
    
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
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1), p_r
    type(si_pressure)     :: p_in, p_out
    type(si_temperature)  :: temp
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    type(si_time)         :: t
    
    type(si_mass_flow_rate) :: m_dot_con
    real(WP) :: m_dot_con_
    
    call t%v%init_const(0.0_WP, 1)
    
    con%active = .true.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    call con%t_opening%v%init_const(0.0_WP, 1)
    call con%alpha_0%v%init_const(1.0_WP, 1)
    call con%alpha_dot_0%v%init_const(0.0_WP, 1)
    call con%m_dot_0%v%init_const(0.0_WP, 1)
    
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
    call delta_pre%v%init_const(0.0_WP, 1)
    
    p_out = p_in*p_r
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(t, cv_from, cv_to)
    
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
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_pressure)     :: p_in, p_out
    type(si_temperature)  :: temp
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    type(si_time)         :: t
    
    type(si_mass_flow_rate) :: m_dot_con
    
    call t%v%init_const(0.0_WP, 1)
    
    con%active = .false.
    call con%a_e%v%init_const(0.25_WP, 1)
    call con%b%v%init_const(0.5_WP, 1)
    call con%t_opening%v%init_const(0.0_WP, 1)
    call con%alpha_0%v%init_const(1.0_WP, 1)
    call con%alpha_dot_0%v%init_const(0.0_WP, 1)
    call con%m_dot_0%v%init_const(0.0_WP, 1)
    
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
    call delta_pre%v%init_const(0.0_WP, 1)
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    
    m_dot_con = con%m_dot(t, cv_from, cv_to)
    
    call tests%real_eq(m_dot_con%v%v,    0.0_WP, "m_dot, inactive, value")
    call tests%real_eq(m_dot_con%v%d(1), 0.0_WP, "m_dot, inactive, derivative")
end subroutine test_m_dot_4

subroutine test_alpha_m_dot(tests)
    use cva, only: cv_type, con_type
    use checks, only: assert
    use gasdata, only: DRY_AIR
    
    type(test_results_type), intent(in out) :: tests
    
    type(con_type) :: con
    type(si_time)  :: t
    type(unitless) :: alpha, alpha_expected
    type(cv_type)  :: cv_from, cv_to
    
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_pressure)     :: p_in, p_out
    type(si_temperature)  :: temp
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_fs, p_fd, p_atm
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    
    type(si_mass_flow_rate) :: m_dot_con
    
    con%active = .true.
    call con%a_e%v%init_const(0.0_WP, 1)
    call con%b%v%init_const(0.0_WP, 1)
    call con%t_opening%v%init_const(5.0_WP, 1)
    call con%alpha_0%v%init_const(0.2_WP, 1)
    con%alpha_dot_0 = 2.0_WP*(1.0_WP - con%alpha_0) ! to set the cubic term to zero
    call con%m_dot_0%v%init_const(2.0_WP, 1)
    
    call x%v%init_const(1.0_WP, 1)
    call x_dot%v%init_const(0.0_WP, 1)
    call y(1)%v%init_const(1.0_WP, 1)
    call p_in%v%init_const(20.0e5_WP, 1)
    call p_out%v%init_const(2.0e5_WP, 1)
    call temp%v%init_const(400.0_WP, 1)
    call csa%v%init_const(1.0_WP, 1)
    call rm_p%v%init_const(0.0_WP, 1)
    call p_fs%v%init_const(0.0_WP, 1)
    call p_fd%v%init_const(0.0_WP, 1)
    call p_atm%v%init_const(1.0e5_WP, 1)
    call k%v%init_const(0.0_WP, 1)
    call delta_pre%v%init_const(0.0_WP, 1)
    
    call cv_from%set(x, x_dot, y, p_in, temp, "from", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    call cv_to%set(x, x_dot, y, p_out, temp, "to", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 2)
    
    call assert(con%alpha_dot_0%v%v > 0.0_WP, "test_cva (test_alpha_m_dot): alpha_dot_0 > 0 violated")
    
    call t%v%init(0.0_WP, 1, 1)
    alpha     = con%alpha(t)
    m_dot_con = con%m_dot(t, cv_from, cv_to)
    call tests%real_eq(alpha%v%v,     con%alpha_0%v%v,                       "alpha_m_dot, t = 0.0, value")
    call tests%real_eq(alpha%v%d(1),  con%alpha_dot_0%v%v/con%t_opening%v%v, "alpha_m_dot, t = 0.0, derivative")
    call tests%real_eq(m_dot_con%v%v, con%alpha_0%v%v*con%m_dot_0%v%v,       "alpha_m_dot, m_dot, t = 0.0")
    
    t = con%t_opening/2.0_WP
    alpha          = con%alpha(t)
    alpha_expected = con%alpha_0 + con%alpha_dot_0*(t/con%t_opening) &
                        + (3.0_WP - 3.0_WP*con%alpha_0 - 2.0_WP*con%alpha_dot_0)*square(t/con%t_opening)
    m_dot_con      = con%m_dot(t, cv_from, cv_to)
    call tests%real_eq(alpha%v%v,    alpha_expected%v%v,    "alpha_m_dot, t = t_opening/2, value")
    call tests%real_eq(alpha%v%d(1), alpha_expected%v%d(1), "alpha_m_dot, t = t_opening/2, derivative")
    call tests%real_gt(alpha%v%v,    0.0_WP, "alpha_m_dot, t = t_opening/2, value, lower bound")
    call tests%real_lt(alpha%v%v,    1.0_WP, "alpha_m_dot, t = t_opening/2, value, upper bound")
    call tests%real_eq(alpha%v%d(1), 0.0_WP, "alpha_m_dot, t = t_opening/2, derivative, lower bound")
    call tests%real_eq(m_dot_con%v%v, alpha_expected%v%v*con%m_dot_0%v%v, "alpha_m_dot, m_dot, t = t_opening/2")
    
    t         = con%t_opening
    alpha     = con%alpha(t)
    m_dot_con = con%m_dot(t, cv_from, cv_to)
    call tests%real_eq(alpha%v%v,    1.0_WP, "alpha_m_dot, t = t_opening, value")
    call tests%real_eq(alpha%v%d(1), 0.0_WP, "alpha_m_dot, t = t_opening, derivative")
    call tests%real_eq(m_dot_con%v%v, con%m_dot_0%v%v, "alpha_m_dot, m_dot, t = t_opening")
end subroutine test_alpha_m_dot

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
    type(si_length)          :: delta_pre
    real(WP)                 :: m_dot_12
    type(si_specific_energy) :: h_1
    type(si_time)            :: t
    
    type(si_mass_flow_rate), allocatable   :: m_dot(:, :)
    type(si_energy_flow_rate), allocatable :: h_dot(:, :)
    
    call t%v%init_const(0.0_WP, 0)
    
    allocate(sys%cv(2))
    allocate(sys%con(2, 2))
    
    sys%con(1, 1)%active = .false.
    sys%con(2, 1)%active = .false.
    sys%con(2, 2)%active = .false.
    
    sys%con(1, 2)%active = .true.
    call sys%con(1, 2)%a_e%v%init_const(0.25_WP, 0)
    call sys%con(1, 2)%b%v%init_const(0.5_WP, 0)
    call sys%con(1, 2)%t_opening%v%init_const(0.0_WP, 0)
    call sys%con(1, 2)%alpha_0%v%init_const(1.0_WP, 0)
    call sys%con(1, 2)%alpha_dot_0%v%init_const(0.0_WP, 0)
    call sys%con(1, 2)%m_dot_0%v%init_const(0.0_WP, 0)
    
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
    call delta_pre%v%init_const(0.0_WP, 0)
    
    p_out = p_in*p_r
    
    call sys%cv(1)%set(x, x_dot, y, p_in, temp, "1", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 3)
    call sys%cv(2)%set(x, x_dot, y, p_out, temp, "2", csa, rm_p, p_fs, p_fd, k, delta_pre, [DRY_AIR], 3)
    
    m_dot_12 = sys%con(1, 2)%a_e%v%v * p_in%v%v * sqrt((1.0_WP - sys%con(1, 2)%b%v%v)/((R_BAR/DRY_AIR%mm) * temp%v%v))
    
    call sys%calculate_flows(t, m_dot, h_dot)
    
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

subroutine test_conservation_1(tests)
    use convert
    use gasdata, only: DRY_AIR
    use prec, only: PI
    use cva, only: MIRROR_CV_TYPE, SUCCESS_RUN_RC, cv_system_type, run_config_type, run_status_type, &
                    run
    
    type(test_results_type), intent(in out) :: tests
    
    type(run_config_type)             :: config
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    integer               :: n_d
    type(si_length)       :: d_e, x_3, x_4, d_3, d_4, x_stop_4
    type(si_velocity)     :: x_dot
    type(unitless)        :: y(1)
    type(si_pressure)     :: p_atm, p_3, p_4, p_fs_3, p_fd_3, p_fs_4, p_fd_4
    type(si_temperature)  :: temp_atm
    type(si_area)         :: csa_3, csa_4
    type(si_mass)         :: m_p_3, m_p_4, m_start, m_end
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    type(si_energy)       :: e_start, e_end, e_s_3_start, e_p_3_start, e_start_3, &
                                e_chamber_atm, e_barrel_atm
    
    n_d = 1
    
    allocate(sys_start)
    allocate(sys_start%cv(4))
    allocate(sys_start%con(4, 4))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(1, 2)%active = .false.
    sys_start%con(1, 3)%active = .false.
    sys_start%con(1, 4)%active = .false.
    
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    sys_start%con(2, 3)%active = .false.
    sys_start%con(2, 4)%active = .false.
    
    sys_start%con(3, 1)%active = .false.
    sys_start%con(3, 2)%active = .false.
    sys_start%con(3, 3)%active = .false.
    sys_start%con(3, 4)%active = .true.
    d_e = inch_const(0.1_WP, n_d)
    sys_start%con(3, 4)%a_e = (PI/4.0_WP)*square(d_e)
    call sys_start%con(3, 4)%b%v%init_const(0.5_WP, n_d)
    call sys_start%con(3, 4)%t_opening%v%init_const(0.0_WP, n_d)
    call sys_start%con(3, 4)%alpha_0%v%init_const(1.0_WP, n_d)
    call sys_start%con(3, 4)%alpha_dot_0%v%init_const(0.0_WP, n_d)
    call sys_start%con(3, 4)%m_dot_0%v%init_const(0.0_WP, n_d)
    
    sys_start%con(4, 1)%active = .false.
    sys_start%con(4, 2)%active = .false.
    sys_start%con(4, 3) = sys_start%con(3, 4)
    sys_start%con(4, 4)%active = .false.
    
    ! The same for every control volume.
    call y(1)%v%init_const(1.0_WP, n_d)
    call temp_atm%v%init_const(300.0_WP, n_d)
    
    ! 1: atmosphere for chamber
    call p_atm%v%init_const(1.0e5_WP, n_d)
    call d_3%v%init_const(2.0e-2_WP, n_d)
    csa_3 = (PI/4.0_WP)*square(d_3)
    call x_dot%v%init_const(-2.0_WP, n_d)
    
    call sys_start%cv(1)%set_const("atmosphere for chamber", csa_3, p_atm, temp_atm, [DRY_AIR], y, 3, x_dot=-x_dot)
    
    call tests%integer_eq(sys_start%cv(1)%type, MIRROR_CV_TYPE, "test_conservation, sys_start%cv(1)%type")
    call tests%integer_eq(sys_start%cv(1)%i_cv_mirror, 3, "test_conservation, sys_start%cv(1)%i_cv_mirror")
    
    ! 2: atmosphere for barrel
    call p_atm%v%init_const(1.0e5_WP, n_d)
    call d_4%v%init_const(2.0e-2_WP, n_d)
    csa_4 = (PI/4.0_WP)*square(d_4)
    
    call sys_start%cv(2)%set_const("atmosphere for barrel", csa_4, p_atm, temp_atm, [DRY_AIR], y, 4)
    
    call tests%integer_eq(sys_start%cv(2)%type, MIRROR_CV_TYPE, "test_conservation, sys_start%cv(2)%type")
    call tests%integer_eq(sys_start%cv(2)%i_cv_mirror, 4, "test_conservation, sys_start%cv(2)%i_cv_mirror")
    
    ! 3: chamber
    call x_3%v%init_const(10.0e-2_WP, n_d)
    call p_3%v%init(5.0e5_WP, 1, n_d) ! This will also test the derivatives a bit.
    call m_p_3%v%init_const(30.0e-3_WP, n_d)
    call p_fs_3%v%init_const(0.2e5_WP, n_d)
    call p_fd_3%v%init_const(0.1e5_WP, n_d)
    call k%v%init_const(700.0_WP, n_d)
    call delta_pre%v%init_const(-1.0e-2_WP, n_d)
    
    call sys_start%cv(3)%set(x_3, x_dot, y, p_3, temp_atm, "pressure chamber", csa_3, 1.0_WP/m_p_3, p_fs_3, p_fd_3, k, &
                                    delta_pre, [DRY_AIR], 1)
    ! `isentropic_filling=.true.` requires that `p_atm > 0`, so it's not used here.
    
    ! 4: barrel
    
    call x_dot%v%init_const(0.0_WP, n_d)
    call x_4%v%init_const(10.0e-2_WP, n_d)
    call p_4%v%init_const(1.0e5_WP, n_d)
    call m_p_4%v%init_const(1.0e-3_WP, n_d)
    call p_fs_4%v%init_const(0.2e5_WP, n_d)
    call p_fd_4%v%init_const(0.1e5_WP, n_d)
    x_stop_4 = x_4 + inch_const(12.0_WP, n_d)
    call k%v%init_const(0.0_WP, n_d)
    call delta_pre%v%init_const(0.0_WP, n_d)
    
    call sys_start%cv(4)%set(x_4, x_dot, y, p_4, temp_atm, "barrel", csa_4, 1.0_WP/m_p_4, p_fs_4, p_fd_4, k, &
                                delta_pre, [DRY_AIR], 2, x_stop=x_stop_4)
    
    call config%set("test_conservation", 1, csv_output=.true., csv_frequency=100)
    call run(config, sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, SUCCESS_RUN_RC, "test_conservation, status%rc")
    
    if (status%rc /= SUCCESS_RUN_RC) then
        print *, status%data(1), status%i_cv(1)
    end if
    
    e_s_3_start = sys_start%cv(3)%e_s()
    call tests%real_eq(e_s_3_start%v%v, 0.5_WP*(700.0_WP)*(9.0e-2_WP)**2, "test_conservation, e_s_3_start")
    
    e_p_3_start = sys_start%cv(3)%e_p()
    call tests%real_eq(e_p_3_start%v%v, 0.5_WP*(30.0e-3_WP)*(2.0_WP)**2, "test_conservation, e_p_3_start")
    
    e_start_3 = sys_start%cv(3)%e_total()
    call tests%real_eq(e_start_3%v%v, sys_start%cv(3)%e%v%v + e_s_3_start%v%v + e_p_3_start%v%v, &
                            "test_conservation, sys_start%cv(3)%e_total()")
    
    m_start = sys_start%m_total()
    m_end   = sys_end%m_total()
    call tests%real_eq(m_start%v%v, m_end%v%v, "test_conservation, m_start == m_end")
    call tests%real_gt(m_start%v%d(1), 0.0_WP, "test_conservation, m_start derivative is positive")
    call tests%real_eq(m_start%v%d(1) - m_end%v%d(1), 0.0_WP, "test_conservation, delta m derivative is zero", abs_tol=1.0e-22_WP)
    
    e_start = sys_start%e_total()
    e_end   = sys_end%e_total()
    call tests%real_eq(e_start%v%v, e_end%v%v, "test_conservation, e_start == e_end", abs_tol=1.0e-3_WP)
    call tests%real_gt(e_start%v%d(1), 0.0_WP, "test_conservation, e_start derivative is positive")
    call tests%real_eq(e_start%v%d(1) - e_end%v%d(1), 0.0_WP, "test_conservation, delta e derivative is zero", abs_tol=1.0e-5_WP)
    
    e_chamber_atm = sys_end%cv(1)%e_total()
    call tests%character_eq(sys_end%cv(1)%label, "atmosphere for chamber", "test_conservation, chamber atmosphere label")
    call tests%real_lt(e_chamber_atm%v%v, 0.0_WP, "test_conservation, chamber atmosphere energy is non-zero")
    call tests%real_gt(abs(sys_end%cv(1)%x_dot%v%v), 0.0_WP, &
                        "test_conservation, chamber atmosphere plunger velocity is non-zero")
    call tests%real_eq(sys_end%cv(1)%x_dot%v%v, -sys_end%cv(3)%x_dot%v%v, &
                        "test_conservation, chamber atmosphere velocity matches plunger velocity")
    call tests%real_eq(e_chamber_atm%v%v, sys_end%cv(1)%e%v%v, &
                            "test_conservation, chamber atmosphere energy is internal energy")
    
    e_barrel_atm = sys_end%cv(2)%e_total()
    call tests%character_eq(sys_end%cv(2)%label, "atmosphere for barrel", "test_conservation, barrel atmosphere label")
    call tests%real_gt(e_barrel_atm%v%v, 0.0_WP, "test_conservation, barrel atmosphere energy is non-zero")
    call tests%real_gt(abs(sys_end%cv(2)%x_dot%v%v), 0.0_WP, &
                        "test_conservation, barrel atmosphere projectile/plunger velocity is non-zero")
    call tests%real_eq(sys_end%cv(2)%x_dot%v%v, -sys_end%cv(4)%x_dot%v%v, &
                        "test_conservation, barrel atmosphere velocity matches projectile velocity")
    call tests%real_eq(e_barrel_atm%v%v, sys_end%cv(2)%e%v%v, &
                            "test_conservation, barrel atmosphere energy is internal energy")
    
    call tests%real_gt(sys_end%cv(3)%e_f%v%v, 0.0_WP, "test_conservation, chamber friction loss is non-zero")
    call tests%real_gt(sys_end%cv(4)%e_f%v%v, 0.0_WP, "test_conservation, barrel friction loss is non-zero")
end subroutine test_conservation_1

subroutine test_conservation_2(tests)
    ! Testing gas species mass conservation.
    
    use gasdata, only: DRY_AIR, H2O
    use cva, only: TIMEOUT_RUN_RC, cv_system_type, run_config_type, run_status_type, run
    
    type(test_results_type), intent(in out) :: tests
    
    type(run_config_type)             :: config
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    integer               :: n_d
    type(si_length)       :: x
    type(si_velocity)     :: x_dot
    type(unitless)        :: one, y_1(2), y_2(2)
    type(si_pressure)     :: p_1, p_2, p_fs, p_fd
    type(si_temperature)  :: temp
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_stiffness)    :: k
    type(si_length)       :: delta_pre
    type(si_time)         :: t_stop
    type(si_mass)         :: m_dry_air_start, m_dry_air_end, m_h2o_start, m_h2o_end
    
    n_d = 2 ! Look at derivatives of mass fraction.
    
    call one%v%init_const(1.0_WP, n_d)
    call y_1(1)%v%init(1.0_WP, 1, n_d) ! CV 1 is `DRY_AIR`
    y_1(2) = one - y_1(1)
    call y_2(2)%v%init(1.0_WP, 2, n_d) ! CV 2 is `H2O`
    y_2(1) = one - y_2(2)
    
    call x%v%init_const(10.0e-2_WP, n_d)
    call x_dot%v%init_const(0.0_WP, n_d)
    call p_1%v%init_const(5.0e5_WP, n_d)
    call p_2%v%init_const(1.0e5_WP, n_d)
    call temp%v%init_const(300.0_WP, n_d)
    call csa%v%init_const(0.0025_WP, n_d)
    call rm_p%v%init_const(0.0_WP, n_d)
    call p_fs%v%init_const(0.0_WP, n_d)
    call p_fd%v%init_const(0.0_WP, n_d)
    call k%v%init_const(0.0_WP, n_d)
    call delta_pre%v%init_const(0.0_WP, n_d)
    call t_stop%v%init_const(0.1_WP, n_d)
    
    allocate(sys_start)
    allocate(sys_start%cv(2))
    allocate(sys_start%con(2, 2))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(1, 2)%active = .true.
    sys_start%con(1, 2)%a_e = 0.25_WP*csa
    call sys_start%con(1, 2)%b%v%init_const(0.5_WP, n_d)
    call sys_start%con(1, 2)%t_opening%v%init_const(0.0_WP, n_d)
    call sys_start%con(1, 2)%alpha_0%v%init_const(1.0_WP, n_d)
    call sys_start%con(1, 2)%alpha_dot_0%v%init_const(0.0_WP, n_d)
    call sys_start%con(1, 2)%m_dot_0%v%init_const(0.0_WP, n_d)
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    
    call sys_start%cv(1)%set(x, x_dot, y_1, p_1, temp, "chamber 1", csa, rm_p, p_fs, p_fd, k, &
                                    delta_pre, [DRY_AIR, H2O], 0)
    call sys_start%cv(2)%set(x, x_dot, y_2, p_2, temp, "chamber 2", csa, rm_p, p_fs, p_fd, k, &
                                    delta_pre, [DRY_AIR, H2O], 0)
    
    call config%set("test_conservation_2", n_d, t_stop=t_stop)
    call run(config, sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, TIMEOUT_RUN_RC, "test_conservation_2, status%rc")
    
    m_dry_air_start = sys_start%cv(1)%m(1) + sys_start%cv(2)%m(1)
    m_dry_air_end   = sys_end%cv(1)%m(1)   + sys_end%cv(2)%m(1)
    m_h2o_start     = sys_start%cv(1)%m(2) + sys_start%cv(2)%m(2)
    m_h2o_end       = sys_end%cv(1)%m(2)   + sys_end%cv(2)%m(2)
    
    call tests%real_eq(sys_start%cv(1)%m(2)%v%v, 0.0_WP, "test_conservation_2, CV 1 at start has no H2O")
    call tests%real_eq(sys_start%cv(2)%m(1)%v%v, 0.0_WP, "test_conservation_2, CV 2 at start has no DRY_AIR")
    
    call tests%real_eq(m_dry_air_start%v%v, m_dry_air_end%v%v, "test_conservation_2, DRY_AIR mass is conserved", &
                        abs_tol=1.0e-14_WP)
    call tests%real_eq(m_h2o_start%v%v, m_h2o_end%v%v, "test_conservation_2, H2O mass is conserved", &
                        abs_tol=1.0e-14_WP)
end subroutine test_conservation_2

subroutine test_mirror_1(tests)
    ! equilibrium test where velocities will be small
    
    use convert
    use gasdata, only: DRY_AIR
    use cva, only: MIRROR_CV_TYPE, TIMEOUT_RUN_RC, cv_system_type, run_config_type, run_status_type, &
                    run
    
    type(test_results_type), intent(in out) :: tests
    
    type(run_config_type)             :: config
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    type(si_length)      :: x_1, x_2, x_1end, x_2end, x_start_sum, x_end_sum
    type(si_velocity)    :: x_dot
    type(unitless)       :: y(1)
    type(si_pressure)    :: p_1, p_2, p_fs, p_fd, p_1end, p_2end
    type(si_temperature) :: temp
    type(si_area)        :: csa
    type(si_mass)        :: m_p
    type(si_stiffness)   :: k
    type(si_length)      :: delta_pre
    type(si_time)        :: t_stop
    
    allocate(sys_start)
    allocate(sys_start%cv(2))
    allocate(sys_start%con(2, 2))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(1, 2)%active = .false.
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    
    call x_1%v%init_const(10.0e-2_WP, 0)
    call x_2%v%init_const(10.0e-2_WP, 0)
    call x_dot%v%init_const(0.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p_1%v%init_const(5.0e5_WP, 0)
    call p_2%v%init_const(1.0e5_WP, 0)
    call temp%v%init_const(300.0_WP, 0)
    call csa%v%init_const(0.0025_WP, 0)
    call m_p%v%init_const(1.0e-3_WP, 0)
    call p_fs%v%init_const(0.0_WP, 0)
    call p_fd%v%init_const(0.1e5_WP, 0)
    call k%v%init_const(0.0_WP, 0)
    call delta_pre%v%init_const(0.0_WP, 0)
    call t_stop%v%init_const(0.1_WP, 0)
    
    call sys_start%cv(1)%set(x_1, x_dot, y, p_1, temp, "chamber 1", csa, 1.0_WP/m_p, p_fs, p_fd, k, &
                                    delta_pre, [DRY_AIR], 2)
    call sys_start%cv(2)%set(x_2, x_dot, y, p_2, temp, "chamber 2", csa, 1.0_WP/m_p, p_fs, p_fd, k, &
                                    delta_pre, [DRY_AIR], 1, type=MIRROR_CV_TYPE)
    
    call config%set("test_mirror_1", 0, t_stop=t_stop)
    call run(config, sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, TIMEOUT_RUN_RC, "test_mirror_1, status%rc")
    
    ! The pressures don't need to equilibrate for the plunger motion test.
    ! But checking that they do is a decent test, so I will.
    p_1end = sys_end%cv(1)%p()
    p_2end = sys_end%cv(2)%p()
    call tests%real_eq(p_1end%v%v, p_2end%v%v, "test_mirror_1, p_1 == p_2 at end", abs_tol=1.0_WP)
    
    x_1end = sys_end%cv(1)%x
    x_2end = sys_end%cv(2)%x
    x_start_sum = x_1    + x_2
    x_end_sum   = x_1end + x_2end
    call tests%real_eq(x_start_sum%v%v, x_end_sum%v%v, "test_mirror_1, plunger motion is as expected")
    
    call tests%real_eq(sys_end%cv(1)%x_dot%v%v, 0.0_WP, "test_mirror_1, CV 1 x_dot", abs_tol=1.0e-4_WP)
    call tests%real_eq(sys_end%cv(2)%x_dot%v%v, 0.0_WP, "test_mirror_1, CV 2 x_dot", abs_tol=1.0e-4_WP)
    call tests%real_eq(sys_end%cv(1)%x_dot%v%v, -sys_end%cv(2)%x_dot%v%v, &
                            "test_mirror_1, CV 1 x_dot = - CV 2 x_dot")
end subroutine test_mirror_1

subroutine test_mirror_2(tests)
    ! short time test where velocities won't be small
    
    use convert
    use gasdata, only: DRY_AIR
    use cva, only: MIRROR_CV_TYPE, TIMEOUT_RUN_RC, cv_system_type, run_config_type, run_status_type, &
                    run
    
    type(test_results_type), intent(in out) :: tests
    
    type(run_config_type)             :: config
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    type(si_length)      :: x_1, x_2, x_1end, x_2end, x_start_sum, x_end_sum
    type(si_velocity)    :: x_dot
    type(unitless)       :: y(1)
    type(si_pressure)    :: p_1, p_2, p_fs, p_fd
    type(si_temperature) :: temp
    type(si_area)        :: csa
    type(si_mass)        :: m_p
    type(si_stiffness)   :: k
    type(si_length)      :: delta_pre
    type(si_time)        :: t_stop
    
    allocate(sys_start)
    allocate(sys_start%cv(2))
    allocate(sys_start%con(2, 2))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(1, 2)%active = .false.
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    
    call x_1%v%init_const(10.0e-2_WP, 0)
    call x_2%v%init_const(10.0e-2_WP, 0)
    call x_dot%v%init_const(0.0_WP, 0)
    call y(1)%v%init_const(1.0_WP, 0)
    call p_1%v%init_const(5.0e5_WP, 0)
    call p_2%v%init_const(1.0e5_WP, 0)
    call temp%v%init_const(300.0_WP, 0)
    call csa%v%init_const(0.0025_WP, 0)
    call m_p%v%init_const(1.0e-3_WP, 0)
    call p_fs%v%init_const(0.0_WP, 0)
    call p_fd%v%init_const(0.1e5_WP, 0)
    call k%v%init_const(0.0_WP, 0)
    call delta_pre%v%init_const(0.0_WP, 0)
    call t_stop%v%init_const(1.0e-3_WP, 0)
    
    call sys_start%cv(1)%set(x_1, x_dot, y, p_1, temp, "chamber 1", csa, 1.0_WP/m_p, p_fs, p_fd, k, &
                                    delta_pre, [DRY_AIR], 2)
    call sys_start%cv(2)%set(x_2, x_dot, y, p_2, temp, "chamber 2", csa, 1.0_WP/m_p, p_fs, p_fd, k, &
                                    delta_pre, [DRY_AIR], 1, type=MIRROR_CV_TYPE)
    
    call config%set("test_mirror_2", 0, t_stop=t_stop)
    call run(config, sys_start, sys_end, status)
    
    call tests%integer_eq(status%rc, TIMEOUT_RUN_RC, "test_mirror_2, status%rc")
    
    x_1end = sys_end%cv(1)%x
    x_2end = sys_end%cv(2)%x
    x_start_sum = x_1    + x_2
    x_end_sum   = x_1end + x_2end
    call tests%real_eq(x_start_sum%v%v, x_end_sum%v%v, "test_mirror_2, plunger motion is as expected")
    call tests%real_eq(sys_end%cv(1)%x_dot%v%v, -sys_end%cv(2)%x_dot%v%v, &
                            "test_mirror_2, CV 1 x_dot = - CV 2 x_dot")
end subroutine test_mirror_2

subroutine test_check_sys(tests)
    use gasdata, only: DRY_AIR, H2O
    use cva, only: IDEAL_EOS, NORMAL_CV_TYPE, CONTINUE_RUN_RC, SUCCESS_RUN_RC, TIMEOUT_RUN_RC, NEGATIVE_CV_M_TOTAL_RUN_RC, &
                    NEGATIVE_CV_TEMP_RUN_RC, MASS_TOLERANCE_RUN_RC, ENERGY_TOLERANCE_RUN_RC, &
                    MASS_DERIV_TOLERANCE_RUN_RC, ENERGY_DERIV_TOLERANCE_RUN_RC, IDEAL_EOS_RUN_RC, &
                    MIRROR_X_TOLERANCE_RUN_RC, MIRROR_CV_TYPE, &
                    !X_BLOW_UP_RUN_RC, M_BLOW_UP_RUN_RC, E_BLOW_UP_RUN_RC, E_F_BLOW_UP_RUN_RC, X_DOT_BLOW_UP_RUN_RC, &
                    run_config_type, cv_system_type, run_status_type, check_sys
    
    type(test_results_type), intent(in out) :: tests
    
    integer                           :: n_d
    type(run_config_type)             :: config
    type(cv_system_type), allocatable :: sys, sys_start
    type(si_time)                     :: t_stop, t
    type(si_pressure)                 :: p
    type(run_status_type)             :: status
    type(si_temperature)              :: temp
    
    n_d = 2
    
    allocate(sys)
    allocate(sys%cv(2))
    ! `sys%con` doesn't need to be allocated for this test.
    allocate(sys%cv(1)%m(2))
    allocate(sys%cv(2)%m(2))
    
    call t_stop%v%init_const(1.0_WP, n_d)
    call config%set("test_check_sys", n_d, t_stop=t_stop)
    
    ! `m_start` is 1.0
    ! `e_start` is 2.0
    call t%v%init_const(0.01_WP, n_d)
    
    ! `CONTINUE_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    sys_start = sys
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, CONTINUE_RUN_RC, "test_check_sys, CONTINUE_RUN_RC, status%rc")
    
    ! `SUCCESS_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(4.0_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, SUCCESS_RUN_RC, "test_check_sys, SUCCESS_RUN_RC, status%rc")
    call tests%integer_eq(size(status%i_cv), 1, "test_check_sys, SUCCESS_RUN_RC, size(status%i_cv)")
    call tests%integer_eq(status%i_cv(1), 1, "test_check_sys, SUCCESS_RUN_RC, status%i_cv(1)")
    
    ! `TIMEOUT_RUN_RC`
    
    call t_stop%v%init_const(0.0_WP, n_d)
    call config%set("test_check_sys", n_d, t_stop=t_stop)
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, TIMEOUT_RUN_RC, "test_check_sys, TIMEOUT_RUN_RC, status%rc")
    
    ! `NEGATIVE_CV_M_TOTAL_RUN_RC`
    
    call t_stop%v%init_const(1.0_WP, n_d)
    call config%set("test_check_sys", n_d, t_stop=t_stop)
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(-0.5_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(1.5_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, NEGATIVE_CV_M_TOTAL_RUN_RC, "test_check_sys, NEGATIVE_CV_M_TOTAL_RUN_RC, status%rc")
    call tests%integer_eq(size(status%i_cv), 1, "test_check_sys, NEGATIVE_CV_M_TOTAL_RUN_RC, size(status%i_cv)")
    call tests%integer_eq(status%i_cv(1), 1, "test_check_sys, NEGATIVE_CV_M_TOTAL_RUN_RC, status%i_cv(1)")
    call tests%integer_eq(size(status%data), 2, "test_check_sys, NEGATIVE_CV_M_TOTAL_RUN_RC, size(status%data)")
    call tests%real_eq(status%data(1), -0.5_WP, "test_check_sys, NEGATIVE_CV_M_TOTAL_RUN_RC, status%data(1)")
    call tests%real_eq(status%data(2), 1.5_WP, "test_check_sys, NEGATIVE_CV_M_TOTAL_RUN_RC, status%data(2)")
    
    ! `NEGATIVE_CV_TEMP_RUN_RC`
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.5_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.0_WP, n_d)
    call temp%v%init_const(-100.0_WP, n_d)
    sys%cv(2)%e = sys%cv(2)%m(1)*DRY_AIR%u(temp)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    sys%cv(1)%e = sys_start%e_total() - sys%cv(2)%e
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, NEGATIVE_CV_TEMP_RUN_RC, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, status%rc")
    call tests%integer_eq(size(status%i_cv), 1, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, size(status%i_cv)")
    call tests%integer_eq(status%i_cv(1), 2, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, status%i_cv(1)")
    call tests%integer_eq(size(status%data), 2, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, size(status%data)")
    call tests%real_gt(status%data(1), 0.0_WP, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, status%data(1) sign")
    call tests%real_lt(status%data(2), 0.0_WP, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, status%data(2) sign")
    call tests%real_eq(status%data(2), -100.0_WP, "test_check_sys, NEGATIVE_CV_TEMP_RUN_RC, status%data(2)")
    
    ! `MASS_TOLERANCE_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.5_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, MASS_TOLERANCE_RUN_RC, "test_check_sys, MASS_TOLERANCE_RUN_RC, status%rc")
    call tests%integer_eq(size(status%data), 1, "test_check_sys, MASS_TOLERANCE_RUN_RC, size(status%data)")
    call tests%real_eq(status%data(1), 0.25_WP, "test_check_sys, MASS_TOLERANCE_RUN_RC, status%data(1)")
    
    ! `ENERGY_TOLERANCE_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, ENERGY_TOLERANCE_RUN_RC, "test_check_sys, ENERGY_TOLERANCE_RUN_RC, status%rc")
    call tests%integer_eq(size(status%data), 1, "test_check_sys, ENERGY_TOLERANCE_RUN_RC, size(status%data)")
    call tests%real_eq(status%data(1), 0.5_WP, "test_check_sys, ENERGY_TOLERANCE_RUN_RC, status%data(1)")
    
    ! `MASS_DERIV_TOLERANCE_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init(0.25_WP, 1, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call tests%integer_eq(size(sys%cv(1)%x%v%d), 2, "test_check_sys, MASS_DERIV_TOLERANCE_RUN_RC, n_d")
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, MASS_DERIV_TOLERANCE_RUN_RC, "test_check_sys, MASS_DERIV_TOLERANCE_RUN_RC, status%rc")
    call tests%integer_eq(size(status%data), 2, "test_check_sys, MASS_DERIV_TOLERANCE_RUN_RC, size(status%data)")
    call tests%real_eq(status%data(1), 1.0_WP, "test_check_sys, MASS_DERIV_TOLERANCE_RUN_RC, status%data(1)")
    call tests%real_eq(status%data(2), 1.0_WP, "test_check_sys, MASS_DERIV_TOLERANCE_RUN_RC, status%data(2)")
    
    ! `ENERGY_DERIV_TOLERANCE_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(1)%e%v%init(1.0_WP, 1, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.25_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call tests%integer_eq(size(sys%cv(1)%x%v%d), 2, "test_check_sys, MASS_DERIV_TOLERANCE_RUN_RC, n_d")
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, ENERGY_DERIV_TOLERANCE_RUN_RC, "test_check_sys, ENERGY_DERIV_TOLERANCE_RUN_RC, status%rc")
    call tests%integer_eq(size(status%data), 2, "test_check_sys, ENERGY_DERIV_TOLERANCE_RUN_RC, size(status%data)")
    call tests%real_eq(status%data(1), 1.0_WP, "test_check_sys, ENERGY_DERIV_TOLERANCE_RUN_RC, status%data(1)")
    call tests%real_eq(status%data(2), 1.0_WP, "test_check_sys, ENERGY_DERIV_TOLERANCE_RUN_RC, status%data(2)")
    
    ! `IDEAL_EOS_RUN_RC`
    
    call sys%cv(1)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(1)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%m(1)%v%init_const(0.5_WP, n_d)
    call sys%cv(1)%m(2)%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(1)%label = "CV1"
    sys%cv(1)%eos   = IDEAL_EOS
    sys%cv(1)%type  = NORMAL_CV_TYPE
    sys%cv(1)%gas   = [DRY_AIR, H2O]
    call sys%cv(1)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(1)%csa%v%init_const(1.0e-4_WP, n_d)
    call sys%cv(1)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(1)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(1)%k%v%init_const(10.0_WP, n_d)
    sys%cv(1)%delta_pre = -sys%cv(1)%x
    sys%cv(1)%i_cv_mirror = 0
    
    call sys%cv(2)%x%v%init_const(0.1_WP, n_d)
    call sys%cv(2)%x_dot%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%m(1)%v%init_const(0.5_WP, n_d)
    call sys%cv(2)%m(2)%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%e%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%e_f%v%init_const(0.0_WP, n_d)
    sys%cv(2)%label = "CV1"
    sys%cv(2)%eos   = IDEAL_EOS
    sys%cv(2)%type  = NORMAL_CV_TYPE
    sys%cv(2)%gas   = [DRY_AIR, H2O]
    call sys%cv(2)%x_stop%v%init_const(2.0_WP, n_d)
    call sys%cv(2)%csa%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%rm_p%v%init_const(1.0_WP, n_d)
    call sys%cv(2)%m_spring%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fs%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%p_fd%v%init_const(0.0_WP, n_d)
    call sys%cv(2)%k%v%init_const(10.0_WP, n_d)
    sys%cv(2)%delta_pre = -sys%cv(2)%x
    sys%cv(2)%i_cv_mirror = 0
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, IDEAL_EOS_RUN_RC, "test_check_sys, IDEAL_EOS_RUN_RC, status%rc")
    call tests%integer_eq(size(status%i_cv), 1, "test_check_sys, IDEAL_EOS_RUN_RC, size(status%i_cv)")
    call tests%integer_eq(status%i_cv(1), 1, "test_check_sys, IDEAL_EOS_RUN_RC, status%i_cv(1)")
    call tests%integer_eq(size(status%data), 2, "test_check_sys, IDEAL_EOS_RUN_RC, size(status%data)")
    p = sys%cv(1)%p()
    call tests%real_gt(p%v%v, DRY_AIR%p_c, "test_check_sys, IDEAL_EOS_RUN_RC, p >= p_c")
    call tests%real_gt(status%data(1), DRY_AIR%p_c, "test_check_sys, IDEAL_EOS_RUN_RC, status%data(1) sign")
    call tests%real_lt(status%data(2), DRY_AIR%p_c, "test_check_sys, IDEAL_EOS_RUN_RC, status%data(2) sign")
    
    ! `MIRROR_X_TOLERANCE_RUN_RC`
    
    sys = sys_start
    sys%cv(1)%i_cv_mirror = 2
    sys%cv(2)%i_cv_mirror = 1
    sys%cv(2)%type        = MIRROR_CV_TYPE
    call sys%cv(2)%x%v%init_const(0.2_WP, n_d)
    
    call check_sys(config, sys, sys_start, t, status)
    call tests%integer_eq(status%rc, MIRROR_X_TOLERANCE_RUN_RC, "test_check_sys, MIRROR_X_TOLERANCE_RUN_RC, status%rc")
    
    ! TODO: `X_BLOW_UP_RUN_RC`
    ! TODO: `X_DOT_BLOW_UP_RUN_RC`
    ! TODO: `M_BLOW_UP_RUN_RC`
    ! TODO: `E_BLOW_UP_RUN_RC`
    ! TODO: `E_F_BLOW_UP_RUN_RC`
end subroutine test_check_sys

!tripwire$ begin C75F7EAC Update `\secref{exact-solution}` of verval.tex when changing the exact solution test if necessary.
pure function exact_x_dot(sys_0, x)
    use checks, only: assert, is_close
    use cva, only: IDEAL_EOS, CONST_EOS, NORMAL_CV_TYPE, MIRROR_CV_TYPE, cv_system_type
    
    type(cv_system_type), intent(in), allocatable :: sys_0
    type(si_length), intent(in)                   :: x
    
    type(si_velocity) :: exact_x_dot
    
    type(si_area)         :: csa
    type(si_inverse_mass) :: rm_p
    type(si_pressure)     :: p_0, p_atm, p_f
    real(WP)              :: gamma
    type(si_length)       :: x_0
    
    call assert(size(sys_0%cv) == 2, "test_cva (exact_x_dot): sys_0%cv has the incorrect size")
    call assert(sys_0%cv(1)%type == MIRROR_CV_TYPE, "test_cva (exact_x_dot): CV 1 must be MIRROR_CV_TYPE")
    call assert(sys_0%cv(1)%eos  == CONST_EOS,      "test_cva (exact_x_dot): CV 1 must be CONST_EOS")
    call assert(size(sys_0%cv(1)%gas) == 1,         "test_cva (exact_x_dot): CV 1 must have only 1 gas")
    call assert(sys_0%cv(2)%type == NORMAL_CV_TYPE, "test_cva (exact_x_dot): CV 2 must be NORMAL_CV_TYPE")
    call assert(sys_0%cv(2)%eos  == IDEAL_EOS,      "test_cva (exact_x_dot): CV 2 must be IDEAL_EOS")
    call assert(size(sys_0%cv(2)%gas) == 1,         "test_cva (exact_x_dot): CV 2 must have only 1 gas")
    
    csa   = sys_0%cv(2)%csa
    rm_p  = sys_0%cv(2)%rm_p
    p_0   = sys_0%cv(2)%p()
    x_0   = sys_0%cv(2)%x
    gamma = sys_0%cv(2)%gas(1)%gamma
    p_atm = sys_0%cv(1)%p()
    p_f   = sys_0%cv(2)%p_fd
    
    call assert(is_close(sys_0%cv(2)%p_fd%v%v, sys_0%cv(2)%p_fs%v%v), &
                    "test_cva (exact_x_dot): p_fs must equal p_fd")
    
    exact_x_dot = sqrt(2.0_WP*csa*rm_p * (p_0*(x*((x_0/x)**gamma) - x_0)/(1.0_WP - gamma) &
                                        - (p_atm + p_f)*(x - x_0)))
end function exact_x_dot

subroutine exact_x_dot_de(n, ne, ne_d)
    use fmad, only: ad
    use gasdata, only: DRY_AIR
    use cva, only: TIMEOUT_RUN_RC, cv_system_type, run_config_type, run_status_type, run, &
                    ENERGY_DERIV_TOLERANCE_RUN_RC, ENERGY_DERIV_TOLERANCE
    use checks, only: assert
    
    integer, intent(in)                :: n
    type(ad), intent(out), allocatable :: ne(:)
    real(WP), intent(out), allocatable :: ne_d(:, :)
    
    type(run_config_type)             :: config
    type(cv_system_type), allocatable :: sys_start, sys_end
    type(run_status_type)             :: status
    
    integer, parameter   :: N_VAR = 1, N_D = 6
    integer              :: i_var, i_d
    type(si_area)        :: csa
    type(si_pressure)    :: p_atm, p_0, p_fs, p_fd
    type(si_temperature) :: temp_atm
    type(si_velocity)    :: x_dot, x_dot_exact
    type(si_length)      :: x_0, delta_pre
    type(si_mass)        :: m_p
    type(si_stiffness)   :: k
    type(unitless)       :: y(1)
    type(si_time)        :: t_stop, dt
    
    allocate(ne(N_VAR))
    allocate(ne_d(N_VAR, N_D))
    
    allocate(sys_start)
    allocate(sys_start%cv(2))
    allocate(sys_start%con(2, 2))
    
    sys_start%con(1, 1)%active = .false.
    sys_start%con(1, 2)%active = .false.
    sys_start%con(2, 1)%active = .false.
    sys_start%con(2, 2)%active = .false.
    
    ! 1: atmosphere for barrel
    
    call csa%v%init(2.0e-2_WP, 1, N_D)
    call p_atm%v%init(1.0e5_WP, 2, N_D)
    call temp_atm%v%init_const(300.0_WP, N_D)
    call y(1)%v%init_const(1.0_WP, N_D)
    
    call sys_start%cv(1)%set_const("atmosphere for barrel", csa, p_atm, temp_atm, [DRY_AIR], y, 2)
    
    ! 2: barrel
    
    call x_dot%v%init_const(0.0_WP, N_D) ! The derivative of this is unstable? I made it constant.
    call x_0%v%init(0.1_WP, 3, N_D)
    call p_0%v%init(10.0e5_WP, 4, N_D)
    call m_p%v%init(2.0_WP, 5, N_D)
    call p_fs%v%init(0.1e5_WP, 6, N_D)
    p_fd = p_fs
    call k%v%init_const(0.0_WP, N_D)
    call delta_pre%v%init_const(0.0_WP, N_D)
    
    call sys_start%cv(2)%set(x_0, x_dot, y, p_0, temp_atm, "barrel", csa, 1.0_WP/m_p, p_fs, p_fd, k, &
                                delta_pre, [DRY_AIR], 1, isentropic_filling=.true., p_atm=p_atm, constant_friction=.true.)
    
    call assert(sys_start%cv(2)%constant_friction, "test_exact, cv%constant_friction")
    
    call t_stop%v%init_const(TEST_EXACT_T_STOP, N_D)
    
    ! At first I thought that the number of time steps needed to be an integer.
    ! That was to avoid a reduction in order-of-accuracy from `sys_interp`.
    ! However, `sys_interp` is not called if `rc == TIMEOUT_RUN_RC` as it is here..
    dt = t_stop / real(n, WP)
    
    call config%set("test_exact", N_D, t_stop=t_stop, dt=dt, tolerance_checks=.false.)
    call run(config, sys_start, sys_end, status)
    
    x_dot_exact = exact_x_dot(sys_start, sys_end%cv(2)%x)
    
    if (status%rc == ENERGY_DERIV_TOLERANCE_RUN_RC) then
        print *, "ENERGY_DERIV_TOLERANCE_RUN_RC"
        print *, status%data(1), int(status%data(2)), ENERGY_DERIV_TOLERANCE
    end if
    call assert(status%rc == TIMEOUT_RUN_RC, "test_exact, status%rc")
    
    do i_var = 1, N_VAR
        select case (i_var)
            case (1)
                ne(i_var) = abs(sys_end%cv(2)%x_dot%v - x_dot_exact%v)
            case default
                error stop "test_cva (exact_x_dot_de): invalid i_var (1)"
        end select
        
        do i_d = 1, N_D
            select case (i_var)
                case (1)
                    ne_d(i_var, i_d) = abs(sys_end%cv(2)%x_dot%v%d(i_d) - x_dot_exact%v%d(i_d))
                case default
                    error stop "test_cva (exact_x_dot_de): invalid i_var (2)"
            end select
        end do
    end do
end subroutine exact_x_dot_de

subroutine test_exact(tests)
    ! Tests using exact solution with projectile and one internal control volume.
    ! There is also an additional constant control volume for the atmosphere, but no real calculations are done by that.
    ! Also tests that `cv%constant_friction = .true.` works.
    
    ! TODO: Can additionally test the following: `e`, `e_f` (easy), `p`, `temp`
    ! I'm not sure all of these will show 4th order accuracy, however, as not all are solved for directly via RK4.
    
    use fmad, only: ad
    use convergence, only: convergence_test
    use checks, only: assert
    use io, only: write_latex_engineering
    
    type(test_results_type), intent(in out) :: tests
    
    type(ad), allocatable :: ne_fine(:)
    real(WP), allocatable :: ne_d_fine(:, :), p(:), ne(:)
    integer  :: n(2), tex_unit
    real(WP) :: dt
    
    n  = [500, 1000]
    dt = TEST_EXACT_T_STOP / real(n(2), WP)
    
    call exact_x_dot_de(10000, ne_fine, ne_d_fine)
    
    call tests%real_eq(ne_fine(1)%v, 0.0_WP, "test_exact, x_dot numerical error", abs_tol=1.0e-12_WP)
    call tests%real_eq(ne_d_fine(1, 1), 0.0_WP, "test_exact, d(x_dot)/d(csa) numerical error", abs_tol=1.0e-11_WP)
    call tests%real_eq(ne_d_fine(1, 2), 0.0_WP, "test_exact, d(x_dot)/d(p_atm) numerical error", abs_tol=1.0e-18_WP)
    call tests%real_eq(ne_d_fine(1, 3), 0.0_WP, "test_exact, d(x_dot)/d(x_0) numerical error", abs_tol=1.0e-11_WP)
    call tests%real_eq(ne_d_fine(1, 4), 0.0_WP, "test_exact, d(x_dot)/d(p_0) numerical error", abs_tol=1.0e-18_WP)
    call tests%real_eq(ne_d_fine(1, 5), 0.0_WP, "test_exact, d(x_dot)/d(m_p) numerical error", abs_tol=1.0e-13_WP)
    call tests%real_eq(ne_d_fine(1, 6), 0.0_WP, "test_exact, d(x_dot)/d(p_f) numerical error", abs_tol=1.0e-18_WP)
    
    ! I guess that I'm running into floating point error if I make the time step smaller than around the default.
    ! 4th order accuracy sure convergences fast!
    call convergence_test(n, exact_x_dot_de, [4.0_WP], "test_exact, passing", tests, &
                            p_tol=[0.03_WP], p_d_tol=[0.14_WP], p=p, ne=ne)
    
    call assert(size(p) == 1, "test_cva (test_exact): size(p) == 1 violated", print_integer=[size(p)])
    open(newunit=tex_unit, action="write", status="replace", position="rewind", file="test_exact.tex", delim="quote")
    write(unit=tex_unit, fmt="(a)") "% auto-generated"
    call write_latex_engineering(tex_unit, dt, "testexactdt", "f4.1")
    call write_latex_engineering(tex_unit, ne(1), "xdoterror", "f5.3")
    write(unit=tex_unit, fmt="(a, f5.3, a)") "\newcommand*{\xdotorder}{", p(1), "}"
    close(tex_unit)
end subroutine test_exact
!tripwire$ end

end program test_cva
