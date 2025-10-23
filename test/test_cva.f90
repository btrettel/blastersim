! tests for control volume analysis
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

program test_cva

use prec, only: WP
use unittest, only: test_results_type
implicit none

type(test_results_type) :: tests

call tests%start_tests("cva.nml")

call test_p_eos(tests)
call test_rho_eos(tests)
call test_u_h(tests)
call test_p_f_1(tests)
call test_p_f_2(tests)
call test_p_f0_1(tests)
call test_p_f0_2(tests)
call test_temp_cv(tests)
call test_set(tests)
call test_smooth_min(tests)
call test_f_m_dot(tests)
call test_g_m_dot(tests)
! TODO: call test_m_dot(tests)

call tests%end_tests()

contains

subroutine test_p_eos(tests)
    use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                     si_temperature  => unit_p00_p00_p00_p10, &
                     si_pressure     => unit_m10_p10_m20_p00
    use prec, only: WP
    use cva, only: P_ATM, T_ATM, RHO_ATM, AIR, cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    type(cv_type)         :: cv
    
    call rho%v%init_const(RHO_ATM, 0)
    call temp%v%init_const(T_ATM, 0)
    call cv%m_1%v%init_const(1.0_WP, 0)
    call cv%m_2%v%init_const(0.0_WP, 0)
    
    cv%gas_1 = AIR
    cv%gas_2 = AIR
    
    p = cv%p_eos(rho, temp)
    
    call tests%real_eq(p%v%v, P_ATM, "p_eos, atmospheric", abs_tol=1.0_WP)
end subroutine test_p_eos

subroutine test_rho_eos(tests)
    use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                     si_temperature  => unit_p00_p00_p00_p10, &
                     si_pressure     => unit_m10_p10_m20_p00
    use prec, only: WP
    use cva, only: P_ATM, T_ATM, RHO_ATM, AIR, cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    type(cv_type)         :: cv
    
    call p%v%init_const(P_ATM, 0)
    call temp%v%init_const(T_ATM, 0)
    call cv%m_1%v%init_const(1.0_WP, 0)
    call cv%m_2%v%init_const(0.0_WP, 0)
    
    cv%gas_1 = AIR
    cv%gas_2 = AIR
    
    rho = cv%rho_eos(p, temp)
    
    call tests%real_eq(rho%v%v, RHO_ATM, "rho_eos, atmospheric", abs_tol=1.0e-4_WP)
end subroutine test_rho_eos

subroutine test_u_h(tests)
    use units, only: si_temperature     => unit_p00_p00_p00_p10, &
                     si_specific_energy => unit_p20_p00_m20_p00
    use prec, only: WP
    use cva, only: AIR
    
    type(test_results_type), intent(in out) :: tests

    type(si_temperature)     :: temp
    type(si_specific_energy) :: u, h
    
    call temp%v%init_const(400.0_WP, 0)
    
    u = AIR%u(temp)
    h = AIR%h(temp)
    
    ! Data from moran_fundamentals_2008 table A-22.
    call tests%real_eq(u%v%v, 286.16e3_WP, "u, 400 K", abs_tol=1.0e3_WP)
    call tests%real_eq(h%v%v, 400.98e3_WP, "h, 400 K", abs_tol=1.0e3_WP)
end subroutine test_u_h

subroutine test_p_f_1(tests)
    use units, only: si_pressure => unit_m10_p10_m20_p00
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
    use units, only: si_pressure => unit_m10_p10_m20_p00
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
    use units, only: si_pressure => unit_m10_p10_m20_p00
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
    use units, only: si_pressure => unit_m10_p10_m20_p00
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
    use units, only: si_temperature => unit_p00_p00_p00_p10
    use cva, only: AIR, cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)        :: cv
    type(si_temperature) :: temp
    
    ! 1 kg of mass
    call cv%m_1%v%init_const(1.0_WP, 0)
    call cv%m_2%v%init_const(0.0_WP, 0)
    
    ! air
    cv%gas_1 = AIR
    cv%gas_2 = AIR
    
    ! moran_fundamentals_2008 table A-22, 300 K with 1 kg of mass
    call cv%e%v%init_const(214.07e3_WP, 0)
    
    temp = cv%temp()
    call tests%real_eq(temp%v%v, 300.0_WP, "temp_cv (qualitative)", abs_tol=5.0_WP)
end subroutine test_temp_cv

subroutine test_set(tests)
    use units, only: si_length          => unit_p10_p00_p00_p00, &
                     si_velocity        => unit_p10_p00_m10_p00, &
                     unitless           => unit_p00_p00_p00_p00, &
                     si_mass            => unit_p00_p10_p00_p00, &
                     si_energy          => unit_p20_p10_m20_p00, &
                     si_area            => unit_p20_p00_p00_p00, &
                     si_pressure        => unit_m10_p10_m20_p00, &
                     si_stiffness       => unit_p00_p10_m20_p00, &
                     si_volume          => unit_p30_p00_p00_p00, &
                     si_mass_density    => unit_m30_p10_p00_p00, &
                     si_temperature     => unit_p00_p00_p00_p10, &
                     si_specific_energy => unit_p20_p00_m20_p00
    use cva, only: P_ATM, T_ATM, RHO_ATM, AIR, cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type) :: cv
    
    type(si_length)          :: x
    type(si_velocity)        :: x_dot
    type(unitless)           :: y_1
    type(si_pressure)        :: p, p_cv
    type(si_temperature)     :: temp, temp_cv
    type(si_area)            :: csa
    type(si_mass)            :: m_p
    type(si_pressure)        :: p_fs, p_fd
    type(si_stiffness)       :: k
    type(si_length)          :: x_z
    type(si_specific_energy) :: u
    
    type(si_volume)       :: vol_cv
    type(si_mass_density) :: rho_cv
    
    call x%v%init_const(0.5_WP, 0)
    call x_dot%v%init_const(0.5_WP, 0)
    call y_1%v%init_const(1.0_WP, 0)
    call p%v%init_const(2.0_WP*P_ATM, 0)
    call temp%v%init_const(sqrt(2.0_WP)*T_ATM, 0)
    call csa%v%init_const(0.1_WP, 0)
    call m_p%v%init_const(0.5_WP, 0)
    call p_fs%v%init_const(0.2e5_WP, 0)
    call p_fd%v%init_const(0.1e5_WP, 0)
    call k%v%init_const(10.0_WP, 0)
    call x_z%v%init_const(3.0_WP, 0)
    
    call cv%set(x, x_dot, y_1, p, temp, csa, m_p, p_fs, p_fd, k, x_z, AIR, AIR)
    
    call tests%real_eq(cv%x%v%v, x%v%v, "set, x")
    call tests%real_eq(cv%x_dot%v%v, x_dot%v%v, "set, x_dot")
    ! no `p` or `temp` member variables
    call tests%real_eq(cv%csa%v%v, csa%v%v, "set, csa")
    call tests%real_eq(cv%rm_p%v%v, 1.0_WP/m_p%v%v, "set, rm_p")
    call tests%real_eq(cv%p_fs%v%v, p_fs%v%v, "set, p_fs")
    call tests%real_eq(cv%p_fd%v%v, p_fd%v%v, "set, p_fd")
    call tests%real_eq(cv%k%v%v, k%v%v, "set, k")
    call tests%real_eq(cv%x_z%v%v, x_z%v%v, "set, x_z")
    
    temp_cv = cv%temp()
    call tests%real_eq(temp_cv%v%v, temp%v%v, "set, temp")
    
    vol_cv = cv%vol()
    call tests%real_eq(vol_cv%v%v, 0.05_WP, "set, vol")
    
    rho_cv = cv%rho()
    call tests%real_eq(rho_cv%v%v, sqrt(2.0_WP)*RHO_ATM, "set, rho", abs_tol=1.0e-4_WP)
    
    p_cv = cv%p()
    call tests%real_eq(p_cv%v%v, p%v%v, "set, p")
    
    u = AIR%u(temp)
    call tests%real_eq(cv%m_1%v%v, sqrt(2.0_WP)*RHO_ATM*0.05_WP, "set, m_1", abs_tol=1.0e-6_WP)
    call tests%real_eq(cv%m_2%v%v, 0.0_WP, "set, m_2")
    call tests%real_eq(cv%e%v%v, u%v%v*sqrt(2.0_WP)*RHO_ATM*0.05_WP, "set, e", abs_tol=1.0_WP)
end subroutine test_set

subroutine test_smooth_min(tests)
    use units, only: unitless => unit_p00_p00_p00_p00
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
    use units, only: unitless => unit_p00_p00_p00_p00
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
    use units, only: unitless => unit_p00_p00_p00_p00
    use cva, only: g_m_dot
    
    type(test_results_type), intent(in out) :: tests
    
    type(unitless) :: p_r, g
    
    call p_r%v%init_const(0.0_WP, 0)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 0.0_WP, "g_m_dot (1)")
    
    call p_r%v%init_const(0.9_WP, 0)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 0.0_WP, "g_m_dot (2)")
    
    call p_r%v%init_const(1.0_WP, 0)
    g = g_m_dot(p_r)
    call tests%real_eq(g%v%v, 1.0_WP, "g_m_dot (3)", abs_tol=1.0e-6_WP)
end subroutine test_g_m_dot

! TODO: Plot `g_m_dot` to test it.

!subroutine test_m_dot(tests)
!    use cva, only: m_dot
    
!    type(test_results_type), intent(in out) :: tests
    
!    ! TODO: test to check that m_dot \varpropto \Delta p at small \Delta p
!end subroutine test_m_dot

end program test_cva
