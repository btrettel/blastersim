! tests for the port module
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
call test_p_f_1(tests)
call test_p_f_2(tests)
call test_p_f0_1(tests)
call test_p_f0_2(tests)
call test_temp_cv(tests)

call tests%end_tests()

contains

subroutine test_p_eos(tests)
    use units, only: si_mass_density => unit_m30_p10_p00_p00, &
                     si_temperature  => unit_p00_p00_p00_p10, &
                     si_pressure     => unit_m10_p10_m20_p00
    use prec, only: WP
    use cva, only: P_ATM, T_ATM, RHO_ATM, p_eos
    
    type(test_results_type), intent(in out) :: tests

    type(si_mass_density) :: rho
    type(si_temperature)  :: temp
    type(si_pressure)     :: p
    
    call rho%v%init_const(RHO_ATM, 0)
    call temp%v%init_const(T_ATM, 0)
    
    p = p_eos(rho, temp)
    
    call tests%real_eq(p%v%v, P_ATM, "p_eos, atmospheric", abs_tol=1.0_WP)
end subroutine test_p_eos

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
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(-10.0_WP, 0)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, -cv%p_fd%v%v, "p_f, x_dot < 0 (dynamic friction), p_fe == 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(-10.0_WP, 0)
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, -cv%p_fd%v%v, "p_f, x_dot < 0 (dynamic friction), p_fe > 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(0.0_WP, 0)
    p_fe = -cv%p_fs/10.0_WP
    p_f  = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, p_fe%v%v, "p_f, x_dot == 0 (static friction), p_fe < 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(0.0_WP, 0)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, 0.0_WP, "p_f, x_dot == 0 (static friction), p_fe == 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(0.0_WP, 0)
    p_fe = cv%p_fs/10.0_WP
    p_f  = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, p_fe%v%v, "p_f, x_dot == 0 (static friction), p_fe > 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(10.0_WP, 0)
    call p_fe%v%init_const(-1.0e5_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, cv%p_fd%v%v, "p_f, x_dot > 0 (dynamic friction), p_fe < 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
    call cv%x_dot%v%init_const(10.0_WP, 0)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f = cv%p_f(p_fe)
    call tests%real_eq(p_f%v%v, cv%p_fd%v%v, "p_f, x_dot > 0 (dynamic friction), p_fe == 0")
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
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
    
    deallocate(cv%x_dot%v%d)
    deallocate(p_fe%v%d)
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
    
    deallocate(p_fe%v%d)
    p_fe = -cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, p_fe%v%v, "p_f0 (2)")
    
    deallocate(p_fe%v%d)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (3)")
    
    deallocate(p_fe%v%d)
    p_fe = cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, p_fe%v%v, "p_f0 (4)")
    
    deallocate(p_fe%v%d)
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
    
    deallocate(p_fe%v%d)
    p_fe = -cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (2)")
    
    deallocate(p_fe%v%d)
    call p_fe%v%init_const(0.0_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (3)")
    
    deallocate(p_fe%v%d)
    p_fe = cv%p_fs/10.0_WP
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (4)")
    
    deallocate(p_fe%v%d)
    call p_fe%v%init_const(1.0e5_WP, 0)
    p_f0 = cv%p_f0(p_fe)
    call tests%real_eq(p_f0%v%v, 0.0_WP, "p_f0 (5)")
end subroutine test_p_f0_2

! TODO: Plot `p_f0` to test it.

subroutine test_temp_cv(tests)
    use units, only: si_temperature => unit_p00_p00_p00_p10
    use cva, only: cv_type
    
    type(test_results_type), intent(in out) :: tests

    type(cv_type)        :: cv
    type(si_temperature) :: temp
    
    ! 1 kg of mass
    call cv%m%v%init_const(1.0_WP, 0)
    
    ! moran_fundamentals_2008 table A-22, 300 K with 1 kg of mass
    call cv%e%v%init_const(214.07e3_WP, 0)
    
    temp = cv%temp()
    call tests%real_eq(temp%v%v, 300.0_WP, "temp_cv (qualitative)", abs_tol=5.0_WP)
end subroutine test_temp_cv

end program test_cva
