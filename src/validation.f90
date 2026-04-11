! Module for validation testing procedures.
! Standard: Fortran 2018
! Preprocessor: none
! Author: Ben Trettel (<http://trettel.us/>)
! Project: [BlasterSim](https://github.com/btrettel/blastersim)
! License: [AGPLv3](https://www.gnu.org/licenses/agpl-3.0.en.html)

module validation

implicit none
private

public :: z_tail_cdf
public :: t_tail_cdf, student_t

contains

pure function z_tail_cdf(z)
    use prec, only: WP
    use checks, only: assert
    
    real(WP), intent(in) :: z
    
    real(WP) :: z_tail_cdf
    
    z_tail_cdf = 0.5_WP*(1.0_WP - erf(z/sqrt(2.0_WP)))
    
    call assert(z_tail_cdf <= 1.0_WP, "validation (z_tail_cdf): z_tail_cdf <= 1 violated", &
                    print_real=[z, z_tail_cdf])
end function z_tail_cdf

pure function t_tail_cdf(student_t, dof)
    use prec, only: WP
    use checks, only: assert
    
    ! ALGORITHM AS 27  APPL. STATIST. VOL.19, NO.1
    ! Calculate the upper tail area under Student's t-distribution.
    ! Translated from Algol by Alan Miller.
    
    ! Modified from <https://jblevins.org/mirror/amiller/as27.f90>.
    ! The constants are single precision, but this is fine for my purposes.
    
    ! Copyright:
    ! Arguably ineligible as it doesn't contain any creative expression.
    ! This is fair use otherwise.
    
    real(WP), intent(in) :: student_t
    integer, intent(in)  :: dof
    
    real(WP) :: t_tail_cdf

    real(WP) :: v, tt
    real(WP), parameter  :: A1 = 0.09979441_WP, A2 = -0.581821_WP, A3 = 1.390993_WP, &
                                    A4 = -1.222452_WP, A5 = 2.151185_WP
    real(WP), parameter  :: B1 = 5.537409_WP, B2 = 11.42343_WP
    real(WP), parameter  :: C1 = 0.04431742_WP, C2 = -0.2206018_WP, C3 = -0.03317253_WP, &
                                    C4 = 5.679969_WP, C5 = -12.96519_WP
    real(WP), parameter  :: D1 = 5.166733_WP, D2 = 13.49862_WP
    real(WP), parameter  :: E1 = 0.009694901_WP, E2 = -0.1408854_WP, E3 = 1.88993_WP, &
                                    E4 = -12.75532_WP, E5 = 25.77532_WP
    real(WP), parameter  :: F1 = 4.233736_WP, F2 = 14.3963_WP
    real(WP), parameter  :: G1 = -9.187228e-5_WP, G2 = 0.03789901_WP, G3 = -1.280346_WP, &
                                    G4 = 9.249528_WP, G5 = -19.08115_WP
    real(WP), parameter  :: H1 = 2.777816_WP, H2 = 16.46132_WP
    real(WP), parameter  :: I1 = 5.79602e-4_WP, I2 = -0.02763334_WP, I3 = 0.4517029_WP, &
                                    I4 = -2.657697_WP, I5 = 5.127212_WP
    real(WP), parameter  :: J1 = 0.5657187_WP, J2 = 21.83269_WP

    ! The original version said `dof > 4`, but my own testing shows that below 10, this is often inaccurate.
    call assert(dof >= 10, "validation (t_tail_cdf): dof >= 10 violated")

    v = 1.0_WP / real(dof, WP)
    tt = abs(student_t)
    t_tail_cdf = 0.5_WP*(1.0_WP + &
                tt*(((A1 + v*(A2 + v*(A3 + v*(A4 + v*A5)))) / (1.0_WP - v*(B1 - v*B2))) + &
                tt*(((C1 + v*(C2 + v*(C3 + v*(C4 + v*C5)))) / (1.0_WP - v*(D1 - v*D2))) + &
                tt*(((E1 + v*(E2 + v*(E3 + v*(E4 + v*E5)))) / (1.0_WP - v*(F1 - v*F2))) + &
                tt*(((G1 + v*(G2 + v*(G3 + v*(G4 + v*G5)))) / (1.0_WP - v*(H1 - v*H2))) + &
                tt*((I1 + v*(I2 + v*(I3 + v*(I4 + v*I5)))) / (1.0_WP - v*(J1 - v*J2))) )))))**(-8)
    
    if (student_t < 0.0_WP) then
        t_tail_cdf = 1.0_WP - t_tail_cdf
    end if
    
    call assert(t_tail_cdf <= 1.0_WP, "validation (t_tail_cdf): t_tail_cdf <= 1 violated", &
                    print_real=[student_t, t_tail_cdf], print_integer=[dof])
end function t_tail_cdf

pure function student_t(t_tail_cdf_, dof)
    use prec, only: WP
    use checks, only: assert, is_close
    
    real(WP), intent(in) :: t_tail_cdf_
    integer, intent(in)  :: dof
    
    real(WP) :: student_t
    
    integer, parameter :: MAX_ITERS = 100
    real(WP) :: student_t_i, student_t_im1, student_t_im2, &
                t_tail_cdf_i, t_tail_cdf_im1, t_tail_cdf_im2
    integer  :: i
    
    ! `student_t_im2` is based on the upper bound from Chebyshev's inequality.
    !student_t_im2  = 1.0_WP/sqrt(2.0_WP*abs(t_tail_cdf_))
    ! This one caused numerical problems.
    
    student_t_im2  = 3.0_WP
    t_tail_cdf_im2 = t_tail_cdf(student_t_im2, dof)
    
    student_t_im1  = 1.96_WP
    t_tail_cdf_im1 = t_tail_cdf(student_t_im1, dof)
    
    if (t_tail_cdf_ > 0.5_WP) then
        student_t_im2  = -student_t_im2
        t_tail_cdf_im2 = 1.0_WP - t_tail_cdf_im2
        
        student_t_im1  = -student_t_im1
        t_tail_cdf_im1 = 1.0_WP - t_tail_cdf_im1
    end if
    
    call assert(.not. is_close(student_t_im1, student_t_im2), &
                    "validation (student_t): student_t_im1 /= student_t_im2 violated")
    do i = 1, MAX_ITERS
        !print *, i, student_t_im1, t_tail_cdf_im1, student_t_im2, t_tail_cdf_im2, t_tail_cdf_
        if (abs(t_tail_cdf_im1 - t_tail_cdf_im2) < 100.0_WP*spacing(t_tail_cdf_)) exit
        
        student_t_i  = student_t_im1 + (student_t_im1 - student_t_im2) * (t_tail_cdf_ - t_tail_cdf_im1) &
                                            / (t_tail_cdf_im1 - t_tail_cdf_im2)
        t_tail_cdf_i = t_tail_cdf(student_t_i, dof)
        
        student_t_im2 = student_t_im1
        student_t_im1 = student_t_i
        
        t_tail_cdf_im2 = t_tail_cdf_im1
        t_tail_cdf_im1 = t_tail_cdf_i
    end do
    
    student_t = student_t_i
end function student_t

end module validation
