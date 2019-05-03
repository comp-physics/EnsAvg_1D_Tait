!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module containing miscellaneous subroutines
!
!  Last update: January 21, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_misc

  USE m_globalvar
  USE mpi_setup
  IMPLICIT NONE

  CONTAINS

  !========================================================================

  SUBROUTINE s_quad( func,mom )

    REAL(KIND(0.D0)), DIMENSION(NR0), INTENT(IN) :: func
    REAL(KIND(0.D0)), INTENT(OUT) :: mom

    mom = DOT_PRODUCT( weight,func )

  END SUBROUTINE s_quad

  !========================================================================

  SUBROUTINE s_comp_n( vftmp,nRtmp,ntmp )

    REAL(KIND(0.D0)), INTENT(IN) :: vftmp
    REAL(KIND(0.D0)), DIMENSION(NR0), INTENT(IN) :: nRtmp
    REAL(KIND(0.D0)), INTENT(OUT) :: ntmp
    REAL(KIND(0.D0)) :: nR3

    CALL s_quad( nRtmp**3,nR3 )
    ntmp = DSQRT( pi43*nR3/vftmp )

  END SUBROUTINE s_comp_n

  !========================================================================

  SUBROUTINE s_NaN_rhs( qgvn )

    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: qgvn
    INTEGER :: i
    INTEGER :: iv

    DO i = i_nopad_mn,i_nopad_mx
       DO iv = 1,Nv

          IF ( qgvn(iv)%f(i)/=qgvn(iv)%f(i) ) THEN
             PRINT*, 'NaN in flux (w/o ptilde) at (i,it,iv):', i, it, iv
             PRINT*, 'x =', xgrid(i)
             PRINT*, 'mpi_rank =', mpi_rank
             PRINT*, 'Computation stopped.'
             STOP
          END IF

       END DO
    END DO

  END SUBROUTINE s_NaN_rhs

  !========================================================================

  SUBROUTINE s_NaN_src( qgvn )

    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: qgvn
    INTEGER :: i
    INTEGER :: iv

    DO i = i_nopad_mn,i_nopad_mx
       DO iv = 1,Nv

          IF ( qgvn(iv)%f(i)/=qgvn(iv)%f(i) ) THEN
             PRINT*, 'NaN in sources (w/ ptilde) at (i,it,iv):', i, it, iv
             PRINT*, 'x =', xgrid(i)
             PRINT*, 'mpi_rank =', mpi_rank
             PRINT*, 'Computation stopped.'
             STOP
          END IF

       END DO
    END DO

  END SUBROUTINE s_NaN_src

  !========================================================================

  SUBROUTINE s_check_sign( qgvn )

    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: qgvn
    INTEGER :: i
    INTEGER :: iv

    DO i = i_nopad_mn,i_nopad_mx
       DO iv = Nveul,Nveul+1

          IF ( qgvn(iv)%f(i)<0.D0 ) THEN
             PRINT*, 'Negative sign at (i,it,iv):', i, it, iv
             PRINT*, 'x =', xgrid(i)
             PRINT*, 'mpi_rank =', mpi_rank
          END IF

       END DO
    END DO

  END SUBROUTINE s_check_sign

  !========================================================================

END MODULE m_misc
