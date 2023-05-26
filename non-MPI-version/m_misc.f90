!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module containing miscellaneous subroutines
!
!  Last update: March 19, 2009
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_misc

  USE m_globalvar
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
    ! vftemp is \alpha, nRtemp is n*R, ntmp is n
    ! compute n from \alpha and nR
    REAL(KIND(0.D0)), INTENT(IN) :: vftmp
    REAL(KIND(0.D0)), DIMENSION(NR0), INTENT(IN) :: nRtmp
    REAL(KIND(0.D0)), INTENT(OUT) :: ntmp
    REAL(KIND(0.D0)) :: nR3

    CALL s_quad( nRtmp**3,nR3 )  !returns itself if NR0 = 1
    ntmp = DSQRT( pi43*nR3/vftmp )
    
    !print*, 'nRtmp^3, nR3', nRtmp(:)**3., nR3

  END SUBROUTINE s_comp_n

  !========================================================================

END MODULE m_misc
