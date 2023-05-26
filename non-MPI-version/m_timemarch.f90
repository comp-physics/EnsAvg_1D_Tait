!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module solving ODE (semi-discrete form of PDE), given RHS
!
!  Last update: February 11, 2009
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_timemarch

  USE m_globalvar
  USE m_rhs
  USE m_timesplit
  IMPLICIT NONE

  TYPE(timelevel), DIMENSION(:), ALLOCATABLE :: rk
  REAL(KIND(0.D0)), DIMENSION(4,0:4) :: ark
  REAL(KIND(0.D0)), DIMENSION(4,0:4) :: brk

  CONTAINS

  !========================================================================

  SUBROUTINE s_start_timevar

    INTEGER :: iv
    INTEGER :: ir
    INTEGER :: irk

    ! allocate storage
    IF ( timescheme=='euler' ) THEN
       ALLOCATE( rk(Nv) )
       DO iv = 1,Nv
          ALLOCATE( rk(iv)%l(1)%f(1:Nx) )
       END DO
    ELSE IF ( timescheme=='rk3tvd' ) THEN
       ALLOCATE( rk(Nv) )
       DO irk = 1,2
          DO iv = 1,Nv
             ALLOCATE( rk(iv)%l(irk)%f(1:Nx) )
          END DO
       END DO
    ELSE IF ( timescheme=='rk4tvd' ) THEN
       ALLOCATE( rk(Nv) )
       DO irk = 1,5
          DO iv = 1,Nv
             ALLOCATE( rk(iv)%l(irk)%f(1:Nx) )
          END DO
       END DO
    END IF

    ! time increment used for adaptive RK
    IF ( timesplit=='y' ) THEN
       ALLOCATE( dttry(NR0) )
       DO ir = 1,NR0
          ALLOCATE( dttry(ir)%f(1:Nx) )
       END DO
    END IF

    ! coefficients for RK4TVD
    ark(2,0) = 649.D0/1600.D0
    ark(2,1) = 951.D0/1600.D0
    ark(3,0) = 53989.D0/2500000.D0
    ark(3,1) = 4806213.D0/20000000.D0
    ark(3,2) = 23619.D0/32000.D0
    ark(4,0) = 2.D-1
    ark(4,1) = 6127.D0/30000.D0
    ark(4,2) = 7873.D0/30000.D0
    ark(4,3) = third
    brk(1,0) = 5.D-1
    brk(2,0) = -10890423.D0/25193600.D0
    brk(2,1) = 5000.D0/7873.D0
    brk(3,0) = -102261.D0/5000000.D0
    brk(3,1) = -5121.D0/20000.D0
    brk(3,2) = 7873.D0/10000.D0
    brk(4,0) = 1.D-1
    brk(4,1) = sixth
    brk(4,2) = 0.D0
    brk(4,3) = sixth

  END SUBROUTINE s_start_timevar

  !========================================================================

  SUBROUTINE s_euler

    INTEGER :: iv

    tvd = 1.D0
    coeff = 1.D0
    CALL s_rhs( rk%l(1) )
    DO iv = 1,Nv
       rk(iv)%l(1)%f = rk(iv)%l(1)%f + dt*rhs(iv)%f + dqbub(iv)%f
    END DO

  END SUBROUTINE s_euler

  !========================================================================

  SUBROUTINE s_rk3tvd

    INTEGER :: iv

    ! step 1
    tvd = 1.D0
    coeff = 1.D0
    CALL s_rhs( rk%l(1) )
    DO iv = 1,Nv
       rk(iv)%l(2)%f = rk(iv)%l(1)%f + dt*rhs(iv)%f + dqbub(iv)%f
    END DO
    ! step 2
    coeff = 0.25D0
    CALL s_rhs( rk%l(2) )
    DO iv = 1,Nv
       rk(iv)%l(2)%f = 0.75D0*rk(iv)%l(1)%f &
                     + 0.25D0*( rk(iv)%l(2)%f+dt*rhs(iv)%f ) + dqbub(iv)%f
    END DO
    ! step 3
    coeff = 2.D0*third
    CALL s_rhs( rk%l(2) )
    DO iv = 1,Nv
       rk(iv)%l(1)%f = ( rk(iv)%l(1)%f+2.D0*(rk(iv)%l(2)%f+dt*rhs(iv)%f) )&
                     * third + dqbub(iv)%f
    END DO

  END SUBROUTINE s_rk3tvd

  !========================================================================

  SUBROUTINE s_rk4tvd

    ! Timestep splitting is not included.
    INTEGER :: iv

    ! step 1a
    tvd = 1.D0
    CALL s_rhs( rk%l(1) )
    DO iv = 1,Nv
       rk(iv)%l(2)%f = rk(iv)%l(1)%f + brk(1,0)*dt*rhs(iv)%f
       rk(iv)%l(3)%f = ark(2,0)*rk(iv)%l(1)%f + ark(2,1)*rk(iv)%l(2)%f
       rk(iv)%l(4)%f = ark(3,0)*rk(iv)%l(1)%f + ark(3,1)*rk(iv)%l(2)%f
       rk(iv)%l(5)%f = ark(4,0)*rk(iv)%l(1)%f + ark(4,1)*rk(iv)%l(2)%f &
                     + brk(4,0)*dt*rhs(iv)%f
    END DO
    ! step 1b (adjoint)
    tvd = -1.D0
    CALL s_rhs( rk%l(1) )
    DO iv = 1,Nv
       rk(iv)%l(3)%f = rk(iv)%l(3)%f + brk(2,0)*dt*rhs(iv)%f
       rk(iv)%l(4)%f = rk(iv)%l(4)%f + brk(3,0)*dt*rhs(iv)%f
    END DO
    ! step 2a
    tvd = 1.D0
    CALL s_rhs( rk%l(2) )
    DO iv = 1,Nv
       rk(iv)%l(3)%f = rk(iv)%l(3)%f + brk(2,1)*dt*rhs(iv)%f
       rk(iv)%l(5)%f = rk(iv)%l(5)%f + brk(4,1)*dt*rhs(iv)%f &
                     + ark(4,2)*rk(iv)%l(3)%f
    END DO
    ! step 2b (adjoint)
    tvd = -1.D0
    CALL s_rhs( rk%l(2) )
    DO iv = 1,Nv
       rk(iv)%l(4)%f = rk(iv)%l(4)%f + brk(3,1)*dt*rhs(iv)%f &
                     + ark(3,2)*rk(iv)%l(3)%f
    END DO
    ! step 3
    tvd = 1.D0
    CALL s_rhs( rk%l(3) )
    DO iv = 1,Nv
       rk(iv)%l(4)%f = rk(iv)%l(4)%f + brk(3,2)*dt*rhs(iv)%f
    END DO
    ! step 4
    tvd = 1.D0
    CALL s_rhs( rk%l(4) )
    DO iv = 1,Nv
       rk(iv)%l(1)%f = rk(iv)%l(5)%f + ark(4,3)*rk(iv)%l(4)%f &
                     + brk(4,3)*dt*rhs(iv)%f
    END DO

  END SUBROUTINE s_rk4tvd

  !========================================================================

  SUBROUTINE s_stop_timevar

    DEALLOCATE( rk )
    IF ( timesplit=='y' ) DEALLOCATE( dttry )

  END SUBROUTINE s_stop_timevar

  !========================================================================

END MODULE m_timemarch
