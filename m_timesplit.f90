!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module for adaptive Runge-Kutta solving bubble-dynamic equations
!
!  Last update: May 13, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_timesplit

  USE mpi_setup
  USE mpi_transfer
  USE m_globalvar
  USE m_timemarch
  USE m_rhs
  USE m_rhsvar
  USE m_bubbles
  USE m_misc
  IMPLICIT NONE

  CONTAINS

  !========================================================================

  SUBROUTINE s_integsrc

    INTEGER :: i
    INTEGER :: ib
    INTEGER :: ir
    REAL(KIND(0.D0)) :: R3bar
    REAL(KIND(0.D0)), DIMENSION(NR0) :: R3
    REAL(KIND(0.D0)), DIMENSION(Nb) :: btmp

    ! returns dpldt
    CALL s_compdpldt( rk%l(1) )

    ! integrate bubble-dynamic sources with adaptive RK
    DO i = 1,Nx

       DO ir = 1,NR0
          ! substitute the initial values
          iR0 = ir
          DO ib = 1,Nb
             btmp(ib) = bub(ib,ir)%f(i)
          END DO
          CALL s_odeint( &
                btmp,Nb,pres(i),dpldt(i),dtsplit,accuracy,dttry(ir),0.D0 )
          ! unpdate
          DO ib = 1,Nb
             rk(ibub(ib,ir))%l(1)%f(i) = nbub(i)*btmp(ib)
          END DO
          R3(ir) = btmp(1)**3
       END DO
       CALL s_quad( R3,R3bar )
       rk(Nveul)%l(1)%f(i) = pi43*nbub(i)*R3bar

    END DO

  END SUBROUTINE s_integsrc

  !========================================================================

  SUBROUTINE s_compdpldt( qgvn )

    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: qgvn
    INTEGER :: iv

    ! transfer neighboring data
    DO iv = 1,Nv
       qval(iv)%f = qgvn(iv)%f
    END DO
    CALL s_mpi_transfer( qval )
    ! compute physical properties & characteristic speeds
    CALL s_compvar
    ! Roe-averaged right & left eigenvectors at celledges
    CALL s_roeaverage
    ! numerical flux
    CALL s_rhsfv
    ! dpl/dt
    CALL s_compsrc

  END SUBROUTINE s_compdpldt

  !========================================================================

  SUBROUTINE s_odeint( ystart,Nvar,pl,pldot,hsplit,acc,htry,hmin )

    ! Runge-Kutta driver:
    ! Integrate initial value ystart from x1=0 to x2=hsplit.
    ! Replace ystart with the new value at x2.
    ! htry is set as a guessed first stepsize.
    ! hmin as the minimum allowed stepsize (can be zero)
    INTEGER, INTENT(IN) :: Nvar
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), INTENT(IN) :: hsplit
    REAL(KIND(0.D0)), INTENT(IN) :: acc
    REAL(KIND(0.D0)), INTENT(IN) :: htry
    REAL(KIND(0.D0)), INTENT(IN) :: hmin
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(INOUT) :: ystart
    ! local variables
    INTEGER :: nstp
    INTEGER :: ivar
    INTEGER, PARAMETER :: MAXSTP = 10000
    REAL(KIND(0.D0)) :: x
    REAL(KIND(0.D0)) :: x1, x2
    REAL(KIND(0.D0)) :: h
    REAL(KIND(0.D0)) :: h1
    REAL(KIND(0.D0)) :: hdid
    REAL(KIND(0.D0)) :: hnext
    REAL(KIND(0.D0)), PARAMETER :: small = 1.D-30
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: y
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: dydx
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: yscal

    ! ICs
    x1 = 0.D0
    x2 = hsplit
    x = x1
    h = MIN( htry,hsplit )
    h1 = h
    y = ystart

    ! integrate
    DO nstp = 1,MAXSTP

       ! derivatives
       CALL s_derivs( y,pl,pldot,dydx,Nvar )
       ! scaling used for monitoring accuracy
       yscal = DABS( y ) + DABS( h*dydx ) + small
       ! correct stepsize for overshooting
       IF ( (x+h-x2)*(x+h-x1)>0.D0 ) h = x2 - x
       !!! ODE solvers !!!
       CALL s_rkqs( y,dydx,Nvar,pl,pldot,h,acc,yscal,hdid,hnext )
       ! NaN
       DO ivar = 1,Nvar
          IF ( y(ivar)/=y(ivar) ) THEN
             y = ystart
             x = x1
             hdid = 0.D0
             h1 = 0.1D0*h1
             hnext = h1
          END IF
       END DO
       ! judgement to exit
       IF ( (x-x2)*(x2-x1)>=0.D0 ) THEN
          ystart = y
          EXIT
       END IF
       IF ( DABS(hnext)<hmin ) THEN
          PRINT*, 'Stepsize smaller than minimum in s_odeint'
          STOP
       END IF
       x = x + hdid
       h = hnext

    END DO

  END SUBROUTINE s_odeint

  !========================================================================

  SUBROUTINE s_derivs( y,pl,pldot,ak,Nvar )

    INTEGER, INTENT(IN) :: Nvar
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(IN) :: y
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(OUT) :: ak
    REAL(KIND(0.D0)) :: vflux

    IF ( polytropic=='y' ) THEN

       ak(1) = y(2)
       ak(2) = f_gilmore( y(1),y(2),0.D0,0.D0,pl,pldot )

    ELSE IF ( polytropic=='n'.AND.vapor=='y' ) THEN

       CALL s_bwproperty( y(3) )
       ak(1) = y(2)
       vflux = f_vflux( y(1),y(2),y(4) )
       ak(4) = vflux*4.D0*pi*y(1)**2
       ak(3) = f_bpres_dot( vflux,y(1),y(2),y(3),y(4) )
       ak(2) = f_gilmore( y(1),y(2),y(3),ak(3),pl,pldot )

    ELSE IF ( polytropic=='n'.AND.vapor=='n' ) THEN

       CALL s_bwproperty( y(3) )
       ak(1) = y(2)
       vflux = f_vflux( y(1),y(2),0.D0 )
       ak(3) = f_bpres_dot( vflux,y(1),y(2),y(3),0.D0 )
       ak(2) = f_gilmore( y(1),y(2),y(3),ak(3),pl,pldot )

    END IF

  END SUBROUTINE s_derivs

  !========================================================================

  SUBROUTINE s_rkqs( y,dydx,Nvar,pl,pldot,htry,acc,yscal,hdid,hnext )

    ! fifth-order RK step with monitoring of local truncation error
    ! to ensure accuracy and adjust stepsize
    INTEGER, INTENT(IN) :: Nvar
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), INTENT(IN) :: htry
    REAL(KIND(0.D0)), INTENT(IN) :: acc
    REAL(KIND(0.D0)), INTENT(OUT) :: hdid
    REAL(KIND(0.D0)), INTENT(OUT) :: hnext
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(INOUT) :: y
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(IN) :: dydx
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(IN) :: yscal
    ! local variables
    INTEGER :: ivar
    REAL(KIND(0.D0)) :: h
    REAL(KIND(0.D0)) :: htmp
    REAL(KIND(0.D0)) :: errmax
    REAL(KIND(0.D0)), PARAMETER :: safety = 0.9D0
    REAL(KIND(0.D0)), PARAMETER :: pgrow = -0.2D0
    REAL(KIND(0.D0)), PARAMETER :: pshrnk = -0.25D0
    REAL(KIND(0.D0)), PARAMETER :: small = 1.D-30
    REAL(KIND(0.D0)), PARAMETER :: errcon = 1.89D-4
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: yerr
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ytmp

    h = htry
    DO

       CALL s_rkck( y,dydx,Nvar,pl,pldot,h,ytmp,yerr )
       errmax = 0.D0
       DO ivar = 1,Nvar
          errmax = MAX( errmax,DABS(yerr(ivar)/yscal(ivar)) )
       END DO
       errmax = errmax/acc
       IF ( errmax>1.D0 ) THEN
          htmp = safety*h*( errmax**pshrnk )
          h = DSIGN( MAX(DABS(htmp),0.1D0*DABS(h)),h )
          IF ( h==0.D0 ) THEN
             PRINT*, 'Stepsize underflow in s_rkqs'
             STOP
          END IF
       ELSE
          IF ( errmax>errcon ) THEN
             hnext = safety*h*( errmax**pgrow )
          ELSE
             hnext = 5.D0*h
          END IF
          hdid = h
          y = ytmp
          EXIT
       END IF

    END DO

  END SUBROUTINE s_rkqs

  !========================================================================

  SUBROUTINE s_rkck( y,dydx,Nvar,pl,pldot,h,yout,yerr )

    ! fifth-order Cash-Karp RK
    INTEGER, INTENT(IN) :: Nvar
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), INTENT(IN) :: h
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(IN) :: y
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(IN) :: dydx
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(OUT) :: yout
    REAL(KIND(0.D0)), DIMENSION(Nvar), INTENT(OUT) :: yerr
    ! local variables
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ytmp
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ak2
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ak3
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ak4
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ak5
    REAL(KIND(0.D0)), DIMENSION(Nvar) :: ak6
    REAL(KIND(0.D0)), PARAMETER :: &
      a2=0.2D0, a3=0.3D0, a4=0.6D0, a5=1.D0, a6=0.875D0,                &
      b21=0.2D0, b31=3.D0/40.D0, b32=9.D0/40.D0, b41=0.3D0, b42=-0.9D0, &
      b43=1.2D0, b51=-11.D0/54.D0, b52=2.5D0, b53=-70.D0/27.D0,         &
      b54=35.D0/27.D0, b61=1631.D0/55296.D0, b62=175.D0/512.D0,         &
      b63=575.D0/13824.D0, b64=44275.D0/110592.D0, b65=253.D0/4096.D0,  &
      c1=37.D0/378.D0, c3=250.D0/621.D0, c4=125.D0/594.D0,              &
      c6=512.D0/1771.D0, dc1=c1-2825.D0/27648.D0,                       &
      dc3=c3-18575.D0/48384.D0, dc4=c4-13525.D0/55296.D0,               &
      dc5=-277.D0/14336.D0, dc6=c6-0.25D0

    ytmp = y + b21*h*dydx
    CALL s_derivs( ytmp,pl,pldot,ak2,Nvar )
    ytmp = y + h*( b31*dydx+b32*ak2 )
    CALL s_derivs( ytmp,pl,pldot,ak3,Nvar )
    ytmp = y + h*( b41*dydx+b42*ak2+b43*ak3 )
    CALL s_derivs( ytmp,pl,pldot,ak4,Nvar )
    ytmp = y + h*( b51*dydx+b52*ak2+b53*ak3+b54*ak4 )
    CALL s_derivs( ytmp,pl,pldot,ak5,Nvar )
    ytmp = y + h*( b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5 )
    CALL s_derivs( ytmp,pl,pldot,ak6,Nvar )
    ! updated values
    yout = y + h*( c1*dydx+c3*ak3+c4*ak4+c6*ak6 )
    ! estimated error
    yerr = h*( dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6 )

  END SUBROUTINE s_rkck

  !========================================================================

END MODULE m_timesplit
