!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module for adaptive Runge-Kutta solving bubble-dynamic equations
!
!  Last update: March 31, 2009
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_timesplit

  USE m_globalvar
  USE m_bubbles
  IMPLICIT NONE

  REAL(KIND(0.D0)) :: coeff

  CONTAINS

  !========================================================================

  SUBROUTINE s_adaptiveRK( bubold,pl,pldot,dt_try,bubnew )

    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(IN) :: bubold
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), INTENT(INOUT) :: dt_try
    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(OUT) :: bubnew

    REAL(KIND(0.D0)), DIMENSION(Nb) :: dbdt
    REAL(KIND(0.D0)), DIMENSION(Nb) :: btmp
    REAL(KIND(0.D0)), DIMENSION(Nb) :: bubtmp
    REAL(KIND(0.D0)) :: t_each
    REAL(KIND(0.D0)) :: dt_done

    t_each = 0.D0
    bubtmp = bubold
    DO

       btmp = bubtmp
       CALL s_derivs( btmp,pl,pldot,dbdt )
       CALL s_rkqs( btmp,dbdt,pl,pldot,dt_try,dt_done )
       t_each = t_each + dt_done
       ! compute values at t_each = coeff*dt
       IF ( t_each>coeff*dt ) THEN
          dt_try = coeff*dt + dt_done - t_each
          CALL s_rkqs( bubtmp,dbdt,pl,pldot,dt_try,dt_done )
          bubnew = bubtmp
          EXIT
       ELSE IF ( t_each==coeff*dt ) THEN
          bubnew = btmp
          EXIT
       ELSE
          bubtmp = btmp
       END IF

    END DO

  END SUBROUTINE s_adaptiveRK

  !========================================================================

  SUBROUTINE s_derivs( y,pl,pldot,ak )

    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(IN) :: y
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(OUT) :: ak
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

  SUBROUTINE s_rkqs( y,dydx,pl,pldot,htry,hdid )

    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(IN) :: dydx
    REAL(KIND(0.D0)), INTENT(OUT) :: hdid
    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(INOUT) :: y
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), INTENT(INOUT) :: htry

    INTEGER :: ib
    REAL(KIND(0.D0)) :: errmax
    REAL(KIND(0.D0)) :: h
    REAL(KIND(0.D0)) :: htmp
    REAL(KIND(0.D0)), PARAMETER :: safety = 0.9D0
    REAL(KIND(0.D0)), PARAMETER :: pgrow = -0.2D0
    REAL(KIND(0.D0)), PARAMETER :: pshrnk = -0.25D0
    REAL(KIND(0.D0)), PARAMETER :: small = 1.D-30
    REAL(KIND(0.D0)), PARAMETER :: errcon = 1.89D-4
    REAL(KIND(0.D0)), PARAMETER :: acc = 1.D-5
    REAL(KIND(0.D0)), DIMENSION(Nb) :: yscal
    REAL(KIND(0.D0)), DIMENSION(Nb) :: yerr
    REAL(KIND(0.D0)), DIMENSION(Nb) :: ytmp

    h = htry
    yscal = DABS( y ) + DABS( h*dydx ) + small
    DO

       CALL s_rkck( y,dydx,pl,pldot,h,ytmp,yerr )
       errmax = 0.D0
       DO ib = 1,Nb
          errmax = MAX( errmax,DABS(yerr(ib)/yscal(ib)) )
       END DO
       errmax = errmax/acc
       IF ( errmax>1.D0 ) THEN
          htmp = safety*h*( errmax**pshrnk )
          h = DSIGN( MAX(DABS(htmp),0.1D0*DABS(h)),h )
       ELSE
          IF ( errmax>errcon ) THEN
             htry = safety*h*( errmax**pgrow )
          ELSE
             htry = 5.D0*h
          END IF
          hdid = h
          y = ytmp
          EXIT
       END IF

    END DO

  END SUBROUTINE s_rkqs

  !========================================================================

  SUBROUTINE s_rkck( y,dydx,pl,pldot,h,yout,yerr )

    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(IN) :: y
    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(IN) :: dydx
    REAL(KIND(0.D0)), INTENT(IN) :: pl
    REAL(KIND(0.D0)), INTENT(IN) :: pldot
    REAL(KIND(0.D0)), INTENT(IN) :: h
    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(OUT) :: yout
    REAL(KIND(0.D0)), DIMENSION(Nb), INTENT(OUT) :: yerr

    REAL(KIND(0.D0)), DIMENSION(Nb) :: ytmp
    REAL(KIND(0.D0)), DIMENSION(Nb) :: ak2
    REAL(KIND(0.D0)), DIMENSION(Nb) :: ak3
    REAL(KIND(0.D0)), DIMENSION(Nb) :: ak4
    REAL(KIND(0.D0)), DIMENSION(Nb) :: ak5
    REAL(KIND(0.D0)), DIMENSION(Nb) :: ak6
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
    CALL s_derivs( ytmp,pl,pldot,ak2 )
    ytmp = y + h*( b31*dydx+b32*ak2 )
    CALL s_derivs( ytmp,pl,pldot,ak3 )
    ytmp = y + h*( b41*dydx+b42*ak2+b43*ak3 )
    CALL s_derivs( ytmp,pl,pldot,ak4 )
    ytmp = y + h*( b51*dydx+b52*ak2+b53*ak3+b54*ak4 )
    CALL s_derivs( ytmp,pl,pldot,ak5 )
    ytmp = y + h*( b61*dydx+b62*ak2+b63*ak3+b64*ak4+b65*ak5 )
    CALL s_derivs( ytmp,pl,pldot,ak6 )
    yout = y + h*( c1*dydx+c3*ak3+c4*ak4+c6*ak6 )
    yerr = h*( dc1*dydx+dc3*ak3+dc4*ak4+dc5*ak5+dc6*ak6 )

  END SUBROUTINE s_rkck

  !========================================================================

END MODULE m_timesplit
