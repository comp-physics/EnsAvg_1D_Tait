!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module computing RHS of the PDE in semi-discrete form
!
!  Last update: April 22, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_rhs

  USE mpi_setup
  USE mpi_transfer
  USE m_globalvar
  USE m_rhsvar
  USE m_bubbles
  USE m_misc

  IMPLICIT NONE
  REAL(KIND(0.D0)) :: duwall

  CONTAINS

  !========================================================================

  SUBROUTINE s_rhs( qgvn )

    INTEGER :: iv
    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: qgvn

    ! transfer neighboring data
    DO iv = 1,Nv
       qval(iv)%f = qgvn(iv)%f
    END DO
    CALL s_mpi_transfer( qval )

    ! compute primitive variables & characteristic speeds
    CALL s_compvar
    ! Roe-averaged right & left eigenvectors at cell edges
    CALL s_roeaverage

    ! numerical flux
    CALL s_rhsfv
    ! source flux for momentum equation & sources in bubble dynamics
    DO iv = 1,Nv
       rhsrc(iv)%f = 0.D0
    END DO
    IF ( source=='y' ) THEN
       CALL s_compsrc
    END IF

    ! RHS including sources
    DO iv = 1,Nv
       rhs(iv)%f = rhsfx(iv)%f + rhsrc(iv)%f
    END DO

    ! check Not-a-Number
    CALL s_NaN_rhs( rhsfx )
    CALL s_NaN_src( rhsrc )

  END SUBROUTINE s_rhs

  !========================================================================

  SUBROUTINE s_rhsfv

    INTEGER :: i
    INTEGER :: iv

    ! substitute conserved variables for WENO reconstruction
    DO iv = 1,Nv
       DO i = 1,N1
          inweno(iv)%f(i) = qval(iv)%f(i)
       END DO
    END DO
    IF ( bound%beg=='reflect'.AND.mpi_rank==0 ) THEN
       DO iv = 1,Nv
          IF ( iv==dir1 ) THEN
             DO i = 1,wenonum
                ! reflected
                ! moving wall BC incorporated
                inweno(iv)%f(1-i) = 2.D0*uwall - inweno(iv)%f(i)
             END DO
          ELSE
             DO i = 1,wenonum
                ! mirrored
                inweno(iv)%f(1-i) = inweno(iv)%f(i)
             END DO
          END IF
       END DO
    END IF
    IF ( bound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN
       DO iv = 1,Nv
          IF ( iv==dir1 ) THEN
             DO i = 1,wenonum
                ! reflected
                inweno(iv)%f(N1+i) = -inweno(iv)%f(N1+1-i)
             END DO
          ELSE
             DO i = 1,wenonum
                ! mirrored
                inweno(iv)%f(N1+i) = inweno(iv)%f(N1+1-i)
             END DO
          END IF
       END DO
    END IF

    ! WENO reconstruction with monotonicity preserving
    CALL s_weno

    ! characteristic fields => physical space
    IF ( chardecomp=='y' ) THEN
       CALL s_multiply_rr( lweno,0,-1 )
       lweno = v_out ! "R" state
       CALL s_multiply_rr( rweno,0,0 )
       rweno = v_out ! "L" state
    END IF

    ! correct negative values
    IF ( negative=='y' ) CALL s_reconstruct1

    ! approximate Riemann solver (HLLC)
    ! boundary conditions implemented
    CALL s_approx_riemann

  END SUBROUTINE s_rhsfv

  !========================================================================

  SUBROUTINE s_reconstruct1

    INTEGER :: i
    INTEGER :: iv
    INTEGER :: ir
    REAL(KIND(0.D0)), DIMENSION(NR0) :: bval

    !!! correct lweno (1st-order reconstruction) !!!
    ! mixture density
!    DO i = 1,N1
!       IF ( lweno(1)%f(i)/=lweno(1)%f(i) ) THEN
!          DO iv = 1,Nv
!             lweno(iv)%f(i) = qval(iv)%f(i)
!          END DO
!       END IF
!    END DO
    ! void fraction
    DO i = 1,N1
       IF ( lweno(Nveul)%f(i)<=0.D0 ) THEN
          DO iv = 1,Nv
             lweno(iv)%f(i) = qval(iv)%f(i)
          END DO
       END IF
    END DO
    ! radius
    DO i = 1,N1
       DO ir = 1,NR0
          bval(ir) = lweno(ibub(1,ir))%f(i)
       END DO
       IF ( ANY(bval<=0.D0) ) THEN
          DO iv = 1,Nv
             lweno(iv)%f(i) = qval(iv)%f(i)
          END DO
       END IF
    END DO
    ! bubble pressure
    IF ( polytropic=='n' ) THEN
       DO i = 1,N1
          DO ir = 1,NR0
             bval(ir) = lweno(ibub(3,ir))%f(i)
          END DO
          IF ( ANY(bval<=0.D0) ) THEN
             DO iv = 1,Nv
                lweno(iv)%f(i) = qval(iv)%f(i)
             END DO
          END IF
       END DO
    END IF
    ! mass of vapor
    IF ( polytropic=='n'.AND.vapor=='y' ) THEN
       DO i = 1,N1
          DO ir = 1,NR0
             bval(ir) = lweno(ibub(4,ir))%f(i)
          END DO
          IF ( ANY(bval<=0.D0) ) THEN
             DO iv = 1,Nv
                lweno(iv)%f(i) = qval(iv)%f(i)
             END DO
          END IF
       END DO
    END IF

    !!! correct rweno (1st-order reconstruction) !!!
    ! mixture density
!    DO i = 1,N1
!       IF ( rweno(1)%f(i)/=rweno(1)%f(i) ) THEN
!          DO iv = 1,Nv
!             rweno(iv)%f(i) = qval(iv)%f(i)
!          END DO
!       END IF
!    END DO
    ! void fraction
    DO i = 1,N1
       IF ( rweno(Nveul)%f(i)<=0.D0 ) THEN
          DO iv = 1,Nv
             rweno(iv)%f(i) = qval(iv)%f(i)
          END DO
       END IF
    END DO
    ! radius
    DO i = 1,N1
       DO ir = 1,NR0
          bval(ir) = rweno(ibub(1,ir))%f(i)
       END DO
       IF ( ANY(bval<=0.D0) ) THEN
          DO iv = 1,Nv
             rweno(iv)%f(i) = qval(iv)%f(i)
          END DO
       END IF
    END DO
    ! bubble pressure
    IF ( polytropic=='n' ) THEN
       DO i = 1,N1
          DO ir = 1,NR0
             bval(ir) = rweno(ibub(3,ir))%f(i)
          END DO
          IF ( ANY(bval<=0.D0) ) THEN
             DO iv = 1,Nv
                rweno(iv)%f(i) = qval(iv)%f(i)
             END DO
          END IF
       END DO
    END IF
    ! mass of vapor
    IF ( polytropic=='n'.AND.vapor=='y' ) THEN
       DO i = 1,N1
          DO ir = 1,NR0
             bval(ir) = rweno(ibub(4,ir))%f(i)
          END DO
          IF ( ANY(bval<=0.D0) ) THEN
             DO iv = 1,Nv
                rweno(iv)%f(i) = qval(iv)%f(i)
             END DO
          END IF
       END DO
    END IF

  END SUBROUTINE s_reconstruct1

  !========================================================================

  SUBROUTINE s_compsrc

    INTEGER :: i
    INTEGER :: ir
    REAL(KIND(0.D0)) :: ra
    REAL(KIND(0.D0)) :: vv
    REAL(KIND(0.D0)) :: pb
    REAL(KIND(0.D0)) :: mv
    REAL(KIND(0.D0)) :: vflux
    REAL(KIND(0.D0)) :: pbdot
    REAL(KIND(0.D0)) :: pldot
    REAL(KIND(0.D0)) :: srcvf
    REAL(KIND(0.D0)) :: divum
    REAL(KIND(0.D0)) :: R2Vbar
    REAL(KIND(0.D0)), DIMENSION(NR0) :: R2V

    ! compute numerical source flux in momentum equation
    ! returns nsfx(Nbeg:N1) out of nsfx(0:N1)
    CALL s_hllc_src
    ! for nonreflective BCs
    IF ( xincoming%end=='undex' .or. xincoming%beg=='undex' ) THEN
       IF ( bound%beg=='nonreflect'.AND.mpi_rank==0 ) THEN
          nsfx(0) = nsfx(1)
          velhf(0) = velhf(1)
       END IF
       IF ( bound%end=='nonreflect'.AND.mpi_rank==mpi_size-1 ) THEN
          nsfx(N1) = nsfx(N1-1)
          velhf(N1) = velhf(N1-1)
       END IF
    ELSE
       IF ( bound%beg=='nonreflect'.AND.mpi_rank==0 ) THEN
          ! returns nsfx(0) and nsfx(1)
          CALL s_srcbeg
       END IF
       IF ( bound%end=='nonreflect'.AND.mpi_rank==mpi_size-1 ) THEN
          ! returns nsfx(N1) and nsfx(N1-1)
          CALL s_srcend
       END IF
    END IF
    ! source in momentum equation
    DO i = 1,N1
       rhsrc(dir1)%f(i) = ( nsfx(i)-nsfx(i-1) )/ds(i)
    END DO

    ! sources for bubble-dynamic equations
    IF ( timesplit=='unsplit' ) THEN

       IF ( polytropic=='n'.AND.vapor=='y' ) THEN
          DO i = 1,N1
             DO ir = 1,NR0
                R2V(ir) = bub(1,ir)%f(i)**2*bub(2,ir)%f(i)
             END DO
             CALL s_quad( R2V,R2Vbar )
             srcvf = 4.D0*pi*nbub(i)*R2Vbar
             divum = ( velhf(i)-velhf(i-1) )/ds(i)
             pldot = qval(1)%f(i)*sound(i)**2*( srcvf-divum )
             rhsrc(Nveul)%f(i) = srcvf
             DO ir = 1,NR0
                iR0 = ir
                ra = bub(1,ir)%f(i)
                vv = bub(2,ir)%f(i)
                pb = bub(3,ir)%f(i)
                mv = bub(4,ir)%f(i)
                CALL s_bwproperty( pb )
                rhsrc(ibub(1,ir))%f(i) = qval(ibub(2,ir))%f(i)
                vflux = f_vflux( ra,vv,mv )
                rhsrc(ibub(4,ir))%f(i) = nbub(i)*vflux*4.D0*pi*ra**2
                pbdot = f_bpres_dot( vflux,ra,vv,pb,mv )
                rhsrc(ibub(3,ir))%f(i) = nbub(i)*pbdot
                rhsrc(ibub(2,ir))%f(i) = nbub(i) &
                                * f_gilmore( ra,vv,pb,pbdot,pres(i),pldot )
             END DO
          END DO
       ELSE IF ( polytropic=='n'.AND.vapor=='n' ) THEN
          DO i = 1,N1
             DO ir = 1,NR0
                R2V(ir) = bub(1,ir)%f(i)**2*bub(2,ir)%f(i)
             END DO
             CALL s_quad( R2V,R2Vbar )
             srcvf = 4.D0*pi*nbub(i)*R2Vbar
             divum = ( velhf(i)-velhf(i-1) )/ds(i)
             pldot = qval(1)%f(i)*sound(i)**2*( srcvf-divum )
             rhsrc(Nveul)%f(i) = srcvf
             DO ir = 1,NR0
                iR0 = ir
                ra = bub(1,ir)%f(i)
                vv = bub(2,ir)%f(i)
                pb = bub(3,ir)%f(i)
                CALL s_bwproperty( pb )
                rhsrc(ibub(1,ir))%f(i) = qval(ibub(2,ir))%f(i)
                vflux = f_vflux( ra,vv,0.D0 )
                pbdot = f_bpres_dot( vflux,ra,vv,pb,0.D0 )
                rhsrc(ibub(3,ir))%f(i) = nbub(i)*pbdot
                rhsrc(ibub(2,ir))%f(i) = nbub(i) &
                                * f_gilmore( ra,vv,pb,pbdot,pres(i),pldot )
             END DO
          END DO
       ELSE IF ( polytropic=='y' ) THEN
          DO i = 1,N1
             DO ir = 1,NR0
                R2V(ir) = bub(1,ir)%f(i)**2*bub(2,ir)%f(i)
             END DO
             CALL s_quad( R2V,R2Vbar )
             srcvf = 4.D0*pi*nbub(i)*R2Vbar
             divum = ( velhf(i)-velhf(i-1) )/ds(i)
             pldot = qval(1)%f(i)*sound(i)**2*( srcvf-divum )
             rhsrc(Nveul)%f(i) = srcvf
             DO ir = 1,NR0
                iR0 = ir
                ra = bub(1,ir)%f(i)
                vv = bub(2,ir)%f(i)
                rhsrc(ibub(1,ir))%f(i) = qval(ibub(2,ir))%f(i)
                rhsrc(ibub(2,ir))%f(i) = nbub(i) &
                               * f_gilmore( ra,vv,0.D0,0.D0,pres(i),pldot )
             END DO
          END DO
       END IF

    ELSE

       DO i = 1,N1
          DO ir = 1,NR0
             R2V(ir) = bub(1,ir)%f(i)**2*bub(2,ir)%f(i)
          END DO
          CALL s_quad( R2V,R2Vbar )
          srcvf = 4.D0*pi*nbub(i)*R2Vbar
          divum = ( velhf(i)-velhf(i-1) )/ds(i)
          dpldt(i) = qval(1)%f(i)*sound(i)**2*( srcvf-divum )
       END DO

    END IF

  END SUBROUTINE s_compsrc

  !========================================================================

  SUBROUTINE s_weno

    INTEGER :: i
    INTEGER :: m
    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(0:2) :: beta
    REAL(KIND(0.D0)), DIMENSION(0:2) :: alpha
    REAL(KIND(0.D0)), DIMENSION(0:2) :: poly
    REAL(KIND(0.D0)), DIMENSION(0:2) :: d
    REAL(KIND(0.D0)), DIMENSION(-2:1) :: dvd
    REAL(KIND(0.D0)) :: tmpmin
    REAL(KIND(0.D0)) :: tmpmax
    REAL(KIND(0.D0)), DIMENSION(0:2) :: curv
    REAL(KIND(0.D0)) :: ulc
    REAL(KIND(0.D0)) :: umd
    REAL(KIND(0.D0)) :: uul
    REAL(KIND(0.D0)) :: minmod
    REAL(KIND(0.D0)) :: diff1
    REAL(KIND(0.D0)) :: diff2
    REAL(KIND(0.D0)) :: wenosum

    ! characteristic decomposition
    IF ( chardecomp=='y' ) THEN

       DO m = -wenonum,wenonum
          !!! left-going wave at i-1/2 (left celledge) !!!
          CALL s_multiply_ll( inweno,m,-1 )
          linw%w(m) = v_out
          !!! right-going wave at i+1/2 (right celledge) !!!
          CALL s_multiply_ll( inweno,m,0 )
          rinw%w(m) = v_out
       END DO

    ELSE ! component-wise

       IF ( bound%beg=='reflect'.AND.mpi_rank==0 ) THEN
          DO m = -wenonum,0
             DO iv = 1,Nv
                DO i = 1,N1
                   linw(iv)%w(m)%f(i) = inweno(iv)%f(i+m)
                END DO
             END DO
          END DO
       ELSE
          DO m = -wenonum,0
             DO iv = 1,Nv
                DO i = 1,N1
                   linw(iv)%w(m)%f(i) = inweno(iv)%f(ip(i,m))
                END DO
             END DO
          END DO
       END IF
       IF ( bound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN
          DO m = 0,wenonum
             DO iv = 1,Nv
                DO i = 1,N1
                   linw(iv)%w(m)%f(i) = inweno(iv)%f(i+m)
                END DO
             END DO
          END DO
       ELSE
          DO m = 0,wenonum
             DO iv = 1,Nv
                DO i = 0,N1
                   linw(iv)%w(m)%f(i) = inweno(iv)%f(ip(i,m))
                END DO
             END DO
          END DO
       END IF
       rinw = linw

    END IF

    ! build interpolation polynomial with WENO reconstruction
    IF ( wenoord==1 ) THEN

       DO iv = 1,Nv
          DO i = 1,N1
             rweno(iv)%f(i) = rinw(iv)%w(0)%f(i)
             lweno(iv)%f(i) = linw(iv)%w(0)%f(i)
          END DO
       END DO

    ELSE IF ( wenoord==3 ) THEN

       DO iv = 1,Nv
          DO i = 1,N1

             !!! right-going wave at i+1/2 ("L" state) !!!
             ! differences
             dvd(-1) = rinw(iv)%w(0)%f(i) - rinw(iv)%w(-1)%f(i)
             dvd(0) = rinw(iv)%w(1)%f(i) - rinw(iv)%w(0)%f(i)
             ! smoothness indicators
             beta(0) = w3beta(0,i)*dvd(-1)*dvd(-1) + eps
             beta(1) = w3beta(1,i)*dvd(0)*dvd(0) + eps
             ! weights
             d(0) = w3dr(0,i)
             d(1) = w3dr(1,i)
             DO m = 0,1
                alpha(m) = d(m)/( beta(m)*beta(m) )
             END DO
             wenosum = alpha(0) + alpha(1)
             alpha = alpha/wenosum
             ! polynomials
             poly(0) = rinw(iv)%w(0)%f(i) + w3polyr(0,i)*dvd(-1)
             poly(1) = rinw(iv)%w(0)%f(i) + w3polyr(1,i)*dvd(0)
             ! WENO value
             rweno(iv)%f(i) = alpha(0)*poly(0) + alpha(1)*poly(1)

             !!! left-going wave at i-1/2 ("R" state) !!!
             ! differences
             dvd(-1) = linw(iv)%w(0)%f(i) - linw(iv)%w(-1)%f(i)
             dvd(0) = linw(iv)%w(1)%f(i) - linw(iv)%w(0)%f(i)
             ! smoothness indicators
             beta(0) = w3beta(0,i)*dvd(-1)*dvd(-1) + eps
             beta(1) = w3beta(1,i)*dvd(0)*dvd(0) + eps
             ! weights
             d(0) = w3dl(0,i)
             d(1) = w3dl(1,i)
             DO m = 0,1
                alpha(m) = d(m)/( beta(m)*beta(m) )
             END DO
             wenosum = alpha(0) + alpha(1)
             alpha = alpha/wenosum
             ! polynomials
             poly(0) = linw(iv)%w(0)%f(i) + w3polyl(0,i)*dvd(-1)
             poly(1) = linw(iv)%w(0)%f(i) + w3polyl(1,i)*dvd(0)
             ! WENO value
             lweno(iv)%f(i) = alpha(0)*poly(0) + alpha(1)*poly(1)

          END DO
       END DO

    ELSE IF ( wenoord==5 ) THEN

       DO iv = 1,Nv
          DO i = 1,N1 ! within computational domain

             !!! right-going wave at i+1/2 ("L" state) !!!
             ! differences
             DO m = -2,1
                dvd(m) = rinw(iv)%w(m+1)%f(i) - rinw(iv)%w(m)%f(i)
             END DO
             ! smoothness indicators
             beta(0) = eps + betaw0(0,i)*dvd(-1)*dvd(-1) &
                     + betaw0(1,i)*dvd(-1)*dvd(-2) &
                     + betaw0(2,i)*dvd(-2)*dvd(-2)
             beta(1) = eps + betaw1(0,i)*dvd(0)*dvd(0) &
                     + betaw1(1,i)*dvd(0)*dvd(-1) &
                     + betaw1(2,i)*dvd(-1)*dvd(-1)
             beta(2) = eps + betaw2(0,i)*dvd(1)*dvd(1) &
                     + betaw2(1,i)*dvd(1)*dvd(0) &
                     + betaw2(2,i)*dvd(0)*dvd(0)
             ! weights
             DO m = 0,2
                d(m) = dwenor(m,i)
             END DO
             DO m = 0,2
                alpha(m) = d(m)/( beta(m)*beta(m) )
             END DO
             wenosum = alpha(0) + alpha(1) + alpha(2)
             alpha = alpha/wenosum
             ! polynomials
             poly(0) = rinw(iv)%w(0)%f(i) + dvd(-1)*polyr0(0,i) &
                     + dvd(-2)*polyr0(1,i)
             poly(1) = rinw(iv)%w(0)%f(i) + dvd(0)*polyr1(0,i) &
                     + dvd(-1)*polyr1(1,i)
             poly(2) = rinw(iv)%w(0)%f(i) + dvd(1)*polyr2(0,i) &
                     + dvd(0)*polyr2(1,i)
             ! WENO value
             rweno(iv)%f(i) = alpha(0)*poly(0) + alpha(1)*poly(1) &
                            + alpha(2)*poly(2)
             ! for nonreflective BCs
             IF ( ((bound%beg=='nonreflect'.AND.i==2).AND.mpi_rank==0) &
                  .OR. &
                  ((bound%end=='nonreflect'.AND.i==N1-1).AND. &
                  mpi_rank==mpi_size-1) ) THEN
                ! WENO3
                dvd(-1) = rinw(iv)%w(0)%f(i) - rinw(iv)%w(-1)%f(i)
                dvd(0) = rinw(iv)%w(1)%f(i) - rinw(iv)%w(0)%f(i)
                beta(0) = w3beta(0,i)*dvd(-1)*dvd(-1) + eps
                beta(1) = w3beta(1,i)*dvd(0)*dvd(0) + eps
                d(0) = w3dr(0,i)
                d(1) = w3dr(1,i)
                DO m = 0,1
                   alpha(m) = d(m)/( beta(m)*beta(m) )
                END DO
                wenosum = alpha(0) + alpha(1)
                alpha = alpha/wenosum
                poly(0) = rinw(iv)%w(0)%f(i) + w3polyr(0,i)*dvd(-1)
                poly(1) = rinw(iv)%w(0)%f(i) + w3polyr(1,i)*dvd(0)
                rweno(iv)%f(i) = alpha(0)*poly(0) + alpha(1)*poly(1)
             END IF
             ! monotonicity preserving (Balsara 2000)
             IF ( mpweno=='y' ) THEN
                curv(0) = dvd(-1) - dvd(-2)
                curv(1) = dvd(0) - dvd(-1)
                minmod = 0.5D0 &
                       * ( SIGN(1.D0,curv(0))+SIGN(1.D0,curv(1)) ) &
                       * MIN( DABS(curv(0)),DABS(curv(1)) )
                ulc = rinw(iv)%w(0)%f(i) + 0.5D0*dvd(-1) &
                    + four3rd*minmod
                curv(0) = dvd(0) - dvd(-1)
                curv(1) = dvd(1) - dvd(0)
                minmod = 0.5D0 &
                       * ( SIGN(1.D0,curv(0))+SIGN(1.D0,curv(1)) ) &
                       * MIN( DABS(curv(0)),DABS(curv(1)) )
                umd = 0.5D0*( rinw(iv)%w(0)%f(i)+rinw(iv)%w(1)%f(i) &
                    - minmod )
                uul = rinw(iv)%w(0)%f(i) + 2.D0*dvd(-1)
                tmpmin = MAX( &
                         MIN(rinw(iv)%w(0)%f(i),rinw(iv)%w(1)%f(i),umd), &
                         MIN(rinw(iv)%w(0)%f(i),uul,ulc) )
                tmpmax = MIN( &
                         MAX(rinw(iv)%w(0)%f(i),rinw(iv)%w(1)%f(i),umd), &
                         MAX(rinw(iv)%w(0)%f(i),uul,ulc) )
                diff1 = tmpmin - rweno(iv)%f(i)
                diff2 = tmpmax - rweno(iv)%f(i)
                rweno(iv)%f(i) = rweno(iv)%f(i) + 0.5D0*( &
                                 SIGN(1.D0,diff1)+SIGN(1.D0,diff2) ) &
                               * MIN( DABS(diff1),DABS(diff2) )
             END IF

             !!! left-going wave at i-1/2 ("R" state) !!!
             ! differences
             DO m = -2,1
                dvd(m) = linw(iv)%w(m+1)%f(i) - linw(iv)%w(m)%f(i)
             END DO
             ! smoothness indicatiors
             beta(0) = eps + betaw0(0,i)*dvd(-1)*dvd(-1) &
                     + betaw0(1,i)*dvd(-1)*dvd(-2) &
                     + betaw0(2,i)*dvd(-2)*dvd(-2)
             beta(1) = eps + betaw1(0,i) * dvd(0) * dvd(0) &
                     + betaw1(1,i)*dvd(0)*dvd(-1) &
                     + betaw1(2,i)*dvd(-1)*dvd(-1)
             beta(2) = eps + betaw2(0,i)*dvd(1)*dvd(1) &
                     + betaw2(1,i)*dvd(1)*dvd(0) &
                     + betaw2(2,i)*dvd(0)*dvd(0)
             ! weights
             DO m = 0,2
                d(m) = dwenol(m,i)
             END DO
             DO m = 0,2
                alpha(m) = d(m)/( beta(m)*beta(m) )
             END DO
             wenosum = alpha(0) + alpha(1) + alpha(2)
             alpha = alpha/wenosum
             ! polynomials
             poly(0) = linw(iv)%w(0)%f(i) + dvd(-1)*polyl0(0,i) &
                     + dvd(-2) * polyl0(1,i)
             poly(1) = linw(iv)%w(0)%f(i) + dvd(0)*polyl1(0,i) &
                     + dvd(-1)*polyl1(1,i)
             poly(2) = linw(iv)%w(0)%f(i) + dvd(1)*polyl2(0,i) &
                     + dvd(0)*polyl2(1,i)
             ! WENO value
             lweno(iv)%f(i) = alpha(0)*poly(0) + alpha(1)*poly(1) &
                            + alpha(2)*poly(2)
             ! for nonreflective BCs
             IF ( ((bound%beg=='nonreflect'.AND.i==2).AND.mpi_rank==0) &
                  .OR. &
                  ((bound%end=='nonreflect'.AND.i==N1-1).AND. &
                  mpi_rank==mpi_size-1) ) THEN
                ! WENO3
                dvd(-1) = linw(iv)%w(0)%f(i) - linw(iv)%w(-1)%f(i)
                dvd(0) = linw(iv)%w(1)%f(i) - linw(iv)%w(0)%f(i)
                beta(0) = w3beta(0,i)*dvd(-1)*dvd(-1) + eps
                beta(1) = w3beta(1,i)*dvd(0)*dvd(0) + eps
                d(0) = w3dl(0,i)
                d(1) = w3dl(1,i)
                DO m = 0,1
                   alpha(m) = d(m)/( beta(m)*beta(m) )
                END DO
                wenosum = alpha(0) + alpha(1)
                alpha = alpha/wenosum
                poly(0) = linw(iv)%w(0)%f(i) + w3polyl(0,i)*dvd(-1)
                poly(1) = linw(iv)%w(0)%f(i) + w3polyl(1,i)*dvd(0)
                lweno(iv)%f(i) = alpha(0)*poly(0) + alpha(1)*poly(1)
             END IF
             ! monotonicity preserving (Balsara 2000)
             IF ( mpweno=='y' ) THEN
                curv(0) = dvd(0) - dvd(-1)
                curv(1) = dvd(1) - dvd(0)
                minmod = 0.5D0*( SIGN(1.D0,curv(0))+SIGN(1.D0,curv(1)) ) &
                       * MIN( DABS(curv(0)),DABS(curv(1)) )
                ulc = linw(iv)%w(0)%f(i) - 0.5D0*dvd(0) + four3rd*minmod
                curv(0) = dvd(-1) - dvd(-2)
                curv(1) = dvd(0) - dvd(-1)
                minmod = 0.5D0*( SIGN(1.D0,curv(0))+SIGN(1.D0,curv(1)) ) &
                       * MIN( DABS(curv(0)),DABS(curv(1)) )
                umd = 0.5D0*( linw(iv)%w(0)%f(i)+linw(iv)%w(-1)%f(i) &
                    - minmod )
                uul = linw(iv)%w(0)%f(i) - 2.D0*dvd(0)
                tmpmin = MAX( &
                         MIN(linw(iv)%w(0)%f(i),linw(iv)%w(-1)%f(i),umd), &
                         MIN(linw(iv)%w(0)%f(i),uul,ulc) )
                tmpmax = MIN( &
                         MAX(linw(iv)%w(0)%f(i),linw(iv)%w(-1)%f(i),umd), &
                         MAX(linw(iv)%w(0)%f(i),uul,ulc) )
                diff1 = tmpmin - lweno(iv)%f(i)
                diff2 = tmpmax - lweno(iv)%f(i)
                lweno(iv)%f(i) = lweno(iv)%f(i) + 0.5D0*( &
                                 SIGN(1.D0,diff1)+SIGN(1.D0,diff2) ) &
                               * MIN( DABS(diff1),DABS(diff2) )
             END IF

          END DO
       END DO

    END IF

  END SUBROUTINE s_weno

  !========================================================================

  SUBROUTINE s_approx_riemann

    INTEGER :: i
    INTEGER :: iv

    ! compute conserved variables & flux at cell edges
    ! reflective BCs implmented
    CALL s_celledgevalue

    ! HLLC Riemann solvers (returns nfx(:)%f)
    CALL s_hllc

    ! compute wall pressure
    duwall = 0.D0
    IF ( xincoming%beg=='freeplate'.AND.mpi_rank==0 ) CALL s_wallpres

    ! modify nfx(:)%f(1) or nfx(:)%f(N1-1) for nonreflective BCs
    IF ( wenoord/=1 ) THEN
       IF ( bound%beg=='nonreflect'.AND.mpi_rank==0 ) THEN
          CALL s_modifynfx_beg
       END IF
       IF ( bound%end=='nonreflect'.AND.mpi_rank==mpi_size-1 ) THEN
          CALL s_modifynfx_end
       END IF
    END IF

    ! RHS for the hyperbolic parts
    DO iv = 1,Nv
       DO i = 1,N1 ! within computational domain
          rhsfx(iv)%f(i) = ( nfx(iv)%f(ip(i,-1))-nfx(iv)%f(i) )/ds(i)
       END DO
    END DO

    ! modify rhsfx(:)%f(1) or rhsfx(:)%f(N1) for Thompson BCs
    IF ( bound%beg=='nonreflect'.AND.mpi_rank==0 ) THEN
       CALL s_thompsonbc_beg
    END IF
    IF ( bound%end=='nonreflect'.AND.mpi_rank==mpi_size-1 ) THEN
       CALL s_thompsonbc_end
    END IF

  END SUBROUTINE s_approx_riemann

  !========================================================================

  SUBROUTINE s_hllc

    INTEGER :: i
    INTEGER :: iv
    REAL(KIND(0.D0)) :: tmp, tmpm
    REAL(KIND(0.D0)) :: Sminus, Splus

    ! returns sl, sr & sstar at celledges
    CALL s_hllwavespeed

    ! intermediate state
    DO i = Nbeg,N1

       !!! "R*" state !!!
       tmp = ( sr(i)-rvel(ip(i,1)) )/( sr(i)-sstar(i) )
       rqstar(1)%f(i) = tmp*rqval(1)%f(ip(i,1))
       rqstar(dir1)%f(i) = rqstar(1)%f(i)*sstar(i)
       DO iv = Nveul,Nv
          rqstar(iv)%f(i) = tmp*rqval(iv)%f(ip(i,1))
       END DO
       !!! "L*" state !!!
       tmp = ( sl(i)-lvel(i) )/( sl(i)-sstar(i) )
       lqstar(1)%f(i) = tmp*lqval(1)%f(i)
       lqstar(dir1)%f(i) = lqstar(1)%f(i)*sstar(i)
       DO iv = Nveul,Nv
          lqstar(iv)%f(i) = tmp*lqval(iv)%f(i)
       END DO
       !!! HLLC flux !!!
       tmpm = SIGN( 0.5D0,sstar(i) )
       Sminus = MIN( 0.D0,sl(i) )
       Splus = MAX( 0.D0,sr(i) )
       DO iv = 1,Nv
          nfx(iv)%f(i) = ( 0.5D0+tmpm )*( lflux(iv)%f(i) &
                       + Sminus*(lqstar(iv)%f(i)-lqval(iv)%f(i)) ) &
                       + ( 0.5D0-tmpm )*( rflux(iv)%f(ip(i,1)) &
                       + Splus*(rqstar(iv)%f(i)-rqval(iv)%f(ip(i,1))) )
       END DO
!       IF ( nfx(1)%f(i)/=nfx(1)%f(i) ) THEN
!          PRINT*, 'NaN in nfx(1)%f(i) at (i,it,mpi_rank): ', i, it, mpi_rank
!          PRINT*, 'lqval(iv)%f(i) =', (lqval(iv)%f(i),iv=1,1)
!          PRINT*, 'lqstar(iv)%f(i) =', (lqstar(iv)%f(i),iv=1,1)
!          PRINT*, 'rqstar(iv)%f(i) =', (rqstar(iv)%f(i),iv=1,1)
!          PRINT*, 'rqval(iv)%f(i+1) =', (rqval(iv)%f(ip(i,1)),iv=1,1)
!       END IF

    END DO

  END SUBROUTINE s_hllc

  !========================================================================

  SUBROUTINE s_wallpres

    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(Nv) :: qtmp
    REAL(KIND(0.D0)) :: tmpm
    REAL(KIND(0.D0)) :: tmpl
    REAL(KIND(0.D0)) :: tmpr

    ! choose one state out of four states
    ! IF ( sstar(0)>=0.D0.AND.sl(0)>=0.D0 ) qtmp = lqval
    ! IF ( sstar(0)>=0.D0.AND.sl(0)<0.D0 ) qtmp = lqstar
    ! IF ( sstar(0)<0.D0.AND.sr(0)>=0.D0 ) qtmp = rqstar
    ! IF ( sstar(0)<0.D0.AND.sr(0)<0.D0 ) qtmp = rqval
    tmpm = SIGN( 0.5D0,sstar(0) )
    tmpl = SIGN( 0.5D0,sl(0) )
    tmpr = SIGN( 0.5D0,sr(0) )
    DO iv = 1,Nv
       qtmp(iv) = ( 0.5D0+tmpm )*( lqstar(iv)%f(0) &
                + (0.5D0+tmpl)*(lqval(iv)%f(0)-lqstar(iv)%f(0)) ) &
                + ( 0.5D0-tmpm )*( rqstar(iv)%f(0) &
                + (0.5D0-tmpr)*(rqval(iv)%f(ip(0,1))-rqstar(iv)%f(0)) )
    END DO

    ! wall pressure
    pwall = qtmp(1)/( 1.D0-qtmp(Nveul) )
    pwall = ( pl0+B_tait )*pwall**n_tait - B_tait
    ! wall velocity increment
    duwall = mp_inv*( pl0-pwall )

  END SUBROUTINE s_wallpres

  !========================================================================

  SUBROUTINE s_hllc_src

    INTEGER :: i
    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(Nv) :: qtmp
    REAL(KIND(0.D0)) :: tmpm
    REAL(KIND(0.D0)) :: tmpl
    REAL(KIND(0.D0)) :: tmpr

    DO i = Nbeg,N1

       ! choose one state out of four states
       ! IF ( sstar(i)>=0.D0.AND.sl(i)>=0.D0 ) qtmp = lqval
       ! IF ( sstar(i)>=0.D0.AND.sl(i)<0.D0 ) qtmp = lqstar
       ! IF ( sstar(i)<0.D0.AND.sr(i)>=0.D0 ) qtmp = rqstar
       ! IF ( sstar(i)<0.D0.AND.sr(i)<0.D0 ) qtmp = rqval
       tmpm = SIGN( 0.5D0,sstar(i) )
       tmpl = SIGN( 0.5D0,sl(i) )
       tmpr = SIGN( 0.5D0,sr(i) )
       DO iv = 1,Nv
          qtmp(iv) = ( 0.5D0+tmpm )*( lqstar(iv)%f(i) &
                   + (0.5D0+tmpl)*(lqval(iv)%f(i)-lqstar(iv)%f(i)) ) &
                   + ( 0.5D0-tmpm )*( rqstar(iv)%f(i) &
                   + (0.5D0-tmpr)*(rqval(iv)%f(ip(i,1))-rqstar(iv)%f(i)) )
       END DO
       ! numerical source flux in momentum equation
       nsfx(i) = f_nsfx( qtmp )
!       IF ( nsfx(i)/=nsfx(i) ) THEN
!          PRINT*, 'NaN in pltilde at (i,it,mpi_rank): ', i, it, mpi_rank
!          PRINT*, 'lqval(iv)%f(i) =', (lqval(iv)%f(i),iv=1,1)
!          PRINT*, 'lqstar(iv)%f(i) =', (lqstar(iv)%f(i),iv=1,1)
!          PRINT*, 'rqstar(iv)%f(i) =', (rqstar(iv)%f(i),iv=1,1)
!          PRINT*, 'rqval(iv)%f(i+1) =', (rqval(iv)%f(ip(i,1)),iv=1,1)
!          PRINT*, 'Computation stopped.'
!          STOP
!       END IF
       ! cell-edge velocity
       velhf(i) = qtmp(dir1)/qtmp(1)

    END DO

  END SUBROUTINE s_hllc_src

  !========================================================================

  FUNCTION f_nsfx( q_in )

    REAL(KIND(0.D0)), DIMENSION(Nv), INTENT(IN) :: q_in
    REAL(KIND(0.D0)) :: f_nsfx
    INTEGER :: ir
    REAL(KIND(0.D0)) :: pltmp
    REAL(KIND(0.D0)) :: ntmp
    REAL(KIND(0.D0)) :: R3av
    REAL(KIND(0.D0)) :: R3pbwav
    REAL(KIND(0.D0)) :: R3V2av
    REAL(KIND(0.D0)), DIMENSION(NR0) :: ra
    REAL(KIND(0.D0)), DIMENSION(NR0) :: vv
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pbw
    REAL(KIND(0.D0)), DIMENSION(NR0) :: nRtmp
    REAL(KIND(0.D0)) :: pbtmp
    REAL(KIND(0.D0)) :: vf1

    DO ir = 1,NR0
       nRtmp(ir) = q_in(ibub(1,ir))
    END DO
    CALL s_comp_n( q_in(Nveul),nRtmp,ntmp )
    pbtmp = 0.D0
    DO ir = 1,NR0
       iR0 = ir
       ra(ir) = q_in(ibub(1,ir))/ntmp
       vv(ir) = q_in(ibub(2,ir))/ntmp
       IF ( polytropic=='n' ) pbtmp = q_in(ibub(3,ir))/ntmp
       pbw(ir) = f_bwpres_1( ra(ir),vv(ir),pbtmp )
    END DO
    vf1 = 1.D0 - q_in(Nveul)
    pltmp = ( pl0+B_tait )*( q_in(1)/vf1 )**n_tait - B_tait
    CALL s_quad( ra**3,R3av )
    CALL s_quad( ra**3*pbw,R3pbwav )
    CALL s_quad( ra**3*vv**2,R3V2av )
    ! numerical source flux in momentum equation
    f_nsfx = q_in(Nveul)*( pltmp-(R3pbwav+q_in(1)*R3V2av)/R3av )

  END FUNCTION f_nsfx

  !========================================================================

  SUBROUTINE s_modifynfx_beg

    INTEGER :: i
    INTEGER :: iv

    ! numerical flux at 3/2 (Pirozzoli 2002) used for rhsfx(:)%f(2)
    DO i = 1,4
       fxtmp(1)%f(i) = qval(dir1)%f(i)
       fxtmp(dir1)%f(i) = qval(dir1)%f(i)*vel1(i) + pres(i)
       DO iv = Nveul,Nv
          fxtmp(iv)%f(i) = qval(iv)%f(i)*vel1(i)
       END DO
    END DO
    DO iv = 1,Nv
       nfx(iv)%f(1) = ( 3.D0*fxtmp(iv)%f(1)+13.D0*fxtmp(iv)%f(2) &
                    - 5.D0*fxtmp(iv)%f(3)+fxtmp(iv)%f(4) )*twelfth
    END DO

  END SUBROUTINE s_modifynfx_beg

  !========================================================================

  SUBROUTINE s_modifynfx_end

    INTEGER :: i
    INTEGER :: iv

    ! numerical flux at N1-1/2 used for rhsfx(:)%f(N1-1)
    DO i = N1-3,N1
       fxtmp(1)%f(i) = qval(dir1)%f(i)
       fxtmp(dir1)%f(i) = qval(dir1)%f(i)*vel1(i) + pres(i)
       DO iv = Nveul,Nv
          fxtmp(iv)%f(i) = qval(iv)%f(i)*vel1(i)
       END DO
    END DO
    DO iv = 1,Nv
       nfx(iv)%f(N1-1) = ( fxtmp(iv)%f(N1-3)-5.D0*fxtmp(iv)%f(N1-2) &
                       + 13.D0*fxtmp(iv)%f(N1-1) &
                       + 3.D0*fxtmp(iv)%f(N1) )*twelfth
    END DO

  END SUBROUTINE s_modifynfx_end

  !========================================================================

  SUBROUTINE s_thompsonbc_beg

    INTEGER :: ib
    INTEGER :: ir
    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(Nv) :: diff
    REAL(KIND(0.D0)), DIMENSION(Nv) :: lch
    REAL(KIND(0.D0)), DIMENSION(Nv) :: charsptmp
    REAL(KIND(0.D0)) :: vftmp
    REAL(KIND(0.D0)) :: vftmp1
    REAL(KIND(0.D0)) :: cltmp
    REAL(KIND(0.D0)) :: utmp
    REAL(KIND(0.D0)) :: pltmp
    REAL(KIND(0.D0)) :: rhotmp
    REAL(KIND(0.D0)), DIMENSION(Nb,NR0) :: nbubtmp
    REAL(KIND(0.D0)) :: tmp1
    REAL(KIND(0.D0)) :: tmp2

    !!! one-sided differences !!!
    DO iv = 1,Nv
       CALL s_onesidebc( qval(iv)%f,'beg' )
       diff(iv) = oneside
    END DO

    !!! Substitution of values at i=1 !!!
    rhotmp = qval(1)%f(1)
    utmp = vel1(1)
    cltmp = sound(1)
    pltmp = pres(1)
    vftmp = qval(Nveul)%f(1)
    vftmp1 = 1.D0 - qval(Nveul)%f(1)
    DO ib = 1,Nb
       DO ir = 1,NR0
          nbubtmp(ib,ir) = qval(ibub(ib,ir))%f(1)
       END DO
    END DO
    DO iv = 1,Nv
       charsptmp(iv) = charsp(iv)%f(1)
    END DO

    !!! modification on sonic speed !!!
    IF ( modifycl=='y' ) THEN
       IF ( thermal=='adiabat' ) THEN
          tmp1 = vftmp1*gamma_b*pltmp
          tmp2 = vftmp*n_tait*( pltmp+B_tait )
          cltmp = cltmp*DSQRT( tmp1/(tmp1+tmp2) )
          charsptmp(1) = utmp - cltmp
          charsptmp(2) = utmp + cltmp
       ELSE
          tmp1 = vftmp1*pltmp
          tmp2 = vftmp*n_tait*( pltmp+B_tait )
          cltmp = cltmp*DSQRT( tmp1/(tmp1+tmp2) )
          charsptmp(1) = utmp - cltmp
          charsptmp(2) = utmp + cltmp
       END IF
    END IF

    !!! Thompson-type BCs (1987) !!!
    lch = 0.D0
    ! L1
    IF ( charsptmp(1)<0.D0 ) THEN
       tmp2 = 0.5D0/rhotmp/cltmp
       tmp1 = ( utmp+vftmp1*cltmp )*tmp2
       lch(1) = ( tmp1*diff(1)-tmp2*diff(dir1)+0.5D0*diff(Nveul) ) &
              * charsptmp(1)
    END IF
    ! L2
    IF ( charsptmp(2)<0.D0 ) THEN
       tmp2 = 0.5D0/rhotmp/cltmp
       tmp1 = ( -utmp+vftmp1*cltmp )*tmp2
       lch(2) = ( tmp1*diff(1)+tmp2*diff(dir1)+0.5D0*diff(Nveul) ) &
              * charsptmp(2)
    END IF
    ! L3
    IF ( charsptmp(3)<0.D0 ) THEN
       tmp1 = vftmp/rhotmp
       lch(3) = charsptmp(3)*( tmp1*diff(1)-diff(Nveul) )
    END IF
    ! L_{nR}, L_{n\dot{R}}, L_{np_b}, L_{nm_v}
    DO ib = 1,Nb
       DO ir = 1,NR0
          IF ( charsptmp(ibub(ib,ir))<0.D0 ) THEN
             tmp2 = nbubtmp(ib,ir)
             tmp1 = -vftmp1*tmp2/rhotmp
             lch(ibub(ib,ir)) = charsptmp(ibub(ib,ir)) &
                        *( tmp1*diff(1)-tmp2*diff(3)+diff(ibub(ib,ir)) )
          END IF
       END DO
    END DO

    !!! RHS for i=1 !!!
    tmp1 = lch(1) + lch(2) + lch(3)
    rhsfx(1)%f(1) = -rhotmp*tmp1
    rhsfx(dir1)%f(1) = rhotmp*( cltmp*(lch(1)-lch(2))-utmp*tmp1 )
    rhsfx(Nveul)%f(1) = lch(3) - vftmp*tmp1
    DO ib = 1,Nb
       DO ir = 1,NR0
          rhsfx(ibub(ib,ir))%f(1) = -nbubtmp(ib,ir)*( lch(1)+lch(2) ) &
                                  - lch(ibub(ib,ir))
       END DO
    END DO

  END SUBROUTINE s_thompsonbc_beg

  !========================================================================

  SUBROUTINE s_thompsonbc_end

    INTEGER :: ib
    INTEGER :: ir
    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(Nv) :: diff
    REAL(KIND(0.D0)), DIMENSION(Nv) :: lch
    REAL(KIND(0.D0)), DIMENSION(Nv) :: charsptmp
    REAL(KIND(0.D0)) :: vftmp
    REAL(KIND(0.D0)) :: vftmp1
    REAL(KIND(0.D0)) :: cltmp
    REAL(KIND(0.D0)) :: utmp
    REAL(KIND(0.D0)) :: pltmp
    REAL(KIND(0.D0)) :: rhotmp
    REAL(KIND(0.D0)), DIMENSION(Nb,NR0) :: nbubtmp
    REAL(KIND(0.D0)) :: tmp1
    REAL(KIND(0.D0)) :: tmp2
!    REAL(KIND(0.D0)) :: dpdx
!    REAL(KIND(0.D0)) :: drhodx
!    REAL(KIND(0.D0)) :: dudx
!    REAL(KIND(0.D0)) :: drhoudx

    !!! one-sided differences !!!
    DO iv = 1,Nv
       CALL s_onesidebc( qval(iv)%f,'end' )
       diff(iv) = oneside
    END DO

    !!! Substitution of values at i=N1 !!!
    rhotmp = qval(1)%f(N1)
    utmp = vel1(N1)
    cltmp = sound(N1)
    pltmp = pres(N1)
    vftmp = qval(Nveul)%f(N1)
    vftmp1 = 1.D0 - qval(Nveul)%f(N1)
    DO ib = 1,Nb
       DO ir = 1,NR0
          nbubtmp(ib,ir) = qval(ibub(ib,ir))%f(N1)
       END DO
    END DO
    DO iv = 1,Nv
       charsptmp(iv) = charsp(iv)%f(N1)
    END DO

    !!! modification on sonic speed !!!
    IF ( modifycl=='y' ) THEN
       IF ( thermal=='adiabat' ) THEN
          tmp1 = vftmp1*gamma_b*pltmp
          tmp2 = vftmp*n_tait*( pltmp+B_tait )
          cltmp = cltmp*DSQRT( tmp1/(tmp1+tmp2) )
          charsptmp(1) = utmp - cltmp
          charsptmp(2) = utmp + cltmp
       ELSE
          tmp1 = vftmp1*pltmp
          tmp2 = vftmp*n_tait*( pltmp+B_tait )
          cltmp = cltmp*DSQRT( tmp1/(tmp1+tmp2) )
          charsptmp(1) = utmp - cltmp
          charsptmp(2) = utmp + cltmp
       END IF
    END IF

    !!! incoming wave !!!
    ! Assume very small void fraction and
    ! acoustic relations of pure liquid holds.
!    IF ( xincoming%end=='undex' ) THEN
!       dpdx = ( pl0-pltmp )/xdecay
!       drhodx = dpdx/cl0/cl0
!       dudx = -dpdx/cl0
!       drhoudx = dudx + utmp*drhodx
!    END IF

    !!! Thompson-type BCs (1987) !!!
    lch = 0.D0
    ! L1
    IF ( charsptmp(1)>0.D0 ) THEN ! outgoing
       tmp2 = 0.5D0/rhotmp/cltmp
       tmp1 = ( utmp+vftmp1*cltmp )*tmp2
       lch(1) = ( tmp1*diff(1)-tmp2*diff(dir1)+0.5D0*diff(Nveul) ) &
              * charsptmp(1)
!    ELSE ! incoming
!       IF ( xincoming%end=='undex' ) THEN
!          tmp2 = 0.5D0/rhotmp/cltmp
!          tmp1 = ( utmp+vftmp1*cltmp )*tmp2
!          lch(1) = ( tmp1*drhodx-tmp2*drhoudx+0.5D0*diff(Nveul) ) &
!                 * charsptmp(1)
!       END IF
    END IF
    ! L2
    IF ( charsptmp(2)>0.D0 ) THEN ! outgoing
       tmp2 = 0.5D0/rhotmp/cltmp
       tmp1 = ( -utmp+vftmp1*cltmp )*tmp2
       lch(2) = ( tmp1*diff(1)+tmp2*diff(dir1)+0.5D0*diff(Nveul) ) &
              * charsptmp(2)
!    ELSE ! incoming
!       IF ( xincoming%end=='undex' ) THEN
!          tmp2 = 0.5D0/rhotmp/cltmp
!          tmp1 = ( -utmp+vftmp1*cltmp )*tmp2
!          lch(2) = ( tmp1*drhodx+tmp2*drhoudx+0.5D0*diff(Nveul) ) &
!                 * charsptmp(2)
!       END IF
    END IF
    ! L3
    IF ( charsptmp(3)>0.D0 ) THEN ! outgoing
       tmp1 = vftmp/rhotmp
       lch(3) = charsptmp(3)*( tmp1*diff(1)-diff(Nveul) )
!    ELSE ! incoming
!       IF ( xincoming%end=='undex' ) THEN
!          tmp1 = vftmp/rhotmp
!          lch(3) = charsptmp(3)*( tmp1*drhodx-diff(Nveul) )
!       END IF
    END IF
    ! L_{nR}, L_{n\dot{R}}, L_{np_b}, L_{nm_v}
    DO ib = 1,Nb
       DO ir = 1,NR0
          IF ( charsptmp(ibub(ib,ir))>0.D0 ) THEN ! outgoing
             tmp2 = nbubtmp(ib,ir)
             tmp1 = -vftmp1*tmp2/rhotmp
             lch(ibub(ib,ir)) = charsptmp(ibub(ib,ir)) &
                        *( tmp1*diff(1)-tmp2*diff(3)+diff(ibub(ib,ir)) )
!          ELSE ! incoming
!             IF ( xincoming%end=='undex' ) THEN
!                tmp2 = nbubtmp(ib,ir)
!                tmp1 = -vftmp1*tmp2/rhotmp
!                lch(ibub(ib,ir)) = charsptmp(ibub(ib,ir)) &
!                        *( tmp1*drhodx-tmp2*diff(3)+diff(ibub(ib,ir)) )
!             END IF
          END IF
       END DO
    END DO

    !!! RHS for i=N1 !!!
    tmp1 = lch(1) + lch(2) + lch(3)
    rhsfx(1)%f(N1) = -rhotmp*tmp1
    rhsfx(dir1)%f(N1) = rhotmp*( cltmp*(lch(1)-lch(2))-utmp*tmp1 )
    rhsfx(Nveul)%f(N1) = lch(3) - vftmp*tmp1
    DO ib = 1,Nb
       DO ir = 1,NR0
          rhsfx(ibub(ib,ir))%f(N1) = -nbubtmp(ib,ir)*( lch(1)+lch(2) ) &
                                   - lch(ibub(ib,ir))
       END DO
    END DO

    !!! check NaN !!!
!    DO iv = 1,Nv
!       IF ( rhsfx(iv)%f(N1)/=rhsfx(iv)%f(N1) ) THEN
!          PRINT*, 'NaN at s_thompsonbc_end with (it,iv) =', it, iv
!          PRINT*, 'diff:'
!          PRINT*, ( diff(ib),ib=1,Nv )
!          PRINT*, 'lch:'
!          PRINT*, ( lch(ib),ib=1,Nv )
!          PRINT*, 'rhotmp, utmp, cltmp, pltmp, vftmp:'
!          PRINT*, rhotmp, utmp, cltmp, pltmp, vftmp
!          PRINT*, 'Computation stopped.'
!          STOP
!       END IF
!    END DO

  END SUBROUTINE s_thompsonbc_end

  !========================================================================

  SUBROUTINE s_srcbeg

    INTEGER :: i
    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(Nv) :: qtmp

    ! q_{1/2} extrapolated by 4-th order formula (Pirozzoli 2002)
    DO iv = 1,Nv
       qtmp(iv) = ( 25.D0*qval(iv)%f(1)-23.D0*qval(iv)%f(2) &
                + 13.D0*qval(iv)%f(3)-3.D0*qval(iv)%f(4) )*twelfth
    END DO
    nsfx(0) = f_nsfx( qtmp )
    velhf(0) = qtmp(dir1)/qtmp(1)

    ! q_{3/2} interpolated by 4-th order formula (Pirozzoli 2002)
    DO iv = 1,Nv
       qtmp(iv) = ( 3.D0*qval(iv)%f(1)+13.D0*qval(iv)%f(2) &
                - 5.D0*qval(iv)%f(3)+qval(iv)%f(4) )*twelfth
    END DO
    nsfx(1) = f_nsfx( qtmp )
    velhf(1) = qtmp(dir1)/qtmp(1)

  END SUBROUTINE s_srcbeg

  !========================================================================

  SUBROUTINE s_srcend

    INTEGER :: i
    INTEGER :: iv
    REAL(KIND(0.D0)), DIMENSION(Nv) :: qtmp

    ! q_{N1+1/2} extrapolated by 4-th order formula (Pirozzoli 2002)
    DO iv = 1,Nv
       qtmp(iv) = ( -3.D0*qval(iv)%f(N1-3)+13.D0*qval(iv)%f(N1-2) &
                - 23.D0*qval(iv)%f(N1-1)+25.D0*qval(iv)%f(N1) )*twelfth
    END DO
    nsfx(N1) = f_nsfx( qtmp )
    velhf(N1) = qtmp(dir1)/qtmp(1)
!    IF ( nsfx(N1)/=nsfx(N1) ) THEN
!       PRINT*, 'NaN at s_srcend with (it,iv) =', it, iv
!       DO iv = 1,Nv
!          PRINT*, 'iv, qtmp(iv) =', iv, qtmp(iv)
!       END DO
!       STOP
!    END IF

    ! q_{N1-1/2} interpolated by 4-th order formula (Pirozzoli 2002)
    DO iv = 1,Nv
       qtmp(iv) = ( qval(iv)%f(N1-3)-5.D0*qval(iv)%f(N1-2) &
                + 13.D0*qval(iv)%f(N1-1)+3.D0*qval(iv)%f(N1) )*twelfth
    END DO
    nsfx(N1-1) = f_nsfx( qtmp )
    velhf(N1-1) = qtmp(dir1)/qtmp(1)
!    IF ( nsfx(N1-1)/=nsfx(N1-1) ) THEN
!       PRINT*, 'NaN at s_srcend with (it,iv) =', it, iv
!       DO iv = 1,Nv
!          PRINT*, 'iv, qtmp(iv) =', iv, qtmp(iv)
!       END DO
!       STOP
!    END IF

  END SUBROUTINE s_srcend

  !========================================================================

  SUBROUTINE s_onesidebc( input,boundary )

    REAL(KIND(0.D0)), DIMENSION(1:N1), INTENT(IN) :: input
    CHARACTER(LEN=3), INTENT(IN) :: boundary

    IF ( boundary=='beg' ) THEN
       oneside = onesidebeg(1)*input(1) + onesidebeg(2)*input(2) &
               + onesidebeg(3)*input(3) + onesidebeg(4)*input(4) &
               + onesidebeg(5)*input(5)
    ELSE IF ( boundary=='end' ) THEN
       oneside = onesideend(1)*input(N1) + onesideend(2)*input(N1-1) &
               + onesideend(3)*input(N1-2) + onesideend(4)*input(N1-3) &
               + onesideend(5)*input(N1-4)
    END IF

  END SUBROUTINE s_onesidebc

  !========================================================================

END MODULE m_rhs
