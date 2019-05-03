!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module that assigns values used in m_rhs
!
!  Last update: January 24, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_rhsvar

  USE mpi_setup
  USE m_globalvar
  USE m_misc
  IMPLICIT NONE

  ! for direction 1 to be solved
  INTEGER :: Nbeg
  INTEGER :: Nend
  INTEGER :: Nroe
  INTEGER :: N1
  INTEGER :: dir1

  ! index matrix
  INTEGER, DIMENSION(:,:), POINTER :: ip
  ! grid spacing
  REAL(KIND(0.D0)), DIMENSION(:), POINTER :: ds

  ! WENO coefficients
  ! smoothness indicators
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: betaw0
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: betaw1
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: betaw2
  ! ideal weights
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: dwenor
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: dwenol
  ! polynomials
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: polyr0
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: polyr1
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: polyr2
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: polyl0
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: polyl1
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: polyl2
  ! WENO3
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: w3beta
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: w3dr
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: w3dl
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: w3polyr
  REAL(KIND(0.D0)), DIMENSION(:,:), POINTER :: w3polyl

  ! velocity in x direction
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE, TARGET :: xvel
  ! velocity in direction 1 to be solved
  REAL(KIND(0.D0)), DIMENSION(:), POINTER :: vel1

  ! physical quantities for cells
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: pres
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: sound
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: nbub
  ! "L" & "R" quantities
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: lpres
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: rpres
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: lvel
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: rvel

  ! characteristic speeds for cells
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: charsp
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: charsproe

  ! bubble dynamics variables for cells
  ! bub(1,:): R, bub(2,:): \dot{R}, bub(3,:): p_b, bub(4,:): m_v
  TYPE(coordinate), DIMENSION(:,:), ALLOCATABLE :: bub

  ! BCs
  TYPE(bc) :: bound
  ! used for nonreflective BCs
  REAL(KIND(0.D0)) :: oneside
  REAL(KIND(0.D0)), DIMENSION(:), POINTER :: onesidebeg
  REAL(KIND(0.D0)), DIMENSION(:), POINTER :: onesideend

  ! conserved variables
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: qval
  ! numerical flux at cell edges
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: nfx
  ! numerical source flux in momentum equation at cell edges
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: nsfx
  ! cell-edge velocity
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: velhf
  ! flux used in subroutines for nonreflective BCs
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: fxtmp
  ! rhs for each cell
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rhsfx
  ! sources for each cell
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rhsrc

  ! quantities for WENO
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: inweno
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: lweno
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rweno
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: lflux
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rflux
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: lqval
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rqval
  TYPE(wenochar), DIMENSION(:), ALLOCATABLE :: linw
  TYPE(wenochar), DIMENSION(:), ALLOCATABLE :: rinw

  ! for characteristic decomposition
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: v_out
  ! left & right eigenvectors at cell edges
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rr
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: ll

  ! HLLC
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: sr
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: sl
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: sstar
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: lqstar
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rqstar

  CONTAINS

  !========================================================================

  SUBROUTINE s_start_rhsvar

    INTEGER :: n
    INTEGER :: iv
    INTEGER :: ir
    INTEGER :: ib
    INTEGER :: ieig
    INTEGER :: istart
    INTEGER :: istop

    ! pointers for direction to be solved
    N1 = Nx
    bound = xbound
    ip => ix
    ds => dxs
    ! for WENO coefficients
    betaw0 => betax0
    betaw1 => betax1
    betaw2 => betax2
    dwenor => dwxr
    dwenol => dwxl
    polyr0 => polyrx0
    polyr1 => polyrx1
    polyr2 => polyrx2
    polyl0 => polylx0
    polyl1 => polylx1
    polyl2 => polylx2
    onesidebeg => onesidebegx
    onesideend => onesideendx
    w3beta => w3betax
    w3dr => w3drx
    w3dl => w3dlx
    w3polyr => w3polyrx
    w3polyl => w3polylx
    dir1 = 2
    IF ( bound%beg=='reflect'.AND.mpi_rank==0 ) THEN
       istart = 1 - wenonum
       Nbeg = 0
    ELSE
       istart = 1
       Nbeg = 1
    END IF
    IF ( bound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN
       istop = N1 + wenonum
       Nend = N1 + 1
       Nroe = N1 - 1
    ELSE
       istop = N1
       Nend = N1
       Nroe = N1
    END IF

    ! velocity
    ALLOCATE( xvel(1:N1) )
    vel1 => xvel

    ! cell-averaged quantities
    ALLOCATE( pres(1:N1) )
    ALLOCATE( sound(1:N1) )
    ALLOCATE( nbub(1:N1) )
    ! "L" & "R" quantities
    ALLOCATE( lpres(Nbeg:N1) )
    ALLOCATE( rpres(1:Nend) )
    ALLOCATE( lvel(Nbeg:N1) )
    ALLOCATE( rvel(1:Nend) )
    ! HLLC wavespeeds
    ALLOCATE( sl(Nbeg:N1) )
    ALLOCATE( sr(Nbeg:N1) )
    ALLOCATE( sstar(Nbeg:N1) )
    ! numerical source flux in momentum equation
    ALLOCATE( nsfx(0:N1) )
    ! cell-edge velocity
    ALLOCATE( velhf(0:N1) )

    ALLOCATE( charsp(Nv),charsproe(Nv) )
    ALLOCATE( qval(Nv) )
    ALLOCATE( nfx(Nv) )
    ALLOCATE( fxtmp(Nv) )
    ALLOCATE( rhsfx(Nv) )
    ALLOCATE( rhsrc(Nv) )
    ALLOCATE( inweno(Nv) )
    ALLOCATE( lweno(Nv) )
    ALLOCATE( rweno(Nv) )
    ALLOCATE( lflux(Nv) )
    ALLOCATE( rflux(Nv) )
    ALLOCATE( lqval(Nv) )
    ALLOCATE( rqval(Nv) )
    ALLOCATE( linw(Nv) )
    ALLOCATE( rinw(Nv) )
    ALLOCATE( v_out(Nv) )
    ALLOCATE( lqstar(Nv) )
    ALLOCATE( rqstar(Nv) )
    DO iv = 1,Nv
       ! numerical flux at cell edges
       ALLOCATE( nfx(iv)%f(Nbeg:N1) )
       ! RHS for cells
       ALLOCATE( rhsfx(iv)%f(1:N1) )
       ALLOCATE( rhsrc(iv)%f(1:N1) )
       ALLOCATE( qval(iv)%f(1:N1) )
       ALLOCATE( fxtmp(iv)%f(1:N1) )
       ALLOCATE( inweno(iv)%f(istart:istop) )
       ! note that lweno is "R" state
       ALLOCATE( lweno(iv)%f(1:N1) )
       ! note that rweno is "L" state
       ALLOCATE( rweno(iv)%f(1:N1) )
       ! eigenvalues for cells
       ALLOCATE( charsp(iv)%f(1:N1) )
       ! Roe-averaged eigenvalues at cell edges
       ALLOCATE( charsproe(iv)%f(Nbeg:N1) )
       DO n = -wenonum,wenonum
          ALLOCATE( linw(iv)%w(n)%f(1:N1) )
          ALLOCATE( rinw(iv)%w(n)%f(1:N1) )
       END DO
       ! used for local characteristic decomposition
       ALLOCATE( v_out(iv)%f(1:N1) )
       ! lqval is "L" state <== rweno
       ALLOCATE( lqval(iv)%f(Nbeg:N1) )
       ! rqval is "R" state <== lweno
       ALLOCATE( rqval(iv)%f(1:Nend) )
       ! lflux "L" state, rflux "R" state
       ALLOCATE( lflux(iv)%f(Nbeg:N1) )
       ALLOCATE( rflux(iv)%f(1:Nend) )
       ALLOCATE( lqstar(iv)%f(Nbeg:N1) )
       ALLOCATE( rqstar(iv)%f(Nbeg:N1) )
    END DO
    ALLOCATE( bub(Nb,NR0) )
    DO ib = 1,Nb
       DO ir = 1,NR0
          ! bubble dynamics variables
          ALLOCATE( bub(ib,ir)%f(1:N1) )
       END DO
    END DO

    ! right (rr) & left (ll) eigenvectors
    ALLOCATE( rr(Neig) )
    ALLOCATE( ll(Neig) )
    DO ieig = 1,Neig
       ALLOCATE( rr(ieig)%f(Nbeg:N1) )
       ALLOCATE( ll(ieig)%f(Nbeg:N1) )
    END DO
    
  END SUBROUTINE s_start_rhsvar

  !========================================================================

  SUBROUTINE s_compvar

    INTEGER :: i
    INTEGER :: ir
    INTEGER :: ib
    INTEGER :: iv
    REAL(KIND(0.D0)) :: plpB
    REAL(KIND(0.D0)) :: vf1
    REAL(KIND(0.D0)), DIMENSION(NR0) :: nRtmp

    ! velocity
    xvel = qval(2)%f/qval(1)%f

    DO i = 1,N1

       DO ir = 1,NR0
          nRtmp(ir) = qval(ibub(1,ir))%f(i)
       END DO
       CALL s_comp_n( qval(Nveul)%f(i),nRtmp,nbub(i) )
       DO ir = 1,NR0
          ! bubble-dynamic variables
          DO ib = 1,Nb
             bub(ib,ir)%f(i) = qval(ibub(ib,ir))%f(i)/nbub(i)
          END DO
       END DO
       ! averaged liquid pressure
       vf1 = 1.D0 - qval(Nveul)%f(i)
       plpB = ( pl0+B_tait )*( qval(1)%f(i)/vf1 )**n_tait
       pres(i) = plpB - B_tait
       ! \tilde{c}_l
       sound(i) = DSQRT( n_tait*plpB/vf1/qval(1)%f(i) )

    END DO

    ! characteristic speeds (eigenvalues) used for nonreflective BCs
    charsp(1)%f = vel1 - sound
    charsp(2)%f = vel1 + sound
    DO iv = Nveul,Nv
       charsp(iv)%f = vel1
    END DO

  END SUBROUTINE s_compvar

  !========================================================================

  SUBROUTINE s_roeaverage

    INTEGER :: i
    INTEGER :: iv
    INTEGER :: ir
    INTEGER :: ib
    REAL(KIND(0.D0)) :: rhoroe
    REAL(KIND(0.D0)) :: ratio
    REAL(KIND(0.D0)) :: denom
    REAL(KIND(0.D0)) :: umroe
    REAL(KIND(0.D0)) :: vfroe
    REAL(KIND(0.D0)) :: vfroe1
    REAL(KIND(0.D0)) :: nbubroe
    REAL(KIND(0.D0)) :: plpBroe
    REAL(KIND(0.D0)) :: cmroe

    DO i = 1,Nroe

       !!! Euler part !!!
       ! Roe averages
       rhoroe = DSQRT( qval(1)%f(i)*qval(1)%f(ip(i,1)) )
       ratio = DSQRT( qval(1)%f(ip(i,1))/qval(1)%f(i) )
       denom = 1.D0 + ratio
       umroe = ( vel1(i)+ratio*vel1(ip(i,1)) )/denom
       DO iv = Nveul,Nv
          charsproe(iv)%f(i) = umroe
       END DO
       vfroe = ( ratio*qval(Nveul)%f(i)+qval(Nveul)%f(ip(i,1)) )/denom
       vfroe1 = 1.D0 - vfroe
       plpBroe = ( pl0+B_tait )*( rhoroe/vfroe1 )**n_tait
       cmroe = DSQRT( n_tait*plpBroe/vfroe1/rhoroe )
       charsproe(1)%f(i) = umroe - cmroe
       charsproe(2)%f(i) = umroe + cmroe
       ! right eigenvectors at i+1/2
       rr(1)%f(i) = rhoroe
       rr(2)%f(i) = rhoroe*charsproe(1)%f(i)
       rr(3)%f(i) = rhoroe*charsproe(2)%f(i)
       rr(4)%f(i) = rhoroe*umroe
       rr(5)%f(i) = vfroe
       ! left eigenvectors at i+1/2
       ll(2)%f(i) = 0.5D0/rhoroe/cmroe
       ll(1)%f(i) = ( umroe+vfroe1*cmroe )*ll(2)%f(i)
       ll(3)%f(i) = ( -umroe+vfroe1*cmroe )*ll(2)%f(i)
       ll(4)%f(i) = vfroe/rhoroe
       ll(5)%f(i) = vfroe1/rhoroe
       !!! bubble-dynamic part !!!
       DO ib = 1,Nb
          DO ir = 1,NR0
             ! Roe averages
             nbubroe = ( ratio*qval(ibub(ib,ir))%f(i) &
                     + qval(ibub(ib,ir))%f(ip(i,1)) )/denom
             ! right eigenvectors at i+1/2
             rr(ieig(ib,ir))%f(i) = nbubroe
             ! left eigenvectors at i+1/2
             ll(ieig(ib,ir))%f(i) = nbubroe
          END DO
       END DO

    END DO

    !!! reflective BCs !!!
    ! for left boundary (i=0), moving wall BC incoporated
    IF ( bound%beg=='reflect'.AND.mpi_rank==0 ) THEN

       !!! Euler part !!!
       ! Roe averages
       rhoroe = qval(1)%f(1)
       umroe = uwall ! wall velocity
       DO iv = Nveul,Nv
          charsproe(iv)%f(0) = umroe
       END DO
       vfroe = qval(Nveul)%f(1)
       vfroe1 = 1.D0 - vfroe
       plpBroe = ( pl0+B_tait )*( rhoroe/vfroe1 )**n_tait
       cmroe = DSQRT( n_tait*plpBroe/vfroe1/rhoroe )
       charsproe(1)%f(0) = umroe - cmroe
       charsproe(2)%f(0) = umroe + cmroe
       ! right eigenvectors at 1/2
       rr(1)%f(0) = rhoroe
       rr(2)%f(0) = rhoroe*charsproe(1)%f(0)
       rr(3)%f(0) = rhoroe*charsproe(2)%f(0)
       rr(4)%f(0) = rhoroe*umroe
       rr(5)%f(0) = vfroe
       ! left eigenvectors at 1/2
       ll(2)%f(0) = 0.5D0/rhoroe/cmroe
       ll(1)%f(0) = ( umroe+vfroe1*cmroe )*ll(2)%f(0)
       ll(3)%f(0) = ( -umroe+vfroe1*cmroe )*ll(2)%f(0)
       ll(4)%f(0) = vfroe/rhoroe
       ll(5)%f(0) = vfroe1/rhoroe
       !!! bubble-dynamic part !!!
       DO ib = 1,Nb
          DO ir = 1,NR0
             ! right eigenvectors at 1/2
             rr(ieig(ib,ir))%f(0) = qval(ibub(ib,ir))%f(1)
             ! left eigenvectors at 1/2
             ll(ieig(ib,ir))%f(0) = qval(ibub(ib,ir))%f(1)
          END DO
       END DO

    END IF
    ! for right boundary (i=N1)
    IF ( bound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN

       !!! Euler part !!!
       ! Roe averages
       rhoroe = qval(1)%f(N1)
       umroe = 0.D0
       DO iv = Nveul,Nv
          charsproe(iv)%f(N1) = umroe
       END DO
       vfroe = qval(Nveul)%f(N1)
       vfroe1 = 1.D0 - vfroe
       plpBroe = ( pl0+B_tait )*( rhoroe/vfroe1 )**n_tait
       cmroe = DSQRT( n_tait*plpBroe/vfroe1/rhoroe )
       charsproe(1)%f(N1) = umroe - cmroe
       charsproe(2)%f(N1) = umroe + cmroe
       ! right eigenvectors at N1+1/2
       rr(1)%f(N1) = rhoroe
       rr(2)%f(N1) = rhoroe*charsproe(1)%f(N1)
       rr(3)%f(N1) = rhoroe*charsproe(2)%f(N1)
       rr(4)%f(N1) = rhoroe*umroe
       rr(5)%f(N1) = vfroe
       ! left eigenvectors at N1+1/2
       ll(2)%f(N1) = 0.5D0/rhoroe/cmroe
       ll(1)%f(N1) = ( umroe+vfroe1*cmroe )*ll(2)%f(N1)
       ll(3)%f(N1) = ( -umroe+vfroe1*cmroe )*ll(2)%f(N1)
       ll(4)%f(N1) = vfroe/rhoroe
       ll(5)%f(N1) = vfroe1/rhoroe
       !!! bubble-dynamic part !!!
       DO ib = 1,Nb
          DO ir = 1,NR0
             ! right eigenvectors at N1+1/2
             rr(ieig(ib,ir))%f(N1) = qval(ibub(ib,ir))%f(N1)
             ! left eigenvectors at N1+1/2
             ll(ieig(ib,ir))%f(N1) = qval(ibub(ib,ir))%f(N1)
          END DO
       END DO

    END IF

  END SUBROUTINE s_roeaverage

  !========================================================================

  SUBROUTINE s_multiply_rr( v_in,st,lr )

    INTEGER, INTENT(IN) :: st ! 0
    INTEGER, INTENT(IN) :: lr ! -1, 0 (left or right celledge)
    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: v_in
    INTEGER :: i
    INTEGER :: ir
    INTEGER :: ib

    DO i = 1,N1

       v_out(1)%f(i) = rr(1)%f(ip(i,lr))*v_in(1)%f(ip(i,st))     &
                     + rr(1)%f(ip(i,lr))*v_in(dir1)%f(ip(i,st))  &
                     + rr(1)%f(ip(i,lr))*v_in(Nveul)%f(ip(i,st))
       v_out(dir1)%f(i) = rr(2)%f(ip(i,lr))*v_in(1)%f(ip(i,st))     &
                        + rr(3)%f(ip(i,lr))*v_in(dir1)%f(ip(i,st))  &
                        + rr(4)%f(ip(i,lr))*v_in(Nveul)%f(ip(i,st))
       v_out(Nveul)%f(i) = rr(5)%f(ip(i,lr))*v_in(1)%f(ip(i,st))    &
                         + rr(5)%f(ip(i,lr))*v_in(dir1)%f(ip(i,st)) &
                         + ( rr(5)%f(ip(i,lr))-1.D0 ) &
                         * v_in(Nveul)%f(ip(i,st))
       DO ib = 1,Nb
          DO ir = 1,NR0
             v_out(ibub(ib,ir))%f(i) = &
                  rr(ieig(ib,ir))%f(ip(i,lr))*v_in(1)%f(ip(i,st))    &
                + rr(ieig(ib,ir))%f(ip(i,lr))*v_in(dir1)%f(ip(i,st)) &
                + v_in(ibub(ib,ir))%f(ip(i,st))
          END DO
       END DO

    END DO

  END SUBROUTINE s_multiply_rr

  !========================================================================

  SUBROUTINE s_multiply_ll( v_in,st,lr )

    INTEGER, INTENT(IN) :: st ! -2, -1, 0, 1, 2 (five stencils for WENO5)
    INTEGER, INTENT(IN) :: lr ! -1, 0 (left or right celledge)
    TYPE(coordinate), DIMENSION(Nv), INTENT(IN) :: v_in
    INTEGER :: i
    INTEGER :: ir
    INTEGER :: ib
    REAL(KIND(0.D0)) :: tmp

    DO i = 1,N1

       v_out(1)%f(i) = ll(1)%f(ip(i,lr))*v_in(1)%f(ip(i,st))     &
                     - ll(2)%f(ip(i,lr))*v_in(dir1)%f(ip(i,st))  &
                     + 0.5D0*v_in(Nveul)%f(ip(i,st))
       v_out(dir1)%f(i) = ll(3)%f(ip(i,lr))*v_in(1)%f(ip(i,st))     &
                        + ll(2)%f(ip(i,lr))*v_in(dir1)%f(ip(i,st))  &
                        + 0.5D0*v_in(Nveul)%f(ip(i,st))
       v_out(Nveul)%f(i) = ll(4)%f(ip(i,lr))*v_in(1)%f(ip(i,st)) &
                         - v_in(Nveul)%f(ip(i,st))
       tmp = ll(5)%f(ip(i,lr))
       DO ib = 1,Nb
          DO ir = 1,NR0
                v_out(ibub(ib,ir))%f(i) = &
                   - tmp*ll(ieig(ib,ir))%f(ip(i,lr))*v_in(1)%f(ip(i,st)) &
                   - ll(ieig(ib,ir))%f(ip(i,lr))*v_in(Nveul)%f(ip(i,st)) &
                   + v_in(ibub(ib,ir))%f(ip(i,st))
          END DO
       END DO

    END DO

  END SUBROUTINE s_multiply_ll

  !========================================================================

  SUBROUTINE s_celledgevalue

    ! creates celledge values used in HLLC
    INTEGER :: i
    INTEGER :: iv
    INTEGER :: ir
    REAL(KIND(0.D0)) :: vf1

    ! conserved variables
    DO iv = 1,Nv
       DO i = 1,N1
          lqval(iv)%f(i) = rweno(iv)%f(i) ! "L" state
          rqval(iv)%f(i) = lweno(iv)%f(i) ! "R" state
       END DO
    END DO
    IF ( bound%beg=='reflect'.AND.mpi_rank==0 ) THEN
       ! moving wall BC incorporated
       lqval(1)%f(0) = lweno(1)%f(1)
       lqval(dir1)%f(0) = 2.D0*uwall - lweno(dir1)%f(1)
       DO iv = Nveul,Nv
          lqval(iv)%f(0) = lweno(iv)%f(1)
       END DO
    END IF
    IF ( bound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN
       rqval(1)%f(N1+1) = rweno(1)%f(N1)
       rqval(dir1)%f(N1+1) = -rweno(dir1)%f(N1)
       DO iv = Nveul,Nv
          rqval(iv)%f(N1+1) = rweno(iv)%f(N1)
       END DO
    END IF

    ! velocity, pressure & flux
    DO i = Nbeg,N1
       lvel(i) = lqval(dir1)%f(i)/lqval(1)%f(i)
       vf1 = 1.D0 - lqval(Nveul)%f(i)
       lpres(i) = ( pl0+B_tait )*( lqval(1)%f(i)/vf1 )**n_tait - B_tait
       lflux(1)%f(i) = lqval(dir1)%f(i)
       lflux(dir1)%f(i) = lqval(dir1)%f(i)*lvel(i) + lpres(i)
       DO iv = Nveul,Nv
          lflux(iv)%f(i) = lqval(iv)%f(i)*lvel(i)
       END DO
    END DO
    DO i = 1,Nend
       rvel(i) = rqval(dir1)%f(i)/rqval(1)%f(i)
       vf1 = 1.D0 - rqval(Nveul)%f(i)
       rpres(i) = ( pl0+B_tait )*( rqval(1)%f(i)/vf1 )**n_tait - B_tait
       rflux(1)%f(i) = rqval(dir1)%f(i)
       rflux(dir1)%f(i) = rqval(dir1)%f(i)*rvel(i) + rpres(i)
       DO iv = Nveul,Nv
          rflux(iv)%f(i) = rqval(iv)%f(i)*rvel(i)
       END DO
    END DO

  END SUBROUTINE s_celledgevalue

  !========================================================================

  SUBROUTINE s_hllwavespeed

    INTEGER :: i
    REAL(KIND(0.D0)) :: soundsqr
    REAL(KIND(0.D0)) :: soundsql
    REAL(KIND(0.D0)) :: tmpl, tmpr
    REAL(KIND(0.D0)) :: lp, lrho, lu, lvf ! "L" state
    REAL(KIND(0.D0)) :: rp, rrho, ru, rvf ! "R" state

    DO i = Nbeg,N1
       rp = rpres(ip(i,1))
       lp = lpres(i)
       ru = rvel(ip(i,1))
       lu = lvel(i)
       rrho = rqval(1)%f(ip(i,1))
       lrho = lqval(1)%f(i)
       rvf = rqval(Nveul)%f(ip(i,1))
       lvf = lqval(Nveul)%f(i)
       soundsqr = n_tait*( rp+B_tait )/( 1.D0-rvf )/rrho
       soundsql = n_tait*( lp+B_tait )/( 1.D0-lvf )/lrho
       sr(i) = MAX( charsproe(2)%f(i),ru+DSQRT(soundsqr) )
       sl(i) = MIN( charsproe(1)%f(i),lu-DSQRT(soundsql) )
       ! middle wavespeed
       tmpl = sl(i) - lu
       tmpr = sr(i) - ru
       sstar(i) = ( rp-lp-rqval(dir1)%f(ip(i,1))*tmpr &
                + lqval(dir1)%f(i)*tmpl )/( lrho*tmpl-rrho*tmpr )
    END DO

  END SUBROUTINE s_hllwavespeed

  !========================================================================

  SUBROUTINE s_stop_rhsvar

    DEALLOCATE( xvel )
    NULLIFY( ip )
    NULLIFY( ds )
    NULLIFY( betaw0 )
    NULLIFY( betaw1 )
    NULLIFY( betaw2 )
    NULLIFY( dwenor )
    NULLIFY( dwenol )
    NULLIFY( polyr0 )
    NULLIFY( polyr1 )
    NULLIFY( polyr2 )
    NULLIFY( polyl0 )
    NULLIFY( polyl1 )
    NULLIFY( polyl2 )
    NULLIFY( w3beta )
    NULLIFY( w3dr )
    NULLIFY( w3dl )
    NULLIFY( w3polyr )
    NULLIFY( w3polyl )
    NULLIFY( vel1 )
    NULLIFY( onesidebeg )
    NULLIFY( onesideend )
    DEALLOCATE( pres )
    DEALLOCATE( sound )
    DEALLOCATE( nbub )
    DEALLOCATE( lpres )
    DEALLOCATE( rpres )
    DEALLOCATE( lvel )
    DEALLOCATE( rvel )
    DEALLOCATE( sl )
    DEALLOCATE( sr )
    DEALLOCATE( sstar )
    DEALLOCATE( nsfx )
    DEALLOCATE( velhf )
    DEALLOCATE( nfx )
    DEALLOCATE( rhsfx )
    DEALLOCATE( rhsrc )
    DEALLOCATE( fxtmp )
    DEALLOCATE( qval )
    DEALLOCATE( lflux )
    DEALLOCATE( rflux )
    DEALLOCATE( inweno )
    DEALLOCATE( lweno )
    DEALLOCATE( rweno )
    DEALLOCATE( charsp )
    DEALLOCATE( charsproe )
    DEALLOCATE( linw )
    DEALLOCATE( rinw )
    DEALLOCATE( v_out )
    DEALLOCATE( lqval )
    DEALLOCATE( rqval )
    DEALLOCATE( lqstar )
    DEALLOCATE( rqstar )
    DEALLOCATE( bub )
    DEALLOCATE( rr )
    DEALLOCATE( ll )

  END SUBROUTINE s_stop_rhsvar

  !========================================================================

END MODULE m_rhsvar
