!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module computing pre-/post-processing data
!
!  Last update: April 23, 2009
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_inoutvar

  USE m_globalvar
  USE m_misc
  IMPLICIT NONE

  REAL(KIND(0.D0)) :: R0ref
  REAL(KIND(0.D0)) :: sd
  REAL(KIND(0.D0)) :: vf0
  REAL(KIND(0.D0)) :: plnew
  REAL(KIND(0.D0)) :: CFLmx
  REAL(KIND(0.D0)), dimension(:),allocatable :: myRho, myVf, myN
  REAL(KIND(0.D0)), dimension(:,:),allocatable :: myR
  INTEGER :: i_probe1
  INTEGER :: i_probe2
  INTEGER :: Nout

  CONTAINS

  !========================================================================

  SUBROUTINE s_parameter

    !!! computational conditions for Euler equations !!!
    timescheme = 'rk3tvd' ! euler, rk3tvd
    timesplit = 'n'       ! y or n
    mpweno = 'y'          ! y or n (monotonicity preserving: Balsara 2000)
    artcomp = 'n'         ! y or n (artificial compression: Yang 1990)
    source = 'y'          ! y or n (static bubbles, zero sources)

    !!! computational conditions for bubble dynamics !!!
    viscous = 'y' ! n for Re = \infty
    polytropic = 'n' ! 'n' for Preston's or Sugiyama's model
    thermal = 'transfer' !'adiabat' ! transfer, adiabat, isotherm
    model = 'Preston' ! Preston or Sugiyama
    nquad = 'g' ! s (Simpson) or g (Gauss-Hermite)
    ! nquad = 's' ! s (Simpson) or g (Gauss-Hermite)

    !!! liquid and gas !!!
    liquid = 'water' ! test_liquid, water, silicone, glycerol
    gas = 'air' ! air, nitrogen, SF6
    vapor = 'n'

    !!! lognormal bubble size distribution !!!
    R0ref = 10.D-6
    !sd = 0.05D0
    sd = 0.3D0
    disperse = 'poly' ! mono or poly
    ! disperse = 'mono' ! mono or poly

    !!! initial void fraction at pl0 !!!
    vf0 = 1d-3

    !!! ICs !!!
    ! for linear computation
    ! 1: narrow perturbation
    ! 2: bubble screen
    ! for nonlinear computation
    ! 1: shock tube
    ! 2: steady shock
    ! 3: Kameda (1998), nitrogen/silicone
    ! 4: Kameda (1998), SF6/silicone
    ! 5: Beylich (1990), SF6/glycerine
    ! 6: bubble screen
    nonlin = 'n'
    ictype = '1' ! '97'

    !!! no. of cells, domain size & CFL !!!
    !!! BCs (periodic, reflect, nonreflect) !!!
    ! for reflecting BCs, solid walls placed at celledges (s_celledgevalue)
    ! but possibly placed at cell centers (see Johnsen 2007)
    ! for nonreflecting BCs, Thompson-type BCs (1987) implemented
    Nx = 800 !300
    length%x = 8000.D0
    cfl = 0.1D0
    xbound%beg = 'periodic' !'nonreflect'
    xbound%end = 'periodic' !'nonreflect'

    !!! probe measurement !!!
    timehis = 'n'

    !!! no. of flow files (excluding IC) !!!
    Nout = 200

    !!! given setups !!!
    CALL s_given_setup

    !!! allocate global variables !!!
    CALL s_start_inoutvar

    !!! nondimensional parameters !!!
    ! Bubble-dynamics-related parameters are computed, assuming that
    ! all the bubbles are initially at equilibrium (pl0,Tl0).
    IF ( polytropic=='y' ) CALL s_physprop_poly
    IF ( polytropic=='n' ) CALL s_physprop_trans

  END SUBROUTINE s_parameter

  !========================================================================

  SUBROUTINE s_given_setup

    IF ( nonlin=='y'.AND.ictype=='3' ) THEN
       liquid = 'silicone'
       gas = 'nitrogen'
       vapor = 'n'
       R0ref = 573.D-6
       vf0 = 0.18D-2
       Nx = 7000
       length%x = 8930.191972D0
       cfl = 0.8D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis = 'y'
    ELSE IF ( nonlin=='y'.AND.ictype=='97' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'y'
       viscous = 'y'
       R0ref = 50.e-6
       vf0 = 0.005d0
       Nx = 300 
       length%x = 1.d0
       cfl = 0.1d0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis = 'n'
   ELSE IF ( nonlin=='y'.AND.ictype=='99'.or.ictype=='98' ) THEN
       liquid = 'test_liquid' !water'
       gas = 'air'
       vapor = 'n'
       viscous = 'n'
       R0ref = 1.d0 !10.e-6
       vf0 = 0.005d0
       Nx = 300 !601
       length%x = 1.d0
       cfl = 0.1d0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='4' ) THEN
       liquid = 'silicone'
       gas = 'SF6'
       vapor = 'n'
       R0ref = 613.D-6
       vf0 = 0.24D-2
       Nx = 7000
       length%x = 8347.5D0
!       cfl = 0.8D0
       cfl = 0.1D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis = 'y'
    ELSE IF ( nonlin=='y'.AND.ictype=='5' ) THEN
       liquid = 'glycerol'
       gas = 'SF6'
       vapor = 'n'
       R0ref = 1.15D-3
       vf0 = 0.25D-2
       Nx = 2400
       length%x = 2543.5D0
       cfl = 0.8D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis= 'y'
   ELSE IF ( nonlin=='n'.AND.ictype=='1' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'n'
       R0ref = 10.D-6
       vf0 = 0.1D-2
       Nx = 1500
       length%x = 1500.D0
       cfl = 0.1D0
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
       timehis = 'y'
    ELSE IF ( nonlin=='n'.AND.ictype=='2' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'n'
       R0ref = 10.D-6
       vf0 = 0.1D-2
       Nx = 5000
       length%x = 2500.D0
       cfl = 0.1D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis = 'y'
    END IF

  END SUBROUTINE s_given_setup

  !========================================================================

  SUBROUTINE s_start_inoutvar

    INTEGER :: ir

    ! no. of conserved variables for the Euler equations
    Nveul = 3
    Nveul1 = Nveul + 1
    ! no. of bubble-dynamics-related variables
    IF ( polytropic=='y' ) Nb = 2
    IF ( polytropic=='n'.AND.vapor=='n' ) Nb = 3
    IF ( polytropic=='n'.AND.vapor=='y' ) Nb = 4
    ! odd no. of descretization points in R0
    IF ( disperse=='mono' ) NR0 = 1
    IF ( disperse=='poly' ) NR0 = 21 !201 
    ! no. of conserved variables for the entire system
    Nv = Nveul + Nb*NR0
    ! no. for artificial compression
    Nart = 3
    ! used for left & right eigenvectors
    Neig1 = 5
    Neig = Neig1 + Nb*NR0

    ! allocation
    ALLOCATE( ibub(Nb,NR0) )
    ALLOCATE( ieig(Nb,NR0) )
    IF ( polytropic=='n' ) THEN
       ALLOCATE( k_n(NR0) )
       ALLOCATE( k_v(NR0) )
       ALLOCATE( pb0(NR0) )
       ALLOCATE( mass_n0(NR0) )
       ALLOCATE( mass_v0(NR0) )
       ALLOCATE( Pe_T(NR0) )
       ALLOCATE( Re_trans_T(NR0) )
       ALLOCATE( Re_trans_c(NR0) )
       ALLOCATE( Im_trans_T(NR0) )
       ALLOCATE( Im_trans_c(NR0) )
       ALLOCATE( omegaN(NR0) )
    END IF
    ALLOCATE( R0(NR0) )
    ALLOCATE( weight(NR0) )

    ! index matrices
    DO ir = 1,NR0
       ! for bubble-dynamic variables
       ! (1,:): nR, (2,:): n\dot{R}, (3,:): np_b, (4,:): nm_v
       ibub(1,ir) = ir + Nveul
       ibub(2,ir) = ibub(1,ir) + NR0
       IF ( polytropic=='n' ) THEN
          ibub(3,ir) = ibub(2,ir) + NR0
          IF ( vapor=='y' ) THEN
             ibub(4,ir) = ibub(3,ir) + NR0
          END IF
       END IF
       ! for eigenvectors
       ieig(1,ir) = ir + Neig1
       ieig(2,ir) = ieig(1,ir) + NR0
       IF ( polytropic=='n' ) THEN
          ieig(3,ir) = ieig(2,ir) + NR0
          IF ( vapor=='y' ) THEN
             ieig(4,ir) = ieig(3,ir) + NR0
          END IF
       END IF
    END DO

  END SUBROUTINE s_start_inoutvar

  !========================================================================

  SUBROUTINE s_physprop_poly

    INTEGER :: ir
    !REAL(KIND(0.D0)) :: rhol0
    REAL(KIND(0.D0)) :: mul0
    REAL(KIND(0.D0)) :: ss
    REAL(KIND(0.D0)) :: uu
    REAL(KIND(0.D0)) :: omega_ref

    !!! dimensional, physical properties !!!
    ! liquid at STP
    
    
    IF ( liquid=='water' ) THEN
       n_tait = 7.15D0
       B_tait = 304.9D6
       rhol0 = 998.2063D0
       mul0 = 1.002D-3
       ss = 0.07275D0
       ! vapor
       pv = 2.3388D3
    ELSE IF ( liquid=='test_liquid' ) THEN
       n_tait = 4.4D0
       B_tait = 1.d0
       rhol0 = 1.D0
       mul0 = 0.d0    !no viscosity
       ss = 0.d0  !no surface tension !1.25D0
       ! vapor
       pv = 0.D0
   ELSE IF ( liquid=='silicone' ) THEN
       n_tait = 10.D0
       B_tait = 92.4D6
       rhol0 = 960.D0
       mul0 = 0.048D0
       ss = 20.8D-3
       ! vapor
       pv = 0.D0
    ELSE IF ( liquid=='glycerol' ) THEN
       n_tait = 8.75D0
       B_tait = 494.6D6
       rhol0 = 1220.D0
       mul0 = 0.0802D0
       ss = 69.D-3
       ! vapor
       pv = 0.D0
    END IF
    pl0 = 101325.D0
    IF ( ictype=='3' ) pl0 = 114.4D3
    IF ( ictype=='4' ) pl0 = 112.9D3
    IF ( ictype=='5' ) pl0 = 111.0D3
    IF ( ictype=='99' ) pl0 = 1.D0
    IF ( ictype=='98' ) pl0 = 1.D0
    IF ( vapor=='n' ) pv = 0.D0
    IF ( viscous=='n' ) mul0 = 0.D0
    ! gas
    IF ( gas=='air' ) gamma_m = 1.4D0
    IF ( gas=='nitrogen' ) gamma_m = 1.399D0
    IF ( gas=='SF6' ) gamma_m = 1.09D0
    IF ( thermal=='isotherm' ) gamma_m = 1.D0
    ! characteristic velocity
    uu = DSQRT( pl0/rhol0 ) ! = 1 for my test case

    !!! nondimensional parameters !!!
    ! initial radius
    IF ( NR0==1 ) THEN
       ! monodisperse
       R0(1) = 1.D0
       weight(1) = 1.D0
       sd = 0.D0
       NR0beg = 1
    ELSE
       ! polydisperse (lognormal distributions)
       IF ( nquad=='g' ) THEN
          ! Gauss-Hermite quadrature
          CALL s_gauher( NR0 )
       ELSE IF ( nquad=='s' ) THEN
          ! Simpson's rule
          CALL s_simpson( NR0 )
       END IF
    END IF
    ! cavitation, Weber & Reynolds no. based on R0ref
    Ca = ( pl0-pv )/( rhol0*uu**2 ) ! = pl0 = 1 in my test case
    print*, 'Ca = ', Ca
    We = rhol0*uu**2*R0ref/ss
    Re_inv = mul0/( rhol0*uu*R0ref )
    B_tait = B_tait/pl0

    !!! isothermal natural frequency associated with R0ref !!!
    omega_ref = 3.D0*gamma_m*Ca + 2.D0*( 3.D0*gamma_m-1.D0 )/We
    IF ( omega_ref<=0.D0 ) THEN
       WRITE(*,*) 'Unstable initial conditions. Calculation stopped.'
       STOP
    END IF

    !!! normalized vapor, ambient and internal bubble pressure !!!
    pv = pv/pl0
    pl0 = 1.D0

  END SUBROUTINE s_physprop_poly

  !========================================================================

  SUBROUTINE s_physprop_trans

    INTEGER :: ir
    REAL(KIND(0.D0)) :: rhol0
    REAL(KIND(0.D0)) :: mul0
    REAL(KIND(0.D0)) :: ss
    REAL(KIND(0.D0)) :: uu
    REAL(KIND(0.D0)) :: mu_n, mu_v
    REAL(KIND(0.D0)) :: D_m
    REAL(KIND(0.D0)) :: gamma_n, gamma_v
    REAL(KIND(0.D0)) :: temp
    REAL(KIND(0.D0)) :: omega_ref
    REAL(KIND(0.D0)), DIMENSION(NR0) :: chi_vw0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: cp_m0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: k_m0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: rho_m0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: x_vw
    ! polytropic index used to compute isothermal natural frequency
    REAL(KIND(0.D0)), PARAMETER :: k_poly = 1.D0
    ! universal gas constant
    REAL(KIND(0.D0)), PARAMETER :: Ru = 8314.D0

    !!! dimensional, physical properties !!!
    ! liquid
    IF ( liquid=='water' ) THEN
       n_tait = 7.15D0
       B_tait = 304.9D6
       rhol0 = 998.2063D0
       mul0 = 1.002D-3
       ss = 0.07275D0
       ! vapor
       pv = 2.3388D3
       gamma_v = 1.33D0
       M_v = 18.02D0
       mu_v = 0.8816D-5
       k_v = 0.019426D0
   elseIF ( liquid=='test_liquid' ) THEN
       n_tait = 4.4d0
       B_tait = 0.d0
       rhol0 = 1.d0
       mul0 = 0.d0
       ss = 1.d0
       ! vapor
       pv = 1.d0
       gamma_v = 1.4d0
       M_v = 18.02D0
       mu_v = 0.d0
       k_v = 0.019426D0
    ELSE IF ( liquid=='silicone' ) THEN
       n_tait = 10.D0
       B_tait = 92.4D6
       rhol0 = 960.D0
       mul0 = 0.048D0
       ss = 20.8D-3
       ! vapor (irrelevant)
       pv = 0.D0
       gamma_v = 1.33D0
       M_v = 18.02D0
       mu_v = 0.8816D-5
       k_v = 0.019426D0
    ELSE IF ( liquid=='glycerol' ) THEN
       n_tait = 8.75D0
       B_tait = 494.6D6
       rhol0 = 1220.D0
       mul0 = 0.0802D0
       ss = 69.D-3
       ! vapor (irrelevant)
       pv = 0.D0
       gamma_v = 1.33D0
       M_v = 18.02D0
       mu_v = 0.8816D-5
       k_v = 0.019426D0
    END IF
    pl0 = 101325.D0
    IF ( ictype=='3' ) pl0 = 114.4D3
    IF ( ictype=='4' ) pl0 = 112.9D3
    IF ( ictype=='5' ) pl0 = 111.0D3
    IF ( ictype=='99' ) pl0 = 1.
    IF ( ictype=='98' ) pl0 = 1.
    temp = 293.15D0
    IF ( ictype=='3'.OR.ictype=='4' ) temp = 298.15D0
    IF ( ictype=='5' ) temp = 298.D0
    IF ( vapor=='n' ) pv = 0.D0
    IF ( viscous=='n' ) mul0 = 0.D0
    ! gas
    IF ( gas=='air' ) THEN
       gamma_n = 1.4D0
       gamma_m = gamma_n
       M_n = 28.97D0
       mu_n = 1.8D-5
       k_n = 0.02556D0
    ELSE IF ( gas=='nitrogen' ) THEN
       gamma_n = 1.399D0
       gamma_m = gamma_n
       M_n = 28.013D0
       mu_n = 1.76D-5
       k_n = 0.02534D0
    ELSE IF ( gas=='SF6' ) THEN
       gamma_n = 1.09D0
       gamma_m = gamma_n
       M_n = 146.06D0
       mu_n = 1.76D-5
       k_n = 14.D-3
    END IF
    IF ( thermal=='isotherm' ) gamma_m = 1.D0
    ! mass diffusion
    IF ( gas=='air'.AND.liquid=='water' ) THEN
       D_m = 0.242D-4
    ELSE
       ! temporary value
       D_m = 0.242D-4
    END IF
    ! characteristic velocity
    uu = DSQRT( pl0/rhol0 ) ! = 1 in my test case 99

    !!! nondimensional parameters !!!
    ! initial radius
    IF ( NR0==1 ) THEN
       ! monodisperse
       R0(1) = 1.D0
       weight(1) = 1.D0
       sd = 0.D0
       NR0beg = 1
    ELSE
       ! polydisperse (lognormal distributions)
       IF ( nquad=='g' ) THEN
          ! Gauss-Hermite quadrature
          CALL s_gauher( NR0 )
       ELSE IF ( nquad=='s' ) THEN
          ! Simpson's rule
          CALL s_simpson( NR0 )
       END IF
    END IF

    print*, 'rhol0, u^2, R0ref, ss', rhol0, uu*uu, R0ref, ss
    ! cavitation, Weber & Reynolds no. based on R0ref
    Ca = ( pl0-pv )/( rhol0*uu**2 )
    We = rhol0*uu**2*R0ref/ss
    Re_inv = mul0/( rhol0*uu*R0ref )
    B_tait = B_tait/pl0

    !!! isothermal natural frequency associated with R0ref !!!
    omega_ref = 3.D0*k_poly*Ca + 2.D0*( 3.D0*k_poly-1.D0 )/We
    IF ( omega_ref<=0.D0 ) THEN
       WRITE(*,*) 'Unstable initial conditions. Calculation stopped.'
       STOP
    END IF

    !!! thermal properties !!!
    ! gas constants
    R_n = Ru/M_n
    R_v = Ru/M_v
    ! phi_vn & phi_nv (phi_nn = phi_vv = 1)
    phi_vn = ( 1.D0+DSQRT(mu_v/mu_n)*(M_n/M_v)**(0.25D0) )**2 &
           / ( DSQRT(8.D0)*DSQRT(1.D0+M_v/M_n) )
    phi_nv = ( 1.D0+DSQRT(mu_n/mu_v)*(M_v/M_n)**(0.25D0) )**2 &
           / ( DSQRT(8.D0)*DSQRT(1.D0+M_n/M_v) )
    ! internal bubble pressure
    pb0 = pl0 + 2.D0*ss/( R0ref*R0 )
    IF ( vapor=='n' ) THEN
       ! mass fraction of vapor
       chi_vw0 = 0.D0
       ! specific heat for gas
       cp_m0 = R_n*gamma_n/( gamma_n-1.D0 )
       ! mole fraction of vapor
       x_vw = 0.D0
       ! thermal conductivity for gas
       k_m0 = k_n
       ! gas density
       rho_m0 = pb0/( R_n*temp )
    ELSE
       ! mass fraction of vapor
       chi_vw0 = 1.D0/( 1.D0+R_v/R_n*(pb0/pv-1.D0) )
       ! specific heat for gas/vapor mixture
       cp_m0 = chi_vw0*R_v*gamma_v/( gamma_v-1.D0 ) &
             + ( 1.D0-chi_vw0 )*R_n*gamma_n/( gamma_n-1.D0 )
       ! mole fraction of vapor
       x_vw = M_n*chi_vw0/( M_v+(M_n-M_v)*chi_vw0 )
       ! thermal conductivity for gas/vapor mixture
       k_m0 = x_vw*k_v/( x_vw+(1.D0-x_vw)*phi_vn ) &
            + ( 1.D0-x_vw )*k_n/( x_vw*phi_nv+1.D0-x_vw )
       ! mixture density
       rho_m0 = pv/( chi_vw0*R_v*temp )
    END IF
    ! mass of gas/vapor computed using dimensional quantities
    mass_n0 = 4.D0*( pb0-pv )*pi/( 3.D0*R_n*temp*rhol0 )*R0**3
    mass_v0 = 4.D0*pv*pi/( 3.D0*R_v*temp*rhol0 )*R0**3
    ! Peclet numbers
    Pe_T = rho_m0*cp_m0*uu*R0ref/k_m0
    Pe_c = uu*R0ref/D_m
    ! nondimensional properties
    R_n = rhol0*R_n*temp/pl0
    R_v = rhol0*R_v*temp/pl0
    k_n = k_n/k_m0
    k_v = k_v/k_m0
    pb0 = pb0/pl0
    pv = pv/pl0
    ! bubble wall temperature, normalized by T0, in the liquid
    ! keeps a constant (cold liquid assumption)
    Tw = 1.D0
    ! natural frequencies
    omegaN = DSQRT( 3.D0*k_poly*Ca+2.D0*(3.D0*k_poly-1.D0)/(We*R0) )/R0

    !!! normalized ambient pressure !!!
    pl0 = 1.D0

    !!! constant transfer coefficients !!!
    ! use of natural frequency is recommended by Preston (2004).
    DO ir = 1,NR0
       CALL s_transcoeff( omegaN(ir)*R0(ir),Pe_T(ir)*R0(ir), &
                          Re_trans_T(ir),Im_trans_T(ir) )
       CALL s_transcoeff( omegaN(ir)*R0(ir),Pe_c*R0(ir), &
                          Re_trans_c(ir),Im_trans_c(ir) )
    END DO
    IF ( model=='Preston' ) THEN
       Im_trans_T = 0.D0
       Im_trans_c = 0.D0
    END IF

    print*, 'pb0,mass_v0', pb0(:), mass_v0(:)

  END SUBROUTINE s_physprop_trans

  !========================================================================

  SUBROUTINE s_transcoeff( omega,peclet,Re_trans,Im_trans )

    REAL(KIND(0.D0)), INTENT(IN) :: omega
    REAL(KIND(0.D0)), INTENT(IN) :: peclet
    REAL(KIND(0.D0)), INTENT(OUT) :: Re_trans
    REAL(KIND(0.D0)), INTENT(OUT) :: Im_trans
    COMPLEX :: trans, c1, c2, c3
    COMPLEX :: imag = ( 0.,1. )
    REAL(KIND(0.D0)) :: f_transcoeff

    c1 = imag*omega*peclet
    c2 = CSQRT( c1 )
    c3 = ( CEXP(c2)-CEXP(-c2) )/( CEXP(c2)+CEXP(-c2) ) ! TANH(c2)
    trans = ( (c2/c3-1.D0)**(-1)-3.D0/c1 )**( -1 ) ! transfer function

    Re_trans = DBLE( trans )
    Im_trans = AIMAG( trans )

  END SUBROUTINE s_transcoeff

  !========================================================================

  SUBROUTINE s_gauher( Npt )

    REAL(KIND(0.D0)), PARAMETER :: psmall = 3.D-14
    REAL(KIND(0.D0)), PARAMETER :: pim4 = 0.7511255444649425D0 ! pi^(-1/4)
    REAL(KIND(0.D0)), PARAMETER :: R0mn = 0.05D0
    INTEGER, PARAMETER :: mxit = 10

    INTEGER :: i, j
    INTEGER :: its, m
    INTEGER, INTENT(IN) :: Npt
    REAL(KIND(0.D0)) :: p1, p2, p3, pp
    REAL(KIND(0.D0)) :: z, z1
    REAL(KIND(0.D0)), DIMENSION(Npt) :: phi_tld

    ! routine for Gauss-Hermite abscissas and weights (Numerical Recipe)
    ! Roots are symmetric about the origin, then find only half of them
    m = ( Npt+1 )/2
    DO i = 1,m

       IF ( i==1 ) THEN
          z = DSQRT( DBLE(2*Npt+1) ) - 1.85575D0*( DBLE(2*Npt+1) ) &
            **( -0.16667D0 )
       ELSE IF ( i==2 ) THEN
          z = z - 1.14D0*( DBLE(Npt) )**0.426D0/z
       ELSE IF ( i==3 ) THEN
          z = 1.86D0*z - 0.86D0*phi_tld(1)
       ELSE IF ( i==4 ) THEN
          z = 1.91D0*z - 0.91D0*phi_tld(2)
       ELSE
          z = 2.D0*z - phi_tld(i-2)
       END IF
       its = 1
       DO
          IF ( its>mxit.OR.DABS(z-z1)<=psmall ) EXIT
          p1 = pim4
          p2 = 0.D0
          DO j = 1,Npt
             p3 = p2
             p2 = p1
             p1 = z*DSQRT( 2.D0/DBLE(j) )*p2 &
                - DSQRT( DBLE(j-1)/DBLE(j) )*p3
          END DO
          pp = DSQRT( 2.D0*DBLE(Npt) )*p2
          z1 = z
          z = z1 - p1/pp
          its = its + 1
       END DO
       ! assign the root
       phi_tld(i) = z
       phi_tld(Npt+1-i) = -z
       ! assign the weight
       weight(i) = 2.D0/( pp*pp )
       weight(Npt+1-i) = weight(i)

    END DO

    ! normalize the weights
    weight = weight/DSQRT( pi )
    ! transform phi_tld into R0, R0(1) = R0mx, R0(Npt) = R0mn
    R0 = DEXP( DSQRT(2.D0)*sd*phi_tld )

    ! sorting R0 s.t. R0(1) = R0mn, R0(Npt) = R0mx
    ! phi_tld is dummy
    DO i = 1,Npt
       phi_tld(Npt+1-i) = R0(i)
    END DO
    R0 = phi_tld
    ! compute NR0beg
    NR0beg = 1
    DO
       IF ( R0(NR0beg)>=R0mn ) EXIT
       NR0beg = NR0beg + 1
    END DO

    print*, 'R0', R0(:)
    print*, 'weight', weight(:)
    stop

  END SUBROUTINE s_gauher

  !========================================================================

  SUBROUTINE s_simpson( Npt )

    INTEGER, INTENT(IN) :: Npt
    INTEGER :: ir
    REAL(KIND(0.D0)) :: R0mn
    REAL(KIND(0.D0)) :: R0mx
    REAL(KIND(0.D0)) :: dphi
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)), DIMENSION(Npt) :: phi

    ! nondiml. min. & max. initial radii for numerical quadrature
    IF ( sd==0.05D0 ) THEN
       R0mn = 0.75D0
       R0mx = 1.3D0
    ELSE IF ( sd==0.1D0 ) THEN
       R0mn = 0.6D0
       R0mx = 1.7D0
    ELSE IF ( sd==0.3D0 ) THEN
       R0mn = 0.3D0
       R0mx = 6.D0
    ELSE IF ( sd==0.5D0 ) THEN
       R0mn = 0.17D0
       R0mx = 25.D0
    ELSE IF ( sd==0.7D0 ) THEN
       R0mn = 0.12D0
       R0mx = 150.D0
    END IF

    ! phi = ln( R0 ) & return R0
    DO ir = 1,Npt
       phi(ir) = DLOG( R0mn ) &
               + DBLE( ir-1 )*DLOG( R0mx/R0mn )/DBLE( Npt-1 )
       R0(ir) = DEXP( phi(ir) )
    END DO
    dphi = phi(2) - phi(1)

    ! weights for quadrature using Simpson's rule
    DO ir = 2,Npt-1
       ! Gaussian
       tmp = DEXP( -0.5D0*(phi(ir)/sd)**2 )/DSQRT( 2.D0*pi )/sd
       IF ( MOD(ir,2)==0 ) THEN
          weight(ir) = tmp*4.D0*dphi/3.D0
       ELSE
          weight(ir) = tmp*2.D0*dphi/3.D0
       END IF
    END DO
    tmp = DEXP( -0.5D0*(phi(1)/sd)**2 )/DSQRT( 2.D0*pi )/sd
    weight(1) = tmp*dphi/3.D0
    tmp = DEXP( -0.5D0*(phi(Npt)/sd)**2 )/DSQRT( 2.D0*pi )/sd
    weight(Npt) = tmp*dphi/3.D0
    NR0beg = 1

    print*, 'weight = ', weight(:)
    print*, 'R0s = ', R0(:)

    open(101,file='D/R0s')
    do ir = 1,nR0
    write(101,*) R0(ir)
    end do
    close(101)
  END SUBROUTINE s_simpson

  !========================================================================

  SUBROUTINE s_newequilibrium( plnew_tmp,Rnew,pbnew,mvnew )

    ! compute new equliribrium states corresponding to plnew
    ! The new state is achieved in an isothermal process.
    ! subscripts "i" denote initial conditions (0) associated with pl0
    REAL(KIND(0.D0)), INTENT(IN) :: plnew_tmp
    REAL(KIND(0.D0)), INTENT(OUT) :: Rnew
    REAL(KIND(0.D0)), INTENT(OUT) :: pbnew
    REAL(KIND(0.D0)), INTENT(OUT) :: mvnew
    REAL(KIND(0.D0)), PARAMETER :: acc = 1.D-8
    REAL(KIND(0.D0)), PARAMETER :: Rmn = 1.D-30
    REAL(KIND(0.D0)) :: Rmx

    ! find Rnew using Newton-Raphson
    plnew = plnew_tmp
    IF ( plnew>= pl0 ) THEN
       Rmx = R0(iR0) + acc
    ELSE
       Rmx = R0(iR0)*1.D6
    END IF
    Rnew = f_rtsafe( Rmn,Rmx,acc )

    ! new internal bubble pressure
    pbnew = plnew + 2.D0/We/Rnew

    ! new mass of vapor
    mvnew = 4.D0*pv*pi/( 3.D0*R_v )*Rnew**3

  END SUBROUTINE s_newequilibrium

  !========================================================================

  SUBROUTINE s_funcd( x,y,dydx )

    REAL(KIND(0.D0)), INTENT(IN) :: x
    REAL(KIND(0.D0)), INTENT(OUT) :: y
    REAL(KIND(0.D0)), INTENT(OUT) :: dydx
    REAL(KIND(0.D0)) :: pni
    REAL(KIND(0.D0)) :: power

    ! initial partial pressure for noncondensible gas
    IF ( polytropic=='n' ) pni = pb0(iR0) - pv
    IF ( polytropic=='y' ) pni = Ca + 2.D0/We/R0(iR0)

    ! returns both the function value and the first derivative
    IF ( thermal=='adiabat' ) THEN
       power = 3.D0*gamma_m
       y = pni*( R0(iR0)/x )**power - 2.D0/We/x + pv - plnew
       dydx = 2.D0/We/x**2 - power*pni*R0(iR0)**power/x**( power+1 )
    ELSE ! isothermal
       y = pni*( R0(iR0)/x )**3 - 2.D0/We/x + pv - plnew
       dydx = 2.D0/We/x**2 - 3.D0*pni*R0(iR0)**3/x**4
    END IF

  END SUBROUTINE s_funcd

  !========================================================================

  FUNCTION f_rtsafe( x1,x2,xacc )

    REAL(KIND(0.D0)), INTENT(IN) :: x1
    REAL(KIND(0.D0)), INTENT(IN) :: x2
    REAL(KIND(0.D0)), INTENT(IN) :: xacc
    REAL(KIND(0.D0)) :: f_rtsafe ! return the root

    INTEGER, PARAMETER :: maxit = 100 ! max. iterations
    INTEGER :: iter
    REAL(KIND(0.D0)) :: df
    REAL(KIND(0.D0)) :: delx
    REAL(KIND(0.D0)) :: delxold
    REAL(KIND(0.D0)) :: f
    REAL(KIND(0.D0)) :: fh, fl
    REAL(KIND(0.D0)) :: xh, xl
    REAL(KIND(0.D0)) :: tmp
    LOGICAL :: judge

    ! using a combination of Newton-Raphson and bisection, find the root
    ! of a function bracketed between x1 and x2. The root will be refined
    ! until its accuracy is known within xacc. 's_funcd' is a user-supplied
    ! subroutine which returns both the function value and the first
    ! derivative of the function.

    CALL s_funcd( x1,fl,df )
    CALL s_funcd( x2,fh,df )

    judge = .TRUE.
    IF ( (fl>0.D0.AND.fh>0.D0).OR.(fl<0.D0.AND.fh<0.D0) ) THEN
       WRITE(*,*) 'The root is not bracketed in f_rtsafe'
       judge = .FALSE.
    ELSE IF ( fl==0.D0 ) THEN
       f_rtsafe = x1
       judge = .FALSE.
    ELSE IF ( fh==0.D0 ) THEN
       f_rtsafe = x2
       judge = .FALSE.
    ELSE IF ( fl<0.D0.AND.fh>0.D0 ) THEN
       xl = x1
       xh = x2
    ELSE IF ( fl>0.D0.AND.fh<0.D0 ) THEN
       xh = x1
       xl = x2
    END IF


    IF ( judge ) THEN

       f_rtsafe = 0.5D0*( x1+x2 )
       delxold = DABS( x2-x1 )
       delx = delxold
       CALL s_funcd( f_rtsafe,f,df )
       iter = 1
       DO
          IF ( ((f_rtsafe-xh)*df-f)*((f_rtsafe-xl)*df-f)>0.D0.OR. &
             DABS(2.D0*f)>DABS(delxold*df) ) THEN
             delxold = delx
             delx = 0.5D0*( xh-xl )
             f_rtsafe = xl + delx
             IF ( xl==f_rtsafe ) EXIT
          ELSE
             delxold = delx
             delx = f/df
             tmp = f_rtsafe
             f_rtsafe = f_rtsafe - delx
             IF ( tmp==f_rtsafe ) EXIT
          END IF
          IF ( DABS(delx)<xacc ) EXIT
          CALL s_funcd( f_rtsafe,f,df )
          IF ( f<0.D0 ) THEN
             xl = f_rtsafe
          ELSE
             xh = f_rtsafe
          END IF
          iter = iter + 1
          IF ( iter>maxit ) THEN
             WRITE(*,*) 'f_rtsafe exceeding max. iterations'
             EXIT
          END IF
       END DO

    END IF

  END FUNCTION f_rtsafe

  !========================================================================

  SUBROUTINE s_gridgen

    INTEGER :: i
    INTEGER :: j
    INTEGER :: iv
    INTEGER :: ibeg
    INTEGER :: ixbeg
    INTEGER :: iend

    ! set ends
    IF ( xbound%beg=='reflect' ) THEN
       ibeg = 1 - wenonum
       ixbeg = 0
    ELSE
       ibeg = 1
       ixbeg = 1
    END IF
    IF ( xbound%end=='reflect' ) THEN
       iend = Nx + wenonum
    ELSE
       iend = Nx
    END IF

    ! allocate rhs, primitive and conserved variables
    ALLOCATE( rhs(Nv),prim(Nv),cons(Nv),dqbub(Nv) )
    DO iv = 1,Nv
       ALLOCATE( rhs(iv)%f(1:Nx) )
       ALLOCATE( prim(iv)%f(1:Nx) )
       ALLOCATE( cons(iv)%f(1:Nx) )
       ALLOCATE( dqbub(iv)%f(1:Nx) )
    END DO
    ! allocate grids, WENO coefficients
    ALLOCATE( ix(0:Nx,-wenonum:wenonum) )
    ALLOCATE( xgrid(1:Nx) )
    ALLOCATE( dxs(ibeg:iend) )
    ALLOCATE( betax0(0:2,1:Nx) )
    ALLOCATE( betax1(0:2,1:Nx) )
    ALLOCATE( betax2(0:2,1:Nx) )
    ALLOCATE( polyrx0(0:1,1:Nx) )
    ALLOCATE( polyrx1(0:1,1:Nx) )
    ALLOCATE( polyrx2(0:1,1:Nx) )
    ALLOCATE( polylx0(0:1,1:Nx) )
    ALLOCATE( polylx1(0:1,1:Nx) )
    ALLOCATE( polylx2(0:1,1:Nx) )
    ALLOCATE( dwxr(0:2,1:Nx) )
    ALLOCATE( dwxl(0:2,1:Nx) )
    ! WENO3 used for nonreflecting BCs
    ALLOCATE( w3betax(0:1,1:Nx) )
    ALLOCATE( w3drx(0:1,1:Nx) )
    ALLOCATE( w3dlx(0:1,1:Nx) )
    ALLOCATE( w3polyrx(0:1,1:Nx) )
    ALLOCATE( w3polylx(0:1,1:Nx) )
    ALLOCATE( onesidebegx(1:5) )
    ALLOCATE( onesideendx(1:5) )

    ! uniform grid
    dx = length%x/DBLE( Nx )
    dxs = dx
    DO i = 1,Nx
       xgrid(i) = ( DBLE(i)-0.5D0 )*dx
    END DO

    ! index matrix with appropriate boundary treatment
    ! for WENO5, five stincils (-2,-1,0,1,2) used
    DO i = ixbeg,Nx
       DO j = -wenonum,wenonum
          ix(i,j) = i + j
       END DO
    END DO
    IF ( xbound%beg=='periodic'.AND.xbound%end=='periodic' ) THEN
       ix(2,-2) = Nx
       ix(1,-2) = Nx - 1
       ix(1,-1) = Nx
       ix(Nx-1,2) = 1
       ix(Nx,2) = 2
       ix(Nx,1) = 1
    END IF
    IF ( xbound%beg=='nonreflect' ) THEN
       ix(2,-2) = 1
       ix(1,-2) = 1
       ix(1,-1) = 1
    END IF
    IF ( xbound%end=='nonreflect' ) THEN
       ix(Nx-1,2) = Nx
       ix(Nx,2) = Nx
       ix(Nx,1) = Nx
    END IF

  END SUBROUTINE s_gridgen

  !========================================================================

  SUBROUTINE s_probe_setup

    IF ( timehis=='n' ) THEN

       i_probe1 = Nx
       i_probe2 = Nx

    ELSE IF ( timehis=='y'.AND.nonlin=='y' ) THEN

       IF ( ictype=='3' ) THEN
          i_probe1 = 3*Nx/7
          i_probe2 = i_probe1 + 1
       ELSE IF ( ictype=='4' ) THEN
          i_probe1 = 3*Nx/7
          i_probe2 = i_probe1 + 1
       ELSE IF ( ictype=='5' ) THEN
          i_probe1 = Nx/2
          i_probe2 = i_probe1 + 1
       END IF

    ELSE IF ( timehis=='y'.AND.nonlin=='n' ) THEN

       IF ( ictype=='1' ) THEN
          i_probe1 = NINT( 25.D0/dx )
          i_probe2 = i_probe1 + 1
       ELSE IF ( ictype=='2' ) THEN
          i_probe1 = Nx*9/25
          i_probe2 = Nx*16/25
       END IF

    END IF

  END SUBROUTINE s_probe_setup

  !========================================================================

  SUBROUTINE s_ic

    INTEGER :: i
    INTEGER :: ir
    INTEGER :: iv
    INTEGER :: ib
    INTEGER :: diaph
    REAL(KIND(0.D0)) :: fx
    REAL(KIND(0.D0)) :: small
    REAL(KIND(0.D0)) :: R3bar
    REAL(KIND(0.D0)) :: rhor, rhol
    REAL(KIND(0.D0)) :: vel1r, vel1l
    REAL(KIND(0.D0)) :: nr, nl, n0
    REAL(KIND(0.D0)) :: vfr, vfl
    REAL(KIND(0.D0)) :: cmr, cml, cm0
    REAL(KIND(0.D0)) :: Smx
    REAL(KIND(0.D0)) :: vftiny, ntiny
    REAL(KIND(0.D0)) :: vftmp
    REAL(KIND(0.D0)) :: pltmp
    REAL(KIND(0.D0)) :: ratio
    REAL(KIND(0.D0)), DIMENSION(NR0) :: Rtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pbtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: mvtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pb0tmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: mv0tmp

    ! substitute initial pb0 and mass_v0
    pb0tmp = 0.D0!/n
    mv0tmp = 0.D0
    IF ( polytropic=='n' ) THEN
       pb0tmp = pb0
       mv0tmp = mass_v0
    END IF

    IF ( nonlin=='y' ) THEN

       IF ( ictype=='1' ) THEN
          ! sonic speed at 0 Hz for pl0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2
          cm0 = DSQRT( cm0 )
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          
          !!! right state associated with pl0 !!!
          ! at equilibrium
          ! made dimensionless by rhol0, ratio = 1 so pl = pl0
          rhor = 1.D0 - vf0
          vel1r = 0.D0
          nr = n0
          vfr = vf0
          
          !!! left state associated with pltmp !!!
          pltmp = 5.D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl )
          vel1l = 0.D0
          
          diaph = 2*Nx/4
          finaltime = 8.D0*length%x/cm0
          dt = cfl*dx/cm0
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
      ELSEIF ( ictype=='99' ) THEN
          ! SOS for CFL condition. sonic speed at 0 Hz for pl0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2
          cm0 = DSQRT(cm0)
          cm0 = 1.23251662148 !1.18916176545
          !cm0 = 1.25377
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          
          !!! right state associated with pl0 !!!
          pltmp = 1.3*pl0
          nr = n0
          vfr = vf0
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          rhor = ratio*(1-vf0)
          vel1r = 0.D0
          
          !!! left state associated with pltmp !!!
          pltmp = 1.5*pl0
          !DO ir = 1,NR0
          !   iR0 = ir
          !   CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          !END DO
          Rtmp(:) = 2*R0; pbtmp(:) = 0.; mvtmp(:) = 0.
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0 !*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl ) !1.33591698254          
          vel1l = 0.D0
          
          diaph = ceiling(Nx/2.)
          dt = cfl*dx/cm0
          Nt = 500 !NINT( finaltime/dt )
          finaltime = Nt*dt
          interval = 100 !Nt/Nout
          CFLmx = CFL

          print*, 'rho_L/R = ', rhol, rhoR
          !print*, 'p_L/R = ', pltmp, pl0
          print*, 'vf_L/R = ', vfl, vfr
          print*, 'n_L/R = ', nl, nr
          print*, 'vel_L/R = ', vel1l, vel1r
          print*, 'R_L/R = ', Rtmp(1), R0(1)
      ELSEIF ( ictype=='98' ) THEN
          !smooth version of 99 SHB
          ! SOS for CFL condition. sonic speed at 0 Hz for pl0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2d0
          cm0 = DSQRT(cm0)
          print*, 'SOS: ', cm0

          allocate( myRho(Nx), myN(Nx), myVf(Nx) )
          allocate( myR(Nx,NR0) )

          pbtmp(:) = 0.d0; mvtmp(:) = 0.d0
          Rtmp(:) = R0; 
          CALL s_quad( R0**3,R3bar ) 
          myN(:) = vf0/( pi43*R3bar )
          myVf(:) = vf0
          !vfl = pi43*nl*(Rtmp**3.)

          do i = 1,NR0
            myR(:,i) = R0(i)
          end do
          
          vel1l = 0.D0
          
          do i = 1,Nx 
            !myR(i,1)  = R0(1)*  ( 1.d0 + 0.2d0*exp(-1d0*(i*dx-0.5d0)**2.d0/(2.d0*0.005d0)) )
            !myVf(i) = vf0*    ( 1.d0 + 0.2d0*exp(-1d0*(i*dx-0.5d0)**2.d0/(2.d0*0.005d0)) )  
            !myN(i)  = myVf(i)/( pi43*(myR(i)**3.d0) )

            pltmp = 1.5d0*pl0 * ( 1d0 + 0.2d0*exp(-1d0*(i*dx-0.5d0)**2.d0/(2.d0*0.005d0)) )
            ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
            !myRho(i) = ratio*(1.d0-myVf(i))
            myRho(i) = ratio*(1.d0-vf0)
          end do

          dt = cfl*dx/cm0
          Nt = 10 !500*2 !NINT( finaltime/dt )
          finaltime = Nt*dt
          interval = 10 !Nt/Nout
          CFLmx = CFL

          print*, 'dt = ', dt
          print*, 'Re_inv, We, Ca'
          print*, Re_inv, We, Ca
      ELSEIF ( ictype=='97' ) THEN
          !smooth version of 99 SHB
          ! SOS for CFL condition. sonic speed at 0 Hz for pl0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2d0
          cm0 = DSQRT(cm0)
          print*, 'SOS: ', cm0

          allocate( myRho(Nx), myN(Nx), myVf(Nx) )
          allocate( myR(Nx,NR0) )

          pbtmp(:) = 0.d0; mvtmp(:) = 0.d0
          Rtmp(:) = R0; 
          CALL s_quad( R0**3,R3bar ) 
          myN(:) = vf0/( pi43*R3bar )
          myVf(:) = vf0
          !vfl = pi43*nl*(Rtmp**3.)

          do i = 1,NR0
            myR(:,i) = R0(i)
          end do
          
          vel1l = 0.D0
          
          do i = 1,Nx 
            !myR(i,1)  = R0(1)*  ( 1.d0 + 0.2d0*exp(-1d0*(i*dx-0.5d0)**2.d0/(2.d0*0.005d0)) )
            !myVf(i) = vf0*    ( 1.d0 + 0.2d0*exp(-1d0*(i*dx-0.5d0)**2.d0/(2.d0*0.005d0)) )  
            !myN(i)  = myVf(i)/( pi43*(myR(i)**3.d0) )

            pltmp = 1.5d0*pl0 * ( 1d0 + 0.2d0*exp(-1d0*(i*dx-0.5d0)**2.d0/(2.d0*0.005d0)) )
            ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
            !myRho(i) = ratio*(1.d0-myVf(i))
            myRho(i) = ratio*(1.d0-vf0)
          end do

          dt = cfl*dx/cm0
          Nt = 500 !500*2 !NINT( finaltime/dt )
          finaltime = Nt*dt
          interval = 10 !Nt/Nout
          CFLmx = CFL

          print*, 'dt = ', dt
          print*, 'Re_inv, We, Ca'
          print*, Re_inv, We, Ca
      ELSE IF ( ictype=='2' ) THEN
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )

          !!! right state associated with pl0 !!!
          rhor = 1.D0 - vf0
          vel1r = 0.D0
          nr = n0
          vfr = vf0
          cmr = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          
          !!! left state associated with pltmp !!!
          pltmp = 2.D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl )
          vel1l = ( pltmp-pl0 )/rhor/( 1.D0-(1.D0-vf0)/(1.D0-vfl)/ratio )
          vel1l = DSQRT( vel1l )*( 1.D0-nr/nl )
          cml = DSQRT( n_tait*(pltmp+B_tait) )/( 1.D0-vfl )
          
          ! Numerical setup
          diaph = Nx/2
          Smx = MAX( DABS(vel1r)+cmr,DABS(vel1l)+cml )
          finaltime = 0.4*length%x/Smx
          dt = cfl*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL

          print*, 'rho_L/R = ', rhol, rhoR
          print*, 'p_L/R = ', pltmp, pl0
          print*, 'vf_L/R = ', vfl, vfr
          print*, 'n_L/R = ', nl, nr
          print*, 'vel_L/R = ', vel1l, vel1r
          print*, 'R_L/R = ', Rtmp(1), R0(1)
       ELSE IF ( ictype=='3' ) THEN
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state associated with pl0 !!!
          rhor = 1.D0 - vf0
          vel1r = -0.012275166D0
          nr = n0
          vfr = vf0
          cmr = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state associated with pltmp !!!
          pltmp = 1.749125874D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl )
          vel1l = ( pltmp-pl0 )/rhor/( 1.D0-(1.D0-vf0)/(1.D0-vfl)/ratio )
          vel1l = DSQRT( vel1l )*( 1.D0-nr/nl ) + vel1r
          cml = DSQRT( n_tait*(pltmp+B_tait) )/( 1.D0-vfl )
          ! Numerical setup
          diaph = Nx/7
          finaltime = 300.D0
          Smx = MAX( DABS(vel1r)+cmr,DABS(vel1l)+cml )
          dt = cfl*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='4' ) THEN
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state associated with pl0 !!!
          rhor = 1.D0 - vf0
          vel1r = -0.01236D0
          nr = n0
          vfr = vf0
          cmr = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state associated with pltmp !!!
          pltmp = 2.157D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl )
          vel1l = ( pltmp-pl0 )/rhor/( 1.D0-(1.D0-vf0)/(1.D0-vfl)/ratio )
          vel1l = DSQRT( vel1l )*( 1.D0-nr/nl ) + vel1r
          cml = DSQRT( n_tait*(pltmp+B_tait) )/( 1.D0-vfl )
          ! Numerical setup
          diaph = Nx/7
          finaltime = 300.D0
          Smx = MAX( DABS(vel1r)+cmr,DABS(vel1l)+cml )
          dt = cfl*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='5' ) THEN
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state associated with pl0 !!!
          rhor = 1.D0 - vf0
          vel1r = 0.D0
          nr = n0
          vfr = vf0
          cmr = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state associated with pltmp !!!
          pltmp = 1.81D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl )
          vel1l = ( pltmp-pl0 )/rhor/( 1.D0-(1.D0-vf0)/(1.D0-vfl)/ratio )
          vel1l = DSQRT( vel1l )*( 1.D0-nr/nl ) + vel1r
          cml = DSQRT( n_tait*(pltmp+B_tait) )/( 1.D0-vfl )
          ! numerical setup
          diaph = Nx/6
          finaltime = 100.D0
          Smx = MAX( DABS(vel1r)+cmr,DABS(vel1l)+cml )
          dt = cfl*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='6' ) THEN
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          vftiny = 1.D-7
          ntiny = ( vftiny/vf0 )*n0
          !!! right state associated with pl0 !!!
          rhor = 1.D0 - vftiny
          vel1r = 0.D0
          nr = ntiny
          vfr = vftiny
          cmr = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vftiny )
          !!! left state associated with pltmp !!!
          pltmp = 5.D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = ntiny*ratio/( 1.D0-vftiny+pi43*ntiny*R3bar*ratio )
          vfl = pi43*nl*R3bar
          rhol = ratio*( 1.D0-vfl )
          vel1l = ( pltmp-pl0 )/rhor/( 1.D0-(1.D0-vftiny)/(1.D0-vfl)/ratio )
          vel1l = DSQRT( vel1l )*( 1.D0-nr/nl )
          cml = DSQRT( n_tait*(pltmp+B_tait) )/( 1.D0-vfl )
          ! Numerical setup
          diaph = Nx/8
          Smx = MAX( DABS(vel1r)+cmr,DABS(vel1l)+cml )
          finaltime = 150. !length%x/Smx
          dt = cfl*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       END IF
       ! assign ICs for all grid points
       DO i = 1,Nx
          if ( (ictype .ne. '98') .and. (ictype .ne. '97') ) then
            IF ( i<=diaph ) THEN
             CALL s_assignIC( rhol,vel1l,vfl,Rtmp,pbtmp,mvtmp,nl,i )
            ELSE
             CALL s_assignIC( rhor,vel1r,vfr,R0,pb0tmp,mv0tmp,nr,i )
            END If
          else if (  (ictype == '98') .or. (ictype == '97') ) then
            call s_assignIC_smooth( myRho, vel1l, myvf, myR, pb0tmp, mv0tmp, myN, i )
          end if
       END DO

       ! bubble screen
       IF ( ictype=='6' ) THEN
          DO i = 1*Nx/4+1,3*Nx/4
             CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
          END DO
       END IF

    ELSE IF ( nonlin=='n' ) THEN

       IF ( ictype=='1' ) THEN
          ! perturbation in liquid pressure
          small = 1.D-4
          ! sonic speed at 0 Hz for pl0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2
          cm0 = DSQRT( cm0 )
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          DO i = 1,Nx
             fx = DEXP( -(xgrid(i)/4.D0)**2 )
             pltmp = pl0*( 1.D0+small*fx )
             DO ir = 1,NR0
                iR0 = ir
                CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
             END DO
             CALL s_quad( Rtmp**3,R3bar )
             ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
             nl = n0*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
             vfl = pi43*nl*R3bar
             rhol = ratio*( 1.D0-vfl )
             CALL s_assignIC( rhol,0.D0,vfl,Rtmp,pbtmp,mvtmp,nl,i )
          END DO
          finaltime = 2.D0*length%x/cm0
          dt = cfl*dx/cm0
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
          print*, 'Nt = ', Nt
       ELSE IF ( ictype=='2' ) THEN
          ! right-going wave
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          small = 1.D-4
          vftiny = 1.D-7
          ntiny = ( vftiny/vf0 )*n0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vftiny )**2
          cm0 = DSQRT( cm0 )
          DO i = 1,Nx
             fx = xgrid(i) - 0.5*( xgrid(19*Nx/50)+xgrid(19*Nx/50+1) )
             fx = DEXP( -(fx/4.D0)**2 )
             pltmp = pl0*( 1.D0+small*fx )
             DO ir = 1,NR0
                iR0 = ir
                CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
             END DO
             CALL s_quad( Rtmp**3,R3bar )
             ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
             nl = ntiny*ratio/( 1.D0-vftiny+pi43*ntiny*R3bar*ratio )
             vfl = pi43*nl*R3bar
             rhol = ratio*( 1.D0-vfl )
             vel1l = cm0*( rhol-1.D0+vftiny )
             CALL s_assignIC( rhol,vel1l,vfl,Rtmp,pbtmp,mvtmp,nl,i )
          END DO
          ! bubble screen
          DO i = 2*Nx/5+1,3*Nx/5
             CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
          END DO
          ! redifine sonic speed
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2
          cm0 = DSQRT( cm0 )
          finaltime = length%x/cm0
          dt = cfl*dx/cm0
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       END IF

    END IF

    ! initial time increment for adaptive RK
    IF ( timesplit=='y' ) THEN
       DO ir = 1,NR0
          dttry(ir)%f = 0.1D0*dt
       END DO
    END IF
    DO iv = 1,Nv
       dqbub(iv)%f = 0.D0
    END DO

  END SUBROUTINE s_ic

  !========================================================================

  SUBROUTINE s_assignIC( rho,vel1,vf,rr,pb,mv,nn,i )

    INTEGER, INTENT(IN) :: i
    REAL(KIND(0.D0)), INTENT(IN) :: rho, vel1, vf, nn
    REAL(KIND(0.D0)), DIMENSION(NR0), INTENT(IN) :: rr, pb, mv
    INTEGER :: ib, ir
    REAL(KIND(0.D0)) :: bubtmp

    ! Euler part
    prim(1)%f(i) = rho
    prim(2)%f(i) = vel1
    prim(Nveul)%f(i) = vf
    ! bubble-dynamic part
    DO ib = 1,Nb
       DO ir = 1,NR0
          IF ( ib==1 ) bubtmp = rr(ir)
          IF ( ib==2 ) bubtmp = 0.D0
          IF ( ib==3 ) bubtmp = pb(ir)
          IF ( ib==4 ) bubtmp = mv(ir)
          prim(ibub(ib,ir))%f(i) = bubtmp
          cons(ibub(ib,ir))%f(i) = nn*bubtmp
       END DO
    END DO

  END SUBROUTINE s_assignIC


  SUBROUTINE s_assignIC_smooth( rho2,vel1,vf,rr,pb,mv,nn,i )

    INTEGER, INTENT(IN) :: i
    REAL(KIND(0.D0)), INTENT(IN) :: vel1
    REAL(KIND(0.D0)), Dimension(Nx), INTENT(IN) :: rho2, nn, vf
    REAL(KIND(0.D0)), DIMENSION(NR0), INTENT(IN) :: pb, mv 
    REAL(KIND(0.D0)), DIMENSION(Nx,NR0), INTENT(IN) :: rr
    INTEGER :: ib, ir
    REAL(KIND(0.D0)) :: bubtmp

    ! Euler part
    prim(1)%f(i) = rho2(i)
    prim(2)%f(i) = vel1
    prim(Nveul)%f(i) = vf(i)
    ! bubble-dynamic part
    DO ib = 1,Nb
       DO ir = 1,NR0
          IF ( ib==1 ) bubtmp = rr(i,ir)
          IF ( ib==2 ) bubtmp = 0.d0
          IF ( ib==3 ) bubtmp = pb(ir)
          IF ( ib==4 ) bubtmp = mv(ir)
          prim(ibub(ib,ir))%f(i) = bubtmp
          cons(ibub(ib,ir))%f(i) = nn(i)*bubtmp
       END DO
    END DO

  END SUBROUTINE s_assignIC_smooth

  !========================================================================

  SUBROUTINE s_wenocoeff

    INTEGER :: i
    REAL(KIND(0.D0)) :: coeff
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: tmp1
    REAL(KIND(0.D0)) :: tmp2
    REAL(KIND(0.D0)) :: tmp10
    REAL(KIND(0.D0)) :: tmp20
    REAL(KIND(0.D0)) :: m1m5, p1m5, p1m3, p3m3, p3m1, p5m1, p5p1


    ! see appendix of Johnsen (2007)
    ! coefficients for polynomials
    DO i = 1,Nx
       m1m5 = dxs(ix(i,-1)) + dxs(ix(i,-2))
       p1m5 = m1m5 + dxs(i)
       p1m3 = dxs(ix(i,-1)) + dxs(i)
       p3m3 = p1m3 + dxs(ix(i,1))
       p3m1 = dxs(i) + dxs(ix(i,1))
       p5m1 = p3m1 + dxs(ix(i,2))
       p5p1 = p5m1 - dxs(i)
       polyrx0(0,i) = dxs(i)/p1m5 + dxs(i)/p1m3
       polyrx0(1,i) = -dxs(i)*p1m3/m1m5/ p1m5
       polyrx1(0,i) = dxs(i)*p1m3/p3m3/p3m1
       polyrx1(1,i) = dxs(i)*dxs(ix(i,1))/p1m3/p3m3
       polyrx2(0,i) = -dxs(i)*dxs(ix(i,1))/p5m1/p5p1
       polyrx2(1,i) = dxs(i)*( p3m1+p5p1 )/p3m1/p5m1
       polylx0(0,i) = -dxs(i)*( p1m3+m1m5 )/p1m5/p1m3
       polylx0(1,i) = dxs(i)*dxs(ix(i,-1))/m1m5/p1m5
       polylx1(0,i) = -dxs(i)*dxs(ix(i,-1))/p3m3/p3m1
       polylx1(1,i) = -dxs(i)*p3m1/p1m3/p3m3
       polylx2(0,i) = dxs(i)*p3m1/p5m1/p5p1
       polylx2(1,i) = -dxs(i)/p3m1 - dxs(i)/p5m1
    END DO

    ! ideal weights
    DO i = 1,Nx
       tmp = dxs(ix(i,-2)) + dxs(ix(i,-1)) + dxs(i) + dxs(ix(i,1)) &
           + dxs(ix(i,2))
       tmp1 = dxs(ix(i,-2)) + dxs(ix(i,-1)) + dxs(i) + dxs(ix(i,1))
       tmp2 = dxs(ix(i,-1)) + dxs(i) + dxs(ix(i,1)) + dxs(ix(i,2))
       dwxr(0,i) = dxs(ix(i,1))*(dxs(ix(i,1)) + dxs(ix(i,2)))/tmp1/tmp
       dwxr(2,i) = ( dxs(ix(i,-2))+dxs(ix(i,-1))+dxs(i) )*(dxs(ix(i,-1)) &
                 + dxs(i) )/tmp2/tmp
       dwxr(1,i) = 1.D0 - dwxr(2,i) - dwxr(0,i)
       dwxl(0,i) = ( dxs(i)+dxs(ix(i,1)) )*( dxs(i)+dxs(ix(i,1) ) &
                 + dxs(ix(i,2)) )/tmp1/tmp
       dwxl(2,i) = ( dxs(ix(i,-2))+dxs(ix(i,-1)) )*dxs(ix(i,-1))/tmp2/tmp
       dwxl(1,i) = 1.D0 - dwxl(2,i) - dwxl(0,i)
    END DO

    ! smoothness indicators
    DO i = 1,Nx
       tmp10 = 10.D0*dxs(i)*dxs(i)
       tmp20 = 2.D0*tmp10
       tmp = dxs(i) + dxs(ix(i,1))
       coeff = 4.D0*dxs(i)*dxs(i)/( tmp*tmp )
       tmp = tmp + dxs(ix(i,2))
       coeff = coeff/( tmp*tmp )
       tmp = dxs(ix(i,1)) + dxs(ix(i,2))
       coeff = coeff/( tmp*tmp )
       tmp = dxs(i) + dxs(ix(i,1))
       tmp1 = dxs(ix(i,1)) + dxs(ix(i,2))
       tmp2 = dxs(ix(i,2)) + 2.D0*tmp
       betax2(0,i) = coeff*tmp*tmp*( tmp10+dxs(ix(i,1))*dxs(ix(i,1)) &
                   + dxs(i)*dxs(ix(i,1)) )
       betax2(1,i) = -coeff*tmp*tmp1*( tmp20+dxs(ix(i,1))*(dxs(ix(i,2)) &
                   + 3.D0*tmp ) + tmp*(tmp1+dxs(i)) )
       betax2(2,i) = coeff*tmp1*tmp1*( tmp10+tmp2*tmp2-dxs(i)*tmp2 )
       tmp = dxs(ix(i,-1)) + dxs(i)
       coeff = 4.D0*dxs(i)*dxs(i)/( tmp*tmp )
       tmp = tmp + dxs(ix(i,1))
       coeff = coeff/( tmp*tmp )
       tmp = tmp - dxs(ix(i,-1))
       coeff = coeff/( tmp*tmp )
       tmp = dxs(ix(i,-1)) + dxs(i)
       tmp1 = dxs(i) + dxs(ix(i,1))
       betax1(0,i) = coeff*tmp*tmp*( tmp10+tmp*dxs(ix(i,-1)) )
       betax1(1,i) = coeff*tmp*tmp1*( tmp*tmp1 &
                   + dxs(ix(i,-1))*dxs(ix(i,1))-tmp20 )
       betax1(2,i) = coeff*tmp1*tmp1*( tmp10+tmp1*dxs(ix(i,1)) )
       tmp = dxs(ix(i,-1)) + dxs(i)
       coeff = 4.D0*dxs(i)*dxs(i)/( tmp*tmp )
       tmp = tmp + dxs(ix(i,-2))
       coeff = coeff/( tmp*tmp )
       tmp = tmp - dxs(i)
       coeff = coeff/( tmp*tmp )
       tmp = dxs(ix(i,-2)) + dxs(ix(i,-1))
       tmp1 = dxs(ix(i,-1)) + dxs(i)
       tmp2 = tmp + tmp1
       betax0(0,i) = coeff*tmp*tmp*( tmp10+tmp2*tmp2+dxs(i)*tmp2 )
       betax0(1,i) = -coeff*tmp*tmp1*( tmp20+2.D0*dxs(ix(i,-1))*tmp1 &
                   + (tmp+dxs(i))*(tmp1+dxs(ix(i,-1))) )
       betax0(2,i) = coeff*tmp1*tmp1*( tmp10+dxs(ix(i,-1))*dxs(ix(i,-1)) &
                   + dxs(i)*dxs(ix(i,-1)) )
    END DO

    ! coefficients for fourth-order one-sided difference
    IF ( xbound%beg=='nonreflect' ) THEN
       onesidebegx(1) = -( 1.D0/(xgrid(2)-xgrid(1)) &
                      + 1.D0/(xgrid(3)-xgrid(1))+1.D0/(xgrid(4)-xgrid(1)) &
                      + 1.D0/(xgrid(5)-xgrid(1)) ) ! -25/(12*dx)
       onesidebegx(2) = ( xgrid(3)-xgrid(1) )*( xgrid(4)-xgrid(1) ) &
                      * ( xgrid(5)-xgrid(1) )/( xgrid(2)-xgrid(1) ) &
                      / ( xgrid(3)-xgrid(2) )/( xgrid(4)-xgrid(2) ) &
                      / ( xgrid(5)-xgrid(2) ) ! 4/dx
       onesidebegx(3) = -( xgrid(2)-xgrid(1) )*( xgrid(4)-xgrid(1) ) &
                      * ( xgrid(5)-xgrid(1) )/( xgrid(3)-xgrid(1) ) &
                      / ( xgrid(3)-xgrid(2) )/( xgrid(4)-xgrid(3) ) &
                      / ( xgrid(5)-xgrid(3) ) ! -3/dx
       onesidebegx(4) = ( xgrid(2)-xgrid(1) )*( xgrid(3)-xgrid(1) ) &
                      * ( xgrid(5)-xgrid(1) )/( xgrid(4)-xgrid(1) ) &
                      / ( xgrid(4)-xgrid(2) )/( xgrid(4)-xgrid(3) ) &
                      / ( xgrid(5)-xgrid(4) ) ! 4/(3*dx)
       onesidebegx(5) = -( xgrid(2)-xgrid(1) )*( xgrid(3)-xgrid(1) ) &
                      * ( xgrid(4)-xgrid(1) )/( xgrid(5)-xgrid(1) ) &
                      / ( xgrid(5)-xgrid(2) )/( xgrid(5)-xgrid(3) ) &
                      / ( xgrid(5)-xgrid(4) ) ! -1/(4*dx)
    END IF
    IF ( xbound%end=='nonreflect' ) THEN
       onesideendx(1) = -( 1.D0/(xgrid(Nx-1)-xgrid(Nx)) &
                      + 1.D0/(xgrid(Nx-2)-xgrid(Nx)) &
                      + 1.D0/(xgrid(Nx-3)-xgrid(Nx)) &
                      + 1.D0/(xgrid(Nx-4)-xgrid(Nx)) ) ! 25/(12*dx)
       onesideendx(2) = ( xgrid(Nx-2)-xgrid(Nx) ) &
                      * ( xgrid(Nx-3)-xgrid(Nx) ) &
                      * ( xgrid(Nx-4)-xgrid(Nx) ) &
                      / ( xgrid(Nx-1)-xgrid(Nx) ) &
                      / ( xgrid(Nx-2)-xgrid(Nx-1) ) &
                      / ( xgrid(Nx-3)-xgrid(Nx-1) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx-1) ) ! -4/dx
       onesideendx(3) = -( xgrid(Nx-1)-xgrid(Nx) ) &
                      * ( xgrid(Nx-3)-xgrid(Nx) ) &
                      * ( xgrid(Nx-4)-xgrid(Nx) ) &
                      / ( xgrid(Nx-2)-xgrid(Nx) ) &
                      / ( xgrid(Nx-2)-xgrid(Nx-1) ) &
                      / ( xgrid(Nx-3)-xgrid(Nx-2) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx-2) ) ! 3/dx
       onesideendx(4) = ( xgrid(Nx-1)-xgrid(Nx) ) &
                      * ( xgrid(Nx-2)-xgrid(Nx) ) &
                      * ( xgrid(Nx-4)-xgrid(Nx) ) &
                      / ( xgrid(Nx-3)-xgrid(Nx) ) &
                      / ( xgrid(Nx-3)-xgrid(Nx-1) ) &
                      / ( xgrid(Nx-3)-xgrid(Nx-2) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx-3) ) ! -4/(3*dx)
       onesideendx(5) = -( xgrid(Nx-1)-xgrid(Nx) ) &
                      * ( xgrid(Nx-2)-xgrid(Nx) ) &
                      * ( xgrid(Nx-3)-xgrid(Nx) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx-1) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx-2) ) &
                      / ( xgrid(Nx-4)-xgrid(Nx-3) ) ! 1/(4*dx)
    END IF

    ! coefficients for WENO3 used for nonreflecting BCs
    DO i = 1,Nx
       tmp = dxs(i) + dxs(ix(i,1))
       tmp1 = dxs(ix(i,-1)) + dxs(i)
       tmp2 = dxs(ix(i,-1)) + dxs(i) + dxs(ix(i,1))
       w3betax(0,i) = 4.D0*dxs(i)*dxs(i)/tmp/tmp
       w3betax(1,i) = 4.D0*dxs(i)*dxs(i)/tmp1/tmp1
       w3drx(0,i) = dxs(ix(i,1))/tmp2
       w3drx(1,i) = tmp1/tmp2
       w3dlx(0,i) = tmp/tmp2
       w3dlx(1,i) = dxs(ix(i,-1))/tmp2
       w3polyrx(0,i) = dxs(i)/tmp1
       w3polyrx(1,i) = dxs(i)/tmp
       w3polylx(0,i) = -w3polyrx(0,i)
       w3polylx(1,i) = -w3polyrx(1,i)
    END DO

  END SUBROUTINE s_wenocoeff

  !========================================================================

  SUBROUTINE s_prim2cons

    !!! Euler equations !!!
    ! density
    cons(1)%f = prim(1)%f
    ! momentum
    cons(2)%f = prim(1)%f*prim(2)%f
    ! void fraction
    cons(Nveul)%f = prim(Nveul)%f

  END SUBROUTINE s_prim2cons

  !========================================================================

  SUBROUTINE s_cons2prim

    INTEGER :: i
    INTEGER :: ir
    INTEGER :: ib
    REAL(KIND(0.D0)) :: ntmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: nRtmp

    !!! Euler equations !!!
    ! density
    prim(1)%f = cons(1)%f
    ! velocity
    prim(2)%f = cons(2)%f/cons(1)%f
    ! void fraction
    prim(Nveul)%f = cons(Nveul)%f

    !!! bubble-dynamic equations !!!
    DO i = 1,Nx

       ! compute total bubble number density
       DO ir = 1,NR0
          nRtmp(ir) = cons(ibub(1,ir))%f(i)
       END DO
       CALL s_comp_n( cons(Nveul)%f(i),nRtmp,ntmp )
       DO ib = 1,Nb
          DO ir = 1,NR0
             prim(ibub(ib,ir))%f(i) = cons(ibub(ib,ir))%f(i)/ntmp
          END DO
       END DO

    END DO

  END SUBROUTINE s_cons2prim

  !========================================================================

  SUBROUTINE s_output

    INTEGER :: i,jj
    INTEGER :: iv
    INTEGER :: ir
    INTEGER :: NR0mid
    CHARACTER(8) :: today_date
    CHARACTER(10) :: today_time
    CHARACTER(10) :: notmp
    CHARACTER(3) :: no
    REAL(KIND(0.D0)) :: vftmp
    REAL(KIND(0.D0)) :: vftmp1
    REAL(KIND(0.D0)) :: pl
    REAL(KIND(0.D0)) :: cltmp
    REAL(KIND(0.D0)) :: Stmp
    REAL(KIND(0.D0)) :: CFLtmp
    REAL(KIND(0.D0)) :: CFLmx1

    IF ( it==0 ) THEN

       ! numerical specification
       NR0mid = ( NR0+1 )/2
       CALL DATE_AND_TIME( today_date,today_time )
       OPEN(UNIT=1,FILE='D/setup.dat')
       WRITE(1,*) 'date: ', today_date, ', time: ', today_time
       WRITE(1,*) 'nonlinear computation: ', nonlin
       WRITE(1,*) 'IC type: ', ictype
       WRITE(1,*) 'time marching: ', timescheme
       WRITE(1,*) 'CFL =', cfl
       WRITE(1,*) 'final time =', finaltime
       WRITE(1,*) 'timestep splitting: ', timesplit
       WRITE(1,*) 'WENO order: ', wenoord
       WRITE(1,*) 'monotonicity preserving: ', mpweno
       WRITE(1,*) 'artificial compression: ', artcomp
       WRITE(1,*) 'source terms: ', source
       WRITE(1,*) 'viscous computation: ', viscous
       WRITE(1,*) 'polytropic gas: ', polytropic
       WRITE(1,*) 'thermal behavior: ', thermal
       WRITE(1,*) 'transfer model: ', model
       WRITE(1,*) 'gas/vapor bubbles: ', vapor
       WRITE(1,*) 'liquid: ', liquid
       WRITE(1,*) 'noncondensible gas: ', gas
       WRITE(1,*) 'quadrature (Simpson or Gauss): ', nquad
       WRITE(1,*) 'Nx =', Nx,', Nt =', Nt
       WRITE(1,*) 'flow file output at intervals of ', interval
       WRITE(1,*) 'length in x: ',length%x
       WRITE(1,*) 'dt =', dt
       WRITE(1,*) 'left BC: ', xbound%beg, ', right BC: ', xbound%end
       WRITE(1,*) 'n_tait =', n_tait
       WRITE(1,*) 'B_tait =', B_tait
       WRITE(1,*) 'gamma_m =', gamma_m
       IF ( polytropic=='n' ) THEN
          WRITE(1,*) 'Pe_T(NR0mid) =', Pe_T(NR0mid)
          WRITE(1,*) 'Pe_c =', Pe_c
       END IF
       WRITE(1,*) 'Ca =', Ca
       WRITE(1,*) 'We =', We
       IF ( viscous=='y' ) THEN
          WRITE(1,*) 'Re =', 1.D0/Re_inv
       END IF
       WRITE(1,*) 'R0ref for f(R0) [m]:', R0ref
       WRITE(1,*) 'standard deviation for f(R0):', sd
       WRITE(1,*) 'initial void fraction =', vf0
       WRITE(1,*) 'NR0 =', NR0
       WRITE(1,*) 'NR0beg =', NR0beg
       CLOSE(1)
       ! weights
       OPEN(UNIT=2,FILE='D/weights.dat',STATUS='UNKNOWN',FORM='FORMATTED')
       DO ir = 1,NR0
          WRITE(2,*) R0(ir), weight(ir)
       END DO
       CLOSE(2)
       ! grid
       OPEN(UNIT=3,FILE='D/xgrid.dat',STATUS='UNKNOWN',FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(3,*) xgrid(i)
       END DO
       CLOSE(3)

    END IF


    ! primitive variables
    CALL s_cons2prim

    ! update max. CFL
    CFLmx1 = 0.D0
    DO i = 1,Nx
       vftmp = prim(Nveul)%f(i)
       vftmp1 = 1.D0 - vftmp
       pl = ( pl0+B_tait )*( prim(1)%f(i)/vftmp1 )**n_tait
       cltmp = DSQRT( n_tait*pl/vftmp1/prim(1)%f(i) )
       Stmp = DABS( prim(2)%f(i) ) + cltmp
       CFLtmp = dt*Stmp/dxs(i)
       IF ( CFLmx1<CFLtmp ) CFLmx1 = CFLtmp
    END DO


    IF ( CFLmx<CFLmx1 ) THEN
       CFLmx = CFLmx1
       OPEN(UNIT=4,FILE='D/maxcfl.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(4,*) time, CFLmx
       CLOSE(4)
    END IF

    ! time history
    IF ( timehis=='y' ) THEN

       i = i_probe1
       vftmp = prim(Nveul)%f(i)
       vftmp1 = 1.D0 - vftmp
       pl = ( pl0+B_tait )*( prim(1)%f(i)/vftmp1 )**n_tait - B_tait
       OPEN(UNIT=5,FILE='D/probe1.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(5,*) time, pl
       CLOSE(5)
       i = i_probe2
       vftmp = prim(Nveul)%f(i)
       vftmp1 = 1.D0 - vftmp
       pl = ( pl0+B_tait )*( prim(1)%f(i)/vftmp1 )**n_tait - B_tait
       OPEN(UNIT=6,FILE='D/probe2.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(6,*) time, pl
       CLOSE(6)

    END IF

    do i = 1,5
        do jj = 1,nx
            !print*, 'j,i,',i,jj,prim(i)%f(jj)
            if (isnan(prim(i)%f(jj))) then
                print*, 'Output: Prim nan', it, i, jj
                stop 'nan'
            end if
        end do
    end do

    do i = 1,5
        do jj = 1,nx
            !print*, 'j,i,',i,jj,cons(i)%f(jj)
            if (isnan(cons(i)%f(jj))) then
                print*, 'Output: Cons nan', it, i, jj
                stop 'nan'
            end if
        end do
    end do

    !stop 'done'


    ! flow files
    IF ( it==0.OR.MOD(it,interval)==0 ) THEN

!       WRITE(notmp,*) itout
       WRITE(notmp,'(I7)') itout
       IF ( itout<10 ) THEN
          no = '00'//TRIM( ADJUSTL(notmp) )
       ELSE IF ( itout>=10.AND.itout<100 ) THEN
          no = '0'//TRIM( ADJUSTL(notmp) )
       ELSE
          no = TRIM( ADJUSTL(notmp) )
       END IF
       ! output time
       OPEN(UNIT=7,FILE='D/flowtime.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(7,*) itout, time
       CLOSE(7)
       OPEN(UNIT=8,FILE='D/p.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          vftmp = prim(Nveul)%f(i)
          vftmp1 = 1.D0 - vftmp
          pl = ( pl0+B_tait )*( prim(1)%f(i)/vftmp1 )**n_tait - B_tait
          !WRITE(8,*) xgrid(i), ( prim(iv)%f(i),iv=1,Nveul ), pl
          !WRITE(8,*) xgrid(i), ',', prim(1)%f(i)
          WRITE(8,*) i*dx,pl !/pl0
          !WRITE(8,*) i*dx,prim(1)%f(i) !/pl0
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/alf.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,prim(3)%f(i)
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/rho.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,prim(1)%f(i)
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/v.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,prim(2)%f(i)
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/R.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,prim(4)%f(i)
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/nR.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,cons(4)%f(i)
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/Rdot.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,prim(5)%f(i)
       END DO
       CLOSE(8)


       OPEN(UNIT=8,FILE='D/nRdot.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,cons(5)%f(i)
       END DO
       CLOSE(8)

       OPEN(UNIT=8,FILE='D/nbub.'//no//'.csv',STATUS='UNKNOWN', &
            FORM='FORMATTED')
       DO i = 1,Nx
          WRITE(8,*) i*dx,cons(4)%f(i)/prim(4)%f(i)
       END DO
       CLOSE(8)

       itout = itout + 1

    END IF

  END SUBROUTINE s_output

  !========================================================================

  SUBROUTINE s_stop_inoutvar

    DEALLOCATE( rhs )
    DEALLOCATE( prim )
    DEALLOCATE( cons )
    DEALLOCATE( dqbub )
    DEALLOCATE( ibub )
    DEALLOCATE( ieig )
    IF ( polytropic=='n' ) THEN
       DEALLOCATE( k_n )
       DEALLOCATE( k_v )
       DEALLOCATE( pb0 )
       DEALLOCATE( mass_n0 )
       DEALLOCATE( mass_v0 )
       DEALLOCATE( Pe_T )
       DEALLOCATE( Re_trans_T )
       DEALLOCATE( Re_trans_c )
       DEALLOCATE( Im_trans_T )
       DEALLOCATE( Im_trans_c )
       DEALLOCATE( omegaN )
    END IF
    DEALLOCATE( R0 )
    DEALLOCATE( weight )
    DEALLOCATE( ix )
    DEALLOCATE( xgrid )
    DEALLOCATE( dxs )
    DEALLOCATE( betax0 )
    DEALLOCATE( betax1 )
    DEALLOCATE( betax2 )
    DEALLOCATE( dwxr )
    DEALLOCATE( dwxl )
    DEALLOCATE( polyrx0 )
    DEALLOCATE( polyrx1 )
    DEALLOCATE( polyrx2 )
    DEALLOCATE( polylx0 )
    DEALLOCATE( polylx1 )
    DEALLOCATE( polylx2 )
    DEALLOCATE( w3betax )
    DEALLOCATE( w3drx )
    DEALLOCATE( w3dlx )
    DEALLOCATE( w3polyrx )
    DEALLOCATE( w3polylx )
    DEALLOCATE( onesidebegx )
    DEALLOCATE( onesideendx )

  END SUBROUTINE s_stop_inoutvar

  !========================================================================

END MODULE m_inoutvar
