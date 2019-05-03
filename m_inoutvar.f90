!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module computing pre-/post-processing data
!
!  Last update: July 23, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_inoutvar

  USE m_globalvar
  USE m_misc
  USE mpi_setup
  IMPLICIT NONE

  ! two probes to measure pressure
  ! used to compute dispersion relation of linear wave propagation
  INTEGER :: mpi_rank_probe1, mpi_rank_probe2
  INTEGER :: i_probe1, i_probe2
  ! used to output integrands R(x,t;R0)
  INTEGER :: mpi_rank_integ1, mpi_rank_integ2
  INTEGER :: i_integ1, i_integ2
  ! grid for entire domain
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xgrid_tot
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: dxs_tot
  ! cell-edge location
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xhalf_tot
  ! primitive variables and bubble radius for entire domain
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: prim_tot
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rad_tot
  INTEGER, DIMENSION(:), ALLOCATABLE :: i_Nx_nopad
  ! used for MPI_GATHERV
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: psend
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: precv
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: psend1
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: precv1
  INTEGER, DIMENSION(:), ALLOCATABLE :: ircnt
  INTEGER, DIMENSION(:), ALLOCATABLE :: idisp
  INTEGER, DIMENSION(:), ALLOCATABLE :: ircnt1
  INTEGER, DIMENSION(:), ALLOCATABLE :: idisp1
  ! bubble size
  REAL(KIND(0.D0)) :: R0ref
  REAL(KIND(0.D0)) :: sd
  ! initial void fraction
  REAL(KIND(0.D0)) :: vf0
  ! steady shock speed & shock Mach number
  REAL(KIND(0.D0)) :: Us, Ms
  ! new liquid pressure deviated from pl0
  REAL(KIND(0.D0)) :: plnew
  ! max. CFL
  REAL(KIND(0.D0)) :: CFLmx
  ! UNDEX & FSI properties
  REAL(KIND(0.D0)) :: xs, ps, fsi
  ! no. of flow files
  INTEGER :: Nout
  ! measure computation time
  REAL(KIND(0.D0)) :: elp1, elp2
  ! label for outputing bubble sizes
  INTEGER, DIMENSION(:), ALLOCATABLE :: iR1
  ! no. of bubble sizes for output
  INTEGER :: NR0out

  CONTAINS

  !========================================================================

  SUBROUTINE s_parameter

    !!! computational conditions for Euler equations !!!
    timescheme = 'rk3tvd' ! euler, rk3tvd
    timesplit = 'unsplit' ! unsplit, godunov, strang
    mpweno = 'y'          ! y or n (MP-WENO: Balsara 2000)
    chardecomp = 'y'      ! y or n (characteristic decomposition)
    stretch = 'n'         ! y or n (grid stretching)
    negative = 'y'        ! correct negative void fraction at cell edges
    source = 'y'          ! y or n (static bubbles, zero sources)
    wenoord = 5           ! 1, 3, 5
    wenonum = ( wenoord-1 )/2

    !!! WENO tolerance !!!
    ! 10^{-40} for cavitation cases
    ! Otherwise, 10^{-12} for example.
!    eps = 1.D-40
    eps = 1.D-12

    !!! Accuracy for adaptive RK (godunov, strang) !!!
    accuracy = 1.D-6

    !!! computational conditions for bubble dynamics !!!
    viscous = 'y' ! n for Re = \infty
    polytropic = 'n' ! 'n' for Preston's or Sugiyama's model
    thermal = 'transfer' ! transfer, adiabat, isotherm
    model = 'Preston' ! Preston or Sugiyama
    nquad = 's' ! s (Simpson) or g (Gauss-Hermite)

    !!! liquid and gas !!!
    liquid = 'water' ! water, silicone, glycerol
    gas = 'air' ! air, nitrogen, SF6
    vapor = 'y'

    !!! lognormal bubble size distribution !!!
    R0ref = 50.D-6
    sd = 0.D0
    disperse = 'mono' ! mono or poly

    !!! initial void fraction at pl0 !!!
    ! 10^{-6} to 10^{-5} for distilled water
    vf0 = 0.5D-2

    !!! ICs !!!
    ! for linear computation
    ! 1: narrow pressure perturbation
    ! 2: bubble screen
    ! 3: Leroy et al. (2008), sample 3
    ! 4: narrow pressure perturbation (convergence test)
    ! for nonlinear computation
    ! 1: steady shock, user-defined
    ! 2: Kameda (1998), nitrogen/silicone
    ! 3: Kameda (1998), SF6/silicone
    ! 4: Beylich (1990), SF6/glycerine
    ! 5: bubble screen
    ! 6: cavitation tube
    ! 7: UNDEX
    nonlin = 'y'
    smoothic = 'n'
    ictype = '2'

    !!! no. of cells, domain size & CFL !!!
    Nx_tot = 2500
    length%x = 2500.D0
    CFL = 0.1D0

    !!! BCs (reflect, nonreflect): periodic BC for single-processor ver. !!!
    ! for reflective BCs, solid walls placed at celledges (s_celledgevalue)
    ! but possibly placed at cell centers (see Johnsen 2007)
    ! for nonreflective BCs, Thompson-type BCs (1987) implemented
    xbound%beg = 'nonreflect'
    xbound%end = 'nonreflect'
    ! modification of sonic speed at boundaries
    ! y for shock computation
    ! n for cavitation computation
    modifycl = 'y'

    !!! BCs (freeplate, undex) for FSI & UNDEX !!!
    ! free plate (n, freeplate) at the left boundary
    xincoming%beg = 'n'
    !xincoming%beg = 'undex'
    IF ( xincoming%beg=='freeplate' ) THEN
        xbound%beg = 'reflect'
    ELSE IF ( xincoming%beg=='undex' ) THEN
        xbound%beg = 'nonreflect'
    END IF

    ! incoming waves (n, undex) at the right boundary
    xincoming%end = 'n'
    !xincoming%end = 'freeplate'
    IF ( xincoming%end=='undex' ) THEN
        xbound%end = 'nonreflect'
    ELSE IF ( xincoming%end=='freeplate' ) THEN
        xbound%end = 'reflect'
    END IF


    !!! probe measurement !!!
    timehis = 'y'

    !!! no. of flow files (excluding IC) !!!
    Nout = 100

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

    IF ( nonlin=='y'.AND.ictype=='2' ) THEN
        !fig 6
       liquid = 'silicone'
       gas = 'nitrogen'
       vapor = 'n'
       R0ref = 573.D-6
       vf0 = 0.18D-2
       Nx_tot = 7000
       length%x = 8930.191972D0
       CFL = 0.1D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       smoothic = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='99'.or.ictype=='98' ) THEN
       liquid = 'test_liquid'
       gas = 'air'
       vapor = 'n'
       R0ref = 1.d0 !10.e-6
       vf0 = 0.005d0
       Nx_tot = 300 !601
       length%x = 1.d0
       cfl = 0.1d0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       timehis = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='3' ) THEN
        !fig 7
       liquid = 'silicone'
       gas = 'SF6'
       vapor = 'n'
       R0ref = 613.D-6
       vf0 = 0.24D-2
       Nx_tot = 7000
       length%x = 8347.5D0
       CFL = 0.1D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       smoothic = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='4' ) THEN
       liquid = 'glycerol'
       gas = 'SF6'
       vapor = 'n'
       R0ref = 1.15D-3
       vf0 = 0.25D-2
       Nx_tot = 1200
       length%x = 2543.5D0
       CFL = 0.8D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       smoothic = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='5' ) THEN
        !bubble screen
       liquid = 'water'
       gas = 'air'
       vapor = 'y'
       R0ref = 50.D-6
       vf0 = 0.5D-2
!       Nx_tot = 500
!       length%x = 500.D0
       Nx_tot = 800
       length%x = 800.D0
       CFL = 0.1D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
       smoothic = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='6' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'y'
       R0ref = 50.D-6
       vf0 = 1.D-6
!       Nx_tot = 1200
       Nx_tot = 6000
       length%x = 6000.D0
       CFL = 0.1D0
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
       smoothic = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='7' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'y'
       R0ref = 50.D-6
       vf0 = 1.D-5
       Nx_tot = 3000
       length%x = 18000.D0
!       Nx_tot = 5000
!       length%x = 25000.D0
!       Nx_tot = 2000
!       length%x = 8000.D0
       CFL = 0.1D0
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
       xincoming%beg = 'freeplate'
       xincoming%end = 'undex'
       smoothic = 'n'
    ELSE IF ( nonlin=='y'.AND.ictype=='8' ) THEN
        !SHB: 
    stretch = 'y'         ! y or n (grid stretching)        
       liquid = 'water'
       gas = 'air'
       vapor = 'y'
       R0ref = 50.D-6
       vf0 = 1.D-5
!       Nx_tot = 3000
!       length%x = 18000.D0
!       Nx_tot = 5000
!       length%x = 25000.D0
       Nx_tot = 1500
       length%x = 9000.D0
       CFL = 0.05D0
       !xbound%end = 'reflect'
       !xbound%beg = 'nonreflect'
       !xincoming%end = 'freeplate'
       !xincoming%beg = 'undex'
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
       xincoming%beg = 'freeplate'
       xincoming%end = 'undex'
       smoothic = 'n'
    ELSE IF ( nonlin=='n'.AND.ictype=='1' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'n'
       R0ref = 10.D-6
       vf0 = 0.1D-2
       Nx_tot = 5000
       length%x = 5000.D0
       CFL = 0.1D0
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
    ELSE IF ( nonlin=='n'.AND.ictype=='2' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'n'
       R0ref = 10.D-6
       vf0 = 0.1D-2
       Nx_tot = 3000
       !Nx_tot = 1500
       length%x = 3000.D0
       !length%x = 1500.D0
       CFL = 0.2D0
       xbound%beg = 'nonreflect'
       xbound%end = 'nonreflect'
    ELSE IF ( nonlin=='n'.AND.ictype=='3' ) THEN
       liquid = 'water'
       gas = 'air'
       vapor = 'n'
       R0ref = 100.D-6
       sd = 0.2D0
       disperse = 'poly'
       vf0 = 0.91D-2
       Nx_tot = 5000
       length%x = 5000.D0
       CFL = 0.1D0
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
    ELSE IF ( nonlin=='n'.AND.ictype=='4' ) THEN
        !narrow perturbation (convergence)
       liquid = 'water'
       gas = 'air'
       vapor = 'n'
       R0ref = 10.D-6
       vf0 = 0.1D-2
       Nx_tot = 800
       length%x = 800.D0
       CFL = 0.1D0
       xbound%beg = 'reflect'
       xbound%end = 'nonreflect'
    END IF

    ! time-step splitting
    IF ( source=='n' ) timesplit = 'unsplit'

  END SUBROUTINE s_given_setup

  !========================================================================

  SUBROUTINE s_start_inoutvar

    INTEGER :: ir

    ! no. of conserved variables for the Euler equations
    Nveul = 3
    ! no. of bubble-dynamics-related variables
    IF ( polytropic=='y' ) Nb = 2
    IF ( polytropic=='n'.AND.vapor=='n' ) Nb = 3
    IF ( polytropic=='n'.AND.vapor=='y' ) Nb = 4
    ! odd no. of descretization points in R0
    IF ( disperse=='mono' ) NR0 = 1
!    IF ( disperse=='poly' ) NR0 = 101
    IF ( disperse=='poly' ) NR0 = 401
    ! no. of conserved variables for the entire system
    Nv = Nveul + Nb*NR0
    ! used for left & right eigenvectors
    Neig1 = 5
    Neig = Neig1 + Nb*NR0
    ! no. of bubble sizes for output
    IF ( disperse=='poly'.AND.sd==0.7D0 ) THEN
       NR0out = 5
    ELSE
       NR0out = 1
    END IF

    ! allocation
    ALLOCATE( ibub(Nb,NR0) )
    ALLOCATE( ieig(Nb,NR0) )
    ALLOCATE( iR1(NR0out) )
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
    IF ( timesplit/='unsplit' ) ALLOCATE( dttry(NR0) )

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

    REAL(KIND(0.D0)) :: rhol0
    REAL(KIND(0.D0)) :: mul0
    REAL(KIND(0.D0)) :: ss
    REAL(KIND(0.D0)) :: uu
    REAL(KIND(0.D0)), DIMENSION(NR0) :: omega_R0

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
       print*, 'test_liquid'

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
    IF ( nonlin=='y'.AND.ictype=='2' ) pl0 = 114.4D3
    IF ( nonlin=='y'.AND.ictype=='3' ) pl0 = 112.9D3
    IF ( nonlin=='y'.AND.ictype=='4' ) pl0 = 111.0D3
    IF ( xincoming%end=='undex' .or.  xincoming%beg=='undex'  ) pl0 = 120.9D3 ! at 2m from water surface
    IF ( vapor=='n' ) pv = 0.D0
    IF ( viscous=='n' ) mul0 = 0.D0

    IF ( ictype=='99' ) pl0 = 1.D0
    IF ( ictype=='98' ) pl0 = 1.D0

    ! gas
    IF ( gas=='air' ) gamma_b = 1.4D0
    IF ( gas=='nitrogen' ) gamma_b = 1.399D0
    IF ( gas=='SF6' ) gamma_b = 1.09D0
    IF ( thermal=='isotherm' ) gamma_b = 1.D0
    ! characteristic velocity
    uu = DSQRT( pl0/rhol0 )

    !!! nondimensional parameters !!!
    ! initial radius
    IF ( NR0==1 ) THEN
       ! monodisperse
       R0(1) = 1.D0
       weight(1) = 1.D0
       sd = 0.D0
       iR1(1) = 1
    ELSE
       ! polydisperse (lognormal distributions)
       IF ( nquad=='g' ) THEN
          ! Gauss-Hermite quadrature
          CALL s_gauher( NR0 )
       ELSE IF ( nquad=='s' ) THEN
          ! Simpson's rule
          CALL s_simpson( NR0 )
       END IF
       ! find R0 \approx 1
       CALL s_find_iR1
    END IF
    ! cavitation, Weber & Reynolds no. based on R0ref
    Ca = ( pl0-pv )/( rhol0*uu**2 )
    We = rhol0*uu**2*R0ref/ss
    Re_inv = mul0/( rhol0*uu*R0ref )
    B_tait = B_tait/pl0

    !!! normalized vapor & ambient pressure, and sonic speed of liquid !!!
    pv = pv/pl0
    pl0 = 1.D0
    cl0 = DSQRT( n_tait*(pl0+B_tait) )

    !!! gussed timestep for timestep splitting !!!
    IF ( timesplit/='unsplit' ) THEN
       omega_R0 = ( 3.D0*gamma_b*Ca+2.D0*(3.D0*gamma_b-1.D0)/We/R0 )/R0**2
       omega_R0 = DSQRT( omega_R0 )
       dttry = 1.D-3*( 2.D0*pi/omega_R0 )
    END IF

    
    print*, 'Re_inv, We, Ca'
    print*, Re_inv, We, Ca

  END SUBROUTINE s_physprop_poly

  !========================================================================

  SUBROUTINE s_physprop_trans

    INTEGER :: ir
    REAL(KIND(0.D0)) :: rhol0
    REAL(KIND(0.D0)) :: mul0
    REAL(KIND(0.D0)) :: ss
    REAL(KIND(0.D0)) :: uu
    REAL(KIND(0.D0)) :: mu_n, mu_v
    REAL(KIND(0.D0)) :: D_b
    REAL(KIND(0.D0)) :: gamma_n, gamma_v
    REAL(KIND(0.D0)) :: temp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: chi_vw0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: cp_b0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: k_b0
    REAL(KIND(0.D0)), DIMENSION(NR0) :: rho_b0
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
    IF ( nonlin=='y'.AND.ictype=='2' ) pl0 = 114.4D3
    IF ( nonlin=='y'.AND.ictype=='3' ) pl0 = 112.9D3
    IF ( nonlin=='y'.AND.ictype=='4' ) pl0 = 111.0D3
    IF ( xincoming%end=='undex' .or.  xincoming%beg=='undex' ) pl0 = 120.9D3 ! at 2m from water surface
    temp = 293.15D0
    IF ( nonlin=='y'.AND.ictype=='2' ) temp = 298.15D0
    IF ( nonlin=='y'.AND.ictype=='3' ) temp = 298.15D0
    IF ( nonlin=='y'.AND.ictype=='4' ) temp = 298.D0
    IF ( vapor=='n' ) pv = 0.D0
    IF ( viscous=='n' ) mul0 = 0.D0
    ! gas
    IF ( gas=='air' ) THEN
       gamma_n = 1.4D0
       gamma_b = gamma_n
       M_n = 28.97D0
       mu_n = 1.8D-5
       k_n = 0.02556D0
    ELSE IF ( gas=='nitrogen' ) THEN
       gamma_n = 1.399D0
       gamma_b = gamma_n
       M_n = 28.013D0
       mu_n = 1.76D-5
       k_n = 0.02534D0
    ELSE IF ( gas=='SF6' ) THEN
       gamma_n = 1.09D0
       gamma_b = gamma_n
       M_n = 146.06D0
       mu_n = 1.76D-5
       k_n = 14.D-3
    END IF
    IF ( thermal=='isotherm' ) gamma_b = 1.D0
    ! mass diffusion
    IF ( gas=='air'.AND.liquid=='water' ) THEN
       D_b = 0.242D-4
    ELSE
       ! temporary value
       D_b = 0.242D-4
    END IF
    ! characteristic velocity
    uu = DSQRT( pl0/rhol0 )

    !!! nondimensional parameters !!!
    ! initial radius
    IF ( NR0==1 ) THEN
       ! monodisperse
       R0(1) = 1.D0
       weight(1) = 1.D0
       sd = 0.D0
       iR1(1) = 1
    ELSE
       ! polydisperse (lognormal distributions)
       IF ( nquad=='g' ) THEN
          ! Gauss-Hermite quadrature
          CALL s_gauher( NR0 )
       ELSE IF ( nquad=='s' ) THEN
          ! Simpson's rule
          CALL s_simpson( NR0 )
       END IF
       ! find R0 \approx 1
       CALL s_find_iR1
    END IF
    ! cavitation, Weber & Reynolds no. based on R0ref
    Ca = ( pl0-pv )/( rhol0*uu**2 )
    We = rhol0*uu**2*R0ref/ss
    Re_inv = mul0/( rhol0*uu*R0ref )
    B_tait = B_tait/pl0

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
       cp_b0 = R_n*gamma_n/( gamma_n-1.D0 )
       ! mole fraction of vapor
       x_vw = 0.D0
       ! thermal conductivity for gas
       k_b0 = k_n
       ! gas density
       rho_b0 = pb0/( R_n*temp )
    ELSE
       ! mass fraction of vapor
       chi_vw0 = 1.D0/( 1.D0+R_v/R_n*(pb0/pv-1.D0) )
       print*, chi_vw0, R_v, gamma_v, R_n, gamma_n
       ! specific heat for gas/vapor mixture
       cp_b0 = chi_vw0*R_v*gamma_v/( gamma_v-1.D0 ) &
             + ( 1.D0-chi_vw0 )*R_n*gamma_n/( gamma_n-1.D0 )
       ! mole fraction of vapor
       x_vw = M_n*chi_vw0/( M_v+(M_n-M_v)*chi_vw0 )
       ! thermal conductivity for gas/vapor mixture
       k_b0 = x_vw*k_v/( x_vw+(1.D0-x_vw)*phi_vn ) &
            + ( 1.D0-x_vw )*k_n/( x_vw*phi_nv+1.D0-x_vw )
       ! mixture density
       rho_b0 = pv/( chi_vw0*R_v*temp )
    END IF
    ! mass of gas/vapor computed using dimensional quantities
    mass_n0 = 4.D0*( pb0-pv )*pi/( 3.D0*R_n*temp*rhol0 )*R0**3
    mass_v0 = 4.D0*pv*pi/( 3.D0*R_v*temp*rhol0 )*R0**3
    ! Peclet numbers
    Pe_T = rho_b0*cp_b0*uu*R0ref/k_b0
    Pe_c = uu*R0ref/D_b

    print*, 'rhom0, cpm0, km0', rho_b0, cp_b0, k_b0
    print*, 'Pe_c, Pe_T', Pe_c, Pe_T

    ! nondimensional properties
    R_n = rhol0*R_n*temp/pl0
    R_v = rhol0*R_v*temp/pl0
    k_n = k_n/k_b0
    k_v = k_v/k_b0
    pb0 = pb0/pl0
    pv = pv/pl0
    ! bubble wall temperature, normalized by T0, in the liquid
    ! keeps a constant (cold liquid assumption)
    Tw = 1.D0
    ! natural frequencies
    omegaN = DSQRT( 3.D0*k_poly*Ca+2.D0*(3.D0*k_poly-1.D0)/(We*R0) )/R0

    !!! normalized ambient pressure & sonic speed of liquid !!!
    pl0 = 1.D0
    cl0 = DSQRT( n_tait*(pl0+B_tait) )

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

    !!! gussed timestep for timestep splitting !!!
    IF ( timesplit/='unsplit' ) THEN
       dttry = 1.D-3*( 2.D0*pi/omegaN )
    END IF

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
    ELSE IF ( sd==0.2D0 ) THEN
       R0mn = 0.4D0
       R0mx = 3.D0
    ELSE IF ( sd==0.3D0 ) THEN
       R0mn = 0.3D0
       R0mx = 6.D0
    ELSE IF ( sd==0.5D0 ) THEN
!       R0mn = 0.17D0
       R0mn = 0.15D0
       R0mx = 25.D0
    ELSE IF ( sd==0.7D0 ) THEN
       R0mn = 0.1D0
       R0mx = 200.D0
!       R0mn = 0.16D0
!       R0mx = 25.D0
!       R0mn = 0.12D0
!       R0mx = 150.D0
!       IF ( nonlin=='y'.AND.ictype=='7' ) THEN
!          R0mn = 0.5D0
!          R0mx = 50.D0
!       ELSE
!          R0mn = 0.1D0
!          R0mx = 200.D0
!          R0mn = 0.12D0
!          R0mx = 150.D0
!       END IF
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

  END SUBROUTINE s_simpson

  !========================================================================

  SUBROUTINE s_find_iR1

    INTEGER :: ir
    INTEGER :: iout
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: tmpmin
    REAL(KIND(0.D0)) :: R0tmp

    IF ( sd==0.7D0 ) THEN

       DO iout = 1,NR0out
          tmpmin = 1.D9
          IF ( iout==1 ) R0tmp = 0.25D0
          IF ( iout==2 ) R0tmp = 0.5D0
!          IF ( iout==1 ) R0tmp = 0.5D0
!          IF ( iout==2 ) R0tmp = 0.75D0
          IF ( iout==3 ) R0tmp = 1.D0
          IF ( iout==4 ) R0tmp = 2.D0
          IF ( iout==5 ) R0tmp = 4.D0
          ! find the closest point
          DO ir = 1,NR0
             tmp = DABS( R0(ir)-R0tmp )
             IF ( tmp<tmpmin ) THEN
                tmpmin = tmp
                iR1(iout) = ir
             END IF
          END DO
       END DO

    ELSE

       tmpmin = 1.D9
       R0tmp = 1.D0
       DO ir = 1,NR0
          tmp = DABS( R0(ir)-R0tmp )
          IF ( tmp<tmpmin ) THEN
             tmpmin = tmp
             iR1(1) = ir
          END IF
       END DO

    END IF

  END SUBROUTINE s_find_iR1

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

    ! initial partial pressure for non-condensible gas
    IF ( polytropic=='n' ) pni = pb0(iR0) - pv
    IF ( polytropic=='y' ) pni = Ca + 2.D0/We/R0(iR0)

    ! returns both the function value and the first derivative
    IF ( thermal=='adiabat' ) THEN
       power = 3.D0*gamma_b
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

    ! for entire domain
    INTEGER :: i
    INTEGER :: j
    INTEGER :: iv
    INTEGER :: iout
    INTEGER :: impi
    INTEGER :: ibeg
    INTEGER :: ixbeg
    INTEGER :: iend

    ! uniform grid for entire domain
    ALLOCATE( dxs_tot(1:Nx_tot) )
    ALLOCATE( xgrid_tot(1:Nx_tot) )
    ALLOCATE( xhalf_tot(0:Nx_tot) )
    dx = length%x/DBLE( Nx_tot )
    dxs_tot = dx
    xhalf_tot(0) = 0.D0
    DO i = 1,Nx_tot
       xgrid_tot(i) = ( DBLE(i)-0.5D0 )*dx
       xhalf_tot(i) = DBLE( i )*dx
    END DO

    ! uniform grid for each process
    ! Nx is no. of cells including overlaps
    IF ( xbound%beg=='reflect'.AND.mpi_rank==0 ) THEN
       ibeg = 1 - wenonum
       ixbeg = 0
    ELSE
       ibeg = 1
       ixbeg = 1
    END IF
    IF ( xbound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN
       iend = Nx + wenonum
    ELSE
       iend = Nx
    END IF
    ALLOCATE( dxs(ibeg:iend) )
    ALLOCATE( xgrid(1:Nx) )
    dxs = dx
    DO i = 1,Nx
       xgrid(i) = xgrid_tot(i+mpi_mn-1)
    END DO

    ! allocate variables, WENO coefficients
    ALLOCATE( rhs(Nv),prim(Nv),cons(Nv) )
    DO iv = 1,Nv
       ALLOCATE( rhs(iv)%f(1:Nx) )
       ALLOCATE( prim(iv)%f(1:Nx) )
       ALLOCATE( cons(iv)%f(1:Nx) )
    END DO
    ALLOCATE( prim_tot(Nveul) )
    DO iv = 1,Nveul
       ALLOCATE( prim_tot(iv)%f(1:Nx_tot) )
    END DO
    ALLOCATE( rad_tot(NR0out) )
    DO iout = 1,NR0out
       ALLOCATE( rad_tot(iout)%f(1:Nx_tot) )
    END DO
    IF ( timesplit/='unsplit' ) ALLOCATE( dpldt(1:Nx) )
    ALLOCATE( ix(0:Nx,-2:2) )
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
    ! WENO3 used for nonreflective BCs
    ! only for mpi_rank==0,mpi_size-1
    ALLOCATE( w3betax(0:1,1:Nx) )
    ALLOCATE( w3drx(0:1,1:Nx) )
    ALLOCATE( w3dlx(0:1,1:Nx) )
    ALLOCATE( w3polyrx(0:1,1:Nx) )
    ALLOCATE( w3polylx(0:1,1:Nx) )
    ALLOCATE( onesidebegx(1:5) )
    ALLOCATE( onesideendx(1:5) )
    ! used for MPI_GATHERV
    ALLOCATE( psend(1:Nveul*Nx_nopad) )
    ALLOCATE( precv(1:Nveul*Nx_tot) )
    ALLOCATE( psend1(1:NR0out*Nx_nopad) )
    ALLOCATE( precv1(1:NR0out*Nx_tot) )
    ALLOCATE( ircnt(0:mpi_size-1) )
    ALLOCATE( idisp(0:mpi_size-1) )
    ALLOCATE( ircnt1(0:mpi_size-1) )
    ALLOCATE( idisp1(0:mpi_size-1) )
    ALLOCATE( i_Nx_nopad(0:mpi_size-1) )
    CALL MPI_GATHER( Nx_nopad,1,MPI_INTEGER, &
                     i_Nx_nopad,1,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err )
    IF ( mpi_rank==0 ) THEN
       idisp(0) = 0
       DO impi = 0,mpi_size-1
          ircnt(impi) = Nveul*i_Nx_nopad(impi)
          IF ( impi==mpi_size-1 ) EXIT
          idisp(impi+1) = idisp(impi) + ircnt(impi)
       END DO
    END IF
    CALL MPI_BCAST( ircnt,mpi_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err )
    CALL MPI_BCAST( idisp,mpi_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err )
    IF ( mpi_rank==0 ) THEN
       idisp1(0) = 0
       DO impi = 0,mpi_size-1
          ircnt1(impi) = NR0out*i_Nx_nopad(impi)
          IF ( impi==mpi_size-1 ) EXIT
          idisp1(impi+1) = idisp1(impi) + ircnt1(impi)
       END DO
    END IF
    CALL MPI_BCAST( ircnt1,mpi_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err )
    CALL MPI_BCAST( idisp1,mpi_size,MPI_INTEGER,0,MPI_COMM_WORLD,mpi_err )

    ! index matrix with appropriate boundary treatment
    ! for WENO5, five stincils (-2,-1,0,1,2) used
    DO i = ixbeg,Nx
       DO j = -2,2
          ix(i,j) = i + j
       END DO
    END DO
    IF ( mpi_rank==0 ) THEN
       IF ( xbound%beg=='nonreflect' ) THEN
          ix(2,-2) = 1
          ix(1,-2) = 1
          ix(1,-1) = 1
       END IF
       IF ( mpi_size/=1 ) THEN
          ix(Nx-1,2) = Nx
          ix(Nx,2) = Nx
          ix(Nx,1) = Nx
       END IF
    END IF
    IF ( mpi_rank==mpi_size-1 ) THEN
       IF ( mpi_size/=1 ) THEN
          ix(2,-2) = 1
          ix(1,-2) = 1
          ix(1,-1) = 1
       END IF
       IF ( xbound%end=='nonreflect' ) THEN
          ix(Nx-1,2) = Nx
          ix(Nx,2) = Nx
          ix(Nx,1) = Nx
       END IF
    END IF
    IF ( mpi_rank/=0.AND.mpi_rank/=mpi_size-1 ) THEN
       ix(2,-2) = 1
       ix(1,-2) = 1
       ix(1,-1) = 1
       ix(Nx-1,2) = Nx
       ix(Nx,2) = Nx
       ix(Nx,1) = Nx
    END IF

  END SUBROUTINE s_gridgen

  !========================================================================

  SUBROUTINE s_stretch_grid

    INTEGER :: i
    REAL(KIND(0.D0)) :: c1, c2
    REAL(KIND(0.D0)) :: aa, bb, cc
    REAL(KIND(0.D0)) :: xs
    REAL(KIND(0.D0)) :: xlength

    ! constants for stretching
    IF ( xincoming%end=='undex' ) THEN
       c1 = 0.84D0
       c2 = 0.14D0
       cc = 4.D0
!       c1 = 0.75D0
!       c2 = 0.15D0
!       cc = 5.D0
!       c1 = 0.8D0
!       c2 = 0.2D0
!       cc = 5.D0
       xlength = length%x
       ! position of symmetry
       xs = 0.D0
    else if (xincoming%beg=='undex') then
       c1 = 0.84D0
       c2 = 0.14D0
       cc = 4.D0
!       c1 = 0.75D0
!       c2 = 0.15D0
!       cc = 5.D0
!       c1 = 0.8D0
!       c2 = 0.2D0
!       cc = 5.D0
       xlength = length%x
       ! position of symmetry
       xs = 9000.D0
    END IF
    bb = DLOG( DCOSH(cc*(1.D0-c1)) ) + DLOG( DCOSH(-cc*(1.D0+c1)) ) &
       - 2.D0*DLOG( DCOSH(-c1*cc) )
    bb = bb/cc/( 1.D0-c2/c1 )
    aa = c2/c1 - 2.D0*DLOG( DCOSH(-c1*cc) )/cc/bb

    ! stretching
    ! Note that xhalf_tot(0) = 0
    DO i = 1,Nx_tot
       ! cell edge
       IF ( xgrid_tot(i)<xs ) THEN
          xhalf_tot(i) = &
                   - f_newgrid( xs-xhalf_tot(i),aa,bb,cc,c1,xlength ) + xs
       ELSE
          xhalf_tot(i) = &
                     f_newgrid( xhalf_tot(i)-xs,aa,bb,cc,c1,xlength ) + xs
       END IF
       ! cell center
       xgrid_tot(i) = 0.5D0*( xhalf_tot(i)+xhalf_tot(i-1) )
       ! cell width
       dxs_tot(i) = xhalf_tot(i) - xhalf_tot(i-1)
    END DO
    ! substitution to each processor's
    DO i = 1,Nx
       xgrid(i) = xgrid_tot(i+mpi_mn-1)
       dxs(i) = dxs_tot(i+mpi_mn-1)
    END DO
    ! for reflective BC
    IF ( xbound%beg=='reflect'.AND.mpi_rank==0 ) THEN
       IF ( wenoord/=1 ) THEN
          dxs(0) = dxs(1)
          IF ( wenoord==5 ) dxs(-1) = dxs(2)
       END IF
    END IF
    IF ( xbound%end=='reflect'.AND.mpi_rank==mpi_size-1 ) THEN
       IF ( wenoord/=1 ) THEN
          dxs(Nx+1) = dxs(Nx)
          IF ( wenoord==5 ) dxs(Nx+2) = dxs(Nx-1)
       END IF
    END IF

  END SUBROUTINE s_stretch_grid

  !========================================================================

  FUNCTION f_newgrid( pos,aa,bb,cc,c1,domlen )

    REAL(KIND(0.D0)), INTENT(IN) :: pos
    REAL(KIND(0.D0)), INTENT(IN) :: aa, bb, cc, c1
    REAL(KIND(0.D0)), INTENT(IN) :: domlen
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: f_newgrid

    ! smooth function for grid stretching
    tmp = DLOG( DCOSH(cc/domlen*(pos-c1*domlen)) ) &
        + DLOG( DCOSH(cc/domlen*(pos+c1*domlen)) )
    f_newgrid = pos*( aa+tmp/bb/cc )

  END FUNCTION f_newgrid

  !========================================================================

  SUBROUTINE s_probe_setup

    INTEGER :: i
    INTEGER :: i_probe1_tot
    INTEGER :: i_probe2_tot
    INTEGER :: i_integ1_tot
    INTEGER :: i_integ2_tot

    ! no probe measurement (default)
    mpi_rank_probe1 = mpi_size
    mpi_rank_probe2 = mpi_size
    mpi_rank_integ1 = mpi_size
    mpi_rank_integ2 = mpi_size

    ! set probe locations
    IF ( timehis=='y'.AND.nonlin=='y' ) THEN

       IF ( ictype=='1' ) THEN
          i_integ1_tot = 101
          i_integ2_tot = 151
          DO i = i_nopad_mn,i_nopad_mx
             IF ( xgrid(i)==xgrid_tot(i_integ1_tot) ) THEN
                mpi_rank_integ1 = mpi_rank
                i_integ1 = i
             ELSE IF ( xgrid(i)==xgrid_tot(i_integ2_tot) ) THEN
                mpi_rank_integ2 = mpi_rank
                i_integ2 = i
             END IF
          END DO
       ELSE
          IF ( ictype=='2' ) THEN
             i_probe1_tot = 3*Nx_tot/7
             i_probe2_tot = i_probe1_tot + 1
          ELSE IF ( ictype=='3' ) THEN
             i_probe1_tot = 3*Nx_tot/7
             i_probe2_tot = i_probe1_tot + 1
          ELSE IF ( ictype=='4' ) THEN
             i_probe1_tot = Nx_tot/2
             i_probe2_tot = i_probe1_tot + 1
          ELSE IF ( ictype=='5' ) THEN
!             i_probe1_tot = Nx_tot/20
!             i_probe2_tot = 19*Nx_tot/20
             i_probe1_tot = Nx_tot/4 - 10
             i_probe2_tot = 3*Nx_tot/4 + 11
          ELSE IF ( ictype=='6' ) THEN
             i_probe1_tot = 1
             i_probe2_tot = Nx_tot/2
          END IF
          DO i = i_nopad_mn,i_nopad_mx
             IF ( xgrid(i)==xgrid_tot(i_probe1_tot) ) THEN
                mpi_rank_probe1 = mpi_rank
                i_probe1 = i
             ELSE IF ( xgrid(i)==xgrid_tot(i_probe2_tot) ) THEN
                mpi_rank_probe2 = mpi_rank
                i_probe2 = i
             END IF
          END DO
       END IF

    ELSE IF ( timehis=='y'.AND.nonlin=='n' ) THEN

       IF ( ictype=='1'.OR.ictype=='3' ) THEN
          i_probe1_tot = NINT( 25.D0/dx )
          i_probe2_tot = i_probe1_tot + 1
          i_integ1_tot = NINT( Nx_tot/10.D0 )
          i_integ2_tot = NINT( Nx_tot*111.D0/250.D0 )
          DO i = i_nopad_mn,i_nopad_mx
             IF ( xgrid(i)==xgrid_tot(i_probe1_tot) ) THEN
                mpi_rank_probe1 = mpi_rank
                i_probe1 = i
             ELSE IF ( xgrid(i)==xgrid_tot(i_probe2_tot) ) THEN
                mpi_rank_probe2 = mpi_rank
                i_probe2 = i
             ELSE IF ( xgrid(i)==xgrid_tot(i_integ1_tot) ) THEN
                mpi_rank_integ1 = mpi_rank
                i_integ1 = i
             ELSE IF ( xgrid(i)==xgrid_tot(i_integ2_tot) ) THEN
                mpi_rank_integ2 = mpi_rank
                i_integ2 = i
             END IF
          END DO
       ELSE IF ( ictype=='2' ) THEN
          i_probe1_tot = Nx_tot*1/6
          i_probe2_tot = Nx_tot*5/6
          DO i = i_nopad_mn,i_nopad_mx
             IF ( xgrid(i)==xgrid_tot(i_probe1_tot) ) THEN
                mpi_rank_probe1 = mpi_rank
                i_probe1 = i
             ELSE IF ( xgrid(i)==xgrid_tot(i_probe2_tot) ) THEN
                mpi_rank_probe2 = mpi_rank
                i_probe2 = i
             END IF
          END DO
       ELSE IF ( ictype=='4' ) THEN
          mpi_rank_integ1 = 0
          i_integ1 = 1
       END IF

    END IF

  END SUBROUTINE s_probe_setup

  !========================================================================

  SUBROUTINE s_ic

    INTEGER :: i
    INTEGER :: ir
    INTEGER :: iv
    INTEGER :: diaph
    REAL(KIND(0.D0)) :: fx
    REAL(KIND(0.D0)) :: x0
    REAL(KIND(0.D0)) :: small
    REAL(KIND(0.D0)) :: R3bar
    REAL(KIND(0.D0)) :: rhor, rhol
    REAL(KIND(0.D0)) :: vel1r, vel1l
    REAL(KIND(0.D0)) :: nr, nl, n0
    REAL(KIND(0.D0)) :: vfr, vfl
    REAL(KIND(0.D0)) :: cltld_r, cltld_l, cltld_0
    REAL(KIND(0.D0)) :: Smx, cm0, ratio
    REAL(KIND(0.D0)) :: vftiny, ntiny
    REAL(KIND(0.D0)) :: vftmp
    REAL(KIND(0.D0)) :: pltmp, plH
    REAL(KIND(0.D0)), DIMENSION(NR0) :: Rtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pbtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: mvtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pb0tmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: mv0tmp

    ! substitute initial pb0 and mass_v0
    pb0tmp = 0.D0
    mv0tmp = 0.D0
    IF ( polytropic=='n' ) THEN
       pb0tmp = pb0
       mv0tmp = mass_v0
    END IF

    !!! shock (non-smooth IC) !!!
    IF ( nonlin=='y'.AND.smoothic=='n' ) THEN

       IF ( ictype=='1' ) THEN
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state (pl0) !!!
          rhor = 1.D0 - vf0
          vel1r = 0.D0
          nr = n0
          vfr = vf0
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state (pltmp) !!!
          pltmp = 2.D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
          cltld_l = DSQRT( n_tait*(pltmp+B_tait)/(1.D0-vfl)/rhol )
          ! numerical setup
!          diaph = Nx_tot/5
          diaph = Nx_tot/2
          Smx = MAX( DABS(vel1r)+cltld_r,DABS(vel1l)+cltld_l )
!          finaltime = 50.D0
          finaltime = 10.D0
          dt = CFL*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL

      ELSEIF ( ictype=='99' ) THEN
          ! SOS for CFL condition. sonic speed at 0 Hz for pl0
          cm0 = n_tait*( pl0+B_tait )/( 1.D0-vf0 )**2
          cm0 = DSQRT(cm0)
          cm0 = 1.23251662148 !1.18916176545
          ! total bubble number density
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          
          !!! right state associated with pl0 !!!
          pltmp = 1.3*pl0
          nr = n0
          vfr = vf0
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          print*, 'ratio', ratio, vfr, pl0, B_tait, n_tait
          rhor = ratio*(1-vf0)
          vel1r = 0.D0
          
          !!! left state associated with pltmp !!!
          pltmp = 1.5*pl0
          !DO ir = 1,NR0
          !   iR0 = ir
          !   CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          !END DO
          Rtmp(:) = R0; pbtmp(:) = 0.; mvtmp(:) = 0.
          CALL s_quad( Rtmp**3,R3bar )
          ratio = ( (pltmp+B_tait)/(pl0+B_tait) )**( 1.D0/n_tait )
          nl = n0 !*ratio/( 1.D0-vf0+pi43*n0*R3bar*ratio )
          vfl = vf0 !pi43*nl*R3bar
          print*, 'ratio', ratio, vfl
          rhol = ratio*( 1.D0-vfl ) !1.33591698254          
          vel1l = 0.D0
          
          diaph = ceiling(Nx/2.)
          dt = cfl*dx/cm0
          Nt = 500 !NINT( finaltime/dt )
          finaltime = Nt*dt
          interval = 50 !Nt/Nout
          CFLmx = CFL

          print*, 'rho_L/R = ', rhol, rhoR
          !print*, 'p_L/R = ', pltmp, pl0
          print*, 'vf_L/R = ', vfl, vfr
          print*, 'n_L/R = ', nl, nr
          print*, 'vel_L/R = ', vel1l, vel1r
          print*, 'R_L/R = ', Rtmp(1), R0(1)

       ELSE IF ( ictype=='2' ) THEN
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state (pl0) !!!
          rhor = 1.D0 - vf0
          vel1r = -0.012275166D0
          nr = n0
          vfr = vf0
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state (pltmp) !!!
          pltmp = 1.749125874D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
          vel1l = vel1l + vel1r
          cltld_l = DSQRT( n_tait*(pltmp+B_tait)/(1.D0-vfl)/rhol )
          ! numerical setup
          diaph = Nx_tot/7
          finaltime = 300.D0
          Smx = MAX( DABS(vel1r)+cltld_r,DABS(vel1l)+cltld_l )
          dt = CFL*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='3' ) THEN
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state (pl0) !!!
          rhor = 1.D0 - vf0
          vel1r = -0.01236D0
          nr = n0
          vfr = vf0
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state (pltmp) !!!
          pltmp = 2.157D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
          vel1l = vel1l + vel1r
          cltld_l = DSQRT( n_tait*(pltmp+B_tait)/(1.D0-vfl)/rhol )
          ! numerical setup
          diaph = Nx_tot/7
          finaltime = 300.D0
          Smx = MAX( DABS(vel1r)+cltld_r,DABS(vel1l)+cltld_l )
          dt = CFL*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='4' ) THEN
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state (pl0) !!!
          rhor = 1.D0 - vf0
          vel1r = 0.D0
          nr = n0
          vfr = vf0
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          !!! left state (pltmp) !!!
          pltmp = 1.81D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
          cltld_l = DSQRT( n_tait*(pltmp+B_tait)/(1.D0-vfl)/rhol )
          ! numerical setup
          diaph = Nx_tot/6
          finaltime = 100.D0
          Smx = MAX( DABS(vel1r)+cltld_r,DABS(vel1l)+cltld_l )
          dt = CFL*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='5' ) THEN
          !for bubble screen
          
          ! bubble number densities at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          vftiny = 1.D-7
          ntiny = ( vftiny/vf0 )*n0
          
          !!! right state (pl0) !!!
          rhor = 1.D0 - vftiny
          vel1r = 0.D0
          nr = ntiny
          vfr = vftiny
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vftiny )
          
          !!! left state !!!
          pltmp = 5.D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_steadyshock( &
               pltmp,Rtmp,vftiny,ntiny,rhol,vel1l,vfl,nl,Us,Ms )
          cltld_l = DSQRT( n_tait*(pltmp+B_tait)/(1.D0-vfl)/rhol )

            print*, 'pl', (pl0+B_tait)*(rhol/(1.d0-vfl))**n_tait - B_tait
            print*, 'vfl, R0l, vl, cl', vfl, Rtmp, vel1l !, cml
            
            print*, 'pr', (pl0+B_tait)*(rhor/(1.d0-vfr))**n_tait - B_tait
            print*, 'vfr, R0r, vr, cr', vfr, R0, vel1r !, cmr
            
            print*, 'screen part'
            print*, 'p', (pl0+B_tait)*((1.d0-vf0)/(1.d0-vf0))**n_tait - B_tait
            print*, 'vf, R0, v, c', vf0, R0, 0d0

          ! numerical setup
!          diaph = 2*Nx_tot/25
          diaph = Nx_tot/8
          Smx = MAX( DABS(vel1r)+cltld_r,DABS(vel1l)+cltld_l )

            print*, 'smx', Smx

!          finaltime = 30.D0*length%x/Smx
          finaltime = 150.D0
          dt = CFL*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
       ELSE IF ( ictype=='6' ) THEN
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          !!! right state (pl0) !!!
          rhor = 1.D0 - vf0
          nr = n0
          vfr = vf0
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          vel1r = 0.5D0/DSQRT( 101325.D0/998.2063D0 )
          ! numerical setup
          diaph = Nx_tot ! irrelevant
          Smx = DABS( vel1r ) + cltld_r
          finaltime = length%x/Smx
          dt = CFL*dx/Smx
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          CFLmx = CFL
          ! for moving wall BC
          fsi = 0.D0
          mp_inv = 0.D0


          
       END IF
       ! assign ICs
       IF ( ictype=='6' ) THEN
          DO i = 1,Nx
             CALL s_assignIC( rhor,vel1r,vfr,R0,pb0tmp,mv0tmp,nr,i )
          END DO
       ELSE
          DO i = 1,Nx
             IF ( xgrid(i)<=xgrid_tot(diaph) ) THEN
                CALL s_assignIC( rhol,vel1l,vfl,Rtmp,pbtmp,mvtmp,nl,i )
             ELSE
                CALL s_assignIC( rhor,vel1r,vfr,R0,pb0tmp,mv0tmp,nr,i )
             END If
          END DO
       END IF
       ! bubble screen
       IF ( ictype=='5' ) THEN
          DO i = 1,Nx
!             IF ( xgrid(i)>=xgrid_tot(Nx_tot/10+1).AND. &
!                  xgrid(i)<=xgrid_tot(9*Nx_tot/10) ) THEN
             IF ( xgrid(i)>=xgrid_tot(Nx_tot/4+1).AND. &
                  xgrid(i)<=xgrid_tot(3*Nx_tot/4) ) THEN
                CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
             END IF
          END DO
       END IF

    !!! shock (smooth IC) !!!
    ELSE IF ( nonlin=='y'.AND.smoothic=='y' ) THEN

       IF ( ictype=='1' ) THEN
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          ! initial state
          rhor = 1.D0 - vf0
          vel1r = 0.D0
          nr = n0
          vfr = vf0
          cltld_r = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          ! shock pressure
          plH = 500.D0*pl0
          DO ir = 1,NR0
             iR0 = ir
             CALL s_newequilibrium( plH,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
          END DO
          CALL s_steadyshock( plH,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
          cltld_l = DSQRT( n_tait*(plH+B_tait)/(1.D0-vfl)/rhol )
          Smx = MAX( DABS(vel1r)+cltld_r,DABS(vel1l)+cltld_l )
          ! diaphragm
          diaph = Nx_tot/10
          x0 = 0.5D0*( xgrid_tot(diaph)+xgrid_tot(diaph+1) )
       END IF
       ! using approximate Heaviside function
       DO i = 1,Nx
          CALL s_smoothshock( xgrid(i),x0,dx,fx )
          pltmp = fx*( plH-pl0 ) + pl0
          IF ( fx<1.D-12 ) THEN
             CALL s_assignIC( rhor,vel1r,vfr,R0,pb0tmp,mv0tmp,nr,i )
          ELSE
             DO ir = 1,NR0
                iR0 = ir
                CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
             END DO
             CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
             CALL s_assignIC( rhol,vel1l,vfl,Rtmp,pbtmp,mvtmp,nl,i )
          END IF
       END DO
       finaltime = 0.75D0*length%x/Smx
       dt = CFL*dx/Smx
       Nt = NINT( finaltime/dt )
       interval = Nt/Nout
       CFLmx = CFL

    !!! linear problem !!!
    ELSE IF ( nonlin=='n' ) THEN

       IF ( ictype=='1'.OR.ictype=='3' ) THEN
          ! perturbation in liquid pressure
          small = 1.D-4
          ! sonic speed at pl0
          cltld_0 = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          ! purturbed quantities
          DO i = 1,Nx
             fx = DEXP( -(xgrid(i)/4.D0)**2 )
             pltmp = pl0*( 1.D0+small*fx )
             DO ir = 1,NR0
                iR0 = ir
                CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
             END DO
             CALL s_steadyshock( &
                  pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
             CALL s_assignIC( rhol,0.D0,vfl,Rtmp,pbtmp,mvtmp,nl,i )
          END DO
          finaltime = 1.5D0*length%x/cltld_0
          dt = CFL*dx/cltld_0
          Nt = NINT( finaltime/dt )
          interval = Nt/Nout
          print*, 'Interval = ', interval
          CFLmx = CFL
       ELSE IF ( ictype=='2' ) THEN
          ! right-going wave
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          small = 1.D-4
          vftiny = 1.D-7
          ntiny = ( vftiny/vf0 )*n0
          cltld_0 = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vftiny )
          DO i = 1,Nx
             fx = xgrid(i) - 0.5*( &
                  xgrid_tot(29*Nx_tot/60)+xgrid_tot(29*Nx_tot/60+1) )
             fx = DEXP( -(fx/8.D0)**2 )
             pltmp = pl0*( 1.D0+small*fx )
             DO ir = 1,NR0
                iR0 = ir
                CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
             END DO
             CALL s_steadyshock( &
                  pltmp,Rtmp,vftiny,ntiny,rhol,vel1l,vfl,nl,Us,Ms )
             vel1l = cltld_0*( rhol-1.D0+vftiny )
             CALL s_assignIC( rhol,vel1l,vfl,Rtmp,pbtmp,mvtmp,nl,i )
          END DO
          ! bubble screen
          DO i = 1,Nx
             IF ( xgrid(i)>=xgrid_tot(Nx_tot/2+1).AND. &
                  xgrid(i)<=xgrid_tot(2*Nx_tot/3) ) THEN
                CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
             END IF
          END DO
          ! redifine sonic speed
          cltld_0 = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          finaltime = 0.5D0*length%x/cltld_0
          dt = CFL*dx/cltld_0
          Nt = NINT( finaltime/dt )
          Nt = 562500
          finaltime=dt*Nt
          interval = Nt/Nout
          CFLmx = CFL
          print*,'final time = ', finaltime
          print*, 'Nt' , Nt, 'dt = ', dt
       ELSE IF ( ictype=='4' ) THEN
          ! perturbation in liquid pressure
          small = 1.D-4
          ! sonic speed at pl0
          cltld_0 = DSQRT( n_tait*(pl0+B_tait) )/( 1.D0-vf0 )
          ! bubble number density at pl0
          CALL s_quad( R0**3,R3bar )
          n0 = vf0/( pi43*R3bar )
          ! purturbed quantities
          DO i = 1,Nx
             fx = DEXP( -(xgrid(i)/4.D0)**2 )
             pltmp = pl0*( 1.D0+small*fx )
             DO ir = 1,NR0
                iR0 = ir
                CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
             END DO
             CALL s_steadyshock( &
                  pltmp,Rtmp,vf0,n0,rhol,vel1l,vfl,nl,Us,Ms )
             CALL s_assignIC( rhol,0.D0,vfl,Rtmp,pbtmp,mvtmp,nl,i )
          END DO
          finaltime = 10.D0
          Nt = 15000
          interval = Nt/Nout
          dt = finaltime/DBLE( Nt )
          CFLmx = dt*cltld_0/dx
       END IF

    END IF

    ! time-step splitting
    IF ( timesplit/='unsplit' ) THEN
       ! splitting time step
       IF ( timesplit=='godunov' ) dtsplit = dt
       IF ( timesplit=='strang' ) dtsplit = 0.5D0*dt
    END IF

  END SUBROUTINE s_ic

  !========================================================================

  SUBROUTINE s_undexic

    INTEGER :: i
    INTEGER :: ir
    INTEGER :: ixs, iscreen
    REAL(KIND(0.D0)) :: pltmp, xscreen
    REAL(KIND(0.D0)) :: R3bar
    REAL(KIND(0.D0)) :: Smax
    REAL(KIND(0.D0)) :: n0
    REAL(KIND(0.D0)) :: rhos, vel1s, vfs, ns
    REAL(KIND(0.D0)), DIMENSION(NR0) :: Rtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pbtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: mvtmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: pb0tmp
    REAL(KIND(0.D0)), DIMENSION(NR0) :: mv0tmp
    REAL(KIND(0.D0)), PARAMETER :: xdecaydml = 0.02247D0 ! [m]
!    REAL(KIND(0.D0)), PARAMETER :: xdecaydml = 0.02757D0 ! [m]
!    REAL(KIND(0.D0)), PARAMETER :: xdecaydml = 0.03014D0 ! [m]
!    REAL(KIND(0.D0)), PARAMETER :: xdecaydml = 0.03211D0 ! [m]
!    REAL(KIND(0.D0)), PARAMETER :: xdecaydml = 0.03373D0 ! [m]

    if (ictype=='7') then
        ! UNDEX properties
        ixs = 3  !SHB: sets where the shock starts
        xs = xgrid_tot(ixs) + 0.5D0*dxs_tot(ixs)
        ps = 41.2D0*pl0 !SHB: shock strength
        !    ps = 14.61D0*pl0
        !    ps = 9.24D0*pl0
        !    ps = 6.677D0*pl0
        !    ps = 5.189D0*pl0
        xdecay = xdecaydml/R0ref !SHB: sets how long the shock is
    else if (ictype=='8') then
        xscreen = 1000.d0
        iscreen = minloc(abs(xgrid_tot(:) - xscreen),1)
        xs = 1.1d0*xscreen
        ixs = minloc(abs(xgrid_tot(:) - xs),1)
 
        !iscreen = floor(Nx_tot/2)
        !ixs = iscreen + 3  !SHB: sets where the shock starts
        !xs = xgrid_tot(ixs) + 0.5D0*dxs_tot(ixs)
        ps = 41.D0*pl0 !SHB: shock strength?
        xdecay = 0.1*xdecaydml/R0ref !SHB: sets how long the shock is?
        !xdecay = 100. !SHB: sets how long the shock is?
    end if

    ! FSI properties (fsi=0:stationary wall)
    ! fsi    = \psi
    ! mp_inv = 1/m_p 
    ! xdecay = \rho_l c_l \tau
    ! \tau sets the decay rate, essentially

!    fsi = 5.D0
    fsi = 0.72D0
!    fsi = 0.1D0
    mp_inv = fsi/xdecay

    ! substitute initial pb0 and mass_v0
    pb0tmp = 0.D0
    mv0tmp = 0.D0
    IF ( polytropic=='n' ) THEN
       pb0tmp = pb0
       mv0tmp = mass_v0
    END IF

    ! quantities at pl0
    CALL s_quad( R0**3,R3bar )
    n0 = vf0/( pi43*R3bar )

    if (ictype=='7') then
        ! set ic
        DO i = 1,Nx
           IF ( xgrid(i)<xs ) THEN ! undisturbed
              CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
           ELSE ! shocked (from the right)
              pltmp = pl0 + ps*DEXP( (xgrid(i)-xs)/xdecay )
              DO ir = 1,NR0
                 iR0 = ir
                 CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
              END DO
              CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhos,vel1s,vfs,ns,Us,Ms )
              vel1s = -vel1s ! sign correction
              CALL s_assignIC( rhos,vel1s,vfs,Rtmp,pbtmp,mvtmp,ns,i )
           END IF
        END DO
    else if (ictype=='8') then
         ! set ic
        DO i = 1,Nx
           IF ( xgrid(i)<xs ) THEN ! undisturbed
              CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
           ELSE ! shocked (from the left)
              !pltmp = pl0 + ps*DEXP( -1.*(xs-xgrid(i))/xdecay )
              pltmp = pl0 + ps*DEXP( (xs-xgrid(i))/xdecay )
              DO ir = 1,NR0
                 iR0 = ir
                 CALL s_newequilibrium( pltmp,Rtmp(ir),pbtmp(ir),mvtmp(ir) )
              END DO
              CALL s_steadyshock( pltmp,Rtmp,vf0,n0,rhos,vel1s,vfs,ns,Us,Ms )
              vel1s = -vel1s ! sign correction
              CALL s_assignIC( rhos,vel1s,vfs,Rtmp,pbtmp,mvtmp,ns,i )
           END IF
        END DO   
    end if


    ! bubble screen
    IF ( ictype=='8' ) THEN
        vf0 = 1e-3
        n0 = vf0/( pi43*R3bar )
        DO i = 1,Nx
            !IF ( xgrid(i)>=xgrid_tot(3*Nx_tot/4).AND. &
            !     xgrid(i)<=xgrid_tot(Nx_tot) ) THEN
            IF ( xgrid(i)>=xgrid_tot(1).AND. &
                 xgrid(i)<=xgrid_tot(iscreen) ) THEN                 
                 CALL s_assignIC( 1.D0-vf0,0.D0,vf0,R0,pb0tmp,mv0tmp,n0,i )
            END IF
        END DO
    END IF
    
    finaltime = 80.D0 ! fsi=5
!    finaltime = 300.D0 ! fsi=0.72
!    finaltime = 600.D0 ! fsi=0.1
    Smax = ( pl0+ps+B_tait )/( pl0+B_tait )
    Smax = Smax**( 0.5D0*(n_tait-1.D0)/n_tait )
    Smax = ps/cl0 + cl0*Smax
    dt = CFL*MINVAL( dxs_tot )/Smax
    Nt = NINT( finaltime/dt )
    interval = Nt/Nout
    CFLmx = CFL

    ! free plate
    uwall = 0.D0
    pwall = pl0

    ! time-step splitting
    IF ( timesplit/='unsplit' ) THEN
       ! splitting time step
       IF ( timesplit=='godunov' ) dtsplit = dt
       IF ( timesplit=='strang' ) dtsplit = 0.5D0*dt
    END IF

  END SUBROUTINE s_undexic

  !========================================================================

  SUBROUTINE s_steadyshock( plH,RH,vfi,ni,rhoH,velH,vfH,nH,Ushock,Mshock )

    REAL(KIND(0.D0)), INTENT(IN) :: plH, vfi, ni
    REAL(KIND(0.D0)), DIMENSION(NR0), INTENT(IN) :: RH
    REAL(KIND(0.D0)), INTENT(OUT) :: rhoH, velH, vfH, nH, Ushock, Mshock
    REAL(KIND(0.D0)) :: RH3bar, ratio, c0, vfi1, rho0
    REAL(KIND(0.D0)) :: deno

    ! sonic speed of the mixture in a quasistatic limit
    vfi1 = 1.D0 - vfi
    rho0 = vfi1
    c0 = n_tait*( pl0+B_tait )/rho0
    IF ( thermal=='adiabat' ) THEN
       c0 = c0*gamma_b*pl0/( vfi*n_tait*(pl0+B_tait)+vfi1*gamma_b*pl0 )
    ELSE
       c0 = c0*pl0/( vfi*n_tait*(pl0+B_tait)+vfi1*pl0 )
    END IF
    c0 = DSQRT( c0 )

    ! bubble number density
    CALL s_quad( RH**3,RH3bar )
    ratio = ( (pl0+B_tait)/(plH+B_tait) )**( 1.D0/n_tait )
    nH = ni/( vfi1*ratio+pi43*ni*RH3bar )
    ! void fraction
    vfH = pi43*nH*RH3bar
    ! mixture density
    rhoH = ( 1.D0-vfH )/ratio
    deno = 1.D0 - rho0/rhoH
    IF ( deno==0.D0 ) THEN
       ! induced velocity
       velH = 0.D0
       Ushock = c0
       Mshock = 1.D0
    ELSE
       ! steady shock speed
       Ushock = ( plH-pl0 )/rho0/deno
       Ushock = DSQRT( Ushock )
       ! induced velocity
       velH = Ushock*deno
       ! shoch Mach number
       Mshock = Ushock/c0
    END IF

  END SUBROUTINE s_steadyshock

  !========================================================================

  SUBROUTINE s_smoothshock( dist,diaph,del,prop )

    REAL(KIND(0.D0)), INTENT(IN) :: dist
    REAL(KIND(0.D0)), INTENT(IN) :: diaph
    REAL(KIND(0.D0)), INTENT(IN) :: del
    REAL(KIND(0.D0)), INTENT(OUT) :: prop
    REAL(KIND(0.D0)), PARAMETER :: coeff = 1.5D0
    REAL(KIND(0.D0)) :: tmp

    ! approximate Heaviside function
    tmp = TANH( coeff*(dist-diaph)/del )
    prop = 0.5D0*( 1.D0-tmp )

  END SUBROUTINE s_smoothshock

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

  !========================================================================

  SUBROUTINE s_wenocoeff

    ! coefficients computed in [i_nopad_mn,i_nopad_mx]
    INTEGER :: i
    REAL(KIND(0.D0)) :: coeff
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: tmp1
    REAL(KIND(0.D0)) :: tmp2
    REAL(KIND(0.D0)) :: tmp10
    REAL(KIND(0.D0)) :: tmp20
    REAL(KIND(0.D0)) :: m1m5, p1m5, p1m3, p3m3, p3m1, p5m1, p5p1

    IF ( wenoord==5 ) THEN

       ! WENO5
       ! polynomials
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
          dwxr(2,i) = ( dxs(ix(i,-2))+dxs(ix(i,-1))+dxs(i) ) &
                    * (dxs(ix(i,-1))+dxs(i) )/tmp2/tmp
          dwxr(1,i) = 1.D0 - dwxr(2,i) - dwxr(0,i)
          dwxl(0,i) = ( dxs(i)+dxs(ix(i,1)) )*( dxs(i)+dxs(ix(i,1) ) &
                    + dxs(ix(i,2)) )/tmp1/tmp
          dwxl(2,i) = ( dxs(ix(i,-2))+dxs(ix(i,-1)) ) &
                    * dxs(ix(i,-1))/tmp2/tmp
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
          betax2(1,i) = -coeff*tmp*tmp1*( tmp20+dxs(ix(i,1)) &
                      * (dxs(ix(i,2))+3.D0*tmp ) + tmp*(tmp1+dxs(i)) )
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
          betax0(2,i) = coeff*tmp1*tmp1*( tmp10+dxs(ix(i,-1)) &
                      * dxs(ix(i,-1))+dxs(i)*dxs(ix(i,-1)) )
       END DO

    END IF
    IF ( wenoord/=1 ) THEN

       ! WENO3
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

    END IF

    ! fourth-order one-sided difference for WENO5
    ! second-order one-sided difference for WENO3
    ! first-order one-sided difference for WENO1
    IF ( xbound%beg=='nonreflect'.AND.mpi_rank==0 ) THEN

       IF ( wenoord==5 ) THEN
          onesidebegx(1) = -( 1.D0/(xgrid(2)-xgrid(1)) &
                         + 1.D0/(xgrid(3)-xgrid(1)) &
                         + 1.D0/(xgrid(4)-xgrid(1)) &
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
       ELSE IF ( wenoord==3 ) THEN
          onesidebegx(1) = ( 2.D0*xgrid(1)-xgrid(2)-xgrid(3) ) &
                         / ( xgrid(2)-xgrid(1) )/( xgrid(3)-xgrid(1) )
          onesidebegx(2) = ( xgrid(3)-xgrid(1) ) &
                         / ( xgrid(2)-xgrid(1) )/( xgrid(3)-xgrid(2) )
          onesidebegx(3) = ( xgrid(1)-xgrid(2) ) &
                         / ( xgrid(3)-xgrid(1) )/( xgrid(3)-xgrid(2) )
          onesidebegx(4) = 0.D0
          onesidebegx(5) = 0.D0
       ELSE IF ( wenoord==1 ) THEN
          onesidebegx(2) = 1.D0/( xgrid(2)-xgrid(1) )
          onesidebegx(1) = -onesidebegx(2)
          onesidebegx(3) = 0.D0
          onesidebegx(4) = 0.D0
          onesidebegx(5) = 0.D0
       END IF

    END IF
    IF ( xbound%end=='nonreflect'.AND.mpi_rank==mpi_size-1 ) THEN

       IF ( wenoord==5 ) THEN
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
       ELSE IF ( wenoord==3 ) THEN
          onesideendx(1) = ( 2.D0*xgrid(Nx)-xgrid(Nx-1)-xgrid(Nx-2) ) &
                         / ( xgrid(Nx)-xgrid(Nx-1) ) &
                         / ( xgrid(Nx)-xgrid(Nx-2) )
          onesideendx(2) = ( xgrid(Nx)-xgrid(Nx-2) ) &
                         / ( xgrid(Nx)-xgrid(Nx-1) ) &
                         / ( xgrid(Nx-2)-xgrid(Nx-1) )
          onesideendx(3) = ( xgrid(Nx-1)-xgrid(Nx) ) &
                         / ( xgrid(Nx)-xgrid(Nx-2) ) &
                         / ( xgrid(Nx-2)-xgrid(Nx-1) )
          onesideendx(4) = 0.D0
          onesideendx(5) = 0.D0
       ELSE IF ( wenoord==1 ) THEN
          onesideendx(1) = 1.D0/( xgrid(Nx)-xgrid(Nx-1) )
          onesideendx(2) = -onesideendx(1)
          onesideendx(3) = 0.D0
          onesideendx(4) = 0.D0
          onesideendx(5) = 0.D0
       END IF

    END IF

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

  SUBROUTINE s_prim_tot

    INTEGER :: i
    INTEGER :: iv
    INTEGER :: impi
    INTEGER :: arg
    INTEGER :: pos

    pos = 0
    DO impi = 0,mpi_size-1

       DO iv = 1,Nveul
          DO i = 1,i_Nx_nopad(impi)
             arg = i_Nx_nopad(impi)*( iv-1 ) + i + idisp(impi)
             prim_tot(iv)%f(i+pos) = precv(arg)
          END DO
       END DO
       pos = pos + i_Nx_nopad(impi)

    END DO

  END SUBROUTINE s_prim_tot

  !========================================================================

  SUBROUTINE s_rad_tot

    INTEGER :: i
    INTEGER :: iout
    INTEGER :: impi
    INTEGER :: arg
    INTEGER :: pos

    pos = 0
    DO impi = 0,mpi_size-1

       DO iout = 1,NR0out
          DO i = 1,i_Nx_nopad(impi)
             arg = i_Nx_nopad(impi)*( iout-1 ) + i + idisp1(impi)
             rad_tot(iout)%f(i+pos) = precv1(arg)
          END DO
       END DO
       pos = pos + i_Nx_nopad(impi)

    END DO

  END SUBROUTINE s_rad_tot

  !========================================================================

  SUBROUTINE s_output

    INTEGER :: i
    INTEGER :: iv
    INTEGER :: ir
    INTEGER :: iout
    INTEGER :: arg
    INTEGER :: NRh
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

    IF ( it==0.AND.mpi_rank==0 ) THEN

       ! numerical specification
       NRh = ( NR0+1 )/2
       CALL DATE_AND_TIME( today_date,today_time )
       OPEN(UNIT=1,FILE='./setup.dat')
       WRITE(1,*) 'date: ', today_date, ', time: ', today_time
       WRITE(1,*) 'nonlinear computation: ', nonlin
       IF ( nonlin=='y' ) THEN
          WRITE(1,*) 'smooth IC: ', smoothic
       END IF
       WRITE(1,*) 'IC type: ', ictype
       WRITE(1,*) 'time marching: ', timescheme
       WRITE(1,*) 'CFL =', CFL
       WRITE(1,*) 'final time =', finaltime
       WRITE(1,*) 'timestep splitting: ', timesplit
       WRITE(1,*) 'WENO order: ', wenoord
       WRITE(1,*) 'WENO tolerance =', eps
       WRITE(1,*) 'accuracy for adaptive RK =', accuracy
       WRITE(1,*) 'monotonicity preserving: ', mpweno
       WRITE(1,*) 'characteristic decomposition: ', chardecomp
       WRITE(1,*) 'correct negative values: ', negative
       WRITE(1,*) 'source terms: ', source
       WRITE(1,*) 'viscous computation: ', viscous
       WRITE(1,*) 'polytropic gas: ', polytropic
       WRITE(1,*) 'thermal behavior: ', thermal
       WRITE(1,*) 'transfer model: ', model
       WRITE(1,*) 'gas/vapor bubbles: ', vapor
       WRITE(1,*) 'liquid: ', liquid
       WRITE(1,*) 'noncondensible gas: ', gas
       WRITE(1,*) 'quadrature (Simpson or Gauss): ', nquad
       WRITE(1,*) 'Nx_tot =', Nx_tot,', Nt =', Nt
       WRITE(1,*) 'flow file output at intervals of ', interval
       WRITE(1,*) 'length in x: ',length%x
       WRITE(1,*) 'grid stretching: ', stretch
       WRITE(1,*) 'dt =', dt
       WRITE(1,*) 'left BC: ', xbound%beg, ', right BC: ', xbound%end
       WRITE(1,*) 'modify cl at boundaries: ', modifycl
       WRITE(1,*) 'free plate (left BC): ', xincoming%beg
       WRITE(1,*) 'incoming wave (right BC): ', xincoming%end
       WRITE(1,*) 'n_tait =', n_tait, ', B_tait =', B_tait
       WRITE(1,*) 'cl0 = ', cl0, ', tau =', xdecay/cl0
       WRITE(1,*) 'xs =', xs, ', ps =', ps
       WRITE(1,*) 'fsi =', fsi, ', 1/mp =', mp_inv
       WRITE(1,*) 'gamma_b =', gamma_b
       WRITE(1,*) 'Ca =', Ca
       WRITE(1,*) 'We =', We
       WRITE(1,*) 'pv =', pv
       IF ( viscous=='y' ) THEN
          WRITE(1,*) 'Re =', 1.D0/Re_inv
       END IF
       IF ( polytropic=='n' ) THEN
          WRITE(1,*) 'Pe_T^ref =', Pe_T(NRh)
          WRITE(1,*) 'Pe_c =', Pe_c
          WRITE(1,*) 'Pe_T^ref based on omegaN =', Pe_T(NRh)*omegaN(NRh)
          WRITE(1,*) 'Pe_c^ref based on omegaN =', Pe_c*omegaN(NRh)
       END IF
       WRITE(1,*) 'Us =', Us, ', Ms =', Ms
       WRITE(1,*) 'R0ref [m] =', R0ref
       WRITE(1,*) 'standard deviation of f(R0):', sd
       WRITE(1,*) 'alpha0 =', vf0
       WRITE(1,*) 'NR0 =', NR0
       WRITE(1,*) 'iR1 =', ( iR1(iout),iout=1,NR0out )
       WRITE(1,*) 'mpi_size =', mpi_size
       CLOSE(1)
       ! weights
       IF ( disperse=='poly' ) THEN
          OPEN(UNIT=2,FILE='./weights.dat',STATUS='UNKNOWN',FORM='FORMATTED')
          DO ir = 1,NR0
             WRITE(2,*) R0(ir), weight(ir)
          END DO
          CLOSE(2)
       END IF
       ! grid
       OPEN(UNIT=3,FILE='./xgrid.dat',STATUS='UNKNOWN',FORM='FORMATTED')
       DO i = 1,Nx_tot
          WRITE(3,*) xgrid_tot(i), dxs_tot(i)
       END DO
       CLOSE(3)

    END IF

    ! compute primitive variables at each processor
    CALL s_cons2prim
    ! gather (rho,u,alpha) at mpi_rank==0
    ! need prim(1:Nveul)%f to compute pl
    DO iv = 1,Nveul
       DO i = 1,Nx_nopad
          arg = Nx_nopad*( iv-1 ) + i
          psend(arg) = prim(iv)%f(i+i_nopad_mn-1)
       END DO
    END DO
    CALL MPI_GATHERV( psend(1),Nveul*Nx_nopad,MPI_DOUBLE_PRECISION, &
                      precv(1),ircnt,idisp,MPI_DOUBLE_PRECISION, &
                      0,MPI_COMM_WORLD,mpi_err )
    ! gather bubble radius at mpi_rank==0
    DO iout = 1,NR0out
       DO i = 1,Nx_nopad
          arg = Nx_nopad*( iout-1 ) + i
          psend1(arg) = prim(ibub(1,iR1(iout)))%f(i+i_nopad_mn-1)
       END DO
    END DO
    CALL MPI_GATHERV( psend1(1),NR0out*Nx_nopad,MPI_DOUBLE_PRECISION, &
                      precv1(1),ircnt1,idisp1,MPI_DOUBLE_PRECISION, &
                      0,MPI_COMM_WORLD,mpi_err )

    ! update max. CFL
    IF ( mpi_rank==0 ) THEN
       CALL s_prim_tot
       CALL s_rad_tot
       CFLmx1 = 0.D0
       DO i = 1,Nx_tot
          vftmp = prim_tot(Nveul)%f(i)
          vftmp1 = 1.D0 - vftmp
          pl = ( pl0+B_tait )*( prim_tot(1)%f(i)/vftmp1 )**n_tait
          cltmp = DSQRT( n_tait*pl/vftmp1/prim_tot(1)%f(i) )
          Stmp = DABS( prim_tot(2)%f(i) ) + cltmp
          CFLtmp = dt*Stmp/dxs_tot(i)
          IF ( CFLmx1<CFLtmp ) CFLmx1 = CFLtmp
       END DO
       IF ( CFLmx<CFLmx1 ) THEN
          CFLmx = CFLmx1
          OPEN(UNIT=4,FILE='./maxcfl.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED',POSITION='APPEND')
          WRITE(4,*) time, CFLmx
          CLOSE(4)
       END IF
    END IF

    ! wall velocity & pressure
    IF ( xincoming%beg=='freeplate'.AND.mpi_rank==0 ) THEN
       OPEN(UNIT=5,FILE='wall.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(5,*) time, uwall, pwall
       CLOSE(5)
    ELSE IF ( xincoming%end=='freeplate'.AND.mpi_rank==0 ) THEN
       OPEN(UNIT=5,FILE='wall.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(5,*) time, uwall, pwall
       CLOSE(5)
    END IF

    ! probe measurement
    IF ( mpi_rank==mpi_rank_probe1 ) THEN

       vftmp = prim(Nveul)%f(i_probe1)
       vftmp1 = 1.D0 - vftmp
       pl = ( pl0+B_tait )*( prim(1)%f(i_probe1)/vftmp1 )**n_tait - B_tait
       OPEN(UNIT=7,FILE='./probe1.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(7,*) time, pl, vftmp
       CLOSE(7)

    END IF
    IF ( mpi_rank==mpi_rank_probe2 ) THEN

       vftmp = prim(Nveul)%f(i_probe2)
       vftmp1 = 1.D0 - vftmp
       pl = ( pl0+B_tait )*( prim(1)%f(i_probe2)/vftmp1 )**n_tait - B_tait
       OPEN(UNIT=8,FILE='./probe2.dat',STATUS='UNKNOWN',FORM='FORMATTED', &
            POSITION='APPEND')
       WRITE(8,*) time, pl, vftmp
       CLOSE(8)

    END IF

    ! flow files
    IF ( it==0.OR.MOD(it,interval)==0 ) THEN

       ! output time
       IF ( mpi_rank==0 ) THEN
          OPEN(UNIT=9,FILE='./flowtime.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED',POSITION='APPEND')
          WRITE(9,*) itout, time
          CLOSE(9)
       END IF
       ! label for flow files
!       WRITE(notmp,*) itout
       WRITE(notmp,'(I7)') itout
       IF ( itout<10 ) THEN
          no = '00'//TRIM( ADJUSTL(notmp) )
       ELSE IF ( itout>=10.AND.itout<100 ) THEN
          no = '0'//TRIM( ADJUSTL(notmp) )
       ELSE
          no = TRIM( ADJUSTL(notmp) )
       END IF
       
       IF ( mpi_rank==0 ) THEN
          ! averaged flow field
          OPEN(UNIT=10,FILE='D/p.'//no//'.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED')
          DO i = 1,Nx_tot
             vftmp = prim_tot(Nveul)%f(i)
             vftmp1 = 1.D0 - vftmp
             pl = ( pl0+B_tait )*( prim_tot(1)%f(i)/vftmp1 )**n_tait &
                - B_tait
             WRITE(10,*) xgrid_tot(i),pl/pl0
          END DO
          CLOSE(10)
 
          OPEN(UNIT=10,FILE='D/alf.'//no//'.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED')
          DO i = 1,Nx_tot
             WRITE(10,*) xgrid_tot(i),prim_tot(Nveul)%f(i)
          END DO
          CLOSE(10)

          OPEN(UNIT=10,FILE='D/u.'//no//'.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED')
          DO i = 1,Nx_tot
             WRITE(10,*) xgrid_tot(i),prim_tot(2)%f(i)
          END DO
          CLOSE(10)

          ! bubble radius field
          OPEN(UNIT=11,FILE='D/radi.'//no//'.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED')
          DO i = 1,Nx_tot
             WRITE(11,*) xgrid_tot(i),( rad_tot(iout)%f(i),iout=1,NR0out )
          END DO
          CLOSE(11)
       END IF


       ! integrands R(x,t;R0)
!       IF ( itout==30 ) THEN
!          IF ( mpi_rank==mpi_rank_integ1 ) THEN
!             OPEN(UNIT=12,FILE='./integ1.dat',STATUS='UNKNOWN', &
!                  FORM='FORMATTED')
!             DO ir = 1,NR0
!                WRITE(12,*) prim(ibub(1,ir))%f(i_integ1)
!             END DO
!             CLOSE(12)
!          ELSE IF ( mpi_rank==mpi_rank_integ2 ) THEN
!             OPEN(UNIT=13,FILE='./integ2.dat',STATUS='UNKNOWN', &
!                  FORM='FORMATTED')
!             DO ir = 1,NR0
!                WRITE(13,*) prim(ibub(1,ir))%f(i_integ2)
!             END DO
!             CLOSE(13)
!          END IF
!       END IF
       IF ( mpi_rank==mpi_rank_integ1 ) THEN
          IF ( nonlin=='n'.AND.ictype=='4' ) THEN
             OPEN(UNIT=12,FILE='./integ_'//no//'.dat',STATUS='UNKNOWN', &
                  FORM='FORMATTED')
             DO ir = 1,NR0
                WRITE(12,*) prim(ibub(1,ir))%f(i_integ1) - R0(ir)
             END DO
             CLOSE(12)
          ELSE
             OPEN(UNIT=12,FILE='./integ1_'//no//'.dat',STATUS='UNKNOWN', &
                  FORM='FORMATTED')
             DO ir = 1,NR0
                WRITE(12,*) prim(ibub(1,ir))%f(i_integ1)
             END DO
             CLOSE(12)
          END IF
       ELSE IF ( mpi_rank==mpi_rank_integ2 ) THEN
          OPEN(UNIT=13,FILE='./integ2_'//no//'.dat',STATUS='UNKNOWN', &
               FORM='FORMATTED')
          DO ir = 1,NR0
             WRITE(13,*) prim(ibub(1,ir))%f(i_integ2)
          END DO
          CLOSE(13)
       END IF
       itout = itout + 1

    END IF

  END SUBROUTINE s_output

  !========================================================================

  SUBROUTINE s_start_clock

    CALL MPI_BARRIER( MPI_COMM_WORLD,mpi_err )
    elp1 = MPI_WTIME()

  END SUBROUTINE s_start_clock
  
  !========================================================================

  SUBROUTINE s_stop_clock

    CALL MPI_BARRIER( MPI_COMM_WORLD,mpi_err )
    elp2 = MPI_WTIME()
    IF ( mpi_rank==0 ) THEN
       PRINT *, 'computation time [sec] =', elp2 - elp1
       PRINT *, 'computation time [min] =', ( elp2-elp1 )/60.D0
       PRINT *, 'computation time [hr]  =', ( elp2-elp1 )/3600.D0
    END IF

  END  SUBROUTINE s_stop_clock

  !========================================================================

  SUBROUTINE s_stop_inoutvar

    DEALLOCATE( rhs )
    DEALLOCATE( prim )
    DEALLOCATE( cons )
    DEALLOCATE( ibub )
    DEALLOCATE( ieig )
    DEALLOCATE( prim_tot )
    DEALLOCATE( rad_tot )
    IF ( timesplit/='unsplit' ) THEN
       DEALLOCATE( dttry )
       DEALLOCATE( dpldt )
    END IF
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
    DEALLOCATE( xgrid_tot )
    DEALLOCATE( xhalf_tot )
    DEALLOCATE( dxs )
    DEALLOCATE( dxs_tot )
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
    DEALLOCATE( psend )
    DEALLOCATE( precv )
    DEALLOCATE( psend1 )
    DEALLOCATE( precv1 )
    DEALLOCATE( ircnt )
    DEALLOCATE( idisp )
    DEALLOCATE( ircnt1 )
    DEALLOCATE( idisp1 )
    DEALLOCATE( i_Nx_nopad )
    DEALLOCATE( iR1 )

  END SUBROUTINE s_stop_inoutvar

  !========================================================================

END MODULE m_inoutvar
