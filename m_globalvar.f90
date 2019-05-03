!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module storing global variables
!
!  Last update: April 12, 2010
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_globalvar

  IMPLICIT NONE

  ! no. of conserved variables for the Euler equations
  INTEGER :: Nveul
  ! odd no. of descretization points in R0
  ! 1 for the monodisperse case
  INTEGER :: NR0
  ! no. of bubble-dynamic variables
  INTEGER :: Nb
  ! no. of conserved variables for the entire system
  INTEGER :: Nv
  ! max. no. of storage locations for RK
  INTEGER, PARAMETER :: Nrk = 2
  ! WENO order
  INTEGER :: wenoord
  INTEGER :: wenonum
  ! used for left & right eigenvectors
  INTEGER :: Neig1
  INTEGER :: Neig

  ! global computational parameters
  INTEGER :: it
  INTEGER :: itout
  INTEGER :: interval
  INTEGER :: Nt
  REAL(KIND(0.D0)) :: dx
  REAL(KIND(0.D0)) :: dt
  REAL(KIND(0.D0)) :: dtsplit
  REAL(KIND(0.D0)) :: CFL
  REAL(KIND(0.D0)) :: time
  REAL(KIND(0.D0)) :: finaltime
  REAL(KIND(0.D0)) :: accuracy ! for adaptive RK
  CHARACTER(LEN=6) :: timescheme
  CHARACTER(LEN=7) :: timesplit
  CHARACTER(LEN=1) :: mpweno
  CHARACTER(LEN=1) :: chardecomp
  CHARACTER(LEN=1) :: stretch
  CHARACTER(LEN=1) :: negative
  CHARACTER(LEN=1) :: source
  CHARACTER(LEN=1) :: nonlin
  CHARACTER(LEN=1) :: smoothic
  CHARACTER(LEN=2) :: ictype
  CHARACTER(LEN=1) :: viscous
  CHARACTER(LEN=1) :: polytropic
  CHARACTER(LEN=8) :: thermal
  CHARACTER(LEN=1) :: vapor
  CHARACTER(LEN=8) :: model
  CHARACTER(LEN=1) :: nquad
  CHARACTER(LEN=1) :: timehis
  CHARACTER(LEN=18) :: liquid
  CHARACTER(LEN=18) :: gas
  CHARACTER(LEN=4) :: disperse
  CHARACTER(LEN=1) :: modifycl

  ! algebraic constants
  REAL(KIND(0.D0)), PARAMETER :: pi = 3.141592653589793238462643D0
  REAL(KIND(0.D0)), PARAMETER :: third = 1.D0/3.D0
!  REAL(KIND(0.D0)), PARAMETER :: sixth = 1.D0/6.D0
  REAL(KIND(0.D0)), PARAMETER :: twelfth = 1.D0/12.D0
  REAL(KIND(0.D0)), PARAMETER :: four3rd = 4.D0*third
  REAL(KIND(0.D0)), PARAMETER :: pi43 = pi*four3rd

  ! derived types
  TYPE :: dir
     REAL(KIND(0.D0)) :: x
  END TYPE dir
  TYPE :: coordinate
     REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: f
  END TYPE coordinate
  TYPE :: timelevel
     TYPE(coordinate), DIMENSION(Nrk) :: l
  END TYPE timelevel
  TYPE :: wenochar
     TYPE(coordinate), DIMENSION(-2:2) :: w
  END TYPE wenochar
  TYPE :: bc
     CHARACTER(LEN=10) :: beg
     CHARACTER(LEN=10) :: end
  END TYPE bc

  ! domain size & BCs
  TYPE(dir) :: length
  TYPE(bc) :: xbound
  TYPE(bc) :: xincoming

  ! rhs, primitive & conserved variables for both Euler equations
  ! and bubble dynamics equations
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: rhs
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: prim
  TYPE(coordinate), DIMENSION(:), ALLOCATABLE :: cons

  ! time increment used for adaptive RK
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: dttry
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: dpldt

  ! index matrices
  INTEGER, DIMENSION(:,:), ALLOCATABLE, TARGET :: ix
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ibub
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: ieig

  ! grid
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: xgrid
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE, TARGET :: dxs

  ! WENO coefficients
  ! WENO tolerance
  REAL(KIND(0.D0)) :: eps
  ! smoothness indicators
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: betax0
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: betax1
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: betax2
  ! ideal weights
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: dwxr
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: dwxl
  ! polynomials
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: polyrx0
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: polyrx1
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: polyrx2
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: polylx0
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: polylx1
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: polylx2
  ! WENO3
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: w3betax
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: w3drx
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: w3dlx
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: w3polyrx
  REAL(KIND(0.D0)), DIMENSION(:,:), ALLOCATABLE, TARGET :: w3polylx
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE, TARGET :: onesidebegx
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE, TARGET :: onesideendx

  ! Tait EOS
  REAL(KIND(0.D0)) :: n_tait
  REAL(KIND(0.D0)) :: B_tait
  REAL(KIND(0.D0)) :: pl0
  REAL(KIND(0.D0)) :: cl0
  ! vapor pressure
  REAL(KIND(0.D0)) :: pv
  ! bubble wall temperature in the liquid
  REAL(KIND(0.D0)) :: Tw
  ! ratio of specific heats
  REAL(KIND(0.D0)) :: gamma_b
  ! molecular weight
  REAL(KIND(0.D0)) :: M_n
  REAL(KIND(0.D0)) :: M_v
  ! thermal conductivity
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: k_n
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: k_v
  ! gas constants
  REAL(KIND(0.D0)) :: R_n
  REAL(KIND(0.D0)) :: R_v
  ! phi used for computing thermal conductivity
  REAL(KIND(0.D0)) :: phi_vn
  REAL(KIND(0.D0)) :: phi_nv
  ! initial, internal bubble pressure
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: pb0
  ! mass of gas/vapor
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: mass_n0
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: mass_v0
  ! Peclet numbers
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: Pe_T
  REAL(KIND(0.D0)) :: Pe_c
  ! constant transfer coefficients
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: Re_trans_T
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: Re_trans_c
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: Im_trans_T
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: Im_trans_c
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: omegaN
  ! cavitaion, Weber & Reynolds no.
  REAL(KIND(0.D0)) :: Ca
  REAL(KIND(0.D0)) :: We
  REAL(KIND(0.D0)) :: Re_inv

  ! UNDEX & FSI properties
  REAL(KIND(0.D0)) :: xdecay
  REAL(KIND(0.D0)) :: mp_inv
  REAL(KIND(0.D0)) :: uwall
  REAL(KIND(0.D0)) :: pwall

  ! initial radius
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: R0
  ! weights for numerical quadrature over f(R0)
  REAL(KIND(0.D0)), DIMENSION(:), ALLOCATABLE :: weight
  ! label for descretization points in R0
  INTEGER :: iR0

END MODULE m_globalvar
