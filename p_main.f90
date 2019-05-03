!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Program solving one-dimensional homogenous bubbly flow
!  (MPI parallel computation)
!
!  Last update: January 27, 2010
!  (c) Keita Ando & Tim Colonius
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!  === Model assumptions ===
!  1. Tait liquid
!  2. Spherical gas/vapor bubbles (Gilmore, no fission or coalescence)
!  3. low void fraction (no direct bubble/bubble interaction)
!  4. No relative motion between the phases
!  5. Bubble size distribution f(R0) included
!  6. Cold liquid (liquid-phase energy equation decoupled, pv const.)
!  7. Perfect gas and vapor
!  8. Spatially uniform pressure inside the bubble
!  9. No rectified diffusion
!
!  === Numerical Method ===
!  Time-marching schemes: explicit Euler, TVD-RK3
!  Time-step splitting: Godunov, Strang (w/ adaptive RK)
!  FV scheme: component/characteristic-wise WENO(1,3,5)
!  Approx. Riemann solver: HLLC
!  BCs: nonreflective (Thompson), reflective (stationary wall)
!       incoming wave (UNDEX), free plate (FSI)
!  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM p_main

  USE mpi_setup
  USE m_globalvar
  USE m_inoutvar
  USE m_timemarch
  USE m_timesplit
  USE m_rhsvar

  IMPLICIT NONE
  INTEGER :: iv

  ! setup MPI
  CALL s_initialize_mpi
  ! start clock
  CALL s_start_clock
  ! computational conditions
  CALL s_parameter
  CALL s_mpi_indices
  ! grid generation, allocation & index matrices
  CALL s_gridgen ! uniform grid
  IF ( stretch=='y' ) CALL s_stretch_grid
  ! set probe locations
  CALL s_probe_setup
  ! allocate array used in time marching
  CALL s_start_timevar
  ! ICs (primitive)
  IF ( xincoming%end=='undex' .or. xincoming%beg=='undex' ) THEN
     CALL s_undexic
  ELSE
     CALL s_ic
  END IF
  ! WENO coefficients
  CALL s_wenocoeff
  ! store ICs (conserved)
  CALL s_prim2cons
  DO iv = 1,Nv
     rk(iv)%l(1) = cons(iv)
  END DO
  it = 0
  itout = 0
  time = 0.D0
  CALL s_output

  ! start computation
  CALL s_start_rhsvar
  DO it = 1,Nt
     IF (mod(it,50)==0 .and. mpi_rank==0) print*, 'i = ', it, ' of ', Nt
     IF ( timesplit=='strang' ) CALL s_integsrc
     IF ( timescheme=='euler'  ) CALL s_euler
     IF ( timescheme=='rk3tvd' ) CALL s_rk3tvd
     IF ( timesplit=='strang' ) CALL s_integsrc
     IF ( timesplit=='godunov' ) CALL s_integsrc
     DO iv = 1,Nv
        cons(iv)%f = rk(iv)%l(1)%f
     END DO
     time = time + dt
     CALL s_output

  END DO
  CALL s_stop_rhsvar

  ! deallocate array
  CALL s_stop_timevar
  CALL s_stop_inoutvar
  ! stop clock
  CALL s_stop_clock
  ! stop mpi
  CALL s_finalize_mpi

END PROGRAM p_main
