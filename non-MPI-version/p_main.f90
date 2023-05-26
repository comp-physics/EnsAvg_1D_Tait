!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Program solving one-dimensional homogeneous bubbly flow
!  (Serial code)
!
!  Last update: February 21, 2009
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
!  Time-marching schemes: explicit Euler, TVD-RK3, TVD-RK4
!  Time-step splitting: Godunov, Strang (w/ adaptive RK)
!  FV scheme: component/characteristic-wise WENO5
!  Approx. Riemann solver: HLLC
!  BCs: nonreflective (Thompson), reflective (stationary wall)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM p_main

  USE m_globalvar
  USE m_inoutvar
  USE m_timemarch
  USE m_rhsvar

  IMPLICIT NONE
  INTEGER :: iv

  print*, 'Start up'
  ! computational conditions
  CALL s_parameter
  print*, 'Got parameters'
  ! grid generation, allocation & index matrices
  CALL s_gridgen
  print*, 'Got grid'
  ! set probe locations
  CALL s_probe_setup
  ! allocate array used in time marching
  CALL s_start_timevar
  ! ICs (primitive)
  CALL s_ic
  print*, 'Got IC'
  
  ! WENO coefficient
  CALL s_wenocoeff
  ! store ICs (conserved)
  CALL s_prim2cons
  print*, 'Prim to cons'

  


  DO iv = 1,Nv
     rk(iv)%l(1) = cons(iv)
  END DO
  it = 0
  itout = 0
  time = 0.D0
  CALL s_output
  print*, 'Initial output'

  !stop
  ! start computation
  CALL s_start_rhsvar
  DO it = 1,Nt
     print*, 'Time step: ', it

     IF ( timescheme=='euler'  ) CALL s_euler
     IF ( timescheme=='rk3tvd' ) CALL s_rk3tvd
     IF ( timescheme=='rk4tvd' ) CALL s_rk4tvd
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

END PROGRAM p_main
