!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Module computing bubble dynamics
!
!  Last update: May 19, 2009
!  Author: Keita Ando
!  Department of Mechanical Engineering
!  Division of Engineering and Applied Science
!  California Institute of Technology, Pasadena CA 91125
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

MODULE m_bubbles

  USE m_globalvar
  IMPLICIT NONE

  ! bubble wall properties
  REAL(KIND(0.D0)) :: chi_vw
  REAL(KIND(0.D0)) :: k_bw
  REAL(KIND(0.D0)) :: rho_bw

  CONTAINS

  !========================================================================

  SUBROUTINE s_bwproperty( pb )

    REAL(KIND(0.D0)), INTENT(IN) :: pb
    REAL(KIND(0.D0)) :: x_vw

    IF ( vapor=='n' ) THEN
       ! mass fraction of vapor
       chi_vw = 0.D0
       ! mole fraction of vapor & thermal conductivity of gas mixture
       x_vw = 0.D0
       k_bw = k_n(iR0)
       ! gas mixture density
       rho_bw = pb/( R_n*Tw )
    ELSE
       ! mass fraction of vapor
       chi_vw = 1.D0/( 1.D0+R_v/R_n*(pb/pv-1.D0) )
       ! mole fraction of vapor & thermal conductivity of gas mixture
       x_vw = M_n*chi_vw/( M_v+(M_n-M_v)*chi_vw )
       k_bw = x_vw*k_v(iR0)/( x_vw+(1.D0-x_vw)*phi_vn ) &
            + ( 1.D0-x_vw )*k_n(iR0)/( x_vw*phi_nv+1.D0-x_vw )
       ! gas mixture density
       rho_bw = pv/( chi_vw*R_v*Tw )
    END IF

  END SUBROUTINE s_bwproperty

  !========================================================================

  FUNCTION f_vflux( R,V,mass_v )

    REAL(KIND(0.D0)), INTENT(IN) :: R
    REAL(KIND(0.D0)), INTENT(IN) :: V
    REAL(KIND(0.D0)), INTENT(IN) :: mass_v
    REAL(KIND(0.D0)) :: chi_bar
    REAL(KIND(0.D0)) :: grad_chi
    REAL(KIND(0.D0)) :: f_vflux

    IF ( vapor=='y' ) THEN

       IF ( thermal=='transfer' ) THEN
          ! constant transfer model
          chi_bar = mass_v/( mass_v+mass_n0(iR0) )
          grad_chi = -Re_trans_c(iR0)*( chi_bar-chi_vw )
          f_vflux = rho_bw*grad_chi/Pe_c/( 1.D0-chi_vw )/R
       ELSE
          ! irrerevant procedure (polytropic)
          f_vflux = pv*V/( R_v*Tw )
       END IF

    ELSE IF ( vapor=='n' ) THEN

       f_vflux = 0.D0

    END IF

  END FUNCTION f_vflux

  !========================================================================

  FUNCTION f_bpres_dot( vflux,R,V,pb,mass_v )

    REAL(KIND(0.D0)), INTENT(IN) :: vflux
    REAL(KIND(0.D0)), INTENT(IN) :: R
    REAL(KIND(0.D0)), INTENT(IN) :: V
    REAL(KIND(0.D0)), INTENT(IN) :: pb
    REAL(KIND(0.D0)), INTENT(IN) :: mass_v
    REAL(KIND(0.D0)) :: T_bar
    REAL(KIND(0.D0)) :: grad_T
    REAL(KIND(0.D0)) :: tmp1, tmp2
    REAL(KIND(0.D0)) :: f_bpres_dot

    IF ( thermal=='transfer'.AND.model=='Preston' ) THEN
       ! constant transfer model
       T_bar = Tw*( pb/pb0(iR0) )*( R/R0(iR0) )**3 &
             * ( mass_n0(iR0)+mass_v0(iR0) )/( mass_n0(iR0)+mass_v )
       grad_T = -Re_trans_T(iR0)*( T_bar-Tw )
       f_bpres_dot = 3.D0*gamma_b*( -V*pb+vflux*R_v*Tw &
                   + pb0(iR0)*k_bw*grad_T/Pe_T(iR0)/R )/R
    ELSE IF ( thermal=='transfer'.AND.model=='Sugiyama' ) THEN
       ! only for gas bubble
       T_bar = Tw*( pb/pb0(iR0) )*( R/R0(iR0) )**3
       tmp1 = -V*pb - pb0(iR0)*k_bw/Pe_T(iR0)/R &
            * ( Re_trans_T(iR0)*(T_bar-Tw) &
            + Im_trans_T(iR0)*T_bar/omegaN(iR0)*3.D0*V/R )
       tmp2 = R/3.D0/gamma_b &
            + pb0(iR0)*k_bw*Im_trans_T(iR0)*T_bar &
            / Pe_T(iR0)/R/omegaN(iR0)/pb
       f_bpres_dot = tmp1/tmp2
    ELSE
       f_bpres_dot = -3.D0*gamma_b*V/R*( pb-pv )
    END IF

  END FUNCTION f_bpres_dot

  !========================================================================

  FUNCTION f_bwpres_1( R,V,pb )

    REAL(KIND(0.D0)), INTENT(IN) :: R
    REAL(KIND(0.D0)), INTENT(IN) :: V
    REAL(KIND(0.D0)), INTENT(IN) :: pb
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: f_bwpres_1

    IF ( polytropic=='y' ) THEN
       tmp = ( R0(iR0)/R )**( 3.D0*gamma_b )
       tmp = ( Ca+2.D0/We/R0(iR0) )*tmp - Ca + 1.D0
    ELSE IF ( polytropic=='n' ) THEN
       tmp = pb
    END IF
    f_bwpres_1 = tmp - 4.D0*Re_inv*V/R - 2.D0/( We*R )

  END FUNCTION f_bwpres_1

  !========================================================================

  FUNCTION f_bwpres( R,V,pb )

    REAL(KIND(0.D0)), INTENT(IN) :: R
    REAL(KIND(0.D0)), INTENT(IN) :: V
    REAL(KIND(0.D0)), INTENT(IN) :: pb
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: f_bwpres

    IF ( polytropic=='y' ) THEN
       tmp = ( R0(iR0)/R )**( 3.D0*gamma_b )
       tmp = ( Ca+2.D0/We/R0(iR0) )*tmp - Ca
    ELSE IF ( polytropic=='n' ) THEN
       tmp = pb - 1.D0
    END IF
    f_bwpres = tmp - 4.D0*Re_inv*V/R - 2.D0/( R*We )

  END FUNCTION f_bwpres

  !========================================================================

  FUNCTION f_enthal( Cpinf,Cpbw )

    REAL(KIND(0.D0)), INTENT(IN) :: Cpinf
    REAL(KIND(0.D0)), INTENT(IN) :: Cpbw
    REAL(KIND(0.D0)) :: tmp1, tmp2, tmp3
    REAL(KIND(0.D0)) :: f_enthal

    tmp1 = ( n_tait-1.D0 )/n_tait
    tmp2 = (  Cpbw/(1.D0+B_tait)+1.D0 )**tmp1
    tmp3 = ( Cpinf/(1.D0+B_tait)+1.D0 )**tmp1
    f_enthal = ( tmp2-tmp3 )*n_tait*( 1.D0+B_tait )/( n_tait-1.D0 )

  END FUNCTION f_enthal

  !========================================================================

  FUNCTION f_enthal_dot( R,V,Cpbw,Cpinf,pb_dot,Cpinf_dot )

    REAL(KIND(0.D0)), INTENT(IN) :: R
    REAL(KIND(0.D0)), INTENT(IN) :: V
    REAL(KIND(0.D0)), INTENT(IN) :: Cpbw
    REAL(KIND(0.D0)), INTENT(IN) :: Cpinf
    REAL(KIND(0.D0)), INTENT(IN) :: pb_dot
    REAL(KIND(0.D0)), INTENT(IN) :: Cpinf_dot
    REAL(KIND(0.D0)) :: tmp1, tmp2
    REAL(KIND(0.D0)) :: f_enthal_dot

    IF ( polytropic=='y' ) THEN
       tmp1 = ( R0(iR0)/R )**( 3.D0*gamma_b )
       tmp1 = -3.D0*gamma_b*( Ca+2.D0/We/R0(iR0) )*tmp1*V/R
    ELSE IF ( polytropic=='n' ) THEN
       tmp1 = pb_dot
    END IF
    tmp2 = ( 2.D0/We+4.D0*Re_inv*V )*V/R**2

    f_enthal_dot = &
             ( Cpbw/(1.D0+B_tait)+1.D0 )**( -1.D0/n_tait )*( tmp1+tmp2 ) &
           - ( Cpinf/(1.D0+B_tait)+1.D0 )**( -1.D0/n_tait )*Cpinf_dot

  END FUNCTION f_enthal_dot

  !========================================================================

  FUNCTION f_sound( Cpinf,H )

    REAL(KIND(0.D0)), INTENT(IN) :: Cpinf
    REAL(KIND(0.D0)), INTENT(IN) :: H
    REAL(KIND(0.D0)) :: tmp
    REAL(KIND(0.D0)) :: f_sound

    tmp = ( Cpinf/(1.D0+B_tait)+1.D0 )**( (n_tait-1.D0)/n_tait )
    tmp = n_tait*( 1.D0+B_tait )*tmp
    f_sound = DSQRT( tmp+(n_tait-1.D0)*H )

  END FUNCTION f_sound

  !========================================================================

  FUNCTION f_gilmore( R,V,pb,pb_dot,pinf,pinf_dot )

    REAL(KIND(0.D0)), INTENT(IN) :: R
    REAL(KIND(0.D0)), INTENT(IN) :: V
    REAL(KIND(0.D0)), INTENT(IN) :: pb
    REAL(KIND(0.D0)), INTENT(IN) :: pb_dot
    REAL(KIND(0.D0)), INTENT(IN) :: pinf
    REAL(KIND(0.D0)), INTENT(IN) :: pinf_dot
    REAL(KIND(0.D0)) :: Cpinf
    REAL(KIND(0.D0)) :: Cpbw
    REAL(KIND(0.D0)) :: H
    REAL(KIND(0.D0)) :: H_dot
    REAL(KIND(0.D0)) :: A
    REAL(KIND(0.D0)) :: tmp1, tmp2, tmp3
    REAL(KIND(0.D0)) :: f_gilmore

    Cpinf = pinf - pl0
    Cpbw = f_bwpres( R,V,pb )
    H = f_enthal( Cpinf,Cpbw )
    H_dot = f_enthal_dot( R,V,Cpbw,Cpinf,pb_dot,pinf_dot )
    A = f_sound( Cpinf,H )
    tmp1 = V/A
    tmp2 = 1.D0 + 4.D0*Re_inv/A/R*( Cpbw/(1.D0+B_tait)+1.D0 ) &
         **( -1.D0/n_tait )
    tmp3 = 1.5D0*V**2*( tmp1/3.D0-1.D0 ) + H*( 1.D0+tmp1 ) &
         + R*H_dot*( 1.D0-tmp1 )/A

    f_gilmore = tmp3/( R*(1.D0-tmp1)*tmp2 )

  END FUNCTION f_gilmore

  !========================================================================

END MODULE m_bubbles
