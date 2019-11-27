MODULE pp_psiX_cb_ME
  USE pp_psiX_cb_global
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  ! Eq.(4) in arXiv:1105.4186
  SUBROUTINE crystalball_gg_psiX(wme)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    REAL(KIND(1d0))::KK,pt
    INTEGER::init=0
    !INTEGER,PARAMETER::nn=2
    !REAL(KIND(1d0)),PARAMETER::ptavg=4.5d0
    REAL(KIND(1d0))::shat
    !INCLUDE "pp_psipsi_dps.inc"
    SAVE init,KK
    IF(init.EQ.0)THEN
       KK=cb_lam**2*cb_kapa/mpsi**2
       icol_un(1,1)=1
       icol_un(4,2)=1
       icol_un(4,1)=2
       icol_un(2,2)=2
       icol_un(2,1)=3
       icol_un(1,2)=3
       icol_un(3,1)=0
       icol_un(3,2)=0
       AQEDUP=1d0/137d0
       AQCDUP=0.118d0
       init=1
    ENDIF
    pt=DSQRT(cb_pmom(3,1)**2+cb_pmom(3,2)**2)
    shat=(cb_pmom(1,1)+cb_pmom(2,1))**2&
         +(cb_pmom(1,2)+cb_pmom(2,2))**2&
         +(cb_pmom(1,3)+cb_pmom(2,3))**2
    shat=(cb_pmom(1,4)+cb_pmom(2,4))**2-shat
    IF(pt.LE.cb_ptavg)THEN
       wme=KK*shat*DEXP(-cb_kapa*pt**2/mpsi**2)
    ELSE
       wme=KK*shat*DEXP(-cb_kapa*cb_ptavg**2/mpsi**2)&
            *1d0/(1d0+cb_kapa/cb_n*(pt**2-cb_ptavg**2)/mpsi**2)**cb_n
    ENDIF
    RETURN
  END SUBROUTINE crystalball_gg_psiX
END MODULE pp_psiX_cb_ME
