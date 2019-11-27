MODULE pp_psipsi_DPS_ME
  USE pp_psipsi_dps_global
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  ! Eq.(4) in arXiv:1105.4186
  SUBROUTINE crystalball_gg_psiX(ipsi,wme)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(OUT)::wme
    INTEGER,INTENT(IN)::ipsi
    REAL(KIND(1d0))::KK,pt
    INTEGER::init=0
    INTEGER::nn=2
    REAL(KIND(1d0))::ptavg=4.5d0
    REAL(KIND(1d0))::shat
    !INCLUDE "pp_psipsi_dps.inc"
    SAVE init
    IF(init.EQ.0)THEN
       !KK=dps_lam**2*dps_kapa/mpsi**2
       icol_un(1,1)=1
       icol_un(6,2)=1
       icol_un(6,1)=2
       icol_un(2,2)=2
       icol_un(2,1)=3
       icol_un(1,2)=3
       icol_un(5,1)=0
       icol_un(5,2)=0
       icol_un(3,1)=4
       icol_un(8,2)=4
       icol_un(8,1)=5
       icol_un(4,2)=5
       icol_un(4,1)=6
       icol_un(3,2)=6
       icol_un(7,1)=0
       icol_un(7,2)=0
       AQEDUP=1d0/137d0
       AQCDUP=0.118d0
       init=1
    ENDIF
    KK=dps_lam(ipsi)**2*dps_kapa(ipsi)/mpsi(ipsi)**2
    nn=dps_n(ipsi)
    ptavg=dps_ptavg(ipsi)
    pt=DSQRT(dps_pmom(ipsi,3,1)**2+dps_pmom(ipsi,3,2)**2)
    shat=(dps_pmom(ipsi,1,1)+dps_pmom(ipsi,2,1))**2&
         +(dps_pmom(ipsi,1,2)+dps_pmom(ipsi,2,2))**2&
         +(dps_pmom(ipsi,1,3)+dps_pmom(ipsi,2,3))**2
    shat=(dps_pmom(ipsi,1,4)+dps_pmom(ipsi,2,4))**2-shat
    IF(pt.LE.ptavg)THEN
       wme=KK*shat*DEXP(-dps_kapa(ipsi)*pt**2/mpsi(ipsi)**2)
    ELSE
       wme=KK*shat*DEXP(-dps_kapa(ipsi)*ptavg**2/mpsi(ipsi)**2)&
            *1d0/(1d0+dps_kapa(ipsi)/nn*(pt**2-ptavg**2)/mpsi(ipsi)**2)**nn
    ENDIF
    RETURN
  END SUBROUTINE crystalball_gg_psiX
END MODULE pp_psipsi_DPS_ME
