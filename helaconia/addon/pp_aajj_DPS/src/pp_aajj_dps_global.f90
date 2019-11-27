MODULE pp_aajj_dps_global
  IMPLICIT NONE
  INCLUDE "pp_aajj_dps.inc"
  REAL(KIND(1d0))::s12_aj,t13_aj,s12_jj,t13_jj,s12_aa,t13_aa
  REAL(KIND(1d0))::EL,Gstrong
  INTEGER,DIMENSION(2,4)::ifldps
  INTEGER,PARAMETER::maxproc_dps=217
  REAL(KIND(1d0)),DIMENSION(2,maxproc_dps)::wme_proc
  LOGICAL,DIMENSION(2,maxproc_dps)::wme_calc
  INTEGER,DIMENSION(2,4,2)::icolun_dps
  REAL(KIND(1d0)),DIMENSION(2,2)::xpi
  REAL(KIND(1d0)),DIMENSION(2,2,-5:5)::strf_pdf_save ! pdf*x
  LOGICAL,DIMENSION(2)::strf_pdf_saveq
  REAL(KIND(1d0)),DIMENSION(2,-5:5,-5:5)::strf_dpdf_save ! dpdf only
  LOGICAL,DIMENSION(2,-5:5,-5:5)::strf_dpdf_saveq
  INTEGER,DIMENSION(maxproc_dps,maxproc_dps,8)::iflh_save
  INTEGER,DIMENSION(8,2)::imothup_save
  LOGICAL,DIMENSION(maxproc_dps,maxproc_dps)::nonzero_subprocess
  REAL(KIND(1d0)),DIMENSION(maxproc_dps,maxproc_dps)::wme_subprocess
  REAL(KIND(1d0))::wmetot,gstrong1,gstrong2
  INTEGER::dpsorsps ! 1 DPS, 2 SPS
  LOGICAL,DIMENSION(maxproc_dps)::contained_subprocess
  REAL(KIND(1d0))::ymax,ymin,ptmin
  ! special cut
  REAL(KIND(1d0))::aajj_pta1u,aajj_pta1l,aajj_etaa1,aajj_pta2u,aajj_pta2l,aajj_etaa2
  REAL(KIND(1d0))::aajj_ptj1u,aajj_ptj1l,aajj_etaj1,aajj_ptj2u,aajj_ptj2l,aajj_etaj2
  ! K factors
  REAL(KIND(1d0))::Kfactor_aa=1d0,Kfactor_jj=1d0,Kfactor_aj=1d0
END MODULE pp_aajj_dps_global
