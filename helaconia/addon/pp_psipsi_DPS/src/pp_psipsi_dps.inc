  REAL(KIND(1d0)),DIMENSION(2,4,1:4)::dps_pmom
  REAL(KIND(1d0)),DIMENSION(8,1:5)::dps_hadron_pmom
  REAL(KIND(1d0)),DIMENSION(2)::dps_kapa=(/0.6d0,0.6d0/)
  REAL(KIND(1d0)),DIMENSION(2)::dps_lam=(/0.327d0,0.327d0/)
  REAL(KIND(1d0)),DIMENSION(2)::dps_ptavg=(/4.5d0,4.5d0/)
  INTEGER,DIMENSION(2)::dps_n=(/2,2/)
  REAL(KIND(1d0))::dps_sigmaeff=14.5d6 ! 1 mb = 1*10^6 nb
  REAL(KIND(1d0)),DIMENSION(2)::mpsi=(/3.096d0,3.096d0/)
  REAL(KIND(1d0)),DIMENSION(2)::dps_lambdath=(/0d0,0d0/)
  !COMMON/pp_psipsi_dps/dps_pmom,dps_hadron_pmom,dps_kapa,dps_lam,dps_sigmaeff,mpsi,dps_lambdath