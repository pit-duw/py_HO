MODULE pp_psipsi_dps_global
  IMPLICIT NONE
  INCLUDE "pp_psipsi_dps.inc"
  INTEGER,DIMENSION(2)::dps_istate
  REAL(KIND(1d0))::mjpsi=3.096d0,mpsi2s=3.686d0
  REAL(KIND(1d0))::mY1S=9.46030d0,mY2S=10.02326d0,mY3S=10.3552d0
  INTEGER::dps_decay_index=0
  LOGICAL::dps_includeqq=.TRUE.
END MODULE pp_psipsi_dps_global
