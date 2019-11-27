MODULE pp_QQ_cb_global
  IMPLICIT NONE
  INCLUDE "pp_QQ_cb.inc"
  INTEGER::cb_istate ! 1:ccbar,2:D0+D0bar,3:D*+-,4:D+-
  REAL(KIND(1d0))::mc=1.5d0,mD0=1.8d0,mDstar=2.0d0,mDpm=1.8d0
  LOGICAL::cb_includeqq=.TRUE.
END MODULE pp_QQ_cb_global
