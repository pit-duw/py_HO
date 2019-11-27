MODULE fit_cb_global
  IMPLICIT NONE
  REAL(KIND(1d0))::fit_xp1,fit_xp2,fit_ehat
  REAL(KIND(1d0)),DIMENSION(4,4)::fit_cb_pmom
  INTEGER::fit_iPDFSUP1=10042
  REAL(KIND(1d0))::fit_sqrtS,fit_ptup,fit_ptlow,fit_yup,fit_ylow
  REAL(KIND(1d0))::fit_kapa,fit_lam,fit_nn,fit_ptavg
  REAL(KIND(1d0))::mc=1.5d0,mD0=1.5d0,mDstar=1.5d0,mDpm=1.5d0,mpsi
  INTEGER::fit_istate ! 1:c+cbar,2:D0+D0bar,3:D*+-,4:D+-
  LOGICAL::fit_uselhapdf=.FALSE.
  LOGICAL::fit_includeqq=.TRUE.
END MODULE fit_cb_global
