MODULE fit_cb_global
  IMPLICIT NONE
  REAL(KIND(1d0))::fit_xp1,fit_xp2,fit_ehat
  REAL(KIND(1d0)),DIMENSION(4,4)::fit_cb_pmom
  INTEGER::fit_iPDFSUP1=10042
  REAL(KIND(1d0))::fit_sqrtS,fit_ptup,fit_ptlow,fit_yup,fit_ylow
  REAL(KIND(1d0))::fit_kapa,fit_lam,fit_nn,fit_ptavg
  ! pdg live 2014 values
  REAL(KIND(1d0))::mY1S=9.4603d0,mY2S=10.02326d0,mY3S=10.3552d0,mpsi
  INTEGER::fit_istate ! 1:Y(1S),2:Y(2S),3:Y(3S)
  LOGICAL::fit_uselhapdf=.FALSE.
  !INTEGER::fit_itype ! 1: d^2sigma/dpT/dy, 2: d^2sigma/dpT/dxF
  LOGICAL::fit_includeqq=.TRUE.
END MODULE fit_cb_global
