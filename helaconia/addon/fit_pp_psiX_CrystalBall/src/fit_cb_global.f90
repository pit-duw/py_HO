MODULE fit_cb_global
  IMPLICIT NONE
  REAL(KIND(1d0))::fit_xp1,fit_xp2,fit_ehat
  REAL(KIND(1d0)),DIMENSION(4,4)::fit_cb_pmom
  INTEGER::fit_iPDFSUP1=10042
  REAL(KIND(1d0))::fit_sqrtS,fit_ptup,fit_ptlow,fit_yup,fit_ylow
  REAL(KIND(1d0))::fit_kapa,fit_lam,fit_nn,fit_ptavg
  REAL(KIND(1d0))::mjpsi=3.096d0,mpsi2s=3.686d0,mpsi
  REAL(KIND(1d0))::mchic0=3.41475d0,mchic1=3.51066d0,mchic2=3.55620d0
  INTEGER::fit_istate ! 1:J/psi,2:psi(2S)
  LOGICAL::fit_uselhapdf=.FALSE.
  LOGICAL::fit_includeqq=.TRUE.
END MODULE fit_cb_global
