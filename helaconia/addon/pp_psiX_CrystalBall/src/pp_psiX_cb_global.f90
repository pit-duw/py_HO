MODULE pp_psiX_cb_global
  IMPLICIT NONE
  INCLUDE "pp_psiX_cb.inc"
  INTEGER::cb_istate ! 1:jpsi,2:psi(2S)
  REAL(KIND(1d0))::mjpsi=3.096d0,mpsi2s=3.686d0
  REAL(KIND(1d0))::mY1S=9.46030d0,mY2S=10.02326d0,mY3S=10.3552d0
  REAL(KIND(1d0))::mchic0=3.41475d0,mchic1=3.51066d0,mchic2=3.55620d0
  ! PDG 2014
  REAL(KIND(1d0))::mchib01P=9.85944d0,mchib11P=9.89278d0,mchib21P=9.91221d0
  REAL(KIND(1d0))::mchib02P=10.2325d0,mchib12P=10.25546d0,mchib22P=10.26865d0
  ! update it if experiments can distuiguish chib(3P) states
  REAL(KIND(1d0))::mchib03P=10.534d0,mchib13P=10.534d0,mchib23P=10.534d0
  LOGICAL::cb_includeqq=.TRUE.
END MODULE pp_psiX_cb_global
