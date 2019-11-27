  REAL(KIND(1d0)),DIMENSION(4,1:4)::cb_pmom
  REAL(KIND(1d0)),DIMENSION(4,1:5)::cb_hadron_pmom
  REAL(KIND(1d0))::cb_kapa=0.6d0,cb_lam=0.327d0,cb_ptavg=4.5d0
  INTEGER::cb_n=2
  REAL(KIND(1d0))::mpsi=3.096d0
  REAL(KIND(1d0))::cb_lambdath=0d0
  !COMMON/pp_QQ_cb/cb_pmom,cb_hadron_pmom,cb_kapa,cb_lam,mpsi,cb_lambdath