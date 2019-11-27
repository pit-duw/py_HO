MODULE func_psi
IMPLICIT NONE
REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
REAL(KIND(1d0)),PARAMETER::EulerGamma=0.57721566490153286061d0

INTERFACE unit_step
   MODULE PROCEDURE unit_step_d
   MODULE PROCEDURE unit_step_s
   MODULE PROCEDURE unit_step_i
END INTERFACE unit_step

INTERFACE DLOG2
   MODULE PROCEDURE DLOG2_d
   MODULE PROCEDURE DLOG2_c
END INTERFACE DLOG2

CONTAINS
FUNCTION lamda_func(a,b,c)
REAL(KIND(1d0))::lamda_func
REAL(KIND(1d0)),INTENT(IN)::a,b,c
lamda_func=a**2+b**2+c**2-2*a*b-2*a*c-2*b*c
END FUNCTION lamda_func

FUNCTION unit_step_d(x)
REAL(KIND(1d0))::unit_step_d
REAL(KIND(1d0))::x
IF(x.LT.0d0) THEN
   unit_step_d=0d0
ELSE
   unit_step_d=1d0
END IF
END FUNCTION unit_step_d

FUNCTION unit_step_s(x)
REAL(KIND(1d0))::unit_step_s
REAL(KIND(1.0))::x
IF(x.LT.0.0) THEN
   unit_step_s=0d0
ELSE
   unit_step_s=1d0
END IF
END FUNCTION unit_step_s

FUNCTION unit_step_i(x)
REAL(KIND(1d0))::unit_step_i
INTEGER::x
IF(x.LT.0) THEN
   unit_step_i=0d0
ELSE
   unit_step_i=1d0
END IF
END FUNCTION unit_step_i

FUNCTION ACosh_p(x)
REAL(KIND(1d0)),INTENT(IN)::x
REAL(KIND(1d0))::ACosh_p
IF(DABS(x).LT.1d0) THEN
  WRITE(*,*)"Error in ACosh_p!"
  ACosh_p=0
ELSE
  ACosh_p=DLOG(x+DSQRT(x**2-1d0))
END IF
END FUNCTION ACosh_p

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! DLOG2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION DLOG2_d(x)
COMPLEX(KIND(1d0))::DLOG2_d
REAL(KIND(1d0)),INTENT(IN)::x
DLOG2_d=CDLOG(DCMPLX(x))
END FUNCTION DLOG2_d
FUNCTION DLOG2_c(x)
COMPLEX(KIND(1d0))::DLOG2_c
COMPLEX(KIND(1d0)),INTENT(IN)::x
DLOG2_c=CDLOG(x)
END FUNCTION DLOG2_c

FUNCTION scalarproduction2(pmom)
  IMPLICIT NONE
  REAL(KIND(1d0))::scalarproduction2
  REAL(KIND(1d0)),DIMENSION(1:4),INTENT(IN)::pmom
  scalarproduction2=pmom(4)**2-pmom(1)**2-pmom(2)**2-pmom(3)**2
  scalarproduction2=DSQRT(ABS(scalarproduction2))
  RETURN
END FUNCTION scalarproduction2

FUNCTION scalarproduction(pmom1,pmom2)
  IMPLICIT NONE
  REAL(KIND(1d0))::scalarproduction
  REAL(KIND(1d0)),DIMENSION(1:4),INTENT(IN)::pmom1,pmom2
  scalarproduction=pmom1(4)*pmom2(4)-pmom1(1)*pmom2(1)&
       -pmom1(2)*pmom2(2)-pmom1(3)*pmom2(3)
  RETURN
END FUNCTION scalarproduction

FUNCTION rapidity2(p)
! rapidity
REAL(KIND(1d0)),DIMENSION(4)::p
REAL(KIND(1d0))::rapidity2,c
IF(p(4).EQ.0d0)THEN
   rapidity2 = 0d0
   RETURN
ENDIF
c=p(3)/ABS(p(4))
IF(ABS(c).EQ.1d0)THEN
   rapidity2 =0d0
ELSE
   rapidity2=0.5d0*LOG((1+c)/(1-c))
ENDIF
END FUNCTION rapidity2

FUNCTION transverse2(p)
REAL(KIND(1d0)),DIMENSION(4)::p
REAL(KIND(1d0))::transverse2
transverse2=SQRT(p(1)**2+p(2)**2)
END FUNCTION transverse2

FUNCTION ph42(p)
IMPLICIT NONE
REAL(KIND(1d0)),DIMENSION(4),INTENT(IN)::p
REAL(KIND(1d0))::px,py,pz
REAL(KIND(1d0))::ph42,S2,S
px=p(1)
py=p(2)
pz=p(3)
IF(px**2+py**2.LE.0d0)THEN
   IF(pz.GE.0d0)ph42=0d0
   IF(pz.LT.0)ph42=pi
   RETURN
ENDIF
S2=px**2/(px**2+py**2)
S=DSQRT(S2)
IF(S.GT.1d0)THEN
   WRITE(*,*)'PH42(X) WARNING S=',S
   ph42=0
   IF(px.LT.0d0)ph42=pi
   RETURN
ENDIF
IF(px.LT.0)S=-DSQRT(S2)
ph42=DACOS(S)
IF(py.LT.0)ph42=2*pi-ph42
RETURN
END FUNCTION ph42
END MODULE func_psi
