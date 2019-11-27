MODULE Kinetic_Func
USE Helac_Global
USE Helac_Func_1
IMPLICIT NONE

CONTAINS

FUNCTION flambda(x1,x2,x3)
REAL(KIND=DBL),INTENT(IN)::x1,x2,x3
REAL(KIND=DBL)::flambda,f
f=x1**2+x2**2+x3**2-dnou(2)*x1*x2-dnou(2)*x1*x3-dnou(2)*x2*x3
flambda=DSQRT(f)
END FUNCTION flambda

FUNCTION cosij(a,b)
REAL(KIND=DBL),DIMENSION(3),INTENT(IN)::a,b
REAL(KIND=DBL)::cosij,den
cosij=a(1)*b(1)+a(2)*b(2)+a(3)*b(3)
den=DSQRT(a(1)**2+a(2)**2+a(3)**2)*&
     DSQRT(b(1)**2+b(2)**2+b(3)**2)
IF(den.EQ.0d0)THEN
   cosij=1d0
   RETURN
ENDIF
cosij=cosij/DSQRT(a(1)**2+a(2)**2+a(3)**2)
cosij=cosij/DSQRT(b(1)**2+b(2)**2+b(3)**2)
END FUNCTION cosij

FUNCTION scalar_product(p1,p2)
! momentums are in normal representation with fourth comp is the zero comp
REAL(KIND=DBL),DIMENSION(4),INTENT(IN)::p1,p2
REAL(KIND=DBL)::scalar_product
scalar_product=p1(4)*p2(4)-p1(3)*p2(3)-p1(2)*p2(2)-p1(1)*p2(1)
END FUNCTION scalar_product


SUBROUTINE cmstolab()
! from the cms to the lab frame
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(20,4)::p
REAL(KIND=DBL),DIMENSION(4)::pboo
REAL(KIND=DBL)::exp1,exp2
INTEGER::i,icut
REAL(KIND=DBL)::q,e
DO i=1,n
   p(i,1:4)=phegas_pmom(i,1:4)
ENDDO
icut=0
IF(istruc.EQ.1.OR..NOT.labeqcoll)THEN
! xp1(Sqrt(S)/2,0,0,Sqrt(S)/2); xp2(Sqrt(S)/2,0,0,-Sqrt(S)/2)
! Sqrt(s)=SQRT(xp1*xp2)*Sqrt(S)
! ehat=Sqrt(s);e=Sqrt(S)
  q=ehat                   
  e=q/SQRT(xp1*xp2)
  exp1=ebeam(1)*xp1
  exp2=ebeam(2)*xp2
  IF(fixtarget)THEN
     IF(.NOT.fixtargetrev)THEN
        exp1=exp1+FT_M*xp1/2d0
        exp2=exp2/2d0
     ELSE
        exp1=exp1/2d0
        exp2=exp2+FT_M*xp2/2d0
     ENDIF
  ENDIF
  pboo(4)=(exp1+exp2)   ! the zero component
  pboo(3)=(exp1-exp2)   ! the z component
  pboo(1:2)=0
  ! boost p(i,1:4) via pboo to the lab frame
  DO i=1,n
     CALL BOOSTL(q,pboo,p(i,1:4))
  ENDDO
ENDIF
phegas_pmom(1:n,1:4)=p(1:n,1:4)
END SUBROUTINE cmstolab     

SUBROUTINE labtocms()
! from the lab frame to the cms
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(20,4)::p
REAL(KIND=DBL),DIMENSION(4)::pboo
REAL(KIND=DBL)::exp1,exp2
INTEGER::i,icut
REAL(KIND=DBL)::q,e
DO i=1,n
   p(i,1:4)=phegas_pmom(i,1:4)
ENDDO
icut=0
IF(istruc.EQ.1.OR..NOT.labeqcoll)THEN
! xp1(Sqrt(S)/2,0,0,Sqrt(S)/2); xp2(Sqrt(S)/2,0,0,-Sqrt(S)/2)
! Sqrt(s)=SQRT(xp1*xp2)*Sqrt(S)
! ehat=Sqrt(s);e=Sqrt(S)
  q=ehat                   
  e=q/SQRT(xp1*xp2)
  exp1=ebeam(1)*xp1
  exp2=ebeam(2)*xp2
  IF(fixtarget)THEN
     IF(.NOT.fixtargetrev)THEN
        exp1=exp1+FT_M*xp1/2d0
        exp2=exp2/2d0
     ELSE
        exp1=exp1/2d0
        exp2=exp2+FT_M*xp2/2d0
     ENDIF
  ENDIF
  pboo(4)=(exp1+exp2)   ! the zero component
  pboo(3)=-(exp1-exp2)   ! the z component
  pboo(1:2)=0
  ! boost p(i,1:4) via pboo to the cms frame
  DO i=1,n
     CALL BOOSTL(q,pboo,p(i,1:4))
  ENDDO
ENDIF
phegas_pmom(1:n,1:4)=p(1:n,1:4)
END SUBROUTINE labtocms

FUNCTION prapidity(p)
! psudo-rapidity
REAL(KIND=DBL),DIMENSION(4)::p
REAL(KIND=DBL)::prapidity,c
IF(p(1)**2+p(2)**2+p(3)**2.LE.0d0)THEN
   prapidity = 0d0
   RETURN
ENDIF
c=p(3)/SQRT(p(1)**2+p(2)**2+p(3)**2)
IF(ABS(c).EQ.1d0)THEN
   prapidity = 0d0
ELSE
   prapidity=0.5d0*LOG((1+c)/(1-c))
ENDIF
END FUNCTION prapidity

! vector product
FUNCTION vec_crossprod(vec1,vec2)
  ! vec1 crosstime vec2
  IMPLICIT NONE
  REAL(KIND=DBL),DIMENSION(3),INTENT(IN)::vec1,vec2
  REAL(KIND=DBL),DIMENSION(3)::vec_crossprod
  vec_crossprod(1)=vec1(2)*vec2(3)-vec1(3)*vec2(2)
  vec_crossprod(2)=vec1(3)*vec2(1)-vec1(1)*vec2(3)
  vec_crossprod(3)=vec1(1)*vec2(2)-vec1(2)*vec2(1)
  RETURN
END FUNCTION vec_crossprod

FUNCTION rapidity(p)
! rapidity
REAL(KIND=DBL),DIMENSION(4)::p
REAL(KIND=DBL)::rapidity,c
IF(p(4).EQ.0d0)THEN
   rapidity = 0d0
   RETURN
ENDIF
c=p(3)/ABS(p(4))
IF(ABS(c).GE.1d0)THEN
   rapidity =0d0
ELSE
   rapidity=0.5d0*LOG((1+c)/(1-c))
ENDIF
END FUNCTION rapidity

FUNCTION xFeynman(p,Ecms)
! xF is the Feynman parameter (for example see hep-ph/9702263)
! Ecms is the total energy in cms
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(4)::p
REAL(KIND=DBL)::xFeynman,Ecms
IF(Ecms.LT.0d0)THEN
   xFeynman=0d0
ELSE
   xFeynman=2d0*p(3)/Ecms ! yi = Log[xFi*E/mti/2 + Sqrt[1 + (xFi*E/mti/2)^2]]
ENDIF
END FUNCTION xFeynman
   
FUNCTION transverse(p)
REAL(KIND=DBL),DIMENSION(4)::p
REAL(KIND=DBL)::transverse
transverse=SQRT(p(1)**2+p(2)**2)
END FUNCTION transverse

FUNCTION getrapidity(en,pl)
  IMPLICIT NONE
  REAL(KIND(1d0))::getrapidity
  REAL(KIND(1d0)),INTENT(IN)::en,pl
  REAL(KIND(1d0)),PARAMETER::tiny=1.d-8
  REAL(KIND(1d0))::xplus,xminus,y
  xplus=en+pl
  xminus=en-pl
  IF(xplus.GT.tiny.AND.xminus.GT.tiny)THEN
     IF( (xplus/xminus).GT.tiny.AND.(xminus/xplus).GT.tiny)THEN
        y=0.5d0*LOG( xplus/xminus  )
     ELSE
        y=SIGN(1.d0,pl)*1.d8
     ENDIF
  ELSE
     y=SIGN(1.d0,pl)*1.d8
  ENDIF
  getrapidity=y
  RETURN
END FUNCTION getrapidity

SUBROUTINE BOOST(Q,PBOO,PCM,PLB)
! momentums are in normal representation with fourth comp is the zero comp
! Boost PCM via PBOO(PBOO^2=Q^2) to PLB
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::Q
REAL(KIND=DBL),DIMENSION(4),INTENT(IN)::PBOO,PCM
REAL(KIND=DBL),DIMENSION(4),INTENT(OUT)::PLB
REAL(KIND=DBL)::FACT
INTEGER::J
PLB(4)=(PBOO(4)*PCM(4)+PBOO(3)*PCM(3)+PBOO(2)*PCM(2)+PBOO(1)*PCM(1))/Q
FACT=(PLB(4)+PCM(4))/(Q+PBOO(4))
DO J=1,3
  PLB(J)=PCM(J)+FACT*PBOO(J)
ENDDO
END SUBROUTINE BOOST


SUBROUTINE BOOSTL(Q,PBOO,P)
! momentums are in normal representation with fourth comp is the zero comp
! Boost P via PBOO(PBOO^2=Q^2) to PLB,and set to P
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::Q
REAL(KIND=DBL),DIMENSION(4),INTENT(IN)::PBOO
REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::P
REAL(KIND=DBL),DIMENSION(4)::PCM,PLB
REAL(KIND=DBL)::FACT
INTEGER::J
PCM(1:4)=P(1:4)
PLB(4)=(PBOO(4)*PCM(4)+PBOO(3)*PCM(3)+PBOO(2)*PCM(2)+PBOO(1)*PCM(1))/Q
FACT=(PLB(4)+PCM(4))/(Q+PBOO(4))
DO J=1,3
   PLB(J)=PCM(J)+FACT*PBOO(J)
ENDDO
P(1:4)=PLB(1:4)
END SUBROUTINE BOOSTL

SUBROUTINE ROTATEL(PROT,P)
  IMPLICIT NONE
  REAL(KIND=DBL),DIMENSION(3),INTENT(IN)::PROT
  REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::P
  REAL(KIND=DBL),DIMENSION(3)::PTMP
  REAL(KIND=DBL)::costh,sinth,cosphi,sinphi
  REAL(KIND=DBL)::pp,pt
  pp=PROT(1)**2+PROT(2)**2+PROT(3)**2
  IF(pp.LT.0d0)RETURN
  pp=DSQRT(pp)
  costh=PROT(3)/pp
  sinth=DSQRT(1d0-costh**2)
  pt=PROT(1)**2+PROT(2)**2
  pt=DSQRT(pt)
  IF(pt.EQ.0d0)THEN
     cosphi=1d0
     sinphi=0d0
  ELSE
     cosphi=PROT(1)/pt
     sinphi=PROT(2)/pt
  ENDIF
  PTMP(1:3)=P(1:3)
  ! rotate via y aix by th
  PTMP(1)=costh*P(1)+sinth*P(3)
  PTMP(3)=-sinth*P(1)+costh*P(3)
  P(1:3)=PTMP(1:3)
  ! rotate via z aix by phi
  PTMP(1)=cosphi*P(1)-sinphi*P(2)
  PTMP(2)=sinphi*P(1)+cosphi*P(2)
  P(1:3)=PTMP(1:3)
  RETURN
END SUBROUTINE ROTATEL

FUNCTION ph4(px,py,pz)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::px,py,pz
REAL(KIND=DBL)::ph4,S2,S
INTEGER::init=0
REAL(KIND=DBL)::pi
SAVE pi,init
IF(init.EQ.0)THEN
   CALL Helac_mypi(pi)
   init=1
ENDIF
IF(px**2+py**2.LE.0d0)THEN
   IF(pz.GE.0d0)ph4=0d0
   IF(pz.LT.0)ph4=pi
   RETURN
ENDIF
S2=px**2/(px**2+py**2)
S=DSQRT(S2)
IF(S.GT.1d0)THEN
   WRITE(*,*)'PH4(X) WARNING S=',S
   ph4=0
   IF(px.LT.0d0)ph4=pi
   RETURN
ENDIF
IF(px.LT.0)S=-DSQRT(S2)
ph4=DACOS(S)
IF(py.LT.0)ph4=2*pi-ph4
RETURN
END FUNCTION ph4


SUBROUTINE trans_p_to_zq(n1,p,zq)
IMPLICIT NONE
INTEGER,INTENT(IN)::n1
COMPLEX(KIND=DBL),DIMENSION(n1,5),INTENT(OUT):: zq
REAL(KIND=DBL),DIMENSION(n1,5),INTENT(IN)::p
INTEGER::i
REAL(KIND=DBL)::pt,pq
DO i=1,n1
   pt=DSQRT(p(i,1)**2+p(i,2)**2)
   pq=DSQRT(pt**2+p(i,3)**2)
   zq(i,1)= DCMPLX(p(i,4)+p(i,3),p(i,3))
   zq(i,2)= DCMPLX(p(i,4)-p(i,3),pt)
   zq(i,3)= DCMPLX(p(i,1), p(i,2))
   zq(i,4)= DCMPLX(p(i,1),-p(i,2))
   zq(i,5)= DCMPLX( p(i,5) , pq )
ENDDO
END SUBROUTINE trans_p_to_zq


SUBROUTINE transform(n1,zq,p)
IMPLICIT NONE
INTEGER,INTENT(IN)::n1
COMPLEX(KIND=DBL),DIMENSION(n1,5),INTENT(IN)::zq
REAL(KIND=DBL),DIMENSION(n1,5),INTENT(OUT)::p
INTEGER::i
DO i=1,n1 
   p(i,4)=DREAL(zq(i,1)+zq(i,2))/2
   p(i,3)=DREAL(zq(i,1)-zq(i,2))/2
   p(i,1)=DREAL(zq(i,3))
   p(i,2)=DIMAG(zq(i,3))
   p(i,5)=DREAL(zq(i,5))
ENDDO
END SUBROUTINE transform

END MODULE Kinetic_Func
