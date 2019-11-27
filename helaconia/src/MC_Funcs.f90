MODULE MC_Funcs
USE Helac_Global    
IMPLICIT NONE
SAVE
!INTEGER,PARAMETER,PRIVATE::DBL=SELECTED_REAL_KIND(p=13)
!       include 'mc.h'
INTEGER,PARAMETER,PRIVATE::nx = 140000
!  end of 'mc.h'
REAL(KIND=DBL),DIMENSION(nx),PRIVATE::f0,f1,f2,f3,f4
CONTAINS
SUBROUTINE geti(i,j,w)
!       include 'declare.h'
!       include 'mc.h'
INTEGER,INTENT(IN)::i,j
REAL(KIND=DBL),INTENT(OUT)::w
IF(j.EQ.0)w=f0(i)
IF(j.EQ.1)w=f1(i)
IF(j.EQ.2)w=f2(i)
IF(j.EQ.3)w=f3(i)
IF(j.EQ.4)w=f4(i)
END SUBROUTINE geti

SUBROUTINE clear(i)
! include 'declare.h'
!       include 'mc.h'
INTEGER,INTENT(IN)::i
f0(i)=0
f1(i)=0
f2(i)=0
f3(i)=0
f4(i)=0
END SUBROUTINE clear
       
SUBROUTINE bookin(i,w)
!       include 'declare.h'
!       include 'mc.h'
INTEGER,INTENT(IN)::i
REAL(KIND=DBL),INTENT(IN)::w
INTEGER::init=0
SAVE init
IF(GLOBALINIT_Bookin.EQ.0)THEN
	init=0
ENDIF
IF(init.EQ.0)THEN
	init=1
	f0(1:nx)=0
	f1(1:nx)=0
	f2(1:nx)=0
	f3(1:nx)=0
	f4(1:nx)=0
ENDIF
f0(i)=f0(i)+1
f1(i)=f1(i)+w
f2(i)=f2(i)+w*w
f3(i)=f3(i)+w*w*w
f4(i)=f4(i)+w*w*w*w
END SUBROUTINE bookin
       
SUBROUTINE errest(i,xest,yest,iwr)
! include 'declare.h'
INTEGER,INTENT(IN)::i,iwr
REAL(KIND=DBL)::zest,s0,s1,s2,s3,s4,rn1,rn2,rn3,rn4,v1,v2,v3,v4,v5,v6,v7,v8,v9,v10,v11
REAL(KIND=DBL),INTENT(OUT)::xest,yest
! include 'mc.h'
INTEGER::nout=6
s0=f0(i) ! N
s1=f1(i) ! N*<w/p> ! p is the probability
s2=f2(i) ! N*<w^2/p^2>
s3=f3(i) ! N*<w^3/p^3>
s4=f4(i) ! N*<w^4/p^4>
        
xest=0
yest=0
IF(s0.EQ.0)RETURN
rn1=s0 ! N
rn2=rn1*(s0-1d0) ! N(N-1)
rn3=rn2*(s0-2d0) ! N(N-1)(N-2)
rn4=rn3*(s0-3d0) ! N(N-1)(N-2)(N-3)
v1=s1 ! N<w/p>
v2=s2 ! N<w^2/p^2>
v3=s1**2-v2 ! N^2<w/p>^2-N<w^2/p^2>
v4=s3 ! N<w^3/p^3>
v5=s2*s1-v4 ! N^2<w^2/p^2><w/p>-N<w^3/p^3>
v6=s1**3-v4-3*v5 ! N^3<w/p>^3+2N<w^3/p^3>-3*N^2<w^2/p^2><w/p>
v7=s4 ! N<w^4/p^4>
v8=s3*s1-v7 ! N^2<w^3/p^3><w/p>-N<w^4/p^4>
v9=s2**2-v7 ! N^2<w^2/p^2>**2-N<w^4/p^4>
v10=s2*s1**2-v7-2d0*v8-v9 ! N^3<w^2/p^2>(<w/p>)**2+2N<w^4/p^4>-2N^2<w^3/p^3><w/p>-N^2<w^2/p^2>**2
v11=s1**4-v7-4*v8-3d0*v9-6d0*v10 ! N^4<w/p>**4-6N<w^4/p^4>+8N^2<w^3/p^3><w/p>+3N^2<w^2/p^2>**2-6N^3<w^2/p^2><w/p>**2
xest=v1/rn1 ! <w/p>
IF(iwr.EQ.1)WRITE (nout,3) xest
3   FORMAT('  estimator x: ',d15.6)
IF(rn2.GT.0) THEN
	yest=(v2/rn1-v3/rn2)/rn1 ! (<w^2/p^2>-(<w/p>)^2)/(N-1)
	IF(iwr.EQ.1)WRITE (nout,4) yest
    4     FORMAT('  estimator y: ',d15.6)
    IF(yest.LT.0) THEN
		IF(iwr.EQ.1)WRITE (nout,5) yest
		5     FORMAT('  warning: variance estimator =',d15.6/,'  sign flipped in order to survive')
		yest=DABS(yest)
	ENDIF
ELSE
	IF(iwr.EQ.1)WRITE (nout,6) rn1
    6     FORMAT('  variance estimate not possible, n=',f4.1/,'  estimator put to zero in order to survive')
	yest=0
ENDIF
IF(rn4.GT.0)THEN
	zest=((rn2+ rn3)*v7/ rn1 &
		-4d0*(rn2+rn3)*v8/ rn2 &
		+(rn2-rn3)*v9/ rn2 &
		+4d0*(rn2+2d0*rn3)*v10/rn3 &
		-2d0*(rn2+2d0*rn3)*v11/rn4 &
		)/(rn1**2*rn2**2)
	IF(iwr.EQ.1)WRITE (nout,7) zest
    7     FORMAT('  estimator z: ',d15.6)
	IF(zest.LT.0)THEN
		IF(iwr.EQ.1)WRITE (nout,8) zest
		8       FORMAT('  warning: variance-variance estimator =',d15.6/,'  sign flipped in order to survive')
        zest=DABS(zest)
	ENDIF
ELSE
	IF(iwr.EQ.1)WRITE (nout,9) rn1
    9     FORMAT('  variance-variance estimate not possible, n=',f4.1/, '  estimator put to zero in order to survive')
    zest=0
ENDIF
IF(iwr.EQ.1)THEN
	WRITE (nout,10) xest,DSQRT(yest),yest,DSQRT(zest)
	IF(nunit1.NE.nout)WRITE(nunit1,10)xest,DSQRT(yest),yest,DSQRT(zest)
	10   FORMAT('  average estimate :',d15.6/,16x,'+\- ',d15.6/, '  variance estimate:',d15.6/,16x,'+\- ',d15.6)
ENDIF
IF(DSQRT(zest).GE.0.5*yest.AND.yest.GT.0d0.AND.iwr.EQ.1) WRITE (nout,11)
11   FORMAT('  be aware that the error estimate may be bad!')
!-----------------------------------------------------------
RETURN
END SUBROUTINE errest
END MODULE MC_Funcs
