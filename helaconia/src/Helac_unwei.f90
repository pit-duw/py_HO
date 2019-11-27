! MODULE unwei
MODULE Helac_unwei
! from main_mc -> Main_Program
USE Helac_histo
USE Helac_Global
IMPLICIT NONE

INTEGER,PARAMETER,PRIVATE::n1=500
REAL(KIND(1d0)),DIMENSION(n1),PRIVATE::y
INTEGER,PRIVATE::init
SAVE y,init

CONTAINS
SUBROUTINE Helac_put_unwei(x)
!      include 'declare.h'
REAL(KIND(1d0)),INTENT(IN)::x
INTEGER::i
IF(GLOBALINIT_unwei.EQ.0)THEN
	init=0
ENDIF
! use unweighted events
IF(init.EQ.0)THEN
   y(1:n1)=0
   init=1
ENDIF
 
DO i=1,n1
   IF(x.GT.y(i))THEN
      y(i+1:n1)=y(i:n1-1) 
      y(i)=x
      RETURN
   ENDIF
ENDDO
!      if(n.gt.1)y(2:n)=y(1:n-1)
!      y(1)=x
RETURN
END SUBROUTINE Helac_put_unwei
       
!       entry get_unwei(t)
SUBROUTINE Helac_get_unwei(t) 
IMPLICIT NONE
REAL(KIND(1d0)),INTENT(OUT)::t
REAL(KIND(1d0))::a,b,w,sum
INTEGER::nb,i
REAL(KIND(1d0))::hxx,hyy,hee,hnn  
CALL Helac_histo3(200)
a=y(n1)-(1d-11)
b=y(1)+(1d-11)
nb=100
w=1.d0
DO i=1,n1
   CALL Helac_histo1(200,nb,a,b,y(i),w)
ENDDO
OPEN(198,FILE=TRIM(tmp_dir)//'tmp16')
CALL Helac_histo2(200,198,0,w)
REWIND(198)
sum=0
DO i=1,nb
   READ(198,*)hxx,hyy,hee,hnn
!       print*,hx,hy,he,hn
   sum=sum+hnn
!       print*,sum,hn
!       if(hn.le.2)then
   IF(sum.GT.0.9d0*n1)THEN
       t=hxx
       CLOSE(198)
       RETURN
   ENDIF
ENDDO
END SUBROUTINE Helac_get_unwei
END MODULE Helac_unwei
