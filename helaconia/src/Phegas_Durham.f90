MODULE Phegas_Durham_mod
! Refer to Appendix of hep-ph/0512150
USE Helac_Func_1
USE Kinetic_Func
USE Helac_ranmar_mod
IMPLICIT NONE
INTEGER,PARAMETER::maxn=3628800   ! 10! i.e. at most 10 particles
INTEGER::n1,gamma,gamma1,ndim1,nch
INTEGER,DIMENSION(10,maxn)::ia
REAL(KIND(1d0))::cutrap,cutpt,pi
REAL(KIND(1d0))::qq
INTEGER,DIMENSION(maxn)::ich
REAL(KIND(1d0)),DIMENSION(100,4)::pp_durham
SAVE
CONTAINS
SUBROUTINE Phegas_durham(np,et,xt,p)
IMPLICIT NONE
REAL(KIND(1d0)),DIMENSION(100,4),INTENT(OUT)::p
REAL(KIND(1d0)),DIMENSION(4)::pboo,pin
INTEGER::i
REAL(KIND(1d0)),INTENT(IN)::et,xt
INTEGER,DIMENSION(10)::ib
INTEGER::init=0,j,k,flag
INTEGER,INTENT(IN)::np
REAL(KIND(1d0))::f,ff,qx,qy,p1,p2,p3,p4,q0,xres,ch1,sh1,c1,s1,e0,x0
INTEGER,DIMENSION(maxn)::ich0
REAL(KIND(1d0)),DIMENSION(100)::xveg
REAL(KIND(1d0)),DIMENSION(100)::qt,phi,y
SAVE
IF(init.EQ.0)THEN
  n1=np               ! the number of external particles
  cutrap=0.9d0 !2.5d0
  cutpt=60d0
  pi=DACOS(-1d0)     ! constant pi=3.1415926...

  gamma=1
  DO i=1,n1
     gamma=gamma*i  ! gamma=n!
  ENDDO       
  gamma1=gamma/n1/(n1-1) ! gamma1=n!/(n(n-1))

  WRITE(*,*)'gamma,gamma1',gamma,gamma1
  ia(1:10,1:gamma)=0
  DO i=1,gamma
     CALL Helac_permu(ib(1:n1),n1)
     ia(1:n1,i)=ib(1:n1)   ! the permutations of 1~n
  ENDDO
  ndim1=2*n1-2            ! 
  DO i=1,gamma
     ich0(i)=i
  ENDDO

  nch=gamma
  nch=100

  IF(nch.GT.gamma)nch=gamma
  IF(nch.EQ.gamma)THEN
      DO i=1,nch
         ich(i)=i
      ENDDO
  ENDIF
 
  WRITE(*,*)'phegas_durham params:   n,ndim,cutrap,cutpt,gamma,nch'
  WRITE(*,*)'phegas_durham params: ',n1,ndim1,cutrap,cutpt,gamma,nch
  WRITE(*,*)et,xt

  init=1
ENDIF

IF(nch.LT.gamma)THEN
  DO i=1,nch
! 1       continue
     DO
       j=Helac_rnmy(0)*gamma+1 
       ich(i)=ich0(j)
	   flag=0
       DO k=1,i-1
         IF(ich(i).EQ.ich(k))THEN
		    flag=1
		    EXIT
		 ENDIF !goto1
       ENDDO
	   IF(flag.EQ.0)EXIT
	 ENDDO
  ENDDO
ENDIF
! 100  continue
DO
  DO
    qq=DSQRT(et**2-xt**2)   ! the invariance mass of initial states
    DO k=1,ndim1
       xveg(k)=Helac_rnmy(0)
    ENDDO

    DO i=2,n1
       qt(i)=-cutpt*DLOG(xveg(i-1))   ! the pt of final states
    ENDDO
      
    y(1)=0
    DO i=2,n1
       f=DEXP(4*cutrap*(2*xveg(i-1+(n1-1))-1))
       ff=DSQRT(2*f*DCOSH(4*cutrap)-1-f**2)
       ch1=(1+f)*DSINH(2*cutrap)/ff
       sh1=(f-1)*DCOSH(2*cutrap)/ff
       y(i)=DLOG(ch1+sh1)
       y(i)=y(i)+y(i-1)    ! the rapdity of final states
    ENDDO
    phi(1)=0
    DO i=2,n1
       phi(i)=2*pi*Helac_rnmy(0)-phi(i-1)
       phi(i)=phi(i)+phi(i-1)
    ENDDO

! define qt1,phi1
    qx=0
    qy=0
    DO i=2,n1
       qx=qx+qt(i)*DCOS(phi(i))
       qy=qy+qt(i)*DSIN(phi(i))
    ENDDO
    qt(1)=DSQRT(qx*qx+qy*qy)   ! the conservation of momentum
    c1=-qx/qt(1)
    s1=-qy/qt(1)
    phi(1)=DACOS(c1)
    IF(s1.LT.0)phi(1)=2*pi-phi(1)


! define x and y1 
    e0=0
    x0=0
    DO i=1,n1
       e0=e0+qt(i)*DCOSH(y(i))
       x0=x0+qt(i)*DSINH(y(i))
    ENDDO
    q0=DSQRT(e0**2-x0**2)
 
    IF(q0/e0.GT.0)EXIT             !goto 100
  ENDDO
  xres=qq/q0

  ch1=( e0*et-x0*xt)/q0/qq
  sh1=(-x0*et+e0*xt)/q0/qq
  IF(ch1+sh1.GT.0)EXIT
ENDDO
y(1)=DLOG(ch1+sh1)
 
DO i=2,n1
   y(i)=y(i)+y(1)
ENDDO

! define momenta
k=INT(Helac_rnmy(0)*nch)+1
p4=0
p3=0
p2=0
p1=0
DO i=1,n1
   j=ia(i,ich(k))
   p(j,1)=xres*qt(i)*DCOS(phi(i))
   p(j,2)=xres*qt(i)*DSIN(phi(i))
   p(j,3)=xres*qt(i)*DSINH(y(i))
   p(j,4)=xres*qt(i)*DCOSH(y(i))
       
   p4=p4+p(j,4)
   p3=p3+p(j,3)
   p2=p2+p(j,2)
   p1=p1+p(j,1)
ENDDO

pboo(4)=p4
pboo(3)=-p3
pboo(1:2)=0
DO i=1,n1
   j=ia(i,ich(k))
   pin(1:4)=p(j,1:4)
   CALL BOOSTL(qq,pboo,pin)  ! boost to the cm frame
   p(j,1:4)=pin(1:4)
   pp_durham(j,1:4)=p(j,1:4)
ENDDO
END SUBROUTINE Phegas_durham

! define weight factor
SUBROUTINE Phegas_durham_wei(w)
REAL(KIND(1d0)),INTENT(OUT)::w
INTEGER::k1,k2
INTEGER::j,i1,i2
REAL(KIND(1d0))::pt,pti1,wg,dy
wg=0
DO k1=1,nch
   k2=ich(k1)
   w=1
   pt=0
   DO j=2,n1
      i1=ia(j  ,k2)
      i2=ia(j-1,k2)
      dy=prapidity(pp_durham(i1,1:4))-prapidity(pp_durham(i2,1:4))
      pti1=DSQRT(pp_durham(i1,1)**2+pp_durham(i1,2)**2)
      w=w*pti1/( DTANH(2*cutrap-dy)+DTANH(2*cutrap+dy) )
      pt=pt+pti1
   ENDDO
   w=w/qq**2*pt**(n1-1)/2**n1*(2*pi)**(n1-1)*(8*cutrap)**(n1-1)/gamma1
   wg=wg+1d0/w/nch
ENDDO
w=1d0/wg
END SUBROUTINE Phegas_durham_wei
END MODULE Phegas_Durham_mod
