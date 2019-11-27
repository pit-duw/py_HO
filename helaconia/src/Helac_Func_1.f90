MODULE Helac_Func_1
USE Helac_Global
IMPLICIT NONE
!INTEGER,DIMENSION(20),PUBLIC::ia=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),&
!                             ia1=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
CONTAINS
SUBROUTINE Helac_id(n1,l,ix1,i0)
! level two is finded
! level l finded
IMPLICIT NONE
INTEGER,INTENT(IN)::n1,l
INTEGER,INTENT(INOUT)::ix1
INTEGER,INTENT(OUT)::i0
INTEGER,DIMENSION(20)::ia=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
INTEGER::k,is
SAVE ia
IF(GLOBALINIT_id.EQ.0)THEN
	ia(1:20)=0
ENDIF
CALL Helac_mycombi(n1,l,ia)

IF(ia(1).NE.0)THEN
    is=0
    DO k=1,l
      is=is+2**(ia(k)-1)
    ENDDO
    i0=is
    RETURN
ENDIF 
ix1=1
END SUBROUTINE Helac_id        
      
SUBROUTINE Helac_id2(n1,l,i0,ix2,i1,i2)
! splitting is done
! returns two splitting of particles i1 and i2 to i0       
INTEGER,INTENT(IN)::n1
INTEGER,INTENT(IN)::l,i0
INTEGER,INTENT(INOUT)::ix2
INTEGER,INTENT(OUT)::i1,i2
INTEGER,DIMENSION(n1)::ic,y
INTEGER,DIMENSION(20)::ia=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
INTEGER::j=1,m,i,k,is
SAVE ia,j
m=0
IF(GLOBALINIT_id2.EQ.0)THEN
	ia(1:20)=0
	j=1
ENDIF
CALL Helac_bin(n1,i0,y)

DO i=n1,1,-1
   IF(y(i).EQ.1)THEN
      m=m+1
      ic(m)=2**(n1-i)
   ENDIF
ENDDO

! 1     continue
DO
  CALL Helac_mycombi(l,j,ia)     
  IF(ia(1).NE.0)THEN
    is=0
    DO k=1,j
      is=is+ic(ia(k))
    ENDDO
    i1=is
    i2=i0-is
    RETURN
  ENDIF 
        
  j=j+1
  IF(j.GE.l)EXIT
ENDDO
IF(j.EQ.l)THEN
  j=1
  ix2=1
ENDIF
END SUBROUTINE Helac_id2
       
SUBROUTINE Helac_id3(n1,l,i0,ix3,i1,i2,i3)
! splitting is done
! returns three splitting of particles i1,i2 and i3 to i0 
INTEGER,INTENT(IN)::n1,l,i0
INTEGER,INTENT(INOUT)::ix3
INTEGER,INTENT(OUT)::i1,i2,i3       
INTEGER,DIMENSION(20)::ic,ic1
INTEGER,DIMENSION(n1)::y
INTEGER,DIMENSION(20)::ia=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/),&
                       ia1=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
INTEGER::j=1,j1=1,ix=0,ii0,ii1,m,i,flag,is,k,flag2
SAVE ia,ia1,ic,ic1,j,j1,ix,ii0,ii1
flag2=0
IF(GLOBALINIT_id3.EQ.0)THEN
	ia(1:20)=0
	ia1(1:20)=0
	j=1
	j1=1
	ix=0
	ic(1:20)=0
	ic1(1:20)=0
	ii0=0
	ii1=0
ENDIF
DO  
! The first splitting           
  IF(ix.NE.1.OR.flag2.EQ.1)THEN                         !goto2
    IF(flag2.EQ.0)THEN
      ii0=i0
      m=0
      CALL Helac_bin(n1,i0,y)
      DO i=n1,1,-1
        IF(y(i).EQ.1)THEN
          m=m+1
          ic(m)=2**(n1-i)
        ENDIF
      ENDDO
	ENDIF

!1      continue
    flag=0
    DO
      CALL Helac_mycombi(l,j,ia)
      IF(ia(1).EQ.0)THEN
        j1=1
        j=j+1
	    flag=1
        IF(j.GE.l-1)THEN
	       EXIT
	    ENDIF
      ELSE
	    flag=0
	    EXIT
      ENDIF
    ENDDO
    IF(j.EQ.l-1.AND.flag.EQ.1)THEN
      j=1
      ix3=1
      RETURN
    ENDIF
    is=0
    DO k=1,j
      is=is+ic(ia(k))
    ENDDO
    i1=is
    ii1=i1
    i2=i0-is
! the second splitting          
    m=0
    CALL Helac_bin(n1,i2,y)
    DO i=n1,1,-1
      IF(y(i).EQ.1)THEN 
         m=m+1
         ic1(m)=2**(n1-i)
      ENDIF
    ENDDO
    ix=1
  ENDIF
 
! 2     continue 

  DO
    CALL Helac_mycombi(l-j,j1,ia1)       
    IF(ia1(1).NE.0)THEN
      is=0
      DO k=1,j1
       is=is+ic1(ia1(k))
      ENDDO
      i2=is
      i3=ii0-ii1-i2
      RETURN
    ENDIF 
        
    j1=j1+1
    IF(j1.GE.l-j)EXIT   ! goto 2
  ENDDO
  IF(j1.EQ.l-j)THEN
    ix=0
    j1=1
    flag2=1               !goto 1
  ELSE
    EXIT
  ENDIF
ENDDO
END SUBROUTINE Helac_id3

SUBROUTINE Helac_idoctet(n1,l,ix1,i0)
! level l finded for the octet in Helac_init
IMPLICIT NONE
INTEGER,INTENT(IN)::n1,l
INTEGER,INTENT(INOUT)::ix1
INTEGER,INTENT(OUT)::i0
INTEGER,DIMENSION(20)::ia=(/0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0/)
INTEGER::k,is
SAVE ia
ix1=0
IF(GLOBALINIT_idoctet.EQ.0)THEN
	ia(1:20)=0
ENDIF
CALL Helac_mycombi(n1,l,ia)

IF(ia(1).NE.0)THEN
    is=0
    DO k=1,l
      is=is+2**(ia(k)-1)
    ENDDO
    i0=is
    RETURN
ENDIF 
ix1=1
END SUBROUTINE Helac_idoctet
       
INTEGER FUNCTION Helac_ipsgn(i1,i2)
! Return the Eq.(9,10) in 0002082v2
INTEGER,INTENT(IN)::i1,i2       
INTEGER,DIMENSION(n)::y1,y2
INTEGER::k,is,i,j
CALL Helac_bin(n,i1,y1)
CALL Helac_bin(n,i2,y2)
DO k=1,n
   IF(ifl(k).GE.31)THEN
      y1(n-k+1)=0
      y2(n-k+1)=0
   ENDIF
ENDDO
is=0
DO i=n,2,-1
   IF(y2(i).EQ.1)THEN
     DO j=1,i-1
         is=is+y1(j)
     ENDDO
   ENDIF
ENDDO
       
Helac_ipsgn=(-1)**is
END FUNCTION Helac_ipsgn
       
SUBROUTINE Helac_bin(n1,x,y)
!      converts to a binary number
!      y -result
!      x -input number
!      n -levels     
INTEGER,INTENT(IN)::x,n1
INTEGER,DIMENSION(n1),INTENT(OUT)::y
INTEGER::X1,J
y=0
IF(x.EQ.0)RETURN
X1=X
DO J=n1,1,-1
   y(J)=X1-INT(X1/2)*2
   IF((X1/2).LT.1.) EXIT
   X1=INT(X1/2)
ENDDO
END SUBROUTINE Helac_bin
       
FUNCTION Helac_ifactorial(n1)
INTEGER::Helac_ifactorial
INTEGER,INTENT(IN)::n1
INTEGER::k
Helac_ifactorial=1
DO k=2,n1
   Helac_ifactorial=Helac_ifactorial*k
ENDDO
END FUNCTION Helac_ifactorial
      
SUBROUTINE Helac_BIN3(N1,X,Y)
!      converts to a trinary number
INTEGER,INTENT(IN)::X
INTEGER,INTENT(IN)::N1
INTEGER,DIMENSION(N1)::Y
INTEGER::J,X1
Y=0
X1=X
DO J=N1,1,-1
   Y(J)=X1-INT(X1/3)*3
   IF((X1/3).LT.1.) EXIT
   X1=INT(X1/3)
ENDDO
END SUBROUTINE Helac_BIN3

SUBROUTINE Helac_BIN6(N1,X,Y)
!      converts to a sixnary number
INTEGER,INTENT(IN)::X
INTEGER,INTENT(IN)::N1
INTEGER,DIMENSION(N1)::Y
INTEGER::J,X1
Y=0
X1=X
DO J=N1,1,-1
   Y(J)=X1-INT(X1/6)*6
   IF((X1/6).LT.1.) EXIT
   X1=INT(X1/6)
ENDDO
END SUBROUTINE Helac_BIN6

SUBROUTINE Helac_BIN18(N1,X,Y)
!      converts to a eigthteennary number
INTEGER,INTENT(IN)::X
INTEGER,INTENT(IN)::N1
INTEGER,DIMENSION(N1)::Y
INTEGER::J,X1
Y=0
X1=X
DO J=N1,1,-1
   Y(J)=X1-INT(X1/18)*18
   IF((X1/18).LT.1.) EXIT
   X1=INT(X1/18)
ENDDO
END SUBROUTINE Helac_BIN18

!*********************************************************
!     ipol(k)=helicity,ipol(k)=2:rigth-handed            *
!                      ipol(k)=1:left-handed             *
!                      ipol(k)=3:long.                   *
!*********************************************************      
SUBROUTINE Helac_sethel(iflag1,istop,inu1,inegl)
INTEGER,INTENT(IN)::iflag1,inegl
INTEGER,INTENT(INOUT)::istop
INTEGER,INTENT(INOUT)::inu1
INTEGER,DIMENSION(20)::lp=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
INTEGER,DIMENSION(:,:),ALLOCATABLE::nhl
INTEGER,DIMENSION(:),ALLOCATABLE::y2,y3,y22,y6,y62,y18
INTEGER::init=0,l1,l2,l3,l22,l6,l62,l18,k,istat,inu,i2,i3,ic2,ic3,ltotal,&
i22,i6,i62,i18,ic22,ic6,ic62,ic18
SAVE nhl,y2,y3,y22,y6,y62,y18,init
l1=0
l2=0
l3=0
l22=0  ! 1S0 state
l6=0   ! 3S1 state
l62=0  ! 1P1 state
l18=0  ! 3PJ state

IF(GLOBALINIT_Sethel.EQ.0)THEN
	init=0
	lp(1:20)=1
ENDIF

IF(inegl.EQ.1)RETURN
!------------------------------------------------------------
!          START A
!
! In this part supply your helicity configuration
!
!------------------------------------------------------------
IF(iflag1.EQ.1)THEN
   WRITE(*,*)'You must give to me the required helicities'
   WRITE(*,*)'The number of particles is:  ',nhad
   REWIND(20)
   READ(20,*)(ipol(k),k=1,nhad)
   istop=1
   RETURN
ENDIF
!------------------------------------------------------------
!          END A
!------------------------------------------------------------
!------------------------------------------------------------
!          START B
!
! Define all possible helicity configurations
!
!------------------------------------------------------------
IF(init.EQ.0)THEN
	DO k=1,nhad
		lp(k)=1
		IF(iranhel.EQ.0)THEN
			IF(iflh(k).GE.-12.AND.iflh(k).LE.12)lp(k)=2
			IF(iflh(k).EQ.31.OR.iflh(k).EQ.35)lp(k)=2            !change 5.4.00
			IF(iflh(k).GE.32.AND.iflh(k).LE.34)lp(k)=3
		ENDIF
		IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
			IF(QN1S0F(k))lp(k)=2
			IF(QN3S1F(k).OR.QN1P1F(k))lp(k)=6
			IF(QN3PJF(k))lp(k)=18
		ELSEIF(iranhel.EQ.2)THEN
			IF(QN1S0F(k).OR.QN1P1F(k))lp(k)=2
			IF(QN3S1F(k).OR.QN3PJF(k))lp(k)=6
		ELSEIF(iranhel.EQ.3)THEN
			IF(ABS(iflh(k)).GT.100)lp(k)=2
		ENDIF
	ENDDO
	DO k=1,nhad
		IF(iranhel.EQ.0)THEN
			IF(lp(k).EQ.1)l1=l1+1
			IF(lp(k).EQ.2.AND.(ABS(iflh(k)).LE.100))l2=l2+1
			IF(lp(k).EQ.3)l3=l3+1
		ENDIF
		IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
			IF(QN1S0F(k))l22=l22+1
			IF(QN3S1F(k))l6=l6+1
			IF(QN1P1F(k))l62=l62+1
			IF(QN3PJF(k))l18=l18+1
		ELSEIF(iranhel.EQ.2)THEN
			IF(QN1S0F(k).OR.QN1P1F(k))l22=l22+1
			IF(QN3S1F(k).OR.QN3PJF(k))l6=l6+1
		ELSEIF(iranhel.EQ.3)THEN
			IF(lp(k).EQ.2)l2=l2+1
		ENDIF
	ENDDO
!  WRITE(*,*)l1,l2,l3
	init=1
	IF(ALLOCATED(nhl))THEN
		DEALLOCATE(nhl)
	ENDIF
	IF(ALLOCATED(y2))THEN
		DEALLOCATE(y2)
	ENDIF
	IF(ALLOCATED(y3))THEN
		DEALLOCATE(y3)
	ENDIF
	IF(ALLOCATED(y22))THEN
		DEALLOCATE(y22)
	ENDIF
	IF(ALLOCATED(y6))THEN
		DEALLOCATE(y6)
	ENDIF
	IF(ALLOCATED(y62))THEN
		DEALLOCATE(y62)
	ENDIF
	IF(ALLOCATED(y18))THEN
		DEALLOCATE(y18)
	ENDIF
	ltotal=2**l2*3**l3*2**l22*6**(l6+l62)*18**l18
	ALLOCATE(nhl(ltotal,nhad),y2(l2),y3(l3),&
	y22(l22),y6(l6),y62(l62),y18(l18),STAT=istat)

!  WRITE(*,*)istat
	IF(istat.NE.0)THEN
		WRITE(*,*)'warning: allocation is not working properly'
		STOP
	ENDIF

	DO inu=1,ltotal
		IF(iranhel.EQ.0)THEN
			i3=MOD((inu-1)/2**l2,3**l3)
			i2=MOD(inu-1,2**l2)
			CALL  Helac_bin(l2,i2,y2)
			CALL  Helac_bin3(l3,i3,y3)
		ENDIF
		IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
			i22=MOD((inu-1)/(2**l2*3**l3),2**l22)
			i6=MOD((inu-1)/(2**(l2+l22)*3**l3),6**l6)
			i62=MOD((inu-1)/(2**(l2+l22)*3**l3*6**l6),6**l62)
			i18=(inu-1)/(2**(l2+l22)*3**l3*6**(l6+l62))
			CALL  Helac_bin(l22,i22,y22)
			CALL  Helac_bin6(l6,i6,y6)
			CALL  Helac_bin6(l62,i62,y62)
			CALL  Helac_bin18(l18,i18,y18)
		ELSEIF(iranhel.EQ.2)THEN
			i22=MOD((inu-1)/(2**l2*3**l3),2**l22)
			i6=MOD((inu-1)/(2**(l2+l22)*3**l3),6**l6)
			CALL  Helac_bin(l22,i22,y22)
			CALL  Helac_bin6(l6,i6,y6)
		ELSEIF(iranhel.EQ.3)THEN
			i2=MOD(inu-1,2**l2)
			CALL Helac_bin(l2,i2,y2)
		ENDIF
		ic2=0
		ic3=0
		ic22=0
		ic6=0
		ic62=0
		ic18=0
		DO k=1,nhad
			IF(lp(k).EQ.1)ipol(k)=1
			IF(lp(k).EQ.2.AND.(ABS(iflh(k)).LE.100))THEN
				ic2=ic2+1
				ipol(k)=y2(ic2)+1
			ENDIF
			IF(lp(k).EQ.3)THEN
				ic3=ic3+1
				ipol(k)=y3(ic3)+1
			ENDIF
			IF(iranhel.EQ.0.OR.iranhel.EQ.1)THEN
				IF(QN1S0F(k))THEN
					ic22=ic22+1
					ipol(k)=y22(ic22)+1
				ENDIF
				IF(QN3S1F(k))THEN
					ic6=ic6+1
					ipol(k)=y6(ic6)+1
				ENDIF
				IF(QN1P1F(k))THEN
					ic62=ic62+1
					ipol(k)=y62(ic62)+1
				ENDIF
				IF(QN3PJF(k))THEN
					ic18=ic18+1
					ipol(k)=y18(ic18)+1
				ENDIF
			ELSEIF(iranhel.EQ.2)THEN
				IF(QN1S0F(k).OR.QN1P1F(k))THEN
					ic22=ic22+1
					ipol(k)=y22(ic22)+1
				ENDIF
				IF(QN3S1F(k).OR.QN3PJF(k))THEN
					ic6=ic6+1
					ipol(k)=y6(ic6)+1
				ENDIF
			ELSEIF(iranhel.EQ.3)THEN
				IF(lp(k).EQ.2)THEN
					ic2=ic2+1
					ipol(k)=y2(ic2)+1
				ENDIF
			ENDIF
		ENDDO

		nhl(inu,1:nhad)=ipol(1:nhad)

	ENDDO
	nhc=ltotal
ENDIF
!------------------------------------------------------------
!          END B
!------------------------------------------------------------
!------------------------------------------------------------
!          START C
!
!          Normal operation
!
!------------------------------------------------------------
IF(inegl.EQ.2)THEN
     nhc=nhc-1
     nhl(inu1:nhc,1:nhad)=nhl(inu1+1:nhc+1,1:nhad)
     inu1=inu1-1
     RETURN
ENDIF

inu1=inu1+1
ipol(1:nhad)=nhl(inu1,1:nhad)
IF(inu1.GE.nhc)istop=1
RETURN
!------------------------------------------------------------
!          END C
!------------------------------------------------------------
END SUBROUTINE Helac_sethel

SUBROUTINE Helac_sethel2(iflag1,istop,inu1,inegl,imode1)
INTEGER,INTENT(IN)::iflag1,inegl,imode1
INTEGER,INTENT(INOUT)::istop
INTEGER,INTENT(INOUT)::inu1
INTEGER,DIMENSION(20)::lp=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
!INTEGER,DIMENSION(20)::sp=(/1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1/)
INTEGER,DIMENSION(:,:),ALLOCATABLE::nhl,phyhl,map_phyhl2nhl
INTEGER,DIMENSION(:),ALLOCATABLE::y2,y3,y22,y6,y62,y18,map_nhl2phyhl
INTEGER::init=0,l1,l2,l3,l22,l6,l62,l18,k,istat,inu,i2,i3,ic2,ic3,ltotal,&
i22,i6,i62,i18,ic22,ic6,ic62,ic18,stotal,inum
SAVE nhl,phyhl,map_phyhl2nhl,map_nhl2phyhl,y2,y3,y22,y6,y62,y18,init,ltotal,stotal
l1=0
l2=0
l3=0
l22=0  ! 1S0 state
l6=0   ! 3S1 state
l62=0  ! 1P1 state
l18=0  ! 3PJ state                                                              
IF(inegl.EQ.1)RETURN
IF(init.EQ.0)THEN
   DO k=1,nhad
      lp(k)=1
      IF(iflh(k).GE.-12.AND.iflh(k).LE.12)lp(k)=2
      IF(iflh(k).EQ.31.OR.iflh(k).EQ.35)lp(k)=2
      IF(iflh(k).GE.32.AND.iflh(k).LE.34)lp(k)=3
      IF(QN1S0F(k))lp(k)=2
      IF(QN3S1F(k).OR.QN1P1F(k))lp(k)=6
      IF(QN3PJF(k))lp(k)=18
   ENDDO
   DO k=1,nhad
      IF(lp(k).EQ.1)l1=l1+1
      IF(lp(k).EQ.2.AND.(ABS(iflh(k)).LE.100))l2=l2+1
      IF(lp(k).EQ.3)l3=l3+1
      IF(QN1S0F(k))l22=l22+1
      IF(QN3S1F(k))l6=l6+1
      IF(QN1P1F(k))l62=l62+1
      IF(QN3PJF(k))l18=l18+1
   ENDDO
   init=1
   IF(ALLOCATED(nhl))THEN
      DEALLOCATE(nhl)
   ENDIF
   IF(ALLOCATED(phyhl))THEN
      DEALLOCATE(phyhl)
   ENDIF
   IF(ALLOCATED(map_phyhl2nhl))THEN
      DEALLOCATE(map_phyhl2nhl)
   ENDIF
   IF(ALLOCATED(map_nhl2phyhl))THEN
      DEALLOCATE(map_nhl2phyhl)
   ENDIF
   IF(ALLOCATED(y2))THEN
      DEALLOCATE(y2)
   ENDIF
   IF(ALLOCATED(y3))THEN
      DEALLOCATE(y3)
   ENDIF
   IF(ALLOCATED(y22))THEN
      DEALLOCATE(y22)
   ENDIF
   IF(ALLOCATED(y6))THEN
      DEALLOCATE(y6)
   ENDIF
   IF(ALLOCATED(y62))THEN
      DEALLOCATE(y62)
   ENDIF
   IF(ALLOCATED(y18))THEN
      DEALLOCATE(y18)
   ENDIF
   ltotal=2**l2*3**l3*3**(l6+l62)*9**l18
   stotal=2**(l22+l6+l62+l18)
   ALLOCATE(nhl(ltotal*stotal,nhad),y2(l2),y3(l3),&
        y22(l22),y6(l6),y62(l62),y18(l18),STAT=istat)
   ALLOCATE(phyhl(ltotal,nhad),map_phyhl2nhl(ltotal,0:stotal),&
        map_nhl2phyhl(ltotal*stotal))
   IF(istat.NE.0)THEN
      WRITE(*,*)'warning: allocation is not working properly in sethel2'
      STOP
   ENDIF
   phyhl(1:ltotal,1:nhad)=0
   map_phyhl2nhl(1:ltotal,0:stotal)=0
   map_nhl2phyhl(1:ltotal*stotal)=0
   DO inu=1,ltotal*stotal
      i3=MOD((inu-1)/2**l2,3**l3)
      i2=MOD(inu-1,2**l2)
      CALL  Helac_bin(l2,i2,y2)
      CALL  Helac_bin3(l3,i3,y3)
      i22=MOD((inu-1)/(2**l2*3**l3),2**l22)
      i6=MOD((inu-1)/(2**(l2+l22)*3**l3),6**l6)
      i62=MOD((inu-1)/(2**(l2+l22)*3**l3*6**l6),6**l62)
      i18=(inu-1)/(2**(l2+l22)*3**l3*6**(l6+l62))
      CALL  Helac_bin(l22,i22,y22)
      CALL  Helac_bin6(l6,i6,y6)
      CALL  Helac_bin6(l62,i62,y62)
      CALL  Helac_bin18(l18,i18,y18)
      ic2=0
      ic3=0
      ic22=0
      ic6=0
      ic62=0
      ic18=0
      DO k=1,nhad
         IF(lp(k).EQ.1)ipol(k)=1
         IF(lp(k).EQ.2.AND.(ABS(iflh(k)).LE.100))THEN
            ic2=ic2+1
            ipol(k)=y2(ic2)+1
         ENDIF
         IF(lp(k).EQ.3)THEN
            ic3=ic3+1
            ipol(k)=y3(ic3)+1
         ENDIF
         IF(QN1S0F(k))THEN
            ic22=ic22+1
            ipol(k)=y22(ic22)+1
         ENDIF
         IF(QN3S1F(k))THEN
            ic6=ic6+1
            ipol(k)=y6(ic6)+1
         ENDIF
         IF(QN1P1F(k))THEN
            ic62=ic62+1
            ipol(k)=y62(ic62)+1
         ENDIF
         IF(QN3PJF(k))THEN
            ic18=ic18+1
            ipol(k)=y18(ic18)+1
         ENDIF
      ENDDO
      CALL PhysicalPol(inum)
      GLOBALINIT_physicalpol=1
      map_nhl2phyhl(inu)=inum
      ic18=map_phyhl2nhl(inum,0)
      IF(ic18.EQ.0)THEN
         DO k=1,nhad
            IF(ABS(iflh(k)).GT.100)THEN
               phyhl(inum,k)=(ipol(k)-1)/2+1
            ELSE
               phyhl(inum,k)=ipol(k)
            ENDIF
         ENDDO
      ENDIF
      map_phyhl2nhl(inum,0)=ic18+1
      map_phyhl2nhl(inum,ic18+1)=inu
      nhl(inu,1:nhad)=ipol(1:nhad)
   ENDDO
   nhc=ltotal*stotal
   nphyhc=ltotal
END IF
IF(inegl.EQ.2)THEN
   nhc=nhc-1
   inum=map_nhl2phyhl(inu1)
   DO k=inu1+1,nhc+1
      map_nhl2phyhl(k-1)=map_nhl2phyhl(k)
   ENDDO
   DO k=1,nphyhc
      ic18=map_phyhl2nhl(k,0)
      DO ic6=1,ic18
         IF(map_phyhl2nhl(k,ic6).GT.inu1)THEN
            map_phyhl2nhl(k,ic6)=map_phyhl2nhl(k,ic6)-1
         ENDIF
      ENDDO
   ENDDO
   ic18=map_phyhl2nhl(inum,0)
   DO k=1,ic18
      IF(map_phyhl2nhl(inum,k).EQ.inu1)THEN
         map_phyhl2nhl(inum,0)=map_phyhl2nhl(inum,0)-1
         IF(ic18-1.GT.0)THEN
            map_phyhl2nhl(inum,k:ic18-1)=map_phyhl2nhl(inum,k+1:ic18)
         ENDIF
         EXIT
      ENDIF
   ENDDO
   IF(ic18-1.EQ.0)THEN
      nphyhc=nphyhc-1
      phyhl(inum:nphyhc,1:nhad)=phyhl(inum+1:nphyhc+1,1:nhad)
      DO k=1,nhc
         IF(map_nhl2phyhl(k).GT.inum)THEN
            map_nhl2phyhl(k)=map_nhl2phyhl(k)-1
         ENDIF
      ENDDO
      map_phyhl2nhl(inum:nphyhc,0:stotal)=map_phyhl2nhl(inum+1:nphyhc+1,0:stotal)
   ENDIF
   nhl(inu1:nhc,1:nhad)=nhl(inu1+1:nhc+1,1:nhad)
   inu1=inu1-1
   RETURN
ENDIF

IF(imode1.EQ.0)THEN
   inu1=inu1+1
   ipol(1:nhad)=nhl(inu1,1:nhad)
   IF(inu1.GE.nhc)istop=1
ELSE
   inu1=inu1+1
   inum=INT(RMCoH*nphyhc)+1
   ic18=map_phyhl2nhl(inum,0)
   ic6=map_phyhl2nhl(inum,inu1)
   ipol(1:nhad)=nhl(ic6,1:nhad)
   IF(inu1.GE.ic18)istop=1   
ENDIF
END SUBROUTINE Helac_sethel2

SUBROUTINE Helac_setncc(nq)
! number of color configurations
INTEGER,INTENT(OUT)::nq
INTEGER::k
nq=0
DO k=1,n
   ! incoming
   IF(io(k).EQ.1)THEN
     IF(ifl(k).EQ.3.OR.ifl(k).EQ.4.OR.ifl(k).EQ.7&
	 .OR.ifl(k).EQ.8.OR.ifl(k).EQ.11.OR.ifl(k).EQ.12)nq=nq+1
   ! outgoing 
   ELSEIF(io(k).EQ.-1)THEN
     IF(ifl(k).EQ.-3.OR.ifl(k).EQ.-4.OR.ifl(k).EQ.-7& 
      .OR.ifl(k).EQ.-8.OR.ifl(k).EQ.-11.OR.ifl(k).EQ.-12)nq=nq+1 
   ENDIF
   IF(ifl(k).EQ.35)nq=nq+1
ENDDO
END SUBROUTINE Helac_setncc
       
SUBROUTINE Helac_setcc(iq)
INTEGER,DIMENSION(20),INTENT(IN)::iq
INTEGER::ic,ic1,k
ic=0
ic1=0
icol(1:20,1:2)=0
DO k=1,n
   ! incoming
   IF(io(k).EQ. 1)THEN
       IF(ifl(k).EQ.3.OR.ifl(k).EQ.4.OR.ifl(k).EQ.7& 
        .OR.ifl(k).EQ.8.OR.ifl(k).EQ.11.OR.ifl(k).EQ.12)THEN
          ic=ic+1
          icol(k,1)=iq(ic)
       ELSEIF(ifl(k).EQ.-3.OR.ifl(k).EQ.-4.OR.ifl(k).EQ.-7& 
        .OR.ifl(k).EQ.-8.OR.ifl(k).EQ.-11.OR.ifl(k).EQ.-12)THEN
          ic1=ic1+1
          icol(k,2)=ic1
       ENDIF
   ! outgoing  
   ELSEIF(io(k).EQ.-1)THEN
       IF(ifl(k).EQ.-3.OR.ifl(k).EQ.-4.OR.ifl(k).EQ.-7&
	    .OR.ifl(k).EQ.-8.OR.ifl(k).EQ.-11.OR.ifl(k).EQ.-12)THEN
          ic=ic+1
          icol(k,1)=iq(ic)
       ELSEIF(ifl(k).EQ.3.OR.ifl(k).EQ.4.OR.ifl(k).EQ.7 &
        .OR.ifl(k).EQ.8.OR.ifl(k).EQ.11.OR.ifl(k).EQ.12)THEN
          ic1=ic1+1
          icol(k,2)=ic1
       ENDIF
   ENDIF
       
   IF(ifl(k).EQ.35)THEN
	
	ic=ic+1
	ic1=ic1+1
	icol(k,1)=iq(ic)
    icol(k,2)=ic1

   ENDIF
ENDDO
END SUBROUTINE Helac_setcc
       
SUBROUTINE Helac_PERMU(IA,n1)
INTEGER,INTENT(IN)::n1
INTEGER,DIMENSION(n1),INTENT(INOUT)::IA
INTEGER::i,K1,K,flag,KN,L,IB,L1
IF(n1.LE.0) RETURN
IF(IA(1).EQ.0) THEN
   DO i = 1,n1
    IA(i)=i
   ENDDO
   IF(n1.EQ.1) IA(1)=0
ELSE
   flag=0
   DO K1 = n1,2,-1
     K=K1
     IF(IA(K-1).LT.IA(K))THEN
	   flag=1
	   EXIT
	 ENDIF
   ENDDO
   IF(flag.EQ.0)THEN
     IA(1)=0
     RETURN
   ENDIF
   KN=K+n1
   DO L = K,KN/2
     IB=IA(KN-L)
     IA(KN-L)=IA(L)
     IA(L)=IB
   ENDDO
   DO L1 = K,n1
     L=L1
     IF(IA(L).GT.IA(K-1)) EXIT
   ENDDO
   IB=IA(K-1)
   IA(K-1)=IA(L)
   IA(L)=IB
ENDIF
END SUBROUTINE Helac_PERMU
       
SUBROUTINE Helac_countzeros(n1,iv,iz)
INTEGER,INTENT(IN)::n1
INTEGER,DIMENSION(n1),INTENT(IN)::iv
INTEGER,INTENT(OUT)::iz
INTEGER::k
iz=0
DO k=1,n1
   IF(iv(k).EQ.0)iz=iz+1
ENDDO
END SUBROUTINE Helac_countzeros
      
SUBROUTINE Helac_average(ah,ac,s)
! we don't include the symmetry factor for double quarkonium
REAL(KIND=DBL),INTENT(OUT)::ah,ac,s
INTEGER,DIMENSION(-12:44)::is
INTEGER::init=0,k,k1,k2
INTEGER,DIMENSION(2,6)::iscc,isbb,iscb,isbc
SAVE init
is(-12:44)=0
iscc(1:2,1:6)=0
isbb(1:2,1:6)=0
iscb(1:2,1:6)=0
isbc(1:2,1:6)=0
IF(GLOBALINIT_average.EQ.0)THEN
	init=0
ENDIF
IF(init.EQ.1)RETURN
ah=dnou(1)
ac=dnou(1)
s=dnou(1)
              
DO k=1,n
   ! incoming
   IF(io(k).EQ. 1)THEN
      IF(ifl(k).LE.31)ah=ah/dnou(2)
      IF(ifl(k).EQ.35)ah=ah/dnou(2)                  !change 5.4.00
      IF(ifl(k).GE.32.AND.ifl(k).LE.34)ah=ah/dnou(3)
      IF(IABS(ifl(k)).EQ.3.OR.IABS(ifl(k)).eq.4.OR.&
        IABS(ifl(k)).EQ.7.OR.IABS(ifl(k)).EQ.8.OR.&
        IABS(ifl(k)).EQ.11.OR.IABS(ifl(k)).EQ.12)ac=ac/dnou(3)
      IF(ifl(k).EQ.35)ac=ac/dnou(8)                  !change 5.4.00
   ENDIF
ENDDO
       
DO k=3,nhad
   ! outgoing
   IF(ABS(iflh(k)).LE.100)THEN
		is(iflh(k))=is(iflh(k))+1
   ENDIF
   IF(iflh(k).EQ.441001)iscc(1,1)=iscc(1,1)+1
   IF(iflh(k).EQ.441008)iscc(2,1)=iscc(2,1)+1
   IF(iflh(k).EQ.443011)iscc(1,2)=iscc(1,2)+1
   IF(iflh(k).EQ.443018)iscc(2,2)=iscc(2,2)+1
   IF(iflh(k).EQ.441111)iscc(1,3)=iscc(1,3)+1
   IF(iflh(k).EQ.441118)iscc(2,3)=iscc(2,3)+1
   IF(iflh(k).EQ.443101)iscc(1,4)=iscc(1,4)+1
   IF(iflh(k).EQ.443108)iscc(2,4)=iscc(2,4)+1
   IF(iflh(k).EQ.443111)iscc(1,5)=iscc(1,5)+1
   IF(iflh(k).EQ.443118)iscc(2,5)=iscc(2,5)+1
   IF(iflh(k).EQ.443121)iscc(1,6)=iscc(1,6)+1
   IF(iflh(k).EQ.443128)iscc(2,6)=iscc(2,6)+1
   IF(iflh(k).EQ.551001)isbb(1,1)=isbb(1,1)+1
   IF(iflh(k).EQ.551008)isbb(2,1)=isbb(2,1)+1
   IF(iflh(k).EQ.553011)isbb(1,2)=isbb(1,2)+1
   IF(iflh(k).EQ.553018)isbb(2,2)=isbb(2,2)+1
   IF(iflh(k).EQ.551111)isbb(1,3)=isbb(1,3)+1
   IF(iflh(k).EQ.551118)isbb(2,3)=isbb(2,3)+1
   IF(iflh(k).EQ.553101)isbb(1,4)=isbb(1,4)+1
   IF(iflh(k).EQ.553108)isbb(2,4)=isbb(2,4)+1
   IF(iflh(k).EQ.553111)isbb(1,5)=isbb(1,5)+1
   IF(iflh(k).EQ.553118)isbb(2,5)=isbb(2,5)+1
   IF(iflh(k).EQ.553121)isbb(1,6)=isbb(1,6)+1
   IF(iflh(k).EQ.553128)isbb(2,6)=isbb(2,6)+1
   IF(iflh(k).EQ.451001)iscb(1,1)=iscb(1,1)+1
   IF(iflh(k).EQ.451008)iscb(2,1)=iscb(2,1)+1
   IF(iflh(k).EQ.453011)iscb(1,2)=iscb(1,2)+1
   IF(iflh(k).EQ.453018)iscb(2,2)=iscb(2,2)+1
   IF(iflh(k).EQ.451111)iscb(1,3)=iscb(1,3)+1
   IF(iflh(k).EQ.451118)iscb(2,3)=iscb(2,3)+1
   IF(iflh(k).EQ.453101)iscb(1,4)=iscb(1,4)+1
   IF(iflh(k).EQ.453108)iscb(2,4)=iscb(2,4)+1
   IF(iflh(k).EQ.453111)iscb(1,5)=iscb(1,5)+1
   IF(iflh(k).EQ.453118)iscb(2,5)=iscb(2,5)+1
   IF(iflh(k).EQ.453121)iscb(1,6)=iscb(1,6)+1
   IF(iflh(k).EQ.453128)iscb(2,6)=iscb(2,6)+1
   IF(iflh(k).EQ.-451001)isbc(1,1)=isbc(1,1)+1
   IF(iflh(k).EQ.-451008)isbc(2,1)=isbc(2,1)+1
   IF(iflh(k).EQ.-453011)isbc(1,2)=isbc(1,2)+1
   IF(iflh(k).EQ.-453018)isbc(2,2)=isbc(2,2)+1
   IF(iflh(k).EQ.-451111)isbc(1,3)=isbc(1,3)+1
   IF(iflh(k).EQ.-451118)isbc(2,3)=isbc(2,3)+1
   IF(iflh(k).EQ.-453101)isbc(1,4)=isbc(1,4)+1
   IF(iflh(k).EQ.-453108)isbc(2,4)=isbc(2,4)+1
   IF(iflh(k).EQ.-453111)isbc(1,5)=isbc(1,5)+1
   IF(iflh(k).EQ.-453118)isbc(2,5)=isbc(2,5)+1
   IF(iflh(k).EQ.-453121)isbc(1,6)=isbc(1,6)+1
   IF(iflh(k).EQ.-453128)isbc(2,6)=isbc(2,6)+1
ENDDO
      
DO k=-12,44
   s=s/Helac_ifactorial(is(k))
ENDDO
DO k1=1,2
	DO k2=1,6
		s=s/(Helac_ifactorial(iscc(k1,k2))*Helac_ifactorial(isbb(k1,k2))&
		*Helac_ifactorial(iscb(k1,k2))*Helac_ifactorial(isbc(k1,k2)))
	ENDDO
ENDDO
init=1
END SUBROUTINE Helac_average
       
INTEGER FUNCTION Helac_level(n1,i)
INTEGER,INTENT(IN)::n1,i       
INTEGER,DIMENSION(n1)::y
INTEGER::j
Helac_level=0
IF(i.EQ.0)RETURN
CALL Helac_bin(n1,i,y)
DO j=1,n1
   Helac_level=Helac_level+y(j)
ENDDO
END FUNCTION Helac_level
      
SUBROUTINE Helac_MYCOMBI(N1,L,IA)
INTEGER,INTENT(IN)::L,N1
INTEGER,DIMENSION(*),INTENT(INOUT)::IA
INTEGER::I,I1
IF(IA(1).EQ. 0) THEN
   DO I = 1,L
      IA(I)=I
   ENDDO
   IA(L+1)=0
ELSE
   DO I1 = 1,N1
      I=I1
      IF(IA(I+1) .NE. IA(I)+1) EXIT
      IA(I)=I
   ENDDO
   IA(I)=IA(I)+1
   IF(IA(L) .EQ. N1+1) IA(1)=0
ENDIF
END SUBROUTINE Helac_MYCOMBI 
      
SUBROUTINE Helac_eat(n1,ic1,ic2,icy)
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!    THE EATING ALGORITHM
!
!  k1 is the position of i1 in ic1 array
!  then look for the position of i1 in
!  the ic2 array and put it to k2
!  if k2=k the original position
!  delta's have been formed a closed loop icy-> icy+1
!  otherwise put k1=k2 and continue
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
INTEGER,INTENT(IN)::n1
INTEGER,INTENT(OUT)::icy
INTEGER,DIMENSION(n1),INTENT(IN)::ic1,ic2
INTEGER,DIMENSION(n1)::iflag1
INTEGER::k,k1,i ,k2,i1,flag
k=0
iflag1=0  
icy=0
flag=0
DO
 DO
  DO
    IF(flag.EQ.0)THEN
       k=k+1
       IF(k.GT.n1)RETURN
       k1=k
	ENDIF

    IF(iflag1(k1).EQ.0)THEN
       i1=ic1(k1)
       iflag1(k1)=1
	   EXIT
    ENDIF
	flag=0
  ENDDO
      
  DO i=1,n1
    IF(ic2(i).EQ.i1)EXIT
  ENDDO
  k2=i                 
  IF(k2.EQ.k)THEN
    icy=icy+1
  ELSE
    EXIT
  ENDIF
  flag=0               
 ENDDO
       
 k1=k2
 flag=1

ENDDO      
END SUBROUTINE Helac_eat

SUBROUTINE Helac_mtime()
CHARACTER(len=13)::cd0, ct0, cz0
INTEGER,DIMENSION(8)::iv0
CALL DATE_AND_TIME(cd0,ct0,cz0,iv0)
WRITE(*,*)'TIME:'
WRITE(*,*)'YEAR',iv0(1),', MONTH',iv0(2),', DAY',iv0(3)
WRITE(*,*)'HOUR',iv0(5),', MINUTE',iv0(6),', SECOND',iv0(7)
WRITE(*,*)'MILLISECOND',iv0(8)
END SUBROUTINE Helac_mtime

! Returns the imin,isec,i100th to positive, while ihr can be negative
! e.g.  1 (h) -2 (min) 0 (s) 0 (10ms) -> 0 (h) 58 (min) 0 (s) 0 (10ms) 
SUBROUTINE Helac_OTime(ihr,imin,isec,i100th)
IMPLICIT NONE
INTEGER,INTENT(INOUT)::ihr,imin,isec,i100th
INTEGER::iigg,h1,h2
iigg=i100th+100*(isec+60*(imin+60*ihr))
h1=MOD(iigg,360000)
ihr=(iigg-h1)/360000
iigg=h1
h1=MOD(iigg,6000)
imin=(iigg-h1)/6000
iigg=h1
h1=MOD(iigg,100)
isec=(iigg-h1)/100
i100th=h1
END SUBROUTINE Helac_OTime

SUBROUTINE Vegas_speed(nmc,iday,ihr,imin,isec,i100th)
INTEGER,INTENT(IN)::nmc
INTEGER,INTENT(INOUT)::iday,ihr,imin,isec,i100th
INTEGER::iigg,h1
REAL::te
iigg=i100th+100*(isec+60*(imin+60*(ihr+24*iday)))
te=dnou(iigg*10**4)/dnou(nmc)
WRITE(*,*)" "
IF(te.LT.1000)THEN
	WRITE(*,*)"speed=",te,"us/point"
ELSEIF(te.GE.1000.AND.te.LT.10**6)THEN
	WRITE(*,*)"speed=",te/1000,"ms/point"
ELSEIF(te.GE.10**6.AND.te.LT.10**9)THEN
	WRITE(*,*),"speed=",te/10**6,"s/point"
ENDIF
h1=MOD(iigg,100*60*60*24)
iday=(iigg-h1)/(100*60*60*24)
iigg=h1
h1=MOD(iigg,360000)
ihr=(iigg-h1)/360000
iigg=h1
h1=MOD(iigg,6000)
imin=(iigg-h1)/6000
iigg=h1
h1=MOD(iigg,100)
isec=(iigg-h1)/100
i100th=h1
WRITE(*,*)"total time:"
WRITE(*,*)iday,"d ;",ihr,"h ;",imin,"m ;",isec,"s ;",i100th,"cms "
WRITE(*,*)" "
END SUBROUTINE Vegas_speed


SUBROUTINE Helac_mytime(i)
IMPLICIT NONE
INTEGER,INTENT(IN)::i
INTEGER::time
CHARACTER(len=13)::cd1, ct1, cz1
INTEGER,DIMENSION(8)::iv1
INTEGER::iclock1=0, irate1, imax1
CHARACTER(len=13)::cd0, ct0, cz0
INTEGER,DIMENSION(8)::iv0
INTEGER::iclock0=0, irate0, imax0

cd1='             '
ct1=cd1
cz1=cd1
cd0='             '
ct0=cd0
cz0=cd0
IF(i.EQ.0)THEN
    CALL SYSTEM_CLOCK(iclock0,irate0,imax0)
    WRITE(*,*)'----REST------>',iclock0-iclock1
ELSEIF(i.EQ.1)THEN
    CALL DATE_AND_TIME(cd1,ct1,cz1,iv1)
    CALL SYSTEM_CLOCK(iclock1,irate1,imax1)
    WRITE(*,*)'----ME  ------>',iclock1-iclock0
ENDIF
END SUBROUTINE Helac_mytime

FUNCTION dnou(i)
IMPLICIT NONE
INTEGER,INTENT(IN)::i
REAL(KIND=DBL)::dnou
dnou=DBLE(i)
END FUNCTION dnou

SUBROUTINE Helac_mypi(p)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(OUT)::p
p=DACOS(-dnou(1))
END SUBROUTINE Helac_mypi

SUBROUTINE Helac_zmset(i)
INTEGER,INTENT(IN)::i
RETURN
END SUBROUTINE Helac_zmset

FUNCTION eps(i)
IMPLICIT NONE
INTEGER,INTENT(IN)::i
REAL(KIND=DBL)::eps
!eps=10d0*epsilon(dnou(i))
eps=1.0d-14
END FUNCTION eps

FUNCTION ALPHAS(sc)
USE alfas_functions
!   Refer to hep-ph/9701390
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::sc
REAL(KIND=DBL)::ALPHAS
INTEGER::NF,K,I
REAL(KIND=DBL)::Q2
REAL(KIND=DBL)::pi
!...HEAVY QUARK LAMBDA VALUES
REAL(KIND=DBL),DIMENSION(3:6)::LAMBDAL=(/ 0.2470d0, 0.2150d0, 0.1650d0, 0.0847d0 /),&
                               LAMBDAN=(/ 0.3500d0, 0.3260d0, 0.2262d0, 0.1510d0 /)
!LAMBDAL=(/ 0.2470, 0.2150, 0.1650, 0.0847 /),&
!LAMBDAN=(/ 0.3500, 0.3260, 0.2260, 0.1510 /)
!LAMBDAL=(/ 0.2041, 0.1750, 0.1320, 0.0665 /),&
!LAMBDAN=(/ 0.2994, 0.2460, 0.1677, 0.0678 /)
!...HEAVY QUARK THRESHOLDS m_c=1.4, m_b=4.5
REAL(KIND=DBL),DIMENSION(3)::Q2THR=(/1.690d0,20.25d0,32400d0 /)
REAL(KIND=DBL)::B0,ALP,LAM2,B1,B10,LQ2,XLP,XL,XLM,Y,Y1
INCLUDE "../input/lhapdf/call_alphas"
IF(useMCFMrun)THEN
   ALPHAS=ALPHAS_MCFM(sc)
   RETURN
ENDIF
!...DETERMINATION OF THE APPROPRIATE NUMBER OF FLAVOURS :
Q2=sc**2
NF = 3
DO K = 1, 3
   IF (Q2 .GT. Q2THR (K)) THEN
      NF = NF + 1
   ELSE
      EXIT
   END IF
ENDDO
!...LO ALPHA_S :
B0 = 11.d0- 2.d0/3.d0* NF
IF(iPDFSUP1.NE.10000.AND.iPDFSUP1.NE.10041)THEN
	LAM2 = LAMBDAL (NF) * LAMBDAL (NF)
	ALP  = 1.d0/(B0 * DLOG (Q2/LAM2))
!...NLO ALPHA_S :
ELSE                     ! NLO
    LAM2 = LAMBDAN (NF) * LAMBDAN (NF)
    B1 = 102.- 38./3.* NF
    B10 = B1 / (B0*B0)
	!...START VALUE FOR NLO ITERATION :
    LQ2 = DLOG (Q2 / LAM2)
    ALP = 1./(B0*LQ2) * (1.- B10*DLOG(LQ2)/LQ2)

    !...EXACT NLO VALUE, FOUND VIA NEWTON PROCEDURE : 
	! Find root of f(x)=0, use iteration of x_n=x_(n-1)-f(x_(n-1))/f'(x_(n-1))
	! Here, f(ALP)=Y
    DO I = 1, 6
      XL  = DLOG (1./(B0*ALP) + B10)
      XLP = DLOG (1./(B0*ALP*1.01) + B10)
      XLM = DLOG (1./(B0*ALP*0.99) + B10)
      Y  = LQ2 - 1./ (B0*ALP) + B10 * XL
      Y1 = (- 1./ (B0*ALP*1.01) + B10 * XLP &
        + 1./ (B0*ALP*0.99) - B10 * XLM) / (0.02d0*ALP)  ! Modification:XLP-> XLM
      ALP = ALP - Y/Y1
    ENDDO
ENDIF

CALL Helac_mypi(pi)
!...OUTPUT :
ALPHAS = ALP*4d0*pi

END FUNCTION ALPHAS

 SUBROUTINE ALPSOR(A,N1,K,IOPT)
!-----------------------------------------------------------------------
!     Sort A(N) into ascending order
!     IOPT = 1 : return sorted A and index array K
!     IOPT = 2 : return index array K only
!-----------------------------------------------------------------------
!     DOUBLE PRECISION A(N),B(5000)
!     INTEGER N,I,J,IOPT,K(N),K1(5000),IL(5000),IR(5000)
INTEGER,INTENT(IN)::N1
REAL(KIND(1d0)),DIMENSION(N1),INTENT(INOUT)::A
REAL(KIND(1d0)),DIMENSION(5000)::B
INTEGER,DIMENSION(N1),INTENT(OUT)::K
!INTEGER,DIMENSION(N1)::K1
INTEGER,INTENT(IN)::IOPT
INTEGER::I,J,flag10,flag8,flag81,flag82
INTEGER,DIMENSION(5000)::IL,IR
IF (N1.GT.5000) then
   WRITE(*,*) 'Too many entries to sort in alpsrt, stop'
   STOP
ENDIF
IF(N1.LE.0)RETURN
IL(1)=0
IR(1)=0
DO I=2,N1
   IL(I)=0
   IR(I)=0
   J=1
   DO
    DO  ! find the smallest index j that A(j)<A(I), 
	    ! otherwise find the largest index j that A(j)>=A(I). 
      IF(A(I).GT.A(J)) EXIT   ! goto 5      ! 2
	  ! don't search the index smaller than j
      IF(IL(J).EQ.0) EXIT     ! goto 4      ! 3
      J=IL(J)
    ENDDO                      !GOTO 2
	flag10=0
	! the largest index j that A(j)>=A(I).
	IF(A(I).LE.A(J))THEN
      IR(I)=-J                !  4
      IL(J)=I
	  flag10=1
      EXIT                     ! goto 10
	ENDIF
	! don't search the index smaller than j ! 5
    IF(IR(J).LE.0) EXIT       ! goto 6
    J=IR(J)
   ENDDO                             !GOTO 2
   IF(flag10.EQ.1)THEN
     CYCLE
   ELSE
     IR(I)=IR(J)                 ! 6
     IR(J)=I
   ENDIF
ENDDO
!  10  CONTINUE
I=1
J=1
flag8=0
flag81=0
!      GOTO 8
DO
  DO
    IF(flag8.EQ.1.AND.flag81.EQ.0)THEN
       J=IL(J)                ! 20
    ENDIF
	! when the element A(k)<A(j) (at the right hand of j) and k>j
    IF(IL(J).GT.0)THEN        ! GOTO 20 ! 8
       flag8=1
	   flag81=0
	ELSE
       EXIT
    ENDIF
  ENDDO

  DO
     flag81=0
     flag82=0
     K(I)=J               !9
     B(I)=A(J)
     I=I+1
     IF(IR(J).GT.0)THEN               !12,30,13
        J=IR(J)                        ! 13
        flag81=1
        EXIT                           ! goto 8
	 ! this index is must the nearest to J
     ELSEIF(IR(J).LT.0)THEN
       J=-IR(J)                       ! 12
                                 !   GOTO 9
     ELSE
       flag82=1
       EXIT
     ENDIF
  ENDDO
                                 !  30  CONTINUE
  IF(flag82.EQ.1)EXIT
ENDDO
!DO I=1,N1
!   K1(I)=K(N1+1-I)
!ENDDO
!DO I=1,N1
!   K(I)=K1(I)
!ENDDO
IF(IOPT.EQ.2) RETURN
DO I=1,N1          ! 31
  A(I)=B(I)
ENDDO
END SUBROUTINE ALPSOR

FUNCTION QN1S0F(kk)
IMPLICIT NONE
INTEGER,INTENT(IN)::kk
LOGICAL::QN1S0F
QN1S0F=.FALSE.
IF(kk.GT.nhad.OR.kk.LE.0)RETURN
IF(SubInteger(iflh(kk),3,4).NE.10)RETURN
QN1S0F=.TRUE.
!QN1S0F=(iflh(kk).EQ.441001).OR.(iflh(kk).EQ.441008).OR.&
!(iflh(kk).EQ.551001).OR.(iflh(kk).EQ.551008).OR.(iflh(kk).EQ.451001)&
!.OR.(iflh(kk).EQ.451008)
RETURN
END FUNCTION QN1S0F

FUNCTION QN3S1F(kk)
IMPLICIT NONE
INTEGER,INTENT(IN)::kk
LOGICAL::QN3S1F
QN3S1F=.FALSE.
IF(kk.GT.nhad.OR.kk.LE.0)RETURN
IF(SubInteger(iflh(kk),3,4).NE.30)RETURN
QN3S1F=.TRUE.
!QN3S1F=(iflh(kk).EQ.443011).OR.(iflh(kk).EQ.443018).OR.&
!(iflh(kk).EQ.553011).OR.(iflh(kk).EQ.553018).OR.(iflh(kk).EQ.453011)&
!.OR.(iflh(kk).EQ.453018)
RETURN
END FUNCTION QN3S1F

FUNCTION QN1P1F(kk)
IMPLICIT NONE
INTEGER,INTENT(IN)::kk
LOGICAL::QN1P1F
QN1P1F=.FALSE.
IF(kk.GT.nhad.OR.kk.LE.0)RETURN
IF(SubInteger(iflh(kk),3,4).NE.11)RETURN
QN1P1F=.TRUE.
!QN1P1F=(iflh(kk).EQ.441111).OR.(iflh(kk).EQ.441118).OR.&
!(iflh(kk).EQ.551111).OR.(iflh(kk).EQ.551118).OR.(iflh(kk).EQ.451111)&
!.OR.(iflh(kk).EQ.451118)
RETURN
END FUNCTION QN1P1F

FUNCTION QN3PJF(kk)
IMPLICIT NONE
INTEGER,INTENT(IN)::kk
LOGICAL::QN3PJF
QN3PJF=.FALSE.
IF(kk.GT.nhad.OR.kk.LE.0)RETURN
IF(SubInteger(iflh(kk),3,4).NE.31)RETURN
QN3PJF=.TRUE.
RETURN
END FUNCTION QN3PJF

FUNCTION SubInteger(int,n1,n2)
IMPLICIT NONE
INTEGER,INTENT(IN)::int,n1,n2
INTEGER::nmax,nmin,mod11,int2
INTEGER::SubInteger
int2=ABS(int)
nmax=MAX(n1,n2)
nmin=MIN(n1,n2)
IF(int2.LT.10**(nmin-1))THEN
!	PRINT *,"Warning in SubInteger !"
	SubInteger=0
	RETURN
ELSEIF(nmin.GT.0)THEN
	mod11=MOD(int2,10**nmax)
	SubInteger=mod11/10**(nmin-1)
	RETURN
ELSE
	PRINT *,"Error in SubInteger ! STOP!"
	STOP
ENDIF
END FUNCTION SubInteger

! for every number, it generates the helicity for (U_lamda,eps(L_z),eps(S_z))
!*********************************************************
!     ipol(k)=helicity,ipol(k)=2:rigth-handed            *
!                      ipol(k)=1:left-handed             *
!                      ipol(k)=3:long.                   *
!*********************************************************
SUBROUTINE Qeveryhel(ihela,num,kk)  !num=ipol(kk)
IMPLICIT NONE
INTEGER,DIMENSION(3),INTENT(OUT)::ihela
INTEGER,INTENT(IN)::num,kk
INTEGER::ii1,jj1,kk1
IF(kk.GT.nhad.OR.kk.LE.2)THEN
	PRINT *,"Error 1 in Qeveryhel ! STOP !"
	STOP
ENDIF
IF(ABS(iflh(kk)).LE.100)THEN
	PRINT *,"Error 2 in Qeveryhel ! STOP !"
	STOP
ENDIF
IF(iranhel.EQ.3)THEN
	ihela(1)=num
	ihela(2:3)=0
	RETURN
ENDIF
IF(QN1S0F(kk))THEN
	ihela(1)=num
	ihela(2:3)=0
	RETURN
ENDIF
IF(QN3S1F(kk))THEN
	ihela(2)=0
	ii1=MOD(num-1,2)
	jj1=(num-1)/2
	ihela(1)=ii1+1
	ihela(3)=jj1+1
	RETURN
ENDIF
IF(QN1P1F(kk))THEN
	IF(iranhel.EQ.2)THEN
		ihela(1)=num
		ihela(2:3)=0
		RETURN
	ELSE
		ihela(3)=0
		ii1=MOD(num-1,2)
		jj1=(num-1)/2
		ihela(1)=ii1+1
		ihela(2)=jj1+1
		RETURN
	ENDIF
ENDIF
IF(QN3PJF(kk))THEN
	IF(iranhel.EQ.2)THEN
		ihela(2)=0
		ii1=MOD(num-1,2)
		jj1=(num-1)/2
		ihela(1)=ii1+1
		ihela(3)=jj1+1
		RETURN
	ELSE
		ii1=MOD(num-1,2)
		jj1=MOD((num-1)/2,3)
		kk1=(num-1)/6
		ihela(1)=ii1+1
		ihela(2)=jj1+1
		ihela(3)=kk1+1
		RETURN
	ENDIF
ENDIF
ihela(1:3)=0
RETURN
END SUBROUTINE Qeveryhel

FUNCTION ExterExistQ(n1,x,k1)
IMPLICIT NONE
INTEGER,INTENT(IN)::n1,k1,x
INTEGER,DIMENSION(n1)::y
LOGICAL::ExterExistQ
ExterExistQ=.FALSE.
IF(k1.LE.0.OR.k1.GT.n1)RETURN
CALL Helac_bin(n1,x,y)
IF(y(n1+1-k1).EQ.1)ExterExistQ=.TRUE.
RETURN
END FUNCTION ExterExistQ
!zq(1:n,1)=(p0+pz,pz);zq(1:n,2)=(p0-pz,pt);zq(1:n,3)=(px,py);
!zq(1:n,4)=(px,-py);zq(1:n,5)=(m,p)
! returns the sum zq0 of two momentum zq1 and zq2
SUBROUTINE PlusZq(zq0,zq1,zq2)
IMPLICIT NONE
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(IN)::zq1,zq2
COMPLEX(KIND=DBL),DIMENSION(5),INTENT(OUT)::zq0
REAL(KIND=DBL)::rqrq,rpt,rm,rp,rpz,rp0
zq0(1)=zq1(1)+zq2(1)
zq0(3:4)=zq1(3:4)+zq2(3:4)
rpt=CDABS(zq0(3))
rpz=DIMAG(zq0(1))
rp0=DREAL(zq0(1))-rpz
rp=DSQRT(rpt**2+rpz**2)
rm=DSQRT(rp0**2-rp**2)
zq0(2)=DCMPLX(rp0-rpz,rpt)
zq0(5)=DCMPLX(rm,rp)
END SUBROUTINE PlusZq

SUBROUTINE MomInclude(ii,iparray)
IMPLICIT NONE
INTEGER,INTENT(IN)::ii
INTEGER,DIMENSION(10),INTENT(OUT)::iparray
INTEGER,DIMENSION(n)::iiy
INTEGER::k,k1,k2,k3
CALL Helac_bin(n,ii,iiy)
iparray(1:10)=0
DO k=1,Pwavenum
	k1=Pwave(k)
	k1=Quarkonium2(k1)
	k1=Quarkonium(k1,2)
	k2=iiy(n+1-k1)*io(k1)
	k3=-iiy(n-k1)*io(k1+1)
	iparray(k)=k2+k3
ENDDO
END SUBROUTINE MomInclude

SUBROUTINE Gethelnum(num)
IMPLICIT NONE
INTEGER,INTENT(OUT)::num
INTEGER::k,l1,l2,l3,l9
INTEGER,DIMENSION(nhad)::lp
l1=0
l2=0
l3=0
l9=0
IF(iranhel.EQ.3)THEN
	num=1
	RETURN
ENDIF
DO k=1,nhad
	lp(k)=1
	IF(iranhel.EQ.0)THEN
		IF(iflh(k).GE.-12.AND.iflh(k).LE.12)lp(k)=2
		IF(iflh(k).EQ.31.OR.iflh(k).EQ.35)lp(k)=2            !change 5.4.00
		IF(iflh(k).GE.32.AND.iflh(k).LE.34)lp(k)=3
	ENDIF
	IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
		IF(QN1S0F(k))lp(k)=1
		IF(QN3S1F(k).OR.QN1P1F(k))lp(k)=3
		IF(QN3PJF(k))lp(k)=9
	ELSEIF(iranhel.EQ.2)THEN
		IF(QN1S0F(k).OR.QN1P1F(k))lp(k)=1
		IF(QN3S1F(k).OR.QN3PJF(k))lp(k)=3
	ENDIF
ENDDO
DO k=1,nhad
	IF(iranhel.EQ.0)THEN
		IF(lp(k).EQ.1.AND.(ABS(iflh(k)).LE.100))l1=l1+1
		IF(lp(k).EQ.2.AND.(ABS(iflh(k)).LE.100))l2=l2+1
		IF(lp(k).EQ.3.AND.(ABS(iflh(k)).LE.100))l3=l3+1
	ENDIF
	IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
		IF(QN1S0F(k))l1=l1+1
		IF(QN3S1F(k))l3=l3+1
		IF(QN1P1F(k))l3=l3+1
		IF(QN3PJF(k))l9=l9+1
	ELSEIF(iranhel.EQ.2)THEN
		IF(QN1S0F(k).OR.QN1P1F(k))l1=l1+1
		IF(QN3S1F(k).OR.QN3PJF(k))l3=l3+1
	ENDIF
ENDDO
num=2**(l2)*3**(l3)*9**(l9)
RETURN
END SUBROUTINE Gethelnum
! return the index for the physical polarizations
SUBROUTINE PhysicalPol(num)
IMPLICIT NONE
INTEGER,INTENT(OUT)::num
INTEGER::i,k,ltot,jz,init=0
SAVE init
INTEGER,DIMENSION(20)::lp
SAVE lp
IF(iranhel.EQ.3)THEN
	num=1
	RETURN
ENDIF
IF(GLOBALINIT_physicalpol.EQ.0)THEN
	init=0
	lp(1:20)=1
ENDIF
IF(init.EQ.0)THEN
	DO k=1,nhad
		lp(k)=1
		IF(iranhel.EQ.0)THEN
			IF(iflh(k).GE.-12.AND.iflh(k).LE.12)lp(k)=2
			IF(iflh(k).EQ.31.OR.iflh(k).EQ.35)lp(k)=2            !change 5.4.00
			IF(iflh(k).GE.32.AND.iflh(k).LE.34)lp(k)=3
		ENDIF
		IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
			IF(QN1S0F(k))lp(k)=2
			IF(QN3S1F(k).OR.QN1P1F(k))lp(k)=6
			IF(QN3PJF(k))lp(k)=18
		ELSEIF(iranhel.EQ.2)THEN
			IF(QN1S0F(k).OR.QN1P1F(k))lp(k)=2
			IF(QN3S1F(k).OR.QN3PJF(k))lp(k)=6
		ENDIF
	ENDDO
	init=1
ENDIF
ltot=1
jz=1
DO i=1,nhad
	IF(ABS(iflh(i)).GT.100)THEN
		ltot=ltot+jz*((ipol(i)-1)/2)
		jz=((lp(i)-1)/2+1)*jz
	ELSE
		ltot=ltot+jz*(ipol(i)-1)
		jz=lp(i)*jz
	ENDIF
ENDDO
num=ltot
END SUBROUTINE PhysicalPol

SUBROUTINE PhysicalPol2(num,factor)
IMPLICIT NONE
INTEGER,INTENT(OUT)::num
COMPLEX(KIND=DBL),INTENT(OUT)::factor
INTEGER::i,k,ltot,jz,init=0,sz,lz
INTEGER,DIMENSION(3)::index2
SAVE init
INTEGER,DIMENSION(20)::lp
SAVE lp
IF(GLOBALINIT_physicalpol2.EQ.0)THEN
	init=0
	lp(1:20)=1
ENDIF
FirstHelicityQ=.TRUE.
factor=DCMPLX(1d0)
IF(init.EQ.0)THEN
	DO k=1,nhad
		lp(k)=1
		IF(iranhel.EQ.0)THEN
			IF(iflh(k).GE.-12.AND.iflh(k).LE.12)lp(k)=2
			IF(iflh(k).EQ.31.OR.iflh(k).EQ.35)lp(k)=2            !change 5.4.00
			IF(iflh(k).GE.32.AND.iflh(k).LE.34)lp(k)=3
		ENDIF
		IF(iranhel.EQ.1.OR.iranhel.EQ.0)THEN
			IF(QN1S0F(k))lp(k)=2
			IF(QN3S1F(k).OR.QN1P1F(k))lp(k)=6
			IF(QN3PJF(k))lp(k)=18
		ELSEIF(iranhel.EQ.2)THEN
			IF(QN1S0F(k).OR.QN1P1F(k))lp(k)=2
			IF(QN3S1F(k).OR.QN3PJF(k))lp(k)=6
		ENDIF
	ENDDO
	init=1
ENDIF
ltot=1
jz=1
DO i=1,nhad
	IF(ABS(iflh(i)).GT.100)THEN
		IF(QN3S1F(i).AND.Quarkonium(1,1).EQ.i)THEN
			CALL Qeveryhel(index2,ipol(i),i)
			sz=transfpol(index2(3))
			IF(ihel1.EQ.sz)THEN
				factor=1
			ELSEIF(ihel1.NE.ihel2.AND.ihel2.EQ.sz)THEN
				factor=1
				FirstHelicityQ=.FALSE.
			ELSE
				factor=0
			ENDIF
			ltot=ltot
			jz=((lp(i)-1)/2+1)*jz
		ELSEIF(QN1P1F(i).AND.i.EQ.Quarkonium(1,1))THEN
			CALL Qeveryhel(index2,ipol(i),i)
			lz=transfpol(index2(3))
			IF(ihel1.EQ.lz)THEN
				factor=1
			ELSEIF(ihel1.NE.ihel2.AND.ihel2.EQ.lz)THEN
				factor=1
				FirstHelicityQ=.FALSE.
			ELSE
				factor=0
			ENDIF
			ltot=ltot
			jz=((lp(i)-1)/2+1)*jz
		ELSEIF(QN3PJF(i).AND.i.EQ.Quarkonium(1,1))THEN
			CALL Qeveryhel(index2,ipol(i),i)
			IF(iLSJ.EQ.0)THEN
				sz=transfpol(index2(2))
				lz=index2(3)-1
			ELSEIF(iLSJ.EQ.1)THEN
				sz=transfpol(index2(3))
				IF(iranhel.EQ.2)THEN
					lz=0
				ELSE
					lz=index2(2)-1
				ENDIF
			ENDIF
			IF(ihel1.EQ.sz)THEN
				factor=1
			ELSEIF(ihel1.NE.ihel2.AND.ihel2.EQ.sz)THEN
				factor=1
				FirstHelicityQ=.FALSE.
			ELSE
				factor=0
			ENDIF
			ltot=ltot+jz*lz
			jz=((lp(i)-1)/2+1)*jz
		ELSE
			ltot=ltot+jz*((ipol(i)-1)/2)
			jz=((lp(i)-1)/2+1)*jz
		ENDIF
	ELSE
		ltot=ltot+jz*(ipol(i)-1)
		jz=lp(i)*jz
	ENDIF
ENDDO
num=ltot
END SUBROUTINE PhysicalPol2

SUBROUTINE ls2j(num,factor)
IMPLICIT NONE
INTEGER,INTENT(OUT)::num
COMPLEX(KIND=DBL),INTENT(OUT)::factor
INTEGER::i,k,ltot,jz,init=0
!INTEGER,DIMENSION(:),ALLOCATABLE::Xprime
!COMPLEX(KIND=DBL),DIMENSION(:),ALLOCATABLE::CGarray
! at most 2 3pj
INTEGER,PARAMETER::len=80
INTEGER,DIMENSION(0:len)::Xprime
COMPLEX(KIND=DBL),DIMENSION(0:len)::CGarray
LOGICAL,DIMENSION(0:len)::HArrayQ
SAVE init,Xprime,CGarray,HArrayQ
INTEGER,DIMENSION(20)::lp
INTEGER,DIMENSION(3)::index2
SAVE lp
INTEGER::ISL2JX,ikik,lz,sz
!INTEGER::istat
IF(GLOBALINIT_ls2j.EQ.0)THEN
	init=0
	lp(1:20)=1
	Xprime(0:80)=-1
	CGarray(0:80)=0
	HArrayQ(0:80)=.TRUE.
!	IF(ALLOCATED(Xprime))DEALLOCATE(Xprime)
!	IF(ALLOCATED(CGarray))DEALLOCATE(CGarray)
ENDIF
IF(init.EQ.0)THEN
!	len=9**Num3pj-1
!	ALLOCATE(Xprime(len),STAT=istat)
!	ALLOCATE(CGarray(len),STAT=istat)
!	IF(istat.NE.0)THEN
!		WRITE(*,*)istat, ' warning: allocation is not working properly in ls2j in Helac_Func_1.f90 !'
!		STOP
!	ENDIF
	DO k=1,nhad
		lp(k)=1
		IF(iranhel.EQ.0)THEN
			IF(iflh(k).GE.-12.AND.iflh(k).LE.12)lp(k)=2
			IF(iflh(k).EQ.31.OR.iflh(k).EQ.35)lp(k)=2            !change 5.4.00
			IF(iflh(k).GE.32.AND.iflh(k).LE.34)lp(k)=3
		ENDIF
		IF(QN1S0F(k))lp(k)=2
		IF(QN3S1F(k).OR.QN1P1F(k))lp(k)=6
		IF(QN3PJF(k))lp(k)=18
	ENDDO
	init=1
ENDIF
factor=1
FirstHelicityQ=.TRUE.
ISL2JX=0
DO ikik=1,Num3pj
	ISL2JX=ISL2JX+(ipol(Pwave3pj(ikik))-1)/2*9**(ikik-1)
ENDDO
IF(ISL2JX.GT.len)THEN
	WRITE(*,*)"Warning : X > len in ls2j in Helac_Func_1.f90"
	STOP
ENDIF
IF(Num3pj.EQ.0)THEN
	IF(imode.EQ.0)THEN
		CALL PhysicalPol(num)
		GLOBALINIT_physicalpol=1
		factor=1
	ELSE
		CALL PhysicalPol2(num,factor)
		GLOBALINIT_physicalpol2=1
	ENDIF
	RETURN
ELSE
	IF(Xprime(ISL2JX).EQ.-1)THEN
		factor=1
		num=1
		DO ikik=1,Num3pj
			CALL Qeveryhel(index2,ipol(Pwave3pj(ikik)),Pwave3pj(ikik))
			lz=transfpol(index2(2))
			sz=transfpol(index2(3))
			IF(imode.EQ.0.OR.(imode.EQ.1.AND.Pwave3pj(ikik).NE.Quarkonium(1,1)))THEN
				factor=factor*ClebschGordan2((/1,lz/),(/1,sz/),(/JNum3pj(ikik),lz+sz/))
				num=num+(JNum3pj(ikik)**2+MOD(lz+sz+JNum3pj(ikik),2*JNum3pj(ikik)+1))*9**(ikik-1)
			ELSEIF(imode.EQ.1.AND.Pwave3pj(ikik).EQ.Quarkonium(1,1))THEN
				IF(lz+sz.EQ.ihel1)THEN
					factor=factor*ClebschGordan2((/1,lz/),(/1,sz/),(/JNum3pj(ikik),lz+sz/))
				ELSEIF(lz+sz.EQ.ihel2.AND.ihel1.NE.ihel2)THEN
					factor=factor*ClebschGordan2((/1,lz/),(/1,sz/),(/JNum3pj(ikik),lz+sz/))
					FirstHelicityQ=.FALSE.
				ELSE
					factor=0
				ENDIF
				num=num+JNum3pj(ikik)**2*9**(ikik-1)
			ENDIF
		ENDDO
		Xprime(ISL2JX)=num
		CGarray(ISL2JX)=factor
		HarrayQ(ISL2JX)=FirstHelicityQ
	ELSE
		num=Xprime(ISL2JX)
		factor=CGarray(ISL2JX)
		FirstHelicityQ=HarrayQ(ISL2JX)
	ENDIF
ENDIF
ltot=1
jz=1
DO i=1,nhad
	IF(Pwave3pjQ(i).EQ.0)THEN
		IF(ABS(iflh(i)).GT.100)THEN
			IF(i.EQ.Quarkonium(1,1).AND.QN3S1F(i).AND.imode.EQ.1)THEN
				CALL Qeveryhel(index2,ipol(i),i)
				sz=transfpol(index2(3))
				IF(ihel1.EQ.sz)THEN
					factor=factor
				ELSEIF(ihel1.NE.ihel2.AND.ihel2.EQ.sz)THEN
					factor=factor
					FirstHelicityQ=.FALSE.
				ELSE
					factor=0
				ENDIF
				ltot=ltot
				jz=((lp(i)-1)/2+1)*jz
			ELSEIF(i.EQ.Quarkonium(1,1).AND.QN1P1F(i).AND.imode.EQ.1)THEN
				CALL Qeveryhel(index2,ipol(i),i)
				lz=transfpol(index2(2))
				IF(ihel1.EQ.lz)THEN
					factor=factor
				ELSEIF(ihel1.NE.ihel2.AND.ihel2.EQ.lz)THEN
					factor=factor
					FirstHelicityQ=.FALSE.
				ELSE
					factor=0
				ENDIF
				ltot=ltot
				jz=((lp(i)-1)/2+1)*jz
			ELSE
				ltot=ltot+jz*((ipol(i)-1)/2)
				jz=((lp(i)-1)/2+1)*jz
			ENDIF
		ELSE
			ltot=ltot+jz*(ipol(i)-1)
			jz=lp(i)*jz
		ENDIF
	ENDIF
ENDDO
num=num+(9**Num3pj)*(ltot-1)
END SUBROUTINE ls2j

FUNCTION transfpol(num)
IMPLICIT NONE
INTEGER,INTENT(IN)::num
INTEGER::transfpol
IF(num.EQ.3)THEN
	transfpol=0
ELSEIF(num.EQ.2)THEN
	transfpol=1
ELSEIF(num.EQ.1)THEN
	transfpol=-1
ELSE
	PRINT *,"Wrong in transfpol of Helac_Func_1.f90 ! STOP !"
	STOP
ENDIF
END FUNCTION transfpol

SUBROUTINE SetRandom(lr,sr)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(*),INTENT(IN)::lr,sr
INTEGER::i,kk1,kk2
Qhelran(1:20)=-1
QSzran(1:20)=-1
IF(iranhel.EQ.0)RETURN
kk1=1
kk2=1
DO i=1,nhad
	IF(ABS(iflh(i)).LE.40.AND.iranhel.NE.0)THEN
		Qhelran(i)=lr(kk1)
		kk1=kk1+1
	ENDIF
	IF(ABS(iflh(i)).GT.100.AND.SubInteger(iflh(i),3,3).EQ.1.AND.iranhel.GE.2)THEN
		Qhelran(i)=lr(kk1)
		kk1=kk1+1
	ENDIF
	IF(ABS(iflh(i)).GT.100.AND.SubInteger(iflh(i),4,4).EQ.3.AND.iranhel.EQ.3)THEN
		QSzran(i)=sr(kk2)
		kk2=kk2+1
	ENDIF
ENDDO
END SUBROUTINE SetRandom

FUNCTION ClebschGordan2(LL,SS,JJ)
IMPLICIT NONE
INTEGER,DIMENSION(2),INTENT(IN)::LL,SS,JJ
COMPLEX(KIND(1d0))::ClebschGordan2
IF((LL(2)+SS(2).NE.JJ(2)).OR.(LL(1).LT.0).OR.(SS(1).LT.0).OR.(JJ(1).LT.0) &
.OR.(ABS(LL(2)).GT.LL(1)).OR.(ABS(SS(2)).GT.SS(1)).OR.(ABS(JJ(2)).GT.JJ(1))&
.OR.(LL(1)+SS(1).LT.JJ(1)).OR.(ABS(LL(1)-SS(1)).GT.JJ(1)))THEN
	ClebschGordan2=DCMPLX(0)
	RETURN
END IF
IF(LL(1).EQ.0)THEN
	IF((SS(1).NE.JJ(1)).OR.(SS(2).NE.JJ(2)))THEN
		ClebschGordan2=DCMPLX(0)
		RETURN
	ELSE
		ClebschGordan2=DCMPLX(1)
		RETURN
	ENDIF
ENDIF
IF(SS(1).EQ.0)THEN
	IF((LL(1).NE.JJ(1)).OR.(LL(2).NE.JJ(2)))THEN
		ClebschGordan2=DCMPLX(0)
		RETURN
	ELSE
		ClebschGordan2=DCMPLX(1)
		RETURN
	ENDIF
ENDIF
IF((LL(1).EQ.1).AND.(SS(1).EQ.1))THEN
	SELECT CASE(JJ(1))
	CASE(1)
		SELECT CASE(LL(2))
		CASE(0)
			ClebschGordan2=DCMPLX(-SS(2)*1d0/DSQRT(2d0))
			RETURN
		CASE(1,-1)
			ClebschGordan2=DCMPLX(DBLE(LL(2))*1d0/DSQRT(2d0))
			RETURN
		END SELECT
	CASE(2)
		SELECT CASE(LL(2))
		CASE(1,-1)
			IF(SS(2).EQ.LL(2))THEN
				ClebschGordan2=DCMPLX(1d0)
				RETURN
			ENDIF
			IF(SS(2).EQ.0)THEN
				ClebschGordan2=DCMPLX(1d0/DSQRT(2d0))
				RETURN
			ENDIF
			IF(SS(2)+LL(2).EQ.0)THEN
				ClebschGordan2=DCMPLX(1d0/DSQRT(6d0))
				RETURN
			ENDIF
		CASE(0)
			SELECT CASE(SS(2))
			CASE(1,-1)
				ClebschGordan2=DCMPLX(1d0/DSQRT(2d0))
				RETURN
			CASE(0)
				ClebschGordan2=DCMPLX(DSQRT(2d0/3d0))
				RETURN
			END SELECT
		END SELECT
	CASE(0)
		SELECT CASE(LL(2))
		CASE(1,-1)
			ClebschGordan2=DCMPLX(1d0/DSQRT(3d0))
			RETURN
		CASE(0)
			ClebschGordan2=DCMPLX(-1d0/DSQRT(3d0))
			RETURN
		END SELECT
	END SELECT
ENDIF
PRINT *,"Wrong in ClebschGordan2 ! STOP !"
STOP
END FUNCTION ClebschGordan2

SUBROUTINE DAYTIME(iy,im,id,ih,imin,isec,imilsec)
INTEGER,INTENT(OUT)::iy,im,id,ih,imin,isec,imilsec
INTEGER,DIMENSION(8)::timevalues
! GETDAY and GETTIM can only be used in ifort (not gfortran)
CALL DATE_AND_TIME(VALUES=timevalues)
iy=timevalues(1)
im=timevalues(2)
id=timevalues(3)
ih=timevalues(5)
imin=timevalues(6)
isec=timevalues(7)
imilsec=timevalues(8)
END SUBROUTINE DAYTIME

FUNCTION lepton_pdg(ipdg)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::ipdg
  LOGICAL::lepton_pdg
  lepton_pdg=.FALSE.
  IF(ABS(ipdg).EQ.11.OR.ABS(ipdg).EQ.13.OR.ABS(ipdg).EQ.15)lepton_pdg=.TRUE.
  RETURN
END FUNCTION lepton_pdg

FUNCTION lepton_ho(ifif)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::ifif
  LOGICAL::lepton_ho
  lepton_ho=.FALSE.
  IF(ABS(ifif).EQ.2.OR.ABS(ifif).EQ.6.OR.ABS(ifif).EQ.10)lepton_ho=.TRUE.
  RETURN
END FUNCTION lepton_ho

FUNCTION neutrino_pdg(ipdg)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::ipdg
  LOGICAL::neutrino_pdg
  neutrino_pdg=.FALSE.
  IF(ABS(ipdg).EQ.12.OR.ABS(ipdg).EQ.14.OR.ABS(ipdg).EQ.16)neutrino_pdg=.TRUE.
  RETURN
END FUNCTION neutrino_pdg

FUNCTION neutrino_ho(ifif)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::ifif
  LOGICAL::neutrino_ho
  neutrino_ho=.FALSE.
  IF(ABS(ifif).EQ.1.OR.ABS(ifif).EQ.5.OR.ABS(ifif).EQ.9)neutrino_ho=.TRUE.
  RETURN
END FUNCTION neutrino_ho

FUNCTION pdgt(ifif)
! from HELAC-Onia notation to PDG notation
INTEGER,INTENT(IN)::ifif
INTEGER::pdgt,init=0,ifac
INTEGER,DIMENSION(-100:100)::ipdt
SAVE ipdt,init,ifac
IF(init.EQ.0)THEN
!                include 'pdt.h'  to the representation in PDG  
   ipdt( 1)= 12
   ipdt(-1)=-12
   ipdt( 5)= 14
   ipdt(-5)=-14
   ipdt( 9)= 16
   ipdt(-9)=-16
   ipdt( 2)= 11
   ipdt(-2)=-11
   ipdt( 6)= 13
   ipdt(-6)=-13
   ipdt( 10)= 15
   ipdt(-10)=-15
   ipdt( 3)= 2
   ipdt(-3)=-2
   ipdt( 4)= 1
   ipdt(-4)=-1
   ipdt( 7)= 4
   ipdt(-7)=-4
   ipdt( 8)= 3
   ipdt(-8)=-3
   ipdt( 11)= 6
   ipdt(-11)=-6
   ipdt( 12)= 5
   ipdt(-12)=-5
   ipdt(35)=21
   ipdt(31)=22
   ipdt(32)=23
   ipdt(33)=24
   ipdt(34)=-24
   ipdt(41)=25
   ! end of 'pdt.h'   
   init=1
ENDIF
IF(ABS(ifif).LT.100)THEN
   pdgt=ipdt(ifif)
   RETURN
ENDIF
ifac=1
IF(ifif.LT.0)ifac=-1
SELECT CASE(ABS(ifif))
   ! J/psi   
CASE(443011)
   pdgt=443
   ! etac  
CASE(441001)
   pdgt=441
   ! hc 
CASE(441111)
   pdgt=10443
   ! chic0 
CASE(443101)
   pdgt=10441
   ! chic1 
CASE(443111)
   pdgt=20443
   ! chic2
CASE(443121)
   pdgt=445
   ! upsilon 
CASE(553011)
   pdgt=553
   ! etab 
CASE(551001)
   pdgt=551
   ! hb 
CASE(551111)
   pdgt=10553
   ! chib0 
CASE(553101)
   pdgt=10551
   ! chib1 
CASE(553111)
   pdgt=20553
   ! chib2
CASE(553121)
   pdgt=555
   ! Bc*+-
CASE(453011)
   pdgt=543*ifac
   ! Bc+- 
CASE(451001)
   pdgt=541*ifac
   ! Bc2*+- 
CASE(453121)
   pdgt=545*ifac
   ! Bc0*+-
CASE(453101)
   pdgt=10541*ifac
   ! Bc1(H)+-
CASE(453111)
   pdgt=20543*ifac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!                                                          !  
!     To matching the notations in PYTHIA 8                ! 
!                                                          !  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   ! ccbar(3P0(8))
CASE(443108)
   pdgt=9910441
   ! bbbar(3P0(8))
CASE(553108)
   pdgt=9910551
   ! bbbar(3S1(8)) 
CASE(553018)
   pdgt=9900553
   ! bbbar(1S0(8))
CASE(551008)
   pdgt=9900551
   ! ccbar(3S1(8)) 
CASE(443018)
   pdgt=9900443
   ! ccbar(1S0(8))
CASE(441008)
   pdgt=9900441
CASE DEFAULT
   pdgt=ifif
END SELECT
RETURN
END FUNCTION pdgt

FUNCTION rpdgt(ifif)
! from PDG notation to HELAC-Onia notation
INTEGER,INTENT(IN)::ifif
INTEGER::rpdgt,init=0,ifac
INTEGER,DIMENSION(-100:100)::ipdt
SAVE ipdt,init,ifac
IF(init.EQ.0)THEN
   ipdt(12)= 1
   ipdt(-12)=-1
   ipdt(14)= 5
   ipdt(-14)=-5
   ipdt(16)= 9
   ipdt(-16)=-9
   ipdt(11)= 2
   ipdt(-11)=-2
   ipdt(13)= 6
   ipdt(-13)=-6
   ipdt(15)= 10
   ipdt(-15)=-10
   ipdt(2)= 3
   ipdt(-2)=-3
   ipdt(1)= 4
   ipdt(-1)=-4
   ipdt( 4)= 7
   ipdt(-4)=-7
   ipdt( 3)= 8
   ipdt(-3)=-8
   ipdt( 6)= 11
   ipdt(-6)=-11
   ipdt( 5)= 12
   ipdt(-5)=-12
   ipdt(21)=35
   ipdt(22)=31
   ipdt(23)=32
   ipdt(24)=33
   ipdt(-24)=34
   ipdt(25)=41
   init=1
ENDIF
IF(ABS(ifif).LT.100)THEN
   rpdgt=ipdt(ifif)
   RETURN
ENDIF
ifac=1
IF(ifif.LT.0)ifac=-1
SELECT CASE(ABS(ifif))
   ! J/psi
CASE(443)
   rpdgt=443011
   ! etac
CASE(441)
   rpdgt=441001
   ! hc
CASE(10443)
   rpdgt=441111
   ! chic0
CASE(10441)
   rpdgt=443101
   ! chic1
CASE(20443)
   rpdgt=443111
   ! chic2
CASE(445)
   rpdgt=443121
   ! upsilon
CASE(553)
   rpdgt=553011
   ! etab
CASE(551)
   rpdgt=551001
   ! hb
CASE(10553)
   rpdgt=551111
   ! chib0
CASE(10551)
   rpdgt=553101
   ! chib1
CASE(20553)
   rpdgt=553111
   ! chib2
CASE(555)
   rpdgt=553121
   ! Bc*+-
CASE(543)
   rpdgt=453011*ifac
   ! Bc+-
CASE(541)
   rpdgt=451001*ifac
   ! Bc2*+-
CASE(545)
   rpdgt=453121*ifac
   ! Bc0*+-
CASE(10541)
   rpdgt=453101*ifac
   ! Bc1(H)+-
CASE(20543)
   rpdgt=453111*ifac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
!     To matching the notations in PYTHIA 8                !
!                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! ccbar(3P0(8))
CASE(9910441)
   rpdgt=443108
   ! bbbar(3P0(8))
CASE(9910551)
   rpdgt=553108
   ! bbbar(3S1(8))
CASE(9900553)
   rpdgt=553018
   ! bbbar(1S0(8))
CASE(9900551)
   rpdgt=551008
   ! ccbar(3S1(8))
CASE(9900443)
   rpdgt=443018
   ! ccbar(1S0(8))
CASE(9900441)
   rpdgt=441008
CASE DEFAULT
   rpdgt=ifif
END SELECT
RETURN
END FUNCTION rpdgt
END MODULE Helac_Func_1
