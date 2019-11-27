MODULE Helac_histo
IMPLICIT NONE
INTEGER,PARAMETER::id=200,nbin=160
REAL(KIND(1d0)),DIMENSION(id,nbin)::H,H2
INTEGER,DIMENSION(id,nbin)::NU
REAL(KIND(1d0)),DIMENSION(id)::HX
INTEGER,DIMENSION(id)::IOO,IU,II
REAL(KIND(1d0)),DIMENSION(id)::Y0,Y1,IC
SAVE H,HX,IOO,IU,II,H2,NU,Y0,Y1,IC
CONTAINS
!----------------------------------------------------------------------
SUBROUTINE Helac_HISTO1(IH,IB,X0,X1,X,W)
INTEGER,INTENT(IN)::IH,IB
REAL(KIND(1d0)),INTENT(IN)::X0,X1,X,W
!CHARACTER(len=1),DIMENSION(30)::REGEL=(/' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '&
!                                      ,' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '&
!									  ,' ',' '/)
!CHARACTER(len=1)::BLANK=' ',STAR='*'
INTEGER::flag,IX
!PRINT *,IH,IB,X0,X1,X,W ! Debug
flag=0
Y0(IH)=X0
Y1(IH)=X1
IC(IH)=IB
IF(X.GE.X0)THEN
  IF(X.GT.X1)THEN
     flag=1
  ELSE
	 IF(X1.EQ.X0)RETURN
     IX=IDINT((X-X0)/(X1-X0)*DBLE(IB))+1
     H(IH,IX)=H(IH,IX)+W
     H2(IH,IX)=H2(IH,IX)+W*W
     NU(IH,IX)=NU(IH,IX)+1
     IF(H(IH,IX).GT.HX(IH)) HX(IH)=H(IH,IX)
     II(IH)=II(IH)+1
     RETURN
  ENDIF
ENDIF
IF(flag.EQ.0)THEN
     IU(IH)=IU(IH)+1
     RETURN
ENDIF
IOO(IH)=IOO(IH)+1
END SUBROUTINE Helac_HISTO1
!----------------------------------------------------------------------
SUBROUTINE Helac_HISTO2(IH,NUNIT,IL,W0)
IMPLICIT NONE
INTEGER,INTENT(IN)::IH,IL,NUNIT
REAL(KIND(1d0)),INTENT(IN)::W0
CHARACTER(len=1),DIMENSION(30)::REGEL=(/' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '&
                                      ,' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '&
									  ,' ',' '/)
CHARACTER(len=1)::BLANK=' ',STAR='*'
INTEGER::IB1,I
REAL(KIND(1d0))::X01,X11
REAL(KIND(1d0))::STEP,Z,v1,v2,v3,RN1,RN2,ERR,Y
INTEGER::IV,IZ
IB1=IC(IH)
X01=Y0(IH)
X11=Y1(IH)
HX(IH)=HX(IH)*(1.D0+1.D-06)
IF(IL.EQ.0) WRITE(6,21) IH,II(IH),IU(IH),IOO(IH)
IF(IL.EQ.1) WRITE(6,22) IH,II(IH),IU(IH),IOO(IH)
21 FORMAT('0NO.',I3,' LIN : INSIDE,UNDER,OVER ',3I6)
22 FORMAT('0NO.',I3,' LOG : INSIDE,UNDER,OVER ',3I6)
IF(II(IH).NE.0)THEN !GOTO 28
   WRITE(6,23)
   23 FORMAT(35(1H ),3(10H----+----I))
   STEP=(X11-X01)/DBLE(IB1)
!---------------------------------------------------
   DO IV=1,IB1      ! 27
!---------------------------------------------------
      Z=DBLE(IV)*STEP+X01-STEP*0.5D0
      IF(IL.EQ.1)THEN
	     IZ=-1
         IF(H(IH,IV).GT.0.D0)IZ=IDINT(DLOG(H(IH,IV))/DLOG(HX(IH))*30.d0)+1
	  ELSE
         IZ=IDINT(H(IH,IV)/HX(IH)*30.d0)+1
      ENDIF 
      IF(IZ.GT.0.AND.IZ.LE.30) REGEL(IZ)=STAR
      WRITE(6,26) Z,H(IH,IV),(REGEL(I),I=1,30)
      26 FORMAT(1H ,2G15.6,4H   I,30A1,1HI)
!*************************************************
      v1=H(IH,IV)
      v2=H2(IH,IV)
      v3=v1**2-v2
      RN1=W0          !NU(IH,IV)
      RN2=RN1*(W0-1)  !(NU(IH,IV)-1)
      IF(RN2.GT.0)THEN
        ERR=(v2/RN1-v3/RN2)/RN1
        ERR=DSQRT(ERR)
      ELSE
        ERR=0
      ENDIF
      Y=0
      IF(RN1.GT.0)Y=H(IH,IV)/W0
      WRITE(NUNIT,216)Z,Y,ERR,NU(IH,IV)
      216 FORMAT(3G15.6,I8)
!*************************************************
      IF(IZ.GT.0.AND.IZ.LE.30) REGEL(IZ)=BLANK
   ENDDO
!   27 CONTINUE
   WRITE(6,23)
   RETURN
ENDIF
WRITE(6,29)  !28
29 FORMAT('  NO ENTRIES INSIDE HISTOGRAM')
END SUBROUTINE Helac_HISTO2

SUBROUTINE Helac_HISTO3(IH)
IMPLICIT NONE
INTEGER,INTENT(IN)::IH
INTEGER::I
DO I=1,nbin
   H(IH,I)=0.d0
ENDDO
HX(IH)=0.d0
IOO(IH)=0
IU(IH)=0
II(IH)=0
END SUBROUTINE Helac_HISTO3
!-----------------------------------------------------------      
END MODULE Helac_histo
