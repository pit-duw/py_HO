MODULE HO_t2bw
  USE Helac_Global
  USE Helac_ranmar_mod
  USE Kinetic_Func
  IMPLICIT NONE
CONTAINS
  SUBROUTINE HO_t2bw_hel(Pm,Pd1,Pd2,frame,ihelm,iheld1,iheld2,anti)
    ! t > b w+ when anti = 1
    ! t~ > b~ w- when anti = -1
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::Pm
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pd1,Pd2
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::frame
    INTEGER,INTENT(IN)::ihelm,anti
    INTEGER,INTENT(OUT)::iheld1,iheld2
    REAL(KIND(1d0))::rrr,rrr2,mt,mw,fac,costh,sinth,ppp,phi
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    iheld1=-1*SIGN(1,anti)
    mt=ABS(Pm(0)**2-Pm(1)**2-Pm(2)**2-Pm(3)**2)
    mt=DSQRT(mt)
    mw=parmas(33)
    IF(mw.GE.mt)THEN
       WRITE(*,*)'You can not decay t -> b W when Mw > Mt!'
       STOP
    ENDIF
    rrr=Helac_rnmy(0)
    fac=2*mw**2/(mt**2+2*mw**2)
    rrr2=Helac_rnmy(0)
    IF(rrr.LT.fac)THEN
       iheld2=-1*SIGN(1,anti)
       IF(ihelm*SIGN(1,anti).GT.0)THEN
          rrr2=2*rrr2-0.5d0
          costh=DSQRT(ABS(1d0+2d0*rrr2))-1d0
       ELSE
          rrr2=2*rrr2-1.5d0
          costh=1d0-DSQRT(ABS(1d0-2d0*rrr2))
       ENDIF
    ELSE
       iheld2=0
       IF(ihelm*SIGN(1,anti).LE.0)THEN
          rrr2=2*rrr2-0.5d0
          costh=DSQRT(ABS(1d0+2d0*rrr2))-1d0
       ELSE
          rrr2=2*rrr2-1.5d0
          costh=1d0-DSQRT(ABS(1d0-2d0*rrr2))
       ENDIF
    ENDIF
    sinth=DSQRT(MIN(1d0-costh**2,1d0))
    phi=2d0*pi*Helac_rnmy(0)
    ppp=(mt**2-mw**2)/(2d0*mt)
    Pd2(0)=DSQRT(ppp**2+mw**2)
    Pd2(1)=ppp*DCOS(phi)*sinth
    Pd2(2)=ppp*DSIN(phi)*sinth
    Pd2(3)=ppp*costh
    Pd1(0)=ppp
    Pd1(1:3)=-Pd2(1:3)
    PBOO(1:3)=Pm(1:3)
    PL(1:3)=Pd1(1:3)
    PL(4)=Pd1(0)
    ! Rotate and boost Pd1
    CALL ROTATEL(frame,PL)
    CALL BOOSTL(mt,PBOO,PL)
    Pd1(0)=PL(4)
    Pd1(1:3)=PL(1:3)
    ! Rotate and boost Pd2
    PL(1:3)=Pd2(1:3)
    PL(4)=Pd2(0)
    CALL ROTATEL(frame,PL)
    CALL BOOSTL(mt,PBOO,PL)
    Pd2(0)=PL(4)
    Pd2(1:3)=PL(1:3)
    RETURN
  END SUBROUTINE HO_t2bw_hel
END MODULE HO_t2bw
