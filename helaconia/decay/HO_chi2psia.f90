MODULE HO_chi2psia
  USE Helac_ranmar_mod
  USE Kinetic_Func
  USE HOVll
  IMPLICIT NONE
CONTAINS
  SUBROUTINE HO_chi12psia_hel(PM,Pd1,Pd2,frame,ihelm,iheld1,iheld2,mchi1,mpsi)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::Pm
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pd1,Pd2
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::frame
    !INTEGER,INTENT(IN)::ihelm
    !INTEGER,INTENT(OUT)::iheld1,iheld2
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    REAL(KIND(1d0)),INTENT(IN)::mchi1,mpsi
    REAL(KIND(1d0))::mass,mass2,ppp
    REAL(KIND(1d0)),DIMENSION(4)::newframe
    ! use CLEO measurement arXiv:0910.0046
    REAL(KIND(1d0)),PARAMETER::aJ10=0.998038696644574d0,aJ11=-0.0626d0
    INCLUDE "HO_SDME_chi12psia.inc"
    mass=DSQRT(PM(0)**2-PM(1)**2-PM(2)**2-PM(3)**2)
    mass2=mass*mpsi/mchi1
    ppp=(mass**2-mass2**2)/(2d0*mass)
    Pd1(0)=DSQRT(ppp**2+mass2**2)
    Pd1(1)=-ppp*DCOS(phi)*sinth
    Pd1(2)=-ppp*DSIN(phi)*sinth
    Pd1(3)=-ppp*costh
    Pd2(0)=ppp
    Pd2(1:3)=-Pd1(1:3)
    PBOO(1:3)=PM(1:3)
    PBOO(4)=PM(0)
    PL(1:3)=Pd1(1:3)
    PL(4)=Pd1(0)
    ! Rotate and boost Pd1
    CALL ROTATEL(frame,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pd1(0)=PL(4)
    Pd1(1:3)=PL(1:3)
    ! Rotate and boost Pd2  
    PL(1:3)=Pd2(1:3)
    PL(4)=Pd2(0)
    CALL ROTATEL(frame,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pd2(0)=PL(4)
    Pd2(1:3)=PL(1:3)
    RETURN
  END SUBROUTINE HO_chi12psia_hel

  SUBROUTINE HO_chi22psia_hel(PM,Pd1,Pd2,frame,ihelm,iheld1,iheld2,mchi2,mpsi)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::Pm
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pd1,Pd2
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::frame
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    REAL(KIND(1d0)),INTENT(IN)::mchi2,mpsi
    REAL(KIND(1d0))::mass,mass2,ppp
    ! use CLEO measurement arXiv:0910.0046
    REAL(KIND(1d0)),PARAMETER::aJ21=0.9956661086930698d0,aJ22=-0.093d0,aJ23=0d0
    INCLUDE "HO_SDME_chi22psia.inc"
    mass=DSQRT(PM(0)**2-PM(1)**2-PM(2)**2-PM(3)**2)
    mass2=mass*mpsi/mchi2
    ppp=(mass**2-mass2**2)/(2d0*mass)
    Pd1(0)=DSQRT(ppp**2+mass2**2)
    Pd1(1)=-ppp*DCOS(phi)*sinth
    Pd1(2)=-ppp*DSIN(phi)*sinth
    Pd1(3)=-ppp*costh
    Pd2(0)=ppp
    Pd2(1:3)=-Pd1(1:3)
    PBOO(1:3)=PM(1:3)
    PBOO(4)=PM(0)
    PL(1:3)=Pd1(1:3)
    PL(4)=Pd1(0)
    ! Rotate and boost Pd1
    CALL ROTATEL(frame,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pd1(0)=PL(4)
    Pd1(1:3)=PL(1:3)
    ! Rotate and boost Pd2
    PL(1:3)=Pd2(1:3)
    PL(4)=Pd2(0)
    CALL ROTATEL(frame,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pd2(0)=PL(4)
    Pd2(1:3)=PL(1:3)
    RETURN
  END SUBROUTINE HO_chi22psia_hel
END MODULE HO_chi2psia
