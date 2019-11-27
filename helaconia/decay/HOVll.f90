MODULE HOVll
  USE Helac_ranmar_mod
  USE Kinetic_Func
  IMPLICIT NONE
CONTAINS
  ! Jpsi -> l+ l- with l+,l- massless
  ! It is a vector current
  ! More general decay see EvtGen
  ! Jpsi_+/- -> l+ l- with l+ l- helicity summed
  ! it takes helicity frame
  SUBROUTINE HO_Vll11(PV,squ,Pl1,Pl2)
    ! probability is 3*(1+costh**2)/8 in the rest frame of V
    ! it is azimuthal flat
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PV
    REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::squ
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pl1,Pl2
    REAL(KIND(1d0))::mass,costh,sinth,rrr,phi
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    mass = ABS(PV(0)**2-PV(1)**2-PV(2)**2-PV(3)**2)
    mass = DSQRT(mass)
    rrr=Helac_rnmy(0)-0.5d0
    IF(rrr.LT.-0.5d0)THEN
       rrr=-0.5d0
    ELSEIF(rrr.GT.0.5d0)THEN
       rrr=0.5d0
    ENDIF
    rrr=ABS(4d0*rrr+DSQRT(1d0+16d0*rrr**2))
    rrr=rrr**(1d0/3d0)
    costh=-1d0/rrr+rrr
    sinth=DSQRT(MIN(1d0-costh**2,1d0))
    phi=2d0*pi*Helac_rnmy(0)
    Pl1(0)=mass/2d0
    Pl1(1)=mass/2d0*sinth*DCOS(phi)
    Pl1(2)=mass/2d0*sinth*DSIN(phi)
    Pl1(3)=mass/2d0*costh
    Pl2(0)=mass/2d0
    Pl2(1:3)=-Pl1(1:3)
    PBOO(1:3)=PV(1:3)
    PBOO(4)=PV(0)
    PL(1:3)=Pl1(1:3)
    PL(4)=Pl1(0)
    ! Rotate and boost Pl1
    CALL ROTATEL(squ,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pl1(0)=PL(4)
    Pl1(1:3)=PL(1:3)
    ! Rotate and boost Pl2
    PL(1:3)=Pl2(1:3)
    PL(4)=Pl2(0)
    CALL ROTATEL(squ,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pl2(0)=PL(4)
    Pl2(1:3)=PL(1:3)
    RETURN
  END SUBROUTINE HO_Vll11
  ! W > l nv or q q'
  ! Z > f f~
  ! long.
  SUBROUTINE HO_Vll00(PV,squ,Pl1,Pl2)
    ! probability is 3*(1-costh**2)/4 in the rest frame of V
    ! it is azimuthal flat
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PV
    REAL(KIND(1d0)),DIMENSION(1:3),INTENT(IN)::squ
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pl1,Pl2
    REAL(KIND(1d0))::mass,costh,sinth,rrr,phi
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    mass = ABS(PV(0)**2-PV(1)**2-PV(2)**2-PV(3)**2)
    mass = DSQRT(mass)
    rrr=2d0*(Helac_rnmy(0)-0.5d0)
    IF(rrr.LT.-1d0)THEN
       rrr=-1d0
    ELSEIF(rrr.GT.1d0)THEN
       rrr=1d0
    ENDIF
    costh=SOLVEVll00(rrr)
    sinth=DSQRT(MIN(1d0-costh**2,1d0))
    phi=2d0*pi*Helac_rnmy(0)
    Pl1(0)=mass/2d0
    Pl1(1)=mass/2d0*sinth*DCOS(phi)
    Pl1(2)=mass/2d0*sinth*DSIN(phi)
    Pl1(3)=mass/2d0*costh
    Pl2(0)=mass/2d0
    Pl2(1:3)=-Pl1(1:3)
    PBOO(1:3)=PV(1:3)
    PBOO(4)=PV(0)
    PL(1:3)=Pl1(1:3)
    PL(4)=Pl1(0)
    ! Rotate and boost Pl1
    CALL ROTATEL(squ,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pl1(0)=PL(4)
    Pl1(1:3)=PL(1:3)
    ! Rotate and boost Pl2
    PL(1:3)=Pl2(1:3)
    PL(4)=Pl2(0)
    CALL ROTATEL(squ,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pl2(0)=PL(4)
    Pl2(1:3)=PL(1:3)
    RETURN
  END SUBROUTINE HO_Vll00

  ! transver W > l nue or q q' decay
  SUBROUTINE HO_VllW11(PV,squ,isign,Pl1,Pl2)
    ! probability is 3*(1+isign*costh)**2/8 in the rest frame of V
    ! it is azimuthal flat
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PV
    REAL(KIND(1d0)),DIMENSION(1:3),INTENT(IN)::squ
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pl1,Pl2
    INTEGER,INTENT(IN)::isign
    REAL(KIND(1d0))::mass,costh,sinth,rrr,phi
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    mass = ABS(PV(0)**2-PV(1)**2-PV(2)**2-PV(3)**2)
    mass = DSQRT(mass)
    rrr=Helac_rnmy(0)
    costh=-isign*(1-2d0*rrr**(1d0/3d0))
    sinth=DSQRT(1d0-costh**2)
    phi=2d0*pi*Helac_rnmy(0)
    INCLUDE 'boost.inc'
    RETURN
  END SUBROUTINE HO_VllW11

  FUNCTION SOLVEVll00(rrr)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::rrr
    REAL(KIND(1d0))::SOLVEVll00
    REAL(KIND(1d0))::xr,xi
    COMPLEX(KIND(1d0))::x
    x=DCMPLX(rrr)+DCMPLX(0d0,1d0)*SQRT(1d0-rrr**2)
    x=x**(1d0/3d0)
    xr=DREAL(x)
    xi=DIMAG(x)
    IF(ABS(xr).LE.0.5d0+1d-99)THEN
       SOLVEVll00=-2d0*xr
       SOLVEVll00=DSIGN(MIN(ABS(SOLVEVll00),1d0),SOLVEVll00)
    ELSEIF(ABS(xr+DSQRT(3d0)*xi).LE.1d0+1d-99)THEN
       SOLVEVll00=xr+DSQRT(3d0)*xi
       SOLVEVll00=DSIGN(MIN(ABS(SOLVEVll00),1d0),SOLVEVll00)
    ELSE
       SOLVEVll00=xr-DSQRT(3d0)*xi
       SOLVEVll00=DSIGN(MIN(ABS(SOLVEVll00),1d0),SOLVEVll00)
    ENDIF
    RETURN
  END FUNCTION SOLVEVll00

  ! Jpsi -> l+ l- with l+,l- massless
  SUBROUTINE HO_Vll_unpolarized(PV,Pl1,Pl2)
    ! probability is 1/2 in the rest frame of V
    ! it is an unpolarized decay of V into l+ l-
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PV
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pl1,Pl2
    REAL(KIND(1d0))::mass,costh,sinth,rrr,phi
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    mass = ABS(PV(0)**2-PV(1)**2-PV(2)**2-PV(3)**2)
    mass = DSQRT(mass)
    rrr=2d0*Helac_rnmy(0)-1d0
    costh=rrr
    sinth=DSQRT(MIN(1d0-costh**2,1d0))
    phi=2d0*pi*Helac_rnmy(0)
    Pl1(0)=mass/2d0
    Pl1(1)=mass/2d0*sinth*DCOS(phi)
    Pl1(2)=mass/2d0*sinth*DSIN(phi)
    Pl1(3)=mass/2d0*costh
    Pl2(0)=mass/2d0
    Pl2(1:3)=-Pl1(1:3)
    PBOO(1:3)=PV(1:3)
    PBOO(4)=PV(0)
    PL(1:3)=Pl1(1:3)
    PL(4)=Pl1(0)
    ! boost Pl1 
    CALL BOOSTL(mass,PBOO,PL)
    Pl1(0)=PL(4)
    Pl1(1:3)=PL(1:3)
    ! boost Pl2  
    PL(1:3)=Pl2(1:3)
    PL(4)=Pl2(0)
    CALL BOOSTL(mass,PBOO,PL)
    Pl2(0)=PL(4)
    Pl2(1:3)=PL(1:3)
    RETURN
  END SUBROUTINE HO_Vll_unpolarized
  ! a=int(f1)/(int(f1)+int(f2))
  ! int(f1+f2)=int(f1)+int(f2)=
  ! int(1/a*Theta(a-r)int(f1)+1/(1-a)*Theta(r-a)*int(f2) dr)
  SUBROUTINE HO_Zff11(PV,squ,isign,YL,YR,Pl1,Pl2)
    ! Z > f f
    ! probability is 3*(YL**2*(1+isign*costh)**2
    ! +YR**2*(1-isign*costh)**2)/8/(YL**2+YR**2) in the rest frame of V
    ! it is azimuthal flat
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::YL,YR ! YL=I3-Q*SW**2,YR=-Q*SW**2
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(IN)::PV
    INTEGER,INTENT(IN)::isign
    REAL(KIND(1d0)),DIMENSION(1:3),INTENT(IN)::squ
    REAL(KIND(1d0)),DIMENSION(0:3),INTENT(OUT)::Pl1,Pl2
    REAL(KIND(1d0))::mass,costh,sinth,rrr,phi
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932385d0
    REAL(KIND(1d0)),DIMENSION(4)::PBOO,PL
    REAL(KIND(1d0))::aaa
    mass = ABS(PV(0)**2-PV(1)**2-PV(2)**2-PV(3)**2)
    mass = DSQRT(mass)
    aaa=YR**2/(YL**2+YR**2)
    rrr=Helac_rnmy(0)
    IF(rrr.LE.aaa)THEN
       costh=isign*(1-2d0*Helac_rnmy(0)**(1d0/3d0))
    ELSE
       costh=-isign*(1-2d0*Helac_rnmy(0)**(1d0/3d0))
    ENDIF
    sinth=DSQRT(1d0-costh**2)
    phi=2d0*pi*Helac_rnmy(0)
    INCLUDE 'boost.inc'
  END SUBROUTINE HO_Zff11

  SUBROUTINE Coup_Zff(ifif,YL,YR)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ifif
    REAL(KIND(1d0)),INTENT(OUT)::YL,YR
    REAL(KIND(1d0)),DIMENSION(12)::EWI3,EWQ
    REAL(KIND(1d0))::EWSW,EWCW
    REAL(KIND(1d0)),PARAMETER::mz=91.188d0
    REAL(KIND(1d0)),PARAMETER::mw=80.419d0
    INTEGER::init=0,i,j
    SAVE EWI3,init,EWQ,EWSW,EWCW
    IF(init.EQ.0)THEN
       EWCW=mw/mz
       EWSW=DSQRT(1-EWCW**2)
       DO i=1,12
          j=MOD(i,4)
          IF(j.EQ.0)THEN
             ! down-type quark
             EWI3(i)=-0.5d0
             EWQ(i)=-1d0/3d0
          ELSEIF(j.EQ.3)THEN
             ! up-type quark
             EWI3(i)=0.5d0
             EWQ(i)=2d0/3d0
          ELSEIF(j.EQ.1)THEN
             ! neutrino
             EWI3(i)=0.5d0
             EWQ(i)=0d0
          ELSE
             ! charged lepton
             EWI3(i)=-0.5d0
             EWQ(i)=-1d0
          ENDIF
       ENDDO
       init=0d0
    ENDIF
    YR=-EWQ(ABS(ifif))*EWSW**2
    YL=EWI3(ABS(ifif))+YR
    RETURN
  END SUBROUTINE Coup_Zff

  SUBROUTINE MCviaf(fun,x,fmax_in)
    IMPLICIT NONE
    REAL(KIND(1d0)),EXTERNAL::fun
    REAL(KIND(1d0)),INTENT(OUT)::x
    REAL(KIND(1d0)),INTENT(IN),OPTIONAL::fmax_in
    REAL(KIND(1d0))::fmax,rrr,fprob
    LOGICAL::lexit
    IF(PRESENT(fmax_in))THEN
       fmax=MAX(fmax_in,1d-99)
    ELSE
       fmax=1d0
    ENDIF
    lexit=.FALSE.
    DO WHILE(.NOT.lexit)
       x=Helac_rnmy(0)
       fprob=fun(x)/fmax
       rrr=Helac_rnmy(0)
       IF(fprob.GT.rrr)THEN
          lexit=.TRUE.
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE MCviaf
END MODULE HOVLL
