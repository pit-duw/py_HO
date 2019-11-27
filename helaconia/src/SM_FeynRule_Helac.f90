MODULE Helac_SM_FeynRule
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  This is a package physics from helac_phegas 1.2.2.
!  
!  In this routine, all couplings of the standard electroweak theory are defined.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
USE Helac_Global
USE Helac_Func_1
USE Structf_PDFs
USE Constants
IMPLICIT NONE
CONTAINS
SUBROUTINE Helac_FeynRule_SM()
!  CHARGE and ISOSPIN
REAL(KIND=DBL),DIMENSION(4)::q,t
!  CKM MATRIX
COMPLEX(KIND=DBL),DIMENSION(3,3)::zckm
REAL(KIND=DBL)::root2,gqcd,rsw2,alpha,chi,hm,e
COMPLEX(KIND=DBL)::zi,zero,zwm,zzm,zgw
COMPLEX(KIND=DBL)::zmas,zsinw2,zsinw,zcosw
INTEGER::iuni,ihi,mm,i,l,j,iqed,nunittemp
INTEGER::init=0,iqcdqcd
REAL(KIND=DBL)::pi
SAVE     
!  -----------------------------------------------------------------------
!  start initialize
!  -----------------------------------------------------------------------
CALL Helac_mypi(pi)
IF(init.EQ.0)THEN
! CKM initialize
   zckm(1:3,1:3)=0
   DO i=1,3
     zckm(i,i)=1
   ENDDO

   WRITE(*,*)' '
   WRITE(*,*)'---------------------------------------------------------------------'
   WRITE(*,*)'WE ARE LOADING THE STANDARD MODEL FEYNMAN RULES ...'
   WRITE(*,*)' '
   CALL ReadElem_integer('qcd',iqcdqcd)
   SELECT CASE(iqcdqcd)
   ! only electroweak
   CASE(0)
	 onlyqcd=.FALSE.
	 withqcd=.FALSE.
	 iqed=0
   ! electroweak and qcd
   CASE(1)
     onlyqcd=.FALSE.
	 withqcd=.TRUE.
	 iqed=0
   ! qcd
   CASE(2)
     onlyqcd=.TRUE.
	 withqcd=.FALSE.
	 iqed=0
   ! qed
   CASE(3)
	 onlyqcd=.FALSE.
	 withqcd=.FALSE.
	 iqed=1
   ! qcd and qed
   CASE(4)
     onlyqcd=.FALSE.
	 withqcd=.TRUE.
	 iqed=1
   END SELECT
   WRITE(*,*)'onlyqcd: T(1) denotes only in QCD , F(0) denotes the oppsite side '
!   READ(*,*)onlyqcd
   WRITE(*,*)'withqcd: T(1) denotes with QCD and Electroweak, F(0) denotes the oppsite '
!   READ(*,*)withqcd
   WRITE(*,*)'irun : 1 the qcd running coupling constant, 0 fixed qcd coupling constant '
!   READ(*,*)irun
   CALL ReadElem_integer('alphasrun',irun)
   CALL ReadElem_logic('reweight_Scale',reweight_scale)
   reweight_scale=reweight_scale.AND.(irun.EQ.1)
   IF(reweight_scale)THEN
      CALL ReadElem_real('rw_RScale_down',rw_Rscale_down)
      CALL ReadElem_real('rw_RScale_up',rw_Rscale_up)
      CALL ReadElem_real('rw_FScale_down',rw_Fscale_down)
      CALL ReadElem_real('rw_FScale_up',rw_Fscale_up)
      ho_nscale=8
      IF(rw_Rscale_down.LE.0d0.OR.rw_Rscale_up.LE.0d0.OR.&
           rw_Fscale_down.LE.0d0.OR.rw_Fscale_up.LE.0d0)THEN
         reweight_scale=.FALSE.
         ho_nscale=0
         WRITE(*,*)"WARNING:Scale reweighting is off."
         WRITE(*,*)"WARNING:Please make sure rw_Rscale_up, rw_Rscale_down, rw_Fscale_up, rw_Fscale_down are larger than 0."
      ENDIF
      IF(reweight_scale.AND.rw_Rscale_down.GT.rw_Rscale_up)THEN
         gqcd=rw_Rscale_down
         rw_Rscale_down=rw_Rscale_up
         rw_Rscale_up=gqcd
      ENDIF
      IF(reweight_scale.AND.rw_Fscale_down.GT.rw_Fscale_up)THEN
         gqcd=rw_Fscale_down
         rw_Fscale_down=rw_Fscale_up
         rw_Fscale_up=gqcd
      ENDIF
   ENDIF
   WRITE(*,*)' '
   WRITE(*,*)'onlyqcd value=',onlyqcd,'   withqcd value= ',withqcd,'   irun value=',irun
   WRITE(*,*)' '

   root2=DSQRT(dnou(2))
!CALL mypi(pi)
! the unit of imag and zero
   zi=DCMPLX(dnou(0),dnou(1))
   zero=DCMPLX(dnou(0),dnou(0))
       
   IF(iunitary.EQ.0)iuni=1
   IF(iunitary.EQ.1)iuni=0
   IF(ihiggs.EQ.0)ihi=0
   IF(ihiggs.EQ.1)ihi=1
!   WRITE(*,*)'iunitary,ihiggs,iuni,ihi'
!   WRITE(*,*)iunitary,ihiggs,iuni,ihi
      
   zgv3(31:35,31:35,31:35)=zero
   zgv4(31:35,31:35,31:35,31:35)=zero
   zgvffl(31:35,-12:-1,1:12)=zero
   zgvffr(31:35,-12:-1,1:12)=zero
   zgvvs(31:35,31:35,41:44)=zero
   zgsvv(41:44,31:35,31:35)=zero
   zgvvss(31:35,31:35,41:44,41:44)=zero
   zgssvv(41:44,41:44,31:35,31:35)=zero
   zgvss(31:35,41:44,41:44)=zero
   zgssv(41:44,41:44,31:35)=zero
   zgsffl(41:44,-12:-1,1:12)=zero
   zgsffr(41:44,-12:-1,1:12)=zero
   zgs3(41:44,41:44,41:44)=zero
   zgs4(41:44,41:44,41:44,41:44)=zero
   parmas(-12:44)=0
   parwid(-12:44)=0

	CALL ReadConst
! part from constants.h     
   rsw2=dnou(1)-rwm**2/rzm**2    ! squre of sine of the weinberg angle
!  G_F/Sqrt(2)=g^2/8/mw^2 & g=e/s_w & alpha=e^2/4/pi
   IF(alphaem.EQ.-1)THEN
		alpha=DSQRT(dnou(2))*gfermi*rwm*rwm*rsw2/pi
   ELSE
		alpha=alphaem
   ENDIF
   gqcd=0
   IF(withqcd)gqcd=DSQRT(alphaQCD2*4*pi)
   IF(onlyqcd)gqcd=DSQRT(alphaQCD2*4*pi)
   parmas(32)=rzm                  ! the mass of Z in GeV
   parwid(32)=rzw                  ! the width of Z in GeV
   parmas(33)=rwm                  ! the mass of W in GeV
   parwid(33)=rww                  ! the width of W in GeV
   parmas(34)=parmas(33)           ! the mass of W in GeV
   parwid(34)=parwid(33)           ! the width of W in GeV
   parmas(41)=rhm                  ! the mass of higgs in GeV
   parwid(41)=rhw                  ! the width of higgs in GeV
   parmas(1)=0.0d0                 ! the mass of nu_e in GeV
   parmas(2)=0.0d0                 ! the mass of electron in GeV
   parmas(3)=0.0d0                 ! the mass of up quark in GeV
   parmas(4)=0.0d0                 ! the mass of down quark in GeV
   parmas(5)=0.0d0                 ! the mass of nu_mu in GeV
   parmas(6)=0.0d0                 ! the mass of mu in GeV
   parmas(7)=rcm                  ! the mass of charm quark in GeV
   parmas(8)=0.0d0                 ! the mass of strange quark in GeV
   parmas(9)=0.0d0                 ! the mass of nu_tau in GeV
   parmas(10)=0.0d0                ! the mass of tau in GeV
   parmas(11)=rtm              ! the mass of top quark in GeV
   parwid(11)=rtw                ! the width of top quark in GeV
   parmas(12)=rbm                ! the mass of bottom quark in GeV
   aqedup=alpha
   aqcdup=gqcd**2/4d0/pi

! HHH coupling
   chi=1d0       
   mm=33
! compl_mass.h
   IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
      zmas=parmas(mm)*DCMPLX(DSQRT(dnou(1)/dnou(2)&
           +dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2))&
       , -DSQRT(-dnou(1)/dnou(2)+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)) )
   ELSE
      zmas=DCMPLX(parmas(mm),dnou(0))
   ENDIF
   zwm=zmas
   mm=32
! compl_mass.h
   IF(parmas(mm).GT.dnou(0).AND.iwidth.EQ.1)THEN
      zmas=parmas(mm)*DCMPLX(DSQRT(dnou(1)/dnou(2)&
      +dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2))&
      , -DSQRT(-dnou(1)/dnou(2)+dnou(1)/dnou(2)*DSQRT(dnou(1)+parwid(mm)**2/parmas(mm)**2)) )
   ELSE
      zmas=DCMPLX(parmas(mm),dnou(0))
   ENDIF
   zzm=zmas
! the mass of goldstone bosons
   parmas(42)=parmas(32)
   parwid(42)=parwid(32)
   parmas(43)=parmas(33)
   parwid(43)=parwid(33)
   parmas(44)=parmas(34)
   parwid(44)=parwid(34)
   nunittemp=201021
   CLOSE(nunittemp)
   OPEN(nunittemp,FILE=TRIM(output_dir)//'MassWidth.out')   
   WRITE(nunittemp,*)' '
   WRITE(nunittemp,*)' The Mass and Width parameters of various particles '
   WRITE(nunittemp,*)'part. index','   Mass(GeV)   ','  Width(GeV)'    
   DO i=1,12
      WRITE(nunittemp,'(i8,5x,e10.4,5x,e10.4)')i,parmas(i),parwid(i)
   ENDDO
   DO i=31,35
      WRITE(nunittemp,'(i8,5x,e10.4,5x,e10.4)')i,parmas(i),parwid(i)
   ENDDO
   DO i=41,44
      WRITE(nunittemp,'(i8,5x,e10.4,5x,e10.4)')i,parmas(i),parwid(i)
   ENDDO
   CLOSE(nunittemp,STATUS='KEEP')

   hm=parmas(41)
        
   DO l=1,12
      parmas(-l)=parmas(l)
      parwid(-l)=parwid(l)
   ENDDO
       
   init=1
ENDIF

!  -----------------------------------------------------------------------
!  end initialize
!  -----------------------------------------------------------------------
       
!  -----------------------------------------------------------------------
!  start QCD couplings
!  -----------------------------------------------------------------------
!      include 'qcdscale.h'
IF(init.EQ.1)THEN
  WRITE(*,*)' '
  IF(gqcd.EQ.0d0)WRITE(*,*)'QCD COUPLING SET TO ZERO'
  IF(gqcd.GT.0d0)WRITE(*,*)'QCD INCLUDED, g=',gqcd,gqcd**2/4/pi
  init=2
!      elseif(init.ge.2.and.istruc.eq.1) then
ELSEIF(init.GE.2.AND.irun.EQ.1) THEN
  CALL qcdscale(scale)
!  PRINT *,scale,parmas(32)
!  PRINT *,ALPHAS(scale),ALPHAS(parmas(32))
!  STOP
  aqcdup=ALPHAS(scale)
  IF(reweight_scale_phase)aqcdup=aqcdup/10d0
  gqcd=DSQRT(aqcdup*4*pi)  ! ALPHAS in Helac_Func_1.f90
  init=3
ENDIF
IF(init.EQ.2)THEN
	WRITE(*,*)'---------------------------------------------------------------------------'  
ENDIF 
     
zgvffl(35,-3,3)=DCMPLX(gqcd/root2)
zgvffr(35,-3,3)=DCMPLX(gqcd/root2)
zgvffl(35,-4,4)=DCMPLX(gqcd/root2)
zgvffr(35,-4,4)=DCMPLX(gqcd/root2)
zgvffl(35,-7,7)=DCMPLX(gqcd/root2)
zgvffr(35,-7,7)=DCMPLX(gqcd/root2)
zgvffl(35,-8,8)=DCMPLX(gqcd/root2)
zgvffr(35,-8,8)=DCMPLX(gqcd/root2)
zgvffl(35,-11,11)=DCMPLX(gqcd/root2)
zgvffr(35,-11,11)=DCMPLX(gqcd/root2)
zgvffl(35,-12,12)=DCMPLX(gqcd/root2)
zgvffr(35,-12,12)=DCMPLX(gqcd/root2)
zgv3(35,35,35)=-DCMPLX(gqcd/4d0*DSQRT(dnou(2))**3)
zgv4(35,35,35,35)= DCMPLX(gqcd**2/8d0*DSQRT(dnou(2))**4)
       
! IF ONLY QCD THEN
! Including QED
!e=DSQRT(dnou(4)*pi*alpha)
! charges
!q(1)= dnou(0)
!q(2)=-dnou(1)
!q(3)= dnou(2)/dnou(3)
!q(4)=-dnou(1)/dnou(3)
!zgvffl(31,-2,2)=-DCMPLX(e*q(2))
!zgvffl(31,-3,3)=-DCMPLX(e*q(3))
!zgvffl(31,-4,4)=-DCMPLX(e*q(4))

!zgvffl(31,-6,6)=zgvffl(31,-2,2)
!zgvffl(31,-7,7)=zgvffl(31,-3,3)
!zgvffl(31,-8,8)=zgvffl(31,-4,4)

!zgvffl(31,-10,10)=zgvffl(31,-6,6)
!zgvffl(31,-11,11)=zgvffl(31,-7,7)
!zgvffl(31,-12,12)=zgvffl(31,-8,8)
      
!zgvffr(31,-2,2)=-DCMPLX(e*q(2))
!zgvffr(31,-3,3)=-DCMPLX(e*q(3))
!zgvffr(31,-4,4)=-DCMPLX(e*q(4))

!zgvffr(31,-6,6)=zgvffr(31,-2,2)
!zgvffr(31,-7,7)=zgvffr(31,-3,3)
!zgvffr(31,-8,8)=zgvffr(31,-4,4)

!zgvffr(31,-10,10)=zgvffr(31,-6,6)
!zgvffr(31,-11,11)=zgvffr(31,-7,7)
!zgvffr(31,-12,12)=zgvffr(31,-8,8)
IF(onlyqcd)RETURN
       
!  -----------------------------------------------------------------------
!  end QCD couplings
!  -----------------------------------------------------------------------
!  no running of couplings: next line does only for qcd
!  -----------------------------------------------------------------------
IF(init.EQ.3)RETURN
!  -----------------------------------------------------------------------

zsinw2=dnou(1)-zwm**2/zzm**2
zsinw =CDSQRT(zsinw2)
zcosw =CDSQRT(DCMPLX(dnou(1),dnou(0))-zsinw2)
       
e=DSQRT(dnou(4)*pi*alpha)
zgw=DCMPLX(e)/zsinw
! charges
q(1)= dnou(0)
q(2)=-dnou(1)
q(3)= dnou(2)/dnou(3)
q(4)=-dnou(1)/dnou(3)
! isospins      
t(1)=+dnou(1)/dnou(2)
t(2)=-dnou(1)/dnou(2)
t(3)=+dnou(1)/dnou(2)
t(4)=-dnou(1)/dnou(2)
! V V V        
zgv3(31,33,34)=-DCMPLX(e) 
zgv3(33,31,33)=-DCMPLX(e) 
zgv3(34,34,31)=-DCMPLX(e) 

zgv3(32,33,34)= zgw*zcosw
zgv3(33,32,33)= zgw*zcosw
zgv3(34,34,32)= zgw*zcosw
! V V V V    
zgv4(31,31,33,34)=-DCMPLX(e**2)                  !AAW+W-
zgv4(31,32,33,34)= DCMPLX(e**2)*zcosw/zsinw      !AZW+W-
     
zgv4(32,32,33,34)=-DCMPLX(e**2)*zcosw**2/zsinw**2  !ZZW+W-
zgv4(32,31,33,34)= DCMPLX(e**2)*zcosw/zsinw        !ZAW+W-
       
zgv4(33,34,33,33)= zgw*zgw                        !W+W-W+W+
zgv4(33,33,32,32)=-zgw*zgw*zcosw**2               !W+W+ZZ
zgv4(33,33,31,32)= DCMPLX(e**2)*zcosw/zsinw       !W+W+AZ
zgv4(33,33,31,31)=-DCMPLX(e**2)                   !W+W+AA
       
zgv4(34,33,34,34)= zgw*zgw                        !W-W+W-W-
zgv4(34,34,32,32)=-zgw*zgw*zcosw**2               !W-W-ZZ
zgv4(34,34,31,32)= DCMPLX(e**2)*zcosw/zsinw       !W-W-AZ
zgv4(34,34,31,31)=-DCMPLX(e**2)                   !W-W-AA
! V f fbar       
zgvffl(31,-2,2)=-DCMPLX(e*q(2))
zgvffl(31,-3,3)=-DCMPLX(e*q(3))
zgvffl(31,-4,4)=-DCMPLX(e*q(4))

zgvffl(31,-6,6)=zgvffl(31,-2,2)
zgvffl(31,-7,7)=zgvffl(31,-3,3)
zgvffl(31,-8,8)=zgvffl(31,-4,4)

zgvffl(31,-10,10)=zgvffl(31,-6,6)
zgvffl(31,-11,11)=zgvffl(31,-7,7)
zgvffl(31,-12,12)=zgvffl(31,-8,8)
      
zgvffr(31,-2,2)=-DCMPLX(e*q(2))
zgvffr(31,-3,3)=-DCMPLX(e*q(3))
zgvffr(31,-4,4)=-DCMPLX(e*q(4))

zgvffr(31,-6,6)=zgvffr(31,-2,2)
zgvffr(31,-7,7)=zgvffr(31,-3,3)
zgvffr(31,-8,8)=zgvffr(31,-4,4)

zgvffr(31,-10,10)=zgvffr(31,-6,6)
zgvffr(31,-11,11)=zgvffr(31,-7,7)
zgvffr(31,-12,12)=zgvffr(31,-8,8)
        
zgvffl(32,-1,1)=(DCMPLX(t(1))-DCMPLX(q(1))*zsinw**2)/zcosw*zgw 
zgvffl(32,-2,2)=(DCMPLX(t(2))-DCMPLX(q(2))*zsinw**2)/zcosw*zgw 
zgvffl(32,-3,3)=(DCMPLX(t(3))-DCMPLX(q(3))*zsinw**2)/zcosw*zgw 
zgvffl(32,-4,4)=(DCMPLX(t(4))-DCMPLX(q(4))*zsinw**2)/zcosw*zgw 

zgvffl(32,-5,5)=zgvffl(32,-1,1)
zgvffl(32,-6,6)=zgvffl(32,-2,2)
zgvffl(32,-7,7)=zgvffl(32,-3,3)
zgvffl(32,-8,8)=zgvffl(32,-4,4)
      
zgvffl(32,-9,9)  =zgvffl(32,-1,1)
zgvffl(32,-10,10)=zgvffl(32,-2,2)
zgvffl(32,-11,11)=zgvffl(32,-3,3)
zgvffl(32,-12,12)=zgvffl(32,-4,4)
      
zgvffr(32,-1,1)=-DCMPLX(q(1))*zsinw**2/zcosw*zgw
zgvffr(32,-2,2)=-DCMPLX(q(2))*zsinw**2/zcosw*zgw 
zgvffr(32,-3,3)=-DCMPLX(q(3))*zsinw**2/zcosw*zgw
zgvffr(32,-4,4)=-DCMPLX(q(4))*zsinw**2/zcosw*zgw

zgvffr(32,-5,5)=zgvffr(32,-1,1)
zgvffr(32,-6,6)=zgvffr(32,-2,2)
zgvffr(32,-7,7)=zgvffr(32,-3,3)
zgvffr(32,-8,8)=zgvffr(32,-4,4)
      
zgvffr(32,-9,9)  =zgvffr(32,-1,1)
zgvffr(32,-10,10)=zgvffr(32,-2,2)
zgvffr(32,-11,11)=zgvffr(32,-3,3)
zgvffr(32,-12,12)=zgvffr(32,-4,4)
         
zgvffl(34,-1,2)=zgw/DCMPLX(root2)                    ! no lepton mixing
zgvffl(34,-3,4)=zgw/DCMPLX(root2)*DCONJG(zckm(1,1))
zgvffl(34,-3,8)=zgw/DCMPLX(root2)*DCONJG(zckm(1,2))
zgvffl(34,-3,12)=zgw/DCMPLX(root2)*DCONJG(zckm(1,3))

zgvffl(34,-5,6)=zgvffl(34,-1,2)
zgvffl(34,-7,4)=zgw/DCMPLX(root2)*DCONJG(zckm(2,1))
zgvffl(34,-7,8)=zgw/DCMPLX(root2)*DCONJG(zckm(2,2))
zgvffl(34,-7,12)=zgw/DCMPLX(root2)*DCONJG(zckm(2,3))
       
zgvffl(34,-9,10)=zgvffl(34,-1,2) 
zgvffl(34,-11,4)=zgw/DCMPLX(root2)*DCONJG(zckm(3,1))
zgvffl(34,-11,8)=zgw/DCMPLX(root2)*DCONJG(zckm(3,2))
zgvffl(34,-11,12)=zgw/DCMPLX(root2)*DCONJG(zckm(3,3))
       
zgvffl(33,-2,1)=zgw/DCMPLX(root2)
zgvffl(33,-4,3)=zgw/DCMPLX(root2)*zckm(1,1) 
zgvffl(33,-8,3)=zgw/DCMPLX(root2)*zckm(2,1) 
zgvffl(33,-12,3)=zgw/DCMPLX(root2)*zckm(3,1) 
       
zgvffl(33,-6,5)=zgvffl(33,-2,1)
zgvffl(33,-4,7)=zgw/DCMPLX(root2)*zckm(1,2)
zgvffl(33,-8,7)=zgw/DCMPLX(root2)*zckm(2,2)
zgvffl(33,-12,7)=zgw/DCMPLX(root2)*zckm(3,2)
       
zgvffl(33,-10,9)=zgvffl(33,-2,1)
zgvffl(33,-4,11)=zgw/DCMPLX(root2)*zckm(1,3)
zgvffl(33,-8,11)=zgw/DCMPLX(root2)*zckm(2,3)
zgvffl(33,-12,11)=zgw/DCMPLX(root2)*zckm(3,3)
! S V V      
zgsvv(41,33,34)= zgw*zwm*DCMPLX(ihi)                                    !HW+W-
zgsvv(41,32,32)= zgw*zwm/zcosw**2 *DCMPLX(ihi)                          !HZZ
      
zgsvv(43,33,32)=-DCMPLX(e)*zwm*zsinw/zcosw*DCMPLX(iuni)                 !phi+W+Z
zgsvv(43,33,31)=-DCMPLX(e)*zwm*DCMPLX(iuni)                             !phi+W+gamma
        
zgsvv(44,34,32)=-DCMPLX(e)*zwm*zsinw/zcosw*DCMPLX(iuni)                 !phi-W-Z
zgsvv(44,34,31)=-DCMPLX(e)*zwm*DCMPLX(iuni)                             !phi-W-gamma
     
zgvvs(31,33,44)=-DCMPLX(e)*zwm*DCMPLX(iuni)                             !gamma W+phi-
zgvvs(31,34,43)=-DCMPLX(e)*zwm*DCMPLX(iuni)                             !gamma W-phi+
       
zgvvs(32,32,41)= zgw*zwm/zcosw**2*DCMPLX(ihi)                           !ZZH
zgvvs(32,33,44)=-DCMPLX(e)*zwm*zsinw/zcosw*DCMPLX(iuni)                 !ZW+phi-
zgvvs(32,34,43)=-DCMPLX(e)*zwm*zsinw/zcosw*DCMPLX(iuni)                 !ZW-phi+
      
zgvvs(33,33,41)= zgw*zwm*DCMPLX(ihi)                                    !W+W+H
zgvvs(33,31,43)=-DCMPLX(e)*zwm*DCMPLX(iuni)                             !W+gammaphi+
zgvvs(33,32,43)=-DCMPLX(e)*zwm*zsinw/zcosw*DCMPLX(iuni)                 !W+Zphi+
      
zgvvs(34,34,41)= zgw*zwm*DCMPLX(ihi)                                    !W-W-H
zgvvs(34,31,44)=-DCMPLX(e)*zwm*DCMPLX(iuni)                             !W-gamma phi-
zgvvs(34,32,44)=-DCMPLX(e)*zwm*zsinw/zcosw*DCMPLX(iuni)                 !W-Zphi-
!V V S S      
zgvvss(34,34,41,41)=zgw**2/DCMPLX(dnou(2))*DCMPLX(ihi)                  !W-W-HH
zgvvss(34,34,42,42)=zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)                 !W-W-XX
zgvvss(34,34,43,44)=zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)                 !W-W-phi+phi-
zgvvss(34,31,44,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni)  !W- gamma phi- H
zgvvss(34,31,44,42)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !W- gamma phi- X
zgvvss(34,32,44,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)  !W-Zphi-H
zgvvss(34,32,44,42)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) !W-Zphi-X

zgvvss(33,33,41,41)=zgw**2/DCMPLX(dnou(2))*DCMPLX(ihi)                  !W+W+HH
zgvvss(33,33,42,42)=zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)                 !W+W+XX
zgvvss(33,33,43,44)=zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)                 !W+W+phi+phi-
zgvvss(33,31,43,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni)  !W+ gamma phi+ H
zgvvss(33,31,43,42)=zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni)!W+ gamma phi+ X
zgvvss(33,32,43,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)  !W+ Z phi+ H
zgvvss(33,32,43,42)=zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)!W+ Z phi+ X

zgvvss(32,32,41,41)=zgw**2/(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(ihi)       !ZZHH
zgvvss(32,32,42,42)=zgw**2/(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(iuni)      !ZZXX
zgvvss(32,32,43,44)=zgw**2*(zsinw**2-zcosw**2)**2&                      !ZZphi+phi-
      /(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(iuni)
zgvvss(32,31,43,44)=DCMPLX(e**2)*(zsinw**2-zcosw**2)/(zsinw*zcosw)*DCMPLX(iuni) !Z gamma phi+ phi-
zgvvss(32,33,44,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)  !Z W+ phi- H
zgvvss(32,34,43,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)  !Z W- phi+ H
zgvvss(32,33,44,42)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) !Z W+ phi- X
zgvvss(32,34,43,42)=zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) !Z W- phi+ X

zgvvss(31,31,43,44)=DCMPLX(dnou(2)*e**2*iuni)                          !gamma gamma phi+ phi-                     
zgvvss(31,32,43,44)=DCMPLX(e**2)*(zsinw**2-zcosw**2)/(zsinw*zcosw)*DCMPLX(iuni) !gamma Z phi+ phi-
zgvvss(31,33,44,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !gamma W+ phi- H
zgvvss(31,34,43,41)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !gamma W- phi+ H      
zgvvss(31,33,44,42)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !gamma W+ phi- X    
zgvvss(31,34,43,42)= zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !gamma W- phi+ X          
! V S S
zgvss(31,43,44)=-DCMPLX(e*iuni)                                        ! A phi+ phi- 

zgvss(32,42,41)=-zi*zgw/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)           ! Z X H 
zgvss(32,43,44)=-zgw*(zsinw**2-zcosw**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) ! Z phi+ phi- 
      
zgvss(33,43,41)= zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                      ! W+ phi+ H 
zgvss(33,43,42)=-zi*zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                   ! W+ phi+ X 
     
zgvss(34,44,41)=-zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                      ! W- phi- H
zgvss(34,44,42)=-zi*zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                   ! W- phi- X
       
zgssv(41,42,32)= zi*zgw/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)           ! H X Z  
zgssv(41,43,34)=-zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                      ! H phi+ W-
zgssv(41,44,33)= zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                      ! H phi- W+
 
zgssv(42,41,32)=-zi*zgw/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)           ! X H Z 
zgssv(42,44,33)= zi*zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                   ! X phi- W+  
zgssv(42,43,34)= zi*zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                   ! X phi+ W-
       
zgssv(43,43,31)= DCMPLX(e*iuni)                                        ! phi+ phi+ A
zgssv(43,43,32)= zgw*(zsinw**2-zcosw**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) ! phi+ phi+ Z
zgssv(43,41,33)=-zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                      ! phi+ H W+ 
zgssv(43,42,33)=-zi*zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                   ! phi+ X W+ 

zgssv(44,44,31)=-DCMPLX(e*iuni)                                        ! phi- phi- A 
zgssv(44,44,32)=-zgw*(zsinw**2-zcosw**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) !phi- phi- Z 
zgssv(44,41,34)= zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                      ! phi- H W- 
zgssv(44,42,34)=-zi*zgw/DCMPLX(dnou(2))*DCMPLX(iuni)                   ! phi- X W- 
! S S V V      
zgssvv(41,41,32,32)= zgw**2/(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(ihi)     !HHZZ
zgssvv(41,41,33,34)= zgw**2/DCMPLX(dnou(2))*DCMPLX(ihi)                !HHW+W-
zgssvv(41,43,34,31)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !H phi+ W- A 
zgssvv(41,44,33,31)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) !H phi- W+ A 
zgssvv(41,43,34,32)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) !H phi+ W- Z 
zgssvv(41,44,33,32)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) !H phi- W+ Z 
       
zgssvv(42,42,32,32)= zgw**2/(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(iuni)    !X X Z Z 
zgssvv(42,42,33,34)= zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)               !X X W+ W- 
zgssvv(42,43,34,31)= zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) ! X phi+ W- A 
zgssvv(42,44,33,31)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) ! X phi- W+ A
zgssvv(42,43,34,32)= zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) ! X phi+ W- Z
zgssvv(42,44,33,32)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) ! X phi- W+ Z

zgssvv(43,43,31,31)= DCMPLX(dnou(2)*e**2*iuni)                          ! phi+ phi+ A A 
zgssvv(43,43,32,31)= DCMPLX(e**2)*(zsinw**2-zcosw**2)/(zsinw*zcosw)*DCMPLX(iuni) ! phi+ phi+ Z A 
zgssvv(43,43,32,32)= zgw**2*(zsinw**2-zcosw**2)**2&                     ! phi+ phi+ Z Z
                /(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(iuni) 
zgssvv(43,43,33,34)= zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)                ! phi+ phi+ W+ W- 
zgssvv(43,41,33,31)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni)  ! phi+ H W+ A 
zgssvv(43,42,33,31)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) ! phi+ X W+ A 
zgssvv(43,41,33,32)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)  ! phi+ H W+ Z 
zgssvv(43,42,33,32)=-zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) ! phi+ X W+ Z
       
zgssvv(44,44,31,31)=DCMPLX(dnou(2)*e**2*iuni)                           ! phi- phi- A A 
zgssvv(44,44,32,31)=DCMPLX(e**2)*(zsinw**2-zcosw**2)/(zsinw*zcosw)*DCMPLX(iuni) ! phi+ phi+ Z A 
zgssvv(44,44,32,32)=zgw**2*(zsinw**2-zcosw**2)**2&                      ! phi- phi- Z Z
               /(DCMPLX(dnou(2))*zcosw**2)*DCMPLX(iuni) 
zgssvv(44,44,33,34)=zgw**2/DCMPLX(dnou(2))*DCMPLX(iuni)                 ! phi- phi- W+ W- 
zgssvv(44,41,34,31)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni)  ! phi- H W- A 
zgssvv(44,42,34,31)= zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zsinw)*DCMPLX(iuni) ! phi- X W- A 
zgssvv(44,41,34,32)=-DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni)  ! phi- H W- Z 
zgssvv(44,42,34,32)= zi*DCMPLX(e**2)/(DCMPLX(dnou(2))*zcosw)*DCMPLX(iuni) ! phi- X W- Z 
! S f fbar      
DO i=1,12
   zgsffl(41,-i,i)=-DCMPLX(e)/(DCMPLX(dnou(2))*zsinw)/zwm&              !H f fbar
             *DCMPLX(parmas(i))*DCMPLX(ihi) 
   zgsffr(41,-i,i)=-DCMPLX(e)/(DCMPLX(dnou(2))*zsinw)/zwm&
             *DCMPLX(parmas(i))*DCMPLX(ihi)
ENDDO
       
DO i=1,12
   j=MOD(i-1,4)+1
   zgsffl(42,-i,i)=-zi*DCMPLX(e)/(DCMPLX(dnou(2))*zsinw) &               ! X f fbar
        *(DCMPLX(dnou(2)*t(j)))/zwm*DCMPLX(parmas(i)*iuni) 
   zgsffr(42,-i,i)= zi*DCMPLX(e)/(DCMPLX(dnou(2))*zsinw) &
        *DCMPLX(dnou(2)*t(j))/zwm*DCMPLX(parmas(i)*iuni) 
ENDDO
       
zgsffr(44,-1,2)= -zgw*DCMPLX(parmas(2))/(DCMPLX(root2)*zwm)*DCMPLX(iuni) ! phi- -nu_e e       
zgsffr(44,-5,6)= -zgw*DCMPLX(parmas(6))/(DCMPLX(root2)*zwm)*DCMPLX(iuni) ! phi- -nu_mu mu 
zgsffr(44,-9,10)=-zgw*DCMPLX(parmas(10))/(DCMPLX(root2)*zwm)*DCMPLX(iuni)! phi- -nu_tau tau 

zgsffl(44,-1,2)= zero
zgsffl(44,-5,6)= zero
zgsffl(44,-9,10)=zero
       
zgsffr(44,-3,4)=  -zgw*DCMPLX(parmas(4))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(1,1)    ! phi- -u d 
zgsffr(44,-7,4)=  -zgw*DCMPLX(parmas(4))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(1,2)    ! phi- -c d
zgsffr(44,-11,4)=  -zgw*DCMPLX(parmas(4))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(1,3)    ! phi- -t d
zgsffr(44,-3,8)=  -zgw*DCMPLX(parmas(8))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(2,1)    ! phi- -u s 
zgsffr(44,-7,8)=  -zgw*DCMPLX(parmas(8))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(2,2)    ! phi- -c s 
zgsffr(44,-11,8)=  -zgw*DCMPLX(parmas(8))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(2,3)    ! phi- -t s 
zgsffr(44,-3,12)=-zgw*DCMPLX(parmas(12))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(3,1)    ! phi- -u b 
zgsffr(44,-7,12)=-zgw*DCMPLX(parmas(12))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(3,2)    ! phi- -c b 
zgsffr(44,-11,12)=-zgw*DCMPLX(parmas(12))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(3,3)    ! phi- -t b 
       
zgsffl(44,-3,4)=  zgw*DCMPLX(parmas(3))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(1,1)    ! phi- -u d
zgsffl(44,-3,8)=  zgw*DCMPLX(parmas(3))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(2,1)    ! phi- -u s
zgsffl(44,-3,12)=  zgw*DCMPLX(parmas(3))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(3,1)    ! phi- -u b
zgsffl(44,-7,4)=  zgw*DCMPLX(parmas(7))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(1,2)    ! phi- -c d
zgsffl(44,-7,8)=  zgw*DCMPLX(parmas(7))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(2,2)    ! phi- -c s
zgsffl(44,-7,12)=  zgw*DCMPLX(parmas(7))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(3,2)    ! phi- -c b
zgsffl(44,-11,4)=zgw*DCMPLX(parmas(11))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(1,3)    ! phi- -t d
zgsffl(44,-11,8)=zgw*DCMPLX(parmas(11))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(2,3)    ! phi- -t s
zgsffl(44,-11,12)=zgw*DCMPLX(parmas(11))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*zckm(3,3)    ! phi- -t b 
       
zgsffl(43,-2,1) =-zgw*DCMPLX(parmas(2))/(DCMPLX(root2)*zwm)*DCMPLX(iuni)  ! phi+ -e nu_e
zgsffl(43,-6,5) =-zgw*DCMPLX(parmas(6))/(DCMPLX(root2)*zwm)*DCMPLX(iuni)  ! phi+ -mu nu_mu 
zgsffl(43,-10,9)=-zgw*DCMPLX(parmas(10))/(DCMPLX(root2)*zwm)*DCMPLX(iuni) ! phi+ -tau nu_tau 
       
zgsffr(43,-2,1)=zero
zgsffr(43,-6,5)=zero
zgsffr(43,-10,9)=zero
       
zgsffl(43,-4,3)=  -zgw*DCMPLX(parmas(4))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(1,1)) ! phi+ -d u 
zgsffl(43,-4,7)=  -zgw*DCMPLX(parmas(4))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(2,1)) ! phi+ -d c
zgsffl(43,-4,11)=  -zgw*DCMPLX(parmas(4))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(3,1)) ! phi+ -d t
zgsffl(43,-8,3)=  -zgw*DCMPLX(parmas(8))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(1,2)) ! phi+ -s u 
zgsffl(43,-8,7)=  -zgw*DCMPLX(parmas(8))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(2,2)) ! phi+ -s c 
zgsffl(43,-8,11)=  -zgw*DCMPLX(parmas(8))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(3,2)) ! phi+ -s t 
zgsffl(43,-12,3)=-zgw*DCMPLX(parmas(12))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(1,3)) ! phi+ -b u
zgsffl(43,-12,8)=-zgw*DCMPLX(parmas(12))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(2,3)) ! phi+ -b c
zgsffl(43,-12,11)=-zgw*DCMPLX(parmas(12))/(DCMPLX(root2)*zwm)&
              *DCMPLX(iuni)*DCONJG(zckm(3,3)) ! phi+ -b t
       
zgsffr(43,-4,3)=   zgw*DCMPLX(parmas(3))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(1,1)) ! phi+ -d u
zgsffr(43,-8,3)=   zgw*DCMPLX(parmas(3))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(1,2)) ! phi+ -s u
zgsffr(43,-12,3)=   zgw*DCMPLX(parmas(3))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(1,3)) ! phi+ -b u
zgsffr(43,-4,7)=   zgw*DCMPLX(parmas(7))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(2,1)) ! phi+ -d c 
zgsffr(43,-8,7)=   zgw*DCMPLX(parmas(7))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(2,2)) ! phi+ -s c 
zgsffr(43,-12,7)=   zgw*DCMPLX(parmas(7))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(2,3)) ! phi+ -b c 
zgsffr(43,-4,11)= zgw*DCMPLX(parmas(11))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(3,1)) ! phi+ -d t 
zgsffr(43,-8,11)= zgw*DCMPLX(parmas(11))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(3,2)) ! phi+ -s t 
zgsffr(43,-12,11)= zgw*DCMPLX(parmas(11))/(DCMPLX(root2)*zwm)&
             *DCMPLX(iuni)*DCONJG(zckm(3,3)) ! phi+ -b t 
! S S S       
zgs3(41,41,41)=(-DCMPLX(dnou(3))*zgw/DCMPLX(dnou(2)))*DCMPLX(hm**2)/zwm*DCMPLX(ihi*chi) ! HHH
zgs3(41,42,42)=(-zgw/DCMPLX(dnou(2)))*DCMPLX(hm**2)/zwm*DCMPLX(iuni)            !HXX
zgs3(41,43,44)=(-zgw/DCMPLX(dnou(2)))*DCMPLX(hm**2)/zwm*DCMPLX(iuni)            !H phi+ phi- 
       
zgs3(42,42,41)=(-zgw/DCMPLX(dnou(2)))*DCMPLX(hm**2)/zwm*DCMPLX(iuni)            ! XXH 
       
zgs3(43,43,41)=(-zgw/DCMPLX(dnou(2)))*DCMPLX(hm**2)/zwm*DCMPLX(iuni)            ! phi+ phi+ H 
       
zgs3(44,44,41)=(-zgw/DCMPLX(dnou(2)))*DCMPLX(hm**2)/zwm*DCMPLX(iuni)            ! phi- phi- H 
! S S S S      
zgs4(41,41,41,41)=(-DCMPLX(dnou(3))/DCMPLX(dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(ihi)!HHHH
zgs4(41,41,42,42)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni)!HHXX
zgs4(41,41,43,44)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni) !HHphi+phi- 

zgs4(42,42,41,41)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni)   !XXHH
zgs4(42,42,42,42)=(-DCMPLX(dnou(3)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni)   !XXXX 
zgs4(42,42,43,44)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni) !XXphi-phi+ 

zgs4(43,43,41,41)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni) !phi+phi+HH
zgs4(43,43,42,42)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni) !F+F+XX 
zgs4(43,43,44,44)=(-DCMPLX(dnou(1)/dnou(2)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni) !F+F+F-F-

zgs4(44,44,41,41)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni) !F-F-HH
zgs4(44,44,42,42)=(-DCMPLX(dnou(1)/dnou(4)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni)!F-F-XX  
zgs4(44,44,43,43)=(-DCMPLX(dnou(1)/dnou(2)))*zgw**2*DCMPLX(hm**2)/zwm**2*DCMPLX(iuni)!F-F-F+F+
       
!iqed=0
IF(iqed.EQ.1)THEN
       
   zgv3(32:34,32:34,32:34)=zero
   zgv4(32:34,32:34,32:34,32:34)=zero
   zgvffl(32:34,-12:-1,1:12)=zero
   zgvffr(32:34,-12:-1,1:12)=zero
   zgvvs(31:34,31:34,41:44)=zero
   zgsvv(41:44,31:34,31:34)=zero
   zgvvss(31:34,31:34,41:44,41:44)=zero
   zgssvv(41:44,41:44,31:34,31:34)=zero
   zgvss(31:34,41:44,41:44)=zero
   zgssv(41:44,41:44,31:34)=zero
   zgsffl(41:44,-12:-1,1:12)=zero
   zgsffr(41:44,-12:-1,1:12)=zero
   zgs3(41:44,41:44,41:44)=zero
   zgs4(41:44,41:44,41:44,41:44)=zero
       
ENDIF
END SUBROUTINE Helac_FeynRule_SM
END MODULE Helac_SM_FeynRule
