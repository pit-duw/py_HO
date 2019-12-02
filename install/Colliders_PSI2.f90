MODULE Colliders_PSI_2
USE MC_VEGAS
USE Helac_Global
USE Constants
USE func_psi
USE Helac_SM_FeynRule
USE Structf_PDFs
USE Helac_master
USE Helac_Func_1
USE Cuts_Module
USE Colliders_PSI_1
USE FO_plot
USE QEDPS_interface
USE Decay_interface
USE reweight_xsec
IMPLICIT NONE
!INTEGER::varnum,lnum,snum,nunwei
!LOGICAL::ptQ,lunwei2
!REAL(KIND(1d0))::EBMUP1,EBMUP2
!INTEGER::NPRUP=0,lwmax
!INTEGER::nprint
!INTEGER,PARAMETER::maxprint=8
!INTEGER::emep_ISR_shower_old
SAVE
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Electron-Positron Collider
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EP_1FS(rslt,itmxn,ncalln)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)::rslt   ! rslt(1)=cross section ,rslt(2)=sigma,rslt(3)=chi2
INTEGER,INTENT(IN),OPTIONAL::itmxn,ncalln
REAL(KIND=DBL)::vfes,sd,chi2a
INTEGER::ncalm,nc,ii
CHARACTER(len=4),DIMENSION(20)::chchar
INTEGER::iday0,ihr0,imin0,isec0,i100th0,iday1,ihr1,imin1,isec1,i100th1,iyr0,iyr1,imon0,imon1
INTEGER::IDBMUP1,IDBMUP2,IDWTUP
CALL ReadElem_logic("ptdisQ",ptQ)
IF(ptQ.AND.emep_ISR_shower.NE.0)STOP "PT distribution is not avialabel when there is ISR shower"
CALL ReadElem_logic('unwgt',lunwei2)
CALL ReadElem_integer('preunw',nunwei)
lwmax=0
NPRN=-1
STOP "Delta function in 2 > 1 process without PDF"
IF(iranhel.EQ.0)THEN
   varnum = 1 ! -2
   lnum = 0
   snum = 0
ELSEIF(iranhel.EQ.1)THEN
        varnum=1
        lnum=0
        snum=0
        DO ii=1,nhad
           IF(ABS(iflh(ii)).LT.40)THEN
              varnum=varnum+1
           ENDIF
        ENDDO
ELSEIF(iranhel.EQ.2)THEN
        varnum=1
        lnum=0
        snum=0
        DO ii=1,nhad
           IF(ABS(iflh(ii)).LT.40)THEN
              varnum=varnum+1
              lnum=lnum+1
           ENDIF
           IF(QN1P1F(ii).OR.QN3PJF(ii))THEN
              varnum=varnum+1
              lnum=lnum+1
           ENDIF
        ENDDO
ELSEIF(iranhel.EQ.3)THEN
        varnum=1
        DO ii=1,nhad
           IF(ABS(iflh(ii)).LT.40)THEN
              varnum=varnum+1
              lnum=lnum+1
           ENDIF
           IF(QN1P1F(ii))THEN
              varnum=varnum+1
              lnum=lnum+1
           ENDIF
           IF(QN3S1F(ii))THEN
              varnum=varnum+1
              snum=snum+1
           ENDIF
           IF(QN3PJF(ii))THEN
              varnum=varnum+2
              lnum=lnum+1
              snum=snum+1
           ENDIF
        ENDDO
ENDIF
DO ii=1,varnum
        XL(ii)=0.0d0
        XU(ii)=1.0d0
ENDDO
IF(PRESENT(itmxn))THEN
   ITMX=itmxn
ELSE
   ITMX=5
END IF
ITMX=1
IF(PRESENT(ncalln))THEN
   ncalm=ncalln
ELSE
   ncalm=5000
END IF
chchar(1)="40K"
chchar(2)="80K"
chchar(3)="160K"
chchar(4)="320K"
chchar(5)="640K"
chchar(6)="1M"
chchar(7)="2M"
chchar(8)="4M"
chchar(9)="8M"
chchar(10)="16M"
chchar(11)='32M'
chchar(12)='64M'
chchar(13)='120M'
chchar(14)='240M'
chchar(15)='480M'
chchar(16)='960M'
chchar(17)='2G'
chchar(18)='4G'
chchar(19)='8G'
chchar(20)='16G'
! call vegas and result
nprint=10000
NCALL=20000
CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
IF(measurespeed)THEN
   NCALL=500
   CALL VEGAS(varnum,EP_1FS_fxn,vfes,sd,chi2a)
   CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
   iyr1=iyr1-iyr0
   imon1=imon1-imon0
   iday1=iday1-iday0
   ihr1=ihr1-ihr0
   imin1=imin1-imin0
   isec1=isec1-isec0
   i100th1=i100th1-i100th0
   CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
   NCALL=20000
ENDIF
CALL VEGAS(varnum,EP_1FS_fxn,vfes,sd,chi2a)
WRITE(*,*)' '
WRITE(*,*)"====================NCALL=20K==========================="
WRITE(*,*)" "
CALL Helac_mtime()
WRITE(*,*)" "
ii=1
WRITE(*,*)"ITERATION ",ii,":"
WRITE(*,*)vfes,"+\-",sd
WRITE(*,*)"precision:",sd/vfes
DO ii=2,10
        IF(plot_output)CALL initplot
        CALL VEGAS(varnum,EP_1FS_fxn,vfes,sd,chi2a,1)
        IF(plot_output)CALL plotout
        WRITE(*,*)"ITERATION ",ii,":"
        WRITE(*,*)vfes,"+\-",sd
        WRITE(*,*)"precision:",sd/vfes
ENDDO
CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
iday1=iday1-iday0
ihr1=ihr1-ihr0
imin1=imin1-imin0
isec1=isec1-isec0
i100th1=i100th1-i100th0
CALL Vegas_speed(10*NCALL,iday1,ihr1,imin1,isec1,i100th1)
WRITE(*,*)' '
WRITE(*,*)' '
ii=1
DO
  nc=2*NCALL
  IF(nc.GT.ncalm)EXIT
  IF(2*nc.GT.ncalm.AND.lunwei2.AND.lwmax.EQ.0)THEN
     lwmax=1
  ENDIF
  NCALL=nc
  IF(NCALL/maxprint.GT.nprint)nprint=NCALL/maxprint
  !CALL GETDAT(iyr0,imon0,iday0)                                                
  !CALL GETTIM(ihr0,imin0,isec0,i100th0)                                        
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_1FS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd
  WRITE(*,*)"precision:",sd/vfes
  !CALL GETDAT(iyr1,imon1,iday1)                                                
  !CALL GETTIM(ihr1,imin1,isec1,i100th1)                                        
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
  ii=ii+1
ENDDO
IF(lunwei2.AND.lwmax.EQ.1)THEN
   lwmax=2
   WRITE(*,*)"START UNWEIGHTING"
   NCALL=nc/2
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii-1)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_1FS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd
  WRITE(*,*)"precision:",sd/vfes
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
ENDIF
SELECT CASE(COLL_TYPE)
CASE(1)
   IDBMUP1=2212
   IDBMUP2=2212
CASE(2)
   IDBMUP1=2212
   IDBMUP2=-2212
CASE DEFAULT
   IDBMUP1=11
   IDBMUP2=-11
END SELECT
IF(lunwei2)THEN
   IDWTUP=3
ELSE
   IDWTUP=1
ENDIF
NPRUP=NPRUP+1
WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&
     iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
!    &             0,       0,   10042,   10042,IDWTUP,NPRUP    
! XSECUP,XERRUP,XMAXUP,LPRUP 
WRITE(200,5200) vfes*10d0**3,sd*10d0**3,1d0, 81
rslt(1)=vfes
rslt(2)=sd
rslt(3)=chi2a
IF(lunwei2)PRINT *,"number of events",Nevents 
RETURN
5100 FORMAT(1P,2I8,2E14.6,6I8)
5200 FORMAT(1P,3E20.10,I6)
END SUBROUTINE EP_1FS

FUNCTION EP_1FS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt ! it is the weight from VEGAS, 
! which will be used to generate unweighted events, wgt=VTOT*Vi/Ni=Vi/Ni since VTOT=1
REAL(KIND=DBL)::EP_1FS_fxn,sqs,m1r
REAL(KIND=DBL)::sq,m12
REAL(KIND=DBL)::temp1,ymax,ymin,Jy,y,w1
REAL(KIND=DBL)::wme,wme0,ptp,pqp,pi,wsf,recmax=0,recmin=0,y1cup,y1clow
REAL(KIND=DBL)::exp1,exp2,ycollcm
REAL(KIND=DBL),DIMENSION(4)::pmom
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
INTEGER::init=0,i,kk,j,istat,nnn=0,nnntot=0,pdfnumpdf,icut,nwmax,nwarn,nwri,nwri_tot
SAVE init,sqs,m1r,m12,sq,temp1,ymin,ymax,pi,lr,sr,nnntot,ycollcm
SAVE nnn,recmax,recmin,y1cup,y1clow,Jy,nwmax,nwarn,nwri,nwri_tot
IF(init.EQ.0)THEN
    pi=DACOS(-1d0)
    wjac=1
    nwmax=0
    nwarn=0
    nwri=0
    CALL ReadElem_integer('unwevt',nwri_tot)
    CALL ReadElem_real("energy_beam1",sqs)
    ebeam(1)=ABS(sqs)
    CALL ReadElem_real("energy_beam2",sqs)
    ebeam(2)=ABS(sqs)
    IF(ABS(ebeam(1)-ebeam(2))/MAX(ebeam(1)+ebeam(2),1d-17).LT.1d-8)THEN
       absrap=absrap
       labeqcoll=.TRUE.
    ELSE
       absrap=.FALSE.
       labeqcoll=.FALSE.
    ENDIF
    sqs=2d0*DSQRT(ebeam(1)*ebeam(2)) ! we always neglect the mass of initial states
    EBMUP1=ebeam(1)
    EBMUP2=ebeam(2)
    ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
    CALL ReadElem_real("maxy1c",y1cup)
    CALL ReadElem_real("miny1c",y1clow)
    CALL ReadElem_integer('pdf',pdfnumpdf)
    ptc(3:20)=0                   ! Debug
    maxptc(3:20)=-1d0
    drc(3:20,3:20)=0
    etac(3:20)=20
    ec(3:20)=0
    c1(3:20)=1
    c2(3:20)=1
    cc(3:20,3:20)=1
    yycut(3:20)=1d9
    IF(absrap)THEN
       yycutlow(3:20)=0d0
    ELSE
       yycutlow(3:20)=-1d9
    ENDIF
    gmas(3:20,3:20)=0             ! Debug
    gbeammass(1:2,3:20)=0d0

    parmas(0)=0
    parwid(0)=-1
    iPDFSUP1=pdfnumpdf
    IF(pdfnumpdf.EQ.0.OR.COLL_TYPE.EQ.3)THEN
        istruc=0
    ELSE
        istruc=1
    ENDIF
    IF(istruc.EQ.1)THEN      
       iPDFGUP1=0 
       iPDFGUP2=0 
       iPDFSUP1=pdfnumpdf 
       iPDFSUP2=pdfnumpdf 
    ELSE
       iPDFGUP1=-1
       iPDFGUP2=-1
       iPDFSUP1=-1
       iPDFSUP2=-1
       reweight_pdf=.FALSE.
       ho_npdf=0
    ENDIF
    istat=0
    xp1=1
    xp2=1
    IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
        CALL readcuts_HC
        y1cup=MIN(y1cup,yycut(hadron2parton(3)))
        y1clow=MAX(y1clow,yycutlow(hadron2parton(3)))
        IF(absrap)y1cup = ABS(y1cup)
        IF(absrap)y1clow = ABS(y1clow)
    ELSEIF(COLL_TYPE.EQ.3)THEN
        CALL readcuts_epem
        IF(absrap)y1cup=ABS(y1cup)
        IF(absrap)y1clow=ABS(y1clow)
    ENDIF
    IF(NDecayChains.GT.0)CALL readcuts_Decay
    IF(lnum.NE.0)THEN
        ALLOCATE(lr(lnum))
    ENDIF
    IF(snum.NE.0)THEN
        ALLOCATE(sr(snum))
    ENDIF
    IF(istat.NE.0)THEN
        PRINT *,"Wrong in allocation of array in EP_1FS_fxn !STOP !"
        STOP
    ENDIF
    kk=3
    IF(ABS(iflh(3)).GT.100)THEN
        m1r=parmas(ifl(kk))+parmas(ifl(kk+1))
        kk=kk+2
    ELSE
        m1r=parmas(ifl(kk))
        kk=kk+1
    ENDIF
    sq=sqs*sqs
    m12=m1r*m1r
    IF(m12.LE.0d0)THEN
        PRINT *,"ERROR: 2 > 1 with massless final state ?"
        PRINT *,"m = ", m1r
        PRINT *,"sqrtS = ",sqs 
        STOP
    ENDIF
    IF(absrap)THEN
       IF(DLOG(sq/m12)/2d0.LE.y1clow)THEN
          PRINT *,"ERROR: outside the rapidity cuts !"
          PRINT *,"m = ",m1r
          PRINT *,"sqrtS = ",sqs
          PRINT *,"yclow = ",y1clow
          STOP
       ENDIF
    ELSE
       IF(DLOG(sq/m12)/2d0.LE.y1clow-ycollcm.OR.&
            -DLOG(sq/m12)/2d0.GE.y1cup-ycollcm)THEN
          PRINT *,"ERROR: outside the rapidity cuts !"
          PRINT *,"m = ",m1r
          PRINT *,"sqrtS = ",sqs
          PRINT *,"ycup = ",y1cup
          PRINT *,"yclow = ",y1clow
          PRINT *,"ycollcm = ",ycollcm
          STOP
       ENDIF
    ENDIF
! ============================================================
! tau=x1*x2, y=1/2*log(x1/x2) => dtau dy = dx1 dx2
! 2>1 body phase space
! pi/2/M^3 delta(M-sqrt(S*tau)) dtau dy = pi/2/M^3 * 2M/S dy
! x1 = sqrt(M^2/S) Exp(y), x2 = sqrt(M^2/S) Exp(-y)
! => |y|<1/2*log(S/M^2)
! ============================================================
    temp1=pi/m12/sq ! phase space factor
    ymax=MIN(y1cup-ycollcm,DLOG(sq/m12)/2d0) ! upper limit of y
    IF(absrap)THEN
       ymin=MAX(y1clow,0d0) ! lower limit of y
    ELSE
       ymin=MAX(y1clow-ycollcm,-DLOG(sq/m12)/2d0)
    ENDIF
    Jy=(ymax-ymin)*temp1 ! jacobian* phase space factor
    IF(absrap)Jy=2d0*Jy
    IF(ptQ)THEN
        PRINT *, "WARNING:Pt = 0 in 2 > 1 processes ! Calculating cross section instead !"
    ENDIF
    init=1
ENDIF
nnntot=nnntot+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The physical region of the integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! total cross section
y=(2d0*x(1)-1d0)*(ymax-ymin)+SIGN(1d0,x(1)-0.5d0)*ymin
xp1=DSQRT(m12/sq)*DEXP(y) ! xp1,xp2 are the same in collision frame and lab frame
xp2=DSQRT(m12/sq)*DEXP(-y)
exp1=xp1*ebeam(1)
exp2=xp2*ebeam(2)
ehat=DSQRT(m12)
! Generate the momenta of external legs
Phegas_pmom(1,1:2)=0
Phegas_pmom(1,3)=exp1
Phegas_pmom(1,4)=exp1
Phegas_pmom(1:2,5)=0
Phegas_pmom(2,1:2)=0
Phegas_pmom(2,3)=-exp2
Phegas_pmom(2,4)=exp2
kk=3
pmom(1)=0
pmom(2)=0
pmom(3)=(exp1-exp2)
pmom(4)=(exp1+exp2)
IF(ABS(iflh(3)).GT.100)THEN
    Phegas_pmom(kk,1:4)=parmas(ifl(kk))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
    Phegas_pmom(kk,5)=parmas(ifl(kk))
    Phegas_pmom(kk+1,1:4)=parmas(ifl(kk+1))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
    Phegas_pmom(kk+1,5)=parmas(ifl(kk+1))
    kk=kk+2
ELSE
    Phegas_pmom(kk,1:4)=pmom(1:4)
    Phegas_pmom(kk,5)=parmas(ifl(kk))
    kk=kk+1
ENDIF
! boost from lab frame to cms
pmom(3)=-(exp1-exp2)
pmom(4)=(exp1+exp2)
pmom(1:2)=0
DO j=1,n
    IF(imode.EQ.0.AND.NDecayChains.EQ.0)CALL Boostl2(ehat,pmom,Phegas_pmom(j,1:4))
    ptp=DSQRT(Phegas_pmom(j,1)**2+Phegas_pmom(j,2)**2)
    pqp=DSQRT(ptp**2+Phegas_pmom(j,3)**2)
    zq(j,1)= DCMPLX(Phegas_pmom(j,4)+Phegas_pmom(j,3),Phegas_pmom(j,3))
    zq(j,2)= DCMPLX(Phegas_pmom(j,4)-Phegas_pmom(j,3),ptp)
    zq(j,3)= DCMPLX(Phegas_pmom(j,1), Phegas_pmom(j,2))
    zq(j,4)= DCMPLX(Phegas_pmom(j,1),-Phegas_pmom(j,2))
    zq(j,5)= DCMPLX(parmas(ifl(j)),pqp )
ENDDO
icut=1
IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
    CALL Cuts_HC(icut)
ELSEIF(COLL_TYPE.EQ.3)THEN
    CALL Cuts_epem(icut)
ELSE
    icut=1
ENDIF

IF(icut.EQ.0)THEN
    EP_1FS_fxn=0d0
    RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
kk=2
DO j=1,lnum
    lr(j)=2*pi*x(kk)
    kk=kk+1
ENDDO
DO j=1,snum
    sr(j)=2*pi*x(kk)
    kk=kk+1
ENDDO
CALL SetRandom(lr,sr)
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1)THEN
   reweight_scale_phase=.TRUE.
   CALL Helac_FeynRule_SM()
   CALL Helac_master_f(wme0)
   IF(nnn.EQ.0)THEN
      alphas_power=NINT(DLOG10(wme/wme0))
   ELSE
      IF(NINT(DLOG10(wme/wme0)).NE.alphas_power)THEN
         WRITE(*,*)"WARNING:The process is a mixed order process.Reweighting scale may be wrong !"
      ENDIF
   ENDIF
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   EP_1FS_fxn=0d0
   RETURN
ENDIF
CALL strf_pdf(wsf)
IF(ABS(wsf).LE.0d0)THEN
   EP_1FS_fxn=0d0
   RETURN
ENDIF
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
EP_1FS_fxn=wme*wsf*LDMEwt*Jy*3.8938573d5
! only unweight events are generated -> LHE files                                                                                                                                                         
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------
   w1=EP_1FS_fxn*wgt ! multiply VEGAS weight
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=EP_1FS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.EP_1FS_fxn)recmax=EP_1FS_fxn
IF(recmin.GT.EP_1FS_fxn)recmin=EP_1FS_fxn
IF(MOD(nnn,nprint).EQ.0)THEN
PRINT *,"max=",recmax
PRINT *,"min=",recmin
IF(lunwei2.AND.lwmax.GT.1)THEN
   PRINT *, "      n_event,    n_pass,    n_total"
   PRINT *, nwri,nnn,nnntot
ELSE
   PRINT *, "      n_pass,     n_total"
   PRINT *, nnn,nnntot
ENDIF
ENDIF
END FUNCTION EP_1FS_fxn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Cross Section of Electron Positron -> 2 Final States
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EP_2FS(rslt,itmxn,ncalln)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)::rslt   ! rslt(1)=cross section ,rslt(2)=sigma,rslt(3)=chi2
INTEGER,INTENT(IN),OPTIONAL::itmxn,ncalln
REAL(KIND=DBL)::vfes,sd,chi2a
INTEGER::ncalm,nc,ii
CHARACTER(len=4),DIMENSION(20)::chchar
INTEGER::iday0,ihr0,imin0,isec0,i100th0,iday1,ihr1,imin1,isec1,i100th1,iyr0,iyr1,imon0,imon1
INTEGER::IDBMUP1,IDBMUP2,IDWTUP
CALL ReadElem_logic("ptdisQ",ptQ)
IF(ptQ.AND.emep_ISR_shower.NE.0)STOP "PT distribution is not avialabel when there is ISR shower"
CALL ReadElem_logic('unwgt',lunwei2)
CALL ReadElem_integer('preunw',nunwei)
lwmax=0
NPRN=-1
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
	IF(ptQ)THEN
		varnum=1
	ELSE
		varnum=2
	ENDIF
        IF(lunwei2)varnum=varnum+1 ! the trivial phi should be used
	lnum=0
	snum=0
ELSEIF(iranhel.EQ.1)THEN
	IF(ptQ.AND..NOT.lunwei2)THEN
		varnum=1
	ELSE
		varnum=2
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
	DO ii=1,nhad
		IF(ABS(iflh(ii)).LT.40)THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
	ENDDO
ELSEIF(iranhel.EQ.2)THEN
	IF(ptQ)THEN
		varnum=1
	ELSE
		varnum=2
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
	DO ii=1,nhad
		IF(ABS(iflh(ii)).LT.40)THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
		IF(QN1P1F(ii).OR.QN3PJF(ii))THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
	ENDDO
ELSEIF(iranhel.EQ.3)THEN
	IF(ptQ)THEN
		varnum=1
	ELSE
		varnum=2
	ENDIF
        IF(lunwei2)varnum=varnum+1
	DO ii=1,nhad
		IF(ABS(iflh(ii)).LT.40)THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
		IF(QN1P1F(ii))THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
		IF(QN3S1F(ii))THEN
			varnum=varnum+1
			snum=snum+1
		ENDIF
		IF(QN3PJF(ii))THEN
			varnum=varnum+2
			lnum=lnum+1
			snum=snum+1
		ENDIF
	ENDDO
ENDIF
ELSE
   IF(ptQ)THEN
      varnum=1
   ELSE
      varnum=2
   ENDIF
   IF(lunwei2)varnum=varnum+1 ! the trivial phi should be used                                               
   lnum=0
   snum=0
   varnum=varnum+1+NDecayIflh
ENDIF
IF(varnum.LE.1)THEN
   WRITE(*,*)"WARNING: it is trivial for the remaining phi distribution in 2 > 1 at e+e-"
ENDIF
DO ii=1,varnum
	XL(ii)=0.0d0
	XU(ii)=1.0d0
ENDDO
IF(PRESENT(itmxn))THEN
   ITMX=itmxn
ELSE
   ITMX=5
END IF
ITMX=1
IF(PRESENT(ncalln))THEN
   ncalm=ncalln
ELSE
   ncalm=5000
END IF
chchar(1)="40K"
chchar(2)="80K"
chchar(3)="160K"
chchar(4)="320K"
chchar(5)="640K"
chchar(6)="1M"
chchar(7)="2M"
chchar(8)="4M"
chchar(9)="8M"
chchar(10)="16M"
chchar(11)='32M'
chchar(12)='64M'
chchar(13)='120M'
chchar(14)='240M'
chchar(15)='480M'
chchar(16)='960M'
chchar(17)='2G'
chchar(18)='4G'
chchar(19)='8G'
chchar(20)='16G'
! call vegas and result
nprint=10000
NCALL=20000
!CALL GETDAT(iyr0,imon0,iday0)
!CALL GETTIM(ihr0,imin0,isec0,i100th0)
CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
IF(measurespeed)THEN
   NCALL=500
   CALL VEGAS(varnum,EP_2FS_fxn,vfes,sd,chi2a)
   CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
   iyr1=iyr1-iyr0
   imon1=imon1-imon0
   iday1=iday1-iday0
   ihr1=ihr1-ihr0
   imin1=imin1-imin0
   isec1=isec1-isec0
   i100th1=i100th1-i100th0
   CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
   NCALL=20000
ENDIF
CALL VEGAS(varnum,EP_2FS_fxn,vfes,sd,chi2a)
WRITE(*,*)' '
WRITE(*,*)"====================NCALL=20K==========================="
WRITE(*,*)" "
CALL Helac_mtime()
WRITE(*,*)" "
ii=1
WRITE(*,*)"ITERATION ",ii,":"
WRITE(*,*)vfes,"+\-",sd
WRITE(*,*)"precision:",sd/vfes
DO ii=2,10
        IF(plot_output)CALL initplot
	CALL VEGAS(varnum,EP_2FS_fxn,vfes,sd,chi2a,1)
        IF(plot_output)CALL plotout
	WRITE(*,*)"ITERATION ",ii,":"
	WRITE(*,*)vfes,"+\-",sd
	WRITE(*,*)"precision:",sd/vfes
ENDDO
!CALL GETDAT(iyr1,imon1,iday1)
!CALL GETTIM(ihr1,imin1,isec1,i100th1)
CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
iday1=iday1-iday0
ihr1=ihr1-ihr0
imin1=imin1-imin0
isec1=isec1-isec0
i100th1=i100th1-i100th0
CALL Vegas_speed(10*NCALL,iday1,ihr1,imin1,isec1,i100th1)
WRITE(*,*)' '
WRITE(*,*)' '
!CALL VEGAS(3,EP_2FS_fxn,vfes,sd,chi2a,1)
!ITMX=1
ii=1
DO
  nc=2*NCALL
  IF(nc.GT.ncalm)EXIT
  IF(2*nc.GT.ncalm.AND.lunwei2.AND.lwmax.EQ.0)THEN
     lwmax=1
  ENDIF
  NCALL=nc
  IF(NCALL/maxprint.GT.nprint)nprint=NCALL/maxprint
  !CALL GETDAT(iyr0,imon0,iday0)
  !CALL GETTIM(ihr0,imin0,isec0,i100th0)
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  !IF(2*nc.GT.ncalm.AND.emep_ISR_shower_old.EQ.1&
  !     .AND..NOT.lunwei2)start_shower=.TRUE.
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_2FS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd
  WRITE(*,*)"precision:",sd/vfes
  !CALL GETDAT(iyr1,imon1,iday1)
  !CALL GETTIM(ihr1,imin1,isec1,i100th1)
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
  ii=ii+1
ENDDO
IF(lunwei2.AND.lwmax.EQ.1)THEN
   lwmax=2
   WRITE(*,*)"START UNWEIGHTING"
   NCALL=nc/2
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii-1)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  !IF(emep_ISR_shower_old.EQ.1)start_shower=.TRUE.
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_2FS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd
  WRITE(*,*)"precision:",sd/vfes
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
ENDIF
SELECT CASE(COLL_TYPE)
CASE(1)
   IDBMUP1=2212
   IDBMUP2=2212
CASE(2)
   IDBMUP1=2212
   IDBMUP2=-2212
CASE DEFAULT
   IDBMUP1=11
   IDBMUP2=-11
END SELECT
IF(lunwei2)THEN
   IDWTUP=3
ELSE
   IDWTUP=1
ENDIF
NPRUP=NPRUP+1
WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&
     iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
!    &             0,       0,   10042,   10042,IDWTUP,NPRUP  
! XSECUP,XERRUP,XMAXUP,LPRUP     
WRITE(200,5200) vfes*10d0**3,sd*10d0**3,1d0, 81
rslt(1)=vfes
rslt(2)=sd
rslt(3)=chi2a
IF(lunwei2)PRINT *,"number of events",Nevents
IF(emep_ISR_shower.EQ.1)CALL QEDPSEND
RETURN
5100 FORMAT(1P,2I8,2E14.6,6I8)
5200 FORMAT(1P,3E20.10,I6)
END SUBROUTINE EP_2FS

FUNCTION EP_2FS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt
REAL(KIND=DBL)::EP_2FS_fxn,sqs,m1r,m2r,w1
REAL(KIND=DBL)::pt,y1,y2,mt1,mt2,sq,m12,m22,phi
REAL(KIND=DBL)::temp1,temp2,temp21,temp3,temp4,Jpt2,Jx1x2
REAL(KIND=DBL)::wme,wme0,wgtbr,ptp,pqp,pi,wsf,recmax=0,recmin=0,pt1c,y1cup,y1clow
REAL(KIND=DBL)::exp1,exp2,ycollcm,y1cup_old,y1clow_old,pt1c_old
REAL(KIND=DBL),DIMENSION(4)::pmom
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
REAL(KIND=DBL),DIMENSION(3:20)::ptc_old,yycut_old,yycutlow_old
INTEGER::init=0,kk,j,istat,nnn=0,nnntot=0,pdfnumpdf,icut,nwri,nwri_tot,nwarn,nwmax
SAVE init,sqs,m1r,m2r,m12,m22,sq,temp1,pi,lr,sr,pt,nnntot,ycollcm
SAVE nnn,recmax,recmin,pt1c,y1cup,y1clow,nwri,nwri_tot,nwarn,nwmax
SAVE ptc_old,yycut_old,yycutlow_old,y1cup_old,y1clow_old,pt1c_old
IF(init.EQ.0)THEN
	pi=DACOS(-1d0)
	wjac=1
        nwmax=0
        nwarn=0
        nwri=0
        CALL ReadElem_integer('unwevt',nwri_tot)
        CALL ReadElem_real("energy_beam1",sqs)
        ebeam(1)=ABS(sqs)
        CALL ReadElem_real("energy_beam2",sqs)
        ebeam(2)=ABS(sqs)
        IF(ABS(ebeam(1)-ebeam(2))/MAX(ebeam(1)+ebeam(2),1d-17).LT.1d-8)THEN
           absrap=absrap
           labeqcoll=.TRUE.
        ELSE
           IF(emep_ISR_shower.EQ.0)THEN
              absrap=.FALSE.
              labeqcoll=.FALSE.
           ELSE
              absrap=.FALSE.
              labeqcoll=.TRUE. ! run in c.m. frame
           ENDIF
        ENDIF
        sqs=2d0*DSQRT(ebeam(1)*ebeam(2)) ! we always neglect the mass of initial states
        EBMUP1=ebeam(1)
        EBMUP2=ebeam(2)
        IF(emep_ISR_shower.EQ.0)THEN
           ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
        ELSE
           ycollcm=0d0
        ENDIF
	CALL ReadElem_real("minpt1c",pt1c)
	CALL ReadElem_real("maxy1c",y1cup)
	CALL ReadElem_real("miny1c",y1clow)
	CALL ReadElem_integer('pdf',pdfnumpdf)
	ptc(3:20)=0                   ! Debug
        maxptc(3:20)=-1d0
	drc(3:20,3:20)=0
	etac(3:20)=20
	ec(3:20)=0
	c1(3:20)=1
	c2(3:20)=1
	cc(3:20,3:20)=1
        yycut(3:20)=1d9
        IF(absrap)THEN
           yycutlow(3:20)=0d0
        ELSE
           yycutlow(3:20)=-1d9
        ENDIF
	gmas(3:20,3:20)=0             ! Debug
        gbeammass(1:2,3:20)=0

        parmas(0)=0
        parwid(0)=-1
	iPDFSUP1=pdfnumpdf
	IF(pdfnumpdf.EQ.0.OR.COLL_TYPE.EQ.3)THEN
		istruc=0
	ELSE
		istruc=1
	ENDIF
        IF(istruc.EQ.1)THEN
                iPDFGUP1=0
                iPDFGUP2=0
                iPDFSUP1=pdfnumpdf
                iPDFSUP2=pdfnumpdf
        ELSE
                iPDFGUP1=-1
                iPDFGUP2=-1
                iPDFSUP1=-1
                iPDFSUP2=-1
                reweight_pdf=.FALSE.
                ho_npdf=0
        ENDIF
	istat=0
	xp1=1
        xp2=1
	IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
		CALL readcuts_HC
                y1cup=MIN(y1cup,yycut(hadron2parton(3)))
                y1clow=MAX(y1clow,yycutlow(hadron2parton(3)))
                IF(absrap)y1cup=ABS(y1cup)
                IF(absrap)y1clow=ABS(y1clow)
	ELSEIF(COLL_TYPE.EQ.3)THEN
		CALL readcuts_epem
                IF(absrap)y1cup=ABS(y1cup)
                IF(absrap)y1clow=ABS(y1clow)
	ENDIF
        IF(NDecayChains.GT.0)CALL readcuts_Decay
	IF(lnum.NE.0)THEN
		ALLOCATE(lr(lnum))
	ENDIF
	IF(snum.NE.0)THEN
		ALLOCATE(sr(snum))
	ENDIF
	IF(istat.NE.0)THEN
		PRINT *,"Wrong in allocation of array in EP_2FS_fxn !STOP !"
		STOP
	ENDIF
	kk=3
	IF(ABS(iflh(3)).GT.100)THEN
		m1r=parmas(ifl(kk))+parmas(ifl(kk+1))
		kk=kk+2
	ELSE
		m1r=parmas(ifl(kk))
		kk=kk+1
	ENDIF
	IF(ABS(iflh(4)).GT.100)THEN
		m2r=parmas(ifl(kk))+parmas(ifl(kk+1))
		kk=kk+2
	ELSE
		m2r=parmas(ifl(kk))
		kk=kk+1
	ENDIF
	sq=sqs*sqs
	m12=m1r*m1r
	m22=m2r*m2r
	temp1=DSQRT(lamda_func(sq,m12,m22))/(2d0*sqs)
        IF(emep_ISR_shower.EQ.1)CALL QEDPSINIT
        ptc_old(3:20)=ptc(3:20)
        yycut_old(3:20)=yycut(3:20)
        yycutlow_old(3:20)=yycutlow(3:20)
        pt1c_old=pt1c
        y1cup_old=y1cup
        y1clow_old=y1clow
        !emep_ISR_shower_old=emep_ISR_shower
	IF(ptQ)THEN
		CALL ReadElem_real("Pt1",pt)
	ENDIF
	init=1
ENDIF
nnntot=nnntot+1
IF(emep_ISR_shower.EQ.1)THEN
   !IF(start_shower)THEN
   CALL QEDPS_ISR_SHOWER
   !   emep_ISR_shower=emep_ISR_shower_old
   !ELSE
   !   Q2OUT_QEDPS=Q2MAX_QEDPS
   !   emep_ISR_shower=0
   !ENDIF
   IF(Q2OUT_QEDPS.LT.(m1r+m2r)**2)THEN
      EP_2FS_fxn=0d0
      RETURN
   ENDIF
   ebeam(1)=DSQRT(Q2OUT_QEDPS)/2d0
   ebeam(2)=ebeam(1)
   sq=4d0*ebeam(1)*ebeam(2)
   sqs=DSQRT(sq)
   temp1=DSQRT(lamda_func(sq,m12,m22))/(2d0*sqs)
   ptc(3:20)=0d0 ! a small nonzero value to aviod IR divergence
   maxptc(3:20)=-1d0
   pt1c=0d0
   yycut(3:20)=1d9
   y1cup=1d9
   IF(absrap)THEN
      yycutlow(3:20)=0d0
      y1clow=0d0
   ELSE
      yycutlow(3:20)=-1d9
      y1clow=-1d9
   ENDIF
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The physical region of the integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(.NOT.ptQ)THEN
! total cross section
	pt=(temp1-pt1c)*x(1)+pt1c
	mt1=DSQRT(pt*pt+m1r*m1r)
	mt2=DSQRT(pt*pt+m2r*m2r)
        temp2=ACosh_p((sq+m12-m22)/(2d0*sqs*mt1))
        IF(temp2.LT.y1clow-ycollcm)THEN
           EP_2FS_fxn=0d0
           RETURN
        ENDIF
        y1=SIGN(1d0,x(2)-0.5d0)*temp2
! in this case -DLOG((-DEXP(-y1)*mt1 +sqs)/mt2)=DLOG((-DEXP(+y1)*mt1 +sqs)/mt2)
	temp3=DLOG((-DEXP(-y1)*mt1 +sqs)/mt2)
	y2= -temp3
!phi=x(4)*2*pi
!Jacobi Determinations
	Jpt2=2d0*pt*(temp1-pt1c)
        temp4=ABS(DSINH(y1-y2))
! a factor of 2 is comming from SIGN(1d0,x(2)-0.5d0)=(2*theta(x(2))-1)
        Jx1x2=sq/(mt1*mt2*temp4)
ELSE
! pt distribution
	mt1=DSQRT(pt*pt+m1r*m1r)
	mt2=DSQRT(pt*pt+m2r*m2r)
	temp2=ACosh_p((sq+m12-m22)/(2d0*sqs*mt1))
	IF(temp2.LT.y1clow-ycollcm)THEN
              EP_2FS_fxn=0d0
              RETURN
	ENDIF
	y1=SIGN(1d0,x(1)-0.5d0)*temp2
	temp3=DLOG((-DEXP(-y1)*mt1 +sqs)/mt2)
	y2= -temp3
!phi=x(4)*2*pi
!Jacobi Determinations
	Jpt2=2d0*pt
! dy1 dy2 -> dx1 dx2
        temp4=ABS(DSINH(y1-y2))
! a factor of 2 is comming from SIGN(1d0,x(2)-0.5d0)=(2*theta(x(2))-1) 
        Jx1x2=sq/(mt1*mt2*temp4)
ENDIF
! The substitution of original variations
xp1=1d0
xp2=1d0
exp1=xp1*ebeam(1)
exp2=xp2*ebeam(2)
ehat=DSQRT(xp1*xp2*sq)
! Generate the momenta of external legs
Phegas_pmom(1,1:2)=0
Phegas_pmom(1,3)=xp1*sqs/2d0
Phegas_pmom(1,4)=xp1*sqs/2d0
Phegas_pmom(1:2,5)=0
Phegas_pmom(2,1:2)=0
Phegas_pmom(2,3)=-xp2*sqs/2d0
Phegas_pmom(2,4)=xp2*sqs/2d0
! we choose phi=0
IF(lunwei2)THEN
   IF(ptQ)THEN
      phi=2*pi*x(2)
   ELSE
      phi=2*pi*x(3)
   ENDIF
ELSE
   phi=0d0
ENDIF
kk=3
pmom(1)=pt*DCOS(phi)
pmom(2)=pt*DSIN(phi)
pmom(3)=mt1*DSINH(y1)
pmom(4)=mt1*DCOSH(y1)
IF(ABS(iflh(3)).GT.100)THEN
	Phegas_pmom(kk,1:4)=parmas(ifl(kk))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	Phegas_pmom(kk+1,1:4)=parmas(ifl(kk+1))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk+1,5)=parmas(ifl(kk+1))
	kk=kk+2
ELSE
	Phegas_pmom(kk,1:4)=pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	kk=kk+1
ENDIF
pmom(1)=-pt*DCOS(phi)
pmom(2)=-pt*DSIN(phi)
pmom(3)=mt2*DSINH(y2)
pmom(4)=mt2*DCOSH(y2)
IF(ABS(iflh(4)).GT.100)THEN
	Phegas_pmom(kk,1:4)=parmas(ifl(kk))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	Phegas_pmom(kk+1,1:4)=parmas(ifl(kk+1))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk+1,5)=parmas(ifl(kk+1))
	kk=kk+2
ELSE
	Phegas_pmom(kk,1:4)=pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	kk=kk+1
ENDIF
IF(imode.EQ.0.AND.NDecayChains.EQ.0)THEN
   ! boost from collision frame to cms
   pmom(3)=-(xp1-xp2)*sqs/2d0
   pmom(4)=(xp1+xp2)*sqs/2d0
   pmom(1:2)=0
ELSE
   pmom(3)=ebeam(1)-ebeam(2)
   pmom(4)=ebeam(1)+ebeam(2)
   pmom(1:2)=0
ENDIF
DO j=1,n
	IF(imode.EQ.0.AND.NDecayChains.EQ.0)THEN
           CALL Boostl2(ehat,pmom,Phegas_pmom(j,1:4))
        ELSEIF(.NOT.labeqcoll)THEN
           CALL Boostl2(sqs,pmom,Phegas_pmom(j,1:4))
        ENDIF
	ptp=DSQRT(Phegas_pmom(j,1)**2+Phegas_pmom(j,2)**2)
	pqp=DSQRT(ptp**2+Phegas_pmom(j,3)**2)
	zq(j,1)= DCMPLX(Phegas_pmom(j,4)+Phegas_pmom(j,3),Phegas_pmom(j,3))
	zq(j,2)= DCMPLX(Phegas_pmom(j,4)-Phegas_pmom(j,3),ptp)
	zq(j,3)= DCMPLX(Phegas_pmom(j,1), Phegas_pmom(j,2))
	zq(j,4)= DCMPLX(Phegas_pmom(j,1),-Phegas_pmom(j,2))
	zq(j,5)= DCMPLX(parmas(ifl(j)),pqp )
ENDDO
! recover the cuts
IF(emep_ISR_shower.EQ.1)THEN
   pt1c=pt1c_old
   ptc(3:20)=ptc_old(3:20)
   yycut(3:20)=yycut_old(3:20)
   yycutlow(3:20)=yycutlow_old(3:20)
   y1cup=y1cup_old
   y1clow=y1clow_old
ENDIF
icut=1
! cuts on the hard particles
IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
	CALL Cuts_HC(icut)
ELSEIF(COLL_TYPE.EQ.3)THEN
        IF(emep_ISR_shower.EQ.0)THEN
           CALL Cuts_epem(icut)
        ELSE
           CALL QEDPSEVNT(icut) ! imode.EQ.1,phegas_pmom has been boost to lab frame
        ENDIF
ELSE
	icut=1
ENDIF

IF(icut.EQ.0)THEN
   EP_2FS_fxn=0d0
   RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
IF(ptQ)THEN
   kk=2
ELSE
   kk=3
ENDIF
IF(lunwei2)kk=kk+1
IF(.NOT.MCoHelicity)THEN
   DO j=1,lnum
      lr(j)=2*pi*x(kk)
      kk=kk+1
   ENDDO
   DO j=1,snum
      sr(j)=2*pi*x(kk)
      kk=kk+1
   ENDDO
   CALL SetRandom(lr,sr)
ELSE
   RMCoH=x(kk)
   DO j=kk+1,varnum
      Decayran(j-kk)=x(j)
   ENDDO
ENDIF
!IF(NDecayChains.GT.0)THEN
!   CALL HO_Decay(wgtbr)
   ! cut on the decay particles
!   CALL Cuts_Decay(icut)
!   IF(icut.EQ.0)THEN
!      EP_2FS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1)THEN
   reweight_scale_phase=.TRUE.
   CALL Helac_FeynRule_SM()
   CALL Helac_master_f(wme0)
   IF(nnn.EQ.0)THEN
      alphas_power=NINT(DLOG10(wme/wme0))
   ELSE
      IF(NINT(DLOG10(wme/wme0)).NE.alphas_power)THEN
         WRITE(*,*)"WARNING:The process is a mixed order process. Reweighting scale may be wrong !"
      ENDIF
   ENDIF
   reweight_scale_phase=.FALSE.
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   EP_2FS_fxn=0d0
   RETURN
ENDIF
!CALL strf_pdf(wsf)
wsf=1d0
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
EP_2FS_fxn=wme*wsf*LDMEwt/(16d0*pi*ehat**4)&
           *Jpt2*Jx1x2*3.8938573d5
! only unweight events are generated -> LHE files  
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------
   w1=EP_2FS_fxn*wgt ! multiply VEGAS weight
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=EP_2FS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.EP_2FS_fxn)recmax=EP_2FS_fxn
IF(recmin.GT.EP_2FS_fxn)recmin=EP_2FS_fxn
IF(MOD(nnn,nprint).EQ.0)THEN
   PRINT *,"max=",recmax
   PRINT *,"min=",recmin
 IF(lunwei2.AND.lwmax.GT.1)THEN
    PRINT *, "      n_event,    n_pass,    n_total"
    PRINT *, nwri,nnn,nnntot
 ELSE
    PRINT *, "      n_pass,     n_total"
    PRINT *, nnn,nnntot
 ENDIF
ENDIF
END FUNCTION EP_2FS_fxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Cross Section of Electron Positron -> 3 Final States
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EP_3FS(rslt,itmxn,ncalln)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)::rslt   ! rslt(1)=cross section ,rslt(2)=sigma,rslt(3)=chi2
!REAL(KIND=DBL)::sqs,m1,m2,m3                    ! sqs is the total energy of the collider;
                                                ! m1,m2,m3 are the masses of the final states;
!REAL(KIND=DBL),EXTERNAL::ME_Func               ! ME_Func should include |M|^2
                                                ! (with average initial states) and pdf1*pdf2
INTEGER,INTENT(IN),OPTIONAL::itmxn,ncalln
!REAL(KIND=DBL),DIMENSION(2),INTENT(IN)::mu     ! mu(1)=renormalization scale
                                                ! mu(2)=factorization scale
REAL(KIND=DBL)::vfes,sd,chi2a
INTEGER::ncalm,nc,ii
CHARACTER(len=4),DIMENSION(20)::chchar
INTEGER::iday0,ihr0,imin0,isec0,i100th0,iday1,ihr1,imin1,isec1,i100th1,iyr0,iyr1,imon0,imon1
INTEGER::IDBMUP1,IDBMUP2,IDWTUP
CALL ReadElem_logic("ptdisQ",ptQ)
IF(ptQ.AND.emep_ISR_shower.NE.0)STOP "PT distribution is not avialabel when there is ISR shower"
CALL ReadElem_logic('unwgt',lunwei2)
CALL ReadElem_integer('preunw',nunwei)
lwmax=0
NPRN=-1
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
	IF(ptQ)THEN
		varnum=4
	ELSE
		varnum=5
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
ELSEIF(iranhel.EQ.1)THEN
	IF(ptQ)THEN
		varnum=4
	ELSE
		varnum=5
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
	DO ii=1,nhad
		IF(ABS(iflh(ii)).LT.40)THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
	ENDDO
ELSEIF(iranhel.EQ.2)THEN
	IF(ptQ)THEN
		varnum=4
	ELSE
		varnum=5
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
	DO ii=1,nhad
		IF(ABS(iflh(ii)).LT.40)THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
		IF(QN1P1F(ii).OR.QN3PJF(ii))THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
	ENDDO
ELSEIF(iranhel.EQ.3)THEN
	IF(ptQ)THEN
		varnum=4
	ELSE
		varnum=5
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
	DO ii=1,nhad
		IF(ABS(iflh(ii)).LT.40)THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
		IF(QN1P1F(ii))THEN
			varnum=varnum+1
			lnum=lnum+1
		ENDIF
		IF(QN3S1F(ii))THEN
			varnum=varnum+1
			snum=snum+1
		ENDIF
		IF(QN3PJF(ii))THEN
			varnum=varnum+2
			lnum=lnum+1
			snum=snum+1
		ENDIF
	ENDDO
ENDIF
ELSE
   IF(ptQ)THEN
      varnum=4
   ELSE
      varnum=5
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
   varnum=varnum+1+NDecayIflh
ENDIF
DO ii=1,varnum
	XL(ii)=0.0d0
	XU(ii)=1.0d0
ENDDO
IF(PRESENT(itmxn))THEN
   ITMX=itmxn
ELSE
   ITMX=5
END IF
ITMX=1
IF(PRESENT(ncalln))THEN
   ncalm=ncalln
ELSE
   ncalm=5000
END IF

! call vegas and result
nprint=10000
NCALL=20000
!CALL GETDAT(iyr0,imon0,iday0)
!CALL GETTIM(ihr0,imin0,isec0,i100th0)
CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
IF(measurespeed)THEN
   NCALL=500
   CALL VEGAS(varnum,EP_3FS_fxn,vfes,sd,chi2a)
   CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
   iyr1=iyr1-iyr0
   imon1=imon1-imon0
   iday1=iday1-iday0
   ihr1=ihr1-ihr0
   imin1=imin1-imin0
   isec1=isec1-isec0
   i100th1=i100th1-i100th0
   CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
   NCALL=20000
ENDIF
CALL VEGAS(varnum,EP_3FS_fxn,vfes,sd,chi2a)
WRITE(*,*)"====================NCALL=20K==========================="
WRITE(*,*)" "
CALL Helac_mtime()
WRITE(*,*)" "
ii=1
WRITE(*,*)"ITERATION ",ii,":"
WRITE(*,*)vfes,"+\-",sd
WRITE(*,*)"precision:",sd/vfes
DO ii=2,10
        IF(plot_output)CALL initplot
	CALL VEGAS(varnum,EP_3FS_fxn,vfes,sd,chi2a,1)
        IF(plot_output)CALL plotout
	WRITE(*,*)"ITERATION ",ii,":"
	WRITE(*,*)vfes,"+\-",sd
	WRITE(*,*)"precision:",sd/vfes
ENDDO
!CALL GETDAT(iyr1,imon1,iday1)
!CALL GETTIM(ihr1,imin1,isec1,i100th1)
CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
iday1=iday1-iday0
ihr1=ihr1-ihr0
imin1=imin1-imin0
isec1=isec1-isec0
i100th1=i100th1-i100th0
CALL Vegas_speed(10*NCALL,iday1,ihr1,imin1,isec1,i100th1)
WRITE(*,*)''
WRITE(*,*)' '
chchar(1)="40K"
chchar(2)="80K"
chchar(3)="160K"
chchar(4)="320K"
chchar(5)="640K"
chchar(6)="1M"
chchar(7)="2M"
chchar(8)="4M"
chchar(9)="8M"
chchar(10)="16M"
chchar(11)='32M'
chchar(12)='64M'
chchar(13)='120M'
chchar(14)='240M'
chchar(15)='480M'
chchar(16)='960M'
chchar(17)='2G'
chchar(18)='4G'
chchar(19)='8G'
chchar(20)='16G'
ii=1
DO
  nc=2*NCALL
  IF(nc.GT.ncalm)EXIT
  IF(2*nc.GT.ncalm.AND.lunwei2.AND.lwmax.EQ.0)THEN
     lwmax=1
  ENDIF
  NCALL=nc
  IF(NCALL/maxprint.GT.nprint)nprint=NCALL/maxprint
  !CALL GETDAT(iyr0,imon0,iday0)
  !CALL GETTIM(ihr0,imin0,isec0,i100th0)
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_3FS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd  
  WRITE(*,*)"precision:",sd/vfes
  !CALL GETDAT(iyr1,imon1,iday1)
  !CALL GETTIM(ihr1,imin1,isec1,i100th1)
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
  ii=ii+1
ENDDO
IF(lunwei2.AND.lwmax.EQ.1)THEN
   lwmax=2
   WRITE(*,*)"START UNWEIGHTING"
   NCALL=nc/2
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii-1)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_3FS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd
  WRITE(*,*)"precision:",sd/vfes
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
ENDIF
SELECT CASE(COLL_TYPE)
CASE(1)
   IDBMUP1=2212
   IDBMUP2=2212
CASE(2)
   IDBMUP1=2212
   IDBMUP2=-2212
CASE DEFAULT
   IDBMUP1=11
   IDBMUP2=-11
END SELECT
IF(lunwei2)THEN
   IDWTUP=3
ELSE
   IDWTUP=1
ENDIF
NPRUP=NPRUP+1
WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&
     iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
!    &             0,       0,   10042,   10042,IDWTUP,NPRUP 
! XSECUP,XERRUP,XMAXUP,LPRUP 
WRITE(200,5200) vfes*10d0**3,sd*10d0**3,1d0, 81
rslt(1)=vfes
rslt(2)=sd
rslt(3)=chi2a
IF(lunwei2)PRINT *,"number of events",Nevents
IF(emep_ISR_shower.EQ.1)CALL QEDPSEND
RETURN
5100 FORMAT(1P,2I8,2E14.6,6I8)
5200 FORMAT(1P,3E20.10,I6)
END SUBROUTINE EP_3FS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the external function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION EP_3FS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt
REAL(KIND=DBL)::EP_3FS_fxn,m1r,m2r,m3r,sqs,phi0,w1
REAL(KIND=DBL)::pt1,pt12,pt2,y1,y2,y3,phi,mt1,mt2,mt3,sq,m12,m22,m32,aa,ff,ma,hh,cosphi,cosphi2
REAL(KIND=DBL)::temp1,temp2,temp21,temp3,temp4,temp5,temp6,Jpt22,Jpt12,Jy1,Jx1x2,Jphi,&
                lh         ! in 2-> 3, we use t=p2k2 u=p1k2
INTEGER::init=0,kk,j,icut,istat,nnn=0,nnntot=0,pdfnumpdf,nwri,nwri_tot,nwmax,nwarn
REAL(KIND=DBL)::exp1,exp2,ycollcm,y1cup_old,y1clow_old,pt1c_old
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
REAL(KIND=DBL),DIMENSION(3:20)::ptc_old,yycut_old,yycutlow_old
REAL(KIND=DBL)::ptp,pqp,wme,wme0,wgtbr,pi,wsf,y1cup,y1clow,pt1c,recmax=0,recmin=0
REAL(KIND=DBL),DIMENSION(4)::pmom
SAVE init,sqs,sq,m1r,m2r,m3r,m12,m22,m32,temp1,pi,y1cup,y1clow,pt1c,lr,sr,pt1,pt12
SAVE nnn,recmax,recmin,nwri,nwri_tot,nwmax,nwarn,nnntot,ycollcm
SAVE ptc_old,yycut_old,yycutlow_old,y1cup_old,y1clow_old,pt1c_old
IF(init.EQ.0)THEN
	pi=DACOS(-1d0)
	wjac=1
        nwri=0
        nwmax=0
        nwarn=0
        CALL ReadElem_integer('unwevt',nwri_tot)
        CALL ReadElem_real("energy_beam1",sqs)
        ebeam(1)=ABS(sqs)
        CALL ReadElem_real("energy_beam2",sqs)
        ebeam(2)=ABS(sqs)
        IF(ABS(ebeam(1)-ebeam(2))/MAX(ebeam(1)+ebeam(2),1d-17).LT.1d-8)THEN
           absrap=absrap
           labeqcoll=.TRUE.
        ELSE
           IF(emep_ISR_shower.EQ.0)THEN
              absrap=.FALSE.
              labeqcoll=.FALSE.
           ELSE
              absrap=.FALSE.
              labeqcoll=.TRUE. ! run in c.m. frame
           ENDIF
        ENDIF
        sqs=2d0*DSQRT(ebeam(1)*ebeam(2)) ! we always neglect the mass of initial states
        EBMUP1=ebeam(1)
        EBMUP2=ebeam(2)
        IF(emep_ISR_shower.EQ.0)THEN
           ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
        ELSE
           ycollcm=0d0
        ENDIF
	CALL ReadElem_real("minpt1c",pt1c)
	CALL ReadElem_real("maxy1c",y1cup)
	CALL ReadElem_real("miny1c",y1clow)
	CALL ReadElem_integer('pdf',pdfnumpdf)
	ptc(3:20)=0                   ! Debug
        maxptc(3:20)=-1d0
	drc(3:20,3:20)=0
	etac(3:20)=20
	ec(3:20)=0
	c1(3:20)=1
	c2(3:20)=1
	cc(3:20,3:20)=1
        yycut(3:20)=1d9
        IF(absrap)THEN
           yycutlow(3:20)=0d0
        ELSE
           yycutlow(3:20)=-1d9
        ENDIF
	gmas(3:20,3:20)=0             ! Debug
        gbeammass(1:2,3:20)=0

    parmas(0)=0
    parwid(0)=-1

	iPDFSUP1=pdfnumpdf
	IF(pdfnumpdf.EQ.0.OR.COLL_TYPE.EQ.3)THEN
		istruc=0
	ELSE
		istruc=1
	ENDIF
        IF(istruc.EQ.1)THEN
                iPDFGUP1=0
                iPDFGUP2=0
                iPDFSUP1=pdfnumpdf
                iPDFSUP2=pdfnumpdf
        ELSE
                iPDFGUP1=-1
                iPDFGUP2=-1
                iPDFSUP1=-1
                iPDFSUP2=-1
                reweight_pdf=.FALSE.
                ho_npdf=0
        ENDIF
	xp1=1
        xp2=1
	IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
		CALL readcuts_HC
                y1cup=MIN(y1cup,yycut(hadron2parton(3)))
                y1clow=MAX(y1clow,yycutlow(hadron2parton(3)))
                IF(absrap)y1cup=ABS(y1cup)
                IF(absrap)y1clow=ABS(y1clow)
	ELSEIF(COLL_TYPE.EQ.3)THEN
		CALL readcuts_epem
                IF(absrap)y1cup=ABS(y1cup)
                IF(absrap)y1clow=ABS(y1clow)                
	ENDIF
        IF(NDecayChains.GT.0)CALL readcuts_Decay
	kk=3
	istat=0
	IF(lnum.NE.0)THEN
		ALLOCATE(lr(lnum),STAT=istat)
	ENDIF
	IF(snum.NE.0)THEN
		ALLOCATE(sr(snum),STAT=istat)
	ENDIF
	IF(istat.NE.0)THEN
		PRINT *,"Wrong in allocation of array in EP_3FS_fxn !STOP !"
		STOP
	ENDIF
	IF(ABS(iflh(3)).GT.100)THEN
		m1r=parmas(ifl(kk))+parmas(ifl(kk+1))
		kk=kk+2
	ELSE
		m1r=parmas(ifl(kk))
		kk=kk+1
	ENDIF
	IF(ABS(iflh(4)).GT.100)THEN
		m2r=parmas(ifl(kk))+parmas(ifl(kk+1))
		kk=kk+2
	ELSE
		m2r=parmas(ifl(kk))
		kk=kk+1
	ENDIF
	IF(ABS(iflh(5)).GT.100)THEN
		m3r=parmas(ifl(kk))+parmas(ifl(kk+1))
		kk=kk+2
	ELSE
		m3r=parmas(ifl(kk))
		kk=kk+1
	ENDIF
	sq=sqs*sqs
	m12=m1r*m1r
	m22=m2r*m2r
	m32=m3r*m3r
	temp1=(DSQRT(lamda_func(sq,m12,(m2r+m3r)**2))/(2d0*sqs))
        IF(emep_ISR_shower.EQ.1)CALL QEDPSINIT
        ptc_old(3:20)=ptc(3:20)
        yycut_old(3:20)=yycut(3:20)
        yycutlow_old(3:20)=yycutlow(3:20)
        pt1c_old=pt1c
        y1cup_old=y1cup
        y1clow_old=y1clow
	IF(ptQ)THEN
		CALL ReadElem_real("Pt1",pt1)
		pt12=pt1*pt1
	ENDIF
	init=1
ENDIF
nnntot=nnntot+1
IF(emep_ISR_shower.EQ.1)THEN
   CALL QEDPS_ISR_SHOWER
   IF(Q2OUT_QEDPS.LT.(m1r+m2r+m3r)**2)THEN
      EP_3FS_fxn=0d0
      RETURN
   ENDIF
   ebeam(1)=DSQRT(Q2OUT_QEDPS)/2d0
   ebeam(2)=ebeam(1)
   sq=4d0*ebeam(1)*ebeam(2)
   sqs=DSQRT(sq)
   temp1=DSQRT(lamda_func(sq,m12,m22))/(2d0*sqs)
   ptc(3:20)=0d0 ! a small nonzero value to aviod IR divergence
   maxptc(3:20)=-1d0                
   pt1c=0d0
   yycut(3:20)=1d9
   y1cup=1d9
   IF(absrap)THEN
      yycutlow(3:20)=0d0
      y1clow=0d0
   ELSE
      yycutlow(3:20)=-1d9
      y1clow=-1d9
   ENDIF
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The physical region of the integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(.NOT.ptQ)THEN
! total cross section
	pt1=(temp1-pt1c)*x(1)+pt1c
	pt12=pt1*pt1
	mt1=DSQRT(pt12+m12)
        IF(absrap)THEN
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup)
           IF(temp2.LT.y1clow)THEN
              EP_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(2d0*x(2)-1d0)*(temp2-y1clow)+SIGN(1d0,x(2)-0.5d0)*y1clow
        ELSE
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup-ycollcm)
           temp21=MAX(-ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1clow-ycollcm)
           IF(temp2.LT.temp21)THEN
              EP_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(temp2-temp21)*x(2)+temp21 ! in collision frame
        ENDIF
	aa=sq+mt1*mt1-2d0*sqs*mt1*DCOSH(y1)
	ff=m32+pt12
	ma=DMAX1(-lamda_func(aa,m22,ff)/(4d0*m22*pt12),0d0)
	ma=DSQRT(ma)
        IF(lunwei2)THEN
           phi0=2*pi*x(6)
        ELSE
           phi0=0d0
        ENDIF
	phi=unit_step(1d0/2d0-x(3))*(4d0*DACOS(ma)*x(3)-DACOS(ma))&
     +unit_step(x(3)-1d0/2d0)*(2d0*(pi-DACOS(-ma))*(2d0*x(3)-1d0)+DACOS(-ma))&
     +pi+phi0
	cosphi=DCOS(phi)
	cosphi2=cosphi*cosphi
	hh=2d0*pt1*cosphi
	temp3=(DSQRT((4d0*aa-hh**2)*(lamda_func(aa,m22,ff))+((aa-ff+m22)*hh)**2) &
      -hh*(aa-ff+m22))/(4d0*aa-hh**2)
	pt2=x(4)*temp3
	mt2=DSQRT(pt2*pt2+m22)
	mt3=DSQRT(pt12+pt2*pt2+2d0*pt1*pt2*cosphi+m32)
	lh=(sq+m12+m22-m32-2d0*sqs*mt1*DCOSH(y1)-2d0*pt1*pt2*cosphi) &
      /(2d0*mt2*DSQRT(aa))
	IF(lh.LT.1d0)THEN
		EP_3FS_fxn=0d0
		RETURN
	ENDIF
	temp4=ACosh_p(lh)
	y2=SIGN(1d0,x(5)-0.5d0)*temp4+DLOG((sqs-mt1*DEXP(y1))/(sqs-mt1*DEXP(-y1)))/2d0
	temp5=DLOG((sqs-DEXP(-y1)*mt1-DEXP(-y2)*mt2)/mt3)
	y3= -temp5
!Jacobi Determinations
	Jpt12=2d0*pt1*(temp1-pt1c)
        IF(absrap)THEN
           Jy1=2d0*(temp2-y1clow)
        ELSE
           Jy1=(temp2-temp21)
        ENDIF
	Jphi=unit_step(1d0/2d0-x(3))*(4d0*DACOS(ma))&
     +unit_step(x(3)-1d0/2d0)*(4d0*(pi-DACOS(-ma)))
	Jpt22=2d0*pt2*temp3
        temp6=ABS(DSINH(y2-y3))
! a factor of 2 is comming from SIGN(1d0,x(5)-0.5d0)=(2*theta(x(5))-1) 
        Jx1x2=sq/(mt2*mt3*temp6)
ELSE
! pt distribution
	mt1=DSQRT(pt12+m12)
        IF(absrap)THEN
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup)
           IF(temp2.LT.y1clow)THEN
              EP_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(2d0*x(1)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1)-0.5d0)*y1clow
        ELSE
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup-ycollcm)
           temp21=MAX(-ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1clow-ycollcm)
           IF(temp2.LT.temp21)THEN
              EP_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(temp2-temp21)*x(1)+temp21 ! in collision frame
        ENDIF
	aa=sq+mt1*mt1-2d0*sqs*mt1*DCOSH(y1)
	ff=m32+pt12
	ma=DMAX1(-lamda_func(aa,m22,ff)/(4d0*m22*pt12),0d0)
	ma=DSQRT(ma)
        IF(lunwei2)THEN
           phi0=2*pi*x(5)
        ELSE
           phi0=0d0
        ENDIF
	phi=unit_step(1d0/2d0-x(2))*(4d0*DACOS(ma)*x(2)-DACOS(ma))&
     +unit_step(x(2)-1d0/2d0)*(2d0*(pi-DACOS(-ma))*(2d0*x(2)-1d0)+DACOS(-ma))&
     +pi+phi0 ! + pi ?
	cosphi=DCOS(phi)
	cosphi2=cosphi*cosphi
	hh=2d0*pt1*cosphi
	temp3=(DSQRT((4d0*aa-hh**2)*(lamda_func(aa,m22,ff))+((aa-ff+m22)*hh)**2) &
      -hh*(aa-ff+m22))/(4d0*aa-hh**2)
	pt2=x(3)*temp3
	mt2=DSQRT(pt2*pt2+m22)
	mt3=DSQRT(pt12+pt2*pt2+2d0*pt1*pt2*cosphi+m32)
	lh=(sq+m12+m22-m32-2d0*sqs*mt1*DCOSH(y1)-2d0*pt1*pt2*cosphi) &
      /(2d0*mt2*DSQRT(aa))
	IF(lh.LT.1d0)THEN
		EP_3FS_fxn=0d0
		RETURN
	ENDIF
	temp4=ACosh_p(lh)
	y2=SIGN(1d0,x(4)-0.5d0)*temp4+DLOG((sqs-mt1*DEXP(y1))/(sqs-mt1*DEXP(-y1)))/2d0
	temp5=DLOG((sqs-DEXP(-y1)*mt1-DEXP(-y2)*mt2)/mt3)
	y3= -temp5
!Jacobi Determinations
	Jpt12=2d0*pt1
        IF(absrap)THEN
           Jy1=2d0*(temp2-y1clow)
        ELSE
           Jy1=temp2-temp21
        ENDIF
	Jphi=unit_step(1d0/2d0-x(2))*(4d0*DACOS(ma))&
     +unit_step(x(2)-1d0/2d0)*(4d0*(pi-DACOS(-ma)))
	Jpt22=2d0*pt2*temp3
        temp6=ABS(DSINH(y2-y3))
! a factor of 2 is comming from SIGN(1d0,x(5)-0.5d0)=(2*theta(x(5))-1) 
        Jx1x2=sq/(mt2*mt3*temp6)
ENDIF
! The substitution of original variations
xp1=1d0
xp2=1d0
ehat=DSQRT(xp1*xp2*sq)
! Generate the momenta of external legs
Phegas_pmom(1,1:2)=0
Phegas_pmom(1,3)=xp1*sqs/2d0
Phegas_pmom(1,4)=xp1*sqs/2d0
Phegas_pmom(1:2,5)=0
Phegas_pmom(2,1:2)=0
Phegas_pmom(2,3)=-xp2*sqs/2d0
Phegas_pmom(2,4)=xp2*sqs/2d0
! we choose phi=0
!IF(lunwei2)THEN
!   IF(ptQ)THEN
!      phi0=2*pi*x(5)
!   ELSE
!      phi0=2*pi*x(6)
!   ENDIF
!ELSE
!   phi0=0d0
!ENDIF
kk=3
pmom(1)=pt1*DCOS(phi0)
pmom(2)=pt1*DSIN(phi0)
pmom(3)=mt1*DSINH(y1)
pmom(4)=mt1*DCOSH(y1)
IF(ABS(iflh(3)).GT.100)THEN
	Phegas_pmom(kk,1:4)=parmas(ifl(kk))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	Phegas_pmom(kk+1,1:4)=parmas(ifl(kk+1))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk+1,5)=parmas(ifl(kk+1))
	kk=kk+2
ELSE
	Phegas_pmom(kk,1:4)=pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	kk=kk+1
ENDIF
pmom(1)=pt2*DCOS(phi)
pmom(2)=pt2*DSIN(phi)
pmom(3)=mt2*DSINH(y2)
pmom(4)=mt2*DCOSH(y2)
IF(ABS(iflh(4)).GT.100)THEN
	Phegas_pmom(kk,1:4)=parmas(ifl(kk))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	Phegas_pmom(kk+1,1:4)=parmas(ifl(kk+1))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk+1,5)=parmas(ifl(kk+1))
	kk=kk+2
ELSE
	Phegas_pmom(kk,1:4)=pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	kk=kk+1
ENDIF
pmom(1)=-pt2*DCOS(phi)-pt1*DCOS(phi0)
pmom(2)=-pt2*DSIN(phi)-pt1*DSIN(phi0)
pmom(3)=mt3*DSINH(y3)
pmom(4)=mt3*DCOSH(y3)
IF(ABS(iflh(5)).GT.100)THEN
	Phegas_pmom(kk,1:4)=parmas(ifl(kk))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	Phegas_pmom(kk+1,1:4)=parmas(ifl(kk+1))/(parmas(ifl(kk))+parmas(ifl(kk+1)))*pmom(1:4)
	Phegas_pmom(kk+1,5)=parmas(ifl(kk+1))
	kk=kk+2
ELSE
	Phegas_pmom(kk,1:4)=pmom(1:4)
	Phegas_pmom(kk,5)=parmas(ifl(kk))
	kk=kk+1
ENDIF
IF(imode.EQ.0.AND.NDecayChains.EQ.0)THEN
   ! boost from collision frame to cms
   pmom(3)=-(xp1-xp2)*sqs/2d0
   pmom(4)=(xp1+xp2)*sqs/2d0
   pmom(1:2)=0
ELSE
   ! boost from collision frame to lab frame in polarisation
   pmom(3)=ebeam(1)-ebeam(2)
   pmom(4)=ebeam(1)+ebeam(2)
   pmom(1:2)=0
ENDIF
DO j=1,n
	IF(imode.EQ.0.AND.NDecayChains.EQ.0)THEN
           CALL Boostl2(ehat,pmom,Phegas_pmom(j,1:4))
        ELSEIF(.NOT.labeqcoll)THEN
           CALL Boostl2(sqs,pmom,Phegas_pmom(j,1:4))
        ENDIF
	ptp=DSQRT(Phegas_pmom(j,1)**2+Phegas_pmom(j,2)**2)
	pqp=DSQRT(ptp**2+Phegas_pmom(j,3)**2)
	zq(j,1)= DCMPLX(Phegas_pmom(j,4)+Phegas_pmom(j,3),Phegas_pmom(j,3))
	zq(j,2)= DCMPLX(Phegas_pmom(j,4)-Phegas_pmom(j,3),ptp)
	zq(j,3)= DCMPLX(Phegas_pmom(j,1), Phegas_pmom(j,2))
	zq(j,4)= DCMPLX(Phegas_pmom(j,1),-Phegas_pmom(j,2))
	zq(j,5)= DCMPLX(parmas(ifl(j)),pqp )
ENDDO
! recover the cuts                                                              
IF(emep_ISR_shower.EQ.1)THEN
   pt1c=pt1c_old
   ptc(3:20)=ptc_old(3:20)
   yycut(3:20)=yycut_old(3:20)
   yycutlow(3:20)=yycutlow_old(3:20)
   y1cup=y1cup_old
   y1clow=y1clow_old
ENDIF
! cut something
IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
	CALL Cuts_HC(icut)
ELSEIF(COLL_TYPE.EQ.3)THEN
        IF(emep_ISR_shower.EQ.0)THEN
           CALL Cuts_epem(icut)
        ELSE
           CALL QEDPSEVNT(icut) ! imode.EQ.1, phegas_pmom has been boost to lab frame
        ENDIF
ELSE
	icut=1
ENDIF
IF(icut.EQ.0)THEN
	EP_3FS_fxn=0d0
	RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
IF(ptQ)THEN
	kk=5
ELSE
	kk=6
ENDIF
IF(lunwei2)kk=kk+1
IF(.NOT.MCoHelicity)THEN
   DO j=1,lnum
      lr(j)=2*pi*x(kk)
      kk=kk+1
   ENDDO
   DO j=1,snum
      sr(j)=2*pi*x(kk)
      kk=kk+1
   ENDDO
   CALL SetRandom(lr,sr)
ELSE
   RMCoH=x(kk)
   DO j=kk+1,varnum
      Decayran(j-kk)=x(j)
   ENDDO
ENDIF
!IF(NDecayChains.GT.0)THEN
!   CALL HO_Decay(wgtbr)
!   CALL cuts_Decay(icut)
!   IF(icut.EQ.0)THEN
!      EP_3FS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1)THEN
   reweight_scale_phase=.TRUE.
   CALL Helac_FeynRule_SM()
   CALL Helac_master_f(wme0)
   IF(nnn.EQ.0)THEN
      alphas_power=NINT(DLOG10(wme/wme0))
   ELSE
      IF(NINT(DLOG10(wme/wme0)).NE.alphas_power)THEN
         WRITE(*,*)"WARNING:The process is a mixed order process.Reweighting scale may be wrong !"
      ENDIF
   ENDIF
   reweight_scale_phase=.FALSE.
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   EP_3FS_fxn=0d0
   RETURN
ENDIF
! Generate weight for pdf
!CALL strf_pdf(wsf)
wsf=1d0
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
EP_3FS_fxn=wme*wsf*LDMEwt/(512d0*(pi**4)*ehat**4)&
           *Jpt12*Jpt22*Jy1*Jx1x2*Jphi*3.8938573d5
! only unweight events are generated -> LHE files                                                                                                                                                           
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------   
   w1=EP_3FS_fxn*wgt ! multiply VEGAS weight  
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=EP_3FS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.EP_3FS_fxn)recmax=EP_3FS_fxn
IF(recmin.GT.EP_3FS_fxn)recmin=EP_3FS_fxn
IF(MOD(nnn,nprint).EQ.0)THEN
   PRINT *,"max=",recmax
   PRINT *,"min=",recmin
 IF(lunwei2.AND.lwmax.GT.1)THEN
    PRINT *, "      n_event,    n_pass,    n_total"
    PRINT *, nwri,nnn,nnntot
 ELSE
    PRINT *, "      n_pass,     n_total"
    PRINT *, nnn,nnntot
 ENDIF
ENDIF
END FUNCTION EP_3FS_fxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Cross Section of Electron Positron -> n Final States with n>=2                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE EP_nFS(rslt,itmxn,ncalln)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(3),INTENT(OUT)::rslt   ! rslt(1)=cross section ,rslt(2)=sigma,rslt(3)=chi2
!REAL(KIND=DBL)::sqs,m1,m2,m3                    ! sqs is the total energy of the collider;
                                                ! m1,m2,m3 are the masses of the final states;
!REAL(KIND=DBL),EXTERNAL::ME_Func               ! ME_Func should include |M|^2
                                                ! (with average initial states) and pdf1*pdf2
INTEGER,INTENT(IN),OPTIONAL::itmxn,ncalln
!REAL(KIND=DBL),DIMENSION(2),INTENT(IN)::mu     ! mu(1)=renormalization scale
                                                ! mu(2)=factorization scale
REAL(KIND=DBL)::vfes,sd,chi2a
INTEGER::ncalm,nc,ii
CHARACTER(len=4),DIMENSION(20)::chchar
INTEGER::iday0,ihr0,imin0,isec0,i100th0,iday1,ihr1,imin1,isec1,i100th1,iyr0,iyr1,imon0,imon1
INTEGER::IDBMUP1,IDBMUP2,IDWTUP
CALL ReadElem_logic("ptdisQ",ptQ)
IF(ptQ.AND.emep_ISR_shower.NE.0)STOP "PT distribution is not avialabel when there is ISR shower"
CALL ReadElem_logic('unwgt',lunwei2)
CALL ReadElem_integer('preunw',nunwei)
lwmax=0
NPRN=-1
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
   IF(ptQ)THEN
      varnum=3*(nhad-3)-2
   ELSE
      varnum=3*(nhad-3)-1
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
ELSEIF(iranhel.EQ.1)THEN
   IF(ptQ)THEN
      varnum=3*(nhad-3)-2
   ELSE
      varnum=3*(nhad-3)-1
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
   DO ii=1,nhad
      IF(ABS(iflh(ii)).LT.40)THEN
         varnum=varnum+1
         lnum=lnum+1
      ENDIF
   ENDDO
ELSEIF(iranhel.EQ.2)THEN
   IF(ptQ)THEN
      varnum=3*(nhad-3)-2
   ELSE
      varnum=3*(nhad-3)-1
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
   DO ii=1,nhad
      IF(ABS(iflh(ii)).LT.40)THEN
         varnum=varnum+1
         lnum=lnum+1
      ENDIF
      IF(QN1P1F(ii).OR.QN3PJF(ii))THEN
         varnum=varnum+1
         lnum=lnum+1
      ENDIF
   ENDDO
ELSEIF(iranhel.EQ.3)THEN
   IF(ptQ)THEN
      varnum=3*(nhad-3)-2
   ELSE
      varnum=3*(nhad-3)-1
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
   DO ii=1,nhad
      IF(ABS(iflh(ii)).LT.40)THEN
         varnum=varnum+1
         lnum=lnum+1
      ENDIF
      IF(QN1P1F(ii))THEN
         varnum=varnum+1
         lnum=lnum+1
      ENDIF
      IF(QN3S1F(ii))THEN
         varnum=varnum+1
         snum=snum+1
      ENDIF
      IF(QN3PJF(ii))THEN
         varnum=varnum+2
         lnum=lnum+1
         snum=snum+1
      ENDIF
   ENDDO
ENDIF
ELSE
   IF(ptQ)THEN
      varnum=3*(nhad-3)-2
   ELSE
      varnum=3*(nhad-3)-1
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
   varnum=varnum+1+NDecayIflh
ENDIF
DO ii=1,varnum
   XL(ii)=0.0d0
   XU(ii)=1.0d0
ENDDO
IF(PRESENT(itmxn))THEN
   ITMX=itmxn
ELSE
   ITMX=5
END IF
ITMX=1
IF(PRESENT(ncalln))THEN
   ncalm=ncalln
ELSE
   ncalm=5000
END IF

! call vegas and result
nprint=10000
NCALL=20000
!CALL GETDAT(iyr0,imon0,iday0)
!CALL GETTIM(ihr0,imin0,isec0,i100th0)
CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
IF(measurespeed)THEN
   NCALL=500
   CALL VEGAS(varnum,EP_nFS_fxn,vfes,sd,chi2a)
   CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
   iyr1=iyr1-iyr0
   imon1=imon1-imon0
   iday1=iday1-iday0
   ihr1=ihr1-ihr0
   imin1=imin1-imin0
   isec1=isec1-isec0
   i100th1=i100th1-i100th0
   CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
   NCALL=20000
ENDIF
CALL VEGAS(varnum,EP_nFS_fxn,vfes,sd,chi2a)
WRITE(*,*)"====================NCALL=20K==========================="
WRITE(*,*)" "
CALL Helac_mtime()
WRITE(*,*)" "
ii=1
WRITE(*,*)"ITERATION ",ii,":"
WRITE(*,*)vfes,"+\-",sd
WRITE(*,*)"precision:",sd/vfes
DO ii=2,10
   IF(plot_output)CALL initplot
   CALL VEGAS(varnum,EP_nFS_fxn,vfes,sd,chi2a,1)
   IF(plot_output)CALL plotout
   WRITE(*,*)"ITERATION ",ii,":"
   WRITE(*,*)vfes,"+\-",sd
   WRITE(*,*)"precision:",sd/vfes
ENDDO
!CALL GETDAT(iyr1,imon1,iday1)
!CALL GETTIM(ihr1,imin1,isec1,i100th1)
CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
iday1=iday1-iday0
ihr1=ihr1-ihr0
imin1=imin1-imin0
isec1=isec1-isec0
i100th1=i100th1-i100th0
CALL Vegas_speed(10*NCALL,iday1,ihr1,imin1,isec1,i100th1)
WRITE(*,*)''
WRITE(*,*)' '
chchar(1)="40K"
chchar(2)="80K"
chchar(3)="160K"
chchar(4)="320K"
chchar(5)="640K"
chchar(6)="1M"
chchar(7)="2M"
chchar(8)="4M"
chchar(9)="8M"
chchar(10)="16M"
chchar(11)='32M'
chchar(12)='64M'
chchar(13)='120M'
chchar(14)='240M'
chchar(15)='480M'
chchar(16)='960M'
chchar(17)='2G'
chchar(18)='4G'
chchar(19)='8G'
chchar(20)='16G'
ii=1
DO
   nc=2*NCALL
   IF(nc.GT.ncalm)EXIT
   IF(2*nc.GT.ncalm.AND.lunwei2.AND.lwmax.EQ.0)THEN
      lwmax=1
   ENDIF
   NCALL=nc
   IF(NCALL/maxprint.GT.nprint)nprint=NCALL/maxprint
   !CALL GETDAT(iyr0,imon0,iday0)
   !CALL GETTIM(ihr0,imin0,isec0,i100th0)
   CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
   WRITE(*,*)"====================NCALL="//chchar(ii)//"==========================="
   WRITE(*,*)" "
   CALL Helac_mtime()
   WRITE(*,*)" "
   IF(plot_output)CALL initplot
   CALL VEGAS(varnum,EP_nFS_fxn,vfes,sd,chi2a,1)
   IF(plot_output)CALL plotout
   WRITE(*,*)vfes,"+\-",sd
   WRITE(*,*)"precision:",sd/vfes
   !CALL GETDAT(iyr1,imon1,iday1)
   !CALL GETTIM(ihr1,imin1,isec1,i100th1)
   CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
   iday1=iday1-iday0
   ihr1=ihr1-ihr0
   imin1=imin1-imin0
   isec1=isec1-isec0
   i100th1=i100th1-i100th0
   CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
   WRITE(*,*)" "
   ii=ii+1
ENDDO
IF(lunwei2.AND.lwmax.EQ.1)THEN
   lwmax=2
   WRITE(*,*)"START UNWEIGHTING"
   NCALL=nc/2
  CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
  WRITE(*,*)"====================NCALL="//chchar(ii-1)//"==========================="
  WRITE(*,*)" "
  CALL Helac_mtime()
  WRITE(*,*)" "
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,EP_nFS_fxn,vfes,sd,chi2a,1)
  IF(plot_output)CALL plotout
  WRITE(*,*)vfes,"+\-",sd
  WRITE(*,*)"precision:",sd/vfes
  CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
  iday1=iday1-iday0
  ihr1=ihr1-ihr0
  imin1=imin1-imin0
  isec1=isec1-isec0
  i100th1=i100th1-i100th0
  CALL Vegas_speed(NCALL,iday1,ihr1,imin1,isec1,i100th1)
  WRITE(*,*)" "
ENDIF
SELECT CASE(COLL_TYPE)
CASE(1)
   IDBMUP1=2212
   IDBMUP2=2212
CASE(2)
   IDBMUP1=2212
   IDBMUP2=-2212
CASE DEFAULT
   IDBMUP1=11
   IDBMUP2=-11
END SELECT
IF(lunwei2)THEN
   IDWTUP=3
ELSE
   IDWTUP=1
ENDIF
NPRUP=NPRUP+1
WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&
     iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
!    &             0,       0,   10042,   10042,IDWTUP,NPRUP
! XSECUP,XERRUP,XMAXUP,LPRUP 
WRITE(200,5200) vfes*10d0**3,sd*10d0**3,1d0, 81
rslt(1)=vfes
rslt(2)=sd
rslt(3)=chi2a
IF(lunwei2)PRINT *,"number of events",Nevents
IF(emep_ISR_shower.EQ.1)CALL QEDPSEND
RETURN
5100 FORMAT(1P,2I8,2E14.6,6I8)
5200 FORMAT(1P,3E20.10,I6)
END SUBROUTINE EP_nFS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the external function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION EP_nFS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt
REAL(KIND=DBL)::EP_nFS_fxn,sqs,sq
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::mr,ptcut,ymincut,ymaxcut
REAL(KIND=DBL),DIMENSION(nhad,1:4)::phadmom
REAL(KIND=DBL)::pt0,pti,phi,phi0,y0,yi,yn,Q,mi,mres,jac,jactemp,mtot,pt1,pt12,mti
REAL(KIND=DBL)::rlogmax,rlogmin,w1,mtnm1,ynm1
REAL(KIND=DBL)::exp1,exp2,ycollcm,y1cup_old,y1clow_old,pt1c_old
LOGICAL::lflag=.TRUE.
INTEGER::init=0,kk,kk2,ii,j,icut,istat,nnn=0,nnntot=0,pdfnumpdf,nwri,nwri_tot,nwarn,nwmax
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
REAL(KIND=DBL),DIMENSION(3:20)::ptc_old,yycut_old,yycutlow_old
REAL(KIND=DBL)::ptp,pqp,wme,wme0,wgtbr,pi,wsf,y1cup,y1clow,pt1c,recmax=0,recmin=0
REAL(KIND=DBL),DIMENSION(4)::pmom
SAVE init,sqs,sq,pi,y1cup,y1clow,pt1c,lr,sr,pt1,pt12,mr,mtot,nnntot,ycollcm
SAVE nnn,recmax,recmin,ptcut,ymincut,ymaxcut,nwri,nwri_tot,nwarn,nwmax
SAVE ptc_old,yycut_old,yycutlow_old,y1cup_old,y1clow_old,pt1c_old
IF(init.EQ.0)THEN
        pi=DACOS(-1d0)
        wjac=1
        nwmax=0
        nwarn=0
        nwri=0
        CALL ReadElem_integer('unwevt',nwri_tot)
        CALL ReadElem_real("energy_beam1",sqs)
        ebeam(1)=ABS(sqs)
        CALL ReadElem_real("energy_beam2",sqs)
        ebeam(2)=ABS(sqs)
        IF(ABS(ebeam(1)-ebeam(2))/MAX(ebeam(1)+ebeam(2),1d-17).LT.1d-8)THEN
           absrap=absrap
           labeqcoll=.TRUE.
        ELSE
           IF(emep_ISR_shower.EQ.0)THEN
              absrap=.FALSE.
              labeqcoll=.FALSE.
           ELSE
              absrap=.FALSE.
              labeqcoll=.TRUE. ! run in c.m. frame
           ENDIF
        ENDIF
        sqs=2d0*DSQRT(ebeam(1)*ebeam(2)) ! we always neglect the mass of initial states
        EBMUP1=ebeam(1)
        EBMUP2=ebeam(2)
        IF(emep_ISR_shower.EQ.0)THEN
           ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
        ELSE
           ycollcm=0d0
        ENDIF
        CALL ReadElem_real("minpt1c",pt1c)
        CALL ReadElem_real("maxy1c",y1cup)
        CALL ReadElem_real("miny1c",y1clow)
        CALL ReadElem_integer('pdf',pdfnumpdf)
        ptc(3:20)=0                   ! Debug
        maxptc(3:20)=-1d0
        drc(3:20,3:20)=0
        etac(3:20)=20
        ec(3:20)=0
        c1(3:20)=1
        c2(3:20)=1
        cc(3:20,3:20)=1
        yycut(3:20)=1d9
        IF(absrap)THEN
           yycutlow(3:20)=0d0
        ELSE
           yycutlow(3:20)=-1d9
        ENDIF
        gmas(3:20,3:20)=0             ! Debug
        gbeammass(1:2,3:20)=0

        parmas(0)=0
        parwid(0)=-1

        iPDFSUP1=pdfnumpdf
        IF(pdfnumpdf.EQ.0.OR.COLL_TYPE.EQ.3)THEN
                istruc=0
        ELSE
                istruc=1
        ENDIF
        IF(istruc.EQ.1)THEN
                iPDFGUP1=0
                iPDFGUP2=0
                iPDFSUP1=pdfnumpdf
                iPDFSUP2=pdfnumpdf
        ELSE
                iPDFGUP1=-1
                iPDFGUP2=-1
                iPDFSUP1=-1
                iPDFSUP2=-1
                reweight_pdf=.FALSE.
                ho_npdf=0
        ENDIF
        xp1=1
        xp2=1
        IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
                CALL readcuts_HC
                y1cup=MIN(y1cup,yycut(hadron2parton(3)))
                y1clow=MAX(y1clow,yycutlow(hadron2parton(3)))
                IF(absrap)y1cup=ABS(y1cup)
                IF(absrap)y1clow=ABS(y1clow)
        ELSEIF(COLL_TYPE.EQ.3)THEN
                CALL readcuts_epem
                IF(absrap)y1cup=ABS(y1cup)
                IF(absrap)y1clow=ABS(y1clow)
        ENDIF
        IF(NDecayChains.GT.0)CALL readcuts_Decay
        istat=0
        IF(lnum.NE.0)THEN
           ALLOCATE(lr(lnum),STAT=istat)
        ENDIF
        IF(snum.NE.0)THEN
           ALLOCATE(sr(snum),STAT=istat)
        ENDIF
        ALLOCATE(mr(nhad),STAT=istat)
        ALLOCATE(ptcut(nhad),STAT=istat)
        ALLOCATE(ymincut(nhad),STAT=istat)
        ALLOCATE(ymaxcut(nhad),STAT=istat)
        IF(istat.NE.0)THEN
                PRINT *,"Wrong in allocation of array in EP_nFS_fxn !STOP !"
                STOP
        ENDIF
        kk=3
        mtot=0d0
        mr(1)=0d0
        mr(2)=0d0
        ptcut(1:2)=0d0
        ymincut(1:2)=0d0
        ymaxcut(1:2)=0d0
        DO ii=3,nhad
           IF(ABS(iflh(ii)).GT.100)THEN
              mr(ii)=parmas(ifl(kk))+parmas(ifl(kk+1))
              IF(absrap)THEN
                 ymincut(ii)=yycutlow(kk)
                 ymaxcut(ii)=yycut(kk)
              ELSE
                 ymincut(ii)=yycutlow(kk)-ycollcm
                 ymaxcut(ii)=yycut(kk)-ycollcm
              ENDIF
              ptcut(ii)=ptc(kk)+ptc(kk+1)
              kk=kk+2
           ELSE
              mr(ii)=parmas(ifl(kk))
              IF(absrap)THEN
                 ymincut(ii)=yycutlow(kk)
                 ymaxcut(ii)=yycut(kk)
              ELSE
                 ymincut(ii)=yycutlow(kk)-ycollcm
                 ymaxcut(ii)=yycut(kk)-ycollcm
              ENDIF
              ptcut(ii)=ptc(kk)
              kk=kk+1
           ENDIF
           IF(ymincut(ii).GT.ymaxcut(ii))THEN
              WRITE(*,*)"All events are cutted off !"
              WRITE(*,*)"k, ycmin, ycmax"
              WRITE(*,*)ii,ymincut(ii),ymaxcut(ii)
              WRITE(*,*)"BYE !"
              STOP
           ENDIF
           mtot=mtot+mr(ii)
        ENDDO
        IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
           ptcut(3)=MAX(ptcut(3),pt1c)
           ymincut(3)=MAX(y1clow,ymincut(3))
           ymaxcut(3)=MIN(y1cup,ymaxcut(3))
        ELSE
           ptcut(3)=MAX(ptcut(3),pt1c)
           ymincut(3)=MAX(y1clow-ycollcm,ymincut(3))
           ymaxcut(3)=MIN(y1cup-ycollcm,ymaxcut(3))
        ENDIF
        IF(ymincut(3).GE.ymaxcut(3))THEN
           WRITE(*,*)"All events are cutted off !"
           WRITE(*,*)"k, ycmin, ycmax"
           WRITE(*,*)3,ymincut(3),ymaxcut(3)
           WRITE(*,*)"BYE !"
           STOP
        ENDIF
        sq=sqs*sqs
        IF(emep_ISR_shower.EQ.1)THEN
           CALL QEDPSINIT
           ptcut(1:nhad)=0d0
           IF(absrap)THEN
              ymincut(1:nhad)=0d0
           ELSE
              ymincut(1:nhad)=-1d9
           ENDIF
           ymaxcut(1:nhad)=1d9
        ENDIF
        ptc_old(3:20)=ptc(3:20)
        yycut_old(3:20)=yycut(3:20)
        yycutlow_old(3:20)=yycutlow(3:20)
        pt1c_old=pt1c
        y1cup_old=y1cup
        y1clow_old=y1clow
        IF(ptQ)THEN
           CALL ReadElem_real("Pt1",pt1)
           pt12=pt1*pt1
        ENDIF
        init=1
ENDIF
nnntot=nnntot+1
IF(emep_ISR_shower.EQ.1)THEN
   CALL QEDPS_ISR_SHOWER
   IF(Q2OUT_QEDPS.LT.mtot**2)THEN
      EP_nFS_fxn=0d0
      RETURN
   ENDIF
   ebeam(1)=DSQRT(Q2OUT_QEDPS)/2d0
   ebeam(2)=ebeam(1)
   sq=4d0*ebeam(1)*ebeam(2)
   sqs=DSQRT(sq)
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The physical region of the integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
jac=1d0
phadmom(1,1:2)=0
phadmom(1,3)=sqs/2d0
phadmom(1,4)=sqs/2d0
phadmom(2,1:2)=0
phadmom(2,3)=-sqs/2d0
phadmom(2,4)=sqs/2d0
pmom(1:4)=phadmom(1,1:4)+phadmom(2,1:4)
rlogmax=sqs
rlogmin=sqs
xp1=0d0
xp2=0d0
kk=1
mres=mtot
DO ii=3,nhad-1
   Q=scalarproduction2(pmom(1:4))
   mi=mr(ii)
   mres=mres-mi
   phi0=ph42(pmom)
   pt0=transverse2(pmom)
   y0=rapidity2(pmom)
   ! we take the first phi as trivial
   ! generate phi
   IF(ii.EQ.3)THEN
      IF(lunwei2)THEN
         phi=2*pi*x(kk)
         kk=kk+1
      ELSE
         phi=0d0
      ENDIF
      jactemp=2*pi
      lflag=.TRUE.
   ELSEIF(ii.EQ.nhad-1)THEN
      CALL generate_phi(1,Q,mi,mres,phi0,pt0,x(kk),phi,jactemp,lflag)
!      phi=2*pi*x(kk) ! Debug
!      jactemp=2*pi
!      lflag=.TRUE.
      kk=kk+1
   ELSE
      CALL generate_phi(1,Q,mi,mres,phi0,pt0,x(kk),phi,jactemp,lflag)
!      phi=2*pi*x(kk) ! Debug
!      jactemp=2*pi
!      lflag=.TRUE.
      kk=kk+1
   ENDIF
   IF(.NOT.lflag)THEN
      EP_nFS_fxn=0d0
      RETURN
   ENDIF
   jac=jac*jactemp
   ! generate pt
   IF(.NOT.ptQ.OR.ii.NE.3)THEN
      CALL generate_pt(Q,mi,mres,phi0,phi,pt0,x(kk),ptcut(ii),pti,jactemp,lflag)
      kk=kk+1
   ELSE
      ! pt distribution for the first final state particle
      pti=pt1
      lflag=.TRUE.
      jactemp=1d0
   ENDIF
   IF(.NOT.lflag)THEN
      EP_nFS_fxn=0d0
      RETURN
   ENDIF
   jac=jac*jactemp*2d0*pti ! from d pti^2 -> d pti
   ! generate yi
   IF(ii.NE.nhad-1)THEN
      CALL generate_yi(Q,mi,mres,phi0,phi,pt0,pti,y0,x(kk),ymaxcut(ii),ymincut(ii),&
           yi,jactemp,lflag)
   ELSE
      CALL assign_ynm1(Q,mi,mres,phi0,phi,pt0,pti,y0,x(kk),ymaxcut(ii),ymincut(ii),&
           yi,jactemp,lflag)
   ENDIF
   kk=kk+1
   IF(.NOT.lflag)THEN
      EP_nFS_fxn=0d0
      RETURN
   ENDIF
   jac=jac*jactemp
   ! momentum assignment
   mti=DSQRT(mi**2+pti**2)
   phadmom(ii,4)=mti*DCOSH(yi)
   phadmom(ii,3)=mti*DSINH(yi)
   phadmom(ii,1)=pti*DCOS(phi)
   phadmom(ii,2)=pti*DSIN(phi)
   pmom(1:4)=pmom(1:4)-phadmom(ii,1:4)
   rlogmax=rlogmax-mti*DEXP(yi)
   rlogmin=rlogmin-mti*DEXP(-yi)
   IF(ii.EQ.nhad-1)THEN
      mtnm1=mti
      ynm1=yi
   ENDIF
!   xp1=xp1+DEXP(yi)*mti
!   xp2=xp2+DEXP(-yi)*mti
ENDDO
! generate yn
IF(rlogmax.LE.0d0.OR.rlogmin.LE.0d0)THEN
   EP_nFS_fxn=0d0
   RETURN
ENDIF
pti=transverse2(pmom)
mi=mr(nhad)
mti=DSQRT(pti**2+mi**2)
rlogmax=DLOG(rlogmax/mti)
rlogmin=-DLOG(rlogmin/mti)
yn=rlogmax
! weight for 1/ABS(Det({{dxp1/dynm1,dxp1/dyn},{dxp2/dynm1,dxp2/dyn}}))
! a factor of 2 is comming from SIGN(1d0,r-0.5d0)=(2*theta(r)-1) 
jac=jac*sq/(mti*mtnm1*ABS(DSINH(ynm1-yn)))
!phadmom(nhad,3)=mti*DSINH(yn)
!phadmom(nhad,4)=mti*DCOSH(yn)
!xp1=xp1+DEXP(yn)*mti
!xp2=xp2+DEXP(-yn)*mti
! The substitution of original variations
xp1=1d0
xp2=1d0
exp1=ebeam(1)*xp1
exp2=ebeam(2)*xp2
ehat=DSQRT(xp1*xp2*sq)
! Generate the momenta of external legs
! initial states
Phegas_pmom(1,1:2)=0
Phegas_pmom(1,3)=xp1*sqs/2d0
Phegas_pmom(1,4)=xp1*sqs/2d0
Phegas_pmom(1:2,5)=0
Phegas_pmom(2,1:2)=0
Phegas_pmom(2,3)=-xp2*sqs/2d0
Phegas_pmom(2,4)=xp2*sqs/2d0
! final states
kk2=3
phadmom(nhad,1:4)=Phegas_pmom(1,1:4)+Phegas_pmom(2,1:4)
DO ii=3,nhad-1
   pmom(1)=phadmom(ii,1)
   pmom(2)=phadmom(ii,2)
   pmom(3)=phadmom(ii,3)
   pmom(4)=phadmom(ii,4)
   phadmom(nhad,1)=phadmom(nhad,1)-phadmom(ii,1)
   phadmom(nhad,2)=phadmom(nhad,2)-phadmom(ii,2)
   phadmom(nhad,3)=phadmom(nhad,3)-phadmom(ii,3)
   phadmom(nhad,4)=phadmom(nhad,4)-phadmom(ii,4)
   IF(ABS(iflh(ii)).GT.100)THEN
      Phegas_pmom(kk2,1:4)=parmas(ifl(kk2))/(parmas(ifl(kk2))+parmas(ifl(kk2+1)))*pmom(1:4)
      Phegas_pmom(kk2,5)=parmas(ifl(kk2))
      Phegas_pmom(kk2+1,1:4)=parmas(ifl(kk2+1))/(parmas(ifl(kk2))+parmas(ifl(kk2+1)))*pmom(1:4)
      Phegas_pmom(kk2+1,5)=parmas(ifl(kk2+1))
      kk2=kk2+2
   ELSE
      Phegas_pmom(kk2,1:4)=pmom(1:4)
      Phegas_pmom(kk2,5)=parmas(ifl(kk2))
      kk2=kk2+1
   ENDIF
ENDDO
! for the nhad-th particle
pmom(1)=phadmom(nhad,1)
pmom(2)=phadmom(nhad,2)
pmom(3)=phadmom(nhad,3)
pmom(4)=phadmom(nhad,4)
IF(ABS(iflh(nhad)).GT.100)THEN
   Phegas_pmom(kk2,1:4)=parmas(ifl(kk2))/(parmas(ifl(kk2))+parmas(ifl(kk2+1)))*pmom(1:4)
   Phegas_pmom(kk2,5)=parmas(ifl(kk2))
   Phegas_pmom(kk2+1,1:4)=parmas(ifl(kk2+1))/(parmas(ifl(kk2))+parmas(ifl(kk2+1)))*pmom(1:4)
   Phegas_pmom(kk2+1,5)=parmas(ifl(kk2+1))
   kk2=kk2+2
ELSE
   Phegas_pmom(kk2,1:4)=pmom(1:4)
   Phegas_pmom(kk2,5)=parmas(ifl(kk2))
   kk2=kk2+1
ENDIF
IF(imode.EQ.0.AND.NDecayChains.EQ.0)THEN
   ! boost from collision frame to cms
   pmom(3)=-(xp1-xp2)*sqs/2d0
   pmom(4)=(xp1+xp2)*sqs/2d0
   pmom(1:2)=0
ELSE
   ! boost from collision frame to lab frame
   pmom(3)=ebeam(1)-ebeam(2)
   pmom(4)=ebeam(1)+ebeam(2)
   pmom(1:2)=0
ENDIF
DO j=1,n
   IF(imode.EQ.0.AND.NDecayChains.EQ.0)THEN
      CALL Boostl2(ehat,pmom,Phegas_pmom(j,1:4))
   ELSEIF(.NOT.labeqcoll)THEN
      CALL Boostl2(sqs,pmom,Phegas_pmom(j,1:4))
   ENDIF
   ptp=DSQRT(Phegas_pmom(j,1)**2+Phegas_pmom(j,2)**2)
   pqp=DSQRT(ptp**2+Phegas_pmom(j,3)**2)
   zq(j,1)= DCMPLX(Phegas_pmom(j,4)+Phegas_pmom(j,3),Phegas_pmom(j,3))
   zq(j,2)= DCMPLX(Phegas_pmom(j,4)-Phegas_pmom(j,3),ptp)
   zq(j,3)= DCMPLX(Phegas_pmom(j,1), Phegas_pmom(j,2))
   zq(j,4)= DCMPLX(Phegas_pmom(j,1),-Phegas_pmom(j,2))
   zq(j,5)= DCMPLX(parmas(ifl(j)),pqp )
ENDDO
! cut something
IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
   CALL Cuts_HC(icut)
ELSEIF(COLL_TYPE.EQ.3)THEN
   IF(emep_ISR_shower.EQ.0)THEN
      CALL Cuts_epem(icut)
   ELSE
      CALL QEDPSEVNT(icut) ! imode.EQ.1,phegas_pmom has been boost to lab frame
   ENDIF
ELSE
   icut=1
ENDIF
IF(icut.EQ.0)THEN
   EP_nFS_fxn=0d0
   RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
!IF(ptQ)THEN
!        kk=6
!ELSE
!        kk=5
!ENDIF
IF(.NOT.MCoHelicity)THEN
   DO j=1,lnum
      lr(j)=2*pi*x(kk)
      kk=kk+1
   ENDDO
   DO j=1,snum
      sr(j)=2*pi*x(kk)
      kk=kk+1
   ENDDO
   CALL SetRandom(lr,sr)
ELSE
   RMCoH=x(kk)
   DO j=kk+1,varnum
      Decayran(j-kk)=x(kk)
   ENDDO
ENDIF
!IF(NDecayChains.GT.0)THEN
!   CALL HO_Decay(wgtbr)
!   CALL cuts_Decay(icut)
!   IF(icut.EQ.0)THEN
!      EP_nFS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1)THEN
   reweight_scale_phase=.TRUE.
   CALL Helac_FeynRule_SM()
   CALL Helac_master_f(wme0)
   IF(nnn.EQ.0)THEN
      alphas_power=NINT(DLOG10(wme/wme0))
   ELSE
      IF(NINT(DLOG10(wme/wme0)).NE.alphas_power)THEN
         WRITE(*,*)"WARING:The process is a mixed order process.Reweighting scale may be wrong !"
      ENDIF
   ENDIF
   reweight_scale_phase=.FALSE.
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   EP_nFS_fxn=0d0
   RETURN
ENDIF
! Generate weight for pdf
!CALL strf_pdf(wsf)
wsf=1d0
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
EP_nFS_fxn=wme*wsf*LDMEwt/(4d0*(4d0*pi)**(2*nhad-7)*(2*pi)**(nhad-3)*ehat**4)&
           *Jac*3.8938573d5
! only unweight events are generated -> LHE files 
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------   
   w1=EP_nFS_fxn*wgt ! multiply VEGAS weight  
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=EP_nFS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.EP_nFS_fxn)recmax=EP_nFS_fxn
IF(recmin.GT.EP_nFS_fxn)recmin=EP_nFS_fxn
IF(MOD(nnn,nprint).EQ.0)THEN
   PRINT *,"max=",recmax
   PRINT *,"min=",recmin
   IF(lunwei2.AND.lwmax.GT.1)THEN
      PRINT *, "      n_event,    n_pass,    n_total"
      PRINT *, nwri,nnn,nnntot
   ELSE
      PRINT *, "      n_pass,     n_total"
      PRINT *, nnn,nnntot
   ENDIF
ENDIF
END FUNCTION EP_nFS_fxn

SUBROUTINE assign_ynm1(Q,mi,mres,phi0,phi,pt,pti,y0,xr,ycmax,ycmin,&
     yi,wgt,lflag)
! assign yi for i=n-1, Q -> ki+Q2, mres<Q2<(Q-mi) and Q2=mn**2
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::Q,mi,mres,phi0,phi,pt,pti,y0,xr,ycmax,ycmin
  REAL(KIND(1d0)),INTENT(OUT)::yi,wgt
  LOGICAL,INTENT(OUT)::lflag
  REAL(KIND(1d0))::lh,mti,ymin,ymax,ypmin,ypmax,ymmin,ymmax
  lflag=.TRUE.
  mti=DSQRT(mi**2+pti**2)
  IF(mti.LE.0d0)THEN
     lflag=.FALSE.
     yi=0d0
     wgt=0d0
     RETURN
  ENDIF
  lh=Q**2+mi**2-mres**2+2d0*pt*pti*DCOS(phi-phi0)
  lh=lh/(2d0*DSQRT(Q**2+pt**2)*mti)
  IF(lh.LT.1d0)THEN
     lflag=.FALSE.
     yi=0d0
     wgt=0d0
     RETURN
  ENDIF
  ymin=ACosh_p(lh)
  yi=SIGN(1d0,xr-0.5d0)*ymin+y0
  wgt=1d0
  RETURN
END SUBROUTINE ASSIGN_YNM1

END MODULE Colliders_PSI_2
