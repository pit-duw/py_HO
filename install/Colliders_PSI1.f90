MODULE Colliders_PSI_1
USE MC_VEGAS
USE Helac_Global
USE Constants
USE func_psi
USE Helac_SM_FeynRule
USE Structf_PDFs
USE Helac_master
USE Helac_Func_1
USE Cuts_Module
USE FO_plot
USE Decay_interface
USE reweight_xsec
USE KT_Clustering
IMPLICIT NONE
INTEGER::varnum,lnum,snum,nunwei
LOGICAL::ptQ,lunwei2
REAL(KIND(1d0))::EBMUP1,EBMUP2
INTEGER::NPRUP=0,lwmax
INTEGER::nprint
INTEGER,PARAMETER::maxprint=8
!LOGICAL::start_shower
SAVE
CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!  Hadronic Collider
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BB_1FS(rslt,itmxn,ncalln)
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
CALL ReadElem_logic('ktrw',ktrw)
IF(ktrw)THEN
   CALL ReadElem_integer("ktmeasure",ktmeasure)
   CALL ReadElem_real("Rparameter",RR)
ENDIF
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
   varnum = 1
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
ELSE
   ! MC over helicities with explicit helicity
   varnum=1
   varnum=varnum+1+NDecayIflh
   lnum=0
   snum=0
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
   CALL VEGAS(varnum,BB_1FS_fxn,vfes,sd,chi2a)
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
CALL VEGAS(varnum,BB_1FS_fxn,vfes,sd,chi2a)
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
        CALL VEGAS(varnum,BB_1FS_fxn,vfes,sd,chi2a,1)
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
  CALL VEGAS(varnum,BB_1FS_fxn,vfes,sd,chi2a,1)
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
  CALL VEGAS(varnum,BB_1FS_fxn,vfes,sd,chi2a,1)
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
END SUBROUTINE BB_1FS

FUNCTION BB_1FS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt ! it is the weight from VEGAS, 
! which will be used to generate unweighted events, wgt=VTOT*Vi/Ni=Vi/Ni since VTOT=1
REAL(KIND=DBL)::BB_1FS_fxn,sqs,m1r
REAL(KIND=DBL)::sq,m12
REAL(KIND=DBL)::temp1,ymax,ymin,Jy,y,w1,wr
REAL(KIND=DBL)::wme,wme0,wgtbr,ptp,pqp,pi,wsf,recmax=0,recmin=0,y1cup,y1clow
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
    IF(fixtarget)sqs=sqs/DSQRT(2d0)
    EBMUP1=ebeam(1)
    EBMUP2=ebeam(2)
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
    xFcut(3:20)=1d0
    xFcutlow(3:20)=-1d0
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
    IF(.NOT.fixtarget)THEN
       ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
    ELSE
       IF(ABS(ebeam(1)).GT.ABS(ebeam(2)))THEN
          fixtargetrev=.FALSE.
          ycollcm=DLOG(1d0+2d0*ABS(ebeam(1))/ABS(ebeam(2)))/2d0
          FT_M=ebeam(2)
          FT_E1=ebeam(1)
          sqs=DSQRT(sqs**2+FT_M**2)
       ELSE
          fixtargetrev=.TRUE.
          ycollcm=-DLOG(1d0+2d0*ABS(ebeam(2))/ABS(ebeam(1)))/2d0
          FT_M=ebeam(1)
          FT_E1=ebeam(2)
          sqs=DSQRT(sqs**2+FT_M**2)
       ENDIF
    ENDIF
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
        PRINT *,"Wrong in allocation of array in BB_1FS_fxn !STOP !"
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
          PRINT *,"y1cup = ",y1cup
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
IF(absrap)THEN
   y=(2d0*x(1)-1d0)*(ymax-ymin)+SIGN(1d0,x(1)-0.5d0)*ymin
ELSE
   y=x(1)*(ymax-ymin)+ymin
ENDIF
xp1=DSQRT(m12/sq)*DEXP(y) ! xp1,xp2 are same in collision frame and lab frame
xp2=DSQRT(m12/sq)*DEXP(-y)
exp1=xp1*ebeam(1)
exp2=xp2*ebeam(2)
IF(fixtarget)THEN
   IF(.NOT.fixtargetrev)THEN
      exp1=exp1+FT_M*xp1/2d0
      exp2=exp2/2d0
   ELSE
      exp1=exp1/2d0
      exp2=exp2+FT_M*xp2/2d0
   ENDIF
ENDIF
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
    BB_1FS_fxn=0d0
    RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
kk=2
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
!      BB_1FS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1.AND..NOT.ktrw)THEN
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
IF(ktrw.AND.ABS(wme).GT.0d0)THEN
   CALL ktreweight(wr)
   wme=wme*wr
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   BB_1FS_fxn=0d0
   RETURN
ENDIF
CALL strf_pdf(wsf)
IF(ABS(wsf).LE.0d0)THEN
   BB_1FS_fxn=0d0
   RETURN
ENDIF
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
BB_1FS_fxn=wme*wsf*LDMEwt*Jy*3.8938573d5
! only unweight events are generated -> LHE files                                                                                                                                                         
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------
   w1=BB_1FS_fxn*wgt ! multiply VEGAS weight
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=BB_1FS_fxn*wgt ! multiply VEGAS weight
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.BB_1FS_fxn)recmax=BB_1FS_fxn
IF(recmin.GT.BB_1FS_fxn)recmin=BB_1FS_fxn
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
END FUNCTION BB_1FS_fxn
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Cross Section of BoundState BoundState -> 2 Final States
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BB_2FS(rslt,itmxn,ncalln)
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
CALL ReadElem_logic('ktrw',ktrw)
IF(ktrw)THEN
   CALL ReadElem_integer("ktmeasure",ktmeasure)
   CALL ReadElem_real("Rparameter",RR)
ENDIF
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
	IF(ptQ)THEN
		varnum=2
	ELSE
		varnum=3
	ENDIF
        IF(lunwei2)varnum=varnum+1 ! the trivial phi should be used
	lnum=0
	snum=0
ELSEIF(iranhel.EQ.1)THEN
	IF(ptQ.AND..NOT.lunwei2)THEN
		varnum=2
	ELSE
		varnum=3
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
		varnum=2
	ELSE
		varnum=3
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
		varnum=2
	ELSE
		varnum=3
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
      varnum=2
   ELSE
      varnum=3
   ENDIF
   IF(lunwei2)varnum=varnum+1 ! the trivial phi should be used                                               
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
   CALL VEGAS(varnum,BB_2FS_fxn,vfes,sd,chi2a)
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
CALL VEGAS(varnum,BB_2FS_fxn,vfes,sd,chi2a)
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
	CALL VEGAS(varnum,BB_2FS_fxn,vfes,sd,chi2a,1)
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
!CALL VEGAS(3,BB_2FS_fxn,vfes,sd,chi2a,1)
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
  IF(plot_output)CALL initplot
  CALL VEGAS(varnum,BB_2FS_fxn,vfes,sd,chi2a,1)
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
  CALL VEGAS(varnum,BB_2FS_fxn,vfes,sd,chi2a,1)
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
END SUBROUTINE BB_2FS

FUNCTION BB_2FS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt
REAL(KIND=DBL)::BB_2FS_fxn,sqs,m1r,m2r,w1,wr
REAL(KIND=DBL)::pt,y1,y2,mt1,mt2,sq,m12,m22,phi
REAL(KIND=DBL)::temp1,temp2,temp21,temp3,temp4,Jpt2,Jy1,Jy2
REAL(KIND=DBL)::wme,wme0,wgtbr,ptp,pqp,pi,wsf,recmax=0,recmin=0,pt1c,maxpt1c,y1cup,y1clow
REAL(KIND=DBL)::exp1,exp2,ycollcm
REAL(KIND=DBL),DIMENSION(4)::pmom
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
INTEGER::init=0,kk,j,istat,nnn=0,nnntot=0,pdfnumpdf,icut,nwri,nwri_tot,nwarn,nwmax
SAVE init,sqs,m1r,m2r,m12,m22,sq,temp1,pi,lr,sr,pt,nnntot,ycollcm
SAVE nnn,recmax,recmin,pt1c,maxpt1c,y1cup,y1clow,nwri,nwri_tot,nwarn,nwmax
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
        IF(fixtarget)sqs=sqs/DSQRT(2d0)
        EBMUP1=ebeam(1)
        EBMUP2=ebeam(2)
        IF(.NOT.fixtarget)THEN
           ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
        ELSE
           IF(ABS(ebeam(1)).GT.ABS(ebeam(2)))THEN
              fixtargetrev=.FALSE.
              ycollcm=DLOG(1d0+2d0*ABS(ebeam(1))/ABS(ebeam(2)))/2d0
              FT_M=ebeam(2)
              FT_E1=ebeam(1)
              sqs=DSQRT(sqs**2+FT_M**2)
           ELSE
              fixtargetrev=.TRUE.
              ycollcm=-DLOG(1d0+2d0*ABS(ebeam(2))/ABS(ebeam(1)))/2d0
              FT_M=ebeam(1)
              FT_E1=ebeam(2)
              sqs=DSQRT(sqs**2+FT_M**2)
           ENDIF
        ENDIF
	!CALL ReadElem_real("energy",sqs)
        !EBMUP1=sqs/2d0
        !EBMUP2=sqs/2d0
	CALL ReadElem_real("minpt1c",pt1c)
        CALL ReadElem_real("maxpt1c",maxpt1c)
        IF(maxpt1c.LT.0d0)THEN
           maxpt1c=-1d0
        ELSE
           IF(maxpt1c.LE.pt1c)THEN
              WRITE(*,*)"ERROR: the first final state was cut off by pt cut (pt1c,maxpt1c) "
              STOP
           ENDIF
        ENDIF
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
        xFcut(3:20)=1d0
        xFcutlow(3:20)=-1d0
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
		PRINT *,"Wrong in allocation of array in BB_2FS_fxn !STOP !"
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
        IF(maxpt1c.GE.0d0)THEN
           temp1=MIN(temp1,maxpt1c)
        ENDIF
	IF(ptQ)THEN
		CALL ReadElem_real("Pt1",pt)
	ENDIF
	init=1
ENDIF
nnntot=nnntot+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The physical region of the integration
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(.NOT.ptQ)THEN
! total cross section
	pt=(temp1-pt1c)*x(3)+pt1c
	mt1=DSQRT(pt*pt+m1r*m1r)
	mt2=DSQRT(pt*pt+m2r*m2r)
        IF(absrap)THEN
           temp2=MIN(ACosh_p((sq+m12-m22)/(2d0*sqs*mt1)),y1cup)
           IF(temp2.LT.y1clow)THEN
              BB_2FS_fxn=0d0
              RETURN
           ENDIF
           y1=(2d0*x(1)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1)-0.5d0)*y1clow
        ELSE
           temp2=MIN(ACosh_p((sq+m12-m22)/(2d0*sqs*mt1)),y1cup-ycollcm)
           temp21=MAX(-ACosh_p((sq+m12-m22)/(2d0*sqs*mt1)),y1clow-ycollcm)
           IF(temp2.LT.temp21)THEN
              BB_2FS_fxn=0d0
              RETURN
           ENDIF
           y1=(temp2-temp21)*x(1)+temp21 ! in collision frame
        ENDIF
	temp3=DLOG((-DEXP(-y1)*mt1 +sqs)/mt2)
	temp4=temp3+DLOG((-DEXP(y1)*mt1 + sqs)/mt2)
	y2= -temp3 +temp4*x(2) ! in collision frame
!phi=x(4)*2*pi
!Jacobi Determinations
	Jpt2=2d0*pt*(temp1-pt1c)
        IF(absrap)THEN
           Jy1=2d0*(temp2-y1clow)
        ELSE
           Jy1=temp2-temp21
        ENDIF
	Jy2=temp4
ELSE
! pt distribution
	mt1=DSQRT(pt*pt+m1r*m1r)
	mt2=DSQRT(pt*pt+m2r*m2r)
        IF(absrap)THEN
           temp2=MIN(ACosh_p((sq+m12-m22)/(2d0*sqs*mt1)),y1cup)
           IF(temp2.LT.y1clow)THEN
              BB_2FS_fxn=0d0
              RETURN
           ENDIF
           y1=(2d0*x(1)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1)-0.5d0)*y1clow
        ELSE
           temp2=MIN(ACosh_p((sq+m12-m22)/(2d0*sqs*mt1)),y1cup-ycollcm)
           temp21=MAX(-ACosh_p((sq+m12-m22)/(2d0*sqs*mt1)),y1clow-ycollcm)
           IF(temp2.LT.temp21)THEN
              BB_2FS_fxn=0d0
              RETURN
           ENDIF
           y1=(temp2-temp21)*x(1)+temp21 ! in collision frame
        ENDIF
	temp3=DLOG((-DEXP(-y1)*mt1 +sqs)/mt2)
	temp4=temp3+DLOG((-DEXP(y1)*mt1 + sqs)/mt2)
	y2= -temp3 +temp4*x(2) ! in collision frame
!phi=x(4)*2*pi
!Jacobi Determinations
	Jpt2=2d0*pt
        IF(absrap)THEN
           Jy1=2d0*(temp2-y1clow)
        ELSE
           Jy1=(temp2-temp21)
        ENDIF
	Jy2=temp4
ENDIF
! The substitution of original variations
xp1=(DEXP(y1)*mt1 + DEXP(y2)*mt2)/sqs ! xp1,xp2 are same in collision frame and lab frame
xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2)/sqs
ehat=DSQRT(xp1*xp2*sq)
exp1=xp1*ebeam(1)
exp2=xp2*ebeam(2)
IF(fixtarget)THEN
   IF(.NOT.fixtargetrev)THEN
      exp1=exp1+FT_M*xp1/2d0
      exp2=exp2/2d0
   ELSE
      exp1=exp1/2d0
      exp2=exp2+FT_M*xp2/2d0
   ENDIF
ENDIF
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
      phi=2*pi*x(3)
   ELSE
      phi=2*pi*x(4)
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
   ! boost from collision frame to lab frame in polarization
   IF(.NOT.fixtarget)THEN
      pmom(3)=ebeam(1)-ebeam(2)
      pmom(4)=ebeam(1)+ebeam(2)
   ELSE
      IF(.NOT.fixtargetrev)THEN
         pmom(3)=ebeam(1)
         pmom(4)=ebeam(1)+ebeam(2)
      ELSE
         pmom(3)=-ebeam(2)
         pmom(4)=ebeam(1)+ebeam(2)
      ENDIF
   ENDIF
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
icut=1
IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
	CALL Cuts_HC(icut)
ELSEIF(COLL_TYPE.EQ.3)THEN
	CALL Cuts_epem(icut)
ELSE
	icut=1
ENDIF

IF(icut.EQ.0)THEN
	BB_2FS_fxn=0d0
	RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
IF(ptQ)THEN
	kk=3
ELSE
	kk=4
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
!      STOP
!      BB_2FS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1.AND..NOT.ktrw)THEN
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
IF(ktrw.AND.ABS(wme).GT.0d0)THEN
   CALL ktreweight(wr)
   wme=wme*wr
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   BB_2FS_fxn=0d0
   RETURN
ENDIF
CALL strf_pdf(wsf)
IF(ABS(wsf).LE.0d0)THEN
   BB_2FS_fxn=0d0
   RETURN
ENDIF
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
BB_2FS_fxn=wme*wsf*LDMEwt/(16d0*pi*ehat**4)&
           *xp1*xp2*Jpt2*Jy1*Jy2*3.8938573d5
! only unweight events are generated -> LHE files  
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------
   w1=BB_2FS_fxn*wgt ! multiply VEGAS weight
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=BB_2FS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.BB_2FS_fxn)recmax=BB_2FS_fxn
IF(recmin.GT.BB_2FS_fxn)recmin=BB_2FS_fxn
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
END FUNCTION BB_2FS_fxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Cross Section of BoundState BoundState -> 3 Final States
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BB_3FS(rslt,itmxn,ncalln)
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
CALL ReadElem_logic('ktrw',ktrw)
IF(ktrw)THEN
   CALL ReadElem_integer("ktmeasure",ktmeasure)
   CALL ReadElem_real("Rparameter",RR)
ENDIF
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
	IF(ptQ)THEN
		varnum=5
	ELSE
		varnum=6
	ENDIF
        IF(lunwei2)varnum=varnum+1
	lnum=0
	snum=0
ELSEIF(iranhel.EQ.1)THEN
	IF(ptQ)THEN
		varnum=5
	ELSE
		varnum=6
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
		varnum=5
	ELSE
		varnum=6
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
		varnum=5
	ELSE
		varnum=6
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
      varnum=5
   ELSE
      varnum=6
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
   CALL VEGAS(varnum,BB_3FS_fxn,vfes,sd,chi2a)
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
CALL VEGAS(varnum,BB_3FS_fxn,vfes,sd,chi2a)
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
	CALL VEGAS(varnum,BB_3FS_fxn,vfes,sd,chi2a,1)
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
  CALL VEGAS(varnum,BB_3FS_fxn,vfes,sd,chi2a,1)
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
  CALL VEGAS(varnum,BB_3FS_fxn,vfes,sd,chi2a,1)
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
END SUBROUTINE BB_3FS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the external function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION BB_3FS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt
REAL(KIND=DBL)::BB_3FS_fxn,m1r,m2r,m3r,sqs,phi0,w1,wr
REAL(KIND=DBL)::pt1,pt12,pt2,y1,y2,y3,phi,mt1,mt2,mt3,sq,m12,m22,m32,aa,ff,ma,hh,cosphi,cosphi2
REAL(KIND=DBL)::temp1,temp2,temp21,temp3,temp4,temp5,temp6,Jpt22,Jpt12,Jy1,Jy2,Jy3,Jphi,&
                lh         ! in 2-> 3, we use t=p2k2 u=p1k2
INTEGER::init=0,kk,j,icut,istat,nnn=0,nnntot=0,pdfnumpdf,nwri,nwri_tot,nwmax,nwarn
REAL(KIND=DBL)::exp1,exp2,ycollcm
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
REAL(KIND=DBL)::ptp,pqp,wme,wme0,wgtbr,pi,wsf,y1cup,y1clow,pt1c,maxpt1c,recmax=0,recmin=0
REAL(KIND=DBL),DIMENSION(4)::pmom
SAVE init,sqs,sq,m1r,m2r,m3r,m12,m22,m32,temp1,pi,y1cup,y1clow,pt1c,maxpt1c,lr,sr,pt1,pt12
SAVE nnn,recmax,recmin,nwri,nwri_tot,nwmax,nwarn,nnntot,ycollcm
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
           absrap=.FALSE.
           labeqcoll=.FALSE.
        ENDIF
        sqs=2d0*DSQRT(ebeam(1)*ebeam(2)) ! we always neglect the mass of initial states
        IF(fixtarget)sqs=sqs/DSQRT(2d0)
        EBMUP1=ebeam(1)
        EBMUP2=ebeam(2)
        IF(.NOT.fixtarget)THEN
           ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
        ELSE
           IF(ABS(ebeam(1)).GT.ABS(ebeam(2)))THEN
              fixtargetrev=.FALSE.
              ycollcm=DLOG(1d0+2d0*ABS(ebeam(1))/ABS(ebeam(2)))/2d0
              FT_M=ebeam(2)
              FT_E1=ebeam(1)
              sqs=DSQRT(sqs**2+FT_M**2)
           ELSE
              fixtargetrev=.TRUE.
              ycollcm=-DLOG(1d0+2d0*ABS(ebeam(2))/ABS(ebeam(1)))/2d0
              FT_M=ebeam(1)
              FT_E1=ebeam(2)
              sqs=DSQRT(sqs**2+FT_M**2)
           ENDIF
        ENDIF
	CALL ReadElem_real("minpt1c",pt1c)
        CALL ReadElem_real("maxpt1c",maxpt1c)
        IF(maxpt1c.LT.0d0)THEN
           maxpt1c=-1d0
        ELSE
           IF(maxpt1c.LE.pt1c)THEN
              WRITE(*,*)"ERROR: the first final state was cut off by pt cut (pt1c,maxpt1c) "
              STOP
           ENDIF
        ENDIF
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
        xFcut(3:20)=1d0
        xFcutlow(3:20)=-1d0
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
		PRINT *,"Wrong in allocation of array in BB_3FS_fxn !STOP !"
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
        IF(maxpt1c.GE.0d0)THEN
           temp1=MIN(temp1,maxpt1c)
        ENDIF
	IF(ptQ)THEN
		CALL ReadElem_real("Pt1",pt1)
		pt12=pt1*pt1
	ENDIF
	init=1
ENDIF
nnntot=nnntot+1
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
              BB_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(2d0*x(2)-1d0)*(temp2-y1clow)+SIGN(1d0,x(2)-0.5d0)*y1clow
        ELSE
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup-ycollcm)
           temp21=MAX(-ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1clow-ycollcm)
           IF(temp2.LT.temp21)THEN
              BB_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(temp2-temp21)*x(2)+temp21 ! in collision frame
        ENDIF
	aa=sq+mt1*mt1-2d0*sqs*mt1*DCOSH(y1)
	ff=m32+pt12
	ma=DMAX1(-lamda_func(aa,m22,ff)/(4d0*m22*pt12),0d0)
	ma=DSQRT(ma)
	phi=unit_step(1d0/2d0-x(3))*(4d0*DACOS(ma)*x(3)-DACOS(ma))&
     +unit_step(x(3)-1d0/2d0)*(2d0*(pi-DACOS(-ma))*(2d0*x(3)-1d0)+DACOS(-ma)) ! in collision frame
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
		BB_3FS_fxn=0d0
		RETURN
	ENDIF
	temp4=ACosh_p(lh)
	y2=(2d0*x(5)-1d0)*temp4+DLOG((sqs-mt1*DEXP(y1))/(sqs-mt1*DEXP(-y1)))/2d0 ! in collision frame
	temp5=DLOG((sqs-DEXP(-y1)*mt1-DEXP(-y2)*mt2)/mt3)
	temp6=temp5+DLOG((sqs-DEXP(y1)*mt1-DEXP(y2)*mt2)/mt3)
	y3= -temp5 +temp6*x(6) ! in collision frame
!Jacobi Determinations
	Jpt12=2d0*pt1*(temp1-pt1c)
        IF(absrap)THEN
           Jy1=2d0*(temp2-y1clow)
        ELSE
           Jy1=temp2-temp21
        ENDIF
	Jphi=unit_step(1d0/2d0-x(3))*(4d0*DACOS(ma))&
     +unit_step(x(3)-1d0/2d0)*(4d0*(pi-DACOS(-ma)))
	Jpt22=2d0*pt2*temp3
	Jy2=2d0*temp4
	Jy3=temp6
ELSE
! pt distribution
	mt1=DSQRT(pt12+m12)
        IF(absrap)THEN
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup)
           IF(temp2.LT.y1clow)THEN
              BB_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(2d0*x(1)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1)-0.5d0)*y1clow
        ELSE
           temp2=MIN(ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1cup-ycollcm)
           temp21=MAX(-ACosh_p((sq+m12-(m2r+m3r)**2)/(2d0*sqs*mt1)),y1clow-ycollcm)
           IF(temp2.LT.temp21)THEN
              BB_3FS_fxn=0d0
              RETURN
           ENDIF
           y1=(temp2-temp21)*x(1)+temp21 ! in collision frame           
        ENDIF
	aa=sq+mt1*mt1-2d0*sqs*mt1*DCOSH(y1)
	ff=m32+pt12
	ma=DMAX1(-lamda_func(aa,m22,ff)/(4d0*m22*pt12),0d0)
	ma=DSQRT(ma)
	phi=unit_step(1d0/2d0-x(2))*(4d0*DACOS(ma)*x(2)-DACOS(ma))&
     +unit_step(x(2)-1d0/2d0)*(2d0*(pi-DACOS(-ma))*(2d0*x(2)-1d0)+DACOS(-ma)) ! + pi ?
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
		BB_3FS_fxn=0d0
		RETURN
	ENDIF
	temp4=ACosh_p(lh)
	y2=(2d0*x(4)-1d0)*temp4+DLOG((sqs-mt1*DEXP(y1))/(sqs-mt1*DEXP(-y1)))/2d0 ! in collision frame
	temp5=DLOG((sqs-DEXP(-y1)*mt1-DEXP(-y2)*mt2)/mt3)
	temp6=temp5+DLOG((sqs-DEXP(y1)*mt1-DEXP(y2)*mt2)/mt3)
	y3= -temp5 +temp6*x(5) ! in collision frame
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
	Jy2=2d0*temp4
	Jy3=temp6
ENDIF
! The substitution of original variations
xp1=(DEXP(y1)*mt1 + DEXP(y2)*mt2+DEXP(y3)*mt3)/sqs ! xp1,xp2 is the same in collision frame and lab frame
xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2+DEXP(-y3)*mt3)/sqs
ehat=DSQRT(xp1*xp2*sq)
exp1=xp1*ebeam(1)
exp2=xp2*ebeam(2)
IF(fixtarget)THEN
   IF(.NOT.fixtargetrev)THEN
      exp1=exp1+FT_M*xp1/2d0
      exp2=exp2/2d0
   ELSE
      exp1=exp1/2d0
      exp2=exp2+FT_M*xp2/2d0
   ENDIF
ENDIF
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
   WRITE(*,*)"ERROR:Please use BB_nFS for unweighted events !"
   STOP
   IF(ptQ)THEN
      phi0=2*pi*x(6)
   ELSE
      phi0=2*pi*x(7)
   ENDIF
ELSE
   phi0=0d0
ENDIF
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
   IF(.NOT.fixtarget)THEN
      pmom(3)=ebeam(1)-ebeam(2)
      pmom(4)=ebeam(1)+ebeam(2)
   ELSE
      IF(.NOT.fixtargetrev)THEN
         pmom(3)=ebeam(1)
         pmom(4)=ebeam(1)+ebeam(2)
      ELSE
         pmom(3)=-ebeam(2)
         pmom(4)=ebeam(1)+ebeam(2)
      ENDIF
   ENDIF
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
	CALL Cuts_epem(icut)
ELSE
	icut=1
ENDIF

IF(icut.EQ.0)THEN
	BB_3FS_fxn=0d0
	RETURN
ENDIF
! Reloading Feynman Rules
IF(irun.EQ.1)CALL Helac_FeynRule_SM()
! set the random number for random helicity method
IF(ptQ)THEN
	kk=6
ELSE
	kk=7
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
!      BB_3FS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1.AND..NOT.ktrw)THEN
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
IF(ktrw.AND.ABS(wme).GT.0d0)THEN
   CALL ktreweight(wr)
   wme=wme*wr
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   BB_3FS_fxn=0d0
   RETURN
ENDIF
! Generate weight for pdf
CALL strf_pdf(wsf)
IF(ABS(wsf).LE.0d0)THEN
   BB_3FS_fxn=0d0
   RETURN
ENDIF
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
BB_3FS_fxn=wme*wsf*LDMEwt/(512d0*(pi**4)*ehat**4)&
           *xp1*xp2*Jpt12*Jpt22*Jy1*Jy2*Jy3*Jphi*3.8938573d5
! only unweight events are generated -> LHE files                                                                                                                                                           
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------   
   w1=BB_3FS_fxn*wgt ! multiply VEGAS weight  
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=BB_3FS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.BB_3FS_fxn)recmax=BB_3FS_fxn
IF(recmin.GT.BB_3FS_fxn)recmin=BB_3FS_fxn
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
END FUNCTION BB_3FS_fxn

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The Cross Section of BoundState BoundState -> n Final States with n>=2                  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE BB_nFS(rslt,itmxn,ncalln)
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
CALL ReadElem_logic('ktrw',ktrw)
IF(ktrw)THEN
   CALL ReadElem_integer("ktmeasure",ktmeasure)
   CALL ReadElem_real("Rparameter",RR)
ENDIF
IF(.NOT.MCoHelicity)THEN
IF(iranhel.EQ.0)THEN
   IF(ptQ)THEN
      varnum=3*(nhad-3)-1
   ELSE
      varnum=3*(nhad-3)
   ENDIF
   IF(lunwei2)varnum=varnum+1
   lnum=0
   snum=0
ELSEIF(iranhel.EQ.1)THEN
   IF(ptQ)THEN
      varnum=3*(nhad-3)-1
   ELSE
      varnum=3*(nhad-3)
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
      varnum=3*(nhad-3)-1
   ELSE
      varnum=3*(nhad-3)
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
      varnum=3*(nhad-3)-1
   ELSE
      varnum=3*(nhad-3)
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
      varnum=3*(nhad-3)-1
   ELSE
      varnum=3*(nhad-3)
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
   CALL VEGAS(varnum,BB_nFS_fxn,vfes,sd,chi2a)
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
CALL VEGAS(varnum,BB_nFS_fxn,vfes,sd,chi2a)
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
   CALL VEGAS(varnum,BB_nFS_fxn,vfes,sd,chi2a,1)
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
   CALL VEGAS(varnum,BB_nFS_fxn,vfes,sd,chi2a,1)
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
  CALL VEGAS(varnum,BB_nFS_fxn,vfes,sd,chi2a,1)
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
END SUBROUTINE BB_nFS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Define the external function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
FUNCTION BB_nFS_fxn(x,wgt)
REAL(KIND=DBL),DIMENSION(varnum),INTENT(IN)::x
REAL(KIND=DBL),INTENT(IN)::wgt
REAL(KIND=DBL)::BB_nFS_fxn,sqs,sq
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::mr,ptcut,ymincut,ymaxcut
REAL(KIND=DBL),DIMENSION(nhad,1:4)::phadmom
REAL(KIND=DBL)::pt0,pti,phi,phi0,y0,yi,yn,Q,mi,mres,jac,jactemp,mtot,pt1,pt12,mti
REAL(KIND=DBL)::rlogmax,rlogmin,w1,wr
LOGICAL::lflag=.TRUE.
INTEGER::init=0,kk,kk2,ii,j,icut,istat,nnn=0,nnntot=0,pdfnumpdf,nwri,nwri_tot,nwarn,nwmax
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::lr,sr
REAL(KIND=DBL)::ptp,pqp,wme,wme0,wgtbr,pi,wsf,y1cup,y1clow,pt1c,maxpt1c,recmax=0,recmin=0
REAL(KIND=DBL),DIMENSION(4)::pmom
REAL(KIND=DBL)::exp1,exp2,ycollcm
SAVE init,sqs,sq,pi,y1cup,y1clow,pt1c,maxpt1c,lr,sr,pt1,pt12,mr,mtot,nnntot,ycollcm
SAVE nnn,recmax,recmin,ptcut,ymincut,ymaxcut,nwri,nwri_tot,nwarn,nwmax
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
        IF(fixtarget)sqs=sqs/DSQRT(2d0)
        EBMUP1=ebeam(1)
        EBMUP2=ebeam(2)
        IF(.NOT.fixtarget)THEN
           ycollcm=DLOG(ABS(ebeam(1))/ABS(ebeam(2)))/2d0
        ELSE
           IF(ABS(ebeam(1)).GT.ABS(ebeam(2)))THEN
              fixtargetrev=.FALSE.
              ycollcm=DLOG(1d0+2d0*ABS(ebeam(1))/ABS(ebeam(2)))/2d0
              FT_M=ebeam(2)
              FT_E1=ebeam(1)
              sqs=DSQRT(sqs**2+FT_M**2)
           ELSE
              fixtargetrev=.TRUE.
              ycollcm=-DLOG(1d0+2d0*ABS(ebeam(2))/ABS(ebeam(1)))/2d0
              FT_M=ebeam(1)
              FT_E1=ebeam(2)
              sqs=DSQRT(sqs**2+FT_M**2)
           ENDIF
        ENDIF
        CALL ReadElem_real("minpt1c",pt1c)
!        CALL ReadElem_real("maxpt1c",maxpt1c)
!        IF(maxpt1c.LT.0d0)THEN
!           maxpt1c=-1d0
!        ELSE
!           IF(maxpt1c.LE.pt1c)THEN
!              WRITE(*,*)"ERROR: the first final state was cut off by pt cut (pt1c,maxpt1c)"
!              STOP
!           ENDIF
!        ENDIF
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
        xFcut(3:20)=1d0
        xFcutlow(3:20)=-1d0
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
                PRINT *,"Wrong in allocation of array in BB_nFS_fxn !STOP !"
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
           IF(absrap)THEN
              ymincut(3)=MAX(y1clow,ymincut(3))
              ymaxcut(3)=MIN(y1cup,ymaxcut(3))
           ELSE
              ymincut(3)=MAX(y1clow-ycollcm,ymincut(3))
              ymaxcut(3)=MIN(y1cup-ycollcm,ymaxcut(3))
           ENDIF
        ENDIF
        IF(ymincut(3).GE.ymaxcut(3))THEN
           WRITE(*,*)"All events are cutted off !"
           WRITE(*,*)"k, ycmin, ycmax"
           WRITE(*,*)3,ymincut(3),ymaxcut(3)
           WRITE(*,*)"BYE !"
           STOP
        ENDIF
        sq=sqs*sqs
        IF(ptQ)THEN
           CALL ReadElem_real("Pt1",pt1)
           pt12=pt1*pt1
        ENDIF
        init=1
ENDIF
nnntot=nnntot+1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The physical region of the integration in collision frame
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
      BB_nFS_fxn=0d0
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
      BB_nFS_fxn=0d0
      RETURN
   ENDIF
   jac=jac*jactemp*2d0*pti ! from d pti^2 -> d pti
   ! generate yi
   CALL generate_yi(Q,mi,mres,phi0,phi,pt0,pti,y0,x(kk),ymaxcut(ii),ymincut(ii),&
        yi,jactemp,lflag)
   kk=kk+1
   IF(.NOT.lflag)THEN
      BB_nFS_fxn=0d0
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
   xp1=xp1+DEXP(yi)*mti
   xp2=xp2+DEXP(-yi)*mti
ENDDO
! generate yn
IF(rlogmax.LE.0d0.OR.rlogmin.LE.0d0)THEN
   BB_nFS_fxn=0d0
   RETURN
ENDIF
pti=transverse2(pmom)
mi=mr(nhad)
mti=DSQRT(pti**2+mi**2)
rlogmax=DLOG(rlogmax/mti)
rlogmin=-DLOG(rlogmin/mti)
CALL generate_yn(rlogmax,rlogmin,x(kk),ymaxcut(nhad),ymincut(nhad),yn,jactemp,lflag)
kk=kk+1
IF(.NOT.lflag)THEN
   BB_nFS_fxn=0d0
   RETURN
ENDIF
jac=jac*jactemp
!phadmom(nhad,3)=mti*DSINH(yn)
!phadmom(nhad,4)=mti*DCOSH(yn)
xp1=xp1+DEXP(yn)*mti
xp2=xp2+DEXP(-yn)*mti
! The substitution of original variations
xp1=xp1/sqs
xp2=xp2/sqs
IF(xp1.GT.1d0.OR.xp1.LE.0d0.OR.xp2.GT.1d0.OR.xp2.LE.0d0)THEN
   BB_nFS_fxn=0d0
   RETURN
ENDIF
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
   ! boost from collision frame to lab frame in polarisation
   IF(.NOT.fixtarget)THEN
      pmom(3)=ebeam(1)-ebeam(2)
      pmom(4)=ebeam(1)+ebeam(2)
   ELSE
      IF(.NOT.fixtargetrev)THEN
         pmom(3)=ebeam(1)
         pmom(4)=ebeam(1)+ebeam(2)
      ELSE
         pmom(3)=-ebeam(2)
         pmom(4)=ebeam(1)+ebeam(2)
      ENDIF
   ENDIF
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
   CALL Cuts_epem(icut)
ELSE
   icut=1
ENDIF
IF(icut.EQ.0)THEN
   BB_nFS_fxn=0d0
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
      Decayran(j-kk)=x(j)
   ENDDO
ENDIF
!IF(NDecayChains.GT.0)THEN
!   CALL HO_Decay(wgtbr)
!   CALL cuts_Decay(icut)
!   IF(icut.EQ.0)THEN
!      BB_nFS_fxn=0d0
!      RETURN
!   ENDIF
!ENDIF
! Generate the square of matrix element
CALL Helac_master_f(wme)
GLOBALINIT_master=1
IF(reweight_scale.AND.nnn.LE.1.AND..NOT.ktrw)THEN
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
IF(ktrw.AND.ABS(wme).GT.0d0)THEN
   CALL ktreweight(wr)
   wme=wme*wr
ENDIF
IF(NDecayChains.GT.0)wme=wme*weight_br
IF(ABS(wme).LE.0d0)THEN
   BB_nFS_fxn=0d0
   RETURN
ENDIF
! Generate weight for pdf
CALL strf_pdf(wsf)
IF(ABS(wsf).LE.0d0)THEN
   BB_nFS_fxn=0d0
   RETURN
ENDIF
IF(reweight_scale)THEN
   CALL reweight_xsec_scale(wsf)
ENDIF
IF(reweight_pdf)THEN
   CALL reweight_xsec_pdf(wsf)
ENDIF
BB_nFS_fxn=wme*wsf*LDMEwt/(4d0*(4d0*pi)**(2*nhad-7)*(2*pi)**(nhad-3)*ehat**4)&
           *xp1*xp2*Jac*3.8938573d5
! only unweight events are generated -> LHE files 
IF(lunwei2.AND.lwmax.GT.0)THEN
! --------------   
   w1=BB_nFS_fxn*wgt ! multiply VEGAS weight  
   CALL unwei_procedure(w1,nwri,nwmax,nwarn)
ENDIF
IF(plot_output)THEN
   w1=BB_nFS_fxn*wgt
   CALL outfun(w1)
ENDIF
nnn=nnn+1
IF(recmax.LT.BB_nFS_fxn)recmax=BB_nFS_fxn
IF(recmin.GT.BB_nFS_fxn)recmin=BB_nFS_fxn
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
END FUNCTION BB_nFS_fxn

SUBROUTINE generate_phi(itype,Q,mi,mres,phi0,pt,xr,phi,wgt,lflag)
  IMPLICIT NONE
  INTEGER,INTENT(IN)::itype
  REAL(KIND(1d0)),INTENT(IN)::Q,pt,mi,mres,phi0,xr
  REAL(KIND(1d0)),INTENT(OUT)::phi,wgt
  LOGICAL,INTENT(OUT)::lflag
  REAL(KIND(1d0)),PARAMETER::pi=3.141592653589793d0
  REAL(KIND(1d0))::cosphimin,acos1,acos2
  lflag=.TRUE.
  IF(Q-mi-mres.LT.0d0)THEN
     lflag=.FALSE.
     wgt=0d0
     phi=0d0
     RETURN
  ENDIF
  IF(itype.EQ.0)THEN
     ! i<n-1, Q -> ki+Q2, mres<Q2<(Q-mi)
     IF(mi*pt.LE.0d0)THEN
        wgt=2*pi
        phi=wgt*xr
     ELSE
        cosphimin=((Q-mi)**2-mres**2)*Q/(2d0*mi*pt**2)
        IF(cosphimin.GT.2d0)THEN
           wgt=2*pi
           phi=wgt*xr
        ELSE
           cosphimin=1d0-cosphimin
           wgt=DACOS(cosphimin)
           phi=wgt*(2d0*xr-1d0)+phi0
           wgt=2d0*wgt
        ENDIF
     ENDIF
  ELSEIF(itype.EQ.1)THEN
     ! i=n-1, Q -> kn-1+kn, mres =mn
     IF(pt*mi.LE.0d0)THEN
        wgt=2*pi
        phi=wgt*xr
     ELSE
        cosphimin=DMAX1(-lamda_func(Q**2+pt**2,mi**2,mres**2+pt**2)/(4d0*mi**2*pt**2),0d0)
        cosphimin=DSQRT(cosphimin)
        acos1=DACOS(cosphimin)
        acos2=DACOS(-cosphimin)
        phi=unit_step(0.5d0-xr)*(4d0*acos1*xr-acos1)&
             +unit_step(xr-0.5d0)*(2d0*(pi-acos2)*(2d0*xr-1d0)+acos2)+phi0
        wgt=unit_step(0.5d0-xr)*(4d0*acos1)&
             +unit_step(xr-0.5d0)*(4d0*(pi-acos2))
     ENDIF
  ELSE
     WRITE(*,*)"Wrong itype(0 or 1) in generate_phi ",itype
     STOP
  ENDIF
  RETURN
END SUBROUTINE generate_phi

SUBROUTINE generate_pt(Q,mi,mres,phi0,phi,pt,xr,ptc,pti,wgt,lflag)
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::Q,pt,mi,mres,phi0,phi,xr,ptc ! ptc is the pt cut
  REAL(KIND(1d0)),INTENT(OUT)::pti,wgt
  LOGICAL,INTENT(OUT)::lflag
  REAL(KIND(1d0))::a,b,c,cosphi,delta,ptmax,ptmin
  lflag=.TRUE.
!  IF(Q-mi-mres.LT.0d0)THEN
!     lflag=.FALSE.
!     pti=0d0
!     wgt=0d0
!     RETURN
!  ENDIF
  cosphi=DCOS(phi-phi0)
  a=4d0*(Q**2+pt**2-cosphi**2*pt**2)
  b=-4d0*cosphi*pt*(Q**2+mi**2-mres**2)
  c=-lamda_func(Q**2+pt**2,mi**2,mres**2+pt**2)
  delta=b**2-4d0*a*c
  IF(delta.LT.0d0)THEN
     lflag=.FALSE.
     pti=0d0
     wgt=0d0
     RETURN
  ENDIF
  ptmax=(-b+DSQRT(delta))/(2d0*a)
  ptmin=DMAX1((-b-DSQRT(delta))/(2d0*a),0d0)
  ptmin=DMAX1(ptmin,ptc) ! cut impose
  IF(ptmax.LT.ptmin)THEN
     lflag=.FALSE.
     pti=0d0
     wgt=0d0
     RETURN
  ENDIF
  wgt=ptmax-ptmin
  pti=wgt*xr+ptmin
  RETURN
END SUBROUTINE generate_pt

SUBROUTINE generate_yi(Q,mi,mres,phi0,phi,pt,pti,y0,xr,ycmax,ycmin,&
     yi,wgt,lflag)
! generate yi for i<=n-1, Q -> ki+Q2, mres<Q2<(Q-mi) and Q2=mn**2
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
  ymax=ymin+y0
  ymin=-ymin+y0
  IF(absrap)THEN
     IF(ymax.LE.0d0)THEN
        ! MAX(ymin,-ABS(ycmax))<yi<MIN(ymax,-ABS(ycmin))
        ymmin=DMAX1(ymin,-ABS(ycmax))
        ymmax=DMIN1(ymax,-ABS(ycmin))
        IF(ymmax.LT.ymmin)THEN
           lflag=.FALSE.
           yi=0d0
           wgt=0d0
           RETURN
        ENDIF
        wgt=ymmax-ymmin
        yi=wgt*xr+ymmin
     ELSEIF(ymin.GE.0d0)THEN
        ! MAX(ymin,ABS(ycmin))<yi<MIN(ymax,ABS(ycmax))
        ypmin=DMAX1(ymin,ABS(ycmin))
        ypmax=DMIN1(ymax,ABS(ycmax))
        IF(ypmax.LT.ypmin)THEN
           lflag=.FALSE.
           yi=0d0
           wgt=0d0
           RETURN
        ENDIF
        wgt=ypmax-ypmin
        yi=wgt*xr+ypmin
     ELSE ! ymin<0,ymax>0
        ! MAX(ymin,-ABS(ycmax))<yi<MIN(0d0,-ABS(ycmin)) or
        ! MAX(0d0,ABS(ycmin))<yi<MIN(ymax,ABS(ycmax))
        ypmin=DMAX1(0d0,ABS(ycmin))
        ypmax=DMIN1(ymax,ABS(ycmax))
        ymmin=DMAX1(ymin,-ABS(ycmax))
        ymmax=DMIN1(0d0,-ABS(ycmin))
        IF(ymmax.LT.ymmin.AND.ypmax.LT.ypmin)THEN
           lflag=.FALSE.
           yi=0d0
           wgt=0d0
           RETURN
        ENDIF
        IF(ymmax.LT.ymmin)THEN
           wgt=ypmax-ypmin
           yi=wgt*xr+ypmin
        ELSEIF(ypmax.LT.ypmin)THEN
           wgt=ymmax-ymmin
           yi=wgt*xr+ymmin
        ELSE
           yi=unit_step(0.5d0-xr)*(2d0*xr*(ymmax-ymmin)+ymmin)&
                +unit_step(xr-0.5d0)*((2d0*xr-1d0)*(ypmax-ypmin)+ypmin)
           wgt=unit_step(0.5d0-xr)*2d0*(ymmax-ymmin)&
                +unit_step(xr-0.5d0)*2d0*(ypmax-ypmin)
        ENDIF
     ENDIF
  ELSE
     ymax=MIN(ymax,ycmax)
     ymin=MAX(ymin,ycmin)
     IF(ymax.LT.ymin)THEN
        lflag=.FALSE.
        yi=0d0
        wgt=0d0
        RETURN
     ENDIF
     wgt=ymax-ymin
     yi=wgt*xr+ymin
  ENDIF
  RETURN
END SUBROUTINE GENERATE_YI

SUBROUTINE generate_yn(rlogmax,rlogmin,xr,ycmax,ycmin,yn,wgt,lflag)
! generate yn for Q -> kn-1,kn within rlogmin<yn<rlogmax
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::rlogmax,rlogmin,xr,ycmax,ycmin ! ycmax,ycmin is the cut of yn
  REAL(KIND(1d0)),INTENT(OUT)::yn,wgt
  LOGICAL,INTENT(OUT)::lflag
  REAL(KIND(1d0))::ypmin,ypmax,ymmin,ymmax
  lflag=.TRUE.
  IF(rlogmax.LT.rlogmin.OR.ycmax.LT.ycmin)THEN
     lflag=.FALSE.
     yn=0d0
     wgt=0d0
     RETURN
  ENDIF
  IF(absrap)THEN
     IF(rlogmax.LE.0d0)THEN
        ! MAX(rlogmin,-ABS(ycmax))<yn<MIN(rlogmax,-ABS(ycmin))
        ymmin=DMAX1(rlogmin,-ABS(ycmax))
        ymmax=DMIN1(rlogmax,-ABS(ycmin))
        IF(ymmax.LT.ymmin)THEN
           lflag=.FALSE.
           yn=0d0
           wgt=0d0
           RETURN
        ENDIF
        wgt=ymmax-ymmin
        yn=wgt*xr+ymmin
     ELSEIF(rlogmin.GE.0d0)THEN
        ! MAX(rlogmin,ABS(ycmin))<yn<MIN(rlogmax,ABS(ycmax))
        ypmin=DMAX1(rlogmin,ABS(ycmin))
        ypmax=DMIN1(rlogmax,ABS(ycmax))
        IF(ypmax.LT.ypmin)THEN
           lflag=.FALSE.
           yn=0d0
           wgt=0d0
           RETURN
        ENDIF
        wgt=ypmax-ypmin
        yn=wgt*xr+ypmin
     ELSE ! rlogmin<0 and rlogmax > 0
        ! MAX(rlogmin,-ABS(ycmax))<yn<MIN(0d0,-ABS(ycmin)) or
        ! MAX(0d0,ABS(ycmin))<yn<MIN(rlogmax,ABS(ycmax))
        ypmin=DMAX1(ABS(ycmin),0d0)
        ypmax=DMIN1(rlogmax,ABS(ycmax))
        ymmin=DMAX1(rlogmin,-ABS(ycmax))
        ymmax=DMIN1(0d0,-ABS(ycmin))
        IF(ymmax.LT.ymmin.AND.ypmax.LT.ypmin)THEN
           lflag=.FALSE.
           yn=0d0
           wgt=0d0
           RETURN
        ENDIF
        IF(ymmax.LT.ymmin)THEN
           wgt=ypmax-ypmin
           yn=wgt*xr+ypmin
        ELSEIF(ypmax.LT.ypmin)THEN
           wgt=ymmax-ymmin
           yn=wgt*xr+ymmin
        ELSE
           yn=unit_step(0.5d0-xr)*(2d0*xr*(ymmax-ymmin)+ymmin)&
                +unit_step(xr-0.5d0)*((2d0*xr-1d0)*(ypmax-ypmin)+ypmin)
           wgt=unit_step(0.5d0-xr)*2d0*(ymmax-ymmin)&
                +unit_step(xr-0.5d0)*2d0*(ypmax-ypmin)
        ENDIF
     ENDIF
  ELSE
     ypmin=DMAX1(rlogmin,ycmin)
     ypmax=DMIN1(rlogmax,ycmax)
     IF(ypmin.GT.ypmax)THEN
        lflag=.FALSE.
        yn=0d0
        wgt=0d0
        RETURN
     ENDIF
     wgt=ypmax-ypmin
     yn=wgt*xr+ypmin
  ENDIF
  RETURN
END SUBROUTINE generate_yn

SUBROUTINE unwei_procedure(w1,nwri,nwmax,nwarn)
  USE QEDPS_interface
  USE Decay_interface
  IMPLICIT NONE
  REAL(KIND(1d0)),INTENT(IN)::w1
  INTEGER,INTENT(INOUT)::nwri,nwmax,nwarn
  REAL(KIND(1d0)),DIMENSION(1)::ranr
  REAL(KIND=DBL)::vtime,vspin,scale1,scalup,xwgtup,px,py,pz,p0,pmass,umax,umax1
  INTEGER::init=0,i
  INTEGER::idup,idprup,istup,imothup1,imothup2,icol1,icol2
  LOGICAL::llwri
  SAVE init,umax,umax1
  IF(init.EQ.0)THEN
     umax=0d0
     umax1=0d0
!     NPRUP=0
     init=1
  ENDIF

  IF(lwmax.EQ.1)THEN
     IF(umax.LT.w1)THEN
        umax=w1
     ENDIF
!     IF(MOD(nwmax,10000).EQ.0)PRINT *, umax,umax1
     nwmax=nwmax+1
     umax1=umax
  ENDIF
  IF(lwmax.EQ.2)THEN
     IF(w1.GT.umax1)THEN
        umax1=w1
        WRITE(*,*)'WARNING:umax1,umax',umax1,umax
     ENDIF
     llwri=.FALSE.
     CALL RANDA(1,ranr)
     IF(umax*ranr(1).LT.w1)llwri=.TRUE.
     IF(umax.LT.w1)nwarn=nwarn+1
     IF(llwri)THEN
        nwri=nwri+1
        IF(NDecayChains.GT.0)THEN
           CALL unwei_writer_Decay
           Nevents=nwri
           RETURN
        ENDIF
        IF(emep_ISR_shower.EQ.1)THEN
           CALL unwei_writer_QEDPS
           Nevents=nwri
           RETURN
        ENDIF
        idprup=81 ! id for the process                                                                                                                                                                     
        xwgtup=1 !w1*10**3 !1                                                                                                                                                                              
        CALL qcdscale(scale1)
        scalup=scale1
        WRITE(nunit3)nhad,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
        IF(imode.EQ.0.AND.NDecayChains.EQ.0)CALL cmstolab2() ! boost to lab                                                                                                                                                      
        DO i=1,nhad
           idup=pdgt(iflh(i))
           istup=-ioh(i)
           imothup1=0
           imothup2=0
           IF(i.GT.2.AND.i.LE.nhad)THEN
              imothup1=1
              imothup2=2
           ENDIF
           IF(ioh(i).EQ.1)THEN
              icol1=icol_un(i,1)+100
              IF(icol1.EQ.100)icol1=0
              icol2=icol_un(i,2)+100
              IF(icol2.EQ.100)icol2=0
           ELSE
              icol1=icol_un(i,2)+100
              IF(icol1.EQ.100)icol1=0
              icol2=icol_un(i,1)+100
              IF(icol2.EQ.100)icol2=0
           ENDIF
           ! from parton to hadrons                                                                                                                                                                        
           px=Phegas_hmom2(i,1)
           py=Phegas_hmom2(i,2)
           pz=Phegas_hmom2(i,3)
           p0=Phegas_hmom2(i,4)
           pmass=Phegas_hmom2(i,5)
           vtime=0
           vspin=9
           
           WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2 &
                ,px,py,pz,p0,pmass,vtime,vspin
        ENDDO
     ENDIF
  ENDIF
  Nevents=nwri
END SUBROUTINE unwei_procedure

SUBROUTINE BOOSTL2(Q,PBOO,P)
! momentums are in normal representation with fourth comp is the zero comp
! Boost P via PBOO(PBOO^2=Q^2) to PLB,and set to P
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::Q
REAL(KIND=DBL),DIMENSION(4),INTENT(IN)::PBOO
REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::P
REAL(KIND=DBL),DIMENSION(4)::PCM,PLB
REAL(KIND=DBL)::FACT
INTEGER::J
PCM(1:4)=P(1:4)
PLB(4)=(PBOO(4)*PCM(4)+PBOO(3)*PCM(3)+PBOO(2)*PCM(2)+PBOO(1)*PCM(1))/Q
FACT=(PLB(4)+PCM(4))/(Q+PBOO(4))
DO J=1,3
   PLB(J)=PCM(J)+FACT*PBOO(J)
ENDDO
P(1:4)=PLB(1:4)
END SUBROUTINE BOOSTL2

SUBROUTINE cmstolab2()
! from the cms to the lab frame 
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(20,4)::p
REAL(KIND=DBL),DIMENSION(4)::pboo
INTEGER::i,icut
REAL(KIND=DBL)::q,e,exp1,exp2
DO i=1,n
   p(i,1:4)=phegas_pmom(i,1:4)
ENDDO
icut=0
IF(istruc.EQ.1.OR..NOT.labeqcoll)THEN
  q=ehat
  e=q/SQRT(xp1*xp2)
  exp1=ebeam(1)*xp1
  exp2=ebeam(2)*xp2
  IF(fixtarget)THEN
     IF(.NOT.fixtargetrev)THEN
        exp1=exp1+FT_M*xp1/2d0
        exp2=exp2/2d0
     ELSE
        exp1=exp1/2d0
        exp2=exp2+FT_M*xp2/2d0
     ENDIF
  ENDIF
  pboo(4)=(exp1+exp2)   ! the zero component 
  pboo(3)=(exp1-exp2)   ! the z component 
  pboo(1:2)=0
  ! boost p(i,1:4) via pboo to the lab frame                                                                                    
  DO i=1,n
     CALL BOOSTL2(q,pboo,p(i,1:4))
  ENDDO
ENDIF
phegas_pmom(1:n,1:4)=p(1:n,1:4)
END SUBROUTINE cmstolab2

FUNCTION Phegas_hmom2(i,j)
USE Helac_Global
INTEGER,INTENT(IN)::i,j
INTEGER::kkkij
REAL(KIND=DBL)::Phegas_hmom2
IF(nhad.EQ.n)THEN
   Phegas_hmom2=Phegas_pmom(MOD(i-1,nhad)+1,MOD(j-1,5)+1)
   RETURN
ENDIF
kkkij=hadron2parton(MOD(i-1,nhad)+1)
IF(Quarkonium3(kkkij).EQ.0)THEN
   Phegas_hmom2=Phegas_pmom(kkkij,MOD(j-1,5)+1)
ELSE
   Phegas_hmom2=Phegas_pmom(kkkij,MOD(j-1,5)+1)+Phegas_pmom(kkkij+1,MOD(j-1,5)+1)
ENDIF
RETURN
END FUNCTION Phegas_hmom2
END MODULE Colliders_PSI_1
