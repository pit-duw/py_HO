MODULE pp_psipsi_DPS
  USE Helac_Global
  USE pp_psipsi_dps_global
  USE Helac_Func_1
  USE Kinetic_Func
  USE MC_VEGAS
  USE Constants
  USE Decay_interface
  USE pp_psipsi_DPS_ME
  USE plot_pp_psipsi_dps
  USE DecayInfo
  USE Func_PSI
  IMPLICIT NONE
  INTEGER::nprint
  INTEGER::varnum,nunwei
  LOGICAL::lunwei2
  REAL(KIND(1d0))::EBMUP1,EBMUP2
  INTEGER::NPRUP=0,lwmax
  INTEGER,PARAMETER::maxprint=8
  SAVE
CONTAINS
  SUBROUTINE calc_pp_psipsi_DPS
    IMPLICIT NONE
    CHARACTER(len=10),DIMENSION(20)::char
    INTEGER::i,i1,j1,k1,s1x
    INTEGER::nh,innum,outnum
    REAL(KIND(1d0))::w0,ptp,pqp,wme
    INTEGER::ioerror,icase=0
    LOGICAL::lunwei=.FALSE.,lhewgt=.FALSE.
    REAL(KIND(1d0)),DIMENSION(3)::rslt
    INTEGER::itmxn,ncalln
    REAL(KIND(1d0)),DIMENSION(4)::pmomtemp
    LOGICAL::lexist
    SAVE lunwei,lhewgt,icase
    WRITE(*,*)'                     THE BEGINNING OF HELAC-Onia'
    WRITE(*,*)'    AddOn Process: '
    WRITE(*,*)'    Double Patron Scattering for p p > psi(Upsilon) psi(Upsilon) + X '
    WRITE(*,*)'======================================================================='
    WRITE(*,*)'======================================================================='
    WRITE(*,*)' '
    CALL ReadElem_integer('colpar',Coll_Type)
    CALL ReadElem_logic('unwgt',lunwei)
    CALL ReadElem_logic('lhewgtup',lhewgt)
    CALL ReadElem_logic('topdrawer_output',topdrawer_output)
    CALL ReadElem_logic('gnuplot_output',gnuplot_output)
    CALL ReadElem_logic('root_output',root_output)
    plot_output=topdrawer_output.OR.gnuplot_output.OR.root_output
    CALL ReadElem_integer('literature',literature_cutoffs)
    CALL ReadElem_logic('fixtarget',fixtarget)
    IF(fixtarget.AND.Coll_Type.EQ.3)THEN
       WRITE(*,*)"ERROR: Cannot treat fixed-target experiment in e+e- collisions"
       STOP
    ENDIF
    IF(lunwei.AND.lhewgt)icase=1
    CALL ReadElem_logic("useMCFMrun",useMCFMrun)
    CALL ReadElem_logic('lhapdf',uselhapdf)
    INQUIRE(FILE=TRIM(input_dir)//"paths/lhapdfpath",EXIST=lexist)
    IF(lexist)THEN
       OPEN(UNIT=30333,FILE=TRIM(input_dir)//"paths/lhapdfpath")
       READ(30333,'(A)')LHAPath
       CLOSE(UNIT=30333)
       i=LEN_TRIM(LHAPath)
       LHAPath=LHAPath(1:i-13)
    ENDIF
    uselhapdf=uselhapdf.AND.lexist
    WRITE(*,*)'Use LHAPDF ?',uselhapdf
    CALL ReadElem_logic("absoluterap",absrap)
    ! some necesarry initialization
    OPEN(UNIT=12321,FILE="./input/states.inp")
    READ(12321,*)dps_istate(1),dps_istate(2)
    CLOSE(UNIT=12321)
    IF((dps_istate(1).LT.1.OR.dps_istate(1).GT.5).OR.&
         (dps_istate(2).LT.1.OR.dps_istate(2).GT.5))THEN
       WRITE(*,*)"ERROR:Unknown the states = ",dps_istate(1),dps_istate(2)
       WRITE(*,*)"INFO:Please set 1(J/psi) or 2(psi(2S)) or 3-4 (Y(1S-3S) in states.inp"
       STOP
    ENDIF
    DO i=1,2
       IF(dps_istate(i).EQ.1)THEN
          mpsi(i)=mjpsi
          OPEN(UNIT=12321,FILE="./input/crystalball_jpsi.inp")
          OPEN(UNIT=12322,FILE="./input/polarization_jpsi.inp")
       ELSEIF(dps_istate(i).EQ.2)THEN
          mpsi(i)=mpsi2s
          OPEN(UNIT=12321,FILE="./input/crystalball_psi2s.inp")
          OPEN(UNIT=12322,FILE="./input/polarization_psi2s.inp")
       ELSEIF(dps_istate(i).EQ.3)THEN
          mpsi(i)=mY1S
          OPEN(UNIT=12321,FILE="./input/crystalball_Y1S.inp")
          OPEN(UNIT=12322,FILE="./input/polarization_Y1S.inp")
       ELSEIF(dps_istate(i).EQ.4)THEN
          mpsi(i)=mY2S
          OPEN(UNIT=12321,FILE="./input/crystalball_Y2S.inp")
          OPEN(UNIT=12322,FILE="./input/polarization_Y2S.inp")
       ELSE
          mpsi(i)=mY3S
          OPEN(UNIT=12321,FILE="./input/crystalball_Y3S.inp")
          OPEN(UNIT=12322,FILE="./input/polarization_Y3S.inp")
       ENDIF
       READ(12321,*)dps_n(i),dps_ptavg(i),dps_kapa(i),dps_lam(i)
       READ(12322,*)dps_lambdath(i)
       CLOSE(UNIT=12321)
       CLOSE(UNIT=12322)
    ENDDO
    OPEN(UNIT=12321,FILE="./input/includeqq.inp")
    READ(12321,*)dps_includeqq
    CLOSE(UNIT=12321)
    OPEN(UNIT=12321,FILE="./input/sigma_eff.inp")
    READ(12321,*)dps_sigmaeff
    CLOSE(UNIT=12321)
    ! convert mb to nb
    dps_sigmaeff=dps_sigmaeff*1d6
    iflh(1)=35
    iflh(2)=35
    iflh(3)=35
    iflh(4)=35
    IF(dps_istate(1).LE.2)THEN
       iflh(5)=443011
    ELSE
       iflh(5)=553011
    ENDIF
    iflh(6)=35
    IF(dps_istate(2).LE.2)THEN
       iflh(7)=443011
    ELSE
       iflh(7)=553011
    ENDIF
    iflh(8)=35
    nhad=8
    CALL ReadDecayInfo
    IF(nunit1.NE.6)THEN
       OPEN(UNIT=nunit1,FILE=TRIM(output_dir)//"RESULT_pp_psipsi_dps.out")
    ENDIF
    nunit2=32
    CLOSE(nunit2)
    OPEN(nunit2,FILE=TRIM(output_dir)//'kine_pp_psipsi_dps.out')
    nunit3=30
    CLOSE(nunit3)
    OPEN(nunit3,FILE=TRIM(tmp_dir)//'even_pp_psipsi_dps.out',FORM='unformatted')
    ! for LHA                                                                                                                              
    CLOSE(200)
    OPEN(200,FILE=TRIM(tmp_dir)//'sample_pp_psipsi_dps.init')
    CALL ReadElem_integer("itmax",itmxn)
    CALL ReadElem_integer("nmc",ncalln)
    WRITE(*,*)' '
    CALL Helac_mtime()
    WRITE(*,*)' '
    CALL pp_psipsi_DPS_VEGAS(rslt,itmxn,ncalln)
    WRITE(nunit1,*)"sigma(nb)                   sd"
    WRITE(nunit1,*)rslt(1),rslt(2)
    WRITE(*,*)' '
    CALL Helac_mtime()
    CLOSE(nunit2)
    CLOSE(nunit3)
    CLOSE(21)
    CLOSE(200)
    CALL Generate_lhe_pp_psipsi_dps(nhad,Nevents,icase)
  END SUBROUTINE calc_pp_psipsi_DPS

  SUBROUTINE pp_psipsi_DPS_VEGAS(rslt,itmxn,ncalln)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(3),INTENT(OUT)::rslt
    INTEGER,INTENT(IN),OPTIONAL::itmxn,ncalln
    REAL(KIND(1d0))::vfes,sd,chi2a
    INTEGER::ncalm,nc,ii
    CHARACTER(len=4),DIMENSION(20)::chchar
    INTEGER::iday0,ihr0,imin0,isec0,i100th0,iday1,ihr1,imin1,isec1,i100th1,iyr0,iyr1,imon0,imon1
    INTEGER::IDBMUP1,IDBMUP2,IDWTUP
    CALL ReadElem_logic('unwgt',lunwei2)
    CALL ReadElem_integer('preunw',nunwei)
    lwmax=0
    NPRN=-1
    varnum=3
    IF(lunwei2)varnum=varnum+1
    varnum=2*varnum+NDecayIflh
    DO ii=1,varnum
       XL(ii)=0.0d0
       XU(ii)=1.0d0
    ENDDO
    IF(PRESENT(itmxn))THEN
       ITMX=itmxn
    ELSE
       ITMX=5
    ENDIF
    ITMX=1
    IF(PRESENT(ncalln))THEN
       ncalm=ncalln
    ELSE
       ncalm=5000
    ENDIF
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
    nprint=10000
    NCALL=20000
    CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
    CALL VEGAS(varnum,pp_psipsi_DPS_fxn,vfes,sd,chi2a)
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
       CALL VEGAS(varnum,pp_psipsi_DPS_fxn,vfes,sd,chi2a,1)
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
         CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
         WRITE(*,*)"====================NCALL="//chchar(ii)//"==========================="
         WRITE(*,*)" "
         CALL Helac_mtime()
         WRITE(*,*)" "
         IF(plot_output)CALL initplot_pp_psipsi_DPS
         CALL VEGAS(varnum,pp_psipsi_DPS_fxn,vfes,sd,chi2a,1)
         IF(plot_output)CALL plotout_pp_psipsi_DPS
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
         IF(plot_output)CALL initplot_pp_psipsi_DPS
         CALL VEGAS(varnum,pp_psipsi_DPS_fxn,vfes,sd,chi2a,1)
         IF(plot_output)CALL plotout_pp_psipsi_DPS
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
      WRITE(200,5200) vfes*10d0**3,sd*10d0**3,1d0, 82
      rslt(1)=vfes
      rslt(2)=sd
      rslt(3)=chi2a
      IF(lunwei2)PRINT *,"number of events",Nevents
      RETURN
5100  FORMAT(1P,2I8,2E14.6,6I8)
5200  FORMAT(1P,3E20.10,I6)
    END SUBROUTINE pp_psipsi_DPS_VEGAS

    FUNCTION pp_psipsi_DPS_fxn(x,wgt)
      REAL(KIND(1d0)),DIMENSION(varnum),INTENT(IN)::x
      REAL(KIND(1d0)),INTENT(IN)::wgt
      REAL(KIND(1d0))::pp_psipsi_DPS_fxn
      REAL(KIND(1d0))::y1,y2,phi,w1
      REAL(KIND(1d0)),DIMENSION(4)::pmom
      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
      INTEGER::init=0,nwarn,nwmax,nwri,nwri_tot,icut,nnn=0,nnntot=0,&
           pdfnumpdf,ivarnum
      SAVE init,nwarn,nwmax,nwri,nwri_tot,nnn,nnntot,pdfnumpdf,ivarnum
      REAL(KIND(1d0))::sqs,sq,pt1c,maxpt1c,y1cup,y1clow,ycollcm,mt1,mt2,pt
      SAVE sqs,sq,pt1c,maxpt1c,y1cup,y1clow,ycollcm
      REAL(KIND(1d0))::wdps,wsf,wme,wgtbr,temp1,temp2,temp21,temp3,temp4
      REAL(KIND(1d0))::exp1,exp2,Jpt2,Jy1,Jy2,xp11,xp21,scale1,scale2
      REAL(KIND(1d0))::recmax=0,recmin=0
      SAVE temp1
      INTEGER::ipip,ioffset,j
      !INCLUDE "pp_psipsi_dps.inc"
      IF(init.EQ.0)THEN
         wjac=1
         nwmax=0
         nwarn=0
         nwri=0
         CALL ReadElem_integer('unwevt',nwri_tot)
         CALL ReadElem_real('energy_beam1',sqs)
         ebeam(1)=ABS(sqs)
         CALL ReadElem_real('energy_beam2',sqs)
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
         CALL ReadElem_real('minpt1c',pt1c)
         CALL ReadElem_real('maxpt1c',maxpt1c)
         IF(maxpt1c.LT.0d0)THEN
            maxpt1c=-1d0
         ELSE
            IF(maxpt1c.LE.pt1c)THEN
               WRITE(*,*)"ERROR: the first final state was cut off by pt cut (pt1c,maxpt1c)"
               STOP
            ENDIF
         ENDIF
         CALL ReadElem_real('maxy1c',y1cup)
         CALL ReadElem_real('miny1c',y1clow)
         CALL ReadElem_integer('pdf',pdfnumpdf)
         ptc(3:20)=0
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
         gmas(3:20,3:20)=0
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
         ENDIF
         xp1=1
         xp2=1
         IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
            CALL readcuts_pp_psipsi_dps
            y1cup=MIN(y1cup,yycut(3))
            y1clow=MAX(y1clow,yycutlow(3))
            IF(absrap)y1cup=ABS(y1cup)
            IF(absrap)y1clow=ABS(y1clow)
         ENDIF
         IF(NDecayChains.GT.0)CALL readcuts_Decay
         sq=sqs*sqs
         !temp1=(sq-mpsi**2)/(2d0*sqs)
         ivarnum=3
         IF(lunwei2)ivarnum=ivarnum+1
         init=1
      ENDIF
      nnntot=nnntot+1
      wdps=1d0
      DO ipip=1,2
         ! the ipip-th SPS
         temp1=(sq-mpsi(ipip)**2)/(2d0*sqs) 
         IF(maxpt1c.GE.0d0)THEN
            temp1=MIN(temp1,maxpt1c)
         ENDIF
         ioffset=ivarnum*(ipip-1)
         pt=(temp1-pt1c)*x(3+ioffset)+pt1c
         mt1=DSQRT(pt**2+mpsi(ipip)**2)
         mt2=pt
         IF(absrap)THEN
            temp2=MIN(ACosh_p((sq+mpsi(ipip)**2)/(2d0*sqs*mt1)),y1cup)
            IF(temp2.LT.y1clow)THEN
               pp_psipsi_DPS_fxn=0d0
               RETURN
            ENDIF
            y1=(2d0*x(1+ioffset)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1+ioffset)-0.5d0)*y1clow
         ELSE
            temp2=MIN(ACosh_p((sq+mpsi(ipip)**2)/(2d0*sqs*mt1)),y1cup-ycollcm)
            temp21=MAX(-ACosh_p((sq+mpsi(ipip)**2)/(2d0*sqs*mt1)),y1clow-ycollcm)
            IF(temp2.LT.temp21)THEN
               pp_psipsi_DPS_fxn=0d0
               RETURN
            ENDIF
            y1=(temp2-temp21)*x(1+ioffset)+temp21 ! in collision frame
         ENDIF
         temp3=DLOG((-DEXP(-y1)*mt1+sqs)/mt2)
         temp4=temp3+DLOG((-DEXP(y1)*mt1+sqs)/mt2)
         y2=-temp3+temp4*x(2+ioffset) ! in collision frame
         Jpt2=2d0*pt*(temp1-pt1c)
         IF(absrap)THEN
            Jy1=2d0*(temp2-y1clow)
         ELSE
            Jy1=temp2-temp21
         ENDIF
         Jy2=temp4
         ! The substitution of original variations
         IF(ipip.EQ.2)THEN
            IF((DEXP(y1)*mt1+DEXP(y2)*mt2)/sqs+xp1.GE.1d0)THEN
               pp_psipsi_DPS_fxn=0d0
               RETURN
            ENDIF
            IF((DEXP(-y1)*mt1+DEXP(-y2)*mt2)/sqs+xp2.GE.1d0)THEN
               pp_psipsi_DPS_fxn=0d0
               RETURN
            ENDIF
         ENDIF
         xp1=(DEXP(y1)*mt1+DEXP(y2)*mt2)/sqs
         xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2)/sqs
         ehat=DSQRT(xp1*xp2*sq)
         wdps=wdps*xp1*xp2*Jpt2*Jy1*Jy2*3.8938573d5/(16d0*pi*ehat**4)
         exp1=xp1*ebeam(1)
         exp2=xp2*ebeam(2)
         IF(fixtarget)THEN
            IF(.NOT.fixtargetrev)THEN
               exp1=exp1+FT_M*xp1/2d0
               exp2=exp2/2d0
            ELSE
               exp2=exp2+FT_M*xp2/2d0
               exp1=exp1/2d0
            ENDIF
         ENDIF
         ! Generate the momenta of external legs
         dps_pmom(ipip,1,1:2)=0
         dps_pmom(ipip,1,3)=xp1*sqs/2d0
         dps_pmom(ipip,1,4)=xp1*sqs/2d0
         dps_pmom(ipip,2,1:2)=0
         dps_pmom(ipip,2,3)=-xp2*sqs/2d0
         dps_pmom(ipip,2,4)=xp2*sqs/2d0
         ! we choose phi=0
         IF(lunwei2)THEN
            phi=2*pi*x(4+ioffset)
         ELSE
            phi=0d0
         ENDIF
         dps_pmom(ipip,3,1)=pt*DCOS(phi)
         dps_pmom(ipip,3,2)=pt*DSIN(phi)
         dps_pmom(ipip,3,3)=mt1*DSINH(y1)
         dps_pmom(ipip,3,4)=mt1*DCOSH(y1)
         dps_pmom(ipip,4,1)=-pt*DCOS(phi)
         dps_pmom(ipip,4,2)=-pt*DSIN(phi)
         dps_pmom(ipip,4,3)=mt2*DSINH(y2)
         dps_pmom(ipip,4,4)=mt2*DCOSH(y2)
         ! boost from collision frame to lab frame
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
         IF(.NOT.labeqcoll)THEN
            DO j=1,4
               CALL Boostl(sqs,pmom,dps_pmom(ipip,j,1:4))
            ENDDO
         ENDIF
         IF(pdfnumpdf.EQ.921000)THEN
            ! GS 09 dPDF
            IF(ipip.EQ.1)THEN
               xp11=xp1
               xp21=xp2
               scale1=dps_pmom(ipip,3,1)**2+dps_pmom(ipip,3,2)**2+mpsi(ipip)**2
               scale1=DSQRT(scale1)
            ELSE
               scale2=dps_pmom(ipip,3,1)**2+dps_pmom(ipip,3,2)**2+mpsi(ipip)**2
               scale2=DSQRT(scale2)
               CALL GS09(xp11,xp1,DSQRT(scale1*scale2),0,0,wsf)
               wdps=wdps*wsf
               CALL GS09(xp21,xp2,DSQRT(scale1*scale2),0,0,wsf)
               wdps=wdps*wsf
               ! CALL GSALPS(Q) to run alpha_S, which is not used here
            ENDIF
         ELSE
            CALL strf_pdf_pp_psipsi_dps(ipip,wsf)
            wdps=wdps*wsf
         ENDIF
      ENDDO
      IF(wdps.LE.0d0)THEN
         pp_psipsi_DPS_fxn=0d0
         RETURN
      ENDIF
      CALL pp_psipsi_dps_hadronmom
      icut=1
      CALL Cuts_pp_psipsi_dps(icut)
      IF(icut.EQ.0)THEN
         pp_psipsi_DPS_fxn=0d0
         RETURN
      ENDIF
      DO j=2*ivarnum+1,varnum
         Decayran(j-2*ivarnum)=x(j)
      ENDDO
      IF(NDecayChains.GT.0)THEN
         CALL HO_Decay_pp_psipsi_dps(weight_br)
         IF(weight_br.LE.0d0)THEN
            pp_psipsi_DPS_fxn=0d0
            RETURN
         ENDIF
         CALL cuts_Decay_pp_psipsi_dps(icut)
         IF(icut.EQ.0)THEN
            pp_psipsi_DPS_fxn=0d0
            RETURN
         ENDIF
         wdps=wdps*weight_br
      ENDIF
      DO ipip=1,2
         CALL crystalball_gg_psiX(ipip,wme)
         wdps=wdps*wme
         IF(wdps.LE.0d0)THEN
            pp_psipsi_DPS_fxn=0d0
            RETURN
         ENDIF
      ENDDO
      IF(dps_istate(1).EQ.dps_istate(2))THEN
         pp_psipsi_DPS_fxn=wdps/2d0/dps_sigmaeff
      ELSE
         pp_psipsi_DPS_fxn=wdps/dps_sigmaeff
      ENDIF
      IF(lunwei2.AND.lwmax.GT.0)THEN
         w1=pp_psipsi_DPS_fxn*wgt ! multiply VEGAS weight
         CALL unwei_procedure_pp_psipsi_dps(w1,nwri,nwmax,nwarn)
      ENDIF
      IF(plot_output)THEN
         w1=pp_psipsi_DPS_fxn*wgt
         CALL outfun_pp_psipsi_dps(w1)
      ENDIF
      nnn=nnn+1
      IF(recmax.LT.pp_psipsi_DPS_fxn)recmax=pp_psipsi_DPS_fxn
      IF(recmin.GT.pp_psipsi_DPS_fxn)recmin=pp_psipsi_DPS_fxn
      IF(MOD(nnn,nprint).EQ.0)THEN
         PRINT *,"max=",recmax
         PRINT *,"min=",recmin
         IF(lunwei2.AND.lwmax.GT.1)THEN
            PRINT *,"      n_event,    n_pass,    n_total"
            PRINT *,nwri,nnn,nnntot
         ELSE
            PRINT *, "      n_pass,     n_total"
            PRINT *,nnn,nnntot
         ENDIF
      ENDIF
    END FUNCTION pp_psipsi_DPS_fxn

    SUBROUTINE strf_pdf_pp_psipsi_dps(ipip,wsf)
      USE CTEQ6PDF
      USE Structf_PDFs
      IMPLICIT NONE
      INTEGER::ipp=1
      INTEGER::ih=1 ! ih=1 no photon PDF, ih=2, photon from proton/anti-proton, ih=3 phton from electron/positron
      REAL(KIND(1d0)),INTENT(OUT)::wsf
      INTEGER,INTENT(IN)::ipip
      INTEGER::init=0
      LOGICAL::use_cteq6_f90=.TRUE.
      REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
      SAVE init,ipp,use_cteq6_f90
      REAL(KIND(1d0))::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
                ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
      INCLUDE "../lhapdf/call_strf_lhapdf"
      IF(init.EQ.0)THEN
         SELECT CASE(iPDFSUP1)
         CASE(10000)
            CALL SetCtq6f90(1)
            pdlabel='cteq6_m'
            nloop=2
            alphaQCD2=0.118d0
            use_cteq6_f90=.TRUE.
         CASE(10041)
            CALL SetCtq6f90(3)
            pdlabel='cteq6_l'
            nloop=2
            alphaQCD2=0.118d0
            use_cteq6_f90=.TRUE.
         CASE(10042)
            CALL SetCtq6f90(4)
            pdlabel='cteq6l1'
            nloop=1
            alphaQCD2=0.130d0
            use_cteq6_f90=.TRUE.
         CASE DEFAULT
            CALL pdfset_internal
            use_cteq6_f90=.FALSE.
         END SELECT
         IF(COLL_TYPE.EQ.1)THEN
            ipp=2
         ELSE
            ipp=1
         ENDIF
         init=1
      ENDIF

      scale=dps_pmom(ipip,3,1)**2+dps_pmom(ipip,3,2)**2+mpsi(ipip)**2
      scale=DSQRT(scale)

      IF(use_cteq6_f90)THEN
         glu1_ct = xp1*Ctq6Pdf_f90(0,xp1,scale)
         glu2_ct = xp2*Ctq6Pdf_f90(0,xp2,scale)
         u1_ct   = xp1*Ctq6Pdf_f90(1,xp1,scale)
         u2_ct   = xp2*Ctq6Pdf_f90(1,xp2,scale)
         d1_ct   = xp1*Ctq6Pdf_f90(2,xp1,scale)
         d2_ct   = xp2*Ctq6Pdf_f90(2,xp2,scale)
         s1_ct   = xp1*Ctq6Pdf_f90(3,xp1,scale)
         s2_ct   = xp2*Ctq6Pdf_f90(3,xp2,scale)
         c1_ct   = xp1*Ctq6Pdf_f90(4,xp1,scale)
         c2_ct   = xp2*Ctq6Pdf_f90(4,xp2,scale)
         b1_ct   = xp1*Ctq6Pdf_f90(5,xp1,scale)
         b2_ct   = xp2*Ctq6Pdf_f90(5,xp2,scale)
         ub1_ct  = xp1*Ctq6Pdf_f90(-1,xp1,scale)
         ub2_ct  = xp2*Ctq6Pdf_f90(-1,xp2,scale)
         db1_ct  = xp1*Ctq6Pdf_f90(-2,xp1,scale)
         db2_ct  = xp2*Ctq6Pdf_f90(-2,xp2,scale)
         sb1_ct  = xp1*Ctq6Pdf_f90(-3,xp1,scale)
         sb2_ct  = xp2*Ctq6Pdf_f90(-3,xp2,scale)
         cb1_ct  = xp1*Ctq6Pdf_f90(-4,xp1,scale)
         cb2_ct  = xp2*Ctq6Pdf_f90(-4,xp2,scale)
         bb1_ct  = xp1*Ctq6Pdf_f90(-5,xp1,scale)
         bb2_ct  = xp2*Ctq6Pdf_f90(-5,xp2,scale)
      ELSE
         CALL fdist(ih,xp1,scale,pdflist(-7:7))
         glu1_ct = xp1*pdflist(0)
         u1_ct   = xp1*pdflist(2)
         d1_ct   = xp1*pdflist(1)
         s1_ct   = xp1*pdflist(3)
         c1_ct   = xp1*pdflist(4)
         b1_ct   = xp1*pdflist(5)
         ub1_ct  = xp1*pdflist(-2)
         db1_ct  = xp1*pdflist(-1)
         sb1_ct  = xp1*pdflist(-3)
         cb1_ct  = xp1*pdflist(-4)
         bb1_ct  = xp1*pdflist(-5)
         CALL fdist(ih,xp2,scale,pdflist(-7:7))
         glu2_ct = xp2*pdflist(0)
         u2_ct   = xp2*pdflist(2)
         d2_ct   = xp2*pdflist(1)
         s2_ct   = xp2*pdflist(3)
         c2_ct   = xp2*pdflist(4)
         b2_ct   = xp2*pdflist(5)
         ub2_ct  = xp2*pdflist(-2)
         db2_ct  = xp2*pdflist(-1)
         sb2_ct  = xp2*pdflist(-3)
         cb2_ct  = xp2*pdflist(-4)
         bb2_ct  = xp2*pdflist(-5)
      ENDIF
      wsf=glu1_ct*glu2_ct*wjac/xp1/xp2
      IF(dps_includeqq)THEN
         wsf=wsf+(u1_ct*ub2_ct+d1_ct*db2_ct+s1_ct*sb2_ct+c1_ct*cb2_ct+&
              ub1_ct*u2_ct+db1_ct*d2_ct+sb1_ct*s2_ct+cb1_ct*c2_ct)*wjac/xp1/xp2
      ENDIF
      IF(init.EQ.0)init=1
    END SUBROUTINE strf_pdf_pp_psipsi_dps

    SUBROUTINE readcuts_pp_psipsi_dps
      IMPLICIT NONE
      CHARACTER(len=24)::file
      LOGICAL::lexist
      INTEGER::iounit,flag=0,i,i1,j,j1
      REAL(KIND(1d0))::ptq,ptcharm,ptconia,ptbonia,etaq,ycq,ycqlow,etaconia,ycconia,ycconialow,ycbonia,ycbonialow
      REAL(KIND(1d0))::cutoff,xFcq,xFcqlow,xFcconia,xFcconialow,xFcbonia,xFcbonialow,etabonia
      REAL(KIND(1d0))::gbeamq,drqq,gqq,maxptq,maxptcharm,maxptconia,maxptbonia
      ! open default input file
      INQUIRE(FILE=TRIM(input_dir)//"default.inp",EXIST=lexist)
      IF(.NOT.lexist)THEN
         PRINT *,"Warning: the file default.inp does not exist ! STOP !"
         STOP
      ENDIF
      INQUIRE(FILE=TRIM(input_dir)//"default.inp",OPENED=lexist)
      IF(lexist)THEN
         INQUIRE(FILE=TRIM(input_dir)//"default.inp",NUMBER=iounit)
         IF(iounit.NE.udefault)THEN
            PRINT *,"WARNING: the default.inp has been linked with another unit ! Close and reopen !"
            CLOSE(UNIT=iounit)
            OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
         ENDIF
      ELSE
         OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
      ENDIF
      ! open user's input file
      IF(TRIM(Input_File)/="default.inp")THEN
         INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
         IF(.NOT.lexist)THEN
            PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
            STOP
         ENDIF
         INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
         IF(lexist)THEN
            INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
            IF(iounit.NE.uinput)THEN
               PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
               CLOSE(UNIT=iounit)
               OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
            ENDIF
         ELSE
            OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
         ENDIF
      ELSE
         flag=1
      ENDIF
      cutoff=readvalue_r("cutoffp",flag)
      PRINT *,"WARNING CUTOFF SET:",cutoff
      ptc(3:20)=cutoff
      maxptc(3:20)=-1d0
      drc(3:20,3:20)=0
      etac(3:20)=20
      yycut(3:20)=1d9
      IF(absrap)THEN
         yycutlow(3:20)=0d0
      ELSE
         yycutlow(3:20)=-1d9
      ENDIF
      xFcut(3:20)=1d0
      xFcutlow(3:20)=1d0
      xFcutflag=.FALSE.
      ec(3:20)=cutoff
      c1(3:20)=1d0
      c2(3:20)=1d0
      cc(3:20,3:20)=1d0
      gmas(3:20,3:20)=cutoff
      gbeammass(1:2,3:20)=0d0
      y1cup=30d0
      y1clow=0d0
      y1cup=readvalue_r("maxy1c",flag)
      y1clow=readvalue_r("miny1c",flag)
      ! minimum quark pt
      ptq=readvalue_r("minptq",flag)
      ! maximum quark pt
      maxptq=readvalue_r("maxptq",flag)
      ! minimum charm pt
      ptcharm=readvalue_r("minptc",flag)
      ! maximum charm pt
      maxptcharm=readvalue_r("maxptc",flag)
      ! minimum charmonia pt
      ptconia=readvalue_r("minptconia",flag)
      ! maximum charmonia pt
      maxptconia=readvalue_r("maxptconia",flag)
      ! minimum bottonia pt
      ptbonia=readvalue_r("minptbonia",flag)
      ! maximum bottonia pt
      maxptbonia=readvalue_r("maxptbonia",flag)
      ! maximum rapidity quark 
      etaq=readvalue_r("maxrapq",flag)
      ! maximum y rapidity quark
      ycq=readvalue_r("maxyrapq",flag)
      ! minimum y rapidity quark
      ycqlow=readvalue_r("minyrapq",flag)
      ! maximum Feynman parameter xF
      xFcq=readvalue_r("maxxFq",flag)
      ! minimum Feynman parameter xF
      xFcqlow=readvalue_r("minxFq",flag)
      ! maximum rapidity charmonium
      etaconia=readvalue_r("maxrapconia",flag)
      ! maximum y rapidity charmonium
      ycconia=readvalue_r("maxyrapconia",flag)
      ! minimum y rapidity charmonium
      ycconialow=readvalue_r("minyrapconia",flag)
      ! maximum Feynman parameter xF
      xFcconia=readvalue_r("maxxFconia",flag)
      ! minimum Feynman parameter xF
      xFcconialow=readvalue_r("minxFconia",flag)
      ! maximum rapidity bottomonia
      etabonia=readvalue_r("maxrapbonia",flag)
      ! maximum y rapidity bottomonia
      ycbonia=readvalue_r("maxyrapbonia",flag)
      ! minimum y rapidity bottomonia
      ycbonialow=readvalue_r("minyrapbonia",flag)
      ! maximum Feynman parameter xF
      xFcbonia=readvalue_r("maxxFbonia",flag)
      ! minimum Feynman parameter xF
      xFcbonialow=readvalue_r("minxFbonia",flag)
      ! minimum DR quark with quark
      drqq=readvalue_r("mindrqq",flag)
      ! minimum mass quark with quark
      gqq=readvalue_r("minmqqp",flag)
      ! minimum mass u,d,s quarks and gluon  with partonic beam
      gbeamq=readvalue_r("minmqbeam",flag)
      CLOSE(UNIT=udefault)
      CLOSE(UNIT=uinput)
      DO i=5,nhad
         i1=iflh(i)
         IF(iflh(i).EQ.35)i1=0
         IF(i1.EQ.0)THEN
            ptc(i)=ptq
            IF(maxptq.GE.0d0)maxptc(i)=maxptq
            etac(i)=etaq
            yycut(i)=ycq
            yycutlow(i)=ycqlow
            xFcut(i)=xFcq
            xFcutlow(i)=xFcqlow
         ELSEIF(i1.GE.440000.AND.i1.LE.449999)THEN
            ptc(i)=ptconia
            IF(maxptconia.GE.0d0)maxptc(i)=maxptconia
            etac(i)=etaconia
            yycut(i)=ycconia
            yycutlow(i)=ycconialow
            xFcut(i)=xFcconia
            xFcutlow(i)=xFcconialow
         ELSE
            ptc(i)=ptbonia
            IF(maxptbonia.GE.0d0)maxptc(i)=maxptbonia
            etac(i)=etabonia
            yycut(i)=ycbonia
            yycutlow(i)=ycbonialow
            xFcut(i)=xFcbonia
            xFcutlow(i)=xFcbonialow
         ENDIF
      ENDDO
      DO i=5,nhad-1
         i1=iflh(i)
         IF(i1.EQ.35)i1=0
         DO j=i+1,nhad
            j1=iflh(j)
            IF(j1.EQ.35)j1=0
            IF(i1.EQ.0.AND.j1.EQ.0)THEN
               drc(i,j)=drqq
               gmas(i,j)=MAX(gqq,gmas(i,j))
            ENDIF
         ENDDO
      ENDDO
      DO i=5,nhad-1
         DO j=i+1,nhad
            gmas(i,j)=MAX(gmas(i,j),DSQRT(2*ptc(i)*ptc(j)*(1-COS(drc(i,j)))))
         ENDDO
      ENDDO
      DO i=6,nhad
         DO j=5,i-1
            drc(i,j)=drc(j,i)
            gmas(i,j)=gmas(j,i)
         ENDDO
      ENDDO
      WRITE(*,*)'---------------------------------------------------'
      WRITE(*,*)'    the cuts for p p > psi(Upsilon) psi(Upsilon) + X '
      WRITE(*,*)'    with double parton scattering (DPS) '
      DO i=5,nhad
         WRITE(*,*)'pt     of  ',i,'   particle   ',ptc(i)
         IF(maxptc(i).GE.0d0)THEN
            WRITE(*,*)'max pt     of  ',i,'   particle   ',maxptc(i)
            IF(maxptc(i).LE.ptc(i))THEN
               WRITE(*,*)"ERROR: One of the final state ",i," was cut off by pt cut"
               STOP
            ENDIF
         ELSE
            WRITE(*,*)'no max pt cut of  ',i,'   particle  '
         ENDIF
         WRITE(*,*)'energy of  ',i,'   particle   ',ec(i)
         WRITE(*,*)'rapidity of  ',i,'   particle   ',etac(i)
         WRITE(*,*)'max y rapidity of ',i,'   particle   ',yycut(i)
         WRITE(*,*)'min y rapidity of ',i,'   particle   ',yycutlow(i)
         IF(yycut(i).LE.yycutlow(i))THEN
            WRITE(*,*)"ERROR: One of the final state ",i," was cut off by y rapidity cut"
            STOP
         ENDIF
         WRITE(*,*)'max Feynman parameter xF of ',i,' particle ',xFcut(i)
         IF(xFcut(i).LT.1d0)xFcutflag=.TRUE.
         WRITE(*,*)'min Feynman parameter xF of ',i,' particle ',xFcutlow(i)
         IF(xFcutlow(i).GT.-1d0)xFcutflag=.TRUE.
      ENDDO
      WRITE(*,*)'The maxrapidity of the first particle',y1cup
      WRITE(*,*)'The minrapidity of the first particle',y1clow
      DO i=5,nhad-1
         DO j=i+1,nhad
            WRITE(*,*)'DR     ',i,'  with  ',j,drc(i,j)
            WRITE(*,*)'mass of ',i,'  with  ',j,gmas(i,j)
         ENDDO
      ENDDO
      WRITE(*,*)'---------------------------------------------------'
     END SUBROUTINE readcuts_pp_psipsi_dps

    SUBROUTINE Cuts_pp_psipsi_dps(icut)
      IMPLICIT NONE
      INTEGER,INTENT(OUT)::icut
      INTEGER::l,l1,l2,flag
      REAL(KIND(1d0))::s,d1,d2,dr,pt,eta,aaa,bbb,ptcut
      REAL(KIND(1d0)),DIMENSION(4)::ponia,pboo2
      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
      REAL(KIND(1d0))::e,q
      !INCLUDE "pp_psipsi_dps.inc"
      icut=0
      ! invariant mass cuts
      DO l1=5,nhad-1
         DO l2=l1+1,nhad
            s=2*scalar_product(dps_hadron_pmom(l1,1:4),dps_hadron_pmom(l2,1:4))&
                 +scalar_product(dps_hadron_pmom(l1,1:4),dps_hadron_pmom(l1,1:4))&
                 +scalar_product(dps_hadron_pmom(l2,1:4),dps_hadron_pmom(l2,1:4))
            IF(s.LT.gmas(l1,l2)**2)RETURN
         ENDDO
      ENDDO
      flag=0
      DO l=5,nhad
         pt=SQRT(dps_hadron_pmom(l,1)**2+dps_hadron_pmom(l,2)**2)
         IF(pt.LT.ptc(l))THEN
            flag=1
            EXIT
         ENDIF
         IF(maxptc(l).GE.0d0)THEN
            IF(pt.GT.maxptc(l))THEN
               flag=1
               EXIT
            ENDIF
         ENDIF
         eta=prapidity(dps_hadron_pmom(l,1:4))
         IF(ABS(eta).GT.etac(l))THEN
            flag=1
            EXIT
         ENDIF
         eta=rapidity(dps_hadron_pmom(l,1:4))
         IF(absrap)eta=ABS(eta)
         IF(eta.GT.yycut(l).OR.eta.LT.yycutlow(l))THEN
            flag=1
            EXIT
         ENDIF
         ! special for the first particle , which can be used to calculate the y distribution  
         IF(l.EQ.5.OR.l.EQ.7)THEN
            IF(eta.GT.y1cup.OR.eta.LT.y1clow)THEN
               flag=1
               EXIT
            ENDIF
         ENDIF
      ENDDO
      IF(flag.EQ.0)THEN
         DO l1=5,nhad-1
            DO l2=l1+1,nhad
               d1=prapidity(dps_hadron_pmom(l1,1:4))-prapidity(dps_hadron_pmom(l2,1:4))
               d2=ph4(dps_hadron_pmom(l1,1),dps_hadron_pmom(l1,2),dps_hadron_pmom(l1,3))&
                    -ph4(dps_hadron_pmom(l2,1),dps_hadron_pmom(l2,2),dps_hadron_pmom(l2,3))
               d2=MIN(DABS(d2),2*pi-DABS(d2))
               IF(d2/pi.GT.1d0)WRITE(*,*)d2/pi
               dr=SQRT(d1**2+d2**2)
               IF(dr.LT.drc(l1,l2))THEN
                  flag=1
                  EXIT
               ENDIF
            ENDDO
            IF(flag.EQ.1)EXIT
         ENDDO
         IF(flag.EQ.0)icut=1
      ENDIF
      ! special cutoffs in the literature
      IF(icut.EQ.1.AND.flag.EQ.0)THEN
         IF(literature_cutoffs.EQ.14060484)THEN
            DO l=5,nhad
               IF(iflh(l).GT.440001.AND.iflh(l).LT.449999)THEN
                  ! charmonia
                  ponia(1:4)=dps_hadron_pmom(l,1:4)
                  eta=rapidity(ponia(1:4))
                  IF(ABS(eta).LT.1.2d0)THEN
                     pt=SQRT(ponia(1)**2+ponia(2)**2)
                     IF(pt.LE.6.5d0)THEN
                        icut=0
                        EXIT
                     ENDIF
                  ELSEIF(ABS(eta).GT.1.2d0.AND.ABS(eta).LT.1.43d0)THEN
                     pt=SQRT(ponia(1)**2+ponia(2)**2)
                     aaa=-8.695652173913045d0
                     bbb=16.934782608695656d0
                     ptcut=aaa*ABS(eta)+bbb
                     IF(pt.LE.ptcut)THEN
                        icut=0
                        EXIT
                     ENDIF
                  ELSEIF(ABS(eta).GT.1.43d0.AND.ABS(eta).LT.2.2d0)THEN
                     pt=SQRT(ponia(1)**2+ponia(2)**2)
                     IF(pt.LE.4.5d0)THEN
                        icut=0
                        EXIT
                     ENDIF
                  ELSE
                     icut=0
                     EXIT
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF
      IF(xFcutflag.AND.icut.EQ.1.AND.flag.EQ.0)THEN
         q=ehat
         e=q/SQRT(xp1*xp2)
         IF(.NOT.labeqcoll)THEN
            IF(.NOT.fixtarget)THEN
               pboo2(4)=(ebeam(1)+ebeam(2))
               pboo2(3)=-(ebeam(1)-ebeam(2))
            ELSE
               IF(.NOT.fixtargetrev)THEN
                  pboo2(4)=(ebeam(1)+ebeam(2))
                  pboo2(3)=-ebeam(1)
               ELSE
                  pboo2(4)=(ebeam(1)+ebeam(2))
                  pboo2(3)=ebeam(2)
               ENDIF
            ENDIF
            pboo2(1:2)=0
            ! boost from the lab frame to the collision frame 
            DO l=5,nhad
               CALL boostl(e,pboo2,dps_hadron_pmom(l,1:4))
            ENDDO
         ENDIF
         DO l=5,nhad
            eta=xFeynman(dps_hadron_pmom(l,1:4),e)
            IF(eta.GT.xFcut(l).OR.eta.LT.xFcutlow(l))THEN
               icut=0
               EXIT
            ENDIF
         ENDDO
         IF(.NOT.labeqcoll)THEN
            pboo2(3)=-pboo2(3)
            ! boost back to the lab frame 
            DO l=5,nhad
               CALL boostl(e,pboo2,dps_hadron_pmom(l,1:4))
            ENDDO
         ENDIF
      ENDIF
    END SUBROUTINE Cuts_pp_psipsi_dps

    SUBROUTINE HO_Decay_pp_psipsi_dps(wgtbr)
      IMPLICIT NONE
      INCLUDE "stdhep_pp_psipsi_dps.inc"
      INTEGER::i,j,k,kk,m,ndecay
      REAL(KIND(1d0)),INTENT(OUT)::wgtbr
      REAL(KIND(1d0))::br,r
      REAL(KIND(1d0)),DIMENSION(MAX_DecayChain)::braccum
      INTEGER::ipip
      !INCLUDE "pp_psipsi_dps.inc"
      wgtbr=1d0
      ndecay=0
      NEVHEP=1
      NHEP=nhad
      DO i=1,nhad
         ISTHEP(i)=0
         IDHEP(i)=pdgt(iflh(i))
         IF(i.GT.4.AND.i.LE.6)THEN
            ISHEP(i)=1
            JMOHEP(1,i)=1
            JMOHEP(2,i)=2
            JDAHEP(1,i)=0
            JDAHEP(2,i)=0
         ELSEIF(i.GT.6.AND.i.LE.nhad)THEN
            ISHEP(i)=1
            JMOHEP(1,i)=3
            JMOHEP(2,i)=4
            JDAHEP(1,i)=0
            JDAHEP(2,i)=0
         ELSEIF(i.GT.2)THEN
            ISHEP(i)=-1
            JMOHEP(1,i)=0
            JMOHEP(2,i)=0
            JDAHEP(1,i)=7
            JDAHEP(2,i)=8
         ELSE
            ISHEP(i)=-1
            JMOHEP(1,i)=0
            JMOHEP(2,i)=0
            JDAHEP(1,i)=5
            JDAHEP(2,i)=6
         ENDIF
         PHEP(1:5,i)=dps_hadron_pmom(i,1:5)
         VHEP(1:4,i)=0d0
      ENDDO
      ipip=0
      dps_decay_index=0
      DO i=1,nhad
         j=iflh2DecayChains(i,0)
         IF(ABS(iflh(i)).GT.100)THEN
            ipip=ipip+1
            dps_decay_index=ipip
         ENDIF
         IF(j.GT.0)THEN
            br=0d0
            braccum(1:MAX_DecayChain)=0d0
            DO k=1,j
               m=iflh2DecayChains(i,k)
               br=br+DecayBR(m)
               braccum(k)=br
            ENDDO
            ! randomly select the decay process
            IF(j.EQ.1)THEN
               kk=1
            ELSE
               ndecay=ndecay+1
               r=Decayran(ndecay)*br
               DO k=1,j
                  IF(braccum(k).GT.r)THEN
                     kk=k
                     EXIT
                  ELSEIF(k.EQ.j)THEN
                     kk=j
                     EXIT
                  ENDIF
               ENDDO
            ENDIF
            m=iflh2DecayChains(i,kk)
            CALL HO_One_Decay_pp_psipsi_dps(i,m)
            wgtbr=wgtbr*br
         ENDIF
         IF(dps_decay_index.GT.2)THEN
            WRITE(*,*)"ERROR:too many psi (>2)"
            STOP
         ELSEIF(dps_decay_index.GT.0)THEN
            dps_decay_index=0
         ENDIF
      ENDDO
      RETURN
    END SUBROUTINE HO_Decay_pp_psipsi_dps

    SUBROUTINE HO_One_Decay_pp_psipsi_dps(imoth,idecay)
      USE HOVll
      USE Helac_ranmar_mod
      IMPLICIT NONE
      INTEGER,INTENT(IN)::imoth,idecay
      INTEGER,PARAMETER::available_ndecay=1
      CHARACTER(len=20),DIMENSION(available_ndecay)::available_decay
      INTEGER::i
      INTEGER::lavail
      REAL(KIND(1d0)),DIMENSION(3)::svec
      REAL(KIND(1d0)),DIMENSION(0:3)::PM,PD1,PD2
      REAL(KIND(1d0))::probT
      INCLUDE "stdhep_pp_psipsi_dps.inc"
      !INCLUDE "pp_psipsi_dps.inc"
      CALL HO_judge_avail_decay(imoth,idecay,lavail)
      IF(lavail.GT.0)THEN
         SELECT CASE(lavail)
            CASE(1)
               ! 3S11 > l+ l-
               IF(dps_decay_index.NE.1.AND.dps_decay_index.NE.2)THEN
                  WRITE(*,*)"ERROR:unknow which psi to decay"
                  STOP
               ENDIF
               PM(0)=dps_hadron_pmom(imoth,4)
               PM(1:3)=dps_hadron_pmom(imoth,1:3)
               svec(1:3)=PM(1:3)
               ! unpolarized in any frame
               IF(dps_lambdath(dps_decay_index).EQ.0d0)THEN
                  CALL HO_Vll_unpolarized(PM,PD1,PD2)
               ELSEIF(dps_lambdath(dps_decay_index).GE.1d0)THEN
                  CALL HO_Vll11(PM,svec,PD1,PD2)
               ELSEIF(dps_lambdath(dps_decay_index).LE.-1d0)THEN
                  CALL HO_Vll00(PM,svec,PD1,PD2)
               ELSE
                  probT=2d0*(1d0+dps_lambdath(dps_decay_index))/(3d0+dps_lambdath(dps_decay_index))
                  IF(Helac_rnmy(0).LT.probT)THEN
                     CALL HO_Vll11(PM,svec,PD1,PD2)
                  ELSE
                     CALL HO_Vll00(PM,svec,PD1,PD2)
                  ENDIF
               ENDIF
               ISTHEP(imoth)=1
               ISHEP(NHEP+1)=1
               ISHEP(NHEP+2)=1
               ISTHEP(NHEP+1)=0
               ISTHEP(NHEP+2)=0
               IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
               IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
               JMOHEP(1,NHEP+1:NHEP+2)=imoth
               JMOHEP(2,NHEP+1:NHEP+2)=0
               JDAHEP(1,imoth)=NHEP+1
               JDAHEP(2,imoth)=NHEP+2
               PHEP(5,NHEP+1:NHEP+2)=0d0
               PHEP(1:3,NHEP+1)=PD1(1:3)
               PHEP(1:3,NHEP+2)=PD2(1:3)
               PHEP(4,NHEP+1)=PD1(0)
               PHEP(4,NHEP+2)=PD2(0)
               NHEP=NHEP+2
            END SELECT
         ELSE
            CALL HO_avail_decays(available_ndecay,available_decay)
            PRINT *,"Only the following decay processes are available in HELAC-Onia"
            DO i=1,available_ndecay
               PRINT *,available_decay(i)
            ENDDO
            STOP
         ENDIF
         RETURN
    END SUBROUTINE HO_One_Decay_pp_psipsi_dps

    SUBROUTINE pp_psipsi_dps_hadronmom
      IMPLICIT NONE
      !INCLUDE "pp_psipsi_dps.inc"
      dps_hadron_pmom(1,1:4)=dps_pmom(1,1,1:4)
      dps_hadron_pmom(1,5)=0d0
      dps_hadron_pmom(2,1:4)=dps_pmom(1,2,1:4)
      dps_hadron_pmom(2,5)=0d0
      dps_hadron_pmom(3,1:4)=dps_pmom(2,1,1:4)
      dps_hadron_pmom(3,5)=0d0
      dps_hadron_pmom(4,1:4)=dps_pmom(2,2,1:4)
      dps_hadron_pmom(4,5)=0d0
      dps_hadron_pmom(5,1:4)=dps_pmom(1,3,1:4)
      dps_hadron_pmom(5,5)=mpsi(1)
      dps_hadron_pmom(6,1:4)=dps_pmom(1,4,1:4)
      dps_hadron_pmom(6,5)=0d0
      dps_hadron_pmom(7,1:4)=dps_pmom(2,3,1:4)
      dps_hadron_pmom(7,5)=mpsi(2)
      dps_hadron_pmom(8,1:4)=dps_pmom(2,4,1:4)
      dps_hadron_pmom(8,5)=0d0
      RETURN
    END SUBROUTINE pp_psipsi_dps_hadronmom

    SUBROUTINE cuts_Decay_pp_psipsi_dps(icut)
      ! only cut on the decay particles                                                                
      IMPLICIT NONE
      INTEGER,INTENT(OUT)::icut
      INTEGER::i,kk
      REAL(KIND(1d0))::eta,pt,c
      REAL(KIND(1d0))::ptminmuon
      iNTEGER::flag
      INTEGER::muon4psiATLAS
      INCLUDE "stdhep_pp_psipsi_dps.inc"
      icut=0
      flag=0
      muon4psiATLAS=0
      DO i=1,NHEP
         IF(ISHEP(i).EQ.-1)CYCLE ! only cut on the final state                                               
         IF(ISTHEP(i).EQ.1)CYCLE ! exclude the mother particles, which have been cutted before               
         IF(JMOHEP(1,i).LE.0.OR.JMOHEP(1,i).GT.NHEP)CYCLE
         IF(ISTHEP(JMOHEP(1,i)).NE.1.OR.ISHEP(JMOHEP(1,i)).EQ.-1)CYCLE
         ! cut is only applied on the decay particles                                                        
         IF(lepton_pdg(IDHEP(i)))THEN
            ! cut on the leptons                                                                             
            pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
            IF(pt.LT.ho_dptcl)THEN
               flag=1
               EXIT
            ENDIF
            eta=prapidity(PHEP(1:4,i))
            IF(ABS(eta).GT.ho_detacl)THEN
               flag=2
               EXIT
            ENDIF
            eta=rapidity(PHEP(1:4,i))
            IF(absrap)eta=ABS(eta)
            IF(eta.GT.ho_dycl.OR.eta.LT.ho_dycllow)THEN
               flag=3
               EXIT
            ENDIF
            IF(PHEP(4,i).LT.ho_decl)THEN
               flag=4
               EXIT
            ENDIF
            c=PHEP(3,i)/PHEP(4,i)
            IF(ABS(c).GT.ho_dcl)THEN
               flag=5
               EXIT
            ENDIF
            ! special cutoffs                                                                                
            IF(literature_cutoffs.EQ.14062380)THEN
               ! DZero (arXiv:1406.2380) special cutoffs                                                     
               eta=prapidity(PHEP(1:4,i))
               IF(ABS(eta).LT.1.35d0)THEN
                  pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                  IF(pt.LE.2d0)THEN
                     flag=6
                     EXIT
                  ENDIF
               ELSEIF(ABS(eta).LT.2d0.AND.ABS(eta).GT.1.35d0)THEN
                  IF(PHEP(4,i).LE.4d0)THEN
                     flag=7
                     EXIT
                  ENDIF
               ELSE
                  flag=8
                  EXIT
               ENDIF
            ELSEIF(literature_cutoffs.EQ.14060000)THEN
               ! ATLAS condition, speical cutoffs                                                            
               ! there is no arXiv number now                                                                
               ! At least one of the Jpsi must have two muons with pT>4 GeV each                             
               ! muon pT>2.5 GeV and muon |eta|<2.3                                                          
               eta=prapidity(PHEP(1:4,i))
               IF(ABS(eta).GE.2.3d0)THEN
                  flag=9
                  EXIT
               ENDIF
               pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
               IF(pt.LE.2.5d0)THEN
                  flag=10
                  EXIT
               ENDIF
               IF(IDHEP(JMOHEP(1,i)).EQ.443.AND.muon4psiATLAS.LT.2)THEN
                  IF(pt.GT.4d0)THEN
                     muon4psiATLAS=muon4psiATLAS+1
                  ELSE
                     muon4psiATLAS=muon4psiATLAS-1
                  ENDIF
               ENDIF
               IF(muon4psiATLAS.EQ.-2)muon4psiATLAS=0
            ELSEIF(literature_cutoffs.EQ.14070000)THEN
               ! CMS condition (from Junquan Tao for psi+psi) 
               eta=prapidity(PHEP(1:4,i))
               IF(ABS(eta).LT.1.2d0)THEN
                  pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                  IF(pt.LE.3.5d0)THEN
                     flag=12
                     EXIT
                  ENDIF
               ELSEIF(ABS(eta).GT.1.2d0.AND.ABS(eta).LT.1.6d0)THEN
                  ptminmuon=3.5d0-(ABS(eta)-1.2d0)*1.5d0/0.4d0
                  pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                  IF(pt.LE.ptminmuon)THEN
                     flag=13
                     EXIT
                  ENDIF
               ELSEIF(ABS(eta).GT.1.6d0.AND.ABS(eta).LT.2.4d0)THEN
                  pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                  IF(pt.LE.2.0d0)THEN
                     flag=14
                     EXIT
                  ENDIF
               ELSE
                  flag=15
                  EXIT
               ENDIF
            ENDIF
         ENDIF
      ENDDO
      IF(flag.EQ.0.AND.literature_cutoffs.EQ.14060000)THEN
         ! ATLAS condition, speical cutoffs                                                                  
         ! there is no arXiv number now                                                                      
         ! At least one of the Jpsi must have two muons with pT>4 GeV each                                   
         IF(muon4psiATLAS.LT.2)THEN
            flag=11
         ENDIF
      ENDIF
      IF(flag.EQ.0)THEN
         icut=1
      ENDIF
      RETURN
    END SUBROUTINE cuts_Decay_pp_psipsi_dps

    SUBROUTINE unwei_procedure_pp_psipsi_dps(w1,nwri,nwmax,nwarn)
      IMPLICIT NONE
      REAL(KIND(1d0)),INTENT(IN)::w1
      INTEGER,INTENT(INOUT)::nwri,nwmax,nwarn
      REAL(KIND(1d0)),DIMENSION(1)::ranr
      REAL(KIND=DBL)::vtime,vspin,scale1,scalup,xwgtup,px,py,pz,p0,pmass,umax,umax1
      INTEGER::init=0,i
      INTEGER::idup,idprup,istup,imothup1,imothup2,icol1,icol2
      LOGICAL::llwri
      !INCLUDE "pp_psipsi_dps.inc"
      SAVE init,umax,umax1
      IF(init.EQ.0)THEN
         umax=0d0
         umax1=0d0
         init=1
      ENDIF
      IF(lwmax.EQ.1)THEN
         IF(umax.LT.w1)THEN
            umax=w1
         ENDIF
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
               CALL unwei_writer_Decay_pp_psipsi_dps
               Nevents=nwri
               RETURN
            ENDIF
            idprup=82 ! id for the process 
            xwgtup=1 !w1*10**3 !1  
            scalup=dps_pmom(1,3,1)**2+dps_pmom(1,3,2)**2+mpsi(1)**2
            scalup=DSQRT(scalup)
            WRITE(nunit3)nhad,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
            DO i=1,nhad
               idup=pdgt(iflh(i))
               IF(i.LE.4)THEN
                  istup=-1
               ELSE
                  istup=1
               ENDIF
               imothup1=0
               imothup2=0
               IF(i.EQ.5.OR.i.EQ.6)THEN
                  imothup1=1
                  imothup2=2
               ELSEIF(i.EQ.7.OR.i.EQ.8)THEN
                  imothup1=3
                  imothup2=4
               ENDIF
               IF(i.LE.4)THEN
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
               px=dps_hadron_pmom(i,1)
               py=dps_hadron_pmom(i,2)
               pz=dps_hadron_pmom(i,3)
               p0=dps_hadron_pmom(i,4)
               pmass=dps_hadron_pmom(i,5)
               vtime=0
               vspin=9

               WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2&
                    ,px,py,pz,p0,pmass,vtime,vspin
            ENDDO
         ENDIF
      ENDIF
      Nevents=nwri
    END SUBROUTINE unwei_procedure_pp_psipsi_dps

    SUBROUTINE unwei_writer_Decay_pp_psipsi_dps
      IMPLICIT NONE
      INCLUDE "stdhep_pp_psipsi_dps.inc"
      INTEGER::icol1,icol2,idup,idprup,istup,imothup1,imothup2
      REAL(KIND(1d0))::px,py,pz,p0,pmass,scalup,vtime,vspin,xwgtup
      INTEGER::i
      !INCLUDE "pp_psipsi_dps.inc"
      idprup = 82
      xwgtup=1
      scalup=dps_pmom(1,3,1)**2+dps_pmom(1,3,2)**2+mpsi(1)**2
      scalup=DSQRT(scalup)
      WRITE(nunit3)NHEP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
      DO i=1,NHEP
         idup=IDHEP(i)
         istup=ISHEP(i)
         imothup1=JMOHEP(1,i)
         imothup2=JMOHEP(2,i)
         IF(i.GT.nhad)THEN
            ! no hadronic decay
            icol1=0
            icol2=0
         ELSE
            IF(istup.LT.0)THEN
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
         ENDIF
         px=PHEP(1,i)
         py=PHEP(2,i)
         pz=PHEP(3,i)
         p0=PHEP(4,i)
         pmass=PHEP(5,i)
         vtime=0
         vspin=9
         WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2,&
              px,py,pz,p0,pmass,vtime,vspin
      ENDDO
    END SUBROUTINE unwei_writer_Decay_pp_psipsi_dps

    SUBROUTINE Generate_lhe_pp_psipsi_dps(n1,nevent,icase)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::n1,nevent,icase
      CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/",output_dir="./output/"
      INTEGER::i,nunit4,nunit5,nunit3
      INTEGER::istop,k
      REAL(KIND(1d0))::p0,px,py,pz,SPINUP,EBMUP1,EBMUP2,XSECUP1,XSECUP2,XERRUP1,XMAXUP1,&
           XWGTUP,VTIMUP,SCALUP,PM0,AQEDUP,AQCDUP,XWGTUP2
      INTEGER::IDBMUP1,IDBMUP2,IDWTUP,NPRUP,IDPRUP,NUP,IDUP,ISTUP,IMOTHUP1,IMOTHUP2,&
           ICOLUP1,ICOLUP2,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,LPRUP1
      INTEGER::ipip
      nunit3=30
      CLOSE(nunit3)
      OPEN(nunit3,FILE=TRIM(tmp_dir)//'even_pp_psipsi_dps.out',FORM='unformatted')
      nunit4=31
      CLOSE(nunit4)
      OPEN(nunit4,FILE=TRIM(tmp_dir)//'sample_pp_psipsi_dps.init')
      nunit5=32
      CLOSE(nunit5)
      OPEN(nunit5,FILE=TRIM(output_dir)//'sample_pp_psipsi_dps.lhe')
      WRITE(nunit5,'(A)') '<LesHouchesEvents version="1.0">'
      WRITE(nunit5,'(A)') '<!--'
      WRITE(nunit5,'(A)') 'File generated with HELAC-ONIA '
      WRITE(nunit5,'(A)') '-->'
      WRITE(nunit5,'(A)') '<init>'
      READ(nunit4,*)IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
      READ(nunit4,*)XSECUP1,XERRUP1,XMAXUP1,LPRUP1
      WRITE(nunit5,5000)IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
      WRITE(nunit5,5100)XSECUP1,XERRUP1,XMAXUP1,LPRUP1
      WRITE(nunit5,'(A)') '</init>'
      istop=1
      k=0
      XWGTUP2=XSECUP1/nevent
      DO WHILE(istop.EQ.1)
         k=k+1
         READ(nunit3,END=100) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
         IF(icase.EQ.1)XWGTUP=XWGTUP2
         WRITE(nunit5,'(A)') '<event>'
         WRITE(nunit5,5200) NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
         ipip=1
         DO i=1,NUP
            READ(nunit3)IDUP,ISTUP,iMOTHUP1,iMOTHUP2,ICOLUP1,ICOLUP2,px,py,pz,p0,pm0,VTIMUP,SPINUP
            IF(IDUP.EQ.443.OR.IDUP.EQ.553)THEN
               IF(ipip.GT.2)THEN
                  WRITE(*,*)"ERROR:Too many psi(Upsilon) found (>2) !"
                  STOP
               ENDIF
               IF(dps_istate(ipip).EQ.2)THEN
                  IDUP=100443 ! psi(2S)
               ELSEIF(dps_istate(ipip).EQ.4)THEN
                  IDUP=100553 ! Y(2S)
               ELSEIF(dps_istate(ipip).EQ.5)THEN
                  IDUP=200553 ! Y(3S)
               ENDIF
               ipip=ipip+1
            ENDIF
            WRITE(nunit5,5300)IDUP,ISTUP,iMOTHUP1,iMOTHUP2,ICOLUP1,ICOLUP2,px,py,pz,p0,pm0,VTIMUP,SPINUP
         ENDDO
         WRITE(nunit5,'(A)') '</event>'
      ENDDO
100   CONTINUE
      IF(icase.EQ.1.AND.k.NE.nevent+1)THEN
         WRITE(*,*)"WARNING:mismatching of the unweighted lhe events number ",k-1,nevent
      ENDIF
      WRITE(nunit5,'(A)') '</LesHouchesEvents>'
      CLOSE(nunit3,STATUS='delete')
      CLOSE(nunit4,STATUS='delete')
5200  FORMAT(1P,2I6,4E14.6)
5300  FORMAT(1P,I8,5I5,5E18.10,E14.6,E12.4)
5000  FORMAT(1P,2I8,2E14.6,6I8)
5100  FORMAT(1P,3E20.10,I6)
      CLOSE(nunit5,STATUS='keep')
    END SUBROUTINE Generate_lhe_pp_psipsi_dps
END MODULE pp_psipsi_DPS
