MODULE pp_aajj_DPS
  USE Helac_Global
  USE pp_aajj_dps_global
  USE Helac_Func_1
  USE Kinetic_Func
  USE MC_VEGAS
  USE Constants
  USE ME2_pp_aj
  USE ME2_pp_aa
  USE ME2_pp_jj
  USE plot_pp_aajj_dps
  USE Func_PSI
  USE pp_aajj_dps_func
  USE pp_aajj_DPS_cuts
  IMPLICIT NONE
  INTEGER::nprint
  INTEGER::varnum,nunwei
  LOGICAL::lunwei2
  REAL(KIND(1d0))::EBMUP1,EBMUP2
  INTEGER::NPRUP=0,lwmax
  INTEGER,PARAMETER::maxprint=8
  !REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
  SAVE
CONTAINS
  SUBROUTINE calc_pp_aajj_DPS
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
    REAL(KIND(1d0))::rwtmp
    REAL(KIND(1d0))::rsw2,alpha
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
    INTEGER::nnproc,iiproc
    LOGICAL::hasjjsubproc=.FALSE.,hasajsubproc=.FALSE.,hasaasubproc=.FALSE.
    INTEGER::irslt
    SAVE lunwei,lhewgt,icase
    WRITE(*,*)'                     THE BEGINNING OF HELAC-Onia'
    WRITE(*,*)'    AddOn Process: '
    WRITE(*,*)'    Double Patron Scattering for p p > diphoton+dijet+X '
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
          rwtmp=rw_Rscale_down
          rw_Rscale_down=rw_Rscale_up
          rw_Rscale_up=rwtmp
       ENDIF
       IF(reweight_scale.AND.rw_Fscale_down.GT.rw_Fscale_up)THEN
          rwtmp=rw_Fscale_down
          rw_Fscale_down=rw_Fscale_up
          rw_Fscale_up=rwtmp
       ENDIF
    ENDIF
    CALL ReadConst
    IF(alphaem.EQ.-1)THEN
       rsw2=1d0-rwm**2/rzm**2
       alpha=DSQRT(dnou(2))*gfermi*rwm*rwm*rsw2/pi
    ELSE
       alpha=alphaem
    ENDIF
    EL=DSQRT(alpha*4*pi)
    Gstrong=DSQRT(alphaQCD2*4*pi)
    OPEN(UNIT=12321,FILE="./input/subprocess.inp")
    nnproc = 0
    !  go through the file to count the number of entries
    DO
       READ(12321,*, iostat = irslt)i
       IF(irslt.NE.0)EXIT
       nnproc = nnproc + 1
    ENDDO
    ! subtract the first line
    nnproc=nnproc-1
    REWIND(12321)
    READ(12321,*)dpsorsps
    IF(dpsorsps.NE.1.AND.dpsorsps.NE.2)THEN
       WRITE(*,*)"ERROR:Incorrect (first) number in the first line of subprocess.inp"
       STOP
    ENDIF
    contained_subprocess(1:maxproc_dps)=.FALSE.
    hasjjsubproc=.FALSE.
    hasajsubproc=.FALSE.
    hasaasubproc=.FALSE.
    DO i=1, nnproc
       READ(12321,*)iiproc
       IF(iiproc.EQ.1001)THEN
          contained_subprocess(1:176)=.TRUE.
          hasjjsubproc=.TRUE.
          CYCLE
       ENDIF
       IF(iiproc.EQ.1002)THEN
          contained_subprocess(177:206)=.TRUE.
          hasajsubproc=.TRUE.
          CYCLE
       ENDIF
       IF(iiproc.EQ.1003)THEN
          contained_subprocess(207:217)=.TRUE.
          hasaasubproc=.TRUE.
          CYCLE
       ENDIF
       IF(iiproc.LE.0.OR.iiproc.GT.maxproc_dps)THEN
          WRITE(*,*)"WARNING:Unknown processs with number = ",iiproc
          CYCLE
       ENDIF
       contained_subprocess(iiproc)=.TRUE.
       IF(iiproc.GE.1.AND.iiproc.LE.176)hasjjsubproc=.TRUE.
       IF(iiproc.GE.177.AND.iiproc.LE.206)hasajsubproc=.TRUE.
       IF(iiproc.GE.207.AND.iiproc.LE.217)hasaasubproc=.TRUE.
    ENDDO
    CLOSE(UNIT=12321)
    IF(dpsorsps.EQ.1)THEN
       ! DPS
       IF(.NOT.hasajsubproc.AND..NOT.(hasjjsubproc.AND.hasaasubproc))THEN
          WRITE(*,*)"No valid subprocess in subprocess.inp for p p > a a j j"
          STOP
       ENDIF
    ELSE
       ! SPS
       IF(.NOT.hasjjsubproc.AND..NOT.hasajsubproc.AND..NOT.hasaasubproc)THEN
          WRITE(*,*)"No valid subprocess in subprocess.inp for p p > a a or j j or a j"
          STOP
       ENDIF
    ENDIF
   ! some necesarry initialization
    OPEN(UNIT=12321,FILE="./input/sigma_eff.inp")
    READ(12321,*)dps_sigmaeff
    CLOSE(UNIT=12321)
    ! convert mb to nb
    dps_sigmaeff=dps_sigmaeff*1d6
    ! read K factors
    OPEN(UNIT=12321,FILE="./input/Kfactors.inp")
    READ(12321,*)Kfactor_jj,Kfactor_aj,Kfactor_aa
    CLOSE(UNIT=12321)
    IF(dpsorsps.EQ.1)THEN
       ! DPS
       nhad=8
    ELSE
       ! SPS
       nhad=4
    ENDIF
    IF(nunit1.NE.6)THEN
       OPEN(UNIT=nunit1,FILE=TRIM(output_dir)//"RESULT_pp_aajj_dps.out")
    ENDIF
    nunit2=32
    CLOSE(nunit2)
    OPEN(nunit2,FILE=TRIM(output_dir)//'kine_pp_aajj_dps.out')
    nunit3=30
    CLOSE(nunit3)
    OPEN(nunit3,FILE=TRIM(tmp_dir)//'even_pp_aajj_dps.out',FORM='unformatted')
    ! for LHA
    CLOSE(200)
    OPEN(200,FILE=TRIM(tmp_dir)//'sample_pp_aajj_dps.init')
    CALL ReadElem_integer("itmax",itmxn)
    CALL ReadElem_integer("nmc",ncalln)
    WRITE(*,*)' '
    CALL Helac_mtime()
    WRITE(*,*)' '
    IF(dpsorsps.EQ.1)THEN
       ! DPS
       CALL pp_aajj_DPS_VEGAS(rslt,itmxn,ncalln)
    ELSE
       ! SPS
       CALL pp_jj_aj_aa_VEGAS(rslt,itmxn,ncalln)
    ENDIF
    WRITE(nunit1,*)"sigma(nb)                   sd"
    WRITE(nunit1,*)rslt(1),rslt(2)
    WRITE(*,*)' '
    CALL Helac_mtime()
    CLOSE(nunit2)
    CLOSE(nunit3)
    CLOSE(21)
    CLOSE(200)
    CALL Generate_lhe_pp_aajj_dps(nhad,Nevents,icase)
  END SUBROUTINE calc_pp_aajj_DPS

  SUBROUTINE pp_aajj_DPS_VEGAS(rslt,itmxn,ncalln)
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
    varnum=2*varnum
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
    CALL VEGAS(varnum,pp_aajj_DPS_fxn,vfes,sd,chi2a)
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
       CALL VEGAS(varnum,pp_aajj_DPS_fxn,vfes,sd,chi2a,1)
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
         IF(plot_output)CALL initplot_pp_aajj_DPS
         CALL VEGAS(varnum,pp_aajj_DPS_fxn,vfes,sd,chi2a,1)
         IF(plot_output)CALL plotout_pp_aajj_DPS
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
         IF(plot_output)CALL initplot_pp_aajj_DPS
         CALL VEGAS(varnum,pp_aajj_DPS_fxn,vfes,sd,chi2a,1)
         IF(plot_output)CALL plotout_pp_aajj_DPS
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
    END SUBROUTINE pp_aajj_DPS_VEGAS

    SUBROUTINE pp_jj_aj_aa_VEGAS(rslt,itmxn,ncalln)
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
      CALL VEGAS(varnum,pp_jj_aj_aa_fxn,vfes,sd,chi2a)
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
         CALL VEGAS(varnum,pp_jj_aj_aa_fxn,vfes,sd,chi2a,1)
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
         !IF(plot_output)CALL initplot_pp_aajj_DPS
         CALL VEGAS(varnum,pp_jj_aj_aa_fxn,vfes,sd,chi2a,1)
         !IF(plot_output)CALL plotout_pp_aajj_DPS
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
         !IF(plot_output)CALL initplot_pp_aajj_DPS
         CALL VEGAS(varnum,pp_jj_aj_aa_fxn,vfes,sd,chi2a,1)
         !IF(plot_output)CALL plotout_pp_aajj_DPS
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
      WRITE(200,6100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&
           iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
      WRITE(200,6200) vfes*10d0**3,sd*10d0**3,1d0, 82
      rslt(1)=vfes
      rslt(2)=sd
      rslt(3)=chi2a
      IF(lunwei2)PRINT *,"number of events",Nevents
      RETURN
6100  FORMAT(1P,2I8,2E14.6,6I8)
6200  FORMAT(1P,3E20.10,I6)
    END SUBROUTINE pp_jj_aj_aa_VEGAS

    FUNCTION pp_aajj_DPS_fxn(x,wgt)
      REAL(KIND(1d0)),DIMENSION(varnum),INTENT(IN)::x
      REAL(KIND(1d0)),INTENT(IN)::wgt
      REAL(KIND(1d0))::pp_aajj_DPS_fxn
      REAL(KIND(1d0))::y1,y2,phi,w1
      REAL(KIND(1d0)),DIMENSION(4)::pmom
      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
      INTEGER::init=0,nwarn,nwmax,nwri,nwri_tot,icut,nnn=0,nnntot=0,&
           pdfnumpdf,ivarnum
      SAVE init,nwarn,nwmax,nwri,nwri_tot,nnn,nnntot,pdfnumpdf,ivarnum
      REAL(KIND(1d0))::sqs,sq,pt1c,y1cup,y1clow,ycollcm,mt1,mt2,pt
      SAVE sqs,sq,pt1c,y1cup,y1clow,ycollcm
      REAL(KIND(1d0))::wdps,wsf,wme,wgtbr,temp1,temp2,temp21,temp3,temp4
      REAL(KIND(1d0))::exp1,exp2,Jpt2,Jy1,Jy2,xp11,xp21,scale1,scale2
      REAL(KIND(1d0))::recmax=0,recmin=0
      SAVE temp1
      INTEGER::ipip,ioffset,j,iproc1,iproc2
      LOGICAL::allowed
      REAL(KIND(1d0))::wme1,wme2,wsf1,wsf2
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
         CALL ReadElem_real('maxy1c',y1cup)
         CALL ReadElem_real('miny1c',y1clow)
         CALL ReadElem_integer('pdf',pdfnumpdf)
         ptc(3:20)=0
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
            CALL readcuts_pp_aajj_dps
            !y1cup=MIN(y1cup,yycut(3))
            !y1clow=MAX(y1clow,yycutlow(3))
            IF(absrap)y1cup=ABS(y1cup)
            IF(absrap)y1clow=ABS(y1clow)
            CALL readcuts_pp_aajj_dps_cuts
         ENDIF
         sq=sqs*sqs
         !temp1=(sq-mpsi**2)/(2d0*sqs)
         ivarnum=3
         IF(lunwei2)ivarnum=ivarnum+1
         ! determine the process
         DO iproc1=1,maxproc_dps
            CALL determine_subprocess(1,iproc1)
            DO iproc2=iproc1,maxproc_dps
               CALL determine_subprocess(2,iproc2)
               CALL subprocess_dps(allowed)
               IF(.NOT.allowed)THEN
                  iflh_save(iproc1,iproc2,1:8)=0
               ELSE
                  iflh_save(iproc1,iproc2,1:8)=iflh(1:8)
               ENDIF
            ENDDO
         ENDDO
         init=1
      ENDIF
      nnntot=nnntot+1
      wdps=1d0
      DO ipip=1,2
         ! the ipip-th SPS
         temp1=(sq)/(2d0*sqs) 
         ioffset=ivarnum*(ipip-1)
         pt=(temp1-pt1c)*x(3+ioffset)+pt1c
         mt1=pt
         mt2=pt
         IF(absrap)THEN
            temp2=MIN(ACosh_p((sq)/(2d0*sqs*mt1)),y1cup)
            IF(temp2.LT.y1clow)THEN
               pp_aajj_DPS_fxn=0d0
               RETURN
            ENDIF
            y1=(2d0*x(1+ioffset)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1+ioffset)-0.5d0)*y1clow
         ELSE
            temp2=MIN(ACosh_p((sq)/(2d0*sqs*mt1)),y1cup-ycollcm)
            temp21=MAX(-ACosh_p((sq)/(2d0*sqs*mt1)),y1clow-ycollcm)
            IF(temp2.LT.temp21)THEN
               pp_aajj_DPS_fxn=0d0
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
               pp_aajj_DPS_fxn=0d0
               RETURN
            ENDIF
            IF((DEXP(-y1)*mt1+DEXP(-y2)*mt2)/sqs+xp2.GE.1d0)THEN
               pp_aajj_DPS_fxn=0d0
               RETURN
            ENDIF
         ENDIF
         xp1=(DEXP(y1)*mt1+DEXP(y2)*mt2)/sqs
         xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2)/sqs
         ehat=DSQRT(xp1*xp2*sq)
         wdps=wdps*xp1*xp2*Jpt2*Jy1*Jy2*3.8938573d5/(16d0*pi*ehat**4)
         xpi(ipip,1)=xp1
         xpi(ipip,2)=xp2
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
         !IF(pdfnumpdf.EQ.921000)THEN
            ! GS 09 dPDF
         !   IF(ipip.EQ.1)THEN
         !      xp11=xp1
         !      xp21=xp2
         !      scale1=dps_pmom(ipip,3,1)**2+dps_pmom(ipip,3,2)**2
         !      scale1=DSQRT(scale1)
         !   ELSE
         !      scale2=dps_pmom(ipip,3,1)**2+dps_pmom(ipip,3,2)**2
         !      scale2=DSQRT(scale2)
         !      CALL GS09(xp11,xp1,DSQRT(scale1*scale2),0,0,wsf)
         !      wdps=wdps*wsf
         !      CALL GS09(xp21,xp2,DSQRT(scale1*scale2),0,0,wsf)
         !      wdps=wdps*wsf
               ! CALL GSALPS(Q) to run alpha_S, which is not used here
         !   ENDIF
         !ELSE
         !   CALL strf_pdf_pp_aajj_dps(ipip,wsf)
         !   wdps=wdps*wsf
         !ENDIF
      ENDDO
      IF(wdps.LE.0d0)THEN
         pp_aajj_DPS_fxn=0d0
         RETURN
      ENDIF
      ! exhaust all possible subprocesses
      nonzero_subprocess(1:maxproc_dps,1:maxproc_dps)=.FALSE.
      wme_subprocess(1:maxproc_dps,1:maxproc_dps)=0d0
      wme_calc(1:2,1:maxproc_dps)=.FALSE.
      wme_proc(1:2,1:maxproc_dps)=0d0
      IF(pdfnumpdf.EQ.921000)THEN
         ! GS 09 dPDF
         strf_dpdf_save(1:2,-5:5,-5:5)=0d0 ! dpdf
         strf_dpdf_saveq(1:2,-5:5,-5:5)=.FALSE.
      ELSE
         ! normal sPDF
         strf_pdf_save(1:2,1:2,-5:5)=0d0 ! pdf*x
         strf_pdf_saveq(1:2)=.FALSE.
      ENDIF
      IF(irun.EQ.1)THEN
         ! running of alpha_s
         CALL RunAlphas(pdfnumpdf)
      ENDIF
      wmetot=0d0
      ! wme is pdf*ME2 and symmetry factor for the two dps subprocesses
      DO iproc1=1,maxproc_dps
         IF(.NOT.contained_subprocess(iproc1))CYCLE
         CALL determine_subprocess(1,iproc1)
         ! determine s and t
         CALL calc_st(1)
         CALL ME2_subprocess(1,iproc1,wme1)
         IF(wme1.LE.0d0)CYCLE
         ! multiply K factor
         IF(iproc1.GE.1.AND.iproc1.LE.176)THEN
            wme1=wme1*Kfactor_jj
         ELSEIF(iproc1.GE.177.AND.iproc1.LE.206)THEN
            wme1=wme1*Kfactor_aj
         ELSEIF(iproc1.GE.207.AND.iproc1.LE.217)THEN
            wme1=wme1*Kfactor_aa
         ENDIF
         IF(wme1.LE.0d0)CYCLE
         !IF(pdfnumpdf.NE.921000)THEN
            ! it uses normal sPDF
         !   CALL strf_pdf_subprocess(pdfnumpdf,wsf1)
         !   wme1=wme1*wsf1
         !   IF(wme1.LE.0d0)CYCLE
         !ENDIF
         DO iproc2=iproc1,maxproc_dps
            IF(.NOT.contained_subprocess(iproc2))CYCLE
            CALL determine_subprocess(2,iproc2)
            CALL subprocess_dps(allowed)
            IF(.NOT.allowed)CYCLE
            ! now we have a+a+j+j
            CALL pp_aajj_dps_hadronmom
            icut=1
            CALL Cuts_pp_aajj_dps(icut)
            IF(icut.EQ.0)CYCLE
            CALL Cuts_pp_aajj_dps_cuts(icut)
            IF(icut.EQ.0)CYCLE
            ! determine s and t
            CALL calc_st(2)
            CALL ME2_subprocess(2,iproc2,wme2)
            IF(wme2.LE.0d0)CYCLE
            CALL strf_pdf_subprocess(pdfnumpdf,wsf2)
            wme2=wme2*wsf2
            IF(wme2.LE.0d0)CYCLE
            ! multiply K factor
            IF(iproc2.GE.1.AND.iproc2.LE.176)THEN
               wme2=wme2*Kfactor_jj
            ELSEIF(iproc2.GE.177.AND.iproc2.LE.206)THEN
               wme2=wme2*Kfactor_aj
            ELSEIF(iproc2.GE.207.AND.iproc2.LE.217)THEN
               wme2=wme2*Kfactor_aa
            ENDIF
            IF(wme2.LE.0d0)CYCLE
            wme=wme1*wme2
            ! the symmetry factor
            IF(iproc1.EQ.iproc2)wme=wme/2d0
            nonzero_subprocess(iproc1,iproc2)=.TRUE.
            wme_subprocess(iproc1,iproc2)=wme
            wmetot=wmetot+wme
         ENDDO
      ENDDO
      IF(wmetot.EQ.0d0)THEN
         pp_aajj_DPS_fxn=0d0
         RETURN
      ENDIF
      pp_aajj_DPS_fxn=wdps*wmetot/dps_sigmaeff
      IF(lunwei2.AND.lwmax.GT.0)THEN
         w1=pp_aajj_DPS_fxn*wgt ! multiply VEGAS weight 
         CALL unwei_procedure_pp_aajj_dps(w1,nwri,nwmax,nwarn)
      ENDIF
      IF(plot_output)THEN
         w1=pp_aajj_DPS_fxn*wgt
         CALL outfun_pp_aajj_dps(w1)
      ENDIF
      nnn=nnn+1
      IF(recmax.LT.pp_aajj_DPS_fxn)recmax=pp_aajj_DPS_fxn
      IF(recmin.GT.pp_aajj_DPS_fxn)recmin=pp_aajj_DPS_fxn
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
    END FUNCTION pp_aajj_DPS_fxn

    FUNCTION pp_jj_aj_aa_fxn(x,wgt)
      REAL(KIND(1d0)),DIMENSION(varnum),INTENT(IN)::x
      REAL(KIND(1d0)),INTENT(IN)::wgt
      REAL(KIND(1d0))::pp_jj_aj_aa_fxn
      REAL(KIND(1d0))::y1,y2,phi,w1
      REAL(KIND(1d0)),DIMENSION(4)::pmom
      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
      INTEGER::init=0,nwarn,nwmax,nwri,nwri_tot,icut,nnn=0,nnntot=0,&
           pdfnumpdf,ivarnum
      SAVE init,nwarn,nwmax,nwri,nwri_tot,nnn,nnntot,pdfnumpdf,ivarnum
      REAL(KIND(1d0))::sqs,sq,pt1c,y1cup,y1clow,ycollcm,mt1,mt2,pt
      SAVE sqs,sq,pt1c,y1cup,y1clow,ycollcm
      REAL(KIND(1d0))::wdps,wsf,wme,wgtbr,temp1,temp2,temp21,temp3,temp4
      REAL(KIND(1d0))::exp1,exp2,Jpt2,Jy1,Jy2,xp11,xp21,scale1,scale2
      REAL(KIND(1d0))::recmax=0,recmin=0
      SAVE temp1
      INTEGER::ipip,ioffset,j,iproc1,iproc2
      LOGICAL::allowed
      REAL(KIND(1d0))::wme1,wme2,wsf1,wsf2
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
         CALL ReadElem_real('maxy1c',y1cup)
         CALL ReadElem_real('miny1c',y1clow)
         CALL ReadElem_integer('pdf',pdfnumpdf)
         IF(pdfnumpdf.EQ.921000)THEN
            WRITE(*,*)"ERROR:Cannot use DPS when calculating 2 > 2 SPS"
            STOP
         ENDIF
         ptc(3:20)=0
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
            CALL readcuts_pp_aajj_dps
            !y1cup=MIN(y1cup,yycut(3))
            !y1clow=MAX(y1clow,yycutlow(3))
            IF(absrap)y1cup=ABS(y1cup)
            IF(absrap)y1clow=ABS(y1clow)
         ENDIF
         sq=sqs*sqs
         init=1
      ENDIF
      nnntot=nnntot+1
      wdps=1d0
      temp1=(sq)/(2d0*sqs)
      ioffset=0
      pt=(temp1-pt1c)*x(3+ioffset)+pt1c
      mt1=pt
      mt2=pt
      IF(absrap)THEN
         temp2=MIN(ACosh_p((sq)/(2d0*sqs*mt1)),y1cup)
         IF(temp2.LT.y1clow)THEN
            pp_jj_aj_aa_fxn=0d0
            RETURN
         ENDIF
         y1=(2d0*x(1+ioffset)-1d0)*(temp2-y1clow)+SIGN(1d0,x(1+ioffset)-0.5d0)*y1clow
      ELSE
         temp2=MIN(ACosh_p((sq)/(2d0*sqs*mt1)),y1cup-ycollcm)
         temp21=MAX(-ACosh_p((sq)/(2d0*sqs*mt1)),y1clow-ycollcm)
         IF(temp2.LT.temp21)THEN
            pp_jj_aj_aa_fxn=0d0
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
      xp1=(DEXP(y1)*mt1+DEXP(y2)*mt2)/sqs
      xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2)/sqs
      ehat=DSQRT(xp1*xp2*sq)
      wdps=wdps*xp1*xp2*Jpt2*Jy1*Jy2*3.8938573d5/(16d0*pi*ehat**4)
      xpi(1,1)=xp1
      xpi(1,2)=xp2
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
      dps_pmom(1,1,1:2)=0
      dps_pmom(1,1,3)=xp1*sqs/2d0
      dps_pmom(1,1,4)=xp1*sqs/2d0
      dps_pmom(1,2,1:2)=0
      dps_pmom(1,2,3)=-xp2*sqs/2d0
      dps_pmom(1,2,4)=xp2*sqs/2d0
      ! we choose phi=0 
      IF(lunwei2)THEN
         phi=2*pi*x(4+ioffset)
      ELSE
         phi=0d0
      ENDIF
      dps_pmom(1,3,1)=pt*DCOS(phi)
      dps_pmom(1,3,2)=pt*DSIN(phi)
      dps_pmom(1,3,3)=mt1*DSINH(y1)
      dps_pmom(1,3,4)=mt1*DCOSH(y1)
      dps_pmom(1,4,1)=-pt*DCOS(phi)
      dps_pmom(1,4,2)=-pt*DSIN(phi)
      dps_pmom(1,4,3)=mt2*DSINH(y2)
      dps_pmom(1,4,4)=mt2*DCOSH(y2)
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
            CALL Boostl(sqs,pmom,dps_pmom(1,j,1:4))
         ENDDO
      ENDIF
      dps_pmom(2,1:4,1:4)=dps_pmom(1,1:4,1:4)
      IF(wdps.LE.0d0)THEN
         pp_jj_aj_aa_fxn=0d0
         RETURN
      ENDIF
      ! exhaust all possible subprocesses
      nonzero_subprocess(1:maxproc_dps,1:maxproc_dps)=.FALSE.
      wme_subprocess(1:maxproc_dps,1:maxproc_dps)=0d0
      wme_calc(1:2,1:maxproc_dps)=.FALSE.
      wme_proc(1:2,1:maxproc_dps)=0d0
      strf_pdf_save(1:2,1:2,-5:5)=0d0 ! pdf*x
      strf_pdf_saveq(1:2)=.FALSE.
      IF(irun.EQ.1)THEN
         ! running of alpha_s
         CALL RunAlphas(pdfnumpdf)
      ENDIF
      wmetot=0d0
      ! wme is pdf*ME2 and symmetry factor for the two dps subprocesses
      DO iproc1=1,maxproc_dps
         IF(.NOT.contained_subprocess(iproc1))CYCLE
         CALL determine_subprocess(1,iproc1)
         icut=1
         CALL Cuts_pp_jj_aj_aa(icut)
         IF(icut.EQ.0)CYCLE
         ! determine s and t
         CALL calc_st(1)
         CALL ME2_subprocess(1,iproc1,wme1)
         IF(wme1.LE.0d0)CYCLE
         ! uses nomrla sPDF
         CALL strf_pdf_subprocess(pdfnumpdf,wsf1)
         wme1=wme1*wsf1
         IF(wme1.LE.0d0)CYCLE
         ! multiply K factor
         IF(iproc1.GE.1.AND.iproc1.LE.176)THEN
            wme1=wme1*Kfactor_jj
         ELSEIF(iproc1.GE.177.AND.iproc1.LE.206)THEN
            wme1=wme1*Kfactor_aj
         ELSEIF(iproc1.GE.207.AND.iproc1.LE.217)THEN
            wme1=wme1*Kfactor_aa
         ENDIF
         IF(wme1.LE.0d0)CYCLE
         nonzero_subprocess(iproc1,iproc1)=.TRUE.
         wme_subprocess(iproc1,iproc1)=wme1
         wmetot=wmetot+wme1
      ENDDO
      IF(wmetot.EQ.0d0)THEN
         pp_jj_aj_aa_fxn=0d0
         RETURN
      ENDIF
      pp_jj_aj_aa_fxn=wdps*wmetot
      IF(lunwei2.AND.lwmax.GT.0)THEN
         w1=pp_jj_aj_aa_fxn*wgt ! multiply VEGAS weight
         CALL unwei_procedure_pp_aajj_dps(w1,nwri,nwmax,nwarn)
      ENDIF
      nnn=nnn+1
      IF(recmax.LT.pp_jj_aj_aa_fxn)recmax=pp_jj_aj_aa_fxn
      IF(recmin.GT.pp_jj_aj_aa_fxn)recmin=pp_jj_aj_aa_fxn
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
    END FUNCTION pp_jj_aj_aa_fxn

    SUBROUTINE RunAlphas(pdfnumpdf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::pdfnumpdf
      REAL(KIND(1d0))::scale1,scale2
      REAL(KIND(1d0))::alphas00
      REAL(KIND(1d0))::scale
      REAL(KIND(1d0)),EXTERNAL::GSALPS
      CALL setscale_dps(1,scale1)
      CALL setscale_dps(2,scale2)
      IF(pdfnumpdf.EQ.921000)THEN
         ! use the running of GS09 dPDF
         CALL Combine_scale_dps(scale1,scale2,scale)
         alphas00=GSALPS(scale)
         gstrong1=DSQRT(4*pi*alphas00)
         gstrong2=gstrong1
      ELSE
         ! use the normal sPDF
         alphas00=ALPHAS(scale1)
         gstrong1=DSQRT(4*pi*alphas00)
         alphas00=ALPHAS(scale2)
         gstrong2=DSQRT(4*pi*alphas00)
      ENDIF
      RETURN
    END SUBROUTINE RunAlphas

    SUBROUTINE select_subprocess(r1,r2,r3,found)
      IMPLICIT NONE
      !REAL(KIND(1d0)),INTENT(IN)::wmetot
      REAL(KIND(1d0)),INTENT(IN)::r1,r2,r3
      REAL(KIND(1d0))::wmeacc
      LOGICAL,INTENT(OUT)::found
      LOGICAL::allowed
      INTEGER::iproc1,iproc2,iproc1_choose,iproc2_choose
      ! according to wme to select the subprocess to output
      wmeacc=0d0
      found=.FALSE.
      DO iproc1=1,maxproc_dps
         DO iproc2=iproc1,maxproc_dps
            IF(.NOT.nonzero_subprocess(iproc1,iproc2))CYCLE
            wmeacc=wmeacc+wme_subprocess(iproc1,iproc2)
            IF(r1.LE.wmeacc/wmetot)THEN
               found=.TRUE.
               iproc1_choose=iproc1
               iproc2_choose=iproc2
               CALL determine_subprocess(1,iproc1)
               CALL determine_subprocess(2,iproc2)
               CALL subprocess_dps(allowed)
               IF(.NOT.allowed)THEN
                  WRITE(*,*)"ERROR 3: Inconsistent of allowed processes"
                  STOP
               ENDIF
               ! now we have a+a+j+j
               CALL pp_aajj_dps_hadronmom
               EXIT
            ENDIF
         ENDDO
         IF(found)EXIT
      ENDDO
      IF(found)THEN
         ! choose the corresponding color flow
         CALL calc_st(1)
         CALL CF_subprocess(1,iproc1_choose,r2)
         CALL calc_st(2)
         CALL CF_subprocess(2,iproc2_choose,r3)
         ! define the icol_un and imothup_save
         ! in the new sorted order
         CALL pp_aajj_dps_colmoth
      ENDIF
      RETURN
    END SUBROUTINE select_subprocess

    SUBROUTINE select_subprocess_sps(r1,r2,found)
      IMPLICIT NONE
      REAL(KIND(1d0)),INTENT(IN)::r1,r2
      REAL(KIND(1d0))::wmeacc
      LOGICAL,INTENT(OUT)::found
      LOGICAL::allowed
      INTEGER::iproc1,iproc1_choose
      wmeacc=0d0
      found=.FALSE.
      DO iproc1=1,maxproc_dps
         IF(.NOT.nonzero_subprocess(iproc1,iproc1))CYCLE
         wmeacc=wmeacc+wme_subprocess(iproc1,iproc1)
         IF(r1.LE.wmeacc/wmetot)THEN
            found=.TRUE.
            iproc1_choose=iproc1
            CALL determine_subprocess(1,iproc1)
            EXIT
         ENDIF
      ENDDO
      IF(found)THEN
         ! choose the corresponding color flow 
         CALL calc_st(1)
         CALL CF_subprocess(1,iproc1_choose,r2)
         ! define the icol_un and imothup_save
         ! in the new sorted order         
      ENDIF
      RETURN
    END SUBROUTINE select_subprocess_sps

    SUBROUTINE calc_st(ipip)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::ipip
      REAL(KIND(1d0))::s12,t13
      REAL(KIND(1d0)),DIMENSION(4)::p12,p13
      p12(1:4)=dps_pmom(ipip,1,1:4)+dps_pmom(ipip,2,1:4)
      p13(1:4)=dps_pmom(ipip,1,1:4)-dps_pmom(ipip,3,1:4)
      s12=p12(4)**2-p12(1)**2-p12(2)**2-p12(3)**2
      t13=p13(4)**2-p13(1)**2-p13(2)**2-p13(3)**2
      s12_aj=s12
      s12_jj=s12
      s12_aa=s12
      t13_aj=t13
      t13_aa=t13
      t13_jj=t13
      RETURN
    END SUBROUTINE calc_st

    SUBROUTINE subprocess_dps(allowed)
      IMPLICIT NONE
      LOGICAL,INTENT(OUT)::allowed ! judge if it is aa+jj
      INTEGER::na,ii,ipip,nj
      iflh(1)=ifldps(1,1)
      iflh(2)=ifldps(1,2)
      iflh(3)=ifldps(2,1)
      iflh(4)=ifldps(2,2)
      na=0
      DO ipip=1,2
         DO ii=3,4
            IF(ifldps(ipip,ii).EQ.22)THEN
               na=na+1
               iflh(na+4)=ifldps(ipip,ii)
            ENDIF
         ENDDO
      ENDDO
      IF(na.NE.2)THEN
         allowed=.FALSE.
         RETURN
      ENDIF
      nj=0
      DO ipip=1,2
         DO ii=3,4
            IF(ifldps(ipip,ii).NE.22)THEN
               nj=nj+1
               iflh(nj+6)=ifldps(ipip,ii)
            ENDIF
         ENDDO
      ENDDO
      allowed=.TRUE.
      RETURN
    END SUBROUTINE subprocess_dps

    SUBROUTINE determine_subprocess(ipip,iproc)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::ipip,iproc
      INTEGER::i1,i2,f1,f2,ii,jj,kk,ll
      IF(iproc.EQ.1)THEN
         ! g g > g g
         ifldps(ipip,1)=21
         ifldps(ipip,2)=21
         ifldps(ipip,3)=21
         ifldps(ipip,4)=21
      ELSEIF(iproc.GE.2.AND.iproc.LE.6)THEN
         ! g g > q q~
         ifldps(ipip,1)=21
         ifldps(ipip,2)=21
         f1=iproc-1
         f2=-f1
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.7.AND.iproc.LE.16)THEN
         ! q q~ > g g or q~ q > g g
         ii=iproc-6
         IF(MOD(ii,2).EQ.0)THEN
            ! q~ q > g g
            i1=-(ii+1)/2
         ELSE
            ! q q~ > g g
            i1=(ii+1)/2
         ENDIF
         i2=-i1
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=21
         ifldps(ipip,4)=21
      ELSEIF(iproc.GE.17.AND.iproc.LE.36)THEN
         ii=iproc-16
         IF(MOD(ii,4).EQ.1)THEN
            ! q g > g q
            i1=(ii+3)/4
            i2=21
            f1=21
            f2=i1
         ELSEIF(MOD(ii,4).EQ.2)THEN
            ! q~ g > g q~
            i1=-(ii+3)/4
            i2=21
            f1=21
            f2=i1
         ELSEIF(MOD(ii,4).EQ.3)THEN
            ! g q > g q
            i1=21
            i2=(ii+3)/4
            f1=21
            f2=i2
         ELSE
            ! g q~ > g q~
            i1=21
            i2=-(ii+3)/4
            f1=21
            f2=i2
         ENDIF
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.37.AND.iproc.LE.76)THEN
         ! q1 q1~ > q2 q2~ or q1~ q1 > q2 q2~ with q1!=q2
         ii=iproc-36
         jj=(ii+7)/8
         kk=MOD(ii-1,8)
         kk=kk/2+1
         IF(MOD(ii,2).EQ.0)THEN
            ! q1~ q1 > q2 q2~
            i1=-jj
            i2=jj
            IF(kk.GE.jj)THEN
               f1=kk+1
               f2=-f1
            ELSE
               f1=kk
               f2=-f1
            ENDIF
         ELSE
            ! q1 q1~ > q2 q2~
            i1=jj
            i2=-jj
            IF(kk.GE.jj)THEN
               f1=kk+1
               f2=-f1
            ELSE
               f1=kk
               f2=-f1
            ENDIF
         ENDIF
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.77.AND.iproc.LE.86)THEN
         ! q1 q1~ > q1 q1~  or q1~ q1 > q1 q1~
         ii=iproc-76
         IF(MOD(ii,2).EQ.2)THEN
            ! q1~ q1 > q1 q1~
            i1=-(ii+1)/2
         ELSE
            ! q1 q1~ > q1 q1~
            i1=(ii+1)/2
         ENDIF
         i2=-i1
         f1=ABS(i1)
         f2=-f2
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.87.AND.iproc.LE.96)THEN
         ! q q > q q or q~ q~ > q~ q~
         ii=iproc-86
         IF(MOD(ii,2).EQ.0)THEN
            ! q~ q~ > q~ q~
            i1=-(ii+1)/2
         ELSE
            ! q q > q q
            i1=(ii+1)/2
         ENDIF
         i2=i1
         f1=i1
         f2=i1
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.97.AND.iproc.LE.176)THEN
         ! q1 q2 > q1 q2 or q2 q1 > q1 q2
         ! where q1 and q2 can be quark or antiquark
         ! q1 != q2
         ii=iproc-96
         jj=(ii+7)/8
         kk=MOD(ii-1,8)+1
         IF(jj.LE.4)THEN
            f1=1
            f2=jj+1
         ELSEIF(jj.GE.5.AND.jj.LE.7)THEN
            f1=2
            f2=jj-4+2
         ELSEIF(jj.GE.8.AND.jj.LE.9)THEN
            f1=3
            f2=jj-7+3
         ELSEIF(jj.EQ.10)THEN
            f1=4
            f2=jj-9+4
         ENDIF
         IF(kk.GT.4)THEN
            f1=-f1
         ENDIF
         IF(MOD(kk,4).EQ.3.OR.MOD(kk,4).EQ.0)THEN
            f2=-f2
         ENDIF
         IF(MOD(kk,2).EQ.0)THEN
            i1=f2
            i2=f1
         ELSE
            i1=f1
            i2=f2
         ENDIF
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.177.AND.iproc.LE.186)THEN
         ! q q~ > a g or q~ q > a g
         ii=iproc-176
         jj=(ii+1)/2
         IF(MOD(ii,2).EQ.0)THEN
            ! q~ q > a g
            i1=-jj
         ELSE
            ! q q~ > a g
            i1=jj
         ENDIF
         i2=-i1
         f1=22
         f2=21
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.187.AND.iproc.LE.206)THEN
         ! q g > a q, q~ g > a q~
         ! g q > a q, g q~ > a q~
         ii=iproc-186
         jj=(ii+3)/4
         IF(MOD(ii,4).EQ.1)THEN
            ! q g > a q
            i1=jj
            i2=21
            f1=22
            f2=i1
         ELSEIF(MOD(ii,4).EQ.2)THEN
            ! q~ g > a q~
            i1=-jj
            i2=21
            f1=22
            f2=i1
         ELSEIF(MOD(ii,4).EQ.3)THEN
            ! g q > a q
            i1=21
            i2=jj
            f1=22
            f2=i2
         ELSE
            ! g q~ > a q~
            i1=21
            i2=-jj
            f1=22
            f2=i2
         ENDIF
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.GE.207.AND.iproc.LE.216)THEN
         ! q q~ > a a or q~ q > a a
         ii=iproc-206
         jj=(ii+1)/2
         IF(MOD(ii,2).EQ.0)THEN
            ! q~ q > a a
            i1=-jj
         ELSE
            ! q q~ > a a
            i1=jj
         ENDIF
         i2=-i1
         f1=22
         f2=22
         ifldps(ipip,1)=i1
         ifldps(ipip,2)=i2
         ifldps(ipip,3)=f1
         ifldps(ipip,4)=f2
      ELSEIF(iproc.EQ.217)THEN
         ! g g > a a
         ifldps(ipip,1)=21
         ifldps(ipip,2)=21
         ifldps(ipip,3)=22
         ifldps(ipip,4)=22
      ELSE
         WRITE(*,*)"ERROR:Unknow process number = ",iproc
         STOP
      ENDIF
      RETURN
    END SUBROUTINE determine_subprocess

    SUBROUTINE ME2_subprocess(ipip,iproc,wme)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::ipip,iproc
      REAL(KIND(1d0)),INTENT(OUT)::wme
      INTEGER::iqqx,ii,idq,jj,kk,ll
      IF(iproc.GT.0.AND.iproc.LE.maxproc_dps)THEN
         IF(wme_calc(ipip,iproc))THEN
            ! recycling
            wme=wme_proc(ipip,iproc)
            RETURN
         ENDIF
      ENDIF
      IF(irun.EQ.1)THEN
         IF(ipip.EQ.1)THEN
            Gstrong=gstrong1
         ELSE
            Gstrong=gstrong2
         ENDIF
      ENDIF
      IF(iproc.EQ.1)THEN
         ! g g > g g
         CALL ME2_gg_gg(wme)
         ! save it
         wme_calc(ipip,iproc)=.TRUE.
         wme_proc(ipip,iproc)=wme
      ELSEIF(iproc.GE.2.AND.iproc.LE.6)THEN
         ! g g > q q~
         CALL ME2_gg_qq(wme)
         ! save it
         wme_calc(ipip,2:6)=.TRUE.
         wme_proc(ipip,2:6)=wme
      ELSEIF(iproc.GE.7.AND.iproc.LE.16)THEN
         ! q q~ > g g or q~ q > g g
         CALL ME2_qq_gg(wme)
         ! save it
         wme_calc(ipip,7:16)=.TRUE.
         wme_proc(ipip,7:16)=wme
      ELSEIF(iproc.GE.17.AND.iproc.LE.36)THEN
         ii=iproc-16
         iqqx=MOD(ii-1,4)+1
         CALL ME2_gq_gq(iqqx,wme)
         ! save it
         DO jj=17,36
            kk=jj-16
            kk=MOD(kk-1,4)+1
            IF(iqqx.EQ.kk)THEN
               wme_calc(ipip,jj)=.TRUE.
               wme_proc(ipip,jj)=wme
            ENDIF
         ENDDO
      ELSEIF(iproc.GE.37.AND.iproc.LE.76)THEN
         ! q1 q1~ > q2 q2~ or q1~ q1 > q2 q2~ with q1!=q2
         CALL ME2_q1q1_q2q2(wme)
         ! save it
         wme_calc(ipip,37:76)=.TRUE.
         wme_proc(ipip,37:76)=wme
      ELSEIF(iproc.GE.77.AND.iproc.LE.86)THEN
         ! q1 q1~ > q1 q1~  or q1~ q1 > q1 q1~
         ii=iproc-76
         iqqx=MOD(ii-1,2)+1
         CALL ME2_q1q1x_q1q1x(iqqx,wme)
         ! save it
         DO jj=77,86
            kk=jj-76
            kk=MOD(kk-1,2)+1
            IF(iqqx.EQ.kk)THEN
               wme_calc(ipip,jj)=.TRUE.
               wme_proc(ipip,jj)=wme
            ENDIF
         ENDDO
      ELSEIF(iproc.GE.87.AND.iproc.LE.96)THEN
         ! q q > q q or q~ q~ > q~ q~
         CALL ME2_q1q1_q1q1(wme)
         ! save it
         wme_calc(ipip,87:96)=.TRUE.
         wme_proc(ipip,87:96)=wme
      ELSEIF(iproc.GE.97.AND.iproc.LE.176)THEN
         ! q1 q2 > q1 q2 or q2 q1 > q1 q2
         ! where q1 and q2 can be quark or antiquark 
         ! q1 != q2
         ii=iproc-96
         iqqx=MOD(ii-1,8)+1
         CALL ME2_q1q2_q1q2(iqqx,wme)
         ! save it
         DO jj=97,176
            kk=jj-96
            kk=MOD(kk-1,8)+1
            IF(iqqx.EQ.kk)THEN
               wme_calc(ipip,jj)=.TRUE.
               wme_proc(ipip,jj)=wme
            ENDIF
         ENDDO
      ELSEIF(iproc.GE.177.AND.iproc.LE.186)THEN
         ! q q~ > a g or q~ q > a g
         ii=iproc-176
         idq=(ii+1)/2
         CALL ME2_qq_ag(idq,wme)
         ! save it
         DO jj=177,186
            kk=jj-176
            kk=(kk+1)/2
            IF(idq.EQ.kk)THEN
               wme_calc(ipip,jj)=.TRUE.
               wme_proc(ipip,jj)=wme
            ENDIF
         ENDDO
      ELSEIF(iproc.GE.187.AND.iproc.LE.206)THEN
         ! q g > a q, q~ g > a q~ 
         ! g q > a q, g q~ > a q~
         ii=iproc-186
         idq=(ii+3)/4
         iqqx=MOD(ii-1,4)+1
         CALL ME2_qg_aq(idq,iqqx,wme)
         ! save it
         DO jj=187,206
            kk=jj-186
            ll=(kk+3)/4
            kk=MOD(kk-1,4)+1
            IF(iqqx.EQ.kk.AND.idq.EQ.ll)THEN
               wme_calc(ipip,jj)=.TRUE.
               wme_proc(ipip,jj)=wme
            ENDIF
         ENDDO
      ELSEIF(iproc.GE.207.AND.iproc.LE.216)THEN
         ! q q~ > a a or q~ q > a a
         ii=iproc-206
         idq=(ii+1)/2
         CALL ME2_qq_aa(idq,wme)
         ! save it
         DO jj=207,216
            kk=jj-206
            kk=(kk+1)/2
            IF(idq.EQ.kk)THEN
               wme_calc(ipip,jj)=.TRUE.
               wme_proc(ipip,jj)=wme
            ENDIF
         ENDDO
      ELSEIF(iproc.EQ.217)THEN
         ! g g > a a
         CALL ME2_gg_aa(wme)
         ! save it
         wme_calc(ipip,iproc)=.TRUE.
         wme_proc(ipip,iproc)=wme
      ELSE
         WRITE(*,*)"ERROR 2:Unknow process number = ",iproc
         STOP
      ENDIF
      RETURN
    END SUBROUTINE ME2_subprocess

    SUBROUTINE CF_subprocess(ipip,iproc,rcol)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::ipip
      INTEGER,INTENT(IN)::iproc
      REAL(KIND(1d0)),INTENT(IN)::rcol
      INTEGER::iqqx,ii
      IF(irun.EQ.1)THEN
         IF(ipip.EQ.1)THEN
            Gstrong=gstrong1
         ELSE
            Gstrong=gstrong2
         ENDIF
      ENDIF
      IF(iproc.EQ.1)THEN
         ! g g > g g
         CALL CF_gg_gg(ipip,rcol)
      ELSEIF(iproc.GE.2.AND.iproc.LE.6)THEN
         ! g g > q q~ 
         CALL CF_gg_qq(ipip,rcol)
      ELSEIF(iproc.GE.7.AND.iproc.LE.16)THEN
         ! q q~ > g g or q~ q > g g
         ii=iproc-6
         iqqx=MOD(ii-1,2)+1
         CALL CF_qq_gg(iqqx,ipip,rcol)
      ELSEIF(iproc.GE.17.AND.iproc.LE.36)THEN
         ii=iproc-16
         iqqx=MOD(ii-1,4)+1
         CALL CF_gq_gq(iqqx,ipip,rcol)
      ELSEIF(iproc.GE.37.AND.iproc.LE.76)THEN
         ! q1 q1~ > q2 q2~ or q1~ q1 > q2 q2~ with q1!=q2
         ii=iproc-36
         iqqx=MOD(ii-1,2)+1
         CALL CF_q1q1_q2q2(iqqx,ipip)
      ELSEIF(iproc.GE.77.AND.iproc.LE.86)THEN
         ! q1 q1~ > q1 q1~  or q1~ q1 > q1 q1~
         ii=iproc-76
         iqqx=MOD(ii-1,2)+1
         CALL CF_q1q1x_q1q1x(iqqx,ipip,rcol)
      ELSEIF(iproc.GE.87.AND.iproc.LE.96)THEN
         ! q q > q q or q~ q~ > q~ q~
         ii=iproc-86
         iqqx=MOD(ii-1,2)+1
         CALL CF_q1q1_q1q1(iqqx,ipip,rcol)
      ELSEIF(iproc.GE.97.AND.iproc.LE.176)THEN
         ! q1 q2 > q1 q2 or q2 q1 > q1 q2 
         ! where q1 and q2 can be quark or antiquark
         ! q1 != q2
         ii=iproc-96
         iqqx=MOD(ii-1,8)+1
         CALL CF_q1q2_q1q2(iqqx,ipip)
      ELSEIF(iproc.GE.177.AND.iproc.LE.186)THEN
         ! q q~ > a g or q~ q > a g
         ii=iproc-176
         iqqx=MOD(ii-1,2)+1
         CALL CF_qq_ag(iqqx,ipip)
      ELSEIF(iproc.GE.187.AND.iproc.LE.206)THEN
         ! q g > a q, q~ g > a q~ 
         ! g q > a q, g q~ > a q~
         ii=iproc-186
         iqqx=MOD(ii-1,4)+1
         CALL CF_qg_aq(iqqx,ipip)
      ELSEIF(iproc.GE.207.AND.iproc.LE.216)THEN
         ! q q~ > a a or q~ q > a a
         ii=iproc-206
         iqqx=MOD(ii-1,2)+1
         CALL CF_qq_aa(iqqx,ipip)
      ELSEIF(iproc.EQ.217)THEN
         ! g g > a a
         CALL CF_gg_aa(ipip)
      ENDIF
      RETURN
    END SUBROUTINE CF_subprocess

    SUBROUTINE strf_pdf_subprocess(pdfnumpdf,wsf)
      IMPLICIT NONE
      INTEGER,INTENT(IN)::pdfnumpdf
      REAL(KIND(1d0)),INTENT(OUT)::wsf
      REAL(KIND(1d0))::scale,scale1,scale2
      REAL(KIND(1d0))::xp11,xp21,xp12,xp22
      REAL(KIND(1d0))::wsft
      INTEGER::id11,id21,id12,id22,nid21,nid22
      xp11=xpi(1,1)
      xp21=xpi(1,2)
      xp12=xpi(2,1)
      xp22=xpi(2,2)
      id11=ifldps(1,1)
      IF(id11.EQ.21)id11=0
      id21=ifldps(1,2)
      IF(id21.EQ.21)id21=0
      id12=ifldps(2,1)
      IF(id12.EQ.21)id12=0
      id22=ifldps(2,2)
      IF(id22.EQ.21)id22=0
      wsf=1d0
      IF(pdfnumpdf.EQ.921000)THEN
         ! GS 09 dPDF
         IF(strf_dpdf_saveq(1,id11,id12))THEN
            wsf=wsf*strf_dpdf_save(1,id11,id12)
         ENDIF
         IF(strf_dpdf_saveq(2,id21,id22))THEN
            wsf=wsf*strf_dpdf_save(2,id21,id22)
         ENDIF
         IF(strf_dpdf_saveq(1,id11,id12).AND.strf_dpdf_saveq(2,id21,id22))RETURN
         CALL setscale_dps(1,scale1)
         CALL setscale_dps(2,scale2)
         CALL Combine_scale_dps(scale1,scale2,scale)
         IF(COLL_TYPE.EQ.2)THEN
            ! p p~
            nid21=-id21
            nid22=-id22
         ELSE
            ! p p
            nid21=id21
            nid22=id22
         ENDIF
         ! iproton, f1 from iproton in ipip=1, f2 from iproton in ipip=2
         IF(.NOT.strf_dpdf_saveq(1,id11,id12))THEN
            strf_dpdf_saveq(1,id11,id12)=.TRUE.
            CALL GS09(xp11,xp12,scale,id11,id12,wsft)
            strf_dpdf_save(1,id11,id12)=wsft
            wsf=wsf*wsft
         ENDIF
         IF(.NOT.strf_dpdf_saveq(2,id21,id22))THEN
            strf_dpdf_saveq(2,id21,id22)=.TRUE.
            CALL GS09(xp21,xp22,scale,nid21,nid22,wsft)
            strf_dpdf_save(2,id21,id22)=wsft
            wsf=wsf*wsft
         ENDIF
         ! CALL GSALPS(Q) to run alpha_S, which is not used here
      ELSE
         CALL strf_pdf_pp_aajj_dps(1,id11,id21,wsft)
         wsf=wsf*wsft
         IF(dpsorsps.EQ.1)THEN
            ! dps process
            CALL strf_pdf_pp_aajj_dps(2,id12,id22,wsft)
            wsf=wsf*wsft
         ENDIF
      ENDIF
      RETURN
    END SUBROUTINE strf_pdf_subprocess

    SUBROUTINE strf_pdf_pp_aajj_dps(ipip,id1,id2,wsf)
      USE CTEQ6PDF
      USE Structf_PDFs
      IMPLICIT NONE
      INTEGER::ipp=1
      INTEGER::ih=1 ! ih=1 no photon PDF, ih=2, photon from proton/anti-proton, ih=3 phton from electron/positron
      REAL(KIND(1d0)),INTENT(OUT)::wsf
      INTEGER,INTENT(IN)::ipip,id1,id2
      INTEGER::init=0,icut
      LOGICAL::use_cteq6_f90=.TRUE.
      REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
      SAVE init,ipp,use_cteq6_f90
      REAL(KIND(1d0))::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
                ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
      REAL(KIND(1d0))::xpp1,xpp2
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

      CALL setscale_dps(ipip,scale)
      xpp1=xpi(ipip,1)
      xpp2=xpi(ipip,2)
      IF(strf_pdf_saveq(ipip))THEN
         wsf=strf_pdf_save(ipip,1,id1)*strf_pdf_save(ipip,2,id2)/xpp1/xpp2
         IF(init.EQ.0)init=1
         RETURN
      ENDIF

      IF(use_cteq6_f90)THEN
         glu1_ct = xpp1*Ctq6Pdf_f90(0,xpp1,scale)
         glu2_ct = xpp2*Ctq6Pdf_f90(0,xpp2,scale)
         u1_ct   = xpp1*Ctq6Pdf_f90(1,xpp1,scale)
         u2_ct   = xpp2*Ctq6Pdf_f90(1,xpp2,scale)
         d1_ct   = xpp1*Ctq6Pdf_f90(2,xpp1,scale)
         d2_ct   = xpp2*Ctq6Pdf_f90(2,xpp2,scale)
         s1_ct   = xpp1*Ctq6Pdf_f90(3,xpp1,scale)
         s2_ct   = xpp2*Ctq6Pdf_f90(3,xpp2,scale)
         c1_ct   = xpp1*Ctq6Pdf_f90(4,xpp1,scale)
         c2_ct   = xpp2*Ctq6Pdf_f90(4,xpp2,scale)
         b1_ct   = xpp1*Ctq6Pdf_f90(5,xpp1,scale)
         b2_ct   = xpp2*Ctq6Pdf_f90(5,xpp2,scale)
         ub1_ct  = xpp1*Ctq6Pdf_f90(-1,xpp1,scale)
         ub2_ct  = xpp2*Ctq6Pdf_f90(-1,xpp2,scale)
         db1_ct  = xpp1*Ctq6Pdf_f90(-2,xpp1,scale)
         db2_ct  = xpp2*Ctq6Pdf_f90(-2,xpp2,scale)
         sb1_ct  = xpp1*Ctq6Pdf_f90(-3,xpp1,scale)
         sb2_ct  = xpp2*Ctq6Pdf_f90(-3,xpp2,scale)
         cb1_ct  = xpp1*Ctq6Pdf_f90(-4,xpp1,scale)
         cb2_ct  = xpp2*Ctq6Pdf_f90(-4,xpp2,scale)
         bb1_ct  = xpp1*Ctq6Pdf_f90(-5,xpp1,scale)
         bb2_ct  = xpp2*Ctq6Pdf_f90(-5,xpp2,scale)
      ELSE
         CALL fdist(ih,xpp1,scale,pdflist(-7:7))
         glu1_ct = xpp1*pdflist(0)
         u1_ct   = xpp1*pdflist(2)
         d1_ct   = xpp1*pdflist(1)
         s1_ct   = xpp1*pdflist(3)
         c1_ct   = xpp1*pdflist(4)
         b1_ct   = xpp1*pdflist(5)
         ub1_ct  = xpp1*pdflist(-2)
         db1_ct  = xpp1*pdflist(-1)
         sb1_ct  = xpp1*pdflist(-3)
         cb1_ct  = xpp1*pdflist(-4)
         bb1_ct  = xpp1*pdflist(-5)
         CALL fdist(ih,xpp2,scale,pdflist(-7:7))
         glu2_ct = xpp2*pdflist(0)
         u2_ct   = xpp2*pdflist(2)
         d2_ct   = xpp2*pdflist(1)
         s2_ct   = xpp2*pdflist(3)
         c2_ct   = xpp2*pdflist(4)
         b2_ct   = xpp2*pdflist(5)
         ub2_ct  = xpp2*pdflist(-2)
         db2_ct  = xpp2*pdflist(-1)
         sb2_ct  = xpp2*pdflist(-3)
         cb2_ct  = xpp2*pdflist(-4)
         bb2_ct  = xpp2*pdflist(-5)
      ENDIF
      strf_pdf_saveq(ipip)=.TRUE.
      strf_pdf_save(ipip,1,0)=glu1_ct
      strf_pdf_save(ipip,1,1)=d1_ct
      strf_pdf_save(ipip,1,2)=u1_ct
      strf_pdf_save(ipip,1,3)=s1_ct
      strf_pdf_save(ipip,1,4)=c1_ct
      strf_pdf_save(ipip,1,5)=b1_ct
      strf_pdf_save(ipip,1,-1)=db1_ct
      strf_pdf_save(ipip,1,-2)=ub1_ct
      strf_pdf_save(ipip,1,-3)=sb1_ct
      strf_pdf_save(ipip,1,-4)=cb1_ct
      strf_pdf_save(ipip,1,-5)=bb1_ct
      IF(ipp.EQ.2)THEN
         ! p p
         strf_pdf_save(ipip,2,0)=glu2_ct
         strf_pdf_save(ipip,2,1)=d2_ct
         strf_pdf_save(ipip,2,2)=u2_ct
         strf_pdf_save(ipip,2,3)=s2_ct
         strf_pdf_save(ipip,2,4)=c2_ct
         strf_pdf_save(ipip,2,5)=b2_ct
         strf_pdf_save(ipip,2,-1)=db2_ct
         strf_pdf_save(ipip,2,-2)=ub2_ct
         strf_pdf_save(ipip,2,-3)=sb2_ct
         strf_pdf_save(ipip,2,-4)=cb2_ct
         strf_pdf_save(ipip,2,-5)=bb2_ct
      ELSE
         ! p p~
         strf_pdf_save(ipip,2,0)=glu2_ct
         strf_pdf_save(ipip,2,-1)=d2_ct
         strf_pdf_save(ipip,2,-2)=u2_ct
         strf_pdf_save(ipip,2,-3)=s2_ct
         strf_pdf_save(ipip,2,-4)=c2_ct
         strf_pdf_save(ipip,2,-5)=b2_ct
         strf_pdf_save(ipip,2,1)=db2_ct
         strf_pdf_save(ipip,2,2)=ub2_ct
         strf_pdf_save(ipip,2,3)=sb2_ct
         strf_pdf_save(ipip,2,4)=cb2_ct
         strf_pdf_save(ipip,2,5)=bb2_ct
      ENDIF
      wsf=strf_pdf_save(ipip,1,id1)*strf_pdf_save(ipip,2,id2)/xpp1/xpp2
      IF(init.EQ.0)init=1
      RETURN
    END SUBROUTINE strf_pdf_pp_aajj_dps

!    SUBROUTINE readcuts_pp_aajj_dps
!      IMPLICIT NONE
!      CHARACTER(len=24)::file
!      LOGICAL::lexist
!      INTEGER::iounit,flag=0,i,i1,j,j1
!      REAL(KIND(1d0))::ptq,ptg,etaq,ycq,ycqlow,etag,ycg,ycglow
!      REAL(KIND(1d0))::cutoff,xFcq,xFcqlow,xFcg,xFcglow
!      REAL(KIND(1d0))::gbeamq,drqq,gqq
!      INTEGER::nnhad
!      ! open default input file
!      INQUIRE(FILE=TRIM(input_dir)//"default.inp",EXIST=lexist)
!      IF(.NOT.lexist)THEN
!         PRINT *,"Warning: the file default.inp does not exist ! STOP !"
!         STOP
!      ENDIF
!      INQUIRE(FILE=TRIM(input_dir)//"default.inp",OPENED=lexist)
!      IF(lexist)THEN
!         INQUIRE(FILE=TRIM(input_dir)//"default.inp",NUMBER=iounit)
!         IF(iounit.NE.udefault)THEN
!            PRINT *,"WARNING: the default.inp has been linked with another unit ! Close and reopen !"
!            CLOSE(UNIT=iounit)
!            OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
!         ENDIF
!      ELSE
!         OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
!      ENDIF
!      ! open user's input file
!      IF(TRIM(Input_File)/="default.inp")THEN
!         INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
!         IF(.NOT.lexist)THEN
!            PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
!            STOP
!         ENDIF
!         INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
!         IF(lexist)THEN
!            INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
!            IF(iounit.NE.uinput)THEN
!               PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
!               CLOSE(UNIT=iounit)
!               OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
!            ENDIF
!         ELSE
!            OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
!         ENDIF
!      ELSE
!         flag=1
!      ENDIF
!      cutoff=readvalue_r("cutoffp",flag)
!      PRINT *,"WARNING CUTOFF SET:",cutoff
!      ptc(3:20)=cutoff
!      drc(3:20,3:20)=0d0
!      etac(3:20)=20
!      yycut(3:20)=1d9
!      IF(absrap)THEN
!         yycutlow(3:20)=0d0
!      ELSE
!         yycutlow(3:20)=-1d9
!      ENDIF
!      xFcut(3:20)=1d0
!      xFcutlow(3:20)=1d0
!      xFcutflag=.FALSE.
!      ec(3:20)=cutoff
!      c1(3:20)=1d0
!      c2(3:20)=1d0
!      cc(3:20,3:20)=1d0
!      gmas(3:20,3:20)=cutoff
!      gbeammass(1:2,3:20)=0d0
!      y1cup=30d0
!      y1clow=0d0
!      y1cup=readvalue_r("maxy1c",flag)
!      y1clow=readvalue_r("miny1c",flag)
!      ! minimum quark pt (5-flavor quarks and gluon)
!      ptq=readvalue_r("minptq",flag)
!      ! minimum photon pt
!      ptg=readvalue_r("minptp",flag)
!      ! maximum rapidity quark 
!      etaq=readvalue_r("maxrapq",flag)
!      ! maximum y rapidity quark
!      ycq=readvalue_r("maxyrapq",flag)
!      ! minimum y rapidity quark
!      ycqlow=readvalue_r("minyrapq",flag)
!      ! maximum Feynman parameter xF
!      xFcq=readvalue_r("maxxFq",flag)
!      ! minimum Feynman parameter xF
!     xFcqlow=readvalue_r("minxFq",flag)
!      ! maximum rapidity photon
!      etag=readvalue_r("maxrapp",flag)
!      ! maximum y rapidity photon
!      ycg=readvalue_r("maxyrapp",flag)
!      ! minimum y rapidity photon
!      ycglow=readvalue_r("minyrapp",flag)
!     ! maximum Feynman parameter xF
!      xFcg=readvalue_r("maxxFp",flag)
!      ! minimum Feynman parameter xF
!      xFcglow=readvalue_r("minxFp",flag)
!      ! minimum mass quark with quark
!      gqq=readvalue_r("minmqqp",flag)
!      ! minimum mass u,d,s quarks and gluon  with partonic beam
!      gbeamq=readvalue_r("minmqbeam",flag)
!      ! minimum dr quark with quark
!      drqq=readvalue_r("mindrqq",flag)
!      CLOSE(UNIT=udefault)
!      CLOSE(UNIT=uinput)
!      IF(dpsorsps.EQ.1)THEN
!         nnhad=nhad
!     ELSE
!         nnhad=2*nhad
!      ENDIF
!      DO i=5,nnhad
!         ! the first two final states are photon
!         ! while the last two final states are jets
!         IF(i.GT.6)THEN
!            ! jets
!            ptc(i)=ptq
!            etac(i)=etaq
!            yycut(i)=ycq
!            yycutlow(i)=ycqlow
!            xFcut(i)=xFcq
!            xFcutlow(i)=xFcqlow
!         ELSE
!            ! photon
!            ptc(i)=ptg
!            etac(i)=etag
!            yycut(i)=ycg
!            yycutlow(i)=ycglow
!            xFcut(i)=xFcg
!            xFcutlow(i)=xFcglow
!         ENDIF
!      ENDDO
!      DO i=5,nnhad
!         DO j=i+1,nnhad
!            IF(i.EQ.7.AND.j.EQ.8)THEN
!               ! dijets
!               drc(i,j)=drqq
!               gmas(i,j)=MAX(gqq,gmas(i,j))
!            ENDIF
!         ENDDO
!      ENDDO
!      DO i=5,nnhad-1
!         DO j=i+1,nnhad
!            gmas(i,j)=MAX(gmas(i,j),DSQRT(2*ptc(i)*ptc(j)*(1-COS(drc(i,j)))))
!         ENDDO
!      ENDDO
!      DO i=6,nnhad
!         DO j=5,i-1
!            drc(i,j)=drc(j,i)
!            gmas(i,j)=gmas(j,i)
!         ENDDO
!      ENDDO
!      WRITE(*,*)'---------------------------------------------------'
!      WRITE(*,*)'    the cuts for p p > diphoton+dijet + X '
!      WRITE(*,*)'    with double parton scattering (DPS) '
!      DO i=5,nnhad
!         WRITE(*,*)'pt     of  ',i,'   particle   ',ptc(i)
!         WRITE(*,*)'energy of  ',i,'   particle   ',ec(i)
!         WRITE(*,*)'rapidity of  ',i,'   particle   ',etac(i)
!         WRITE(*,*)'max y rapidity of ',i,'   particle   ',yycut(i)
!         WRITE(*,*)'min y rapidity of ',i,'   particle   ',yycutlow(i)
!         WRITE(*,*)'max Feynman parameter xF of ',i,' particle ',xFcut(i)
!         IF(xFcut(i).LT.1d0)xFcutflag=.TRUE.
!         WRITE(*,*)'min Feynman parameter xF of ',i,' particle ',xFcutlow(i)
!         IF(xFcutlow(i).GT.-1d0)xFcutflag=.TRUE.
!      ENDDO
!      WRITE(*,*)'The maxrapidity of the first particle',y1cup
!      WRITE(*,*)'The minrapidity of the first particle',y1clow
!      DO i=5,nnhad-1
!         DO j=i+1,nnhad
!            WRITE(*,*)'DR     ',i,'  with  ',j,drc(i,j)
!            WRITE(*,*)'mass of ',i,'  with  ',j,gmas(i,j)
!         ENDDO
!      ENDDO
!      WRITE(*,*)'---------------------------------------------------'
!    END SUBROUTINE readcuts_pp_aajj_dps

!    SUBROUTINE Cuts_pp_aajj_dps(icut)
!      IMPLICIT NONE
!      INTEGER,INTENT(OUT)::icut
!      INTEGER::l,l1,l2,flag
!      REAL(KIND(1d0))::s,d1,d2,dr,pt,eta,aaa,bbb,ptcut
!      REAL(KIND(1d0)),DIMENSION(4)::ponia,pboo2
!      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
!      REAL(KIND(1d0))::e,q
!      !INCLUDE "pp_psipsi_dps.inc"
!      icut=0
!      ! invariant mass cuts
!      DO l1=5,nhad-1
!         DO l2=l1+1,nhad
!            s=2*scalar_product(dps_hadron_pmom(l1,1:4),dps_hadron_pmom(l2,1:4))&
!                 +scalar_product(dps_hadron_pmom(l1,1:4),dps_hadron_pmom(l1,1:4))&
!                 +scalar_product(dps_hadron_pmom(l2,1:4),dps_hadron_pmom(l2,1:4))
!            IF(s.LT.gmas(l1,l2)**2)RETURN
!         ENDDO
!      ENDDO
!      flag=0
!      DO l=5,nhad
!         pt=DSQRT(dps_hadron_pmom(l,1)**2+dps_hadron_pmom(l,2)**2)
!         IF(pt.LT.ptc(l))THEN
!            flag=1
!            EXIT
!         ENDIF
!         eta=prapidity(dps_hadron_pmom(l,1:4))
!         IF(ABS(eta).GT.etac(l))THEN
!           flag=1
!            EXIT
!         ENDIF
!         eta=rapidity(dps_hadron_pmom(l,1:4))
!         IF(absrap)eta=ABS(eta)
!         IF(eta.GT.yycut(l).OR.eta.LT.yycutlow(l))THEN
!            flag=1
!            EXIT
!         ENDIF
!      ENDDO
!      IF(flag.EQ.0)THEN
!         DO l1=5,nhad-1
!            DO l2=l1+1,nhad
!               d1=prapidity(dps_hadron_pmom(l1,1:4))-prapidity(dps_hadron_pmom(l2,1:4))
!               d2=ph4(dps_hadron_pmom(l1,1),dps_hadron_pmom(l1,2),dps_hadron_pmom(l1,3))&
!                    -ph4(dps_hadron_pmom(l2,1),dps_hadron_pmom(l2,2),dps_hadron_pmom(l2,3))
!              d2=MIN(DABS(d2),2*pi-DABS(d2))
!               IF(d2/pi.GT.1d0)WRITE(*,*)d2/pi
!               dr=SQRT(d1**2+d2**2)
!               IF(dr.LT.drc(l1,l2))THEN
!                  flag=1
!                  EXIT
!               ENDIF
!            ENDDO
!            IF(flag.EQ.1)EXIT
!         ENDDO
!         IF(flag.EQ.0)icut=1
!      ENDIF
!      ! special cutoffs in the literature
!      IF(xFcutflag.AND.icut.EQ.1.AND.flag.EQ.0)THEN
!         q=ehat
!         e=q/SQRT(xp1*xp2)
!         IF(.NOT.labeqcoll)THEN
!            IF(.NOT.fixtarget)THEN
!               pboo2(4)=(ebeam(1)+ebeam(2))
!               pboo2(3)=-(ebeam(1)-ebeam(2))
!            ELSE
!               IF(.NOT.fixtargetrev)THEN
!                  pboo2(4)=(ebeam(1)+ebeam(2))
!                  pboo2(3)=-ebeam(1)
!               ELSE
!                  pboo2(4)=(ebeam(1)+ebeam(2))
!                  pboo2(3)=ebeam(2)
!               ENDIF
!            ENDIF
!            pboo2(1:2)=0
!            ! boost from the lab frame to the collision frame 
!            DO l=5,nhad
!               CALL boostl(e,pboo2,dps_hadron_pmom(l,1:4))
!            ENDDO
!         ENDIF
!         DO l=5,nhad
!            eta=xFeynman(dps_hadron_pmom(l,1:4),e)
!            IF(eta.GT.xFcut(l).OR.eta.LT.xFcutlow(l))THEN
!               icut=0
!               EXIT
!            ENDIF
!         ENDDO
!         IF(.NOT.labeqcoll)THEN
!            pboo2(3)=-pboo2(3)
!            ! boost back to the lab frame 
!            DO l=5,nhad
!               CALL boostl(e,pboo2,dps_hadron_pmom(l,1:4))
!            ENDDO
!         ENDIF
!      ENDIF
!    END SUBROUTINE Cuts_pp_aajj_dps

!    SUBROUTINE Cuts_pp_jj_aj_aa(icut)
!      IMPLICIT NONE
!      INTEGER,INTENT(OUT)::icut
!      INTEGER::l,l1,l2,flag,iil,iil1,iil2
!      REAL(KIND(1d0))::s,d1,d2,dr,pt,eta,aaa,bbb,ptcut
!      REAL(KIND(1d0)),DIMENSION(4)::ponia,pboo2
!      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
!      REAL(KIND(1d0))::e,q
!      icut=0
!      ! invariant mass cuts
!      DO l1=3,nhad-1
!         DO l2=l1+1,nhad
!            s=2*scalar_product(dps_pmom(1,l1,1:4),dps_pmom(1,l2,1:4))&
!                 +scalar_product(dps_pmom(1,l1,1:4),dps_pmom(1,l1,1:4))&
!                 +scalar_product(dps_pmom(1,l2,1:4),dps_pmom(1,l2,1:4))
!            IF((ifldps(1,l1).EQ.22.AND.ifldps(1,l2).NE.22).OR.&
!                 (ifldps(1,l1).NE.22.AND.ifldps(1,l2).EQ.22))THEN
!               ! a j
!               IF(s.LT.gmas(5,7)**2)RETURN
!            ELSEIF(ifldps(1,l1).EQ.22.AND.ifldps(1,l2).EQ.22)THEN
!               ! a a
!               IF(s.LT.gmas(5,6)**2)RETURN
!            ELSE
!               ! j j
!               IF(s.LT.gmas(7,8)**2)RETURN
!            ENDIF
!         ENDDO
!      ENDDO
!      flag=0
!      DO l=3,nhad
!         IF(ifldps(1,l).EQ.22)THEN
!            ! a
!            iil=5
!         ELSE
!            ! j
!            iil=7
!         ENDIF
!         pt=SQRT(dps_pmom(1,l,1)**2+dps_pmom(1,l,2)**2)
!         IF(pt.LT.ptc(iil))THEN
!            flag=1
!            EXIT
!         ENDIF
!         eta=prapidity(dps_pmom(1,l,1:4))
!         IF(ABS(eta).GT.etac(iil))THEN
!            flag=1
!            EXIT
!         ENDIF
!         eta=rapidity(dps_pmom(1,l,1:4))
!         IF(absrap)eta=ABS(eta)
!         IF(eta.GT.yycut(iil).OR.eta.LT.yycutlow(iil))THEN
!            flag=1
!            EXIT
!         ENDIF            
!      ENDDO
!      IF(flag.EQ.0)THEN
!         DO l1=3,nhad-1
!            IF(ifldps(1,l1).EQ.22)THEN
!               ! a 
!               iil1=5
!            ELSE
!               ! j
!               iil1=7
!            ENDIF
!            DO l2=l1+1,nhad
!               IF(ifldps(1,l2).EQ.22)THEN
!                  ! a2
!                  iil2=5
!               ELSE
!                  ! j
!                  iil2=7
!               ENDIF
!               d1=prapidity(dps_pmom(1,l1,1:4))-prapidity(dps_pmom(1,l2,1:4))
!               d2=ph4(dps_pmom(1,l1,1),dps_pmom(1,l1,2),dps_pmom(1,l1,3))&
!                    -ph4(dps_pmom(1,l2,1),dps_pmom(1,l2,2),dps_pmom(1,l2,3))
!               d2=MIN(DABS(d2),2*pi-DABS(d2))
!               IF(d2/pi.GT.1d0)WRITE(*,*)d2/pi
!               dr=SQRT(d1**2+d2**2)
!               IF(dr.LT.drc(iil1,iil2))THEN
!                  flag=1
!                  EXIT
!               ENDIF
!            ENDDO
!            IF(flag.EQ.1)EXIT
!         ENDDO
!         IF(flag.EQ.0)icut=1
!      ENDIF
!      ! special cutoffs in the literature 
!      IF(xFcutflag.AND.icut.EQ.1.AND.flag.EQ.0)THEN
!         q=ehat
!         e=q/SQRT(xp1*xp2)
!         IF(.NOT.labeqcoll)THEN
!            IF(.NOT.fixtarget)THEN
!               pboo2(4)=(ebeam(1)+ebeam(2))
!               pboo2(3)=-(ebeam(1)-ebeam(2))
!            ELSE
!               IF(.NOT.fixtargetrev)THEN
!                  pboo2(4)=(ebeam(1)+ebeam(2))
!                  pboo2(3)=-ebeam(1)
!               ELSE
!                  pboo2(4)=(ebeam(1)+ebeam(2))
!                  pboo2(3)=ebeam(2)
!               ENDIF
!            ENDIF
!            pboo2(1:2)=0
!            ! boost from the lab frame to the collision frame
!            DO l=3,nhad
!               CALL boostl(e,pboo2,dps_pmom(1,l,1:4))
!            ENDDO
!         ENDIF
!         DO l=3,nhad
!            IF(ifldps(1,l).EQ.22)THEN
!               ! a
!               iil=5
!            ELSE
!               ! j
!               iil=7
!            ENDIF
!            eta=xFeynman(dps_pmom(1,l,1:4),e)
!            IF(eta.GT.xFcut(iil).OR.eta.LT.xFcutlow(iil))THEN
!               icut=0
!               EXIT
!            ENDIF
!         ENDDO
!         IF(.NOT.labeqcoll)THEN
!            pboo2(3)=-pboo2(3)
!            ! boost back to the lab frame
!            DO l=3,nhad
!               CALL boostl(e,pboo2,dps_pmom(1,l,1:4))
!            ENDDO
!         ENDIF
!      ENDIF
!    END SUBROUTINE Cuts_pp_jj_aj_aa

    SUBROUTINE pp_aajj_dps_colmoth
      IMPLICIT NONE
      INTEGER::ipip,ii,na,nj
      icol_un(1,1:2)=icolun_dps(1,1,1:2)
      icol_un(2,1:2)=icolun_dps(1,2,1:2)
      icol_un(3,1:2)=icolun_dps(2,1,1:2)
      icol_un(4,1:2)=icolun_dps(2,2,1:2)
      imothup_save(1:4,1:2)=0
      na=0
      DO ipip=1,2
         DO ii=3,4
            IF(ifldps(ipip,ii).EQ.22)THEN
               na=na+1
               icol_un(na+4,1:2)=icolun_dps(ipip,ii,1:2)
               imothup_save(na+4,1)=ipip*2-1
               imothup_save(na+4,2)=ipip*2
            ENDIF
         ENDDO
      ENDDO
      nj=0
      DO ipip=1,2
         DO ii=3,4
            IF(ifldps(ipip,ii).NE.22)THEN
               nj=nj+1
               icol_un(nj+6,1:2)=icolun_dps(ipip,ii,1:2)
               imothup_save(nj+6,1)=ipip*2-1
               imothup_save(nj+6,2)=ipip*2
            ENDIF
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE pp_aajj_dps_colmoth

    SUBROUTINE pp_aajj_dps_hadronmom
      IMPLICIT NONE
      INTEGER::ipip,ii,na,nj
      dps_hadron_pmom(1,1:4)=dps_pmom(1,1,1:4)
      dps_hadron_pmom(1,5)=0d0
      dps_hadron_pmom(2,1:4)=dps_pmom(1,2,1:4)
      dps_hadron_pmom(2,5)=0d0
      dps_hadron_pmom(3,1:4)=dps_pmom(2,1,1:4)
      dps_hadron_pmom(3,5)=0d0
      dps_hadron_pmom(4,1:4)=dps_pmom(2,2,1:4)
      dps_hadron_pmom(4,5)=0d0
      na=0
      DO ipip=1,2
         DO ii=3,4
            IF(ifldps(ipip,ii).EQ.22)THEN
               na=na+1
               dps_hadron_pmom(na+4,1:4)=dps_pmom(ipip,ii,1:4)
               dps_hadron_pmom(na+4,5)=0d0
            ENDIF
         ENDDO
      ENDDO
      nj=0
      DO ipip=1,2
         DO ii=3,4
            IF(ifldps(ipip,ii).NE.22)THEN
               nj=nj+1
               dps_hadron_pmom(nj+6,1:4)=dps_pmom(ipip,ii,1:4)
               dps_hadron_pmom(nj+6,5)=0d0
            ENDIF
         ENDDO
      ENDDO
      RETURN
    END SUBROUTINE pp_aajj_dps_hadronmom

    SUBROUTINE unwei_procedure_pp_aajj_dps(w1,nwri,nwmax,nwarn)
      IMPLICIT NONE
      REAL(KIND(1d0)),INTENT(IN)::w1
      INTEGER,INTENT(INOUT)::nwri,nwmax,nwarn
      REAL(KIND(1d0)),DIMENSION(1)::ranr
      REAL(KIND=DBL)::vtime,vspin,scalup,xwgtup,px,py,pz,p0,pmass,umax,umax1
      INTEGER::init=0,i
      INTEGER::idup,idprup,istup,imothup1,imothup2,icol1,icol2
      LOGICAL::llwri
      REAL(KIND(1d0))::scale1,scale2,r1,r2,r3
      LOGICAL::found
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
            CALL RANDA(1,ranr)
            r1=ranr(1)
            CALL RANDA(1,ranr)
            r2=ranr(1)
            IF(dpsorsps.EQ.1)THEN
               ! DPS
               CALL RANDA(1,ranr)
               r3=ranr(1)
               CALL select_subprocess(r1,r2,r3,found)
            ELSE
               ! SPS
               CALL select_subprocess_sps(r1,r2,found)
            ENDIF
            IF(found)THEN
               nwri=nwri+1
               idprup=83 ! id for the process 
               xwgtup=1 !w1*10**3 !1
               CALL setscale_dps(1,scale1)
               IF(dpsorsps.EQ.1)THEN
                  ! DPS
                  CALL setscale_dps(2,scale2)
                  CALL Combine_scale_dps(scale1,scale2,scalup)
               ELSE
                  ! SPS
                  scalup=scale1
               ENDIF
               WRITE(nunit3)nhad,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
               DO i=1,nhad
                  IF(dpsorsps.EQ.1)THEN
                     ! DPS
                     idup=iflh(i)
                     IF(i.LE.4)THEN
                        istup=-1
                     ELSE
                        istup=1
                     ENDIF
                     imothup1=imothup_save(i,1)
                     imothup2=imothup_save(i,2)
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
                  ELSE
                     ! SPS
                     idup=ifldps(1,i)
                     IF(i.LE.2)THEN
                        istup=-1
                        imothup1=0
                        imothup2=0
                     ELSE
                        istup=1
                        imothup1=1
                        imothup2=2
                     ENDIF
                     IF(i.LE.2)THEN
                        icol1=icolun_dps(1,i,1)+100
                        IF(icol1.EQ.100)icol1=0
                        icol2=icolun_dps(1,i,2)+100
                        IF(icol2.EQ.100)icol2=0
                     ELSE
                        icol1=icolun_dps(1,i,2)+100
                        IF(icol1.EQ.100)icol1=0
                        icol2=icolun_dps(1,i,1)+100
                        IF(icol2.EQ.100)icol2=0
                     ENDIF
                     px=dps_pmom(1,i,1)
                     py=dps_pmom(1,i,2)
                     pz=dps_pmom(1,i,3)
                     p0=dps_pmom(1,i,4)
                     pmass=0d0
                  ENDIF
                  vtime=0
                  vspin=9

                  WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2&
                       ,px,py,pz,p0,pmass,vtime,vspin
               ENDDO
            ENDIF
         ENDIF
      ENDIF
      Nevents=nwri
    END SUBROUTINE unwei_procedure_pp_aajj_dps

    SUBROUTINE Generate_lhe_pp_aajj_dps(n1,nevent,icase)
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
      OPEN(nunit3,FILE=TRIM(tmp_dir)//'even_pp_aajj_dps.out',FORM='unformatted')
      nunit4=31
      CLOSE(nunit4)
      OPEN(nunit4,FILE=TRIM(tmp_dir)//'sample_pp_aajj_dps.init')
      nunit5=32
      CLOSE(nunit5)
      OPEN(nunit5,FILE=TRIM(output_dir)//'sample_pp_aajj_dps.lhe')
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
    END SUBROUTINE Generate_lhe_pp_aajj_dps
  END MODULE pp_aajj_DPS
