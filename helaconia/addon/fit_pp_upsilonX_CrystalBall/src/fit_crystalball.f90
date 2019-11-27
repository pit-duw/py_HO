PROGRAM fit_crystalball
  USE Constants
  USE fit_cb_global
  USE Helac_Global
  IMPLICIT NONE
  EXTERNAL FCN_FIT
  REAL(KIND(1d0)),EXTERNAL::CRYSTALBALL
  INTEGER,DIMENSION(4)::NPRM=(/1,2,3,4/)
  REAL(KIND(1d0)),DIMENSION(4)::VSTRT=(/0.5d0,0.35d0,13.5d0,2d0/),&
       STP=(/0.1d0,0.01d0,0d0,0d0/),ARGLIS
  REAL(KIND(1d0)),PARAMETER::ZERO=0d0,ONE=1d0
  REAL(KIND(1d0)),DIMENSION(4)::lower,upper
  CHARACTER*6,DIMENSION(4)::PNAM=(/'kapa  ','lambda','<pt>  ','n     '/)
  INTEGER::i,IERFLG
  LOGICAL::lexist
  WRITE(*,*)'                     THE BEGINNING OF HELAC-Onia'
  WRITE(*,*)'    AddOn Process: '
  WRITE(*,*)'    Fit Crystal Ball Function for p p > Upsilon + X '
  WRITE(*,*)'======================================================================='
  WRITE(*,*)'======================================================================='
  WRITE(*,*)' '
  OPEN(UNIT=1216,FILE="./input/state.inp")
  READ(1216,*)fit_istate
  CLOSE(UNIT=1216)
  IF(fit_istate.NE.1.AND.fit_istate.NE.2.AND.fit_istate.NE.3)THEN
     WRITE(*,*)"ERROR:Unknown the state = ",fit_istate
     WRITE(*,*)"INFO:Please set 1(Y(1S)) or 2(Y(2S)) or 3(Y(3S)) in state.inp"
     STOP
  ENDIF
  OPEN(UNIT=1216,FILE="./input/includeqq.inp")
  READ(1216,*)fit_includeqq
  CLOSE(UNIT=1216)
  OPEN(UNIT=1216,FILE="./input/fit_param_card.inp")
  DO i=1,4
     READ(1216,*)VSTRT(i),STP(i),lower(i),upper(i)
  ENDDO
  CLOSE(UNIT=1216)
  ! Initialize Minuit,define I/O unit numbers
  OPEN(UNIT=1216,FILE="./output/fit_crystalball.out")
  CALL MNINIT(5,1216,7)
  ! Define parameters, set initial values
  DO i=1,4
     ! don't take the limit if lower=upper=0d0
     CALL MNPARM(NPRM(i),PNAM(i),VSTRT(i),STP(i),lower(i),upper(i),IERFLG)
     IF(IERFLG.NE.0)THEN
        WRITE(1216,'(A,I2)') ' UNABLE TO DEFINE PARAMETER NO.',I
        STOP
     ENDIF
  ENDDO

  CALL ReadElem_logic('lhapdf',fit_uselhapdf)
  INQUIRE(FILE=TRIM(input_dir)//"paths/lhapdfpath",EXIST=lexist)
  IF(lexist)THEN
     OPEN(UNIT=30333,FILE=TRIM(input_dir)//"paths/lhapdfpath")
     READ(30333,'(A)')LHAPath
     CLOSE(UNIT=30333)
     i=LEN_TRIM(LHAPath)
     LHAPath=LHAPath(1:i-13)
  ENDIF
  fit_uselhapdf=fit_uselhapdf.AND.lexist
  uselhapdf=fit_uselhapdf
  CALL ReadElem_integer('pdf',fit_iPDFSUP1)
  iPDFSUP1=fit_iPDFSUP1

  CALL MNSETI('fit crystal ball function')
  ARGLIS(1)=1.
  CALL MNEXCM(FCN_FIT,'CALL FCN', ARGLIS, 1, IERFLG, CRYSTALBALL)
  CALL MNEXCM(FCN_FIT,'MIGRAD', ARGLIS, 0, IERFLG, CRYSTALBALL)
  CALL MNEXCM(FCN_FIT,'MINOS', ARGLIS,0,IERFLG, CRYSTALBALL)
  ARGLIS(1)=3.
  CALL MNEXCM(FCN_FIT,'CALL FCN',ARGLIS,1,IERFLG, CRYSTALBALL)
  CALL MNEXCM(FCN_FIT, 'STOP', 0,0,IERFLG, CRYSTALBALL)
  CLOSE(UNIT=1216)
END PROGRAM fit_crystalball

SUBROUTINE FCN_FIT(NPAR,GIN,F,X,IFLAG,FUTIL)
  USE fit_cb_global
  USE plot_fit_cb
  IMPLICIT NONE
  REAL(KIND(1d0)),DIMENSION(*)::X,GIN
  !INTEGER,PARAMETER::MXBIN=2000,MXSET=20,MXFILE=20 ! will be defined in plot_fit_cb
  REAL(KIND(1d0)),DIMENSION(MXBIN,2)::XBIN,EXVAL ! XBIN(1),XBIN(2) is Ptbin, EXVAL(1) is value and EXVAL(2) is err
  INTEGER,DIMENSION(MXSET)::NBIN_SET ! the number of bins in each experimental data set
  INTEGER,DIMENSION(MXSET)::ISTART_SET
  INTEGER,DIMENSION(MXSET)::ITYPE_SET
  REAL(KIND(1d0)),DIMENSION(MXSET)::SQRTS_SET
  REAL(KIND(1d0)),DIMENSION(MXSET,2)::Y_SET ! the y interval
  SAVE XBIN,EXVAL,NBIN_SET,SQRTS_SET,Y_SET,ISTART_SET,ITYPE_SET
  INTEGER::NPAR,IFLAG,ISET,ii,jj,nfile,istart,itype
  SAVE ISET
  REAL(KIND(1d0))::F,THVAL,XX,sqrtS,ptup,ptlow,yup,ylow
  REAL(KIND(1d0))::kapa,lam,nn,ptavg,CHI,CHISQ
  REAL(KIND(1d0)),EXTERNAL::FUTIL
  CHARACTER*100::stringin
  CHARACTER*500,DIMENSION(MXFILE)::stringfile
  INTEGER::IOstatus
  LOGICAL::lexist
  CHARACTER*15::data_filename
  IF(IFLAG.EQ.1)THEN
     IF(fit_istate.EQ.1)THEN
        data_filename="data_list_1.inp"
        mpsi=mY1S
     ELSEIF(fit_istate.EQ.2)THEN
        data_filename="data_list_2.inp"
        mpsi=mY2S
     ELSEIF(fit_istate.EQ.3)THEN
        data_filename="data_list_3.inp"
        mpsi=mY3S
     ELSE
        WRITE(*,*)"ERROR: wrong value in state.inp"
        STOP
     ENDIF
     ! read in the fitted data
     INQUIRE(FILE="./input/"//TRIM(data_filename),EXIST=lexist)
     IF(.NOT.lexist)THEN
        WRITE(*,*)"ERROR:Cannot find the file "//TRIM(data_filename)
        STOP
     ENDIF
     OPEN(UNIT=2231,FILE="./input/"//TRIM(data_filename))
     IOstatus=0
     nfile=0
     DO WHILE(IOstatus.EQ.0)
        READ(2231,*,IOSTAT=IOstatus)stringin
        IF(IOstatus.NE.0)EXIT
        IF(LEN(TRIM(stringin)).EQ.0)CYCLE
        stringin="./input/"//TRIM(stringin)
        INQUIRE(FILE=stringin,EXIST=lexist)
        IF(.NOT.lexist)THEN
           WRITE(*,*)"WARNING:Cannot find the file "//stringin
        ENDIF
        nfile=nfile+1
        IF(nfile.GT.MXFILE)THEN
           WRITE(*,*)"ERROR:Please enlarge MXFILE to be ",nfile
           CLOSE(UNIT=2231)
           STOP
        ENDIF
        stringfile(nfile)=stringin
     ENDDO
     CLOSE(UNIT=2231)
     IF(nfile.EQ.0)THEN
        WRITE(*,*)"ERROR: No valid data files chosen. Bye."
        STOP
     ENDIF
     ! read in the experimental data from files specified by data_list.inp
     ISET=0
     DO ii=1,nfile
        OPEN(UNIT=2231,FILE=TRIM(stringfile(ii)))
        CALL READ_EXPDATA(2231,ITYPE_SET,ISET,SQRTS_SET,Y_SET,NBIN_SET,ISTART_SET,XBIN,EXVAL)
        CLOSE(UNIT=2231)
     ENDDO
     WRITE(*,*)"Number of data to be fitted:",ISTART_SET(ISET)+NBIN_SET(ISET)
     FIT_NMAX=ISET
  ENDIF
  kapa=X(1)
  lam=X(2)
  ptavg=X(3)
  nn=X(4)
  CHISQ=0d0
  DO ii=1, ISET
     itype=ITYPE_SET(ii)
     sqrtS=SQRTS_SET(ii)
     ylow=Y_SET(ii,1)
     yup=Y_SET(ii,2)
     istart=ISTART_SET(ii)
     DO jj=1, NBIN_SET(ii)
        ptlow=XBIN(istart+jj,1)
        ptup=XBIN(istart+jj,2)
        THVAL=FUTIL(itype,sqrtS,ptlow,ptup,ylow,yup,kapa,lam,ptavg,nn)
        THVAL=THVAL/(ptup-ptlow)/(yup-ylow)
        IF(IFLAG.EQ.3)THEN
           FIT_TH_HIST(ii,jj)=THVAL
           FIT_TH_ERR(ii,jj)=0d0
        ENDIF
        CHI=(THVAL-EXVAL(istart+jj,1))**2/EXVAL(istart+jj,2)**2
        CHISQ=CHISQ+CHI
     ENDDO
  ENDDO
  F=CHISQ
  IF(IFLAG.EQ.3)THEN
     ! output the final values
     OPEN(UNIT=2235,FILE="./output/fit_results.out")
     WRITE(2235,*)"number of data = ",ISTART_SET(ISET)+NBIN_SET(ISET)
     WRITE(2235,*)"kapa = ",kapa
     WRITE(2235,*)"lam = ",lam
     WRITE(2235,*)"<pt> = ",ptavg
     WRITE(2235,*)"n = ",nn
     WRITE(2235,*)"chi^2 = ",F
     CLOSE(UNIT=2235)
     DO ii=1,ISET
        CALL FGNUPLOT(ii,ITYPE_SET(ii))
     ENDDO
  ENDIF
  RETURN
END SUBROUTINE FCN_FIT

SUBROUTINE READ_EXPDATA(iunit,ITYPE_SET,ISET,SQRTS_SET,Y_SET,NBIN_SET,ISTART_SET,XBIN,EXVAL)
  USE plot_fit_cb
  IMPLICIT NONE
  INTEGER,INTENT(IN)::iunit
  INTEGER,INTENT(INOUT)::ISET
  INTEGER,DIMENSION(MXSET),INTENT(INOUT)::ITYPE_SET
  !INTEGER,PARAMETER::MXBIN=2000,MXSET=20
  REAL(KIND(1d0)),DIMENSION(MXBIN,2),INTENT(INOUT)::XBIN,EXVAL 
  ! XBIN(1),XBIN(2) is Ptbin, EXVAL(1) is value and EXVAL(2) is err
  INTEGER,DIMENSION(MXSET),INTENT(INOUT)::NBIN_SET 
  ! the number of bins in each experimental data set
  INTEGER,DIMENSION(MXSET),INTENT(INOUT)::ISTART_SET
  REAL(KIND(1d0)),DIMENSION(MXSET),INTENT(INOUT)::SQRTS_SET
  REAL(KIND(1d0)),DIMENSION(MXSET,2),INTENT(INOUT)::Y_SET ! the y interval 
  INTEGER::IOstatus,IOstatus2,ii,ibin,itype=1
  CHARACTER*100::stringin,stringin2
  REAL(KIND(1d0))::value
  CHARACTER*3::unit
  IOstatus=0
  DO WHILE(IOstatus.EQ.0)
     READ(iunit,'(A100)',IOSTAT=IOstatus)stringin
     IF(LEN(TRIM(stringin)).EQ.0)CYCLE
     IF(stringin(1:1).EQ."#")THEN
        ii=INDEX(stringin,"d^2sigma/dpT/dxF")
        IF(ii.NE.0)THEN
           ITYPE=2
        ELSE
           ii=INDEX(stringin,"d^2sigma/dpT/dy")
           IF(ii.NE.0)THEN
              ITYPE=1
           ENDIF
        ENDIF
        CYCLE
     ENDIF
     BACKSPACE(UNIT=iunit)
     ISET=ISET+1
     IF(ISET.GT.MXSET)THEN
        WRITE(*,*)"ERROR:Please enlarge MXSET to be ", ISET
        CLOSE(iunit)
        STOP
     ENDIF
     ITYPE_SET(ISET)=ITYPE
     READ(iunit,*,IOSTAT=IOstatus)SQRTS_SET(ISET),Y_SET(ISET,1),Y_SET(ISET,2),NBIN_SET(ISET)
     IF(IOstatus.NE.0)THEN
        ISET=ISET-1
        EXIT
     ELSE
        IF(ISET.EQ.1)THEN
           ISTART_SET(ISET)=0
        ELSE
           ISTART_SET(ISET)=ISTART_SET(ISET-1)+NBIN_SET(ISET-1)
        ENDIF
     ENDIF
     IF(ISTART_SET(ISET)+NBIN_SET(ISET).GT.MXBIN)THEN
        WRITE(*,*)"ERROR:Please enlarge MXBIN to be ", ISTART_SET(ISET)+NBIN_SET(ISET)
        CLOSE(iunit)
        STOP
     ENDIF
     IOstatus2=0
     ibin=0
     DO WHILE(IOstatus2.EQ.0.AND.ibin.LT.NBIN_SET(ISET))
        READ(iunit,*,IOSTAT=IOstatus2)stringin2
        IF(LEN(TRIM(stringin2)).EQ.0)CYCLE
        IF(stringin2(1:1).EQ."#")CYCLE
        BACKSPACE(UNIT=iunit)
        DO ii=1,NBIN_SET(ISET)
           READ(iunit,*,IOSTAT=IOstatus2)XBIN(ISTART_SET(ISET)+ii,1),&
                XBIN(ISTART_SET(ISET)+ii,2),EXVAL(ISTART_SET(ISET)+ii,1),&
                EXVAL(ISTART_SET(ISET)+ii,2)
           FIT_XHIS(ISET,ii)=(XBIN(ISTART_SET(ISET)+ii,1)+XBIN(ISTART_SET(ISET)+ii,2))/2d0
           FIT_EX_HIST(ISET,ii)=EXVAL(ISTART_SET(ISET)+ii,1)
           FIT_EX_ERR(ISET,ii)=EXVAL(ISTART_SET(ISET)+ii,2)
           IF(ii.EQ.1)FIT_HMIN(ISET)=XBIN(ISTART_SET(ISET)+ii,1)
           IF(ii.EQ.NBIN_SET(ISET))FIT_HMAX(ISET)=XBIN(ISTART_SET(ISET)+ii,2)
        ENDDO
        ibin=NBIN_SET(ISET)
        FIT_NBIN(ISET)=NBIN_SET(ISET)
        IF(SQRTS_SET(ISET).GE.1000.AND.SQRTS_SET(ISET).LT.1d6)THEN
           unit="TeV"
           value=SQRTS_SET(ISET)/1000.
        ELSEIF(SQRTS_SET(ISET).LT.1000.AND.SQRTS_SET(ISET).GE.1d0)THEN
           unit="GeV"
           value=SQRTS_SET(ISET)
        ELSE
           WRITE(*,*)"ERROR:SQRTS is outside the plot range"
           STOP
        ENDIF
        SELECT CASE(ITYPE_SET(ISET))
        CASE(2)
           WRITE(FIT_LABEL(ISET),"(A21,F7.2,A3,A1,F5.2,A5,F5.2)")&
                "{/Symbol @\326\140}S=",value,unit,&
                ",",Y_SET(ISET,1),"<x_F<",Y_SET(ISET,2)
        CASE DEFAULT
           WRITE(FIT_LABEL(ISET),"(A21,F7.2,A3,A1,F6.2,A3,F6.2)")&
                "{/Symbol @\326\140}S=",value,unit,&
                ",",Y_SET(ISET,1),"<y<",Y_SET(ISET,2)
        END SELECT
     ENDDO
     IOstatus=IOstatus2
     ITYPE=1
  ENDDO
  RETURN
END SUBROUTINE READ_EXPDATA

FUNCTION CRYSTALBALL(itype,sqrtS,ptlow,ptup,ylow,yup,kapa,lam,ptavg,nn)
  USE fit_cb_global
  USE MC_VEGAS
  IMPLICIT NONE
  INTEGER::itype
  REAL(KIND(1d0))::sqrtS,ptlow,ptup,ylow,yup,kapa,lam,ptavg,nn
  REAL(KIND(1d0))::CRYSTALBALL
  REAL(KIND(1d0)),EXTERNAL::fit_cb_fxn,fit_cb_fxn2
  INTEGER::ii,varnum
  REAL(KIND(1d0))::vfes,sd,chi2a
  IF(lam.EQ.0d0.OR.kapa.LE.0d0)THEN
     CRYSTALBALL=0d0
     RETURN
  ENDIF
  varnum=3
  DO ii=1,varnum
     XL(ii)=0.0d0
     XU(ii)=1.0d0
  ENDDO
  ITMX=1
  NCALL=20000
  fit_sqrtS=sqrtS
  fit_ptlow=ptlow
  fit_ptup=ptup
  fit_ylow=ylow
  fit_yup=yup
  fit_kapa=kapa
  fit_lam=lam
  fit_ptavg=ptavg
  fit_nn=nn
  SELECT CASE(itype)
  CASE(2)
     ! d^2sigma/dpT/dxF
     CALL VEGAS(varnum,fit_cb_fxn2,vfes,sd,chi2a)
  CASE DEFAULT
     ! d^2sigma/dpT/dy
     CALL VEGAS(varnum,fit_cb_fxn,vfes,sd,chi2a)
  END SELECT
  CRYSTALBALL=vfes
  RETURN
END FUNCTION CRYSTALBALL

FUNCTION fit_cb_fxn(x,wgt)
  USE Func_PSI
  USE fit_cb_global
  IMPLICIT NONE
  ! varnum=3
  REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::x
  REAL(KIND(1d0)),INTENT(IN)::wgt
  REAL(KIND(1d0))::fit_cb_fxn
  REAL(KIND(1d0))::y1,y2,phi,w1
  !REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
  INTEGER::nwarn,nwmax,nwri,nwri_tot,icut,pdfnumpdf,ivarnum
  REAL(KIND(1d0))::sqs,sq,pt1c,y1cup,y1clow,ycollcm,mt1,mt2,pt
  SAVE sqs,sq,pt1c,y1cup,y1clow,ycollcm
  REAL(KIND(1d0))::wdps,wsf,wme,wgtbr,temp1,temp2,temp21,temp3,temp4
  REAL(KIND(1d0))::exp1,exp2,Jpt2,Jy1,Jy2,xp11,xp21,scale1,scale2
  REAL(KIND(1d0))::recmax=0,recmin=0
  INTEGER::j
  REAL(KIND(1d0)),EXTERNAL::CRYSTALBALL_ME
  wdps=1d0
  pt=(fit_ptup-fit_ptlow)*x(3)+fit_ptlow
  mt1=DSQRT(pt**2+mpsi**2)
  mt2=pt
  temp2=MIN(ACosh_p((fit_sqrtS**2+mpsi**2)/(2d0*fit_sqrtS*mt1)),fit_yup)
  temp21=MAX(-ACosh_p((fit_sqrtS**2+mpsi**2)/(2d0*fit_sqrtS*mt1)),fit_ylow)
  IF(temp2.LT.temp21)THEN
     fit_cb_fxn=0d0
     RETURN
  ENDIF
  y1=(temp2-temp21)*x(1)+temp21
  temp3=DLOG((-DEXP(-y1)*mt1+fit_sqrtS)/mt2)
  temp4=temp3+DLOG((-DEXP(y1)*mt1+fit_sqrtS)/mt2)
  y2=-temp3+temp4*x(2) ! in collision frame
  Jpt2=2d0*pt*(fit_ptup-fit_ptlow)
  Jy1=temp2-temp21
  Jy2=temp4
  fit_xp1=(DEXP(y1)*mt1+DEXP(y2)*mt2)/fit_sqrtS
  fit_xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2)/fit_sqrtS
  fit_ehat=DSQRT(fit_xp1*fit_xp2)*fit_sqrtS
  wdps=wdps*fit_xp1*fit_xp2*Jpt2*Jy1*Jy2*3.8938573d5/(16d0*pi*fit_ehat**4)
  phi=0d0
  fit_cb_pmom(1,1:2)=0
  fit_cb_pmom(1,3)=fit_xp1*fit_sqrtS/2d0
  fit_cb_pmom(1,4)=fit_xp1*fit_sqrtS/2d0
  fit_cb_pmom(2,1:2)=0
  fit_cb_pmom(2,3)=-fit_xp2*fit_sqrtS/2d0
  fit_cb_pmom(2,4)=fit_xp2*fit_sqrtS/2d0
  fit_cb_pmom(3,1)=pt*DCOS(phi)
  fit_cb_pmom(3,2)=pt*DSIN(phi)
  fit_cb_pmom(3,3)=mt1*DSINH(y1)
  fit_cb_pmom(3,4)=mt1*DCOSH(y1)
  fit_cb_pmom(4,1)=-pt*DCOS(phi)
  fit_cb_pmom(4,2)=-pt*DSIN(phi)
  fit_cb_pmom(4,3)=mt2*DSINH(y2)
  fit_cb_pmom(4,4)=mt2*DCOSH(y2)
  CALL strf_pdf_fit_cb(wsf)
  wdps=wdps*wsf
  IF(wdps.LE.0d0)THEN
     fit_cb_fxn=0d0
     RETURN
  ENDIF
  icut=1
  CALL Cuts_fit_cb(icut)
  IF(icut.EQ.0)THEN
     fit_cb_fxn=0d0
     RETURN
  ENDIF
  wme=CRYSTALBALL_ME(1)
  wdps=wdps*wme
  IF(wdps.LE.0d0)THEN
     fit_cb_fxn=0d0
     RETURN
  ENDIF
  fit_cb_fxn=wdps
  RETURN
END FUNCTION fit_cb_fxn

FUNCTION fit_cb_fxn2(x,wgt)
  USE Func_PSI
  USE fit_cb_global
  IMPLICIT NONE
  ! varnum=3
  REAL(KIND(1d0)),DIMENSION(3),INTENT(IN)::x
  REAL(KIND(1d0)),INTENT(IN)::wgt
  REAL(KIND(1d0))::fit_cb_fxn2
  REAL(KIND(1d0))::y1,y2,phi,w1
  INTEGER::nwarn,nwmax,nwri,nwri_tot,icut,pdfnumpdf,ivarnum
  REAL(KIND(1d0))::sqs,sq,pt1c,y1cup,y1clow,ycollcm,mt1,mt2,pt,xF1,xF2
  SAVE sqs,sq,pt1c,y1cup,y1clow,ycollcm
  REAL(KIND(1d0))::wdps,wsf,wme,wgtbr,temp1,temp2,temp21,temp3,temp4
  REAL(KIND(1d0))::exp1,exp2,Jpt2,Jy1,Jy2,xp11,xp21,scale1,scale2
  REAL(KIND(1d0))::recmax=0,recmin=0
  INTEGER::j
  REAL(KIND(1d0)),EXTERNAL::CRYSTALBALL_ME
  wdps=1d0
  pt=(fit_ptup-fit_ptlow)*x(3)+fit_ptlow
  mt1=DSQRT(pt**2+mpsi**2)
  mt2=pt
  ! The following are defined in the collision frame
  ! dyi=Sqrt[S]/2/mti/Cosh[yi]*dxFi
  ! xFi=(2*mti*Sinh[yi])/Sqrt[S]
  ! yi = Log[xFi*Sqrt[S]/mti/2 + Sqrt[1 + (xFi*Sqrt[S]/mti/2)^2]]
  temp2=(fit_sqrtS**2+mpsi**2)/(2d0*fit_sqrtS*mt1)
  temp21=-DSQRT(temp2**2-1d0)*2d0*mt1/fit_sqrtS
  temp2=MIN(-temp21,fit_yup)
  temp21=MAX(temp21,fit_ylow)
  IF(temp2.LT.temp21.OR.mt2.LT.1d-3)THEN
     fit_cb_fxn2=0d0
     RETURN
  ENDIF
  xF1=(temp2-temp21)*x(1)+temp21
  y1=DLOG(ABS(xF1*fit_sqrtS/mt1/2d0+DSQRT(1d0+(xF1*fit_sqrtS/mt1/2d0)**2)))
  temp3=DLOG((-DEXP(-y1)*mt1+fit_sqrtS)/mt2)
  temp3=-DSINH(temp3)*2d0*mt2/fit_sqrtS
  temp4=DLOG((-DEXP(y1)*mt1+fit_sqrtS)/mt2)
  temp4=DSINH(temp4)*2d0*mt2/fit_sqrtS
  xF2=(temp4-temp3)*x(2)+temp3 ! in collison frame
  y2=DLOG(ABS(xF2*fit_sqrtS/mt2/2d0+DSQRT(1d0+(xF2*fit_sqrtS/mt2/2d0)**2)))
  Jpt2=2d0*(fit_ptup-fit_ptlow) ! pt has been divided from Jaccobi dy2/dxF2
  Jy1=temp2-temp21
  Jy2=temp4-temp3
  fit_xp1=(DEXP(y1)*mt1+DEXP(y2)*mt2)/fit_sqrtS
  fit_xp2=(DEXP(-y1)*mt1+DEXP(-y2)*mt2)/fit_sqrtS
  fit_ehat=DSQRT(fit_xp1*fit_xp2)*fit_sqrtS
  wdps=wdps*Jpt2*Jy1*Jy2*3.8938573d5/(16d0*pi*fit_ehat**2)
  ! remaining Jacobi from dyi/dxFi
  wdps=wdps/4d0/mt1/DCOSH(y1)/DCOSH(y2)
  IF(wdps.LE.0d0.OR.wdps.NE.wdps)THEN
     fit_cb_fxn2=0d0
     RETURN
  ENDIF
  phi=0d0
  fit_cb_pmom(1,1:2)=0
  fit_cb_pmom(1,3)=fit_xp1*fit_sqrtS/2d0
  fit_cb_pmom(1,4)=fit_xp1*fit_sqrtS/2d0
  fit_cb_pmom(2,1:2)=0
  fit_cb_pmom(2,3)=-fit_xp2*fit_sqrtS/2d0
  fit_cb_pmom(2,4)=fit_xp2*fit_sqrtS/2d0
  fit_cb_pmom(3,1)=pt*DCOS(phi)
  fit_cb_pmom(3,2)=pt*DSIN(phi)
  fit_cb_pmom(3,3)=mt1*DSINH(y1)
  fit_cb_pmom(3,4)=mt1*DCOSH(y1)
  fit_cb_pmom(4,1)=-pt*DCOS(phi)
  fit_cb_pmom(4,2)=-pt*DSIN(phi)
  fit_cb_pmom(4,3)=mt2*DSINH(y2)
  fit_cb_pmom(4,4)=mt2*DCOSH(y2)
  CALL strf_pdf_fit_cb(wsf)
  wdps=wdps*wsf
  IF(wdps.LE.0d0)THEN
     fit_cb_fxn2=0d0
     RETURN
  ENDIF
  icut=1
  CALL Cuts_fit_cb2(icut)
  IF(icut.EQ.0)THEN
     fit_cb_fxn2=0d0
     RETURN
  ENDIF
  wme=CRYSTALBALL_ME(1)
  wdps=wdps*wme
  IF(wdps.LE.0d0)THEN
     fit_cb_fxn2=0d0
     RETURN
  ENDIF
  fit_cb_fxn2=wdps
  RETURN
END FUNCTION fit_cb_fxn2

SUBROUTINE strf_pdf_fit_cb(wsf)
  USE CTEQ6PDF
  USE Structf_PDFs
  USE fit_cb_global
  IMPLICIT NONE
  INTEGER::ipp=1
  INTEGER::ih=1 
  REAL(KIND(1d0)),INTENT(OUT)::wsf
  LOGICAL::use_cteq6_f90=.TRUE.
  INTEGER::init=0
  REAL(KIND=DBL),DIMENSION(-7:7)::pdflist
  SAVE init,ipp,use_cteq6_f90
  REAL(KIND(1d0))::glu1_ct,glu2_ct,u1_ct,u2_ct,d1_ct,d2_ct,s1_ct,s2_ct,c1_ct,c2_ct,b1_ct,b2_ct,&
       ub1_ct,ub2_ct,db1_ct,db2_ct,sb1_ct,sb2_ct,cb1_ct,cb2_ct,bb1_ct,bb2_ct,sf_ct
  INCLUDE "../lhapdf/call_strf_lhapdf"
  IF(init.EQ.0)THEN
     SELECT CASE(fit_iPDFSUP1)
     CASE(10000)
        CALL SetCtq6f90(1)
        use_cteq6_f90=.TRUE.
     CASE(10041)
        CALL SetCtq6f90(3)
        use_cteq6_f90=.TRUE.
     CASE(10042)
        CALL SetCtq6f90(4)
        use_cteq6_f90=.TRUE.
     CASE DEFAULT
        CALL pdfset_internal
        use_cteq6_f90=.FALSE.
     END SELECT
     ipp=2
     init=1
  ENDIF
  scale=fit_cb_pmom(3,1)**2+fit_cb_pmom(3,2)**2+mpsi**2
  scale=DSQRT(scale)

  IF(use_cteq6_f90)THEN
         glu1_ct = fit_xp1*Ctq6Pdf_f90(0,fit_xp1,scale)
         glu2_ct = fit_xp2*Ctq6Pdf_f90(0,fit_xp2,scale)
         u1_ct   = fit_xp1*Ctq6Pdf_f90(1,fit_xp1,scale)
         u2_ct   = fit_xp2*Ctq6Pdf_f90(1,fit_xp2,scale)
         d1_ct   = fit_xp1*Ctq6Pdf_f90(2,fit_xp1,scale)
         d2_ct   = fit_xp2*Ctq6Pdf_f90(2,fit_xp2,scale)
         s1_ct   = fit_xp1*Ctq6Pdf_f90(3,fit_xp1,scale)
         s2_ct   = fit_xp2*Ctq6Pdf_f90(3,fit_xp2,scale)
         c1_ct   = fit_xp1*Ctq6Pdf_f90(4,fit_xp1,scale)
         c2_ct   = fit_xp2*Ctq6Pdf_f90(4,fit_xp2,scale)
         b1_ct   = fit_xp1*Ctq6Pdf_f90(5,fit_xp1,scale)
         b2_ct   = fit_xp2*Ctq6Pdf_f90(5,fit_xp2,scale)
         ub1_ct  = fit_xp1*Ctq6Pdf_f90(-1,fit_xp1,scale)
         ub2_ct  = fit_xp2*Ctq6Pdf_f90(-1,fit_xp2,scale)
         db1_ct  = fit_xp1*Ctq6Pdf_f90(-2,fit_xp1,scale)
         db2_ct  = fit_xp2*Ctq6Pdf_f90(-2,fit_xp2,scale)
         sb1_ct  = fit_xp1*Ctq6Pdf_f90(-3,fit_xp1,scale)
         sb2_ct  = fit_xp2*Ctq6Pdf_f90(-3,fit_xp2,scale)
         cb1_ct  = fit_xp1*Ctq6Pdf_f90(-4,fit_xp1,scale)
         cb2_ct  = fit_xp2*Ctq6Pdf_f90(-4,fit_xp2,scale)
         bb1_ct  = fit_xp1*Ctq6Pdf_f90(-5,fit_xp1,scale)
         bb2_ct  = fit_xp2*Ctq6Pdf_f90(-5,fit_xp2,scale)
      ELSE
         CALL fdist(ih,fit_xp1,scale,pdflist(-7:7))
         glu1_ct = fit_xp1*pdflist(0)
         u1_ct   = fit_xp1*pdflist(2)
         d1_ct   = fit_xp1*pdflist(1)
         s1_ct   = fit_xp1*pdflist(3)
         c1_ct   = fit_xp1*pdflist(4)
         b1_ct   = fit_xp1*pdflist(5)
         ub1_ct  = fit_xp1*pdflist(-2)
         db1_ct  = fit_xp1*pdflist(-1)
         sb1_ct  = fit_xp1*pdflist(-3)
         cb1_ct  = fit_xp1*pdflist(-4)
         bb1_ct  = fit_xp1*pdflist(-5)
         CALL fdist(ih,fit_xp2,scale,pdflist(-7:7))
         glu2_ct = fit_xp2*pdflist(0)
         u2_ct   = fit_xp2*pdflist(2)
         d2_ct   = fit_xp2*pdflist(1)
         s2_ct   = fit_xp2*pdflist(3)
         c2_ct   = fit_xp2*pdflist(4)
         b2_ct   = fit_xp2*pdflist(5)
         ub2_ct  = fit_xp2*pdflist(-2)
         db2_ct  = fit_xp2*pdflist(-1)
         sb2_ct  = fit_xp2*pdflist(-3)
         cb2_ct  = fit_xp2*pdflist(-4)
         bb2_ct  = fit_xp2*pdflist(-5)
      ENDIF
      wsf=glu1_ct*glu2_ct/fit_xp1/fit_xp2
      IF(fit_includeqq)THEN
         wsf=wsf+(u1_ct*ub2_ct+d1_ct*db2_ct+s1_ct*sb2_ct+c1_ct*cb2_ct+&
              ub1_ct*u2_ct+db1_ct*d2_ct+sb1_ct*s2_ct+cb1_ct*c2_ct)/fit_xp1/fit_xp2
      ENDIF
      IF(init.EQ.0)init=1
END SUBROUTINE strf_pdf_fit_cb

SUBROUTINE Cuts_fit_cb(icut)
  USE fit_cb_global
  USE Kinetic_Func
  IMPLICIT NONE
  INTEGER,INTENT(OUT)::icut
  REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
  REAL(KIND(1d0))::pt,y
  icut=0
  pt=SQRT(fit_cb_pmom(3,1)**2+fit_cb_pmom(3,2)**2)
  IF(pt.GT.fit_ptup.OR.pt.LT.fit_ptlow)RETURN
  y=rapidity(fit_cb_pmom(3,1:4))
  IF(y.GT.fit_yup.OR.y.LT.fit_ylow)RETURN
  icut=1
  RETURN
END SUBROUTINE Cuts_fit_cb

SUBROUTINE Cuts_fit_cb2(icut)
  USE fit_cb_global
  USE Kinetic_Func
  IMPLICIT NONE
  INTEGER,INTENT(OUT)::icut
  REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
  REAL(KIND(1d0))::pt,y
  icut=0
  pt=SQRT(fit_cb_pmom(3,1)**2+fit_cb_pmom(3,2)**2)
  IF(pt.GT.fit_ptup.OR.pt.LT.fit_ptlow)RETURN
  y=xFeynman(fit_cb_pmom(3,1:4),fit_sqrtS)
  IF(y.GT.fit_yup.OR.y.LT.fit_ylow)RETURN
  icut=1
  RETURN
END SUBROUTINE Cuts_fit_cb2

FUNCTION CRYSTALBALL_ME(idummy)
  USE fit_cb_global
  IMPLICIT NONE
  INTEGER,INTENT(IN)::idummy ! dummy argument
  REAL(KIND(1d0))::CRYSTALBALL_ME
  REAL(KIND(1d0))::shat,pt,KK
  KK=fit_lam**2*fit_kapa/mpsi**2
  shat=(fit_cb_pmom(1,1)+fit_cb_pmom(2,1))**2&
       +(fit_cb_pmom(1,2)+fit_cb_pmom(2,2))**2&
       +(fit_cb_pmom(1,3)+fit_cb_pmom(2,3))**2
  shat=(fit_cb_pmom(1,4)+fit_cb_pmom(2,4))**2-shat
  pt=DSQRT(fit_cb_pmom(3,1)**2+fit_cb_pmom(3,2)**2)
  IF(pt.LE.fit_ptavg)THEN
     CRYSTALBALL_ME=KK*shat*DEXP(-fit_kapa*pt**2/mpsi**2)
  ELSE
     CRYSTALBALL_ME=KK*shat*DEXP(-fit_kapa*fit_ptavg**2/mpsi**2)&
          *1d0/(1d0+fit_kapa/fit_nn*(pt**2-fit_ptavg**2)/mpsi**2)**fit_nn
  ENDIF
  RETURN
END FUNCTION CRYSTALBALL_ME
