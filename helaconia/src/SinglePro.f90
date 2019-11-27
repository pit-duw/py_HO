MODULE SinglePro
IMPLICIT NONE
CONTAINS
SUBROUTINE SinProcess
USE Helac_Global
USE Helac_master
USE Helac_histo
USE CONSTANTS
USE unweight_lhe
USE Colliders_PSI_1
USE Colliders_PSI_2
USE DecayInfo
IMPLICIT NONE
CHARACTER(len=10),DIMENSION(20)::char
!CHARACTER(len=400)::process
INTEGER::i,i1,j1,k1,s1x    !,gener
INTEGER::nh,innum,outnum
REAL(KIND=DBL)::w0,ptp,pqp,wme
INTEGER::initsig=0,ioerror,icase=0
LOGICAL::lunwei=.FALSE.,lhewgt=.FALSE.
REAL(KIND=DBL),DIMENSION(3)::rslt
INTEGER::itmxn,ncalln
REAL(KIND=DBL),DIMENSION(4)::pmomtemp
LOGICAL::lexist
SAVE initsig,lunwei,lhewgt,icase
CALL Helac_zmset(34)
!nunit1=100                 ! standard output
!IF(nunit1.NE.6)THEN
!	OPEN(UNIT=nunit1,FILE="RESULT.out")
!ENDIF
!WRITE(*,*)'histo file'
!READ(*,*)histo
!WRITE(*,*)'error file'
!READ(*,*)error

IF(initsig.EQ.0)THEN
WRITE(*,*)'                     THE BEGINNING OF HELAC-Onia'
WRITE(*,*)' '
WRITE(*,*)'======================================================================='
WRITE(*,*)'======================================================================='
WRITE(*,*)' '
WRITE(*,*)'histo file: ',TRIM(histo_file),'    error file: ',TRIM(error_file)
WRITE(*,*)' '
WRITE(*,*)'-----------------------------------------------------------------------'
WRITE(*,*)' '
CALL ReadElem_real('repeat',repeat)
WRITE(*,*)'repeat:'
WRITE(*,*)'0 run both phase of Helac_Phegas; 1/2 run only the first/second phase .'
!READ(*,*)repeat
WRITE(*,*)' '
CALL ReadElem_integer('ranhel',iranhel)
WRITE(*,*)'ranhel:'
WRITE(*,*)'0 calculate explicitly all helicity conf; 1 do a MC over helicities .'
!READ(*,*)iranhel
WRITE(*,*)' '
WRITE(*,*)'repeat value=  ',repeat,'  ranhel value=  ',iranhel
CALL ReadElem_integer('colpar',Coll_Type)
CALL ReadElem_logic('exp3pjQ',exp3pjQ)
CALL ReadElem_logic('OctetQ',octetQ)
CALL ReadElem_integer('modes',imode)
CALL ReadElem_logic('MeasureSpeed',measurespeed)
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
WRITE(*,*)' '
WRITE(*,*)'modes=',imode
IF(imode.EQ.1)THEN
	! iLSJ only will be used when the first quarkonium is 3pj
	CALL ReadElem_integer('LSJ',iLSJ)
	WRITE(*,*)'LSJ=',iLSJ
	CALL ReadElem_integer('SDME1',ihel1)
	WRITE(*,*)'SDME1=',ihel1
	CALL ReadElem_integer('SDME2',ihel2)
	WRITE(*,*)'SDME2=',ihel2
	CALL ReadElem_integer('PolarFrame',PolarFrame)
	WRITE(*,*)'PolarFrame=',PolarFrame
	IF(ihel1.NE.ihel2)THEN
		CALL ReadElem_integer('SDPart',SDPart)
		WRITE(*,*)'SDPart=',SDPart
	ENDIF
ENDIF
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
IF(uselhapdf)THEN
   CALL ReadElem_logic('reweight_pdf',reweight_pdf)
   IF(reweight_pdf)THEN
      CALL ReadElem_integer('pdf_min',ipdf_set_min)
      CALL ReadElem_integer('pdf_max',ipdf_set_max)
      ho_npdf=ipdf_set_max-ipdf_set_min+1
   ENDIF
ENDIF
CALL ReadElem_logic('quarksumQ',qsumQ)
WRITE(*,*)'Sum quark ?',qsumQ
IF(qsumQ)THEN
	CALL ReadElem_integer('iqnum',iqnum)
	WRITE(*,*)'number of flavor in initial state:',iqnum
ENDIF
CALL ReadElem_logic("absoluterap",absrap)
IF(.NOT.SUMMATIONQ)THEN
	WRITE(*,*)' '
	WRITE(*,*)'------------------------------------------------------------------------'
	WRITE(*,*)' '
!	WRITE(*,*)'Number of particles :'
	CLOSE(21711)
	OPEN(21711,FILE=TRIM(input_dir)//'process.inp')
	READ(21711,*)nhad
	WRITE(*,*)'The number of hadron particles in the process is :',nhad
	WRITE(*,*)' '
	WRITE(*,*)'Flavour of particles in SM : '
	WRITE(*,*)'nue=1 e=2 u=3 d=4 numu=5 mu=6 c=7 s=8 nutau=9 tau=10 t=11 b=12'
	WRITE(*,*)'A=31 Z=32 W+=33 W-=34 g=35 H=41 G0=42 G+=43 G-=44 Jet=100'
	WRITE(*,*)'all antifermions are corresponding negative integers '
	WRITE(*,*)'Jets are only allowed in the summation mode'
	WRITE(*,*)'etac(1S01)=441001, etac(1S08)=441008, psi(3S11)=443011'
	WRITE(*,*)'psi(3S18)=443018,hc(1P11)=441111,hc(1P18)=441118'
	WRITE(*,*)"chicJ(3PJ1)=4431J1,chicJ(3PJ8)=4431J8,"
	WRITE(*,*)" similar to Bottomnium but with 44 replaced by 55"
	WRITE(*,*)" Bc+(3S11)=453011,Bc+(3S18)=453018,Bc-(3S11)=-453011,.."
	READ(21711,*)(iflh(i),i=1,nhad)
	CLOSE(21711)
	WRITE(*,*)(iflh(i),i=1,nhad)
        ! begin of ISR shower stuff
        CALL ReadElem_integer('emep_ISR_shower',emep_ISR_shower)
        IF((ABS(iflh(1)).NE.2.OR.ABS(iflh(2)).NE.2.OR.iflh(1)+iflh(2).NE.0.OR.COLL_TYPE.NE.3).AND.emep_ISR_shower.NE.0)THEN
           WRITE(*,*)"WARNING:shower on e- e+ will be swich off !"
           emep_ISR_shower=0
        END IF
        ! make sure the first beam is e- while the second beam is e+ when using QEDPS
        IF(emep_ISR_shower.EQ.1)THEN
           IF(iflh(1).EQ.-2.AND.iflh(2).EQ.2)THEN
              STOP "Please make sure e- e+ (not e+ e-) when using QEDPS"
           ENDIF
        ELSE
           emep_ISR_shower=0
        ENDIF
        ! end of ISR shower stuff
	n=nhad
	j1=1
	Quarkonium2(1:20)=0
	Quarkonium3(1:20)=0
	Quarkonium(1:10,1:3)=0
	s1x=1
	Pwavenum=0
	Pwave(1:10)=0
	PwaveQ(1:20)=0
	Num3pj=0
	Pwave3pj(1:10)=0
	JNum3pj(1:10)=-1
	Pwave3pjQ(1:20)=0
	hadron2parton(1:20)=0
	parton2hadron(1:20)=0
        parton2hadrontype(1:20)=0
	octetnum=0
	octetlist(1:10)=0
	DO i=1,nhad
		hadron2parton(i)=j1
		parton2hadron(j1)=i
		IF(ABS(iflh(i)).LE.100)THEN
			ifl(j1)=iflh(i)
			j1=j1+1
		ELSE
			parton2hadron(j1+1)=i
			IF(iflh(i).GE.550000)THEN
			! Bottonium
				ifl(j1)=12
				ifl(j1+1)=-12
				Quarkonium(s1x,1)=i
				Quarkonium(s1x,2)=j1
				Quarkonium(s1x,3)=j1+1
				Quarkonium2(i)=s1x
				Quarkonium3(j1)=s1x
				Quarkonium3(j1+1)=s1x
                                parton2hadrontype(j1)=2
                                parton2hadrontype(j1+1)=2
				IF(SubInteger(iflh(i),3,3).EQ.1)THEN
					Pwavenum=Pwavenum+1
					Pwave(Pwavenum)=i
					PwaveQ(i)=Pwavenum
					IF(SubInteger(iflh(i),4,4).EQ.3)THEN
						Num3pj=Num3pj+1
						Pwave3pj(Num3pj)=i
						Pwave3pjQ(i)=Num3pj
						JNum3pj(Num3pj)=SubInteger(iflh(i),2,2)
					ENDIF
				ENDIF
				IF(SubInteger(iflh(i),1,1).EQ.8)THEN
					octetnum=octetnum+1
					octetlist(octetnum)=i
				ENDIF
				s1x=s1x+1
			ELSEIF(iflh(i).LE.449999.AND.iflh(i).GE.440000)THEN
			! Charmonium
				ifl(j1)=7
				ifl(j1+1)=-7
				Quarkonium(s1x,1)=i
				Quarkonium(s1x,2)=j1
				Quarkonium(s1x,3)=j1+1
				Quarkonium2(i)=s1x
				Quarkonium3(j1)=s1x
				Quarkonium3(j1+1)=s1x
                                parton2hadrontype(j1)=1
                                parton2hadrontype(j1+1)=1
				IF(SubInteger(iflh(i),3,3).EQ.1)THEN
					Pwavenum=Pwavenum+1
					Pwave(Pwavenum)=i
					PwaveQ(i)=Pwavenum
					IF(SubInteger(iflh(i),4,4).EQ.3)THEN
						Num3pj=Num3pj+1
						Pwave3pj(Num3pj)=i
						Pwave3pjQ(i)=Num3pj
						JNum3pj(Num3pj)=SubInteger(iflh(i),2,2)
					ENDIF
				ENDIF
				IF(SubInteger(iflh(i),1,1).EQ.8)THEN
					octetnum=octetnum+1
					octetlist(octetnum)=i
				ENDIF
				s1x=s1x+1
			ELSEIF(ABS(iflh(i)).GE.450000.AND.ABS(iflh(i)).LE.459999)THEN
			! Bc systems
				IF(iflh(i).GT.0)THEN
					ifl(j1)=7
					ifl(j1+1)=-12
				ELSE
					ifl(j1)=12
					ifl(j1+1)=-7
				ENDIF
				Quarkonium(s1x,1)=i
				Quarkonium(s1x,2)=j1
				Quarkonium(s1x,3)=j1+1
				Quarkonium2(i)=s1x
				Quarkonium3(j1)=s1x
				Quarkonium3(j1+1)=s1x
                                parton2hadrontype(j1)=3
                                parton2hadrontype(j1+1)=3
				IF(SubInteger(iflh(i),3,3).EQ.1)THEN
					Pwavenum=Pwavenum+1
					Pwave(Pwavenum)=i
					PwaveQ(i)=Pwavenum
					IF(SubInteger(iflh(i),4,4).EQ.3)THEN
						Num3pj=Num3pj+1
						Pwave3pj(Num3pj)=i
						Pwave3pjQ(i)=Num3pj
						JNum3pj(Num3pj)=SubInteger(iflh(i),2,2)
					ENDIF
				ENDIF
				IF(SubInteger(iflh(i),1,1).EQ.8)THEN
					octetnum=octetnum+1
					octetlist(octetnum)=i
				ENDIF
				s1x=s1x+1
			ELSE
				PRINT *,"Wrong in flavors of hadron in process ! STOP !"
				STOP
			ENDIF
			j1=j1+2
			n=n+1
		ENDIF
	ENDDO
	WRITE(*,*)' '
	WRITE(*,*)'The process at the parton level is'
	WRITE(*,*)(ifl(i),i=1,n)
	WRITE(*,*)' '
ENDIF
IF(imode.EQ.1.AND.iLSJ.EQ.2.AND..NOT.exp3pjQ.AND.QN3PJF(Quarkonium(1,1)))THEN
	PRINT *,"Please set exp3pjQ=T !"
	STOP
ELSEIF(imode.EQ.1.AND.iLSJ.NE.2.AND.exp3pjQ.AND.QN3PJF(Quarkonium(1,1)))THEN
	PRINT *,"Please set iLSJ = 2 or exp3pjQ=.FALSE.!"
	STOP
ENDIF
MCoHelicity=.FALSE.
CALL ReadDecayInfo
IF(NDecayChains.GT.0.OR.iranhel.EQ.4)THEN
   iranhel=0
   MCoHelicity=.TRUE.
ENDIF
IF(exp3pjQ.AND.Num3pj.GE.1.AND.iranhel.GE.2)THEN
	PRINT *,"Please set ranhel < 2 !"
	STOP
ELSEIF(QN3S1F(Quarkonium(1,1)).AND.imode.EQ.1.AND.iranhel.GE.3)THEN
	PRINT *,"Please set ranhel < 3 !"
	STOP
ELSEIF(QN1P1F(Quarkonium(1,1)).AND.imode.EQ.1.AND.iranhel.GE.2)THEN
	PRINT *,"Please set ranhel < 2 !"
	STOP
ELSEIF(QN3PJF(Quarkonium(1,1)).AND.imode.EQ.1.AND.iLSJ.EQ.0.AND.iranhel.GE.2)THEN	
	PRINT *,"Please set ranhel < 2 !"
	STOP
ELSEIF(QN3PJF(Quarkonium(1,1)).AND.imode.EQ.1.AND.iLSJ.EQ.1.AND.iranhel.GE.3)THEN
	PRINT *,"Please set ranhel < 3 !"
	STOP
ELSEIF(QN1S0F(Quarkonium(1,1)).AND.imode.EQ.1)THEN
	IF(ihel1.NE.0.OR.ihel2.NE.0)THEN
		PRINT *,"Please set SDME1=0 and SDME2=0 !"
		STOP
	ENDIF
	imode=0
ENDIF
! ENDIF

!IF(initsig.EQ.0)THEN
WRITE(*,*)'-------------------------------------------------------------------------'
WRITE(*,*)' '
WRITE(*,*)' FLAGS '
WRITE(*,*)' '
!WRITE(*,*)'iflag,iunitary,ihiggs,iwidth'
WRITE(*,*)'iflag: 0 sum over all the helicity confs, 1 choose the specific helicity conf '
CALL ReadElem_integer('iflag',iflag)
!READ(*,*)iflag
WRITE(*,*)'iflag value = ', iflag
WRITE(*,*)' '
WRITE(*,*)"iunitary: 1 use unitary gauge, 0 use 't Hooft-Feynman gauge "
CALL ReadElem_integer('gauge',iunitary)
!READ(*,*)iunitary
WRITE(*,*)'iunitary value = ',iunitary
WRITE(*,*)' '
WRITE(*,*)'ihiggs: include(1) or not(0) the higgs as an intermediate state '
!READ(*,*)ihiggs
CALL ReadElem_integer('ihiggs',ihiggs)
WRITE(*,*)'ihiggs value = ',ihiggs
WRITE(*,*)' '
WRITE(*,*)'iwidth: use(1) or not(0) the complex mass scheme for W and Z '
!READ(*,*)iwidth
CALL ReadElem_integer('widsch',iwidth)
WRITE(*,*)'iwidth value = ',iwidth
!READ(*,*)iflag,iunitary,ihiggs,iwidth
!WRITE(*,*)iflag,iunitary,ihiggs,iwidth
CALL Helac_FeynRule_SM()          ! choose the feynman rules
DO i=1,n
   IF(i.LE.2)THEN
      io(i)=1                     ! the first two are incoming particles
   ELSE
      io(i)=-1                    ! the rest are outgoing particles
   ENDIF
ENDDO
DO i=1,nhad
   IF(i.LE.2)THEN
      ioh(i)=1                     ! the first two are incoming particles
   ELSE
      ioh(i)=-1                    ! the rest are outgoing particles
   ENDIF
ENDDO
!initsig=1
ENDIF
LDMEwt=1d0
char(1:nhad)=' '
DO i=1,nhad
   ! parton
   IF(iflh(i).EQ. 2)char(i)='e'
   IF(iflh(i).EQ.-2)char(i)='ebar'
   IF(iflh(i).EQ. 1)char(i)='nue'
   IF(iflh(i).EQ.-1)char(i)='nuebar'
   IF(iflh(i).EQ. 3)char(i)='u'
   IF(iflh(i).EQ.-3)char(i)='ubar'
   IF(iflh(i).EQ. 4)char(i)='d'
   IF(iflh(i).EQ.-4)char(i)='dbar'
   IF(iflh(i).EQ. 6)char(i)='mu'
   IF(iflh(i).EQ.-6)char(i)='mubar'
   IF(iflh(i).EQ. 5)char(i)='numu'
   IF(iflh(i).EQ.-5)char(i)='numubar'
   IF(iflh(i).EQ. 8)char(i)='s'
   IF(iflh(i).EQ. 7)char(i)='c'
   IF(iflh(i).EQ.-7)char(i)='cbar'
   IF(iflh(i).EQ.-8)char(i)='sbar'
   IF(iflh(i).EQ. 10)char(i)='tau'
   IF(iflh(i).EQ.-10)char(i)='taubar'
   IF(iflh(i).EQ. 9)char(i)='nutau'
   IF(iflh(i).EQ.-9)char(i)='nutaubar'
   IF(iflh(i).EQ. 11)char(i)='t'
   IF(iflh(i).EQ.-11)char(i)='tbar'
   IF(iflh(i).EQ. 12)char(i)='b'
   IF(iflh(i).EQ.-12)char(i)='bbar'
   IF(iflh(i).EQ.31)char(i)='A'
   IF(iflh(i).EQ.32)char(i)='Z'
   IF(iflh(i).EQ.33)char(i)='W+'
   IF(iflh(i).EQ.34)char(i)='W-'
   IF(iflh(i).EQ.41)char(i)='H'
   IF(iflh(i).EQ.35)char(i)='g'
   ! charmonium
   IF(iflh(i).EQ.441001)THEN
		CALL ReadElem_real('LDMEcc1S01',LDME441001)
		char(i)='etac1'
		LDMEwt=LDMEwt*LDME441001/parmas(7)
   ENDIF
   IF(iflh(i).EQ.441008)THEN
   		CALL ReadElem_real('LDMEcc1S08',LDME441008)
		char(i)='etac8'
		LDMEwt=LDMEwt*LDME441008/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443011)THEN
   		CALL ReadElem_real('LDMEcc3S11',LDME443011)
		char(i)='psi1'
		LDMEwt=LDMEwt*LDME443011/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443018)THEN
      	CALL ReadElem_real('LDMEcc3S18',LDME443018)
		char(i)='psi8'
		LDMEwt=LDMEwt*LDME443018/parmas(7)
   ENDIF
   IF(iflh(i).EQ.441111)THEN
      	CALL ReadElem_real('LDMEcc1P11',LDME441111)
		char(i)='hc1'
		LDMEwt=LDMEwt*LDME441111/parmas(7)
   ENDIF
   IF(iflh(i).EQ.441118)THEN
		CALL ReadElem_real('LDMEcc1P18',LDME441118)
		char(i)='hc8'
		LDMEwt=LDMEwt*LDME441118/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443101)THEN
   		CALL ReadElem_real('LDMEcc3P01',LDME443101)
		char(i)='chic01'
		LDMEwt=LDMEwt*LDME443101/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443108)THEN
		CALL ReadElem_real('LDMEcc3P08',LDME443108)
		char(i)='chic08'
		LDMEwt=LDMEwt*LDME443108/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443111)THEN
   		CALL ReadElem_real('LDMEcc3P11',LDME443111)
		char(i)='chic11'
		LDMEwt=LDMEwt*LDME443111/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443118)THEN
		CALL ReadElem_real('LDMEcc3P18',LDME443118)
		char(i)='chic18'
		LDMEwt=LDMEwt*LDME443118/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443121)THEN
   		CALL ReadElem_real('LDMEcc3P21',LDME443121)
		char(i)='chic21'
		LDMEwt=LDMEwt*LDME443121/parmas(7)
   ENDIF
   IF(iflh(i).EQ.443128)THEN
		CALL ReadElem_real('LDMEcc3P28',LDME443128)
		char(i)='chic28'
		LDMEwt=LDMEwt*LDME443128/parmas(7)
   ENDIF
   ! bottonium
   IF(iflh(i).EQ.551001)THEN
		CALL ReadElem_real('LDMEbb1S01',LDME551001)
		char(i)='etab1'
		LDMEwt=LDMEwt*LDME551001/parmas(12)
   ENDIF
   IF(iflh(i).EQ.551008)THEN
		CALL ReadElem_real('LDMEbb1S08',LDME551008)
		char(i)='etab8'
		LDMEwt=LDMEwt*LDME551008/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553011)THEN
   		CALL ReadElem_real('LDMEbb3S11',LDME553011)
		char(i)='upsilon1'
		LDMEwt=LDMEwt*LDME553011/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553018)THEN
		CALL ReadElem_real('LDMEbb3S18',LDME553018)
		char(i)='upsilon8'
		LDMEwt=LDMEwt*LDME553018/parmas(12)
   ENDIF
   IF(iflh(i).EQ.551111)THEN
		CALL ReadElem_real('LDMEbb1P11',LDME551111)
		char(i)='hb1'
		LDMEwt=LDMEwt*LDME551111/parmas(12)
   ENDIF
   IF(iflh(i).EQ.551118)THEN
		CALL ReadElem_real('LDMEbb1P18',LDME551118)
		char(i)='hb8'
		LDMEwt=LDMEwt*LDME551118/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553101)THEN
		CALL ReadElem_real('LDMEbb3P01',LDME553101)
		char(i)='chib01'
		LDMEwt=LDMEwt*LDME553101/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553108)THEN
		CALL ReadElem_real('LDMEbb3P08',LDME553108)
		char(i)='chib08'
		LDMEwt=LDMEwt*LDME553108/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553111)THEN
		CALL ReadElem_real('LDMEbb3P11',LDME553111)
		char(i)='chib11'
		LDMEwt=LDMEwt*LDME553111/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553118)THEN
		CALL ReadElem_real('LDMEbb3P18',LDME553118)
		char(i)='chib18'
		LDMEwt=LDMEwt*LDME553118/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553121)THEN
		CALL ReadElem_real('LDMEbb3P21',LDME553121)
		char(i)='chib21'
		LDMEwt=LDMEwt*LDME553121/parmas(12)
   ENDIF
   IF(iflh(i).EQ.553128)THEN
		CALL ReadElem_real('LDMEbb3P28',LDME553128)
		char(i)='chib28'
		LDMEwt=LDMEwt*LDME553128/parmas(12)
   ENDIF
   ! Bc system
   IF(iflh(i).EQ.453011)THEN
   		CALL ReadElem_real('LDMEbc3S11',LDME453011)
		char(i)='Bc+1'
		LDMEwt=LDMEwt*LDME453011*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453018)THEN
   		CALL ReadElem_real('LDMEbc3S18',LDME453018)
		char(i)='Bc+8'
		LDMEwt=LDMEwt*LDME453018*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453011)THEN
   		CALL ReadElem_real('LDMEbc3S11',LDME453011)
		char(i)='Bc-1'
		LDMEwt=LDMEwt*LDME453011*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453018)THEN
   		CALL ReadElem_real('LDMEbc3S18',LDME453018)
		char(i)='Bc-8'
		LDMEwt=LDMEwt*LDME453018*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.451001)THEN
   		CALL ReadElem_real('LDMEbc1S01',LDME451001)
		char(i)='Bc+1S01'
		LDMEwt=LDMEwt*LDME451001*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-451001)THEN
   		CALL ReadElem_real('LDMEbc1S01',LDME451001)
		char(i)='Bc-1S01'
		LDMEwt=LDMEwt*LDME451001*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.451008)THEN
   		CALL ReadElem_real('LDMEbc1S08',LDME451008)
		char(i)='Bc+1S08'
		LDMEwt=LDMEwt*LDME451008*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-451008)THEN
   		CALL ReadElem_real('LDMEbc1S08',LDME451008)
		char(i)='Bc-1S08'
		LDMEwt=LDMEwt*LDME451008*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.451111)THEN
   		CALL ReadElem_real('LDMEbc1P11',LDME451111)
		char(i)='Bc+1P11'
		LDMEwt=LDMEwt*LDME451111*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-451111)THEN
   		CALL ReadElem_real('LDMEbc1P11',LDME451111)
		char(i)='Bc-1P11'
		LDMEwt=LDMEwt*LDME451111*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.451118)THEN
   		CALL ReadElem_real('LDMEbc1P18',LDME451118)
		char(i)='Bc+1P18'
		LDMEwt=LDMEwt*LDME451118*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-451118)THEN
   		CALL ReadElem_real('LDMEbc1P18',LDME451118)
		char(i)='Bc-1P18'
		LDMEwt=LDMEwt*LDME451118*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453101)THEN
   		CALL ReadElem_real('LDMEbc3P01',LDME453101)
		char(i)='Bc+3P01'
		LDMEwt=LDMEwt*LDME453101*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453101)THEN
   		CALL ReadElem_real('LDMEbc3P01',LDME453101)
		char(i)='Bc-3P01'
		LDMEwt=LDMEwt*LDME453101*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453108)THEN
   		CALL ReadElem_real('LDMEbc3P08',LDME453108)
		char(i)='Bc+3P08'
		LDMEwt=LDMEwt*LDME453108*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453108)THEN
   		CALL ReadElem_real('LDMEbc3P08',LDME453108)
		char(i)='Bc-3P08'
		LDMEwt=LDMEwt*LDME453108*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453111)THEN
   		CALL ReadElem_real('LDMEbc3P11',LDME453111)
		char(i)='Bc+3P11'
		LDMEwt=LDMEwt*LDME453111*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453111)THEN
   		CALL ReadElem_real('LDMEbc3P11',LDME453111)
		char(i)='Bc-3P11'
		LDMEwt=LDMEwt*LDME453111*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453118)THEN
   		CALL ReadElem_real('LDMEbc3P18',LDME453118)
		char(i)='Bc+3P18'
		LDMEwt=LDMEwt*LDME453118*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453118)THEN
   		CALL ReadElem_real('LDMEbc3P18',LDME453118)
		char(i)='Bc-3P18'
		LDMEwt=LDMEwt*LDME453118*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453121)THEN
   		CALL ReadElem_real('LDMEbc3P21',LDME453121)
		char(i)='Bc+3P21'
		LDMEwt=LDMEwt*LDME453121*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453121)THEN
   		CALL ReadElem_real('LDMEbc3P21',LDME453121)
		char(i)='Bc-3P21'
		LDMEwt=LDMEwt*LDME453121*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.453128)THEN
   		CALL ReadElem_real('LDMEbc3P28',LDME453128)
		char(i)='Bc+3P28'
		LDMEwt=LDMEwt*LDME453128*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
   IF(iflh(i).EQ.-453128)THEN
   		CALL ReadElem_real('LDMEbc3P28',LDME453128)
		char(i)='Bc-3P28'
		LDMEwt=LDMEwt*LDME453128*(parmas(7)+parmas(12))/(parmas(7)*parmas(12))/2d0
   ENDIF
ENDDO

WRITE(process,*)(TRIM(char(i)),i=1,nhad)
IF(.NOT.SUMMATIONQ)THEN
	WRITE(*,*)' '
	WRITE(*,*)'The process we are calculating is: ',TRIM(process(2:10*nhad+1))
	WRITE(*,*)' '
ELSE
	WRITE(*,*)' '
	WRITE(*,*)'The subprocess we are calculating is: ',TRIM(process(2:10*nhad+1))
	WRITE(*,*)' '
ENDIF

IF(nunit1.NE.6)THEN
	OPEN(UNIT=nunit1,FILE=TRIM(output_dir)//"RESULT_"//TRIM(process(2:10*nhad+1))//".out")
ENDIF
WRITE(nunit1,*)"LDMEwt=",LDMEwt
nunit2=32
CLOSE(nunit2)
OPEN(nunit2,FILE=TRIM(output_dir)//'kine_'//TRIM(process(2:10*nhad+1))//'.out')
nunit3=30
CLOSE(nunit3)
OPEN(nunit3,FILE=TRIM(tmp_dir)//'even_'//TRIM(process(2:10*nhad+1))//'.out',FORM='unformatted')

! for LHA
CLOSE(200)
OPEN(200,FILE=TRIM(tmp_dir)//'sample'//TRIM(process(2:10*nhad+1))//'.init')
! Scale scheme
CALL ReadElem_integer('Scale',nscheme)
CALL ReadElem_real('FScaleValue',fschemevalue)
CALL ReadELem_real('ScaleFactor',scalefactor)

IF(repeat.EQ.0) THEN
!    CALL CPU_TIME(start_time)
   WRITE(*,*)' '
   CALL Helac_mtime()    ! in Helac_Func_1
   WRITE(*,*)' '
   GLOBALINIT_average=0
   GLOBALINIT_checkgluons=0
   GLOBALINIT_id=0
   GLOBALINIT_id2=0
   GLOBALINIT_id3=0
   GLOBALINIT_redo=0
   CALL Helac_init()     ! in Helac_master
ENDIF
IF(initsig.EQ.0)THEN
	IF(repeat.EQ.1) THEN
!  OPEN(20,FILE='cmat_'//process(2:2*n+1)//'.in')
		OPEN(21,FILE=TRIM(input_dir)//'tree_'//TRIM(process(2:10*nhad+1))//'.in')
!  OPEN(22,FILE='kine_'//process(2:2*n+1)//'.in')
!  OPEN(20,FILE='cmat.in')
!  OPEN(21,FILE='tree.in')
!  OPEN(22,FILE='kine.in')
		WRITE(*,*)' '
		CALL Helac_mtime()
		WRITE(*,*)' '
		GLOBALINIT_average=0
		GLOBALINIT_checkgluons=0
		GLOBALINIT_id=0
		GLOBALINIT_id2=0
		GLOBALINIT_id3=0
		GLOBALINIT_redo=0
		CALL Helac_init()
		WRITE(*,*)' '
		CALL Helac_mtime()
		CLOSE(0)
		STOP
	ENDIF
	IF(repeat.EQ.2)THEN
!   (20,FILE='cmat_'//process(2:2*n+1)//'.in')
		OPEN(21,FILE=TRIM(input_dir)//'tree_'//TRIM(process(2:10*nhad+1))//'.in')
!  OPEN(22,FILE='kine_'//process(2:2*n+1)//'.in')
!  OPEN(20,FILE='cmat.in')
!  OPEN(21,FILE='tree.in')
!  OPEN(22,FILE='kine.in')

		CALL Helac_getlist()  ! in Helac_pan1
		READ(21,*)ng
		WRITE(*,*)'The number of Feynman graphs NG= ',ng
		DO i1=1,ng
			DO j1=1,8
				READ(21,*)(is(i1,j1,k1),k1=1,8)
			ENDDO
		ENDDO
	ENDIF      
! some initialization for histograms       

!IF(initsig.EQ.0)THEN
	OPEN(16,FILE=TRIM(tmp_dir)//TRIM(histo_file))
	OPEN(17,FILE=TRIM(tmp_dir)//TRIM(error_file))
	initsig=1
ENDIF
!      nh=5
!       include 'nh.h'
nh=0
!  end of 'nh.h'
WRITE(16,*)ifl(1:n)
DO i=1,nh
   CALL Helac_histo3(i)   ! in Helac_histo
ENDDO
CALL Helac_histo3(199)
! ---------
CALL ReadElem_integer("gener",gener)
IF(gener.NE.3.AND.gener.GE.0)THEN
	WRITE(*,*)' '
	CALL Helac_mtime()               ! timing
	WRITE(*,*)' '
	CALL Main_integration(w0)
	WRITE(*,*)' '
	CALL Helac_mtime()               ! timing
	CLOSE(nunit2)
	CLOSE(nunit3)
	CLOSE(21)
	CLOSE(200)
        IF(lunwei)CALL Generate_lhe(nhad,Nevents,process,icase)
ELSEIF(gener.LT.0)THEN
   ! read phase space point from PS.input
        OPEN(UNIT=619,FILE=TRIM(input_dir)//'PS.input',STATUS='OLD',ACTION='READ',IOSTAT=ioerror)
        IF(ioerror.NE.0)THEN
           STOP 'Could not read the PS.input phase-space point.'
        ENDIF
        j1 = 1
        DO i=1,nhad
           pmomtemp(1:4)=0d0
           READ(619,*) pmomtemp(4),pmomtemp(1),pmomtemp(2),pmomtemp(3)
           IF(ABS(iflh(i)).LE.100)THEN
              ! it is a parton
              Phegas_pmom(j1,1:4)=pmomtemp(1:4)
              Phegas_pmom(j1,5)=parmas(ifl(j1))
              j1=j1+1
           ELSE
              ! it is a hadron
              Phegas_pmom(j1,1:4)=pmomtemp(1:4)*parmas(ifl(j1))/(parmas(ifl(j1))+parmas(ifl(j1+1)))
              Phegas_pmom(j1,5)=parmas(ifl(j1))
              Phegas_pmom(j1+1,1:4)=pmomtemp(1:4)*parmas(ifl(j1+1))/(parmas(ifl(j1))+parmas(ifl(j1+1)))
              Phegas_pmom(j1+1,5)=parmas(ifl(j1+1))
              j1=j1+2
           ENDIF
        ENDDO
        CLOSE(UNIT=619,STATUS="KEEP")
        DO j1=1,n
           ptp=DSQRT(Phegas_pmom(j1,1)**2+Phegas_pmom(j1,2)**2)
           pqp=DSQRT(ptp**2+Phegas_pmom(j1,3)**2)
           zq(j1,1)= DCMPLX(Phegas_pmom(j1,4)+Phegas_pmom(j1,3),Phegas_pmom(j1,3))
           zq(j1,2)= DCMPLX(Phegas_pmom(j1,4)-Phegas_pmom(j1,3),ptp)
           zq(j1,3)= DCMPLX(Phegas_pmom(j1,1), Phegas_pmom(j1,2))
           zq(j1,4)= DCMPLX(Phegas_pmom(j1,1),-Phegas_pmom(j1,2))
           zq(j1,5)= DCMPLX(parmas(ifl(j1)),pqp )
        ENDDO
        CALL Helac_master_f(wme)
        GLOBALINIT_master=1
        WRITE(*,*)"result=",wme*LDMEwt
ELSE ! VEGAS
	CALL ReadElem_integer("itmax",itmxn)
	CALL ReadElem_integer("nmc",ncalln)
	WRITE(*,*)' '
	CALL Helac_mtime()
	WRITE(*,*)' '
        IF(nhad.EQ.5.AND..NOT.lunwei)THEN
           ! it is wrong for unweighted events since the first final particle's
           ! phi0 is not necessary to be 0
           IF(COLL_TYPE.NE.3)THEN
              CALL BB_3FS(rslt,itmxn,ncalln)
           ELSE
              ! there is still problem with EP_3FS when unwgt=T,ranhel=3 for e- e+ > z z z
              CALL EP_nFS(rslt,itmxn,ncalln)
           ENDIF
           WRITE(nunit1,*)"sigma(nb)                   sd"
           WRITE(nunit1,*)rslt(1),rslt(2)
	ELSEIF(nhad.EQ.4)THEN
           IF(COLL_TYPE.NE.3)THEN
              CALL BB_2FS(rslt,itmxn,ncalln)
           ELSE
              CALL EP_2FS(rslt,itmxn,ncalln)
           ENDIF
           WRITE(nunit1,*)"sigma(nb)                   sd"
           WRITE(nunit1,*)rslt(1),rslt(2)
        ELSEIF(nhad.EQ.3)THEN
           IF(COLL_TYPE.NE.3)THEN
              CALL BB_1FS(rslt,itmxn,ncalln)
           ELSE
              STOP "No integral left for 2 > 1 phase space with x1=x2=1"
           ENDIF
           WRITE(nunit1,*)"sigma(nb)                   sd"
           WRITE(nunit1,*)rslt(1),rslt(2)
        ELSEIF(nhad.GE.4)THEN
           IF(COLL_TYPE.NE.3)THEN
              CALL BB_nFS(rslt,itmxn,ncalln)
           ELSE
              CALL EP_nFS(rslt,itmxn,ncalln)
           ENDIF
           WRITE(nunit1,*)"sigma(nb)                   sd"
           WRITE(nunit1,*)rslt(1),rslt(2)
	ELSE
                PRINT *,"Wrong of nhad in Helac_Vegas ! STOP!"
		STOP
	ENDIF
	WRITE(*,*)' '
	CALL Helac_mtime()
        CLOSE(nunit2)
        CLOSE(nunit3)
        CLOSE(21)
        CLOSE(200)
        IF(lunwei)CALL Generate_lhe(nhad,Nevents,process,icase)
ENDIF
! ---------
! THE MAIN INTEGRATION ROUTINE: see below       
! ---------
!CALL ReadElem_logic('vegasQ',vegasQ)
!WRITE(*,*)' '
!CALL Helac_mtime()               ! timing
!WRITE(*,*)' '
!IF(vegasQ)THEN
!	CALL Collider_PSI(rslt,5,100000)
!	PRINT *,rslt(1:3)
!ELSE
!CALL Main_integration(w0)
!ENDIF
!WRITE(*,*)' '
!CALL Helac_mtime()               ! timing
       
! ---------
! output histograms
!DO i=1,nh
!   CALL Helac_histo2(i,16,0,w0)
!ENDDO
! ---------
!CLOSE(nunit2)
!CLOSE(nunit3)
!CLOSE(21)
!CLOSE(200)
!Generate .lhe file
!IF(.NOT.vegasQ)
!CALL Generate_lhe(nhad,process)
END SUBROUTINE SinProcess
      
SUBROUTINE Main_integration(w0)
USE Helac_Global
USE Helac_Feynman
USE Helac_Func_1
USE Phegas_mod
USE Cuts_Module
USE Helac_SM_FeynRule
USE Helac_ranmar_mod
USE Kinetic_Func
USE MC_Funcs
USE Structf_PDFs
USE KT_Clustering
USE Helac_unwei
USE MC_Helac_GRID
USE Helac_master
USE Helac_histo
USE Constants
USE Adapt
USE FO_plot
USE QEDPS_interface
USE Decay_interface
IMPLICIT NONE
REAL(KIND=DBL),INTENT(OUT)::w0
INTEGER,DIMENSION(-100:100)::ipdt
REAL(KIND=DBL),DIMENSION(ng+1)::alpha,wei,wpsp
REAL(KIND=DBL),DIMENSION(0:ng+1)::beta
REAL(KIND=DBL),DIMENSION(20,4)::p
INTEGER,DIMENSION(ng+1)::igraph
INTEGER,DIMENSION(0:ng+1)::icount0,icount1
REAL(KIND=DBL),DIMENSION(10)::wmax
INTEGER,DIMENSION(ng+1)::ig_m
REAL(KIND=DBL)::hlower,hupper,hvaria,hweigh,w1
CHARACTER(len=2)::file2
CHARACTER(len=3)::file3
REAL(KIND=DBL),DIMENSION(ng+1)::alpha_m,xrr,y
CHARACTER(len=4)::file4
CHARACTER(len=5)::file5
INTEGER::nmc,idup,nwri_tot,ncha,nchch,idprup,iij,l,j,k,icut,number_opt
INTEGER::nopti,nopt_step,maxopt,noptlim,iopt,nopt,imothup1,imothup2,icol1,icol2,iev,nm,&
         nch0,ipass,idir,iweight,ig,IDBMUP1,IDBMUP2,IDWTUP,NPRUP=0,nlim,nlimit,iend,istup,&
		 ndim_adapt,iv
SAVE NPRUP
REAL(KIND=DBL)::optf,xwgtup,px,py,pz,p0,pmass,xx1,yy1,x1,yaa1,wme,wmemax,w,wpsp_t,&
         ss0,alimit,vtime,vspin,sum,scale1,scalup,a,b,wr,wsf,wparni,wo,w2,wgtbr
!LOGICAL::onep       
! WRITE OUT
! START
LOGICAL::lunwei=.FALSE.,llexist,ladapt
! lunwei: for unweightinga   !1/3
!       data lunwei/.false./  !2/3
!      data lunwei/.true./   !3/3
INTEGER::nunwei=1000,flag1,flag2,flag4,flag5
LOGICAL::lwmax=.FALSE.
INTEGER::nwri,nwmax,nwarn
REAL(KIND=DBL)::umax,umax1
!INTEGER::initint=0
!SAVE initint
INTEGER::iunit=80
INTEGER,DIMENSION(ng+1)::ind
INTEGER::i,inonzero
INTEGER::iyr0,iyr1,iday0,iday1,imon0,imon1,ihr0,ihr1,imin0,imin1,isec0,isec1,i100th0,i100th1
umax=0
umax1=0     
nwri=0
nwmax=0
nwarn=0
IF(SUMMATIONQ)THEN
	GLOBALINIT_Phegas=0
	GLOBALINIT_Bookin=0
	GLOBALINIT_master=0
	GLOBALINIT_sethel=0
	GLOBALINIT_optimize=0
	GLOBALINIT_feyn=0
	GLOBALINIT_unwei=0
	GLOBALINIT_setmin=0
	GLOBALINIT_setmax=0
	GLOBALINIT_gen_x=0
	GLOBALINIT_adapt_gen=0
	GLOBALINIT_physicalpol=0
!	GLOBALINIT_findlimits=0
ENDIF
icount0(0:ng+1)=0
icount1(0:ng+1)=0

WRITE(*,*)' '
WRITE(*,*)'BEGINNING THE INTEGRATION SUBROUTINE '
WRITE(*,*)'----------------------------------------------------------------------------'     
WRITE(*,*)'The number of Feynman Graph is ng=',ng
!PAUSE
IF(ng.EQ.0)THEN
   WRITE(*,*)'no Feynman Graph contributions found: BYE'
   STOP
ENDIF
! --------------------------------------------------------------
!         START OF INITIALIZATION
! ---------------------------------------------------------------

!IF(initint.EQ.0)THEN
!    idebug=0
WRITE(*,*)'=========================================================================='
WRITE(*,*)'FOR DEBUGGING',idebug
CALL ReadElem_logic('onep',onep)
!    onep=.FALSE.
! Read adapt
CALL ReadElem_logic("adapt",ladapt)
IF(ladapt)THEN
	CALL ReadElem_integer("maxadaptnum",MaxAdapt)
ENDIF
IF(.NOT.onep)THEN     
! unweighting
!    WRITE(*,*)'Enter the unweighting lunwei(T(1) or F(0)):'
!    READ(*,*,END=1000)lunwei
	CALL ReadElem_logic('unwgt',lunwei)
!	WRITE(*,*)'Enter the unweighting number(begin to unweighting) nunwei='
!	READ(*,*,END=1000)nunwei
	CALL ReadElem_integer('preunw',nunwei)
!	WRITE(*,*)'Enter the total unweighting number(out of MC) nwri_tot='
!	READ(*,*,END=1000)nwri_tot
	CALL ReadElem_integer('unwevt',nwri_tot)
!     read*,lunwei,nunwei,nwri_tot
	WRITE(*,*)'UNWEIGHTING IS lunwei=', lunwei
	CALL ReadElem_logic("ptdisQ",dSigmadPTQ)
	IF(dSigmadPTQ)THEN
		IF(ABS(iflh(3)).LE.100)THEN
			PRINT *,"Main_integration failed to generate the parton differential cross section !"
			STOP
		ENDIF
                IF(emep_ISR_shower.NE.0)STOP "PT distribution is not avialabel when there is ISR shower"
		CALL ReadElem_r("Pt1",PTFirst)
	ENDIF
	beta(0)=0
	igraph(1:ng+1)=0
	inonzero=0
!         read*,ncha
	WRITE(*,*)' '
!    WRITE(*,*)'Enter the types of kinematical chennels'
	WRITE(*,*)'ncha<=1 use RAMBO,ncha=2 use DURHAM'
	WRITE(*,*)'ncha=0/-1,use multichannel(ng+1 or user define),ncha>=1 without multichannel'
!    READ(*,*,END=1000)ncha
	CALL ReadElem_integer('gener',ncha)
	WRITE(*,*)'option for kinematical channels:',ncha
	IF(ncha.EQ.0)THEN
		IF(.NOT.dSigmadPTQ)THEN
			nchch=ng+1
		ELSE
			nchch=ng
			IF((n-2).EQ.2*(nhad-2))nchch=1
		ENDIF
		DO i=1,nchch
			igraph(i)=MOD(i,ng+1)
		ENDDO
	ELSEIF(ncha.GE.1)THEN
		IF(dSigmadPTQ)THEN
			PRINT *,"Wrong of using RAMBO/Durham to calculate dsigma/dpt ! STOP!"
			STOP
		ENDIF
		nchch=1
		igraph(1)=0
	ELSEIF(ncha.EQ.-1)THEN
		WRITE(*,*)' You are using user defined multichannel. Enter the number of channels'
		READ(*,*)nchch
		DO i=1,nchch
			igraph(i)=i
		ENDDO
	ENDIF
	IF(ncha.LE.1)i_psp=1
	IF(ncha.EQ.2)i_psp=2

	WRITE(*,*)' '
	WRITE(*,*)'Number of channels ',nchch
	CALL ReadElem_integer('nmc',nmc)
	WRITE(*,*)'Number of MC points ',nmc
!    READ(*,*)nmc
!    WRITE(*,*)' nmc  = ',nmc
	alpha(1:nchch)=dnou(1)/nchch
	alpha_m(1:nchch)=alpha(1:nchch)
	DO i=1,nchch
		beta(i)=beta(i-1)+alpha(i)
	ENDDO
	WRITE(*,*)' '
	WRITE(*,*)'The  optimization options:'
	WRITE(*,*)'O1=n1,O2=O1+n2*n3,..,On4=O(n4-1)+n2*n3^(n4-1)'
	WRITE(*,*)'n5 is the maximum of MC iterations'
	WRITE(*,*)'n6 is the flag: 0 not perform, 1 perform'
	CALL ReadElem_integer('nopt',nopti)
	WRITE(*,*)'n1=',nopti
!    READ(*,*)nopti
	CALL ReadElem_integer('nopt_step',nopt_step)
	WRITE(*,*)'n2=',nopt_step
!	READ(*,*)nopt_step
	CALL ReadElem_real('optf',optf)
	WRITE(*,*)'n3=',optf
!	READ(*,*)optf
	CALL ReadElem_integer('maxopt',maxopt)
	WRITE(*,*)'n4=',maxopt
!	READ(*,*)maxopt
	CALL ReadElem_integer('noptlim',noptlim)
	WRITE(*,*)'n5=',noptlim
!	READ(*,*)noptlim
	CALL ReadElem_integer('iopt',iopt)
	WRITE(*,*)'n6=',iopt
!	READ(*,*)iopt
	nopt=100
	number_opt=0
	wmax(1:10)=0
	wmemax=0
	wme=1
	IF(nchch.EQ.1)THEN
		iopt=0
		maxopt=-1
	ENDIF
	WRITE(*,*)'nopt,nopt_step(n2),optf(n3),maxopt(n4),iopt(n6)'
	WRITE(*,*) nopt,nopt_step,optf,maxopt,iopt
	WRITE(*,*)'number of channels=',nchch
	WRITE(*,*)'==================================================================='
	WRITE(*,*)'  '
	imc=0
!    IF(iopt.EQ.0)imc=0
	CALL ReadElem_logic('ktrw',ktrw)
	IF(ktrw)THEN
		CALL ReadElem_integer("ktmeasure",ktmeasure)
		CALL ReadElem_real("Rparameter",RR)
	ENDIF
        IF(plot_output)CALL initplot
ELSE
		
	PRINT *,' We use the external phase space points mode '
	PRINT *,'------------------------------------------------------------------'
	WRITE(*,*)' '
	CALL ReadElem_integer('nmc',nmc)
	WRITE(*,*)'Number of MC points ',nmc
	WRITE(*,*)' '
	WRITE(*,*)'Name of the momenta file is ',TRIM(momfile)
    WRITE(*,*)'Name of the output file is ',TRIM(momout)

    OPEN(UNIT=90,FILE=TRIM(tmp_dir)//TRIM(momfile))
    OPEN(UNIT=100,FILE=TRIM(tmp_dir)//TRIM(momout))

	igraph(1:ng+1)=0
	wme=1
	CALL ReadElem_logic('ktrw',ktrw)
	IF(ktrw)THEN
		CALL ReadElem_integer("ktmeasure",ktmeasure)
		CALL ReadElem_real("Rparameter",RR)
	ENDIF
ENDIF
ntotps=0
!    initint=1
!ENDIF

! --------------------------------------------------------------
!         END  OF INITIALIZATION
! ---------------------------------------------------------------

! 4     continue
!CALL GETTIM(ihr0,imin0,isec0,i100th0)
CALL DAYTIME(iyr0,imon0,iday0,ihr0,imin0,isec0,i100th0)
oneponep:IF(.NOT.onep)THEN
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! onep=.FALSE.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
loop1:DO
flag4=0
flag5=0
WRITE(*,*)'NUMBER OF CHANNELS',nchch
ipass=0
!      do iev=1,nmc   ![iev
iev=0
DO WHILE(iev.LE.nmc)  ![iev
        IF(measurespeed.AND.iev.EQ.500)THEN
           CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
           iyr1=iyr1-iyr0
           imon1=imon1-imon0
           iday1=iday1-iday0
           ihr1=ihr1-ihr0
           imin1=imin1-imin0
           isec1=isec1-isec0
           i100th1=i100th1-i100th0
           CALL Vegas_speed(500,iday1,ihr1,imin1,isec1,i100th1)
        ENDIF
	iev=iev+1
	w=0
	wme=0
	wei(1:nchch)=0
	CALL Helac_igen(nchch,beta(1:nchch),ig)
	icount0(igraph(ig))=icount0(igraph(ig))+1
	idir=0
	CALL phegas(wpsp(ig),igraph(ig),idir,iweight)
	GLOBALINIT_Phegas=1
!	WRITE(*,*)' '
!	WRITE(*,*)' Continue?'
!    WRITE(*,*)wpsp(ig)
!	READ(*,*)temp
!	WRITE(*,*)wpsp(ig)  ! Debug
!     call test_gen(wpsp(ig),igraph(ig),0,iweight)
	level21:IF(iweight.NE.0)THEN      !goto 2
		flag2=0
		level1:IF(iweight.NE.-1)THEN     !goto 1
       
			level2:IF(wpsp(ig).GT.0)THEN
				IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
					CALL Cuts_HC(icut)
				ELSEIF(COLL_TYPE.EQ.3)THEN
                                        IF(emep_ISR_shower.EQ.0)THEN
                                           CALL Cuts_epem(icut)
                                        ELSE ! generate momenta with QEDPS after showering
                                           CALL QEDPSEVNT(icut)
                                        ENDIF
				ELSE
				    icut=1
				ENDIF
                                wgtbr=1d0
                                IF(icut.NE.0.AND.MCoHelicity)THEN
                                   RMCoH=Helac_rnmy(0)
                                   DO k=1,NDecayIflh
                                      Decayran(k)=Helac_rnmy(0)
                                   ENDDO
                                   !IF(NDecayChains.GT.0)THEN
                                   !   CALL HO_Decay(wgtbr)
                                   !   CALL cuts_Decay(icut)
                                   !ENDIF                                   
                                ENDIF
				level3:IF(icut.NE.0)THEN   !goto 1
! def include cuts from possible decays
!      call decays(0,icut,w)
!      if(icut.eq.0)goto1

					wpsp_t=0
					flag1=0
					flag2=0
					DO i=1,nchch
						IF(i.NE.ig)THEN
							idir=1
							CALL phegas(wpsp(i),igraph(i),idir,iweight)
							flag2=0
							IF(iweight.EQ.0)THEN      !goto 2
								flag2=1
								EXIT
							ENDIF
							flag1=0
							IF(iweight.EQ.-1)THEN !goto 1
								flag1=1
								EXIT
							ENDIF
							IF(wpsp(i).GT.0)THEN
								wpsp_t=wpsp_t+alpha(i)/wpsp(i)
							ELSE 
								WRITE(*,*)'wpsp2=',wpsp(i),i,wpsp(ig),ig,nchch,igraph(ig)
							ENDIF
						ELSE
							wpsp_t=wpsp_t+alpha(i)/wpsp(i)
						ENDIF
					ENDDO

					level4:IF(flag1.EQ.0.AND.flag2.EQ.0)THEN
						inonzero=inonzero+1
						icount1(igraph(ig))=icount1(igraph(ig))+1
						wpsp_t=wpsp_t**(-1)
						IF(irun.EQ.1)CALL Helac_FeynRule_SM()

						wsf=1
						IF(istruc.EQ.1)CALL strf_pdf(wsf)

						CALL helac_master_f(wme)
						GLOBALINIT_master=1
                                                
						w=wpsp_t*wme*wrest
						w=w*wsf*LDMEwt*weight_br
						ipass=ipass+1

						wr=1
!    include 'adapt1.h'
						ndim_adapt=2
						IF(istruc.EQ.0.OR.(.NOT.ladapt))ndim_adapt=0
						DO iv=1,ndim_adapt
							CALL adapt_weo(wo,iv)
							w=w*wo
						ENDDO
! end of 'adapt1.h'
!    include 'ktreweight.h'
						IF(ktrw.AND.emep_ISR_shower.EQ.0)THEN
							CALL ktreweight(wr)
							w1=w*wr
							w=w1
						END IF
!   end of 'ktreweight.h'

!      if(.not.lwmax)then
!      include 'myparni.h'
!      endif

						w1=w
!                        WRITE(*,*)w1  ! Debug
!  						STOP ! Debug
!def include cuts just for histograms
!      call decays(1,icut,w)
						IF(number_opt.LE.maxopt)THEN
							wei(1:nchch)=0
							DO i=1,nchch
								IF(wpsp(i).GT.0) wei(i)=w*w*wpsp_t/wpsp(i)
							ENDDO
						ENDIF

						IF(wme.GT.wmemax)THEN
							wmemax=wme
						ENDIF
						IF(w.GT.wmax(1))THEN
							DO i=10,2,-1
								wmax(i)=wmax(i-1)
							ENDDO
							IF(wmax(1).GT.0)THEN
								IF(w/wmax(1).GT.10)THEN
									IF(idebug.EQ.1)THEN
										WRITE(nunit2,*)'WARNING ABOUT WEIGHT'
										iij=igraph(ig)
										WRITE(nunit2,*)iev,igraph(ig),alpha(ig)
										WRITE(nunit2,*)w,wme,wpsp_t
										IF(iij.GT.0)THEN
											k=1
											DO WHILE(is(iij,k,1).NE.0)
												WRITE(nunit2,*)is(iij,k,1:8)
												k=k+1
											ENDDO
										ENDIF
										WRITE(nunit2,*)w,wmax
										WRITE(nunit2,*)(i,Phegas_pmom(i,1:4),i=1,n)
										WRITE(nunit2,*)'--------------------'
									ENDIF
								ENDIF
							ENDIF
							wmax(1)=w
						ENDIF

! HISTOGRAMING

!          call tohisto(w)
!      include 'toh.h'

! WRITE OUT EVENTS
! START
! -------------
! only unweight events are generated -> LHE files
						IF(lunwei)THEN
! --------------
							IF(.NOT.lwmax.AND.(number_opt.GT.maxopt.OR.iopt.EQ.0))THEN
								IF(umax.LT.w1)THEN
									umax=w1
								ENDIF
								nwmax=nwmax+1
								umax1=umax
								CALL Helac_put_unwei(w1)
								GLOBALINIT_unwei=1
							ENDIF
!            if(nwmax.eq.nopt_step) then
							IF(nwmax.EQ.nunwei)THEN
								nwmax=nwmax+1
								lwmax=.TRUE.
								CALL Helac_get_unwei(umax)
								WRITE(nunit1,*)'START UNWEIGHTING',iev,umax,nwmax,nunwei
							ENDIF
							IF(lwmax)THEN
								IF(w1.GT.umax1)THEN
									umax1=w1
									WRITE(nunit1,*)'iev,umax1,umax',iev,umax1,umax
								ENDIF
								lwri=.FALSE.
								IF(umax*Helac_rnmy(0).LT.w1)lwri=.TRUE.
								IF(umax.LT.w1)nwarn=nwarn+1
								a=0.d0
								b=10*umax
								CALL Helac_histo1(199,100,a,b,w1,1.d0)
								

                                                                IF(lwri)THEN                   
									nwri=nwri+1
                                                                        IF(NDecayChains.GT.0)THEN
                                                                           CALL unwei_writer_Decay
                                                                        ELSEIF(emep_ISR_shower.EQ.1)THEN
                                                                           CALL unwei_writer_QEDPS
                                                                        ELSE
                                                                           idprup=81 ! id for the process
                                                                           !IF(lunwei.AND..NOT.lhewgtup)THEN
                                                                           xwgtup=1 !w1*10**3 !1
                                                                           !ELSE
                                                                           !xwgtup=w1*10**3 ! in unit of pb for one event
                                                                           !ENDIF
                                                                           CALL qcdscale(scale1)
                                                                           scalup=scale1
                                                                           WRITE(nunit3)nhad,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
                                                                           IF(imode.EQ.0.AND.NDecayChains.EQ.0)CALL cmstolab()
                                                                           
                                                                           DO i=1,nhad
                                                                              !	idup=ifl(i)
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
											

                                                                              px=Phegas_hmom(i,1)
                                                                              py=Phegas_hmom(i,2)
                                                                              pz=Phegas_hmom(i,3)
                                                                              p0=Phegas_hmom(i,4)
                                                                              pmass=Phegas_hmom(i,5)
                                                                              vtime=0
                                                                              vspin=9
                                                                              
                                                                              WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2&
                                                                                   ,px,py,pz,p0,pmass,vtime,vspin

                                                                           ENDDO
                                                                           ! boost back
                                                                           IF(imode.EQ.0.AND.NDecayChains.EQ.0)CALL labtocms()
                                                                        ENDIF
                                                               ENDIF
							ENDIF
!     rewind(nunit3)

							IF(nwri.EQ.nwri_tot)THEN !go to 5
								flag5=1
								EXIT
							ENDIF

! --------------
						ENDIF
					ENDIF level4
				ENDIF level3
			ENDIF level2
		ENDIF level1
! --------------
! END
! 1        continue  
		IF(flag2.EQ.0)THEN
		!       include 'adapt2.h'
			DO iv=1,ndim_adapt
				w2=w**2
				CALL adapt_wei(w2,iv,ndim_adapt,.FALSE.)
				IF(iopt.EQ.0.AND.MOD(inonzero,1000).EQ.0)CALL adapt_wei(w2,iv,ndim_adapt,.TRUE.)
			ENDDO
		! end of 'adapt2.h'
			IF(.NOT.lwmax.AND.COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
!       include 'myparni.h'  
				IF(.NOT.ladapt)THEN 
!    'myparni_yes.h'
					wparni=w
					IF(w.LT.0)wparni=-w
					CALL helac_grid_collect_2( wparni )
				ENDIF
!   end of 'myparni.h'
			ENDIF

			CALL bookin(1,w)  ! in MC_Funcs.f90
			GLOBALINIT_bookin=1
			CALL bookin(2,wme)
!         if(w.gt.0)call bookin(3,dlog(w))
			IF(lwri)CALL bookin(3,w)
			IF(number_opt.LE.maxopt)THEN
				DO i=1,nchch
					CALL bookin(i+3,wei(i))
				ENDDO
			ENDIF
                        ntotps=ntotps+1
                        IF(plot_output)CALL outfun(w)
		ENDIF
	ENDIF level21 
!	WRITE(*,*)iev ! Debug
! 2        continue   
	CALL geti(1,0,w0)
!      print*,'w0=',w0
	IF(MOD(iev,1000).EQ.0)THEN
		CALL errest(1,x1,yaa1,0)
		WRITE(*,'(a6,2d15.6,3i10)')'sigma=',x1,SQRT(yaa1)/x1,inonzero,INT(w0),iev
		IF(nunit1.NE.6)THEN
			WRITE(nunit1,'(a6,2d15.6,3i10)')'sigma=',x1,SQRT(yaa1)/x1,inonzero,INT(w0),iev
		ENDIF
		CALL errest(3,xx1,yy1,0)
		IF(yy1.GT.0)THEN
			WRITE(*,'(a7,2d15.6)')'<w=/=0>',xx1,SQRT(yy1)/DABS(xx1)
			IF(nunit1.NE.6)THEN
				WRITE(nunit1,'(a7,2d15.6)')'<w=/=0>',xx1,SQRT(yy1)/DABS(xx1)
			ENDIF
		ENDIF
		IF(nwri.GT.0)WRITE(*,*)nwri
		WRITE(17,*)iev,x1,SQRT(yaa1)
		WRITE(*,*)'------------------------------------'
		IF(nunit1.NE.6)THEN
			WRITE(nunit1,*)'------------------------------------'
		ENDIF
	ENDIF
	IF(iev.EQ.nopt)THEN
		IF(ipass.EQ.0)nopt=nopt*2
	ENDIF
! optimization   [
!         print*,'nopt,iopt,number_opt,maxopt,nch'
!         print*, nopt,iopt,number_opt,maxopt,nch 
! optimization is performed when the program reaches the MC iterations
! {O1=nopt,O2=O1+nopt_step*optf,...O(maxopt)=O(maxopt-1)+nopt_step*optf^(maxopt-1)}
! It stops when the maximum number of optimizations maxopt is reached,or when the noptlim-th
! MC iteration is reached. iopt is a flag:if set to 0 optimisation is not performed,
! if set to 1 it is performed.
	flag4=0
!    IF(iev.GE.100)WRITE(*,*)'haha1'  ! Debug
	IF(iev.EQ.nopt.AND.iopt.EQ.1.AND.number_opt.LE.maxopt.AND.nchch.GT.1)THEN
		sum=0
		DO i=1,nchch 
			CALL errest(i+3,xrr(i),y(i),0)
			sum=sum+alpha(i)*xrr(i)
!     x(i)=x(i)+2*dsqrt(y(i))
		ENDDO
!       write(1000,*)number_opt+1,sum
! print*,nch,igraph(1:ng)
		nch0=nchch
		CALL Helac_optimize(nchch,igraph(1:nchch),alpha,nm,ig_m,alpha_m,xrr,iend)
		GLOBALINIT_optimize=1
		IF(nchch.EQ.1)iopt=0
		IF(iend.EQ.2)number_opt=maxopt
		IF(number_opt.GE.1)THEN

!       include 'adapt3.h'
			DO iv=1,ndim_adapt
				w2=w**2
				CALL adapt_wei(w2,iv,ndim_adapt,.TRUE.)
			ENDDO
!       end of 'adapt3.h'
!               data iunit/80/
			DO
				INQUIRE(UNIT=iunit,OPENED=llexist)
				IF(llexist)THEN
					iunit=iunit+1
				ELSE
					EXIT
				ENDIF
			ENDDO
			WRITE(nunit1,*)'UNIT=',iunit
			IF(iunit.GE.10.AND.iunit.LE.99)THEN
				WRITE(file2,'(I2)') iunit
				OPEN(iunit,FILE=TRIM(tmp_dir)//'data.'//file2)
			ELSEIF(iunit.GE.100.AND.iunit.LE.999)THEN
				WRITE(file3,'(I3)') iunit
				OPEN(iunit,FILE=TRIM(tmp_dir)//'data.'//file3)
			ELSEIF(iunit.GE.1000.AND.iunit.LE.9999)THEN
				WRITE(file4,'(I4)') iunit
				OPEN(iunit,FILE=TRIM(tmp_dir)//'data.'//file4)
			ELSEIF(iunit.GE.10000.AND.iunit.LE.99999)THEN
				WRITE(file5,'(I5)') iunit
				OPEN(iunit,FILE=TRIM(tmp_dir)//'data.'//file5)
			ELSE
				WRITE(*,*)'The iunit>=100000 in integration ',iunit,' ! STOP !'
				STOP
			ENDIF

			WRITE(iunit,'(i4)')nchch
			DO i=1,nchch
				WRITE(iunit,'(3(e10.4,5x),i8)')alpha(i),xrr(i),DSQRT(y(i)),igraph(i)
			ENDDO
			CLOSE(iunit)
			iunit=iunit+1
		ENDIF
		IF(number_opt.EQ.0)THEN
			nopt=nopti
		ELSE
			nopt=nopt+nopt_step
		ENDIF
		WRITE(nunit1,*)'----------------------------'
		WRITE(nunit1,*)'iev,nopt,number_opt,maxopt,nch,nch0 ',&
						iev,nopt,number_opt,maxopt,nchch,nch0
!       print*,'<w>/w_max,x1,wmax',x1/wmax(1),x1,wmax
		IF(nopt.GT.noptlim)number_opt=maxopt
		nopt_step=nopt_step*optf
		WRITE(nunit1,*)'nopt_step',nopt_step
		number_opt=number_opt+1
		IF(number_opt.GT.maxopt)THEN
			nchch=nm
			igraph(1:nchch)=ig_m(1:nchch)
			WRITE(*,*)'# of graphs , # of channels',ng,nchch
!               print*,(i,igraph(i),i=1,nch)
			DO i=1,nchch
				alpha(i)=alpha_m(i)
				IF(10*alpha(i).GT.1)THEN
					k=1
					iij=igraph(i)
					WRITE(nunit1,*)'the graph',iij,alpha(i)
					IF(iij.EQ.0)CYCLE
					DO WHILE(is(iij,k,1).NE.0)
						WRITE(nunit1,*)is(iij,k,1:8)
						k=k+1
					ENDDO
				ENDIF
			ENDDO
			OPEN(iunit,FILE=TRIM(tmp_dir)//'data.final')
			WRITE(iunit,'(i4)')nchch
			DO i=1,nchch
				WRITE(iunit,'(e10.4,i8)')alpha(i),igraph(i)
			ENDDO
			CLOSE(iunit)
			iunit=iunit+1
		ENDIF
! optimization   ]
		DO i=1,nchch
			beta(i)=beta(i-1)+alpha(i)
		ENDDO
!            print*,'alpha',alpha(1:nch)
!        print*,'----------------------------'
!        print*,nch,beta(nch)
!        eps1= 0 !dnou(1)/100
!        do l=1,nch
!        if(alpha(l).gt.eps1)
!    .	 print'(i4,2d16.8)',igraph(l),alpha(l),beta(l) 
!        enddo
!        print*,'----------------------------'

		IF(number_opt.GT.maxopt.AND.nchch.GT.155)THEN
			CALL SOR(alpha,nchch,ind,1)
			sum=0
			i=0
			DO WHILE(sum.LT.0.90d0)
				i=i+1
				sum=sum+alpha(i)
			ENDDO
			nchch=i
			IF(nchch.GT.155) nchch=155
			ig_m(1:nchch)=igraph(1:nchch)
			DO i=1,nchch
! BUG.22.11.07   igraph(i)=mod(ind(i),ng+1)
				igraph(i)=ig_m(ind(i))
!           WRITE(*,*)i,ind(i),igraph(i)
			ENDDO
			sum=0
			DO i=1,nchch
				sum=sum+alpha(i)
			ENDDO

			WRITE(*,*)'ALPHA(I)',sum
			WRITE(*,'(d12.4)')alpha(1:nchch)
			WRITE(*,*)ind(1:nchch)
			DO i=1,nchch
				WRITE(*,*)igraph(i)
				IF(igraph(i).GT.0)&
				WRITE(*,'(8I6)')(is(igraph(i),j,1:8),j=1,n-2)
			ENDDO

			alpha(1:nchch)=alpha(1:nchch)/sum
			wmax(1:10)=0
			imc=0
			DO i=1,nchch
				beta(i)=beta(i-1)+alpha(i)
			ENDDO
		ENDIF

		IF(number_opt.EQ.1.OR.number_opt.GT.maxopt)THEN
			WRITE(nunit1,*)'clear ',number_opt,' opt'
			DO l=1,nch0
				CALL clear(l+2)
			ENDDO
			inonzero=0
			CALL clear(1)
			CALL clear(2)
!			goto 4
			flag4=1
		ENDIF

! 3           continue       
	ENDIF
	IF(flag4.EQ.1)EXIT
ENDDO ![iev
!       enddo ![idouble
IF(flag4.NE.1.OR.flag5.EQ.1)EXIT
ENDDO  loop1
!      do i=1,nch 
!       call errest(i+2,x(i),y(i),0)
!       print*,'channel',nch,i,igraph(i),x(i),dsqrt(y(i))
!      enddo

! 5      continue
CALL geti(1,0,w0)  ! in MC_Funcs
WRITE(*,*)'out of ',nmc,'  ',INT(w0),' points have been used'
WRITE(*,*)'and  ',inonzero,' points resulted to =/= 0 weight'
WRITE(*,*)'whereas  ',INT(w0)-inonzero,' points to 0 weight'
CALL errest(1,x1,yaa1,1)
WRITE(*,*)'------------------------------------'
WRITE(*,'(a19,1d15.6,a6,1d15.6)')' total sigma (nb) = ',x1,'+/-',SQRT(yaa1)
WRITE(*,*)'------------------------------------'
IF(x1.GT.0)WRITE(*,'(a9,1d15.6,a6)')' % error:',DSQRT(yaa1)/x1*100
WRITE(*,*)'------------------------------------'
IF(nunit1.NE.6)THEN
	WRITE(nunit1,*)'out of ',nmc,'  ',INT(w0),' points have been used'
	WRITE(nunit1,*)'and  ',inonzero,' points resulted to =/= 0 weight'
	WRITE(nunit1,*)'whereas  ',INT(w0)-inonzero,' points to 0 weight'
	WRITE(nunit1,*)'------------------------------------'
	WRITE(nunit1,'(a19,1d15.6,a6,1d15.6)')' total sigma (nb) = ',x1,'+/-',SQRT(yaa1)
	WRITE(nunit1,*)'------------------------------------'
	IF(x1.GT.0)WRITE(nunit1,'(a9,1d15.6,a6)')' % error:',DSQRT(yaa1)/x1*100
	WRITE(nunit1,*)'------------------------------------'	
ENDIF
! plot files
IF(plot_output)CALL plotout
!       if (lunwei) then
!          write(nunit0)n
!          write(nunit0)ifl(1:n)          
!          write(nunit0)x1,dsqrt(y1),wmax(1)
!       endif

CALL geti(3,0,ss0)
WRITE(*,*)'lwri: number of points have used',ss0
IF(nunit1.NE.6)THEN
	WRITE(nunit1,*)'lwri: number of points have used',ss0
ENDIF
!CALL GETTIM(ihr1,imin1,isec1,i100th1)
CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
ihr1=ihr1-ihr0
imin1=imin1-imin0
isec1=isec1-isec0
i100th1=i100th1-i100th0
CALL Helac_OTime(ihr1,imin1,isec1,i100th1)
WRITE(*,*)'The Timing Consuming in Phase Space Integration is:'
WRITE(*,*)ihr1,' h;',imin1,' m;',isec1,' s;',i100th1,' centi s'
IF(nunit1.NE.6)THEN
	WRITE(nunit1,*)'The Timing Consuming in Phase Space Integration is:'
	WRITE(nunit1,*)ihr1,' h;',imin1,' m;',isec1,' s;',i100th1,' centi s'
ENDIF

IF(ncha.NE.0.OR.number_opt.EQ.0)THEN
!	WRITE(*,*)'alimit,nlimit'
!	READ(*,*)alimit,nlimit
   CALL ReadElem_real('alimit',alimit)
   WRITE(*,*)'The lower limit of the alpha(i) alimit=',alimit
!   READ(*,*)alimit
   CALL ReadElem_integer('nlimit',nlim)
   WRITE(*,*)'The lower limit of the number of channels nlim=',nlim
ENDIF

!WRITE(*,*)' '
!WRITE(*,*)'Beam 1/2 identifies (e.g. 2212 for proton) IDBMUP1, IDBMUP2= ' 
!READ(*,*)IDBMUP1,IDBMUP2
!WRITE(*,*)'The energy in GeV of Beam 1/2 EBMUP1, EBMUP2= '
!READ(*,*)EBMUP1,EBMUP2
!WRITE(*,*)'Master switch dictating how the event weights are interpreted IDWTUP='   ! see hep-ph/0109068
!READ(*,*)IDWTUP
!WRITE(*,*)' The number of different user subprocesses NPRUP='
!READ(*,*)NPRUP
!WRITE(*,*)' '
!WRITE(*,*)IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,IDWTUP,NPRUP
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
IF(lunwei)THEN
	IDWTUP=3
ELSE
	IDWTUP=1
ENDIF
NPRUP=NPRUP+1
WRITE(200,5100) IDBMUP1,IDBMUP2,EBMUP1,EBMUP2,&      
                 iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2,IDWTUP,NPRUP
!    &             0,       0,   10042,   10042,IDWTUP,NPRUP
! XSECUP,XERRUP,XMAXUP,LPRUP
WRITE(200,5200) x1*10d0**3,DSQRT(yaa1)*10d0**3,1d0, 81

WRITE(16,*)x1,SQRT(yaa1)
IF(wmax(1).GT.0)WRITE(*,*)'<w>/w_max,w_max',x1/wmax(1),wmax(1)
IF(wmax(10).GT.0)WRITE(*,*)'<w>/w_max,w_max',x1/wmax(10),wmax(10)
WRITE(*,'(10(E14.6,2X))')wmax(1:10)
CALL errest(2,x1,yaa1,0)
IF(wmemax.GT.0)WRITE(*,*)'<me>/memax,memax',x1/wmemax,wmemax

WRITE(*,*)' '     
DO j=1,25
	IF(iwarning(j).GT.0)WRITE(*,*)'iwarning(',j,') = ',iwarning(j)
	IF(iwonders(j).GT.0)WRITE(*,*)'iwonders(',j,') = ',iwonders(j)
ENDDO
WRITE(*,*)' '
WRITE(*,*)'number of w=1 events',nwri
WRITE(*,*)'number of w>1 events',nwarn
WRITE(*,*)'maximum weight used for un:vs',umax,umax1
Nevents=nwri
CALL Helac_histo2(199,199,1,1d0)
CALL Helac_grid_plot_2( 998 )
IF(emep_ISR_shower.EQ.1)CALL QEDPSEND
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  onep=.TRUE. the user supply the momenta for the particles
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ELSE oneponep

! 4     continue in onep
ipass=0
iev=0
DO WHILE(iev.LE.nmc)  ![iev
	iev=iev+1
	w=0
	wme=0

	DO i=1,n
		READ(90,*,END=1000)phegas_pmom(i,1:5),wpsp_t
		WRITE(*,*)' The momenta configurations'
        PRINT *,i,phegas_pmom(i,1:5),wpsp_t
		IF(nunit1.NE.6)THEN
			WRITE(nunit1,*)' The momenta configurations'
			WRITE(nunit1,*)i,phegas_pmom(i,1:5),wpsp_t	
		ENDIF
	ENDDO

    CALL trans_p_to_zq(n,phegas_pmom(1:n,1:5),zq(1:n,1:5))

	IF(irun.EQ.1)CALL Helac_FeynRule_SM()

    CALL helac_master_f(wme)
	GLOBALINIT_master=1

	w=wpsp_t*wme

	ipass=ipass+1

	wr=1
!    include 'adapt1.h'
	ndim_adapt=2
	IF(istruc.EQ.0.OR.(.NOT.ladapt))ndim_adapt=0
	DO iv=1,ndim_adapt
		CALL adapt_weo(wo,iv)
		w=w*wo
	ENDDO
!  end of 'adapt1.h'

!    include 'ktreweight.h'
	IF(ktrw)THEN
		CALL ktreweight(wr)
		w1=w*wr
		w=w1
	END IF
!   end of 'ktreweight.h'

	w1=w

! HISTOGRAMING

!          call tohisto(w)
!      include 'toh.h'

! 1        continue in onep  

	CALL bookin(1,w) ! in MC_Funcs.f90
	GLOBALINIT_bookin=1
	CALL bookin(2,wme)
!         if(w.gt.0)call bookin(3,dlog(w))
! 2        continue  in onep 
	CALL geti(1,0,w0)
    IF(MOD(iev,1000).EQ.0)THEN
		CALL errest(1,x1,yaa1,0)
		PRINT '(A6,2D15.6,3I10)','sigma=',x1,SQRT(yaa1)/x1,INT(w0),iev
		IF(nunit1.NE.6)THEN
			WRITE(nunit1,'(A6,2D15.6,3I10)')'sigma=',x1,SQRT(yaa1)/x1,INT(w0),iev	
		ENDIF
        CALL errest(3,xx1,yy1,0)
        IF(yy1.GT.0)THEN
			PRINT '(A7,2D15.6)','<w=/=0>',xx1,SQRT(yy1)/DABS(xx1)
			IF(nunit1.NE.6)THEN
				WRITE(nunit1,'(A7,2D15.6)')'<w=/=0>',xx1,SQRT(yy1)/DABS(xx1)
			ENDIF
		ENDIF
        WRITE(17,*)iev,x1,SQRT(yaa1)
        PRINT *,'----------------------------'
		IF(nunit1.NE.6)THEN
			WRITE(nunit1,*)'----------------------------'
		ENDIF
	ENDIF

! 3           continue in onep       
ENDDO ![iev
   
! 5      continue in onep
CALL geti(1,0,w0)
PRINT *,'out of ',nmc,'  ',INT(w0),' points have been used'
PRINT *,'whereas  ',INT(w0),' points to 0 weight'
IF(nunit1.NE.6)THEN
	WRITE(nunit1,*)'out of ',nmc,'  ',INT(w0),' points have been used'
	WRITE(nunit1,*)'whereas  ',INT(w0),' points to 0 weight'
ENDIF
CALL errest(1,x1,yaa1,1)

CALL geti(3,0,ss0)

WRITE(16,*)x1,SQRT(yaa1)
CALL errest(2,x1,yaa1,0)

!CALL GETTIM(ihr1,imin1,isec1,i100th1)
CALL DAYTIME(iyr1,imon1,iday1,ihr1,imin1,isec1,i100th1)
ihr1=ihr1-ihr0
imin1=imin1-imin0
isec1=isec1-isec0
i100th1=i100th1-i100th0
CALL Helac_OTime(ihr1,imin1,isec1,i100th1)
WRITE(*,*)'The Timing Consuming in Phase Space Integration is:'
WRITE(*,*)ihr1,' h;',imin1,' m;',isec1,' s;',i100th1,' centi s'
IF(nunit1.NE.6)THEN
	WRITE(nunit1,*)'The Timing Consuming in Phase Space Integration is:'
	WRITE(nunit1,*)ihr1,' h;',imin1,' m;',isec1,' s;',i100th1,' centi s'
ENDIF

ENDIF  oneponep

5100 FORMAT(1P,2I8,2E14.6,6I8)
5200 FORMAT(1P,3E20.10,I6)
!     call test_end

1000   RETURN
END SUBROUTINE Main_integration

SUBROUTINE Helac_optimize(n1,igraph,a,nm,ig_m,am,x,iend)
!       include 'declare.h'
!       include 'common_feyn.h'
!       include 'common_debug.h'
USE Helac_Global
USE Helac_Func_1
USE Helac_Feynman
USE Constants
IMPLICIT NONE
INTEGER,INTENT(INOUT)::n1
REAL(KIND=DBL),DIMENSION(n1),INTENT(INOUT)::a,am,x
INTEGER,DIMENSION(n1),INTENT(INOUT)::igraph
INTEGER,DIMENSION(*),INTENT(INOUT)::ig_m
INTEGER,INTENT(INOUT)::nm
INTEGER::init=0,icount=0,ix=0,n2,i,j,k,l,nlim,istop
INTEGER,INTENT(OUT)::iend
REAL(KIND=DBL)::ex=0.25d0,alimit=0d0,d,a1,xmax,eps1,ratio
INTEGER::flag2=0
!      save d,ex,alimit,nlim
SAVE
IF(GLOBALINIT_optimize.EQ.0)THEN
	init=0
	icount=0
ENDIF
IF(init.EQ.0)THEN
   IF(n1.LE.100)nlim=10
   IF(n1.GT.100.AND.n1.LE.200)nlim=20
   IF(n1.GT.200.AND.n1.LE.500)nlim=40
   IF(n1.GT.500)nlim=60
   d=0
   xmax=0
   iend=0
   CALL ReadElem_real('alimit',alimit)
   WRITE(*,*)'The lower limit of the alpha(i) alimit=',alimit
!   READ(*,*)alimit
   CALL ReadElem_integer('nlimit',nlim)
   WRITE(*,*)'The lower limit of the number of channels nlim=',nlim
!   READ(*,*)nlim
!       if(alimit.gt.0)alimit=dnou(1)/n/10
   WRITE(*,*)'LIMIT=',alimit,nlim
   WRITE(*,*)'DUMPFAC=',ex
!   goto13
!   RETURN
ELSE
     
!      if(alimit.gt.0)alimit=min(dnou(1)/n/10,dnou(1)/nlim/10)
!      if(alimit.gt.0)alimit=dnou(1)/n/10
   WRITE(*,*)'LIMIT=',alimit

!      w_mean=0
!      do i=1,n
!       w_mean=w_mean+a(i)*x(i)
!      enddo
!      xmax=0
!      do i=1,n
!       if(xmax.lt.abs(w_mean-x(i)))xmax=abs(w_mean-x(i))
!      enddo
!      if(xmax/w_mean.gt.d)return
!      d=xmax/w_mean
!      do i=1,n-1
!       do j=i+1,n
!        if(xmax.lt.abs(x(i)-x(j)))xmax=abs(x(i)-x(j))
!       enddo
!      enddo

   a1=0
   DO i=1,n1
       a1=a1+a(i)*x(i)**ex
   ENDDO

   xmax=0
   DO i=1,n1
       IF(i.EQ.1)THEN
           xmax=ABS(a1-x(1)**ex)
       ELSE
           xmax=MAX(xmax,ABS(a1-x(i)**ex))
       ENDIF
   ENDDO
        
   IF(init.EQ.1)THEN
       d=xmax
       init=2
   ENDIF
       
   WRITE(*,*)'This opt. max=',xmax,n1

   IF(xmax.LE.d)THEN
!       d=xmax
!      endif
      ig_m(1:n1)=igraph(1:n1)
      nm=n1
      am(1:n1)=a(1:n1)
      d=xmax
!       iend=max(iend-1,0)
   ELSE
!       iend=iend+1
   ENDIF

   WRITE(*,*)'All opt.    d=',d
   icount=icount+1
!      write(1001,*)icount,xmax
ENDIF
! 13    continue
DO i=1,n1
   a(i)=a(i)*x(i)**ex
ENDDO
       
n2=n1
! 11    continue
CALL Helac_normalize(a(1:n2),n2)

IF(alimit.EQ.0.AND.init.GT.0)RETURN

IF(init.GT.0)THEN
    j=n2
    DO k=1,n2
       IF(a(k).LT.alimit)j=j-1
    ENDDO
    IF(j.LE.nlim)RETURN
ENDIF
       
k=0
j=0
istop=0
DO WHILE(istop.EQ.0)
    k=k+1
    j=j+1
!       if(j.eq.n1-1)istop=1
    IF(j.EQ.n2)istop=1

    flag2=0
	IF(init.EQ.0)THEN
! check for equal gen.
        eps1=dnou(1)/dnou(10)**11
	    DO l=1,n2
            ratio=ABS(a(k)-a(l))/(a(k)+a(l))
            IF(k.NE.l.AND.ratio.LT.eps1.AND.a(k).GT.0)THEN
!Cosats G.Papadopoulos
               IF(idebug.EQ.1)THEN
                  ix=ix+1
                  WRITE(*,*)igraph(k),a(k)
                  WRITE(*,*)igraph(l),a(l)
                  WRITE(*,*)is(igraph(k),1,1:8)
                  WRITE(*,*)is(igraph(k),2,1:8)
                  WRITE(*,*)is(igraph(l),1,1:8)
                  WRITE(*,*)is(igraph(l),2,1:8)
                  WRITE(*,*)ix,'----------------------------------'
               ENDIF
!Costas G.Papadopoulos
               igraph(k:n2-1)=igraph(k+1:n2)
               igraph(n2)=0
               a(k:n2-1)=a(k+1:n2)
               a(n2)=0
               k=k-1
			   flag2=1
			   EXIT
!	  goto 12
            ENDIF
	    ENDDO
	ENDIF

! check for negligible gen.
	IF(flag2.EQ.0)THEN
		IF(a(k).LT.alimit.AND.init.GT.0)THEN
			igraph(k:n2-1)=igraph(k+1:n2)
			igraph(n2)=0
			a(k:n2-1)=a(k+1:n2)
			a(n2)=0
			k=k-1
		ENDIF
	ENDIF

! 12     continue
ENDDO

!k=k+1
!if(a(k).lt.alimit) k=k-1

IF(k.LT.n2)THEN
    n2=k
    CALL Helac_normalize(a(1:n2),n2)
!       goto11
ENDIF
       
IF(n2.LE.n1.AND.init.EQ.0)THEN
	d=0
	init=1
    DO k=1,n2
         a(k)=dnou(1)/n2
    ENDDO
ENDIF
n1=n2

END SUBROUTINE Helac_optimize

SUBROUTINE Helac_igen(n1,b,ig)
!       include 'declare.h'
USE Helac_ranmar_mod
INTEGER,INTENT(IN)::n1
REAL(KIND(1d0)),DIMENSION(1:n1),INTENT(IN)::b
INTEGER,INTENT(OUT)::ig
REAL(KIND(1d0))::r
INTEGER::i
!	ig=n1
!	RETURN
r=Helac_rnmy(0)
!PRINT *,'r=',r
DO i=1,n1
   ig=i
   IF(r-b(i).LT.0)RETURN
ENDDO
END SUBROUTINE Helac_igen
      
SUBROUTINE Helac_normalize(a,n1)
!       include 'declare.h'
INTEGER,INTENT(IN)::n1
REAL(KIND(1d0)),DIMENSION(n1),INTENT(INOUT)::a
REAL(KIND(1d0))::atot
INTEGER::i       
atot=0
DO i=1,n1
   atot=atot+a(i)
ENDDO
DO i=1,n1
   a(i)=a(i)/atot
ENDDO  

END SUBROUTINE Helac_normalize

!
!FUNCTION pdgt(ifif)
!INTEGER,INTENT(IN)::ifif
!INTEGER::pdgt,init=0,ifac
!INTEGER,DIMENSION(-100:100)::ipdt
!SAVE ipdt,init
!IF(init.EQ.0)THEN
!                include 'pdt.h'  to the representation in PDG
!	ipdt( 1)= 12
!	ipdt(-1)=-12
!	ipdt( 5)= 14
!	ipdt(-5)=-14
!	ipdt( 9)= 16
!	ipdt(-9)=-16
!	ipdt( 2)= 11
!	ipdt(-2)=-11
!	ipdt( 6)= 13
!	ipdt(-6)=-13
!	ipdt( 10)= 15
!	ipdt(-10)=-15
!	ipdt( 3)= 2
!	ipdt(-3)=-2
!	ipdt( 4)= 1
!	ipdt(-4)=-1
!	ipdt( 7)= 4
!	ipdt(-7)=-4
!	ipdt( 8)= 3
!	ipdt(-8)=-3
!	ipdt( 11)= 6
!	ipdt(-11)=-6
!	ipdt( 12)= 5
!	ipdt(-12)=-5
!	ipdt(35)=21
!	ipdt(31)=22
!	ipdt(32)=23
!	ipdt(33)=24
!	ipdt(34)=-24
!	ipdt(41)=25
! end of 'pdt.h'
!	init=1
!ENDIF
!IF(ABS(ifif).LT.100)THEN
!	pdgt=ipdt(ifif)
!	RETURN
!ENDIF
!ifac=1
!IF(ifif.LT.0)ifac=-1
!SELECT CASE(ABS(ifif))
! J/psi
!CASE(443011)
!	pdgt=443
! etac
!CASE(441001)
!	pdgt=441
! hc
!CASE(441111)
!	pdgt=10443
! chic0
!CASE(443101)
!	pdgt=10441
! chic1
!CASE(443111)
!	pdgt=20443
! chic2
!CASE(443121)
!	pdgt=445
! upsilon
!CASE(553011)
!	pdgt=553
! etab
!CASE(551001)
!	pdgt=551
! hb
!CASE(551111)
!	pdgt=10553
! chib0
!CASE(553101)
!	pdgt=10551
! chib1
!CASE(553111)
!	pdgt=20553
! chib2
!CASE(553121)
!	pdgt=555
! Bc*+-
!CASE(453011)
!	pdgt=543*ifac
! Bc+-
!CASE(451001)
!	pdgt=541*ifac
! Bc2*+-
!CASE(453121)
!	pdgt=545*ifac
! Bc0*+-
!CASE(453101)
!	pdgt=10541*ifac
! Bc1(H)+-
!CASE(453111)
!	pdgt=20543*ifac
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!                                                          !
!     To matching the notations in PYTHIA 8                !
!                                                          !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ccbar(3P0(8))
!CASE(443108)
!	pdgt=9910441
! bbbar(3P0(8))
!CASE(553108)
!	pdgt=9910551
! bbbar(3S1(8))
!CASE(553018)
!	pdgt=9900553
! bbbar(1S0(8))
!CASE(551008)
!	pdgt=9900551
! ccbar(3S1(8))
!CASE(443018)
!	pdgt=9900443
! ccbar(1S0(8))
!CASE(441008)
!	pdgt=9900441
!CASE DEFAULT
!	pdgt=ifif
!END SELECT
!RETURN
!END FUNCTION pdgt

FUNCTION Phegas_hmom(i,j)
USE Helac_Global
INTEGER,INTENT(IN)::i,j
INTEGER::kkkij
REAL(KIND=DBL)::Phegas_hmom
IF(nhad.EQ.n)THEN
	Phegas_hmom=Phegas_pmom(MOD(i-1,nhad)+1,MOD(j-1,5)+1)
	RETURN
ENDIF
kkkij=hadron2parton(MOD(i-1,nhad)+1)
IF(Quarkonium3(kkkij).EQ.0)THEN
	Phegas_hmom=Phegas_pmom(kkkij,MOD(j-1,5)+1)
ELSE
	Phegas_hmom=Phegas_pmom(kkkij,MOD(j-1,5)+1)+Phegas_pmom(kkkij+1,MOD(j-1,5)+1)
ENDIF
RETURN
END FUNCTION Phegas_hmom

SUBROUTINE SOR(A,N1,K,IOPT)
!-----------------------------------------------------------------------
!     Sort A(N) into ascending order
!     IOPT = 1 : return sorted A and index array K
!     IOPT = 2 : return index array K only
!-----------------------------------------------------------------------
!     DOUBLE PRECISION A(N),B(5000)
!     INTEGER N,I,J,IOPT,K(N),K1(5000),IL(5000),IR(5000)
INTEGER,INTENT(IN)::N1
REAL(KIND(1d0)),DIMENSION(N1),INTENT(INOUT)::A
REAL(KIND(1d0)),DIMENSION(N1)::B
INTEGER,DIMENSION(N1),INTENT(OUT)::K
INTEGER,DIMENSION(N1)::K1
INTEGER,INTENT(IN)::IOPT
INTEGER::I,J,flag10,flag8,flag81,flag82
INTEGER,DIMENSION(N1)::IL,IR
!     IF (N.GT.5000) then
!       write(*,*) 'Too many entries to sort in srt, stop'
!       stop
!     endif
IF(N1.LE.0)RETURN
IL(1)=0
IR(1)=0
DO I=2,N1
   IL(I)=0
   IR(I)=0
   J=1
   flag10=0
   DO
    DO  ! find the smallest index j that A(j)<A(I), 
	    ! otherwise find the largest index j that A(j)>=A(I). 
      IF(A(I).GT.A(J)) EXIT   ! goto 5      ! 2
	  ! don't search the index smaller than j
      IF(IL(J).EQ.0) EXIT     ! goto 4      ! 3
      J=IL(J)
    ENDDO                      !GOTO 2
	flag10=0
	! the largest index j that A(j)>=A(I).
	IF(A(I).LE.A(J))THEN
      IR(I)=-J                !  4
      IL(J)=I
	  flag10=1
      EXIT                     ! goto 10
	ENDIF
	! don't search the index smaller than j ! 5
    IF(IR(J).LE.0) EXIT       ! goto 6
    J=IR(J)
   ENDDO                             !GOTO 2
   IF(flag10.EQ.1)THEN
     CYCLE
   ELSE
     IR(I)=IR(J)                 ! 6
     IR(J)=I
   ENDIF
ENDDO
!  10  CONTINUE
I=1
J=1
flag8=0
flag81=0
!      GOTO 8
DO
  DO
    IF(flag8.EQ.1.AND.flag81.EQ.0)THEN
       J=IL(J)                ! 20
    ENDIF
	! when the element A(k)<A(j) (at the right hand of j) and k>j
    IF(IL(J).GT.0)THEN        ! GOTO 20 ! 8
       flag8=1
	   flag81=0
	ELSE
       EXIT
    ENDIF
  ENDDO

  DO
     flag81=0
     flag82=0
     K(I)=J               !9
     B(I)=A(J)
     I=I+1
     IF(IR(J).GT.0)THEN               !12,30,13
        J=IR(J)                        ! 13
        flag81=1
        EXIT                           ! goto 8
	 ! this index is must the nearest to J
     ELSEIF(IR(J).LT.0)THEN
       J=-IR(J)                       ! 12
                                 !   GOTO 9
     ELSE
       flag82=1
       EXIT
     ENDIF
  ENDDO
                                 !  30  CONTINUE
  IF(flag82.EQ.1)EXIT
ENDDO
DO I=1,N1
   K1(I)=K(N1+1-I)
ENDDO
DO I=1,N1
   K(I)=K1(I)
ENDDO
IF(IOPT.EQ.2) RETURN
DO I=1,N1          ! 31
  A(I)=B(N1+1-I)
ENDDO
END SUBROUTINE SOR
END MODULE SinglePro
