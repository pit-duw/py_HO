MODULE PY6_interface
!**************************************************************
!*              HELAC-Onia - PYTHIA6 interface                *
!*------------------------------------------------------------*
!*  Written by H.-S.Shao, erdissshaw@gmail.com, 08 Jan 2014   *
!*============================================================*
!*                                                            *
!* INSTRUCTION: it will use PYTHIA6/HO2PY6.f                  *
!*                                                            *
!**************************************************************
  USE Helac_Global
  IMPLICIT NONE
CONTAINS
  SUBROUTINE PY6RUN
    IMPLICIT NONE
    ! Pythia parameters.
    INTEGER,DIMENSION(200)::MSTP,MSTI
    REAL(KIND(1d0)),DIMENSION(200)::PARP,PARI
    COMMON/PYPARS/MSTP,PARP,MSTI,PARI
    ! Make sure PYDATA is linked
    EXTERNAL PYDATA
    ! Event record
    INTEGER::N_PY6,NPAD_PY6
    INTEGER,DIMENSION(4000,5)::K_PY6
    REAL(KIND(1d0)),DIMENSION(4000,5)::P_PY6,V_PY6
    COMMON/PYJETS/N_PY6,NPAD_PY6,K_PY6,P_PY6,V_PY6
    ! Add IEVNT to be consistent with HO2PY6.f
    ! Extra commonblock to transfer run info
    INTEGER::LNHIN,LNHOUT,MSCAL,IEVNT
    COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT
    ! Local variables
    INTEGER::NEV,IEV
    CHARACTER(len=5)::CGIVE
    CHARACTER(len=30)::CGIVE0
    ! User process initialization commonblock
    INTEGER,PARAMETER::MAXPUP=100
    INTEGER,DIMENSION(2)::IDBMUP,PDFGUP,PDFSUP
    INTEGER::IDWTUP,NPRUP
    INTEGER,DIMENSION(MAXPUP)::LPRUP
    REAL(KIND(1d0)),DIMENSION(2)::EBMUP
    REAL(KIND(1d0)),DIMENSION(MAXPUP)::XSECUP,XERRUP,XMAXUP
    COMMON/HEPRUP/IDBMUP,EBMUP,PDFGUP,PDFSUP,&
         IDWTUP,NPRUP,XSECUP,XERRUP,XMAXUP,LPRUP
    CHARACTER(len=300)::lheinfile,lheoutfile
    lheinfile=TRIM(output_dir)//'sample'//TRIM(process(2:10*nhad+1))//'.lhe'
    !lheoutfile=TRIM(output_dir)//'sample'//TRIM(process(2:10*nhad+1))//'_PY6.lhe'
    ! Maximum number of events to generate
    ! -1 means all availabel events
    NEV=-1
    ! Initialize HEP logical units
    OPEN(LNHIN,FILE=TRIM(lheinfile))
    ! Set pythia output to LNHOUT
    WRITE(CGIVE,'(I5)')LNHOUT
    
    CGIVE0='MSTU(11)='//CGIVE
    CALL PYGIVE(CGIVE0)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Switches in PY6
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Set pi0 stable to trim event listings
    CALL PYGIVE('MDCY(C111,1)=0')
    CALL PYGIVE('MDCY(C211,1)=0')
    
    ! ISR MSTP(61)
    CALL PYGIVE('MSTP(61)=1')
    ! FSR MSTP(71)
    CALL PYGIVE('MSTP(71)=1')
    ! Multiple interaction
    CALL PYGIVE('MSTP(81)=0')
    ! Hadronization
    CALL PYGIVE('MSTP(111)=0')

    ! Initialize with external process
    CALL PYINIT('USER',' ',' ',0D0)

    ! Read out total number of event generated
    CALL SYSTEM(" grep '<event>' "//lheinfile//" | wc -l > nevents.dat")
    OPEN(UNIT=22210100,FILE='nevents.dat',status='OLD')
    READ(22210100,*)TOTEVNT
    CLOSE(22210100)

    IEV=0
    
    ! Event loop
    DO WHILE(IEV.LT.NEV.OR.NEV.LT.0)
       IEV=IEV+1
       IF(MOD(IEV,TOTEVNT/10).EQ.0)THEN
          WRITE(*,*)(1.0*IEV/TOTEVNT)*100.0,"% done"
       ENDIF
       ! Get event, it will call many subroutines
       CALL PYEVNT
       ! If event generation failed,quit
       IF(MSTI(51).EQ.1)THEN
          CALL PYSTAT(1)
          WRITE(*,*)'IEVNT=',IEVNT
          CLOSE(LNHIN)
          CLOSE(20)
          STOP
       ENDIF
       ! write down the first event two information
       IF(IEV.LT.2)THEN
          CALL PYLIST(2)
       ENDIF
       CALL WRITE_EVENT_LHE
    ENDDO
  END SUBROUTINE PY6RUN

  SUBROUTINE WRITE_EVENT_LHE
    IMPLICIT NONE
    ! The event record.
    INTEGER::N_PY6,NPAD_PY6
    INTEGER,DIMENSION(4000,5)::K_PY6
    REAL(KIND(1d0)),DIMENSION(4000,5)::P_PY6,V_PY6
    COMMON/PYJETS/N_PY6,NPAD_PY6,K_PY6,P_PY6,V_PY6               
    ! User process event common block.
    INTEGER,PARAMETER::MAXNUP=500
    INTEGER::NUP,IDPRUP,
    INTEGER,DIMENSION(MAXNUP)::IDUP,ISTUP
    INTEGER,DIMENSION(2,MAXNUP)::MOTHUP,ICOLUP
    REAL(KIND(1d0))::XWGTUP,SCALUP,AQEDUP,AQCDUP
    REAL(KIND(1d0)),DIMENSION(5,MAXNUP)::PUP
    REAL(KIND(1d0)),DIMENSION(MAXNUP)::VTIMUP,SPINUP
    COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP,&
         ISTUP,MOTHUP,ICOLUP,PUP,VTIMUP,SPINUP
    INTEGER::I,J
    REAL(KIND(1d0))::Pxmiss,Pymiss,Etmiss
    CHARACTER(len=300)::lheoutfile
    INTEGER::init=0
    SAVE lheoutfile,init
    IF(init.EQ.0)THEN
       lheoutfile=TRIM(output_dir)//'sample'//TRIM(process(2:10*nhad+1))//'_PY6.lhe'
       OPEN(UNIT=20, FILE=TRIM(lheoutfile), status="unknown")
       init=1
    ENDIF
    ! Perform event analysis to generate .lhe and .lhco files
    Pxmiss=0
    Pymiss=0
    NUP=0
    ! first four events are beam hadron and initial state without ISR 
    ! thus, we begin at I = 5
    DO I=5,N_PY6
       IF(K_PY6(I,1).GT.20)THEN
          NUP=NUP+1
          IDUP(NUP)=K_PY6(I,2)
          IF(I.LE.6)THEN
             ISTUP(NUP)=-1
             MOTHUP(1,NUP)=0
          ELSE
             ISTUP(NUP)=2
             MOTHUP(1,NUP)=MAX(0,K_PY6(I,3)-4)
          ENDIF
          MOTHUP(2,NUP)=MOTHUP(1,NUP)
          DO J=1,5
             PUP(J,NUP)=P(I,J)
          ENDDO
       ENDIF
       ! lepton and photon final state
       IF(K_PY6(I,1).LT.10.AND.(IABS(K_PY6(I,2)).EQ.11.OR.&
            IABS(K_PY6(I,2)).EQ.13.OR.IABS(K_PY6(I,2)).EQ.15.OR.&
            IABS(K_PY6(I,2)).EQ.22)&
            .AND.K_PY6(K_PY6(I,3),2).EQ.K_PY6(I,2))THEN
          NUP=NUP+1
          IDUP(NUP)=K_PY6(I,2)
          ISTUP(NUP)=1
          MOTHUP(1,NUP)=0
          MOTHUP(2,NUP)=MOTHUP(1,NUP)
          DO J=1,5
             PUP(J,NUP)=P(I,J)
          ENDDO
       ELSEIF(K_PY6(I,1).LT.10.AND.(IABS(K_PY6(I,2)).EQ.12.OR.&
            IABS(K_PY6(I,2)).EQ.14.OR.IABS(K_PY6(I,2)).EQ.16.OR.&
            K_PY6(I,2).EQ.1000022)THEN
          ! neutrino final state
          Pxmiss=Pxmiss+P(I,1)
          Pymiss=Pymiss+P(I,2)
       ELSEIF(K_PY6(I,1).LT.10.AND.&
            (ABS(K_PY6(I,2)).LT.6.OR.K_PY6(I,2).EQ.21))THEN
          NUP=NUP+1
          IDUP(NUP)=K_PY6(I,2)
          ISTUP(NUP)=1
          MOTHUP(1,NUP)=0
          MOTHUP(2,NUP)=MOTHUP(1,NUP)
          DO J=1,5
             PUP(J,NUP)=P(I,J)
          ENDDO
       ENDIF
    ENDDO
    ! Missing transverse energy
    Etmiss=DSQRT(Pxmiss**2+Pymiss**2)
    IF(Etmiss.GT.0d0)THEN
       NUP=NUP+1
       IDUP(NUP)=12
       ISTUP(NUP)=1
       MOTHUP(1,NUP)=0
       MOTHUP(2,NUP)=0
       PUP(1,NUP)=Pxmiss
       PUP(2,NUP)=Pymiss
       PUP(4,NUP)=Etmiss
       PUP(3,NUP)=0
       PUP(5,NUP)=0
    ENDIF

    ! Write event into lhe file
    WRITE(20,'(a)') '<event>'
    WRITE(20,'(I3,I4,4E16.8)')NUP,100,1d0,0d0,0d0,0d0
    DO I=1,NUP
       WRITE(20,'(I8,I3,2I3,2I2,5E16.8,2F3.0)')&
            IDUP(I),ISTUP(I),MOTHUP(1,I),MOTHUP(2,I),0,0,&
            (PUP(J,I),J=1,5),0d0,0d0
    ENDDO
    WRITE(20,'(a)') '</event>'
  END SUBROUTINE WRITE_EVENT_LHE
END MODULE PY6_interface
