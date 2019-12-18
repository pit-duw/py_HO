C******************************************************
C*          HELAC-Onia - Pythia interface.            *
C*    A Modification version from ME2pythia.f         *
C*====================================================*
C* WARNING: remove UPINIT and UPEVNT in pythia 6      *
C*====================================================*
C*            Version 3.9, 3 Dec 2007                 *
C*                                                    *
C*  Version 3.9                                       *
C*                                                    *
C*  - Updated to correspond to Pythia 6.4.14 and new  *
C*    MadEvent header format (MG/ME v. 4.1.47)        *
C*                                                    *
C*  Version 3.8                                       *
C*                                                    *
C*  - Give the event number in the event file in the  *
C*    new variable IEVNT in UPPRIV                    *
C*                                                    *
C*  Version 3.7                                       *
C*                                                    *
C*  - Set mass of massless outgoing particles to      *
C*    Pythia mass (PMAS(I,1))                         *
C*                                                    *
C*  Version 3.6                                       *
C*                                                    *
C*  - Removed the 1st # from the event file header    *
C*                                                    *
C*  Version 3.5                                       *
C*                                                    *
C*  - Reads according to the new LH event file format *
C*  - Now only LNHIN, LNHOUT and MSCAL in UPPRIV      *
C*                                                    *
C*  Version 3.4                                       *
C*                                                    *
C*  - Reads particle masses from event file           *
C*                                                    *
C*  Version 3.3                                       *
C*                                                    *
C*  - Added option MSCAL in common block UPPRIV to    *
C*    choose between fix (0) or event-based (1)       *
C*    scale for Pythia parton showering (SCALUP).     *
C*  - Fixed bug in reading the SLHA file              *
C*                                                    *
C*  Version 3.2                                       *
C*                                                    *
C*  - Reading the SLHA format param_card from the     *
C*    banner                                          *
C*  - Added support for lpp1/lpp2 = 2 or 3            *
C*  - Removed again support for different MadEvent    *
C*    processes in different files (no longer         *
C*    necessary with new multiple processes support   *
C*    in MadGraph/MadEvent                            *
C*                                                    *
C*  Version 3.1                                       *
C*  - Added support for different MadEvent processes  *
C*    in different files                              *
C*  - Fixed bug in e+e- collisions                    *
C*                                                    *
C*   Written by J.Alwall, alwall@slac.stanford.edu    *
C*     Earlier versions by S.Mrenna, M.Kirsanov       *
C*   Modified by H.-S.Shao, erdissshaw@gmail.com      *
C*                                                    *
C******************************************************
C*                                                    *
C* Instructions:                                      *
C* Please use the common block UPPRIV:                *
C* - The logical unit LNHIN must be an opened         *
C*   MadEvent event file                              *
C* - The output unit LNHOUT is by default 6 (std out) *
C* - Set MSCAL to 1 if a dynamical scale is desired   *
C*   for parton showers rather than the one given as  *
C*   factorization scale by HELAC-Onia (otherwise 0)  *
C* - IEVNT gives the number of the event in the event *
C*   file                                             *
C*                                                    *
C* Note that you need to link LHAPDF or PDFLIB        *
C*                                                    *
C******************************************************

C*********************************************************************
C...UPINIT
C...Routine called by PYINIT to set up user-defined processes.
C*********************************************************************
      SUBROUTINE UPINIT
      IMPLICIT NONE
      CALL UPINIT_EVT
      END

C*********************************************************************
C...UPINIT_LHE
C...Routine called by UPINIT for reading LHE
C*********************************************************************      
      SUBROUTINE UPINIT_LHE
      
      IMPLICIT NONE
      CHARACTER*132 CHAR_READ

C...Pythia parameters.
      INTEGER MSTP,MSTI,MRPY
      DOUBLE PRECISION PARP,PARI,RRPY
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &   LPRUP(MAXPUP)

C...Extra commonblock to transfer run info.
      INTEGER LNHIN,LNHOUT,MSCAL,IEVNT
      COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT
      DATA LNHIN,LNHOUT,MSCAL,IEVNT/77,6,1,0/
      SAVE/UPPRIV/

C...Parameter arrays (local)
      integer maxpara
      parameter (maxpara=100)
      integer npara,iseed
      character*20 param(maxpara),value(maxpara)      

C...Lines to read in assumed never longer than 200 characters. 
      INTEGER MAXLEN,IBEG,IPR,I
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING

C...Format for reading lines.
      CHARACTER*6 STRFMT
      STRFMT='(A000)'
      WRITE(STRFMT(3:5),'(I3)') MAXLEN

C...Extract the model parameter card and read it.
c     CALL MODELPAR(LNHIN)

c...Read the <init> block information

C...Loop until finds line beginning with "<init>" or "<init ". 
  100 READ(LNHIN,STRFMT,END=130,ERR=130) STRING ! LNHIN is the logic unit for LHE file
C...Pick out random number seed and use for PYR initialization
      IF(INDEX(STRING,'iseed').NE.0)THEN
         READ(STRING,*) iseed
         IF(iseed.gt.0) THEN
            WRITE(LNHOUT,*) 'Initializing PYR with random seed ',iseed
            MRPY(1) = iseed
            MRPY(2) = 0
         ENDIF
      ENDIF
      IBEG=0
  110 IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-5) GOTO 110 
      IF(STRING(IBEG:IBEG+5).NE.'<init>'.AND.
     &STRING(IBEG:IBEG+5).NE.'<init ') GOTO 100

C...Read first line of initialization info.
      READ(LNHIN,*,END=130,ERR=130) IDBMUP(1),IDBMUP(2),EBMUP(1),
     &EBMUP(2),PDFGUP(1),PDFGUP(2),PDFSUP(1),PDFSUP(2),IDWTUP,NPRUP

C...Read NPRUP subsequent lines with information on each process.
      DO 120 IPR=1,NPRUP
        READ(LNHIN,*,END=130,ERR=130) XSECUP(IPR),XERRUP(IPR),
     &  XMAXUP(IPR),LPRUP(IPR)
  120 CONTINUE

C...Set PDFLIB or LHAPDF pdf number for Pythia
cchen test
c      PDFSUP(1) = 1
c      PDFGUP(1) = 1
       print*,'MSTP(52)',MSTP(52)
       print*,'PDFSUP(1)',PDFSUP(1)

cchen      IF(PDFSUP(1).NE.19070.AND.(PDFSUP(1).NE.0.OR.PDFSUP(2).NE.0))THEN
c     Not CTEQ5L, which is standard in Pythia
cchen         CALL PYGIVE('MSTP(52)=2')
c     The following works for both PDFLIB and LHAPDF (where PDFGUP(1)=0)
c     But note that the MadEvent output uses the LHAPDF numbering scheme
cchen        IF(PDFSUP(1).NE.0)THEN
cchen           MSTP(51)=1000*PDFGUP(1)+PDFSUP(1)
cchen        ELSE
cchen           MSTP(51)=1000*PDFGUP(2)+PDFSUP(2)
cchen        ENDIF
cchen      ENDIF

C...Initialize widths and partial widths for resonances.
      CALL PYINRE
        
C...Calculate xsec reduction due to non-decayed resonances
C...based on first event only!

      CALL BRSUPP

      REWIND(LNHIN)

C...Reset event numbering
      IEVNT=0

      RETURN

C...Error exit: give up if initalization does not work.
  130 WRITE(*,*) ' Failed to read LHEF initialization information.'
      WRITE(*,*) ' Event generation will be stopped.'
      STOP  
      END

C*********************************************************************
C...UPINIT_EVT
C...Routine called by UPINIT
C*********************************************************************
      SUBROUTINE UPINIT_EVT
      USE HELAC_GLOBAL
      IMPLICIT NONE

C...Pythia parameters.                                                          
      INTEGER MSTP,MSTI,MRPY
      DOUBLE PRECISION PARP,PARI,RRPY
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)
      COMMON/PYDATR/MRPY(6),RRPY(100)

C...User process initialization commonblock.                                    
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &   LPRUP(MAXPUP)

C...Extra commonblock to transfer run info.                                     
      INTEGER LNHIN,LNHOUT,MSCAL,IEVNT
      COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT
      DATA LNHIN,LNHOUT,MSCAL,IEVNT/77,6,1,0/
      SAVE/UPPRIV/

      MRPY(1)=7
      MRPY(2)=0
      
      SELECT CASE(COLL_TYPE)
      CASE(1)
         IDBMUP(1)=2212
         IDBMUP(2)=2212
      CASE(2)
         IDBMUP(1)=2212
         IDBMUP(2)=-2212
      CASE DEFAULT
         IDBMUP(1)=11
         IDBMUP(2)=-11
      END SELECT
      
      IDWTUP=1

      EBMUP(1)=ebeam(1)
      EBMUP(2)=ebeam(2)
      
      PDFGUP(1)=iPDFGUP1
      PDFGUP(2)=iPDFGUP2
      PDFSUP(1)=iPDFSUP1
      PDFSUP(2)=iPDFSUP2

      NPRUP=1

C...Initialize widths and partial widths for resonances.
      CALL PYINRE

C...Calculate xsec reduction due to non-decayed resonances                      
C...based on first event only!                                                  

      CALL BRSUPP

C...Reset event numbering                                                       
      IEVNT=0

      RETURN

      END
C*********************************************************************      
C...UPEVNT
C...Routine called by PYEVNT or PYEVNW to get user process event
C*********************************************************************
      SUBROUTINE UPEVNT
      IMPLICIT NONE
      CALL UPEVNT_EVT
      END
C*********************************************************************
C...UPEVNT_LHE
C...Routine called by UPEVNT
C*********************************************************************
      SUBROUTINE UPEVNT_LHE

      IMPLICIT NONE

C...Pythia parameters.
      INTEGER MSTP,MSTI
      DOUBLE PRECISION PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &   LPRUP(MAXPUP)
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP)
C...Pythia common blocks
      INTEGER PYCOMP,KCHG,MINT,NPART,NPARTD,IPART,MAXNUR
      DOUBLE PRECISION PMAS,PARF,VCKM,VINT,PTPART
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYINT1/MINT(400),VINT(400)
      PARAMETER (MAXNUR=1000)
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)

C...Extra commonblock to transfer run info.
      INTEGER LNHIN,LNHOUT,MSCAL,IEVNT
      COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT

C...Local variables
      INTEGER I,J,IBEG,NEX,KP(MAXNUP),MOTH,NUPREAD,IREM
      DOUBLE PRECISION PSUM,ESUM,PM1,PM2,A1,A2,A3,A4,A5
      DOUBLE PRECISION SCALLOW(MAXNUP),PNONJ(4),PMNONJ!,PT2JETS
C...Lines to read in assumed never longer than 200 characters. 
      INTEGER MAXLEN
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING

C...Format for reading lines.
      CHARACTER*6 STRFMT
      CHARACTER*1 CDUM

      STRFMT='(A000)'
      WRITE(STRFMT(3:5),'(I3)') MAXLEN

C...Loop until finds line beginning with "<event>" or "<event ". 
  100 READ(LNHIN,STRFMT,END=130,ERR=130) STRING
      IBEG=0
  110 IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-6) GOTO 110 
      IF(STRING(IBEG:IBEG+6).NE.'<event>'.AND.
     &STRING(IBEG:IBEG+6).NE.'<event ') GOTO 100

C...Read first line of event info.
      READ(LNHIN,*,END=130,ERR=130) NUPREAD,IDPRUP,XWGTUP,SCALUP,
     &AQEDUP,AQCDUP

C...Read NUP subsequent lines with information on each particle.
      ESUM=0d0
      PSUM=0d0
      NEX=2 ! number of external particles (assume two intial particles)
      NUP=1
      IREM=NUPREAD+1

      DO 120 I=1,NUPREAD
        READ(LNHIN,*,END=130,ERR=130) IDUP(NUP),ISTUP(NUP),
     &  MOTHUP(1,NUP),MOTHUP(2,NUP),ICOLUP(1,NUP),ICOLUP(2,NUP),
     &  (PUP(J,NUP),J=1,5),VTIMUP(NUP),SPINUP(NUP)
cchen! read out Ah momentum
          
C...Reset resonance momentum to prepare for mass shifts
        IF(ISTUP(NUP).EQ.2) PUP(3,NUP)=0 ! ISTUP: 2 resonance, 1 final, -1 initial
        IF(ISTUP(NUP).EQ.1)THEN
          NEX=NEX+1
           IF(PUP(5,NUP).EQ.0D0.AND.IABS(IDUP(NUP)).GT.3
     $         .AND.IDUP(NUP).NE.21) THEN ! set masses (c,b,t masses)
C...Set massless particle masses to Pythia default. Adjust z-momentum. 
              PUP(5,NUP)=PMAS(IABS(PYCOMP(IDUP(NUP))),1)
              PUP(3,NUP)=SIGN(SQRT(MAX(0d0,PUP(4,NUP)**2-PUP(5,NUP)**2-
     $           PUP(1,NUP)**2-PUP(2,NUP)**2)),PUP(3,NUP))
           ENDIF
           PSUM=PSUM+PUP(3,NUP)
C...Adjust mother information due to removed mother
           IF(MOTHUP(1,NUP).EQ.IREM)THEN
              MOTHUP(1,NUP)=1
              MOTHUP(2,NUP)=2
           ELSE IF(MOTHUP(1,NUP).GT.IREM)THEN
              MOTHUP(1,NUP)=MOTHUP(1,NUP)-1
              MOTHUP(2,NUP)=MOTHUP(2,NUP)-1
           ENDIF
C...Set mother resonance momenta
           MOTH=MOTHUP(1,NUP)
           DO WHILE (MOTH.GT.2)
             PUP(3,MOTH)=PUP(3,MOTH)+PUP(3,NUP)
             MOTH=MOTHUP(1,MOTH)
           ENDDO
        ENDIF
        NUP=NUP+1
  120 CONTINUE
      NUP=NUP-1

C...Increment event number
      IEVNT=IEVNT+1

C..Adjust mass of resonances
      DO I=1,NUP
         IF(ISTUP(I).EQ.2)THEN
            PUP(5,I)=SQRT(PUP(4,I)**2-PUP(1,I)**2-PUP(2,I)**2-
     $             PUP(3,I)**2)
         ENDIF
      ENDDO

      ESUM=PUP(4,1)+PUP(4,2)

C...Assuming massless incoming particles - otherwise Pythia adjusts
C...the momenta to make them massless
      IF(IDBMUP(1).GT.100.AND.IDBMUP(2).GT.100)THEN
        DO I=1,2
          PUP(3,I)=0.5d0*(PSUM+SIGN(ESUM,PUP(3,I)))
          PUP(5,I)=0d0
        ENDDO
        PUP(4,1)=ABS(PUP(3,1))
        PUP(4,2)=ESUM-PUP(4,1)
      ENDIF
        
C...If you want to use some other scale for parton showering than the 
C...factorisation scale given by HELAC-Onia, please implement the function PYMASC
C...(example function included below) 

      IF(MSCAL.GT.0) CALL PYMASC(SCALUP)
      
      RETURN

C...Error exit, typically when no more events.
  130 WRITE(*,*) ' Failed to read LHEF event information,'
      WRITE(*,*) ' assume end of file has been reached.'
      NUP=0
      MINT(51)=2
      RETURN
      END

C*********************************************************************
C...UPEVNT_EVT
C...Routine called by UPEVNT
C*********************************************************************
      SUBROUTINE UPEVNT_EVT
      USE HELAC_Global
      IMPLICIT NONE
C...Pythia parameters.                                                          
      INTEGER MSTP,MSTI
      DOUBLE PRECISION PARP,PARI
      COMMON/PYPARS/MSTP(200),PARP(200),MSTI(200),PARI(200)

C...User process initialization commonblock.                                    
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &   LPRUP(MAXPUP)

C...User process event common block.                                            
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP)

C...Pythia common blocks
      INTEGER PYCOMP,KCHG,MINT,NPART,NPARTD,IPART,MAXNUR
      DOUBLE PRECISION PMAS,PARF,VCKM,VINT,PTPART
C...Particle properties + some flavour parameters.
      COMMON/PYDAT2/KCHG(500,4),PMAS(500,4),PARF(2000),VCKM(4,4)
      COMMON/PYINT1/MINT(400),VINT(400)
      PARAMETER (MAXNUR=1000)
      COMMON/PYPART/NPART,NPARTD,IPART(MAXNUR),PTPART(MAXNUR)

C...Local variables
      INTEGER I,J,IBEG,NEX,KP(MAXNUP),MOTH,NUPREAD,IREM
      DOUBLE PRECISION PSUM,ESUM,PM1,PM2,A1,A2,A3,A4,A5
      DOUBLE PRECISION SCALLOW(MAXNUP),PNONJ(4),PMNONJ

      NUPREAD = nhad
      IDPRUP = 81
      XWGTUP = 1d0

      END

C*********************************************************************
C   PYMASC
C   Implementation of scale used in Pythia parton showers
C*********************************************************************
      SUBROUTINE PYMASC(scale)
      IMPLICIT NONE

C...Arguments
      REAL*8 scale

C...Functions
      REAL*8 SMDOT5

C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &   LPRUP(MAXPUP)
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP)

C...Extra commonblock to transfer run info.
      INTEGER LNHIN,LNHOUT,MSCAL,IEVNT
      COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT

C...Local variables
      INTEGER ICC1,ICC2,IJ,IDC1,IDC2,IC,IC1,IC2
      REAL*8 QMIN,QTMP

C   Just use the scale read off the event record
      scale=SCALUP

C   Alternatively:

C...  Guesses for the correct scale
C     Assumptions:
C     (1) if the initial state is a color singlet, then
C     use s-hat for the scale
C     
C     (2) if color flow to the final state, use the minimum
C     of the dot products of color connected pairs
C     (times two for consistency with above)

        QMIN=SMDOT5(PUP(1,1),PUP(1,2))
        ICC1=1
        ICC2=2
C     
C     For now, there is no generic way to guarantee the "right"
C     scale choice.  Here, we take the HERWIG pt. of view and
C     choose the dot product of the colored connected "primary"
C     pairs.
C     

        DO 101 IJ=1,NUP
          IF(MOTHUP(2,IJ).GT.2) GOTO 101
          IDC1=ICOLUP(1,IJ)
          IDC2=ICOLUP(2,IJ)
          IF(IDC1.EQ.0) IDC1=-1
          IF(IDC2.EQ.0) IDC2=-2
          
          DO 201 IC=IJ+1,NUP
            IF(MOTHUP(2,IC).GT.2) GOTO 201
            IC1=ICOLUP(1,IC)
            IC2=ICOLUP(2,IC)
            IF(ISTUP(IC)*ISTUP(IJ).GE.1) THEN
              IF(IDC1.EQ.IC2.OR.IDC2.EQ.IC1) THEN
                QTMP=SMDOT5(PUP(1,IJ),PUP(1,IC))
                IF(QTMP.LT.QMIN) THEN
                  QMIN=QTMP
                  ICC1=IJ
                  ICC2=IC
                ENDIF
              ENDIF
            ELSEIF(ISTUP(IC)*ISTUP(IJ).LE.-1) THEN
              IF(IDC1.EQ.IC1.OR.IDC2.EQ.IC2) THEN
                QTMP=SMDOT5(PUP(1,IJ),PUP(1,IC))          
                IF(QTMP.LT.QMIN) THEN
                  QMIN=QTMP
                  ICC1=IJ
                  ICC2=IC
                ENDIF
              ENDIF
            ENDIF
 201      CONTINUE
 101    CONTINUE

        scale=QMIN

      RETURN
      END

C...SMDOT5
C   Helper function

      FUNCTION SMDOT5(V1,V2)
      IMPLICIT NONE
      REAL*8 SMDOT5,TEMP
      REAL*8 V1(5),V2(5)
      INTEGER I

      SMDOT5=0D0
      TEMP=V1(4)*V2(4)
      DO I=1,3
        TEMP=TEMP-V1(I)*V2(I)
      ENDDO

      SMDOT5=SQRT(ABS(TEMP))

      RETURN
      END

C*********************************************************************

      subroutine BRSUPP
c take into account the effect of branching ratios for unstable particles
      IMPLICIT NONE

C...Three Pythia functions return integers, so need declaring.
      INTEGER PYCOMP,MWID
      DOUBLE PRECISION WIDS

C...Resonance width and secondary decay treatment.
      COMMON/PYINT4/MWID(500),WIDS(500,5)

C...User process initialization commonblock.
      INTEGER MAXPUP
      PARAMETER (MAXPUP=100)
      INTEGER IDBMUP,PDFGUP,PDFSUP,IDWTUP,NPRUP,LPRUP
      DOUBLE PRECISION EBMUP,XSECUP,XERRUP,XMAXUP
      COMMON/HEPRUP/IDBMUP(2),EBMUP(2),PDFGUP(2),PDFSUP(2),
     &   IDWTUP,NPRUP,XSECUP(MAXPUP),XERRUP(MAXPUP),XMAXUP(MAXPUP),
     &   LPRUP(MAXPUP)
C...User process event common block.
      INTEGER MAXNUP
      PARAMETER (MAXNUP=500)
      INTEGER NUP,IDPRUP,IDUP,ISTUP,MOTHUP,ICOLUP
      DOUBLE PRECISION XWGTUP,SCALUP,AQEDUP,AQCDUP,PUP,VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,IDUP(MAXNUP),
     &   ISTUP(MAXNUP),MOTHUP(2,MAXNUP),ICOLUP(2,MAXNUP),PUP(5,MAXNUP),
     &   VTIMUP(MAXNUP),SPINUP(MAXNUP)

C...Extra commonblock to transfer run info.
      INTEGER LNHIN,LNHOUT,MSCAL,IEVNT
      COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT

      INTEGER I,J,IBEG
      REAL SUPPCS

C...Lines to read in assumed never longer than 200 characters. 
      INTEGER MAXLEN
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING

C...Format for reading lines.
      CHARACTER*6 STRFMT
      STRFMT='(A000)'
      WRITE(STRFMT(3:5),'(I3)') MAXLEN

C...Loop until finds line beginning with "<event>" or "<event ". 
  100 READ(LNHIN,STRFMT,END=130,ERR=130) STRING
      IBEG=0
  110 IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-6) GOTO 110 
      IF(STRING(IBEG:IBEG+6).NE.'<event>'.AND.
     &STRING(IBEG:IBEG+6).NE.'<event ') GOTO 100

C...Read first line of event info.
      READ(LNHIN,*,END=130,ERR=130) NUP,IDPRUP,XWGTUP,SCALUP,
     &AQEDUP,AQCDUP

C...Read NUP subsequent lines with information on each particle.
      DO 120 I=1,NUP
        READ(LNHIN,*,END=130,ERR=130) IDUP(I),ISTUP(I),
     &  MOTHUP(1,I),MOTHUP(2,I),ICOLUP(1,I),ICOLUP(2,I),
     &  (PUP(J,I),J=1,5),VTIMUP(I),SPINUP(I)
  120 CONTINUE


      SUPPCS=1.
      DO I=3,NUP
        IF (ISTUP(I).EQ.1.AND.(IABS(IDUP(I)).GE.23.OR.
     $     (IABS(IDUP(I)).GE.6.AND.IABS(IDUP(I)).LE.8)))
     $     THEN
          WRITE(LNHOUT,*) 'Resonance ',IDUP(I), ' has BRTOT ',
     $       wids(PYCOMP(IDUP(I)),2)
          if(wids(PYCOMP(IDUP(I)),2).lt.0.95) then ! branching ratio of swiched on( use MDME to see) decay channels over all channels
            write(*,*) 'Decreasing cross section!'
            SUPPCS=SUPPCS*wids(PYCOMP(IDUP(I)),2)
          endif
        endif
      enddo
      if(SUPPCS.gt.0)then
         write(*,*)'Multiplying cross section by ',SUPPCS
         XSECUP(1)=XSECUP(1)*SUPPCS
      else
         write(*,*) 'Warning! Got cross section suppression 0'
         write(*,*) 'No cross section reduction done'
      endif
      RETURN

C...Error exit, typically when no more events.
  130 WRITE(*,*) ' Failed to read LHEF event information.'
      STOP
      END
      
      subroutine case_trap2(name,n)
c**********************************************************
c   change the string to lowercase if the input is not
c**********************************************************
      implicit none
c   
c   ARGUMENT
c   
      character(*) name
      integer n
c   
c   LOCAL
c   
      integer i,k

      do i=1,n
        k=ichar(name(i:i))
        if(k.ge.65.and.k.le.90) then !upper case A-Z
          k=ichar(name(i:i))+32   
          name(i:i)=char(k)        
        endif
      enddo

      return
      end

