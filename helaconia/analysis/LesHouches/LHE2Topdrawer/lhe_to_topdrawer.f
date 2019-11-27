      PROGRAM plot_lhe
      IMPLICIT NONE
      INTEGER LNHIN,LNHOUT,MSCAL,IEVNT
      COMMON/UPPRIV/LNHIN,LNHOUT,MSCAL,IEVNT
      DATA LNHIN,LNHOUT,MSCAL,IEVNT/77,6,1,0/
      SAVE/UPPRIV/
      CHARACTER*1000::LHEFILE
      LOGICAL::lexists
      INTEGER::reason,inindex
      LOGICAL::found
      ! get the files
      call system('ls ./ > fileContents.txt')
      open(2133,FILE='fileContents.txt',action="read")
      ! how many
      found=.FALSE.
      DO
         read(2133,FMT='(a)',iostat=reason)LHEFILE
         IF(reason.NE.0)EXIT
         inindex=0
         inindex=INDEX(LHEFILE,'.lhe')
         IF(inindex.GT.0)THEN
            found=.TRUE.
            EXIT
         ENDIF
      ENDDO      
!      OPEN(UNIT=2133,FILE="lhepath.inp")
!      READ(2133,*)LHEFILE
      CLOSE(UNIT=2133)
      INQUIRE(FILE=TRIM(LHEFILE),EXIST=lexists)
      IF((.NOT.found).OR.(.NOT.lexists))THEN
         WRITE(*,*)"ERROR:Cannot find the lhe file in the current directory"
         STOP
      ENDIF
      call system ("rm fileContents.txt")
      OPEN(UNIT=LNHIN,FILE=TRIM(LHEFILE))
      CALL READ_LHE_EVENT
      CLOSE(UNIT=LNHIN)
      END
      
      SUBROUTINE READ_LHE_EVENT
      
      IMPLICIT NONE
      CHARACTER*132 CHAR_READ

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
C      DATA LNHIN,LNHOUT,MSCAL,IEVNT/77,6,1,0/
      SAVE/UPPRIV/

C...Lines to read in assumed never longer than 200 characters. 
      INTEGER MAXLEN,IBEG,IPR,I,J,NUPREAD
      PARAMETER (MAXLEN=200)
      CHARACTER*(MAXLEN) STRING
      CHARACTER*6,DIMENSION(1)::weights_info
      INTEGER::nwgt
      REAL(KIND(1d0)),DIMENSION(1)::wgts
      REAL(KIND(1d0)),DIMENSION(0:4,MAXNUP)::pmom

C...Format for reading lines.
      CHARACTER*6 STRFMT
      STRFMT='(A000)'
      WRITE(STRFMT(3:5),'(I3)') MAXLEN

      nwgt=1
      weights_info(1)='lhe value'
      CALL analysis_begin(nwgt,weights_info)

C...Loop until finds line beginning with "<init>" or "<init ". 
  100 READ(LNHIN,STRFMT,END=130,ERR=130) STRING ! LNHIN is the logic unit for LHE file
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

      REWIND(LNHIN)

C...Reset event numbering
      IEVNT=0

C...Loop until finds line beginning with "<event>" or "<event ".
 200   READ(LNHIN,STRFMT,END=230,ERR=230) STRING
       IBEG=0
 210   IBEG=IBEG+1
C...Allow indentation.
      IF(STRING(IBEG:IBEG).EQ.' '.AND.IBEG.LT.MAXLEN-6) GOTO 210
      IF(STRING(IBEG:IBEG+6).NE.'<event>'.AND.
     &     STRING(IBEG:IBEG+6).NE.'<event ') GOTO 200

C...Read first line of event info.
      READ(LNHIN,*,END=230,ERR=230) NUPREAD,IDPRUP,XWGTUP,SCALUP,
     &     AQEDUP,AQCDUP

      wgts(1)=XWGTUP ! in unit of pb
      ! convert pb to nb
      wgts(1)=wgts(1)*1d-3

      NUP=1

      DO 220 I=1,NUPREAD
         READ(LNHIN,*,END=230,ERR=230) IDUP(NUP),ISTUP(NUP),
     &        MOTHUP(1,NUP),MOTHUP(2,NUP),ICOLUP(1,NUP),ICOLUP(2,NUP),
     &        (PUP(J,NUP),J=1,5),VTIMUP(NUP),SPINUP(NUP)
         pmom(0,NUP)=PUP(4,NUP)
         pmom(1:3,NUP)=PUP(1:3,NUP)
         pmom(4,NUP)=PUP(5,NUP)
         NUP=NUP+1
 220  CONTINUE
      NUP=NUP-1

      
      CALL analysis_fill(pmom(0:4,1:NUP),NUP,IDUP(1:NUP),wgts)

C...Increment event number
      IEVNT=IEVNT+1

      GO TO 200

      CALL analysis_end(1d0)

      RETURN

C...Error exit: give up if initalization does not work.
  130 WRITE(*,*) ' Failed to read LHEF initialization information.'
      WRITE(*,*) ' Event generation will be stopped.'
      CALL analysis_end(1d0)
      STOP
 230  WRITE(*,*) ' Failed to read LHEF event information,'
      WRITE(*,*) ' assume end of file has been reached.'
      WRITE(*,*) ' number of accepted events = ',IEVNT
      NUP=0
      CALL analysis_end(1d0)
      STOP
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

