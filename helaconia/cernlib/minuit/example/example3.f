      PROGRAM example3
      IMPLICIT NONE
      EXTERNAL FCNMA
      REAL(KIND(1d0)),EXTERNAL::F3
      INTEGER,DIMENSION(3)::NPRM
      REAL(KIND(1d0)),DIMENSION(3)::VSTRT,STP, ARGLIS
      CHARACTER*10,DIMENSION(3)::PNAM
      DATA NPRM / 1 , 2, 3 /
      DATA PNAM /'a','b','c'/
      DATA VSTRT/ 0.,0.,0./
      DATA STP/0.01,0.01,0.01/
      REAL(KIND(1d0)),PARAMETER::ZERO = 0d0
      REAL(KIND(1d0)),DIMENSION(3)::VANS,VERR
      CHARACTER*10,DIMENSION(3)::PCHAR
      INTEGER::i,IERFLG
C      Initialize Minuit, define I/O unit numbers
      OPEN(UNIT=16,FILE="example3.out")
      CALL MNINIT(5,16,7)
C      Define parameters, set initial values
      DO i=1,3
         CALL MNPARM(NPRM(i),PNAM(i),VSTRT(i),STP(i),ZERO,ZERO,IERFLG)
       IF (IERFLG .NE. 0)  THEN
          WRITE (16,'(A,I2)')  ' UNABLE TO DEFINE PARAMETER NO.',I
          STOP
       ENDIF
      ENDDO
      CALL MNSETI('test program for fit')
      ARGLIS(1) = 1.
      CALL MNEXCM(FCNMA, 'CALL FCN', ARGLIS, 1, IERFLG,F3)
      CALL MNEXCM(FCNMA, 'MIGRAD', ARGLIS, 0 , IERFLG,F3)
      CALL MNEXCM(FCNMA, 'MINOS', ARGLIS, 0, IERFLG, F3)
      ARGLIS(1) = 3.
      CALL MNEXCM(FCNMA, 'CALL FCN', ARGLIS, 1, IERFLG,F3)
      CALL MNEXCM(FCNMA,'STOP ', 0,0,IERFLG, F3)
      CLOSE(UNIT=16)
      END

      SUBROUTINE FCNMA(NPAR,GIN,F,X,IFLAG,FUTIL)
      IMPLICIT NONE
      REAL(KIND(1d0)),DIMENSION(*)::X,GIN
      INTEGER,PARAMETER::MXBIN=60
      INTEGER::NPAR,IFLAG
      REAL(KIND(1d0))::F,THVAL,XX
      REAL(KIND(1d0)),EXTERNAL::FUTIL
      REAL(KIND(1d0)),DIMENSION(MXBIN)::XVAL,EXVAL,EXERR
      SAVE XVAL,EXVAL,EXERR
      REAL(KIND(1d0))::CHISQ,CHI
      INTEGER::i
      IF(IFLAG.EQ.1)THEN
         OPEN(UNIT=2231,FILE="example3.dat")
         DO i=1,MXBIN
            READ(2231,*)XVAL(i),EXVAL(i),EXERR(i)
         ENDDO
         CLOSE(UNIT=2231)
      ENDIF
      CHISQ=0d0
      DO i=1,MXBIN
         XX=XVAL(i)
         THVAL=FUTIL(XX,X(1),X(2),X(3))
         CHI=(THVAL-EXVAL(i))**2/EXERR(i)**2
         CHISQ=CHISQ+CHI
      ENDDO
      F=CHISQ
      IF(IFLAG.EQ.3)THEN
         PRINT *,X(1),X(2),X(3)
         PRINT *,F
      ENDIF
      RETURN
      END

      FUNCTION F3(x,a,b,c)
      IMPLICIT NONE
      REAL(KIND(1d0))::x,a,b,c,F3
      F3=a*x**3+b*x+c
      END
