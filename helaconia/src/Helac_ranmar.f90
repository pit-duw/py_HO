MODULE Helac_ranmar_mod
USE Helac_Global  ! we only use the DBL/output_dir definition
IMPLICIT NONE
INTEGER::IJKL,NTOT,NTOT2,I97,J97
REAL(KIND(1d0))::C
REAL(KIND(1d0)),DIMENSION(97)::U
SAVE::IJKL,NTOT,NTOT2,I97,J97,C,U
LOGICAL::FIRST=.TRUE.
SAVE  FIRST
CONTAINS

FUNCTION Helac_rnmy(idummy)
INTEGER,INTENT(IN)::idummy
REAL(KIND(1d0)),DIMENSION(1)::RRR
REAL(KIND(1d0))::Helac_rnmy
INTEGER::L
L=1
CALL Helac_RANMAR(RRR,L)
Helac_rnmy=RRR(1)

!Helac_rnmy=Random_irn(0)*(1.0d-9)
!IF(Helac_rnmy.EQ.0.0d0)Helac_rnmy=Random_irn(0)*(1.0d-9)
!RETURN
END FUNCTION Helac_rnmy
!
! $Id: ranmar.F,v 1.1.1.1 1996/02/15 17:49:53 mclareni Exp $
!
! $Log: ranmar.F,v $
! Revision 1.1.1.1  1996/02/15 17:49:53  mclareni
! Kernlib
!
!
SUBROUTINE Helac_RANMAR(RVEC,LENV)
!
! CERN PROGLIB# V113    RANMAR          .VERSION KERNFOR  4.21  890323
! ORIG. 01/03/89 FCA + FJ
!
IMPLICIT NONE
INTEGER,INTENT(IN)::LENV
REAL(KIND(1d0)),DIMENSION(*),INTENT(OUT)::RVEC
REAL(KIND(1d0)),PARAMETER::TWOM24=2d0**(-24),TWOM48=2d0**(-48)
REAL(KIND(1d0)),PARAMETER::CD=7654321.*TWOM24,CM=16777213.*TWOM24
REAL(KIND(1d0)),PARAMETER::CINT=362436.*TWOM24
INTEGER,PARAMETER::MODCNS=1000000000
INTEGER::iseed,IVEC,IJ,KL,I,J,K,L,II,JJ,M,NITER,LOOP2,IDUM
REAL(KIND(1d0))::UNI,S,T

DO
  IF(FIRST) THEN
     OPEN(77,FILE=TRIM(input_dir)//'seed.input',status='OLD') 
     READ(77,*)iseed
     WRITE(*,*)'Random number seed',iseed
     IJKL = 54217137+10*iseed
     CLOSE(77)
     NTOT = 0
     NTOT2 = 0
!        GO TO 70
   ELSE
!  80 CONTINUE
     DO IVEC= 1, LENV      !100
        UNI = U(I97)-U(J97)
        IF (UNI.LT.0.) UNI=UNI+1.d0
        U(I97) = UNI
        I97 = I97-1
        IF (I97.EQ.0)I97=97
        J97 = J97-1
        IF (J97.EQ.0)J97=97
        C = C - CD
        IF (C.LT.0.)C=C+CM
        UNI = UNI-C
        IF (UNI.LT.0.)UNI=UNI+1.d0
!
!   Replace exact zeroes by uniform distr. *2**-24
!
        IF (UNI.EQ.0.)THEN
           UNI = TWOM24*U(2)
!
!   An exact zero here is very unlikely, but let's be safe.
!
           IF (UNI.EQ.0.) UNI= TWOM48
        ENDIF
        RVEC(IVEC) = UNI
     ENDDO
     NTOT = NTOT + LENV
     IF (NTOT.GE.MODCNS)THEN
        NTOT2 = NTOT2 + 1
        NTOT  = NTOT - MODCNS
     ENDIF
     RETURN
  ENDIF
!   70 CONTINUE
  IJ = IJKL/30082
  KL = IJKL - 30082*IJ
  I = MOD(IJ/177, 177) + 2
  J = MOD(IJ, 177)     + 2
  K = MOD(KL/169, 178) + 1
  L = MOD(KL, 169)
  DO II= 1, 97
     S = 0.0d0
     T = 0.5d0
     DO JJ= 1, 24
         M = MOD(MOD(I*J,179)*K, 179)
         I = J
         J = K
         K = M
         L = MOD(53*L+1, 169)
         IF (MOD(L*M,64) .GE. 32)  S = S+T
         T = 0.5d0*T
      ENDDO
      U(II) = S
  ENDDO
  C   = CINT
  I97 = 97
  J97 = 33
!      Complete initialization by skipping
!            (NTOT2*MODCNS + NTOT) random numbers
  NITER = MODCNS
  DO LOOP2= 1, NTOT2+1
     IF(LOOP2.GT.NTOT2) NITER=NTOT
     DO IDUM = 1, NITER
          UNI = U(I97)-U(J97)
          IF (UNI .LT. 0.) UNI=UNI+1.
          U(I97) = UNI
          I97 = I97-1
          IF (I97 .EQ. 0)  I97=97
          J97 = J97-1
          IF (J97 .EQ. 0)  J97=97
          C = C - CD
          IF (C .LT. 0.)   C=C+CM
      ENDDO
  ENDDO
  NTOT  = 0
  NTOT2 = 0
  IF(FIRST) THEN
      FIRST = .FALSE.
!        GO TO 80
  ELSE
      EXIT
  ENDIF
ENDDO
END SUBROUTINE Helac_RANMAR

SUBROUTINE Helac_RMARIN(IJKLIN,NTOTIN,NTO2IN)
IMPLICIT NONE
INTEGER,INTENT(IN)::IJKLIN,NTOTIN,NTO2IN
REAL(KIND(1d0)),PARAMETER::TWOM24=2d0**(-24),TWOM48=2d0**(-48)
REAL(KIND(1d0)),PARAMETER::CD=7654321.d0*TWOM24,CM=16777213.d0*TWOM24
REAL(KIND(1d0)),PARAMETER::CINT=362436.d0*TWOM24
INTEGER,PARAMETER::MODCNS=1000000000
INTEGER::IVEC,IJ,KL,I,J,K,L,II,JJ,M,NITER,LOOP2,IDUM
REAL(KIND(1d0))::UNI,S,T
FIRST = .FALSE.
IJKL  = IJKLIN
NTOT  = NTOTIN
NTOT2 = NTO2IN
!   70 CONTINUE
IJ = IJKL/30082
KL = IJKL - 30082*IJ
I = MOD(IJ/177, 177) + 2
J = MOD(IJ, 177)     + 2
K = MOD(KL/169, 178) + 1
L = MOD(KL, 169)
DO II= 1, 97
   S = 0.0d0
   T = 0.5d0
   DO JJ= 1, 24
      M = MOD(MOD(I*J,179)*K, 179)
      I = J
      J = K
      K = M
      L = MOD(53*L+1, 169)
      IF (MOD(L*M,64) .GE. 32)  S = S+T
      T = 0.5d0*T
   ENDDO
   U(II) = S
ENDDO
C   = CINT
I97 = 97
J97 = 33
!      Complete initialization by skipping
!            (NTOT2*MODCNS + NTOT) random numbers
NITER = MODCNS
DO LOOP2= 1, NTOT2+1
   IF(LOOP2.GT.NTOT2) NITER=NTOT
   DO IDUM = 1, NITER
      UNI = U(I97)-U(J97)
      IF (UNI .LT. 0.) UNI=UNI+1.d0
      U(I97) = UNI
      I97 = I97-1
      IF (I97 .EQ. 0)  I97=97
      J97 = J97-1
      IF (J97 .EQ. 0)  J97=97
      C = C - CD
      IF (C .LT. 0.)   C=C+CM
  ENDDO
ENDDO
NTOT  = 0
NTOT2 = 0
END SUBROUTINE Helac_RMARIN

SUBROUTINE Helac_RMARUT(IJKLUT,NTOTUT,NTO2UT)
IMPLICIT NONE
INTEGER,INTENT(OUT)::IJKLUT,NTOTUT,NTO2UT
NTOTUT = NTOT
NTO2UT = NTOT2
IJKLUT = IJKL
END SUBROUTINE Helac_RMARUT

!     function rnmy(idummy)
!      include 'declare.h'

!     rnmy=irn(0)*dnou(10)**(-9)
!      if(rnmy.eq.0.0d0)rnmy=irn(0)*dnou(10)**(-9)
!     return
!       end
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Following functions are from Random.f
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     
FUNCTION Random_irn(idummy)
IMPLICIT NONE
INTEGER,INTENT(IN)::idummy
INTEGER,DIMENSION(55)::ia_random
INTEGER::j_random=55,Random_irn
SAVE ia_random,j_random 
      
j_random=j_random+1
IF(j_random.GT.55) j_random=Random_irn55(ia_random)
Random_irn=ia_random(j_random)
RETURN
END FUNCTION Random_irn
      
SUBROUTINE Random_in55(ia_random,ix_random)
IMPLICIT NONE
INTEGER,DIMENSION(55),INTENT(OUT)::ia_random
INTEGER,INTENT(IN)::ix_random
INTEGER::j_random,k_random,i_random,ii_random,iw_random
ia_random(55)=ix_random
j_random=ix_random
k_random=1
DO i_random=1,54
	ii_random=MOD(21*i_random,55)
	ia_random(ii_random)=k_random
	k_random=j_random-k_random
	IF(k_random.LT.0)k_random=k_random+1000000000
	j_random=ia_random(ii_random)
ENDDO
iw_random=Random_irn55(ia_random)
iw_random=Random_irn55(ia_random)
iw_random=Random_irn55(ia_random)
END SUBROUTINE Random_in55
      
FUNCTION Random_irn55(ia_random)
INTEGER,DIMENSION(55),INTENT(OUT)::ia_random
INTEGER::init=0,Random_irn55,iseed_random,ix_random,i_random,j_random
SAVE init
      
IF(init.EQ.0)THEN
	init=1
	OPEN(77,FILE=TRIM(input_dir)//'seed.input',STATUS='OLD') 
	READ(77,*)iseed_random
	PRINT *,'Random number seed',iseed_random
	ix_random=6655070+10*iseed_random
	CLOSE(77)
	CALL Random_in55(ia_random,ix_random)
ENDIF
      
DO i_random=1,24
	j_random=ia_random(i_random)-ia_random(i_random+31)
	IF(j_random.LT.0)j_random=j_random+1000000000
	ia_random(i_random)=j_random
ENDDO
DO i_random=25,55
	j_random=ia_random(i_random)-ia_random(i_random-24)
	IF(j_random.LT.0) j_random=j_random+1000000000
	ia_random(i_random)=j_random
ENDDO
Random_irn55=1
RETURN
END FUNCTION Random_irn55
END MODULE Helac_ranmar_mod
