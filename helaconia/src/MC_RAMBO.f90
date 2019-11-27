MODULE MC_RAMBO
USE Helac_ranmar_mod
IMPLICIT NONE
CONTAINS
SUBROUTINE RAMBO(N,ET,XM,P,WT,IDIR)
!------------------------------------------------------
!
!                       RAMBO
!
!    RA(NDOM)  M(OMENTA)  B(EAUTIFULLY)  O(RGANIZED)
!
!    A DEMOCRATIC MULTI-PARTICLE PHASE SPACE GENERATOR
!    AUTHORS:  S.D. ELLIS,  R. KLEISS,  W.J. STIRLING
!    THIS IS VERSION 1.0 -  WRITTEN BY R. KLEISS
!
!    N  = NUMBER OF PARTICLES (>1, IN THIS VERSION <101)
!    ET = TOTAL CENTRE-OF-MASS ENERGY
!    XM = PARTICLE MASSES ( DIM=100 )
!    P  = PARTICLE MOMENTA ( DIM=(4,100) )
!    WT = WEIGHT OF THE EVENT
!
!------------------------------------------------------
INTEGER,INTENT(IN)::N
REAL(KIND(1d0)),INTENT(IN)::ET
REAL(KIND(1d0)),DIMENSION(100),INTENT(IN)::XM
REAL(KIND(1d0)),DIMENSION(4,100),INTENT(INOUT)::P
REAL(KIND(1d0)),INTENT(OUT)::WT
INTEGER,INTENT(IN)::IDIR
REAL(KIND(1d0)),DIMENSION(4,100)::Q
REAL(KIND(1d0)),DIMENSION(100)::RZ,P2,XM2,E,V
REAL(KIND(1d0)),DIMENSION(4)::R,RN
REAL(KIND(1d0)),DIMENSION(3)::B
INTEGER,DIMENSION(5)::IWARN=(/0,0,0,0,0/)
INTEGER::ITMAX=10,IBEGIN=0,K,I,J
INTEGER::IVERSION=0,NM,ITER
REAL(KIND(1d0))::TWOPI,Costh,Sinth,Fphi,G,A,XMAX,X,WTM,WT2,&
WT3,F0,G0,X2,ACCU,BQ,RMAS,XMT,PI2LOG,ACC
SAVE ACC,ITMAX,IBEGIN,TWOPI,PI2LOG,IWARN,RZ
!------------------------------------------------------
!       VERSION DEFINITION
!------------------------------------------------------
IVERSION=1
!------------------------------------------------------
!         IBEGIN=0
! INITIALIZATION STEP: FACTORIALS FOR THE PHASE SPACE WEIGHT
IF(IBEGIN.EQ.0)THEN                !GOTO 103
   ACC=10d0*epsilon(1d0)   ! the precesion
   IBEGIN=1
   TWOPI=8d0*DATAN(1d0)  ! 2*pi
   PI2LOG=DLOG(TWOPI/4d0)    ! Log[pi/2]
   RZ(2)=PI2LOG
   DO K=3,100
      RZ(K)=RZ(K-1)+PI2LOG-2d0*DLOG(DBLE(K-2))
   ENDDO
   DO K=3,100
      RZ(K)=(RZ(K)-DLOG(DBLE(K-1)))
   ENDDO
ENDIF
! CHECK ON THE NUMBER OF PARTICLES
! 103
IF(N.LE.1.OR.N.GE.101)THEN         !GOTO 104       
   WRITE(*,1001)N
   STOP
ENDIF
!
! CHECK WHETHER TOTAL ENERGY IS SUFFICIENT; COUNT NONZERO MASSES
! 104 
XMT=0d0
NM=0
DO I=1,N
   IF((XM(I).GT.0d0).OR.(XM(I).LT.0d0)) NM=NM+1
   XMT=XMT+DABS(XM(I))
ENDDO
IF(XMT.GT.ET)THEN                     !GOTO 201
   WRITE(*,1002)XMT,ET
   STOP
ENDIF

! THE PARAMETER VALUES ARE NOW ACCEPTED
!
! GENERATE N MASSLESS MOMENTA IN INFINITE PHASE SPACE
! 201  CONTINUE
IF(IDIR.NE.1)THEN             !GOTO 401 !idir
   DO I=1,N
      DO J=1,4
        RN(J)=Helac_rnmy(0)
      ENDDO
      Costh=2d0*RN(1)-1d0   ! Costh
      Sinth=DSQRT(1d0-Costh*Costh)  ! Sinth
      Fphi=TWOPI*RN(2)     ! phi
      Q(4,I)=-DLOG(RN(3)*RN(4))  ! generate by probability q0 exp(-q0) dq0
	                             ! see Appendix(c)  in Comp.Phys.Commun   40 (1986) 359
      Q(3,I)=Q(4,I)*Costh
      Q(2,I)=Q(4,I)*Sinth*DSIN(Fphi)
      Q(1,I)=Q(4,I)*Sinth*DCOS(Fphi)   ! 202
  ENDDO
!
! CALCULATE THE PARAMETERS OF THE CONFORMAL TRANSFORMATION
  DO  I=1,4
      R(I)=0d0
  ENDDO
  DO  I=1,N
      DO  K=1,4
         R(K)=R(K)+Q(K,I)   ! 204
      ENDDO
  ENDDO
  RMAS=DSQRT(R(4)**2-R(3)**2-R(2)**2-R(1)**2)
  DO  K=1,3
      B(K)=-R(K)/RMAS       ! 205
  ENDDO
  G=R(4)/RMAS
  A=1d0/(1d0+G)
  X=ET/RMAS
!
! TRANSFORM THE Q'S CONFORMALLY INTO THE P'S
! Lorentz Boost to the CM of the sum of final states and scale a factor ET/RMAS
  DO  I=1,N
      BQ=B(1)*Q(1,I)+B(2)*Q(2,I)+B(3)*Q(3,I)
      DO K=1,3
         P(K,I)=X*(Q(K,I)+B(K)*(Q(4,I)+A*BQ))  ! 206
      ENDDO
      P(4,I)=X*(G*Q(4,I)+BQ)                   ! 207
  ENDDO
ENDIF
!  401 CONTINUE !idir
!
! CALCULATE WEIGHT AND POSSIBLE WARNINGS
WT=PI2LOG
IF(N.NE.2) WT=(2d0*N-4d0)*DLOG(ET)+RZ(N)
IF(WT.LT.-18d0)THEN               ! GOTO 208
  IF(IWARN(1).LE.5) WRITE(*,1004)WT
  IWARN(1)=IWARN(1)+1
ENDIF
IF(WT.GT.174d0)THEN                 !GOTO 209   ! 208
  IF(IWARN(2).LE.5)WRITE(*,1005)WT
  IWARN(2)=IWARN(2)+1
ENDIF

! RETURN FOR WEIGHTED MASSLESS MOMENTA
IF(NM.EQ.0)THEN  ! all the masses are zero      !GOTO 210        ! 209
  WT=DEXP(WT)  ! uniform weight
  RETURN
ENDIF

! MASSIVE PARTICLES: RESCALE THE MOMENTA BY A FACTOR X
XMAX=DSQRT(1d0-(XMT/ET)**2)    ! 210

IF(IDIR.NE.1)THEN                     !GOTO 402  !idir
  DO I=1,N
     XM2(I)=XM(I)**2 
     P2(I)=P(4,I)**2              !301
  ENDDO
  ITER=0
  X=XMAX
  ACCU=ET*ACC
  DO                        ! 302
     F0=-ET      
     G0=0d0
     X2=X*X
     DO I=1,N
        E(I)=DSQRT(XM2(I)+X2*P2(I))
        F0=F0+E(I)
        G0=G0+P2(I)/E(I)        ! 303
     ENDDO
     IF(DABS(F0).LE.ACCU)EXIT      ! GOTO 305
     ITER=ITER+1
     IF(ITER.GT.ITMAX)THEN        ! GOTO 304
        WRITE(*,1006)ITMAX
        EXIT                      ! GOTO 305
	 ENDIF
     X=X-F0/(X*G0)                ! 304   !  Newton-Raphson iteration 
	                                      ! see Nucl.Phys.B 385 (1992) 413-432
                                  !GOTO 302
  ENDDO
  DO  I=1,N
      V(I)=X*P(4,I)
      DO  K=1,3
          P(K,I)=X*P(K,I)
	  ENDDO
      P(4,I)=E(I)
  ENDDO
ENDIF
!  402 CONTINUE  ! idir
DO I=1,N
   V(I)=SQRT(P(1,I)**2+P(2,I)**2+P(3,I)**2)
   E(I)=P(4,I)
ENDDO
X=0
DO I=1,N
   X=X+V(I)/ET
ENDDO

! CALCULATE THE MASS-EFFECT WEIGHT FACTOR
WT2=1d0
WT3=0d0
DO I=1,N
   WT2=WT2*V(I)/E(I)
   WT3=WT3+V(I)**2/E(I)
ENDDO
WTM=(2d0*N-3d0)*DLOG(X)+DLOG(WT2/WT3*ET)

! RETURN FOR  WEIGHTED MASSIVE MOMENTA
WT=WT+WTM
IF(WT.LT.-18d0)THEN              ! GOTO 309
   IF(IWARN(3).LE.5) WRITE(*,1004)WT,EXP(WT)
   IWARN(3)=IWARN(3)+1
ENDIF
IF(WT.GT. 174d0)THEN       !GOTO 310   ! 309
   IF(IWARN(4).LE.5) WRITE(*,1005)WT,EXP(WT)
   IWARN(4)=IWARN(4)+1
ENDIF
WT=DEXP(WT)     ! 310
RETURN

1001 FORMAT(' RAMBO FAILS: # OF PARTICLES =',I5,' IS NOT ALLOWED')
1002 FORMAT(' RAMBO FAILS: TOTAL MASS =',D15.6,' IS NOT',&
     ' SMALLER THAN TOTAL ENERGY =',D15.6)
1004 FORMAT('RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY UNDERFLOW',D20.9)
1005 FORMAT('RAMBO WARNS: WEIGHT = EXP(',F20.9,') MAY  OVERFLOW',D20.9)
1006 FORMAT('RAMBO WARNS:',I3,' ITERATIONS DID NOT GIVE THE',&
      ' DESIRED ACCURACY =',D15.6)
END SUBROUTINE RAMBO
END MODULE MC_RAMBO
