* QEDPS program library V1.0   1997/Feb./10
* Originally coded by T. Munehisa (Yamanashi UNIV.)
* Modified by Y. Kurihara (KEK)
         SUBROUTINE QPINIT(qq2max,eemass,iseed)
C
	 IMPLICIT REAL*8 (A-H,O-Z)
	 include 'qpincl.f'
         COMMON /QPRND / NSEED
C                                                                               
        PI = acos(-1.d0)
	ALPHA = 1.0D0/137.0359895D0
        IBR  = 0
        NPH = -1                                                               
        NALP = 0                                                               
	if(iseed.gt.0) then
         NSEED = iseed
	else
         NSEED = 1805113725
	end if
                                                                                
        Q2MAX = qq2max
        q2min = 0.d0
                                                                                
        EMASS = eemass
        PHMIN = 2.d-4
        Q02 = (0.1D-5)**2   ! Q0^2                                                    
        Q02M = emass**2*exp(1.d0) ! Qs^2
c       Q02M = (0.1D-2)**2
* PROGRAM TITLE 'QEDPS'.
*=== LOGO
      ISRC = 6

      WRITE(ISRC,100)
  100 FORMAT(
     . '*************************************************************' 
     . '***********'/
     . '*------------------------------------------------------------' 
     . '----------*'/
     . '*------------------------------------------------------------' 
     . '----------*'/
     . '*-----QQQQQQ------EEEEEEEEEEE--DDDDDDDDD------PPPPPPPPP------' 
     . '--SSSS-S--*'/
     . '*---QQQQ--QQQQ-----EEE----EEE---DDD---DDDD-----PPP---PPPP----' 
     . 'SSS--SSS--*')

      WRITE(ISRC,101)
  101 FORMAT(
     . '*--QQQ------QQQ----EEE-----EE---DDD-----DDD----PPP----PPPP--S' 
     . 'SS----SS--*'/
     . '*--QQQ------QQQ----EEE---E--E---DDD-----DDD----PPP----PPPP--S' 
     . 'SS-----S--*'/
     . '*-QQQ--------QQQ---EEE---E--E---DDD------DDD---PPP----PPPP--S' 
     . 'SSS-------*'/
     . '*-QQQ--------QQQ---EEE--EE------DDD------DDD---PPP----PPPP--S' 
     . 'SSSSS-----*'/
     . '*-QQQ--------QQQ---EEEEEEE------DDD------DDD---PPP---PPPP----' 
     . 'SSSSSSS---*')

      WRITE(ISRC,102)
  102 FORMAT(
     . '*-QQQ--------QQQ---EEE--EE------DDD------DDD---PPPPPPPP------' 
     . '--SSSSSS--*'/
     . '*-QQQ--------QQQ---EEE---E------DDD------DDD---PPP-----------' 
     . '----SSSSS-*'/
     . '*-QQQ-QQQQ---QQQ---EEE---E--E---DDD------DDD---PPP----------S' 
     . '-----SSSS-*'/
     . '*--QQQ--QQQ--QQ----EEE------E---DDD-----DDD----PPP----------S' 
     . '------SSS-*'/
     . '*--QQQ---QQ-QQQ----EEE-----EE---DDD-----DDD----PPP----------S' 
     . 'S-----SSS-*')

      WRITE(ISRC,103)
  103 FORMAT(
     . '*---QQQQ-QQQQQ-----EEE----EEE---DDD---DDDD-----PPP----------S' 
     . 'SSS--SSS--*'/
     . '*-----QQQQQQ------EEEEEEEEEEE--DDDDDDDDD------PPPPP---------S' 
     . '--SSSSS---*'/
     . '*---------QQ-------------------------------------------------' 
     . '----------*'/
     . '*---------QQQ-Q----------------------------------------------' 
     . '----------*'/
     . '*---------QQQQQ--------------  Version 1.0  -----  1997 Febur' 
     . 'ary  -----*')

      WRITE(ISRC,104)
  104 FORMAT(
     . '*----------QQQ-------  (c)   Minami Tateya Collaboration,  JA' 
     . 'PAN  -----*'/
     . '*------------------------------------------------------------' 
     . '----------*'/
     . '* Author   :T. Munehisa, J. Fujimoto, Y. Kurihara, and  Y. Sh'
     . 'imizu     *'/
     . '*************************************************************' 
     . '***********'/
     . '* Reference: Y. Kurihara, J. Fujimoto, T. Munehisa, and Y. Sh'
     . 'imizu     *'/
     . '* "Hard Photon Distributions in e+e- Annihilation Processes b'
     . 'y QEDPS"  *'/
     . '* Progress of Theoretical Physics, Vol.96 (1996) 1223,       ' 
     . '          *')

      WRITE(ISRC,105)
  105 FORMAT(
c    . '*       xxxxxxxx, xxxxxxx, xxxxxxxxxxxxxxxxxxxxxxxxxxx       ' 
c    . '          *'/
     . '*************************************************************' 
     . '***********')

*=== END OF LOGO
C                                                                               
        call qpprep
c
        RETURN                                                                  
        END                                                                     
        SUBROUTINE QPPREP
C                                                                               
C        Calculation of Sudakov form factor                                     
C  exp(  -\int^{Q2max}_{Q02} dK^2/K^2\int^{1-Q02/K^2}_0dx alpha/2pi P(x) )      
C  exp(  -\int^{Q2max}_{Q02} dK^2/K^2\int^{1-Q02/K^2}_0dx alpha/2pi P(x) )      
                                                                                
       IMPLICIT REAL*8(A-H,O-Z)                                               
       include 'qpincl.f'
                                                                                
C         P1 =2.0*LOG(K2/Q02)-3.0/2.0= \int^{1-Q02/K2}_0dx P(x),P(x)=(1+x^2)/(1-x)                                           
C         Taking account of effect of the running coupling                      
C         1992.5.7                                                              
C         ALPHA(Q**2)=ALPHA/(1-ALPHA*C*LOG(Q**2/Q02))                           
C         C = 1/pi/3:only electron loop                                         
C         P2 = 1/2/PI/C*(-2*(LOG(Q2MAX/Q02)-1)                                  
C             +(2/ALPHA/C-3/2)*(-LOG(1-ALPHA*C*LOG(Q2MAX/Q02))                  
C                               +LOG(1-ALPHA*C)  )), next order of P(x) with running of ALPHA                             
C                                                                               
	 C= 1.0/3.0/PI                                                         
          NSDK = 1000
          TMAX = LOG(Q2MAX/Q02)                                                 
          TMIN = 0.0                                                            
           TD = (TMAX - TMIN ) / NSDK                                           
	DO 101 I=1,NSDK                                                       
           T = TMIN + TD*I                                                      
      IF( NALP .EQ. 0 )THEN                                                     
          P2 = 2.0*1.0/2.0*T**2-3.0/2.0*T                                       
          P2 = P2 * ALPHA /2.0/PI                                               
                       ELSE                                                     
         P2=1.0/2.0/PI/C*(  -2.0*(T-1)                                          
     & +(2.0/ALPHA/C-3.0/2.0)*(-LOG(1.0-ALPHA*C*T)                              
     & +LOG(1.0-ALPHA*C)  )  )                                                  
                       ENDIF                                                    
                                                                                
          PSDK(I) = -P2                                                         
                                                                                
  101     CONTINUE                                                              
                                                                                
	 RETURN                                                                
          END                                                                   
          SUBROUTINE QPGEN(q2in,  q2out)

         IMPLICIT REAL*8(A-H,O-Z)
        COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP
	include 'qpincl.f'

         do 201 l=1,1000
	 do 201 i =1,10
	 plptn(i,l) = 0.0d0
	 nlptn(i,l) = 0
 201     continue
C
          ECM = sqrt(q2in)
	  s   =      q2in
          NOK = 0
  5       CONTINUE
          CALL QPSINT


          CALL QPANNH(Q2,NOK)
	  q2out=q2
C

C     if NOK is not zero, the event is accepted.
           IF( NOK .EQ. 0 )goto 5
          NVPH = NTOP ! NTOP is total number of particles
C
      IF ( NPH .GE. 0 )THEN
        NNPH = 0 ! number of photon from ISR
         DO 121 I=1,NTOP
      IF( NLPTN(1,I) .EQ. 22   .AND.  NLPTN(3,I).EQ.0)NNPH=NNPH+1
 121   CONTINUE
       IF(  NNPH  .GT. NPH)THEN
                            NOK = 0
                            GOTO 5
                           ENDIF
                       ENDIF

           RETURN
           END

          SUBROUTINE QPSINT
         IMPLICIT REAL*8(A-H,O-Z)

         include 'qpincl.f'
         COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP

C         PLPTN(1, ):  The lightcone fraction
C         PLPTN(2, ):  virtual mass squared
C         PLPTN(3, ):  Pt**2
C         PLPTN(4, ):  x-component of the four momentum
C         PLPTN(5, ):  y-component of the four momentum
C         PLPTN(6, ):  z-component of the four momentum
C         PLPTN(7, ):  E-component of the four momentum
C         PLPTN(8, ):  + component of the lightcone momentum
C         PLPTN(9, ):  - component of the lightcone momentum

C         NLPTN(1, ): particle ID
C         NLPTN(2, ): relative address of the parent
C         NLPTN(3, ): number of the children
C         NLPTN(4, ): relative address of the first child
C         NLPTN(5, ): status of the process
C         NLPTN(6, ): spacelike(-1) or timelike(+1)

C
         ! initialize the beam particles
          NTOP=0
          DO   1  LID =1,2
          IF( LID .EQ. 1) THEN
                            ID= 11 ! make sure the first beam is e-
                          ELSE
                            ID = -11 ! make sure the second beam is e+
                          ENDIF

          NTOP = NTOP + 1
          L = NTOP
          PLPTN(1,L) = 1.0D0
          PLPTN(2,L) = EMASS**2
          PLPTN(3,L) = 0.0D0
          PLPTN(4,L) = 0.0D0
          PLPTN(5,L) = 0.0D0
          IF ( ID. EQ. 11 )THEN
          PLPTN(6,L) = ECM/2.0D0
                           ELSE
          PLPTN(6,L) = -ECM/2.0D0
                           ENDIF
          PLPTN(7,L) = ECM/2.0D0
          PLPTN(8,L) = ECM/SQRT(2.0D0)
          PLPTN(9,L) = 0.0D0

          NLPTN(1,L) = ID
          NLPTN(2,L) = 0
          NLPTN(5,L) = 0
          NLPTN(6,L) = -1 ! -1 spacelike or initial state
  1       CONTINUE

ccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccc    improved QEDPS    1995.8.8
cccccc    start   ccccccccccccccccc
           nccc =1
          if( nccc .eq. 1)then
          DO   4  L =1,2
          call   QPcorr(xc)
        if( xc .ne. 1.0d0)then
ccc         write(6,*)' improved !! '
         lx = ntop+1
         ntop = ntop+2
          PLPTN(1,Lx) = xc*plptn(1,l)
          PLPTN(2,Lx) = 0.0d0
          PLPTN(3,Lx) = 0.0D0
          PLPTN(4,Lx) = 0.0D0
          PLPTN(5,Lx) = 0.0D0
          plptn(6,lx) =plptn(6,l)*xc
          plptn(7,lx) =plptn(7,l)*xc
          plptn(8,lx) =plptn(8,l)*xc
          NLPTN(1,Lx) = nlptn(1,l)
          NLPTN(2,Lx) = lx-l
          NLPTN(5,Lx) = 0
          NLPTN(6,Lx) = -1
ccccc
c
          NLPTN(3,L) = 2
          NLPTN(4,L) = lx-l
          NLPTN(5,L) = 1
ccccc
          lx = lx+1
           xc1 = 1.0d0 - xc
          PLPTN(1,Lx) = xc1*plptn(1,l)
          PLPTN(2,Lx) = 0.0d0
          PLPTN(3,Lx) = 0.0D0
          PLPTN(4,Lx) = 0.0D0
          PLPTN(5,Lx) = 0.0D0
          plptn(6,lx) =plptn(6,l)*xc1
          plptn(7,lx) =plptn(7,l)*xc1
          plptn(8,lx) =plptn(8,l)*xc1
          NLPTN(1,Lx) = 22
          NLPTN(2,Lx) = 0
          NLPTN(5,Lx) = 2
          NLPTN(6,Lx) = 1
                          endif
  4        continue
                                    endif
cccccc    end   ccccccccccccccccc
ccccccccccc    improved QEDPS    1995.8.8 ccccccccccc
cccccccccccccccccccccccccccccccccccccccccccccccc
C    end of  momentum assign for  initial elecrton or positron

C      branching

          NAD = 1

        IF( IBR .NE. 1 )THEN
   2      CONTINUE
          IF( NLPTN(5,NAD) .EQ. 0) THEN
             ! ISR photon shower
          CALL QPGEN1(NAD)
                                   ENDIF
           NAD = NAD + 1
          IF( NAD .LE. NTOP) GOTO 2

                        ENDIF

         DO 3 L=1,NTOP
         IF( NLPTN(5,L) .EQ. 2)GOTO 3
         ! generate momenta
         CALL QPCMSN(L)
 3       CONTINUE

C      WRITE(6,111)
 111   FORMAT('   END OF PSINT ')
         RETURN
         END

          SUBROUTINE QPGEN1(NAD)
         IMPLICIT REAL*8(A-H,O-Z)
         COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP
         include 'qpincl.f'
         REAL*8 K2


C      start of  parton shower

         Q2MX = Q2MAX / 2.0D0
         Q2MN = MAX( -PLPTN(2,NAD) , q02m) ! q02m=Qs^2
         q2mnc = emass**2
         ID   = NLPTN(1,NAD)

  5      CONTINUE
         T1 = LOG(Q2MX / Q02 ) ! Q02=Q0^2,eps=Q02/Q2MAX
         T2 = LOG(Q2MN / Q02 )

         dt =(tmax-tmin)/dfloat(nsdk) ! nsdk=1000,tmax=LOG(Q2MAX/Q02) 
         LT1 = (T1 - TMIN )/dt ! TMIN=0.0
         LT2 = (T2 - TMIN )/dt
         if(lt2.eq.0)lt2=1
         IF( LT1 .LE. LT2 ) GOTO 101
         ! PSDK(LT1)=-\int^{K2}_{Q02} dK^2/K^2\int^{1-Q02/K2}_0dx alpha/2pi P(x),K2 LT1 value in [TMIN,TMAX]
         ! sudakov grid
      psdk1 =  PSDK(LT1)+(psdk(lt1+1)-psdk(lt1))*(t1-(lt1*dt+tmin))/dt
      psdk2 =  PSDK(LT2)+(psdk(lt2+1)-psdk(lt2))*(t2-(lt2*dt+tmin))/dt
          p0 = psdk1 - psdk2
          P1 = LOG( QPRAND(1) )
C
          IF ( P1 .LE. P0 ) GOTO 101 ! there is no emission of photon eta < sudakov
           P2 = ( P1 ) + psdk2

          DO 201 LT =LT2, LT1
          IF( PSDK( LT ) .LE. P2 )GOTO 202
  201     CONTINUE

  202     CONTINUE

C         branching is accepted

          CLT = LT-1 +(P2 -PSDK( LT-1) )/(PSDK(LT)-PSDK(LT-1))
          T = CLT*((TMAX-TMIN)/NSDK ) + TMIN
          K2 = exp( log(q02)+ T ) ! determine k2 from eta=sudakov(k2,Q02)
         IF( K2 .LE. Q02)THEN
            WRITE(6,*)'LT1,LT2         = ',LT1,LT2
            WRITE(6,*)'LT,CLT,T,Q02,P2 = ',LT,CLT,T,Q02,P2
            WRITE(6,*)'PSDK(LT),PSDK(LT-1)=',PSDK(LT),PSDK(LT-1)
                             ENDIF
         IF( ABS(Q02/K2). LE. 1.E-20)THEN
            WRITE(6,*)'Q02,K2          = ',Q02,K2
            WRITE(6,*)'LT1,LT2         = ',LT1,LT2
            WRITE(6,*)'LT,CLT,T,Q02,P2 = ',LT,CLT,T,Q02,P2
            WRITE(6,*)'PSDK(LT),PSDK(LT-1)=',PSDK(LT),PSDK(LT-1)
                             ENDIF

C         determine X according to P(x)=(1+x^2)/(1-x)
C         int^{1-eps}_0dx P(x)=int^{log(eps)}_0dx1l (1+(1-exp(x1l))^2)
          XMAX = 1.0D0 - Q02/K2 ! upper limit of x
          XMIN = 0.0D0 ! lower limit of x
          X1LMX=LOG(q02/k2) ! log(eps)
          X1LMN=LOG(1.0D0 - XMIN) ! log(1)
           X1LD = X1LMX - X1LMN
           PXM =2.0D0
           DO 203 I=1,200

            X1L = X1LMN + X1LD * QPRAND(1)

             X = 1.0D0 - EXP(X1L)
             X1=  EXP(X1L)
         IF( ABS(X) .GE. 1.1)THEN
            WRITE(6,*)'X, XMAX,X1LD,K2=',X,XMAX,X1LD,K2
            WRITE(6,*)'LT,CLT,T,Q02 = ',LT,CLT,T,Q02
            WRITE(6,*)'PSDK(LT),PSDK(LT-1)=',PSDK(LT),PSDK(LT-1)
                             ENDIF

          PX1= 1.0D0 + X**2

          IF( PX1 .GE. PXM*QPRAND(1) )GOTO 204 ! adapt distribution via PX1 form (hit x to P(x))
 203      CONTINUE
c          write(6,*)' fail to hit x to P(x)'

 204      CONTINUE

C         check the constraint is satisfied or not

           eta =1.0d0
           ! x1 is 1-x
           IF(  X1 .LT. eta*K2/Q2MAX ) THEN
                                         Q2MN = K2 ! there is no emission up to K2
                                        GOTO 5
                                        ENDIF

C         write PLPTN and NLPTN of the children

          L = NTOP+1
          PLPTN(1,L) = X*PLPTN(1,NAD)
          PLPTN(2,L) = - K2

          NLPTN(1,L) = ID
          NLPTN(2,L) = L - NAD
          NLPTN(5,L) = 0
          NLPTN(6,L) = -1

          L = NTOP+2
          PLPTN(1,L) =X1*PLPTN(1,NAD)
          PLPTN(2,L) = 0.0D0

C        photon does not branch into electron pair (Now, just appx.)
          NLPTN(1,L) = 22
          NLPTN(2,L) = L - NAD
          NLPTN(3,L) = 0
          NLPTN(5,L) = 1
          NLPTN(6,L) = 1

          NLPTN(3,NAD) = 2
          NLPTN(4,NAD) = (L-1) - NAD
          NLPTN(5,NAD) = 1

           NTOP = NTOP + 2

C      WRITE(6,111)
          RETURN

C         not branch
 101      CONTINUE
C       WRITE(6,*)' NOT BRANCH'

          NLPTN(3,NAD) = 0
          NLPTN(4,NAD) = 0
          NLPTN(5,NAD) = 1
C      WRITE(6,111)
 111   FORMAT('   END OF PSGEN1')

          RETURN
          END

          SUBROUTINE QPCMSN(NAD)
          ! generate momenta
         IMPLICIT REAL*8(A-H,O-Z)
	 include 'qpincl.f'

        COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP

        IF( NAD.EQ.1 .OR. NAD.EQ.2 )THEN
                                   NLPTN(5,NAD)=2
                                   RETURN
                                    ENDIF

        NPR = NAD - NLPTN(2,NAD)
        L1 = NAD
        L2 = L1 + 1
        X= PLPTN(1,L1) / PLPTN(1,NPR)
        PLPTN(8,L1)= PLPTN(8,NPR)*X
        PLPTN(8,L2)= PLPTN(8,NPR)*(1.0D0-X)

        PT2=X*(1.0D0-X)
     &  *(PLPTN(2,NPR)-PLPTN(2,L1)/X-PLPTN(2,L2)/(1.0D0-X))
        IF( PT2 .LT. 0.0)THEN
              WRITE(6,*)' ERROR AT CMSINT, PT2<0, PT2=',PT2
              write(6,*)'l1,l2=',l1,l2,'1-x=',1.0d0-x
              CALL QPDUMP
	      stop
                         ENDIF

         THETA = 2.0*PI*QPRAND(1)
         PT = SQRT( PT2 )
         PX = PT*COS(THETA)
         PY = PT*SIN(THETA)

         PLPTN(3,L1) = PT2
         PLPTN(3,L2) = PT2
         PLPTN(4,L1) = PLPTN(4,NPR)*X         + PX
         PLPTN(4,L2) = PLPTN(4,NPR)*(1.0D0-X) - PX
         PLPTN(5,L1) = PLPTN(5,NPR)*X         + PY
         PLPTN(5,L2) = PLPTN(5,NPR)*(1.0D0-X) - PY
         PLPTN(9,L1)=(PLPTN(2,L1)
     &   +PLPTN(4,L1)**2+PLPTN(5,L1)**2)/PLPTN(8,L1)/2.0D0
         PLPTN(9,L2)=(PLPTN(2,L2)
     &   +PLPTN(4,L2)**2+PLPTN(5,L2)**2)/PLPTN(8,L2)/2.0D0

         PLPTN(6,L1)=( PLPTN(8,L1)-PLPTN(9,L1) )/SQRT(2.0D0)
         PLPTN(6,L2)=( PLPTN(8,L2)-PLPTN(9,L2) )/SQRT(2.0D0)
         IF( NLPTN(1,L1) .LE. 0) THEN
              PLPTN(6,L1) = - PLPTN(6,L1)
              PLPTN(6,L2) = - PLPTN(6,L2)
                                 ENDIF

         PLPTN(7,L1)=( PLPTN(8,L1)+PLPTN(9,L1) )/SQRT(2.0D0)
         PLPTN(7,L2)=( PLPTN(8,L2)+PLPTN(9,L2) )/SQRT(2.0D0)

         NLPTN(5,L1) = 2
         NLPTN(5,L2) = 2

          RETURN
          END

         SUBROUTINE QPDUMP

         IMPLICIT REAL*8(A-H,O-Z)
	 include 'qpincl.f'
        COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP
        DIMENSION PSUM(3)

        WRITE(6,101)NTOP
 101    FORMAT(' LIST OF OUTPUT, # OF LISTS=',I5)
        WRITE(6,102)
 102    FORMAT('     ID',6X,' X ',6X,'Q2',8X,'PT2',8X,'PX',8X,'PY',8X,
     &  'PZ',8X,'E',7X,'NLPTN(I,*)=',
     &  ' I=2',' I=3',' I=4',' I=5',' I=6')


        DO 110 I =1,NTOP
        WRITE(6,103)I,NLPTN(1,I),(PLPTN(K,I),K=1,7),(NLPTN(K,I),K=2,6)
C
 110    CONTINUE

 103    FORMAT(' ',I3,I3,2X,7(G9.2,1X),3X,5I3)
C
        DO 203 I=1,3
 203    PSUM(I)=0.0D0
C
        DO 201  J=1,NTOP
c      modificaton by HSS
c      NLPTN(3,J) - > NLPTN(1,J)
      IF( NLPTN(1,J).EQ.22 .AND. NLPTN(6,J) .GT. 0)THEN
        DO 202 I=1,3
 202    PSUM(I)=PSUM(I) + PLPTN(3+I,J)
                                                  ENDIF
 201   CONTINUE
C
       PSUM2 = 0.0D0
        DO 204 I=1,3
 204    PSUM2  =PSUM(I)**2 + PSUM2
C
        IF( PSUM2 .GT. 1.0D-3)THEN
        WRITE(6,*)' MOMENTUM CONSERVATION ???'
        WRITE(6,*)'  MOMENTUM SUM =',(PSUM(I),I=1,3)
  
                              ENDIF
C

        write(6,*)' end of QPDUMP '
	RETURN
        END



        SUBROUTINE QPANNH(Q2,NOK)
         IMPLICIT REAL*8(A-H,O-Z)
	 include 'qpincl.f'

        COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP

C         determine the virtual mass squared of the virtual photon
C         Then decide to accept or not its value according to
C         the probability (1/Q2)
C
            nok = 0
C

          P0X=0.0
          P0Y=0.0
          P0Z=0.0
          P0E=0.0

          LL = 0
         DO 1 I=1,NTOP
          IF( NLPTN(6,I).LT.0  .AND.  NLPTN(3,I).EQ.0)THEN
          P0X = P0X + PLPTN(4,I)
          P0Y = P0Y + PLPTN(5,I)
          P0Z = P0Z + PLPTN(6,I)
          P0E = P0E + PLPTN(7,I)
          LL = LL + 1
                                                     ENDIF
 1       CONTINUE
         IF( LL.NE.2)WRITE(6,*)'ERROR AT ANNH, LL=',LL
         IF( LL.NE.2)CALL QPDUMP
         IF( ABS(P0E) .GT. 1.0D10 )CALL QPDUMP
         IF(     P0E  .LE. 0.0D0  )RETURN

         Q2 = P0E**2 - P0X**2 - P0Y**2 - P0Z**2
cxxxxx
         IF(  Q2 .LT. Q2MIN) RETURN

c       IF( IBR .EQ. 1 )THEN
c                            NOK = 1
c                            GOTO 110
c                       ENDIF

cxxxxxxxxxxxx
              nok = 1
cxxxxxxxxxxxx

 110    CONTINUE

 120    CONTINUE
          L=NTOP + 1
         PLPTN(1,L)=1.0D0
         PLPTN(2,L)= Q2
         PLPTN(4,L)= P0X
         PLPTN(5,L)= P0Y
         PLPTN(6,L)= P0Z
         PLPTN(7,L)= P0E

         NLPTN(1,L) = 0
         NLPTN(2,L) = 0
         NLPTN(3,L) = 0
         NLPTN(4,L) = 0
         NLPTN(5,L) = 0
         NLPTN(6,L) = 1


         NTOP = NTOP +1


           RETURN
           END
C------------------------------------------------------------------
C         RANDOM
C------------------------------------------------------------------
C         RANDOM NUMBER  ( 0- 1.0 , UNIFORM )
C------------------------------------------------------------------
C
      FUNCTION QPRAND( IDUM  )

      IMPLICIT REAL*8(A-H,O-Z)
 
      COMMON /QPRND / NSEED
      NSEED = NSEED *48828125
      IF(NSEED.LT.0) NSEED=(NSEED+2147483647) + 1
      QPRAND=NSEED*0.4656613D-9
      RETURN
      END
c--------------------------------------------------
        SUBROUTINE QPcorr(xxx)
c see hep-ph/9603322
        IMPLICIT REAL*8(A-H,O-Z)
	include 'qpincl.f'
	data zbmin/1.0d-4/
cccccccccccccccccccccccccccc
c       compensate the finite contribution due to doulble cascade
c
          z = 0.0d0
         fm= QPfint(zbmin,1,0.1d0)
         if( QPRAND(1) .le. fm )then
             z2lmn=log(zbmin)**2
             z2lmx=log(1.0d0)**2
           dz2l = z2lmx - z2lmn
             fmx = 1.0d0/2.0d0/pi*alpha *2.5d0
         do 35 i=1,30
            z2l= z2lmn+QPRAND(1)*dz2l
           fx = QPfint(zbmin,2,z2l)
           if( fmx .le. fx)then
            write(6,*)' fmx < finta ',fmx,fx,' z2l= ',z2l
            z1 =exp(-sqrt(z2l))
            z = 1.0d0 - z1
            zs = sqrt(z)
           fintx =(log(z1/(1.0d0+zs))*(1.0d0+z*z)
     &   -2.0d0*z*log(zs/(1.0d0+zs))+z*z1*log(zs))/log(z1)

c
           ftx =2/sqrt(z)
            write(6,*)'z1,z,zs=',z1,z,zs,' fintx =',fintx,'ftx =',ftx
                           endif
           if( fmx*QPRAND(1) .le. fx)goto 36
 35       continue
          write(6,*)' fail to fix z '
 36       continue
          z = exp(-sqrt(z2l))
          z=1.0d0 -z
cx            write(6,*)' z,fx,fmx=',z,fx,fmx
                           else
           z =1.0d0
                           endif
          xxx = z
          return
          end

      FUNCTION  QPfint(zbmin,nnn,z2l)
      IMPLICIT REAL*8(A-H,O-Z)
      include 'qpincl.f'
ccccccccccccc
c      finite contribution
c
c          z1 = 1 - z
c          z1s = 1 - sqrt(z)
c        f(z) =C*delat(z1)-1/2/pi*alpha*log(z1s)/z1*(2-2z1+z1**2)
c                  -1/2/pi*alpha*(-2z/z1 log(sqrt(z)/(1+sqrt(z))
c                  -1/2/pi*alpha*(zlog(sqrt(z))
c
c        C is fixed by int (f(z)) = 1
c
c         variable z2l =log(z1)**2
ccccccccccc
       fac =1.0/2.0d0/pi*alpha
       if(nnn .eq. 1)then
cx           finta =fac*(log(zbmin))**2
           finta =1.0d0-exp(-fac*(log(zbmin))**2)
	   qpfint=finta
           if( finta .ge. 1.0d0)then
            write(6,*)'  int finta > 1 '
            write(6,*)' zbmin,fac,finta =',zbmin,fac,finta
            stop
                                endif
                     else
            z1 =exp(-sqrt(z2l))
            z = 1.0d0 - z1
            zs = sqrt(z)
           finta =fac*(log(z1/(1.0d0+zs))*(1.0d0+z*z)
     &   +2.0d0*log((1.0+zs)/zs)+z1**2*log(zs)
     &   - (zs-z)*z1)/log(z1)
c
c
        finta=finta*exp(-fac*(log(z1))**2)
	qpfint=finta

c
           if( finta .le. 0.0d0)then
cx            write(6,*)'  finta < 0 '
cx            write(6,*)' z2l,z1,,fac,finta =',z2l,z1,fac,finta
cx              stop
                                endif
                     endif

        return
        end
cccccc    end   ccccccccccccccccc
ccccccccccc    improved QEDPS    1995.8.8
ccccccccccccccccccccccccccccccccccccccccccccccc
        subroutine QPSET(ID,e,px,py,pz,ifirst)
        IMPLICIT REAL*8(A-H,O-Z)
        COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP
c       add by HSS
        SAVE NVPH
*
        if(ifirst.eq.0) then
	  NVPH = NTOP
          NLPTN(3,NVPH) = 2
          NLPTN(4,NVPH) = 1
          NLPTN(5,NVPH) = 2
        end if

          NTOP = NTOP + 1
          PLPTN(1,NTOP) = 1.0
          NLPTN(1,NTOP) = ID
          NLPTN(2,NTOP) = NTOP  - NVPH
          NLPTN(5,NTOP) = 2
          NLPTN(6,NTOP) = 1

          pt2=px**2+py**2
          ppp=sqrt(pt2+pz**2)
          PLPTN(3,NTOP)=0.d0
          PLPTN(4,NTOP)=px
          PLPTN(5,NTOP)=py
          PLPTN(6,NTOP)=pz
          PLPTN(7,NTOP)=e
          ! modification by HSS
          ! eng -> e
          PLPTN(8,NTOP)=( ppp+e)/sqrt(2.d0)
          PLPTN(9,NTOP)=(-ppp+e)/sqrt(2.d0)
*

         if(ifirst.eq.1) then
            ! after generating the final particle in their CMS
            ! boost to the lab frame (i.e. e- e+ CMS before showering)
          CALL QPBSFL(NVPH)
         end if

         RETURN
         END
C**********
C* BOOST$ *
C**********
C MODULE NO.=034
C.....................MAINTENANCE RECORD AREA........................
C   VERSION 84/06/08
C   UPDATE  00/00/00 : PERSON( XXXXXXX )
C
C....................................................................
C
c     SUBROUTINE BOOST$( PPPP , PREF , NNN1 , NNN2 )
      SUBROUTINE QPBOST( PPPP , PREF , NNN1 , NNN2 )
C
      IMPLICIT REAL*8(A-H,O-Z)
C
      DIMENSION PPPP(4,100), PREF(4)
CCC   DATA EPS/1.0E-5/
C
C     BOOST A SET OF 4-VECTOR PPPP(4,I) I=NNN1 TO NNN2
C     BY 4-VECTOR PREF(4)
C
C     BE CAREFUL FOR SIZE OF ARRAY PPPP
C
      SUM =0.
      DO 10 K=1,3
   10 SUM = SUM+PREF(K)**2
      P = SQRT(SUM)
      E0 = PREF(4)
      IF(E0.LE.0.0) THEN
          WRITE(6,601) PREF
  601     FORMAT(' #EPOCS3401(W,BOOST$) ENERGY<=0. PREF=',4G15.5)
          E0=1.E-30
          END IF
      BET = P/E0
C                           NO NEED OF BOOST
CCC   IF(BET.LE.EPS) RETURN
      IF(BET.GE.1.0) THEN
          WRITE(6,602) PREF
  602     FORMAT(' #EPOCS3402(W,BOOST$) VELOCITY>=1. PREF=',4G15.5)
          WRITE(6,612) BET
  612     FORMAT(' #EPOCS3402(W,BOOST$) VELOCITY>=1. BET =',G15.5)
CCC       CALL OUTLST
          BET=0.99999
          END IF
      BE2 =    SQRT((1.-BET)*(1.+BET))
      GAM = 1./BE2
      GA1 = 1./(BE2*(1.+BE2))
C
      DO 100 N = NNN1 , NNN2
        S = 0.
        DO 20 K = 1,3
        S = S + PREF(K)*PPPP(K,N)
   20   CONTINUE
        E = PPPP(4,N)
CCC     PX = (GAM-1.)*S/P**2+GAM*E/E0
        PX =  GA1    *S/E0**2+GAM*E/E0
        DO 30 K = 1,3
        PPPP(K,N) = PPPP(K,N) + PX*PREF(K)
   30   CONTINUE
        PPPP(4,N) = GAM*( E + S/E0 )
  100 CONTINUE
C
      RETURN
      END
C
C------------------------------------------------------------------
C     BOOST PARTICLES IN THE FINAL RADIATION
C------------------------------------------------------------------
C
C     SUBROUTINE BSFNL(NVPH0  ! ORIGINAL NAME
      SUBROUTINE QPBSFL(NVPH)
C
      IMPLICIT REAL*8(A-H,O-Z)
        COMMON/QPLIST/PLPTN(10,1000),NLPTN(10,1000),NTOP
*-----------------------------------------------------------------------
C
      DIMENSION PPPP(4,100), PREF(4)
C
      P2 = 0.0
      DO 11 J=1,3
 11   P2       =PLPTN(J+3,NVPH)**2 +P2
      IF( P2   .LE. 1.0E-5)RETURN
C
      DO 1 I=NVPH+1,NTOP
      DO 1 J=1,4
      PPPP(J,I)=PLPTN(J+3,I)
 1    CONTINUE
C
      DO 3 J=1,4
      PREF(J)=PLPTN(J+3,NVPH)
 3    CONTINUE
      NNN1 = NVPH +1
      NNN2 = NTOP
CCC   IF( PREF(4) .LE.  0.0)CALL OUTLST
C
c     CALL BOOST$( PPPP , PREF , NNN1 , NNN2 )
      CALL QPBOST( PPPP , PREF , NNN1 , NNN2 )
C
      DO 2 I=NVPH+1,NTOP
      DO 2 J=1,4
      PLPTN(J+3,I)  =  PPPP(J,I)
 2    CONTINUE

       RETURN
       END
