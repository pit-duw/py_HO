MODULE Summation_Pro
USE Helac_Global
USE Constants
USE SinglePro
IMPLICIT NONE
LOGICAL::lgluon
INTEGER::nproc,ngluonsi,npro
SAVE lgluon,ngluonsi,nproc
INTEGER,DIMENSION(10000)::mul
INTEGER,DIMENSION(10000,100,15)::ixnew

CONTAINS
SUBROUTINE SUMProcess
IMPLICIT NONE
REAL(KIND(1d0)),DIMENSION(11)::cp=(/35,3,-3,4,-4,8,-8,7,-7,12,-12/)  ! modification
INTEGER,DIMENSION(100)::ip,ix
CHARACTER(1)::last
CHARACTER(10),DIMENSION(20)::cha
CHARACTER(100)::process
INTEGER::iifl,ffl,nj,k,np,i,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,j
!REAL(KIND(1d0))::totmul
INTEGER::n1

nproc=0
mul(1:10000)=1
! Here,n1 is the number of final particles
WRITE(*,*)'-----------------------------------------------------------------------'
WRITE(*,*)' The number of final particles :'
READ *,n1
WRITE(*,*)'The number of final particles in the process is :',n1
n=n1+2
WRITE(*,*)' '
WRITE(*,*)'Flavour of final particles in SM : '
WRITE(*,*)'nue=1 e=2 u=3 d=4 numu=5 mu=6 c=7 s=8 nutau=9 tau=10 t=11 b=12'
WRITE(*,*)'A=31 Z=32 W+=33 W-=34 g=35 H=41 G0=42 G+=43 G-=44 Jet=100'
WRITE(*,*)'all antifermions are corresponding negative integers '
WRITE(*,*)'Jets are only allowed in the summation mode'
READ *,(ip(i),i=1,n1)
WRITE(*,*)(ip(i),i=1,n1)
WRITE(*,*)'------------------------------------------------------------------------'
WRITE(*,*)' '
! the light flavor in initial state and final states
CALL ReadElem_integer('iqnum',iifl)
CALL ReadElem_integer('fqnum',ffl)
WRITE(*,*)' The number of light flavor of initial states and final states are ',iifl,' and ',ffl
iifl=2*iifl+1
ffl=2*ffl+1
! the number of jets      
nj=0
DO k=1,n1
	IF(ip(k).EQ.100)THEN
		nj=nj+1
	ENDIF
ENDDO

np=n1-nj
! always we put the jets at the last  
DO i=1,np
   IF(ip(i).EQ. 2)cha(i)='e'
   IF(ip(i).EQ.-2)cha(i)='ebar'
   IF(ip(i).EQ. 1)cha(i)='nue'
   IF(ip(i).EQ.-1)cha(i)='nuebar'
   IF(ip(i).EQ. 3)cha(i)='u'
   IF(ip(i).EQ.-3)cha(i)='ubar'
   IF(ip(i).EQ. 4)cha(i)='d'
   IF(ip(i).EQ.-4)cha(i)='dbar'
   IF(ip(i).EQ. 6)cha(i)='mu'
   IF(ip(i).EQ.-6)cha(i)='mubar'
   IF(ip(i).EQ. 5)cha(i)='numu'
   IF(ip(i).EQ.-5)cha(i)='numubar'
   IF(ip(i).EQ. 8)cha(i)='s'
   IF(ip(i).EQ. 7)cha(i)='c'
   IF(ip(i).EQ.-7)cha(i)='cbar'
   IF(ip(i).EQ.-8)cha(i)='sbar'
   IF(ip(i).EQ. 10)cha(i)='tau'
   IF(ip(i).EQ.-10)cha(i)='taubar'
   IF(ip(i).EQ. 9)cha(i)='nutau'
   IF(ip(i).EQ.-9)cha(i)='nutaubar'
   IF(ip(i).EQ. 11)cha(i)='t'
   IF(ip(i).EQ.-11)cha(i)='tbar'
   IF(ip(i).EQ. 12)cha(i)='b'
   IF(ip(i).EQ.-12)cha(i)='bbar'
   IF(ip(i).EQ.31)cha(i)='A'
   IF(ip(i).EQ.32)cha(i)='Z'
   IF(ip(i).EQ.33)cha(i)='W+'
   IF(ip(i).EQ.34)cha(i)='W-'
   IF(ip(i).EQ.35)cha(i)='g'
   IF(ip(i).EQ.41)cha(i)='H'
ENDDO   
!DO i=1,np
!	IF(ip(i).EQ.1)cha(i)='ve_'
!    IF(ip(i).EQ.-1)cha(i)='ve~'
!    IF(ip(i).EQ.2)cha(i)='e-'
!	IF(ip(i).EQ.-2)cha(i)='e+'
!    IF(ip(i).EQ.5)cha(i)='vm_'
!    IF(ip(i).EQ.-5)cha(i)='vm~'
!	IF(ip(i).EQ.6)cha(i)='m-'
!	IF(ip(i).EQ.-6)cha(i)='m+'
!	IF(ip(i).EQ.9)cha(i)='vT'
!	IF(ip(i).EQ.-9)cha(i)='vT~'
!	IF(ip(i).EQ.10)cha(i)='T-'
!	IF(ip(i).EQ.-10)cha(i)='T+'
!    IF(ip(i).EQ.11)cha(i)='t_'
!    IF(ip(i).EQ.-11)cha(i)='t~'
!    IF(ip(i).EQ.12)cha(i)='b_'
!	IF(ip(i).EQ.-12)cha(i)='b~'
!	IF(ip(i).EQ.31)cha(i)='A_'
!	IF(ip(i).EQ.32)cha(i)='Z_'
!	IF(ip(i).EQ.33)cha(i)='W+'
!	IF(ip(i).EQ.34)cha(i)='W-'
!	IF(ip(i).EQ.41)cha(i)='H0'
!ENDDO
WRITE(*,*)' '
WRITE(*,*)' The number of final partons(not jets) is ',np
WRITE(process,*)(TRIM(cha(i)),i=1,np)
WRITE(last,'(I1)')nj
!PRINT *,TRIM(process(2:10*np+1))

OPEN(UNIT=10,FILE=TRIM(tmp_dir)//'Pinput'//TRIM(process(2:10*np+1))//'+'//last//'all')
OPEN(UNIT=11,FILE=TRIM(tmp_dir)//'Pinput'//TRIM(process(2:10*np+1))//'+'//last//'red')

IF(iifl.EQ.1)THEN
	ix(1)= -2 
	ix(2)=  2

	ix(3:2+np)=ip(1:np)

	IF(nj.EQ.0) CALL check_sum(n1+2,ix(1:n1+2))

	IF(nj.EQ.1) THEN
		DO l1=1,ffl
			ix(2+np+nj)=cp(l1)
			CALL check_sum(n1+2,ix(1:n1+2))
        ENDDO
    ENDIF

	IF(nj.EQ.2) THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				CALL check_sum(n1+2,ix(1:n1+2))
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.3) THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					CALL check_sum(n1+2,ix(1:n1+2))
				ENDDO
			ENDDO
		ENDDO
	ENDIF

    IF(nj.EQ.4)THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						CALL check_sum(n1+2,ix(1:n1+2))
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.5) THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						DO l5=l4,ffl
							ix(2+np+5)=cp(l5)
							CALL check_sum(n1+2,ix(1:n1+2))
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.6)THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						DO l5=l4,ffl
							ix(2+np+5)=cp(l5)
							DO l6=l5,ffl
								ix(2+np+6)=cp(l6)
								CALL check_sum(n1+2,ix(1:n1+2))
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.7) THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						DO l5=l4,ffl
							ix(2+np+5)=cp(l5)
							DO l6=l5,ffl
								ix(2+np+6)=cp(l6)
								DO l7=l6,ffl
									ix(2+np+7)=cp(l7)
									CALL check_sum(n1+2,ix(1:n1+2))
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.8)THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						DO l5=l4,ffl
							ix(2+np+5)=cp(l5)
							DO l6=l5,ffl
								ix(2+np+6)=cp(l6)
								DO l7=l6,ffl
									ix(2+np+7)=cp(l7)
									DO l8=l7,ffl
										ix(2+np+8)=cp(l8)
										CALL check_sum(n1+2,ix(1:n1+2))
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.9) THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						DO l5=l4,ffl
							ix(2+np+5)=cp(l5)
							DO l6=l5,ffl
								ix(2+np+6)=cp(l6)
								DO l7=l6,ffl
									ix(2+np+7)=cp(l7)
									DO l8=l7,ffl
										ix(2+np+8)=cp(l8)
										DO l9=l8,ffl
											ix(2+np+9)=cp(l9)
											CALL check_sum(n1+2,ix(1:n1+2))
										ENDDO
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

	IF(nj.EQ.10)THEN
		DO l1=1,ffl
			ix(2+np+1)=cp(l1)
			DO l2=l1,ffl
				ix(2+np+2)=cp(l2)
				DO l3=l2,ffl
					ix(2+np+3)=cp(l3)
					DO l4=l3,ffl
						ix(2+np+4)=cp(l4)
						DO l5=l4,ffl
							ix(2+np+5)=cp(l5)
							DO l6=l5,ffl
								ix(2+np+6)=cp(l6)
								DO l7=l6,ffl
									ix(2+np+7)=cp(l7)
									DO l8=l7,ffl
										ix(2+np+8)=cp(l8)
										DO l9=l8,ffl
											ix(2+np+9)=cp(l9)
											DO l10=l9,ffl
												ix(2+np+10)=cp(l10)
												CALL check_sum(n1+2,ix(1:n1+2))
											ENDDO
										ENDDO
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDDO
		ENDDO
	ENDIF

ELSE

	DO j=1,iifl
		ix(1)=cp(j)
		DO k=1,iifl
			ix(2)=cp(k)

			ix(3:2+np)=ip(1:np)

			IF(nj.EQ.0) CALL check_sum(n1+2,ix(1:n1+2))

			IF(nj.EQ.1) THEN
				DO l1=1,ffl
					ix(2+np+nj)=cp(l1)
					CALL check_sum(n1+2,ix(1:n1+2))
				ENDDO
			ENDIF

			IF(nj.EQ.2) THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						CALL check_sum(n1+2,ix(1:n1+2))
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.3)THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							CALL check_sum(n1+2,ix(1:n1+2))
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.4)THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								CALL check_sum(n1+2,ix(1:n1+2))
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.5) THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								DO l5=l4,ffl
									ix(2+np+5)=cp(l5)
									CALL check_sum(n1+2,ix(1:n1+2))
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.6)THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								DO l5=l4,ffl
									ix(2+np+5)=cp(l5)
									DO l6=l5,ffl
										ix(2+np+6)=cp(l6)
										CALL check_sum(n1+2,ix(1:n1+2))
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.7)THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								DO l5=l4,ffl
									ix(2+np+5)=cp(l5)
									DO l6=l5,ffl
										ix(2+np+6)=cp(l6)
										DO l7=l6,ffl
											ix(2+np+7)=cp(l7)
											CALL check_sum(n1+2,ix(1:n1+2))
										ENDDO
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.8)THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								DO l5=l4,ffl
									ix(2+np+5)=cp(l5)
									DO l6=l5,ffl
										ix(2+np+6)=cp(l6)
										DO l7=l6,ffl
											ix(2+np+7)=cp(l7)
											DO l8=l7,ffl
												ix(2+np+8)=cp(l8)
												CALL check_sum(n1+2,ix(1:n1+2))
											ENDDO
										ENDDO
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.9)THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								DO l5=l4,ffl
									ix(2+np+5)=cp(l5)
									DO l6=l5,ffl
										ix(2+np+6)=cp(l6)
										DO l7=l6,ffl
											ix(2+np+7)=cp(l7)
											DO l8=l7,ffl
												ix(2+np+8)=cp(l8)
												DO l9=l8,ffl
													ix(2+np+9)=cp(l9)
													CALL check_sum(n1+2,ix(1:n1+2))
												ENDDO
											ENDDO
										ENDDO
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

			IF(nj.EQ.10) THEN
				DO l1=1,ffl
					ix(2+np+1)=cp(l1)
					DO l2=l1,ffl
						ix(2+np+2)=cp(l2)
						DO l3=l2,ffl
							ix(2+np+3)=cp(l3)
							DO l4=l3,ffl
								ix(2+np+4)=cp(l4)
								DO l5=l4,ffl
									ix(2+np+5)=cp(l5)
									DO l6=l5,ffl
										ix(2+np+6)=cp(l6)
										DO l7=l6,ffl
											ix(2+np+7)=cp(l7)
											DO l8=l7,ffl
												ix(2+np+8)=cp(l8)
												DO l9=l8,ffl
													ix(2+np+9)=cp(l9)
													DO l10=l9,ffl
														ix(2+np+10)=cp(l10)
														CALL check_sum(n1+2,ix(1:n1+2))
													ENDDO
												ENDDO
											ENDDO
										ENDDO
									ENDDO
								ENDDO
							ENDDO
						ENDDO
					ENDDO
				ENDDO
			ENDIF

		ENDDO
	ENDDO

ENDIF

IF(onlyqcd)THEN
!	totmul=0
	DO j=1,npro
!		WRITE(10,*)ixnew(j,1,1:nj+np+2),mul(j)
		ifl(1:n)=ixnew(j,1,1:n)
		CALL SinProcess
!     do k=1,mul(j)
!      write(11,*)'  ',ixnew(j,k,1:nj+np+2)
!     enddo
!		totmul=totmul+mul(j)
	ENDDO
ENDIF


CLOSE(10)
CLOSE(11)
STOP
END SUBROUTINE SUMProcess

SUBROUTINE check_sum(n1,ix1)
IMPLICIT NONE
INTEGER,INTENT(IN)::n1
INTEGER,DIMENSION(n1),INTENT(IN)::ix1
INTEGER,DIMENSION(-12:41)::q
INTEGER::qtot
INTEGER,DIMENSION(-12:41)::f1,f2,f3,f4,f5,f6
INTEGER::f1tot,f2tot,f3tot,f4tot,f5tot,f6tot
INTEGER,DIMENSION(20)::m
INTEGER,DIMENSION(-12:41)::le,lm,lt,l1,l2,l3,lq
INTEGER::letot,lmtot,lttot,l1tot,l2tot,l3tot,lqtot,itot
REAL(KIND(1d0)),DIMENSION(-12:41)::i
REAL(KIND(1d0))::iini,ifin
LOGICAL::lgo1,lgo2,lgo3
INTEGER::kk1,jj1
!LOGICAL::onlyqcd
!      common/onlyqcd/onlyqcd
!      logical lgluon
!INTEGER::ngluonsi,lgluon,nproc
!      common/gluons/lgluon,ngluonsi
!      common/local/nproc
INTEGER::init=0,ngluons
INTEGER::nglu
SAVE

IF(init.EQ.0)THEN
!      print*,'onlyqcd : if T no flavour changing'
!      read*,onlyqcd
!      print*,'onlyqcd=  ',onlyqcd
	CALL ReadElem_integer('qcd',nglu)
	IF(nglu.EQ.2)THEN
		onlyqcd=.TRUE.
		PRINT *,' No Flavor Changing'
	ELSE
		onlyqcd=.FALSE.
	ENDIF
	CALL ReadElem_integer('nglu',nglu)
	IF(nglu.LT.0)THEN
		lgluon=.FALSE.
		ngluonsi=0
	ELSE
		lgluon=.TRUE.
		ngluonsi=nglu
		PRINT *,' The number of gluons is constrained to ',ngluonsi
	ENDIF
!	PRINT *,'Do you like to constrain number of gluons (ngluons=)'
!    READ(*,*)lgluon,ngluonsi
!	PRINT *,'lgluon=  ',lgluon,ngluonsi
! 3*Charges
	q(1)=0
	q(2)=-3
	q(3)=2
	q(4)=-1
	q(5)=q(1)
	q(6)=q(2)
	q(7)=q(3)
	q(8)=q(4)
	q(9)=q(1)
	q(10)=q(2)
	q(11)=q(3)
	q(12)=q(4)
	q(31)=0
	q(32)=0
	q(33)=3
	q(34)=-3
	q(35)=0
	q(41)=0
! 2*isospin
	i(1)=1
	i(2)=-1
	i(3)=1
	i(4)=-1
	i(5)=i(1)
	i(6)=i(2)
	i(7)=i(3)
	i(8)=i(4)
	i(9)=i(1)
	i(10)=i(2)
	i(11)=i(3)
	i(12)=i(4)
	i(31)=0
	i(32)=0
	i(33)=2
	i(34)=-2
	i(35)=0
	i(41)=0

! L-e   i.e. nu_e and e
	le(1:41)=0
	le(1:2)=1
! L-muon i.e. nu_mu and mu
	lm(1:41)=0
	lm(5:6)=1
! L-tau  i.e. nu_tau and tau
	lt(1:41)=0
	lt(9:10)=1

! L-ud  i.e. left-handed of u and d quark
	l1(1:41)=0
	l1(3:4)=1
! L-cs
	l2(1:41)=0
	l2(7:8)=1
! L-tb
	l3(1:41)=0
	l3(11:12)=1

! L-UD  i.e. left-handed of Up-type (u,c,t) and Down-type (d,s,b) quark
	lq(1:41)=0
	lq(3:4)=1
	lq(7:8)=1
	lq(11:12)=1

! f-u  i.e. u quark
	f1(1:41)=0
	f1(3)=1
! f-d
	f2(1:41)=0
	f2(4)=1
! f-c
	f3(1:41)=0
	f3(7)=1
! f-s
	f4(1:41)=0
	f4(8)=1
! f-t
	f5(1:41)=0
	f5(11)=1
! f-b
	f6(1:41)=0
	f6(12)=1

    DO kk1=1,12
       q(-kk1)=-q(kk1)
       i(-kk1)=-i(kk1)
       le(-kk1)=-le(kk1)
       lm(-kk1)=-lm(kk1)
       lt(-kk1)=-lt(kk1)
       l1(-kk1)=-l1(kk1)
       l2(-kk1)=-l2(kk1)
       l3(-kk1)=-l3(kk1)
       lq(-kk1)=-lq(kk1)
       f1(-kk1)=-f1(kk1)
       f2(-kk1)=-f2(kk1)
       f3(-kk1)=-f3(kk1)
       f4(-kk1)=-f4(kk1)
       f5(-kk1)=-f5(kk1)
       f6(-kk1)=-f6(kk1)
	ENDDO

	m(3:20)=1
	m(1:2)=-1

	init=1
ENDIF

! the 3*total Charges
qtot=0
DO kk1=1,n1
	qtot=qtot+q(ix1(kk1))*m(kk1)
ENDDO
! the 2* total isospin
itot=0
DO kk1=1,n1
	itot=itot+i(ix1(kk1))*m(kk1)
ENDDO
! the total electron lepton number
letot=0
DO kk1=1,n1
	letot=letot+le(ix1(kk1))*m(kk1)
ENDDO
! the total muon lepton number
lmtot=0
DO kk1=1,n1
	lmtot=lmtot+lm(ix1(kk1))*m(kk1)
ENDDO
! the total tau lepton number
lttot=0
DO kk1=1,n1
	lttot=lttot+lt(ix1(kk1))*m(kk1)
ENDDO
! the total ud flavor number
l1tot=0
DO kk1=1,n1
	l1tot=l1tot+l1(ix1(kk1))*m(kk1)
ENDDO
! the total cs flavor number
l2tot=0
DO kk1=1,n1
	l2tot=l2tot+l2(ix1(kk1))*m(kk1)
ENDDO
! the total tb flavor number
l3tot=0
DO kk1=1,n1
	l3tot=l3tot+l3(ix1(kk1))*m(kk1)
ENDDO
! the total quark number
lqtot=0
DO kk1=1,n1
	lqtot=lqtot+lq(ix1(kk1))*m(kk1)
ENDDO
! the total u quark flavor number
f1tot=0
DO kk1=1,n1
	f1tot=f1tot+f1(ix1(kk1))*m(kk1)
ENDDO
! the total d quark flavor number
f2tot=0
DO kk1=1,n1
	f2tot=f2tot+f2(ix1(kk1))*m(kk1)
ENDDO
! the total c quark flavor number
f3tot=0
DO kk1=1,n1
	f3tot=f3tot+f3(ix1(kk1))*m(kk1)
ENDDO
! the total s quark flavor number
f4tot=0
DO kk1=1,n1
	f4tot=f4tot+f4(ix1(kk1))*m(kk1)
ENDDO
! the total t quark flavor number
f5tot=0
DO kk1=1,n1
	f5tot=f5tot+f5(ix1(kk1))*m(kk1)
ENDDO
! the total b quark flavor number
f6tot=0
DO kk1=1,n1
	f6tot=f6tot+f6(ix1(kk1))*m(kk1)
ENDDO

ngluons=0
lgo1=(qtot.EQ.0)
lgo2=(itot.EQ.0)
lgo3=(letot.EQ.0.AND.lmtot.EQ.0.AND.lttot.EQ.0)
! We choose the CKM approximate diagonal
lgo3=lgo3.AND.(l1tot.EQ.0.AND.l2tot.EQ.0.AND.l3tot.EQ.0)
!     lgo3=lgo3.and.(lqtot.eq.0)
IF(onlyqcd)lgo3=(f1tot.EQ.0.AND.f2tot.EQ.0.AND.f3tot.EQ.0.AND.&
              f4tot.EQ.0.AND.f5tot.EQ.0.AND.f6tot.EQ.0)
IF(lgo1.AND.lgo2.AND.lgo3)THEN
	DO jj1=1,n1
		IF(ix1(jj1).EQ.35)ngluons=ngluons+1
    ENDDO
    IF(lgluon)THEN
		IF(ngluons.EQ.ngluonsi)THEN
			nproc=nproc+1
			IF(onlyqcd)THEN
				CALL reduce_sum(n1,ix1)
			ELSE
!				WRITE(10,'(20I4)')ix1(1:n1)
				ifl(1:n1)=ix1(1:n1)
				CALL SinProcess
			ENDIF
		ENDIF
	ELSE
		nproc=nproc+1
		IF(onlyqcd)THEN
			CALL reduce_sum(n1,ix1)
		ELSE
!			WRITE(10,'(20I4)')ix1(1:n1)
			ifl(1:n1)=ix1(1:n1)
			CALL SinProcess
		ENDIF
	ENDIF
ENDIF
    
RETURN
END SUBROUTINE check_sum

SUBROUTINE reduce_sum(n1,ix1)
IMPLICIT NONE
INTEGER,INTENT(IN)::n1
INTEGER,DIMENSION(n1),INTENT(IN)::ix1
INTEGER,DIMENSION(n1)::ixn
INTEGER,DIMENSION(5)::f,gf,c=(/3,4,7,8,12/)
INTEGER,DIMENSION(10000,5)::g
INTEGER,DIMENSION(10000,2)::tag
!LOGICAL,DIMENSION(-12:41,-12:41)::log_gg,log_gq,log_qa,log_qg,log_ag,log_aq,log_qq,log_aa,&
!                                  log_qr,log_qb,log_ar,log_ab
!INTEGER::npro
!      common/mul/mul(10000),npro
INTEGER::ipro=0,i,ii,j,ishow
!      common/new/ixnew(10000,100,15)
SAVE c,ipro

!log_gg(i1,i2)=(i1.EQ.35.AND.i2.EQ.35)
!log_gq(i1,i2)=(i1.EQ.35.AND.(i2.GT.0.AND.i2.LE.12))
!log_ga(i1,i2)=(i1.EQ.35.AND.(i2.LT.0.AND.i2.GE.-12))
!log_qg(i1,i2)=(i2.EQ.35.AND.(i1.GT.0.AND.i1.LE.12))
!log_ag(i1,i2)=(i2.EQ.35.AND.(i1.LT.0.AND.i1.GE.-12))
!log_qa(i1,i2)=((i1.GT.0.AND.i1.LE.12).AND.(i2.LT.0.AND.i2.GE.-12) &
!              .AND.ABS(i1).EQ.ABS(i2))
!log_aq(i1,i2)=((i2.GT.0.AND.i2.LE.12).AND.(i1.LT.0.AND.i1.GE.-12) &
!              .AND.ABS(i1).EQ.ABS(i2))
!log_qq(i1,i2)=((i1.GT.0.AND.i1.LE.12).AND.(i2.GT.0.AND.i2.LE.12) &
!              .AND.ABS(i1).EQ.ABS(i2))
!log_aa(i1,i2)=((i1.LT.0.AND.i1.GE.-12).AND.(i2.LT.0.AND.i2.GE.-12) &
!              .AND.ABS(i1).EQ.ABS(i2))
!log_qr(i1,i2)=((i1.GT.0.AND.i1.LE.12).AND.(i2.GT.0.AND.i2.LE.12) &
!              .AND.ABS(i1).NE.ABS(i2))
!log_qb(i1,i2)=((i1.GT.0.AND.i1.LE.12).AND.(i2.LT.0.AND.i2.GE.-12) &
!              .AND.ABS(i1).NE.ABS(i2))
!log_ar(i1,i2)=((i2.GT.0.AND.i2.LE.12).AND.(i1.LT.0.AND.i1.GE.-12) &
!              .AND.ABS(i1).NE.ABS(i2))
!log_ab(i1,i2)=((i1.LT.0.AND.i1.GE.-12).AND.(i2.LT.0.AND.i2.GE.-12) &
!              .AND.ABS(i1).NE.ABS(i2))

ipro=ipro+1

f(1:5)=0
g(ipro,1:5)=0
gf(1:5)=0

DO i=1,n1 
	ii=ix1(i)
	IF(i.LE.2)ii=-ii
	IF(ii.EQ. 3 )f(1)=f(1)+1
	IF(ii.EQ. 4 )f(2)=f(2)+1
	IF(ii.EQ. 7 )f(3)=f(3)+1
	IF(ii.EQ. 8 )f(4)=f(4)+1
	IF(ii.EQ. 12)f(5)=f(5)+1
ENDDO
! sort the number of flavors as the descend array
DO j=1,5
	g(ipro,j)=MAX(f(1),f(2),f(3),f(4),f(5))
	DO i=1,5
		IF(g(ipro,j).EQ.f(i).AND.f(i).GT.0)THEN
			f(i)=0
			gf(j)=c(i)
			EXIT
		ENDIF
    ENDDO
ENDDO

ixn(1:n1)=ix1(1:n1)
DO i=1,n1
	DO j=1,5
		IF(ix1(i).EQ. gf(j))ixn(i)=  99+j
		IF(ix1(i).EQ.-gf(j))ixn(i)=-(99+j)
	ENDDO
ENDDO

IF(log_gg(ix1(1),ix1(2))) tag(ipro,1)=1
IF(log_gq(ix1(1),ix1(2))) tag(ipro,1)=2
IF(log_ga(ix1(1),ix1(2))) tag(ipro,1)=3
IF(log_qg(ix1(1),ix1(2))) tag(ipro,1)=4
IF(log_ag(ix1(1),ix1(2))) tag(ipro,1)=5
IF(log_qa(ix1(1),ix1(2))) tag(ipro,1)=6
IF(log_aq(ix1(1),ix1(2))) tag(ipro,1)=7
IF(log_qq(ix1(1),ix1(2))) tag(ipro,1)=8
IF(log_aa(ix1(1),ix1(2))) tag(ipro,1)=9
IF(log_qr(ix1(1),ix1(2))) tag(ipro,1)=10
IF(log_qb(ix1(1),ix1(2))) tag(ipro,1)=11
IF(log_ar(ix1(1),ix1(2))) tag(ipro,1)=12
IF(log_ab(ix1(1),ix1(2))) tag(ipro,1)=13

ishow=1
ixnew(ipro,1,1:n1)=ix1(1:n1)
! Exclude the same initial configuration (e.g. uubar(udbar) and ddbar(cbbar) are the same)
! and the number of the different flavor of quarks are corresponding the same
DO j=1,ipro-1
	IF(tag(ipro,1).EQ.tag(j,1))THEN
		IF(g(ipro,1).EQ.g(j,1).AND.g(ipro,2).EQ.g(j,2).AND. &
          g(ipro,3).EQ.g(j,3).AND.g(ipro,4).EQ.g(j,4).AND. &
          g(ipro,5).EQ.g(j,5))THEN
			mul(j)=mul(j)+1
			ixnew(j,mul(j),1:n1)=ix1(1:n1)
			ipro=ipro-1
			ishow=0
			EXIT
		ENDIF
	ENDIF
ENDDO

IF(ishow.EQ.1)WRITE(11,*)ipro,ix1(1:n1)
npro=ipro

RETURN
END SUBROUTINE reduce_sum

LOGICAL FUNCTION log_gg(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_gg=(i1.EQ.35.AND.i2.EQ.35)
END FUNCTION log_gg

LOGICAL FUNCTION log_gq(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_gq=(i1.EQ.35.AND.(i2.GT.0.AND.i2.LE.12))
END FUNCTION log_gq

LOGICAL FUNCTION log_ga(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_ga=(i1.EQ.35.AND.(i2.LT.0.AND.i2.GE.-12))
END FUNCTION log_ga

LOGICAL FUNCTION log_qg(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_qg=(i2.EQ.35.AND.(i1.GT.0.AND.i1.LE.12))
END FUNCTION log_qg

LOGICAL FUNCTION log_ag(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_ag=(i2.EQ.35.AND.(i1.LT.0.AND.i1.GE.-12))
END FUNCTION log_ag

LOGICAL FUNCTION log_qa(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_qa=((i1.GT.0.AND.i1.LE.12).AND.(i2.LT.0.AND.i2.GE.-12) &
              .AND.ABS(i1).EQ.ABS(i2))
END FUNCTION log_qa

LOGICAL FUNCTION log_aq(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_aq=((i2.GT.0.AND.i2.LE.12).AND.(i1.LT.0.AND.i1.GE.-12) &
              .AND.ABS(i1).EQ.ABS(i2))
END FUNCTION log_aq

LOGICAL FUNCTION log_qq(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_qq=((i1.GT.0.AND.i1.LE.12).AND.(i2.GT.0.AND.i2.LE.12) &
              .AND.ABS(i1).EQ.ABS(i2))
END FUNCTION log_qq

LOGICAL FUNCTION log_aa(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_aa=((i1.LT.0.AND.i1.GE.-12).AND.(i2.LT.0.AND.i2.GE.-12) &
              .AND.ABS(i1).EQ.ABS(i2))
END FUNCTION log_aa

LOGICAL FUNCTION log_qr(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_qr=((i1.GT.0.AND.i1.LE.12).AND.(i2.GT.0.AND.i2.LE.12) &
              .AND.ABS(i1).NE.ABS(i2))
END FUNCTION log_qr

LOGICAL FUNCTION log_qb(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_qb=((i1.GT.0.AND.i1.LE.12).AND.(i2.LT.0.AND.i2.GE.-12) &
              .AND.ABS(i1).NE.ABS(i2))
END FUNCTION log_qb

LOGICAL FUNCTION log_ar(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_ar=((i2.GT.0.AND.i2.LE.12).AND.(i1.LT.0.AND.i1.GE.-12) &
              .AND.ABS(i1).NE.ABS(i2))
END FUNCTION log_ar

LOGICAL FUNCTION log_ab(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
log_ab=((i1.LT.0.AND.i1.GE.-12).AND.(i2.LT.0.AND.i2.GE.-12) &
              .AND.ABS(i1).NE.ABS(i2))
END FUNCTION log_ab

SUBROUTINE Prodinput
IMPLICIT NONE
CHARACTER(50)::process
CHARACTER(30)::file
CHARACTER(1)::a1
CHARACTER(2)::a2
CHARACTER(3)::a3
CHARACTER(4)::a4
CHARACTER(70)::str
INTEGER::np,i
PRINT *,'file'
READ *,file
OPEN(UNIT=9,FILE= TRIM(tmp_dir)//file)
OPEN(UNIT=8,FILE=TRIM(input_dir)//'infile')
! the total number of external partilces including jets
READ(8,*)np
np=np+2
CLOSE(8)

i=0
! 1    continue
DO
! single processes
	i=i+1
	READ(9,'(A50)') process   ! When one read the endoffile it gets out of the loop
	PRINT *,i,process
	IF(i.LT.10) THEN
		WRITE(a1,'(I1)')i
		OPEN(UNIT=11,FILE=TRIM(tmp_dir)//'cinput'//a1)
	ELSEIF(i.GE.10.AND.i.LE.99)THEN
		WRITE(a2,'(I2)')i
		OPEN(UNIT=11,FILE=TRIM(tmp_dir)//'cinput'//a2)
	ELSEIF(i.GE.100.AND.i.LE.999)THEN
		WRITE(a3,'(I3)')i
		OPEN(UNIT=11,FILE=TRIM(tmp_dir)//'cinput'//a3)
	ELSEIF(i.GE.1000.AND.i.LE.9999)THEN
		WRITE(a4,'(I4)')i
		OPEN(UNIT=11,FILE=TRIM(tmp_dir)//'cinput'//a4)
    ELSEIF(i.GE.10000)THEN
		PRINT *,'number of subprocesses is >= 9999',i
		PRINT *,'please correct in SUBROUTINE Prodinput'
		STOP
    ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! the input file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!hi_file err_file               ! histo and error files
!0 1                            ! repeat ranhel
!4                              ! number of particl
!-8 7 33 35                       ! flavour of particles
!0 1 0 0                        ! iflag,iunitary,ihiggs,iwidth
!F T                            ! onlyqcd,withqcd
!F 1000                         ! unweighting (T/F) - number of pre-unweighted evts
!0                              ! ncha: option for MC generator 1 RAMBO 0 PHEGAS 2 DURHAM
!2000000000                     ! nmc: number of MC iterations
!10000,5000,1,6,70000,1         ! optimization: nopt,nopt_step,optf,maxopt,noptlim,iopt
!1960.                          ! collision energy
!1                              ! structure fucntions (1)
!0                              ! limit for optimization
!2212 2212 7000d0 7000d0 3 1    !IDBMUP1,IDBMUP2,EEBMUP1,EBMUP2,IDWTUP,NPRUP


	OPEN(UNIT=10000,FILE=TRIM(input_dir)//'input')
	READ(10000,'(A30)')str
	WRITE(11,*)str
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)np  ! number of partilces
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)process
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str   ! onlyqcd withqcd
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)'F 1000 1000'  ! unweighting,number of preweighting events
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str             ! ncha: option for MC generator 1 RAMBO 0 PHEGAS 2 DURHAM       
!---------------------------------------
    READ(10000,'(A30)')str
	WRITE(11,*)str            ! nmc
!---------------------------------------
	READ(10000,'(A50)')str
	WRITE(11,*)str            ! optimization: nopt,nopt_step,optf,maxopt,noptlim,iopt
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str            ! collision energy
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str            ! structure function
!---------------------------------------
	READ(10000,'(A30)')str
	WRITE(11,*)str            ! limit for optimization (alimit)
!---------------------------------------
	READ(10000,'(A70)')str
	WRITE(11,*)str            !IDBMUP1,IDBMUP2,EEBMUP1,EBMUP2,IDWTUP,NPRUP
!---------------------------------------
    CLOSE(10000)
	CLOSE(11)
ENDDO 
CLOSE(0)
CLOSE(9)

STOP
END SUBROUTINE Prodinput
END MODULE Summation_Pro
