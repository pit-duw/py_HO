MODULE Helac_master
USE Helac_Global
USE Helac_Func_1
USE Helac_pan1
USE Helac_Feynman
USE Projectors
USE Decay_interface
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,PUBLIC::rmatrix
COMPLEX(KIND=DBL),DIMENSION(:),ALLOCATABLE,PUBLIC::zamp
INTEGER,DIMENSION(:,:,:),ALLOCATABLE,PUBLIC::jcol
COMPLEX(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,PUBLIC::smhelcol
COMPLEX(KIND=DBL),DIMENSION(:,:),ALLOCATABLE,PUBLIC::smhelcol2
INTEGER::init=0,num,nq !,init_col=0 !,iv1,iv2
SAVE rmatrix,zamp,jcol
SAVE init,num !,init_col
!SAVE iv1,iv2
CONTAINS       
! In the begining ---> call SM_FeynRule_Helac first !     
SUBROUTINE Helac_init()
IMPLICIT NONE     
INTEGER,DIMENSION(20)::iquark
INTEGER,DIMENSION(:,:),ALLOCATABLE::iperm
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::divide  ! the factor for swap
INTEGER::ic1,init_unw,istat,l,ix1,ix2,ix3,ijij2,&
         i0,i1,i2,i3,nogo,if1,ino,izeros,ii,jj,&
		 k,kk,j,i,k1,k2,icc,icheck,ijxi,swap,ijxj,colorpos,kkkij,kjkj,ictemp1,ictemp2,ktemp,line
INTEGER,DIMENSION(10)::singletlist
INTEGER,DIMENSION(20,2)::slist
LOGICAL,DIMENSION(10)::used
INTEGER,DIMENSION(10)::ising
!REAL(KIND=DBL)::avhel,avcol,symet
ic1=0
init_unw=0
iquark=0
!*********************************************************
!     ipol(k)=helicity,ipol(k)=2:rigth-handed            *
!                      ipol(k)=1:left-handed             *
!                      ipol(k)=3:long.                   *
!*********************************************************
!*********************************************************
!      start of initialization phase      
!*********************************************************
WRITE(nunit1,*)' '
WRITE(nunit1,*)'--------------------------------------------------------------------------'
CALL Helac_average(avhel,avcol,symet)      ! return the helicity,color and symmetry factors
GLOBALINIT_average=1
WRITE(nunit1,*)'The average over helicity is avhel= ',avhel
WRITE(nunit1,*)'The average over color is avcol= ',avcol
WRITE(nunit1,*)'The symmetry factor symet= ',symet
WRITE(nunit1,*)' '
! WRITE(nunit1,*) avhel,avcol,symet
iv1=31
iv2=34
CALL Helac_setncc(nq)                     ! nq is the number of incoming quarks 
                                          ! or outgoing antiquarks or gluons
       
ncc=Helac_ifactorial(nq)                  ! complete array of nq colors

IF(ALLOCATED(iperm))THEN
	DEALLOCATE(iperm)
ENDIF
IF(ALLOCATED(divide))THEN
	DEALLOCATE(divide)
ENDIF
IF(ALLOCATED(jcol))THEN
	DEALLOCATE(jcol)
ENDIF    
ALLOCATE(iperm(ncc,nq),STAT=istat)
ALLOCATE(divide(ncc),STAT=istat)
ALLOCATE(jcol(ncc,n,0:2),STAT=istat)
IF(istat.NE.0)WRITE(*,*)istat, ' warning: allocation is not working properly in Helac_master !'
WRITE(nunit1,*)'Number of colour configurations ',ncc
WRITE(nunit1,*)' '
CALL Helac_pan(ncc)                       ! initialization in Helac_pan1
IF(ncc.GT.1)iv2=35
DO icc=1,ncc
   IF(nq.NE.0.AND.icc.LE.ncc)THEN
! the arrange of color configurations
      IF(nq.GT.1)CALL Helac_permu(iquark,nq)
      IF(nq.EQ.1)iquark(1)=1
      CALL Helac_setcc(iquark)
      CALL Helac_checkgluons(iquark,nq,icheck)
	  GLOBALINIT_checkgluons=1
      IF(icheck.EQ.0)CYCLE
   ENDIF
   WRITE(nunit1,*)'=============================================================='
   WRITE(nunit1,*)'THE COLOUR OF PARTICLES ONE :'
   WRITE(nunit1,*)'incoming quark,outgoing antiquark or gluon'
   WRITE(nunit1,*)(icol(i,1),i=1,n)
   WRITE(nunit1,*)' '
   WRITE(nunit1,*)'THE COLOUR OF PARTICLES TWO :'
   WRITE(nunit1,*)'incoming antiquark,outgoing quark or gluon'
   WRITE(nunit1,*)(icol(i,2),i=1,n)
   WRITE(nunit1,*)' '
   CALL Helac_pan_ini()             ! initialize the colors,ih and ife of external particles
   ! l means the level l
   DO l=2,n-1
      ix1=0
      DO WHILE(ix1.EQ.0)
       CALL Helac_id(n,l,ix1,i0)
	   GLOBALINIT_id=1
       IF(ix1.EQ.1) EXIT             ! goto 800
	   ! to exclude the first particle
       IF(MOD(i0,2).EQ.0)THEN
          CALL Helac_checki0(i0,l,nogo)
          IF(nogo.EQ.1)CYCLE
          ix2=0
		  ! returns all possible 3-points vertices
          DO WHILE(ix2.EQ.0)
              CALL Helac_id2(n,l,i0,ix2,i1,i2)
			  GLOBALINIT_id2=1
              IF(ix2.EQ.1) EXIT       ! goto801
              CALL Helac_checki0(i1,l,nogo)
              IF(nogo.EQ.1)CYCLE
              CALL Helac_checki0(i2,l,nogo)
              IF(nogo.EQ.1)CYCLE
              i3=0
              DO if1=iv1,iv2
                 CALL Helac_v3(i0,l,if1,i1,i2,i3)
                 CALL Helac_vff(i0,l,if1,i1,i2,i3)
                 CALL Helac_vvs(i0,l,if1,i1,i2,i3)
                 CALL Helac_vss(i0,l,if1,i1,i2,i3)
             ENDDO
             DO if1=41,44
                 CALL Helac_svv(i0,l,if1,i1,i2,i3)
                 CALL Helac_ssv(i0,l,if1,i1,i2,i3)
				 CALL Helac_sff(i0,l,if1,i1,i2,i3)
                 CALL Helac_s3(i0,l,if1,i1,i2,i3)
             ENDDO
             DO if1=1,12
                 CALL Helac_ffv(i0,l, if1,i1,i2,i3)
                 CALL Helac_ffs(i0,l, if1,i1,i2,i3)
                 CALL Helac_affv(i0,l,-if1,i1,i2,i3)
                 CALL Helac_affs(i0,l,-if1,i1,i2,i3)
             ENDDO
          ENDDO                        !801
          IF(l.LT.3)CYCLE             ! goto800.
          ix3=0
          DO WHILE(ix3.EQ.0)
             CALL Helac_id3(n,l,i0,ix3,i1,i2,i3)
			 GLOBALINIT_id3=1
             IF(ix3.EQ.1)EXIT          !goto802
             CALL Helac_checki0(i1,l,nogo)
             IF(nogo.EQ.1)CYCLE
             CALL Helac_checki0(i2,l,nogo)
             IF(nogo.EQ.1)CYCLE
             CALL Helac_checki0(i3,l,nogo)
             IF(nogo.EQ.1)CYCLE
             DO if1=iv1,iv2
                CALL Helac_v4(i0,l,if1,i1,i2,i3)
                CALL Helac_vvss(i0,l,if1,i1,i2,i3)
             ENDDO
             DO if1=41,44
                CALL Helac_s4(i0,l,if1,i1,i2,i3)
                CALL Helac_ssvv(i0,l,if1,i1,i2,i3)
             ENDDO
          ENDDO !802
       ENDIF
     ENDDO      !800
   ENDDO
   CALL Helac_redo(ic1,ino)          ! in Helac_pan1
   GLOBALINIT_redo=1
   IF(nq.NE.0.AND.icc.LE.ncc.AND.ino.EQ.0)THEN
        iperm(ic1,1:nq)=iquark(1:nq)
		singletlist(1:10)=0
		kjkj=0
		DO i=1,nhad
			IF(ABS(iflh(i)).LE.100)THEN
				jcol(ic1,i,1:2)=icol(hadron2parton(i),1:2)
			ELSE
				IF(SubInteger(iflh(i),1,1).NE.1)THEN
					kkkij=hadron2parton(i)
					jcol(ic1,i,1)=icol(kkkij,1)+icol(kkkij+1,1)
					jcol(ic1,i,2)=icol(kkkij,2)+icol(kkkij+1,2)
				ELSE
					jcol(ic1,i,1:2)=0
					kjkj=kjkj+1
					singletlist(kjkj)=i
				ENDIF
			ENDIF
		ENDDO
		slist(1:n,1:2)=icol(1:n,1:2)
		used(1:10)=.FALSE.
		DO i=1,kjkj
			kkkij=singletlist(i)
			IF(used(kkkij))CYCLE
			kkkij=hadron2parton(kkkij)
			ictemp1=slist(kkkij,1)+slist(kkkij+1,1)
			ictemp2=slist(kkkij,2)+slist(kkkij+1,2)
			DO ktemp=1,nhad
				IF(ktemp.NE.singletlist(i).AND..NOT.used(ktemp))THEN
					IF(SubInteger(iflh(ktemp),1,1).NE.1)THEN
						IF(ictemp1.EQ.slist(hadron2parton(ktemp),2))THEN
							jcol(ic1,ktemp,2)=MIN(ictemp1,ictemp2)
							slist(hadron2parton(ktemp),2)=MIN(ictemp1,ictemp2)
						ENDIF
						IF(ictemp2.EQ.slist(hadron2parton(ktemp),1))THEN
							jcol(ic1,ktemp,1)=MIN(ictemp1,ictemp2)
							slist(hadron2parton(ktemp),1)=MIN(ictemp1,ictemp2)
						ENDIF
                                                IF(ABS(iflh(ktemp)).GT.100)THEN
                                                   IF(ictemp1.EQ.slist(hadron2parton(ktemp)+1,2))THEN
                                                      slist(hadron2parton(ktemp)+1,2)=MIN(ictemp1,ictemp2)
                                                   ENDIF
                                                   IF(ictemp2.EQ.slist(hadron2parton(ktemp)+1,1))THEN
                                                      slist(hadron2parton(ktemp)+1,1)=MIN(ictemp1,ictemp2)
                                                   ENDIF
                                                ENDIF
					ELSE
                                                IF(ictemp1.EQ.slist(hadron2parton(ktemp),2))THEN
                                                      slist(hadron2parton(ktemp),2)=MIN(ictemp1,ictemp2)
                                                ENDIF
                                                IF(ictemp2.EQ.slist(hadron2parton(ktemp),1))THEN
                                                      slist(hadron2parton(ktemp),1)=MIN(ictemp1,ictemp2)
                                                ENDIF
                                                IF(ABS(iflh(ktemp)).GT.100)THEN
                                                   IF(ictemp1.EQ.slist(hadron2parton(ktemp)+1,2))THEN
                                                      slist(hadron2parton(ktemp)+1,2)=MIN(ictemp1,ictemp2)
                                                   ENDIF
                                                   IF(ictemp2.EQ.slist(hadron2parton(ktemp)+1,1))THEN
                                                      slist(hadron2parton(ktemp)+1,1)=MIN(ictemp1,ictemp2)
                                                   ENDIF
                                               ENDIF
					ENDIF
				ENDIF
			ENDDO
			used(singletlist(i))=.TRUE.
		ENDDO
        jcol(ic1,1,0)=1
        DO i=1,n
		! exclude the gluon in U(1)
           IF(ifl(i).EQ.35.AND.icol(i,1).EQ.icol(i,2))THEN
              jcol(ic1,1,0)=0
              EXIT             ! goto100
           ENDIF
        ENDDO
        IF(jcol(ic1,1,0).NE.0)THEN
           DO i=1,nhad
              IF((iflh(i).EQ.35.OR.&
                   (ABS(iflh(i)).GT.100.AND.SubInteger(iflh(i),1,1).EQ.8.AND.octetQ)).AND.&
                   jcol(ic1,i,1).EQ.jcol(ic1,i,2))THEN
                 jcol(ic1,1,0)=0
                 EXIT
              ENDIF
           ENDDO
        ENDIF
   !100    continue
   ENDIF
!   WRITE(*,*)icc ! Debug
ENDDO                     !icc
WRITE(nunit1,*)' ' 
WRITE(nunit1,*)'===================================================================='     
WRITE(nunit1,*)'The number of alived colour confs. ic1=',ic1
ncc=ic1                   !-1       
CALL Helac_setlist(ncc)
IF(ALLOCATED(zamp))THEN
	DEALLOCATE(zamp)
ENDIF
IF(ALLOCATED(rmatrix))THEN
	DEALLOCATE(rmatrix)
ENDIF
ALLOCATE(zamp(ncc),rmatrix(ncc,ncc),STAT=istat)
IF(istat.NE.0)WRITE(*,*)istat,&
      'warning: allocation is not working properly in Helac_init of Helac_master'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! reordering for the color siglet of quarkonia
IF(octetQ)THEN
! real color octet
	rmatrix(1:ncc,1:ncc)=0
	octetmsinglet=0
	CALL Helac_rmatrix(iperm(1:ncc,1:nq))
	DO k1=1,octetnum
		GLOBALINIT_idoctet=0
		ijxi=0
		DO WHILE(ijxi.EQ.0)
			CALL Helac_idoctet(octetnum,k1,ijxi,k2)
			GLOBALINIT_idoctet=1
			PRINT *,ijxi,k2
			IF(ijxi.NE.0)EXIT
			CALL Helac_bin(octetnum,k2,ising(1:octetnum))
			octetmsinglet=0
			DO ijxj=1,octetnum
				IF(ising(octetnum+1-ijxj).EQ.1)octetmsinglet=octetmsinglet+2**(octetlist(ijxj)-1)
			ENDDO
			CALL Helac_rmatrix(iperm(1:ncc,1:nq))
		ENDDO
		octetmsinglet=0
	ENDDO

ELSE
! color octet is treated as the color nonet
	DO k1=1,ncc
		ijxi=0
		DO k2=1,nq
			colorpos=QuarkColorPos(k2)
			IF(ColorSinglet(k2,k2).AND.colorpos.NE.0)THEN
				IF(iperm(k1,colorpos).NE.k2)THEN
					DO ijxj=1,nq
						IF(iperm(k1,ijxj).EQ.k2)THEN
							swap=ijxj
							EXIT
						ENDIF
					ENDDO
					iperm(k1,swap)=iperm(k1,colorpos)
					iperm(k1,colorpos)=k2
					ijxi=ijxi+1
				ENDIF
			ENDIF
		ENDDO
		divide(k1)=dnou(NCOL)**(-ijxi)
	ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 3 and 3~ representation has been up<->down
	DO k1=1,ncc
		DO k2=k1,ncc
!returns sum(C_i*C_j), wher C_i is the color configuration,(rmatrix)i,j=C_i*C_j
			CALL Helac_eat(nq,iperm(k1,1:nq),iperm(k2,1:nq),izeros)
!	   CALL Helac_eat(2,(/1,2/),(/1,2/),i)
!	   WRITE(*,*)i
			rmatrix(k1,k2)=dnou(NCOL)**izeros*divide(k1)*divide(k2)
			IF(k2.GT.k1)rmatrix(k2,k1)=rmatrix(k1,k2)
		ENDDO
	ENDDO
ENDIF
       
!WRITE(nunit1,*)'ncc,ic1',ncc,ic1
!*********************************************************
!      end of initialization phase      
!*********************************************************
WRITE(*,*)' '
WRITE(*,*)'The number of colour configurations',ncc
WRITE(*,*)'The colour summation flag(0 specfic color,1 MC selection): imc=',imc
IF(ncc.LE.24)THEN
   WRITE(nunit1,*)'The colour matrix: '
   DO i=1,ncc
       WRITE(nunit1,'(24F8.2)')(rmatrix(i,j),j=1,ncc)
   ENDDO
ENDIF
!   DO k=1,ncc
!	  WRITE(*,*)jcol(k,1,0)
!     WRITE(*,*)jcol(k,1:nhad,1)   ! returns the colors with 1,2
!	  WRITE(*,*)jcol(k,1:nhad,2)
!	  WRITE(*,*)"------------------"
!   ENDDO
!   STOP
!WRITE(*,*)(ifl(i),i=1,n) ! Debug
WRITE(*,*)' '
! there is some problem with ng(public?) and is(allocate?public?) in Feynman_Helac
! ng is the total number of feynman graphs       
IF(repeat.EQ.1)THEN
   WRITE(21,*)ng
   DO ii=1,ng
      DO jj=1,8
          WRITE(21,*)(is(ii,jj,kk),kk=1,8)
      ENDDO
   ENDDO
   WRITE(21,*)avhel,avcol,symet
!  WRITE(21,*)'The colur matrix:'
   WRITE(21,*)rmatrix
   DO k=1,ncc
      WRITE(21,*)jcol(k,1:nhad,0:2)   ! returns the colors with 1,2
   ENDDO
ENDIF  

!IF(ifl(1).EQ.3.AND.ifl(2).EQ.3.AND.ifl(3).EQ.35.AND.ifl(4).EQ.35)PAUSE ! Debug
END SUBROUTINE Helac_init
       
SUBROUTINE Helac_master_f(smel2)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(OUT)::smel2
LOGICAL::lmode
INTEGER::k,ihel,istop,inegl,k1,k2,k_1,k_un,i,istat,inum,icut
REAL(KIND=DBL)::w1,beta0,beta1,wnc,r,alpha
COMPLEX(KIND=DBL)::CGfactor
IF(GLOBALINIT_master.EQ.0)THEN
	init=0
ENDIF
GLOBALINIT_mader=0
IF(init.EQ.0)THEN
	WRITE(*,*)'NCC=',ncc,' imc=',imc
	IF(ALLOCATED(smhelcol))THEN
		DEALLOCATE(smhelcol)
	ENDIF
	CALL Gethelnum(num)
	ALLOCATE(smhelcol(num,ncc),STAT=istat)
	IF(ihel1.NE.ihel2)THEN
		IF(ALLOCATED(smhelcol2))DEALLOCATE(smhelcol2)
		ALLOCATE(smhelcol2(num,ncc),STAT=istat)
	ENDIF
    IF(istat.NE.0)WRITE(*,*)istat,&
      'warning: allocation(1) is not working properly in Helac_master_f'
ENDIF
IF(init.EQ.0.AND.repeat.EQ.2)THEN
	IF(ALLOCATED(zamp))THEN
		DEALLOCATE(zamp)
	ENDIF
	IF(ALLOCATED(rmatrix))THEN
		DEALLOCATE(rmatrix)
	ENDIF
	IF(ALLOCATED(jcol))THEN
		DEALLOCATE(jcol)
	ENDIF
    ALLOCATE(zamp(ncc),rmatrix(ncc,ncc),jcol(ncc,nhad,0:2),STAT=istat)
    IF(istat.NE.0)WRITE(*,*)istat,&
      'warning: allocation(2) is not working properly in Helac_master_f'
    READ(21,*)avhel,avcol,symet
    READ(21,*)rmatrix
    DO k=1,ncc
        READ(21,*)jcol(k,1:nhad,0:2)
    ENDDO
ENDIF
smhelcol(1:num,1:ncc)=DCMPLX(0d0)
IF(ihel1.NE.ihel2)smhelcol2(1:num,1:ncc)=DCMPLX(0d0)
ihel=0
istop=0
IF(gener.NE.3.OR.iranhel.EQ.0)THEN
	Qhelran(1:20)=-1
	QSzran(1:20)=-1
	PhyPolNum=0
ENDIF
IF(.NOT.MCoHelicity)THEN
DO WHILE(istop.EQ.0)
   inegl=0
   CALL Helac_sethel(iflag,istop,ihel,inegl)  ! get every helicity configuration
   GLOBALINIT_Sethel=1
   inegl=1
   CALL Helac_iniqq()   ! initilization and random helicity selection 
                        ! for the external legs in Helac_pan1.f90
   IF(.NOT.exp3pjQ.OR.Num3pj.EQ.0)THEN
		IF(imode.EQ.0)THEN
			CALL PhysicalPol(inum)
			GLOBALINIT_physicalpol=1
			CGfactor=1
		ELSEIF(imode.EQ.1)THEN
			CALL PhysicalPol2(inum,CGfactor)
			GLOBALINIT_physicalpol2=1
		ENDIF
   ELSE
		CALL ls2j(inum,CGfactor)
		GLOBALINIT_ls2j=1
   ENDIF
   IF(imc.EQ.0)THEN
      DO k=1,ncc
         CALL Helac_nextq(k)
         CALL Helac_ampq(k,zamp(k))     ! get the amplitudes
		 zamp(k)=zamp(k)*CGfactor
		 IF(imode.EQ.1.AND..NOT.FirstHelicityQ)THEN
			smhelcol2(inum,k)=smhelcol2(inum,k)+zamp(k)
		 ELSE
			smhelcol(inum,k)=smhelcol(inum,k)+zamp(k)
		 ENDIF
!		 WRITE(*,*)zamp(k)  ! Debug
         IF(onep)THEN 
           WRITE(100,*)'the ',k,'-th color cnfg out of ',ncc
           IF(ncc.GT.1)WRITE(100,*)'the     color',jcol(k,1:nhad,1)
           IF(ncc.GT.1)WRITE(100,*)'the anticolor',jcol(k,1:nhad,2)
           IF(iranhel.EQ.0)WRITE(100,*)'the helicity',ipol(1:n)
           WRITE(100,*)'the amplitude',zamp(k)
         ENDIF
      ENDDO
	  ! color matrix rmatrix is positive definited
      w1=dnou(0)
      DO k1=1,ncc
         DO k2=k1,ncc
            IF(k2.EQ.k1)w1=w1+DREAL( zamp(k1)*DCONJG(zamp(k2)))*rmatrix(k1,k2)!matrix elements 1
			!          include 'largeNclimit1.h'
            IF(k2.GT.k1)w1=w1+2*DREAL( zamp(k1)*DCONJG(zamp(k2)))*rmatrix(k1,k2) 
         ENDDO
      ENDDO
    ELSE
	! color random selection,with Leading Colour Approximation
      k_1=INT(Helac_rnmy(0)*ncc)+1
      CALL Helac_nextq(k_1)
      CALL Helac_ampq(k_1,zamp(k_1))
	  zamp(k_1)=zamp(k_1)*CGfactor
	  smhelcol(inum,k_1)=smhelcol(inum,k_1)+zamp(k_1)*DSQRT(DBLE(ncc))

      w1=DREAL( zamp(k_1)*DCONJG(zamp(k_1)))*rmatrix(k_1,k_1)
      w1=w1*ncc
    ENDIF

! unweight over color
! random color selection
!	IF(n.EQ.nhad)THEN
!		beta0=0d0
!		k_un=1
!		r=Helac_rnmy(0)      ! choose some color configuration
!		wnc=0d0
!		DO k=1,ncc
!			IF(jcol(k,1,0).EQ.1)wnc=wnc+DREAL( zamp(k)*DCONJG(zamp(k) )) ! no U(1) external gluon
!		ENDDO
!		DO k=1,ncc
!			alpha=DREAL( zamp(k)*DCONJG(zamp(k) ))/wnc
!			IF(jcol(k,1,0).EQ.0)alpha=0
!			beta1=beta0+alpha
!			IF(r.GT.beta0.AND.r.LT.beta1)THEN
!				k_un=k
!				EXIT          ! goto 1
!			ENDIF
!			beta0=beta1
!		ENDDO
!		icol_un(1:nhad,1:2)=jcol(k_un,1:nhad,1:2)  ! 1
!	ENDIF
!    smel2=smel2+w1
!	WRITE(*,*)smel2  ! Debug
    IF(w1.EQ.dnou(0))inegl=2             ! delete this helicity configuration
    IF(init.EQ.0)CALL Helac_sethel(iflag,istop,ihel,inegl) ! delete this helicity configuration
    IF(init.EQ.0)THEN
	    WRITE(nunit1,*)' '
	    WRITE(nunit1,*)'polarization states for external legs'
        WRITE(nunit1,*)(ipol(i),i=1,nhad)
		WRITE(nunit1,*)'weighting w1='
        WRITE(nunit1,*)w1
		WRITE(nunit1,*)' '
    ENDIF
ENDDO
ELSE
   IF(init.EQ.0)THEN
      istop=0
      DO WHILE(istop.EQ.0)
         inegl=0
         CALL Helac_sethel2(iflag,istop,ihel,inegl,0)
         inegl=1
         CALL Helac_iniqq()
         IF(imc.EQ.0)THEN
            DO k=1,ncc
               CALL Helac_nextq(k)
               CALL Helac_ampq(k,zamp(k))
            ENDDO
            w1=dnou(0)
            DO k1=1,ncc
               DO k2=k1,ncc
                  IF(k2.EQ.k1)w1=w1+DREAL( zamp(k1)*DCONJG(zamp(k2)))*rmatrix(k1,k2)            
                  IF(k2.GT.k1)w1=w1+2*DREAL( zamp(k1)*DCONJG(zamp(k2)))*rmatrix(k1,k2)
               ENDDO
            ENDDO
         ELSE
            k_1=INT(Helac_rnmy(0)*ncc)+1
            CALL Helac_nextq(k_1)
            CALL Helac_ampq(k_1,zamp(k_1))
            w1=DREAL( zamp(k_1)*DCONJG(zamp(k_1)))*rmatrix(k_1,k_1)
            w1=w1*ncc
         ENDIF
         IF(w1.EQ.dnou(0))inegl=2
         IF(init.EQ.0)CALL Helac_sethel2(iflag,istop,ihel,inegl,0)
         WRITE(nunit1,*)' '
         WRITE(nunit1,*)'polarization states for external legs'
         WRITE(nunit1,*)(ipol(i),i=1,nhad)
         WRITE(nunit1,*)'weighting w1='
         WRITE(nunit1,*)w1
         WRITE(nunit1,*)' '
      ENDDO
      WRITE(*,*)'THE NONVANISHING HELICITY CFGS',nhc,nphyhc
   ENDIF
   !IF(gener.NE.3)THEN
   !   RMCoH=Helac_rnmy(0)
   !   DO k=1,NDecayIflh
   !      Decayran(k)=Helac_rnmy(0)
   !   ENDDO
   !ENDIF
   istop=0
   ihel=0
   weight_br=1d0
   DO WHILE(istop.EQ.0)
      inegl=0
      CALL Helac_sethel2(iflag,istop,ihel,inegl,1)
      inegl=1
      IF(ihel.EQ.1.AND.NDecayChains.GT.0)THEN
         CALL HO_Decay(weight_br)
         IF(weight_br.LE.0d0)THEN
            smel2=dnou(0)
            init=init+1
            RETURN
         ENDIF
         CALL cuts_Decay(icut)
         IF(icut.EQ.0)THEN
            smel2=dnou(0)
            init=init+1
            RETURN
         ENDIF
      ENDIF
      CALL Helac_iniqq()
      IF(.NOT.exp3pjQ.OR.Num3pj.EQ.0)THEN
         CALL PhysicalPol(inum)
         CGfactor=1
      ELSE
         PRINT *,"ERROR:iranhel=4 or decay processes cannot have explicit 3PJ"
         STOP
         CALL ls2j(inum,CGfactor)
         GLOBALINIT_ls2j=1
      ENDIF
      IF(imc.EQ.0)THEN
         DO k=1,ncc
            CALL Helac_nextq(k)
            CALL Helac_ampq(k,zamp(k))     ! get the amplitudes
            zamp(k)=zamp(k)*CGfactor*DSQRT(DBLE(nphyhc))
            smhelcol(inum,k)=smhelcol(inum,k)+zamp(k)
         ENDDO
      ELSE
         k_1=INT(Helac_rnmy(0)*ncc)+1
         CALL Helac_nextq(k_1)
         CALL Helac_ampq(k_1,zamp(k_1))
         zamp(k_1)=zamp(k_1)*CGfactor*DSQRT(DBLE(nphyhc))*DSQRT(DBLE(ncc))
         smhelcol(inum,k_1)=smhelcol(inum,k_1)+zamp(k_1)
      ENDIF
   ENDDO
ENDIF
! sigma
smel2=dnou(0)
zamp(1:ncc)=0
IF(imode.EQ.0.OR.(imode.EQ.1.AND.ihel1.EQ.ihel2))THEN
	DO i=1,num
		DO k1=1,ncc
			DO k2=k1,ncc
				IF(k1.EQ.k2)smel2=smel2+DREAL(smhelcol(i,k1)*DCONJG(smhelcol(i,k2)))*rmatrix(k1,k2)
				IF(k2.GT.k1)smel2=smel2+&
				2d0*DREAL(smhelcol(i,k1)*DCONJG(smhelcol(i,k2)))*rmatrix(k1,k2)
			ENDDO
			zamp(k1)=zamp(k1)+smhelcol(i,k1)*DCONJG(smhelcol(i,k1))
		ENDDO
	ENDDO
ELSEIF(imode.EQ.1.AND.ihel1.NE.ihel2)THEN
	DO i=1,num
		DO k1=1,ncc
			DO k2=1,ncc
				smel2=smel2+DREAL(smhelcol(i,k1)*DCONJG(smhelcol2(i,k2)))*rmatrix(k1,k2)
			ENDDO
			zamp(k1)=zamp(k1)+smhelcol(i,k1)*DCONJG(smhelcol(i,k1))
		ENDDO
	ENDDO
	IF(SDPart.GT.0)THEN
		smel2=MAX(0d0,smel2)
	ELSEIF(SDPart.LT.0)THEN
		smel2=MAX(0d0,-smel2)
	ENDIF
ENDIF
! unweight over color
! random color selection
beta0=0d0
k_un=1
r=Helac_rnmy(0)      ! choose some color configuration
wnc=0d0
DO k=1,ncc
	IF(jcol(k,1,0).EQ.1)wnc=wnc+DREAL(zamp(k)) ! no U(1) external gluon
ENDDO
DO k=1,ncc
	alpha=DREAL(zamp(k))/wnc
	IF(jcol(k,1,0).EQ.0)alpha=0
	beta1=beta0+alpha
	IF(r.GT.beta0.AND.r.LT.beta1)THEN
		k_un=k
		EXIT          ! goto 1
	ENDIF
	beta0=beta1
ENDDO
icol_un(1:nhad,1:2)=jcol(k_un,1:nhad,1:2)  ! 1

!     end of polariazation sum if any
IF(init.EQ.0.AND..NOT.MCoHelicity)WRITE(*,*)'THE NONVANISHING HELICITY CFGS',nhc
!     *****************************************************
!
!WRITE(*,*)smel2,avhel,avcol,symet  ! Debug
smel2=smel2*avhel*avcol*symet
!WRITE(*,*)'total',smel2   ! Debug
!
!     *****************************************************
IF(onep)THEN
    WRITE(100,*)'total amplitude squared = ',smel2
    WRITE(100,*)'----------------------------------------'
ENDIF
      
init=init+1 
IF(onep)init=0
END SUBROUTINE Helac_master_f
        
SUBROUTINE Helac_clean(iflag,zamp,n1)
INTEGER,INTENT(IN)::n1
INTEGER,DIMENSION(n1),INTENT(OUT)::iflag
COMPLEX(KIND=DBL),DIMENSION(n1),INTENT(IN)::zamp
INTEGER::k,l
DO k=2,n1
  DO l=1,k-1
     IF(zamp(l).EQ.zamp(k))iflag(k)=l
     IF(zamp(l).EQ.-zamp(k))iflag(k)=-l
  ENDDO
ENDDO
END SUBROUTINE Helac_clean
      
SUBROUTINE Helac_checkgluons(iq,nq,icheck)
INTEGER,INTENT(IN)::nq
INTEGER,DIMENSION(nq),INTENT(IN)::iq
INTEGER,INTENT(OUT)::icheck
INTEGER::igo=1,init_2=0
INTEGER::i,ix
SAVE igo,init_2
icheck=1
IF(GLOBALINIT_checkgluons.EQ.0)THEN
	igo=1
	init_2=0
ENDIF
IF(igo.EQ.0)RETURN
IF(init_2.EQ.0)THEN
   igo=0
   init_2=1
   DO i=1,n
      IF(ifl(i).NE.35)RETURN
   ENDDO
   igo=1
ENDIF

icheck=0
i=1
DO ix=1,nq-2
   i=iq(i)
   IF(iq(i).EQ.1)RETURN
ENDDO
icheck=1
END SUBROUTINE Helac_checkgluons

SUBROUTINE Helac_rmatrix(iperm)
INTEGER,DIMENSION(ncc,nq),INTENT(IN)::iperm
INTEGER,DIMENSION(:,:),ALLOCATABLE::iperm22
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::divide22
INTEGER::init22=0,istat,k1,k2,ijxi,ijxj,colorpos,swap,izeros
SAVE init22,iperm22,divide22
IF(init22.EQ.0)THEN
	IF(ALLOCATED(iperm22))THEN
		DEALLOCATE(iperm22)
	ENDIF
	IF(ALLOCATED(divide22))THEN
		DEALLOCATE(divide22)
	ENDIF
	ALLOCATE(iperm22(ncc,nq),STAT=istat)
	ALLOCATE(divide22(ncc),STAT=istat)
	IF(istat.NE.0)WRITE(*,*)istat, &
	' warning 2: allocation is not working properly in Helac_rmatrix !'
	init22=1
ENDIF
iperm22(1:ncc,1:nq)=iperm(1:ncc,1:nq)
divide22(1:ncc)=0
DO k1=1,ncc
	ijxi=0
	DO k2=1,nq
		colorpos=QuarkColorPos(k2)
		IF(ColorSinglet(k2,k2).AND.colorpos.NE.0)THEN
			IF(iperm22(k1,colorpos).NE.k2)THEN
				DO ijxj=1,nq
					IF(iperm22(k1,ijxj).EQ.k2)THEN
						swap=ijxj
						EXIT
					ENDIF
				ENDDO
				iperm22(k1,swap)=iperm22(k1,colorpos)
				iperm22(k1,colorpos)=k2
				ijxi=ijxi+1
			ENDIF
		ENDIF
	ENDDO
	divide22(k1)=dnou(NCOL)**(-ijxi)
ENDDO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
! 3 and 3~ representation has been up<->down
ijxi=Helac_level(nhad,octetmsinglet)
ijxj=1
IF(MOD(ijxi,2).NE.0)ijxj=-1
DO k1=1,ncc
    DO k2=k1,ncc
!returns sum(C_i*C_j), wher C_i is the color configuration,(rmatrix)i,j=C_i*C_j
       CALL Helac_eat(nq,iperm22(k1,1:nq),iperm22(k2,1:nq),izeros)
!	   CALL Helac_eat(2,(/1,2/),(/1,2/),i)
!	   WRITE(*,*)i
       rmatrix(k1,k2)=rmatrix(k1,k2)+dnou(NCOL)**izeros*divide22(k1)*divide22(k2)*ijxj
       IF(k2.GT.k1)rmatrix(k2,k1)=rmatrix(k1,k2)
    ENDDO
ENDDO
END SUBROUTINE Helac_rmatrix
   
!      include 'special/master.h'
SUBROUTINE Helac_checki0(i0,l,nogo)
INTEGER,INTENT(IN)::i0,l
INTEGER,INTENT(OUT)::nogo
INTEGER,DIMENSION(n)::y
INTEGER::sum
nogo=0
RETURN

! for t-tbar
nogo=0
IF(l.LE.3)RETURN
nogo=1
CALL Helac_bin(n,i0,y)
IF(y(n-2)+y(n-3)+y(n-4).NE.0.AND.y(n-2)+y(n-3)+y(n-4).NE.3)RETURN
IF(y(n-5)+y(n-6)+y(n-7).NE.0.AND.y(n-5)+y(n-6)+y(n-7).NE.3)RETURN
nogo=0
RETURN
END SUBROUTINE Helac_checki0

END MODULE Helac_master
