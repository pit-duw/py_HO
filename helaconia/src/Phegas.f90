MODULE Phegas_mod
USE Helac_Global
USE Helac_Func_1
USE Helac_pan2
USE Kinetic_Func
USE MC_RAMBO
USE Phegas_Durham_mod
USE Helac_Feynman
USE Phegas_Choice
USE MC_Helac_GRID
USE Cuts_Module
USE Helac_histo
USE Adapt
USE QEDPS_interface
USE Decay_interface
!USE Constants
IMPLICIT NONE
REAL(KIND=DBL)::xx,r2,q,q2,e2,em
REAL(KIND=DBL),DIMENSION(1000)::x
REAL(KIND=DBL)::rq,rq2,re2,rem
REAL(KIND=DBL)::EBMUP1,EBMUP2
CONTAINS
SUBROUTINE Phegas(weight,ig,idir,iweight)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(OUT)::weight
INTEGER,INTENT(IN)::ig
INTEGER,INTENT(INOUT)::idir
INTEGER,INTENT(OUT)::iweight
! for strf use
CHARACTER(len=20),DIMENSION(20)::parm
REAL(KIND=DBL),DIMENSION(20)::val
INTEGER,DIMENSION(n)::y1,y2,y3
LOGICAL::lsym
REAL(KIND=DBL),DIMENSION(4)::p2,pdummy,pboo,pp1,pffp3
REAL(KIND=DBL),DIMENSION(100)::xm
REAL(KIND=DBL),DIMENSION(4,100)::pr
REAL(KIND=DBL),DIMENSION(:,:),ALLOCATABLE::p
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::cmas,tm
INTEGER::init=0,flag,flagif1,pdfnumpdf
REAL(KIND(1d0))::wgt_shower
! durham
REAL(KIND=DBL),DIMENSION(100,4)::p_durham ! durham
REAL(KIND=DBL)::r_0,r_1,gevnb,pi
INTEGER::ig0,i,k,istat,flag12,it1new=0,inum,k1,k2,k3,flagreturn,kk1,kk2
REAL(KIND=DBL)::t1max,t1min,a1,a2,w_q1,gm1,gg1,q1min,q1max,q2min,q2max,&
                q3max,q3min,q1,pq,w_q2,w_costh,w_phi,wtps,eram,&
				pt,tau0,he,weight_s,cha,w1,q3,w2,q23,w,t1l,t1u,cth,&
				t1,r,phi,sth,pm,cth0,phi0,t1new,w_lam,w3,e,tau,taulog,wght,dy,&
				t2max,t2min,phi02,ptcut
INTEGER::ii1,mm1,ii2,mm2,i1,i2,i3,m1,m2,m3,itogo,iq2max,j,i0,k0,j1,j2,j3,m0,mm10,idir0,icase0,&
         icase,temptt,mm20,mm30,flagphi,nchchmax
LOGICAL::adaptlog
INTEGER,DIMENSION(n)::yykk
SAVE
r_0=0
r_1=1
weight=1
iweight=1
IF(idir.EQ.0)ig0=ig
IF(GLOBALINIT_Phegas.EQ.0)THEN
	init=0
	it1new=0
ENDIF
IF(init.EQ.0)THEN
	WRITE(*,*)' '
    WRITE(*,*)'-------------------------------------------------------------------'
    WRITE(*,*)'THE PHEGAS SUBROUTINE'

	ptc(3:20)=0                   ! Debug
	drc(3:20,3:20)=0
	etac(3:20)=20
	ec(3:20)=0
	c1(3:20)=1
	c2(3:20)=1
	cc(3:20,3:20)=1
	gmas(3:20,3:20)=0             ! Debug

    parmas(0)=0
    parwid(0)=-1
    iwarning(1:25)=0
    iwonders(1:25)=0
! ------------------------------------------------------------

    IF(ng.LT.100.AND.ng.GT.0)THEN
	    WRITE(nunit1,*)' '
		WRITE(nunit1,*)'The information for graphs'
        DO i=1,ng
            k=1
            WRITE(nunit1,*)'The graph',i
            DO WHILE(is(i,k,1).NE.0)
			   WRITE(nunit1,*)' '
			   WRITE(nunit1,*)'         i0,         m0'
               WRITE(nunit1,*)is(i,k,1:2)
			   WRITE(nunit1,*)' '
			   WRITE(nunit1,*)'         i1,         m1,         i2,          m2,        i3,         m3'
			   WRITE(nunit1,*)is(i,k,3:8)
			   WRITE(nunit1,*)' '
			   WRITE(nunit1,*)' '
               k=k+1
            ENDDO
        ENDDO
     ENDIF
! ------------------------------------------------------------
     gevnb=dnou(38937966)/100   ! the conservation factor of GeV^-2 to nb
     CALL Helac_mypi(pi)
	 WRITE(*,*)' '
	 WRITE(*,*)'TWO CONSTANTS'
     WRITE(*,*)'The unit GeV^-2 to nb gevnb=',gevnb
     WRITE(*,*)'pi=',pi
!  ---- setup incoming momenta
	 IF(ALLOCATED(p))THEN
		DEALLOCATE(p)
	 ENDIF
	 IF(ALLOCATED(cmas))THEN
		DEALLOCATE(cmas)
	 ENDIF
	 IF(ALLOCATED(tm))THEN
		DEALLOCATE(tm)
	 ENDIF
     ALLOCATE(p(0:2**n-2,4),cmas(0:2**n-2),tm(0:2**n-2),STAT=istat)
     IF(istat.NE.0)WRITE(*,*)'warning: in phegas entry',istat
     cmas(0:2**n-2)=0
     tm(0:2**n-2)=0   
	 WRITE(*,*)' ' 
!	 CALL ReadElem_real('energy',e)
     CALL ReadElem_real('energy_beam1',e)
     ebeam(1)=ABS(e)
     CALL ReadElem_real('energy_beam2',e)
     ebeam(2)=ABS(e)
     IF(ebeam(1).LT.parmas(ifl(1)))STOP "energy of beam 1 is too small"
     IF(ebeam(2).LT.parmas(ifl(2)))STOP "energy of beam 2 is too small"
     IF(ABS(ebeam(1)-ebeam(2))/MAX(ebeam(1)+ebeam(2),1d-17).LT.1d-8)THEN
        absrap=absrap
        labeqcoll=.TRUE.
     ELSE
        IF(emep_ISR_shower.EQ.0)THEN
           absrap=.FALSE.
           labeqcoll=.FALSE.
        ELSE
           absrap=.FALSE.
           labeqcoll=.TRUE. ! run in c.m. frame
        ENDIF
     ENDIF
     e=2d0*ebeam(1)*ebeam(2)+2d0*&
          DSQRT(ebeam(1)**2-parmas(ifl(1))**2)*&
          DSQRT(ebeam(2)**2-parmas(ifl(2))**2)+&
          parmas(ifl(1))**2&
          +parmas(ifl(2))**2
     e=DSQRT(e)
     IF(fixtarget)e=e/DSQRT(2d0)
     IF(fixtarget)THEN
        IF(ABS(ebeam(1)).GT.ABS(ebeam(2)))THEN
           fixtargetrev=.FALSE.
           FT_M=ebeam(2)
           FT_E1=ebeam(1)
           e=DSQRT(e**2+FT_M**2)
        ELSE
           fixtargetrev=.TRUE.
           FT_M=ebeam(1)
           FT_E1=ebeam(2)
           e=DSQRT(e**2+FT_M**2)
        ENDIF
     ENDIF
     WRITE(*,*)'total energy = ',e
     EBMUP1=ebeam(1)
     EBMUP2=ebeam(2)
!     READ(*,*)e   ! Sqrt(S)
	 CALL ReadElem_integer('pdf',pdfnumpdf)
	 IF(pdfnumpdf.EQ.0.OR.COLL_TYPE.EQ.3)THEN
		istruc=0
	 ELSE
		istruc=1
	 ENDIF
     WRITE(*,*)'struc f (1 with structure functions while 0 not) ? ',istruc
 !    READ(*,*)istruc
 !    WRITE(*,*)'istruc = ',istruc
!	 STOP
     IF(emep_ISR_shower.EQ.1)CALL QEDPSINIT
     xp1=1
     xp2=1
     he=e
     ehat=he
     p(0:2**n-2,1:4)=0
	 ! the collision c.m. momenta of incoming two particles
     p(1,4)=(e**2+parmas(ifl(1))**2-parmas(ifl(2))**2)/e/2
     p(1,3)=p(1,4)*DSQRT( 1-(parmas(ifl(1))/p(1,4))**2 )
     p(1,2)=0
     p(1,1)=0
     p(2,4)=(e**2+parmas(ifl(2))**2-parmas(ifl(1))**2)/e/2
     p(2,3)=-p(2,4)*DSQRT( 1-(parmas(ifl(2))/p(2,4))**2 )
     p(2,2)=0
     p(2,1)=0
	 DO i=1,2
       	j=2**(i-1)
       	pt=DSQRT(p(j,1)**2+p(j,2)**2)
       	pq=DSQRT(pt**2+p(j,3)**2)
       	zq(i,1)= DCMPLX(p(j,4)+p(j,3),p(j,3))
       	zq(i,2)= DCMPLX(p(j,4)-p(j,3),pt)
       	zq(i,3)= DCMPLX(p(j,1), p(j,2))
       	zq(i,4)= DCMPLX(p(j,1),-p(j,2))
       	zq(i,5)= DCMPLX( parmas(ifl(i)) , pq )
       	phegas_pmom(i,1:4)=p(j,1:4)
       	phegas_pmom(i,5)=parmas(ifl(i))
	 ENDDO
     wrest=4*DSQRT(scalar_product(p(1,1:4),p(2,1:4))**2&
           -parmas(ifl(1))**2*parmas(ifl(2))**2   )   ! 4*moller's invariant flux factor
     wrest=gevnb/wrest
! the factor is 1/(2*s) and not 1/(2*shat)
! the difference being accounted for in wjac (x1x2gen.h)
     IF(COLL_TYPE.NE.3.AND.istruc.EQ.1)THEN
		CALL readcuts_HC
	 ELSEIF(COLL_TYPE.EQ.3)THEN
		CALL readcuts_epem
	 ENDIF
     IF(NDecayChains.GT.0)CALL readcuts_Decay
     CALL Phegas_setmin(cmas(1:2**n-2))
	 GLOBALINIT_setmin=1
     CALL Phegas_setmax(tm(1:2**n-2),cmas(1:2**n-2),xp1,xp2)
	 GLOBALINIT_setmax=1
     DO i=1,2**n-2
        IF(cmas(i).NE.0)WRITE(nunit2,*)'cmas',i,cmas(i) 
     ENDDO
     DO i=1,2**n-2
        IF(tm(i).NE.0)WRITE(nunit2,*)'tmas',i,tm(i) 
     ENDDO

!  --- 
!  define minimum energy
!  ---
     IF(istruc.EQ.1)THEN      ! istruc=1 means it is a proton (anti-)proton collider
	                          ! istruc=0 means it is a electron-positron collider
!        include 'phegasstrf.h'
		iPDFGUP1=0 !4  ! PDF group according to the Cernlib PDFlib
        iPDFGUP2=0 !4
        iPDFSUP1=pdfnumpdf !48   CTEQ6L1
        iPDFSUP2=pdfnumpdf !48   ! PDF set ID for beam 1 and 2, according to the Cernlib PDFlib
!       end of 'phegasstrf.h'
!       include 'qcdscale.h'
        tau0=(cmas(2**n-4)/e)**2
        IF(cmas(2**n-4).LT.0.1d0)tau0=(0.1d0/e)**2
        WRITE(*,*)'tau0 = ',tau0,cmas(2**n-4)
	 ELSE
		iPDFGUP1=-1
		iPDFGUP2=-1
		iPDFSUP1=-1
		iPDFSUP2=-1
     ENDIF
	 CALL ReadElem_logic("adapt",adaptlog)
     init=1
ENDIF     
! -----
!  structure functions
! -----
IF(istruc.EQ.1.AND.idir.EQ.0)THEN      !start STRF
       
!       r1=rnmy(0)
!       r2=rnmy(0)
	IF(.NOT.adaptlog)THEN

!       include 'x1x2gen.h'
		wjac=1
!		IF(nhad.EQ.4.AND.n.EQ.6.AND.dSigmadPTQ)THEN			
!		ELSE
		IF(init.EQ.1)THEN
			CALL ReadElem_integer('grid_nchmax',nchchmax)
			CALL helac_grid_alloc_2( 2 )
			CALL helac_grid_init_2( 1 ,500,nchchmax,0d0 ,0,1 )
			CALL helac_grid_init_2( 2 ,500,nchchmax,0d0 ,0,2 )
			init=2
		ENDIF

		CALL helac_grid_gnrt_2( 1 ,xx )

!      print*,'tau0=',tau0,dlog(1d0/tau0)
		taulog=xx*DLOG(1d0/tau0)
		tau=DEXP(-taulog)
		CALL helac_grid_wght_2( 1 ,wght ,xx )
		wjac=wjac*wght*DLOG(1d0/tau0)
   

		CALL helac_grid_gnrt_2( 2 ,r2 )
		dy=DLOG(1d0/tau)*(2*r2-1)
		CALL helac_grid_wght_2( 2 ,wght ,r2 )
		wjac=wjac*2*DLOG(1d0/tau)*wght

		xp1=DSQRT(tau*DEXP( dy))
		xp2=DSQRT(tau*DEXP(-dy))
		wjac=wjac/2d0
!		ENDIF


		ehat=DSQRT(xp1*xp2)*e
	ELSE	
		wjac=1

		CALL adapt_gen(15,0d0,DLOG(1d0/tau0),taulog,1)  ! 0<=ln(1/xp1/xp2)<=ln(1/tau0)
		GLOBALINIT_adapt_gen=1
		tau=DEXP(-taulog)                               ! xp1*xp2

		CALL adapt_gen(10,0d0,1d0,r2,2)
		dy=DLOG(1d0/tau)*(2*r2-1)                       ! ln(xp1/xp2)
		wjac=wjac*2*DLOG(1d0/tau)

		xp1=DSQRT(tau*DEXP( dy))
		xp2=DSQRT(tau*DEXP(-dy))
		wjac=wjac/2d0

		ehat=DSQRT(xp1*xp2)*e


!       end of 'x1x2gen.h'
	ENDIF
!      include 'qcdscale.h'

! -----
! in the parton rest frame
   he=DSQRT(xp1*xp2)*e   ! Sqrt(shat)
   p(1,4)=(he**2+parmas(ifl(1))**2-parmas(ifl(2))**2)/he/2
   p(1,3)=p(1,4)*DSQRT( 1-(parmas(ifl(1))/p(1,4))**2 )
   p(1,2)=0
   p(1,1)=0
   p(2,4)=(he**2+parmas(ifl(2))**2-parmas(ifl(1))**2)/he/2
   p(2,3)=-p(2,4)*DSQRT( 1-(parmas(ifl(2))/p(2,4))**2 )
   p(2,2)=0
   p(2,1)=0
   DO i=1,2
      j=2**(i-1)
      pt=DSQRT(p(j,1)**2+p(j,2)**2)
      pq=DSQRT(pt**2+p(j,3)**2)
      zq(i,1)= DCMPLX(p(j,4)+p(j,3),p(j,3))
      zq(i,2)= DCMPLX(p(j,4)-p(j,3),pt)
      zq(i,3)= DCMPLX(p(j,1), p(j,2))
      zq(i,4)= DCMPLX(p(j,1),-p(j,2))
      zq(i,5)= DCMPLX( parmas(ifl(i)) , pq )
      phegas_pmom(i,1:4)=p(j,1:4)
      phegas_pmom(i,5)=parmas(ifl(i))
   ENDDO    
   CALL Phegas_setmax(tm(1:2**n-2),cmas(1:2**n-2),xp1,xp2)
! -----
ELSEIF(emep_ISR_shower.EQ.1.AND.idir.EQ.0)THEN
   CALL QEDPS_ISR_SHOWER
   IF(Q2OUT_QEDPS.LT.cmas(2**n-4)**2)THEN
      iweight=0
      RETURN
   ENDIF
   ebeam(1)=DSQRT(Q2OUT_QEDPS)/2d0
   ebeam(2)=ebeam(1)
   e=2d0*ebeam(1)*ebeam(2)+2d0*&
        DSQRT(ebeam(1)**2-parmas(ifl(1))**2)*&
        DSQRT(ebeam(2)**2-parmas(ifl(2))**2)+&
        parmas(ifl(1))**2&
        +parmas(ifl(2))**2
   e=DSQRT(e)
   he=e
   ehat=e
   p(0:2**n-2,1:4)=0
         ! the collision c.m. momenta of incoming two particles
   p(1,4)=(e**2+parmas(ifl(1))**2-parmas(ifl(2))**2)/e/2
   p(1,3)=p(1,4)*DSQRT( 1-(parmas(ifl(1))/p(1,4))**2 )
   p(1,2)=0
   p(1,1)=0
   p(2,4)=(e**2+parmas(ifl(2))**2-parmas(ifl(1))**2)/e/2
   p(2,3)=-p(2,4)*DSQRT( 1-(parmas(ifl(2))/p(2,4))**2 )
   p(2,2)=0
   p(2,1)=0
   DO i=1,2
      j=2**(i-1)
      pt=DSQRT(p(j,1)**2+p(j,2)**2)
      pq=DSQRT(pt**2+p(j,3)**2)
      zq(i,1)= DCMPLX(p(j,4)+p(j,3),p(j,3))
      zq(i,2)= DCMPLX(p(j,4)-p(j,3),pt)
      zq(i,3)= DCMPLX(p(j,1), p(j,2))
      zq(i,4)= DCMPLX(p(j,1),-p(j,2))
      zq(i,5)= DCMPLX( parmas(ifl(i)) , pq )
      phegas_pmom(i,1:4)=p(j,1:4)
      phegas_pmom(i,5)=parmas(ifl(i))
   ENDDO
   CALL Phegas_setmax(tm(1:2**n-2),cmas(1:2**n-2),xp1,xp2)
   wrest=4*DSQRT(scalar_product(p(1,1:4),p(2,1:4))**2&
        -parmas(ifl(1))**2*parmas(ifl(2))**2   )   ! 4*moller's invariant flux factor
   wrest=gevnb/wrest
ENDIF                             !end the STRF
! -----
p(2**n-2,1:4)=p(1,1:4)  ! the conservation of momentums
IF(idir.EQ.1.AND.(imode.EQ.1.OR.NDecayChains.NE.0))THEN
   IF(emep_ISR_shower.NE.1)THEN
      CALL labtocms()
   ELSE
      CALL labtocms_QEDPS()
   ENDIF
ENDIF

!  ---- choose a graph  : ig

IF(ig.EQ.0)THEN
!IF(ig.EQ.-1)THEN
   IF(dSigmadPTQ)STOP
   IF(istruc.EQ.1)eram=he    !
   IF(istruc.EQ.0)eram=e
   inum=3
   DO i=3,nhad
	  IF(ABS(iflh(i)).GT.100)THEN
		xm(i-2)=parmas(ifl(inum))+parmas(ifl(inum+1))
		IF(idir.EQ.1)pr(1:4,i-2)=phegas_pmom(inum,1:4)+phegas_pmom(inum+1,1:4)
		inum=inum+2
	  ELSE
		xm(i-2)=parmas(ifl(inum))
		IF(idir.EQ.1)pr(1:4,i-2)=phegas_pmom(inum,1:4)
		inum=inum+1
	  ENDIF
   ENDDO
   IF(i_psp.LE.1)CALL RAMBO(nhad-2,eram,xm,pr,wtps,idir)

   IF(idir.EQ.0)THEN
      IF(i_psp.LE.1)THEN
		  inum=3
          DO i=3,nhad
			 IF(ABS(iflh(i)).GT.100)THEN
				j=2**(inum-1)
				p(j,1:4)=parmas(ifl(inum))/(parmas(ifl(inum))+parmas(ifl(inum+1)))&
				*pr(1:4,i-2)   ! generate the momentum randomly
				j=2**(inum)
				p(j,1:4)=parmas(ifl(inum+1))/(parmas(ifl(inum))+parmas(ifl(inum+1)))&
				*pr(1:4,i-2)   ! generate the momentum randomly
				inum=inum+2
			 ELSE
				j=2**(inum-1)
				p(j,1:4)=pr(1:4,i-2)
				inum=inum+1
			 ENDIF
          ENDDO
      ENDIF

      IF(i_psp.EQ.2)THEN
		  IF(n.NE.nhad)THEN
			PRINT *,"Wrong in choosing the Phegas_durham as the exsistences of Quarkonia ! STOP !"
			STOP
		  ENDIF
          CALL Phegas_durham(n-2,p(1,4)+p(2,4),p(1,3)+p(2,3),p_durham)   ! durham 
		                                                !to generate massless momentums
          DO i=3,n
             pr(1:4,i-2)=p_durham(i-2,1:4)       ! durham
             j=2**(i-1)
             p(j,1:4)=pr(1:4,i-2)
          ENDDO
      ENDIF

      DO i=3,n
          j=2**(i-1)
          pt=DSQRT(p(j,1)**2+p(j,2)**2)
          pq=DSQRT(pt**2+p(j,3)**2)
          zq(i,1)= DCMPLX(p(j,4)+p(j,3),p(j,3))
          zq(i,2)= DCMPLX(p(j,4)-p(j,3),pt) 
          zq(i,3)= DCMPLX(p(j,1), p(j,2))
          zq(i,4)= DCMPLX(p(j,1),-p(j,2))
          zq(i,5)= DCMPLX( parmas(ifl(i)) , pq )
          phegas_pmom(i,1:4)=p(j,1:4)
          phegas_pmom(i,5)=parmas(ifl(i))
      ENDDO
    ENDIF
!	OPEN(UNIT=1930,FILE="finalmom.txt")
!	DO i=3,n
!		j=2**(i-1)
!		WRITE(1930,*)i-1,p(j,1:4)
!	ENDDO
!	STOP

    IF(i_psp.EQ.2)CALL Phegas_durham_wei(wtps)         ! durham
    weight=weight*wtps*(2*pi)**(4-3*(nhad-2))                 ! multiply the factor in the phase space
    IF(idir.EQ.0)CALL Phegas_testmom(iweight,xp1,xp2)
    IF((imode.EQ.1.OR.NDecayChains.NE.0).AND.(istruc.EQ.1.OR..NOT.labeqcoll).AND.emep_ISR_shower.EQ.0)THEN
       CALL cmstolab() ! in QEDPS, the boost is done in QEDPSEVNT
       CALL trans_p_to_zq(n,Phegas_pmom(1:n,1:5),zq(1:n,1:5))
    ENDIF
    RETURN
ENDIF
!  ---- start flying over the graph
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! generate the momenta of all hadrons
flagreturn=0       
! Generate q1^2 and q2^2 for hadron line and parton line
IF(istruc.EQ.1)eram=he
IF(istruc.EQ.0)eram=e
k2=0
DO k1=1,n-nhad
	k2=k2+2**(Quarkonium(k1,2)-1)+2**(Quarkonium(k1,3)-1)
ENDDO
IF(k2.GT.0)THEN
	IF(idir.EQ.1)THEN
		p(k2,1:4)=0
		DO k1=1,n-nhad
			p(k2,1:4)=p(k2,1:4)+p(2**(Quarkonium(k1,2)-1),1:4)+p(2**(Quarkonium(k1,3)-1),1:4)
		ENDDO
	ENDIF
	IF(2**n-4-k2.GT.0)THEN
		IF(idir.EQ.1)p(2**n-4-k2,1:4)=p(1,1:4)+p(2,1:4)-p(k2,1:4)
		q1min=cmas(k2)
		q2min=cmas(2**n-4-k2)
		q=eram
		q1max=q-q2min
		IF(Helac_level(n,k2).LE.2)q1max=q1min
		IF(q1max.GT.q1min)THEN
			IF(idir.EQ.1)THEN
				q1=scalar_product(p(k2,1:4),p(k2,1:4))
			ENDIF
			CALL Phegas_gen_had(idir,q1min,q1max,q1,w1)
			w1=w1/(2*pi)
			q1=DSQRT(q1)
		ELSE
			w1=1
			q1=q1min
		ENDIF
		IF(dSigmadPTQ.AND.Helac_level(n,k2).EQ.2)THEN
                       IF((q-DSQRT(q1**2+PTFirst**2))**2-PTFirst**2.LT.q2min**2&
                            .OR.q-q1.LT.q2min)THEN
                           iweight=-1
                           IF(idir.EQ.0)iwarning(25)=iwarning(25)+1
                           IF(idir.EQ.1)iwonders(25)=iwonders(25)+1
                           RETURN
                        ENDIF
			q2max=MIN(DSQRT((q-DSQRT(q1**2+PTFirst**2))**2-PTFirst**2),q-q1)
		ELSE
			q2max=q-q1
		ENDIF
		IF(Helac_level(n,2**n-4-k2).EQ.1)q2max=q2min
		IF(q2max.GT.q2min)THEN
			IF(idir.EQ.1)THEN
				q2=scalar_product(p(2**n-4-k2,1:4),p(2**n-4-k2,1:4))
			ENDIF
			CALL Phegas_gen_had(idir,q2min,q2max,q2,w2)
			w2=w2/(2*pi)
			q2=DSQRT(q2)
		ELSE
			w2=1
			q2=q2min
		ENDIF
! Boost to the initial parton rest frame
!		CALL Phegas_twobody(idir,q,q1,q2,p(k2,1:4),p(2**n-4-k2,1:4),w)
		IF(.NOT.dSigmadPTQ)THEN
			CALL Phegas_twobody(idir,q,q1,q2,p(k2,1:4),p(2**n-4-k2,1:4),w)
		ELSE
			IF(Helac_level(n,k2).EQ.2)THEN
				CALL Phegas_twobody_had_pt(idir,q,q1,q2,PTFirst,p(k2,1:4),p(2**n-4-k2,1:4),w)
			ELSE
				PRINT *,"Wrong in dSigma/dPt ! STOP !"
				STOP
			ENDIF
		ENDIF
		! allowed region but with kinematic cut, so iweight=-1
		IF(w*w1*w2.EQ.0)THEN
			iweight=-1
			IF(idir.EQ.0)iwarning(25)=iwarning(25)+1
			IF(idir.EQ.1)iwonders(25)=iwonders(25)+1
			RETURN
		ENDIF
		IF(idir.EQ.0)THEN
			pboo(1:4)=p(1,1:4)+p(2,1:4)
			pq=eram
			CALL boostl(pq,pboo,p(k2,1:4))
			CALL boostl(pq,pboo,p(2**n-4-k2,1:4))
		ENDIF
		weight=weight*w*w1*w2
	ELSE
		IF(idir.EQ.0)THEN
			p(k2,1:4)=p(1,1:4)+p(2,1:4)
		ENDIF
		q1=eram
		flagreturn=1
	ENDIF
	IF(Helac_level(n,k2).EQ.2)THEN
		IF(idir.EQ.0)THEN
			k1=Quarkonium(1,2)
			k3=Quarkonium(1,3)
			p(2**(k1-1),1:4)=parmas(ifl(k1))/(parmas(ifl(k1))+parmas(ifl(k3)))*p(k2,1:4)
			p(2**(k3-1),1:4)=parmas(ifl(k3))/(parmas(ifl(k1))+parmas(ifl(k3)))*p(k2,1:4)
		ENDIF
	ELSE
		q=q1
		pboo(1:4)=p(k2,1:4)
		k1=k2
		DO i=2,n-nhad
			k3=2**(Quarkonium(i-1,2)-1)+2**(Quarkonium(i-1,3)-1)
			IF(idir.EQ.1)p(k3,1:4)=p(2**(Quarkonium(i-1,2)-1),1:4)+p(2**(Quarkonium(i-1,3)-1),1:4)
			q1=cmas(k3)
			k1=k1-k3
			IF(idir.EQ.1)p(k1,1:4)=p(k1+k3,1:4)-p(k3,1:4)
			q2min=cmas(k1)
			q2max=q-q1
			IF(Helac_level(n,k1).LE.2)q2max=q2min
			IF(q2max.GT.q2min)THEN
				IF(idir.EQ.1)q2=scalar_product(p(k1,1:4),p(k1,1:4))
				CALL Phegas_gen_had(idir,q2min,q2max,q2,w2)
				w2=w2/(2*pi)
				q2=DSQRT(q2)
				IF(i.GT.2.OR.flagreturn.EQ.0.OR.Coll_Type.EQ.3)THEN
					CALL phegas_twobody(idir,q,q1,q2,p(k3,1:4),p(k1,1:4),w)
				ELSE
					IF(.NOT.dSigmadPTQ)THEN
						ptcut=ptc(Quarkonium(1,2))+ptc(Quarkonium(1,3))
						CALL phegas_twobody_had(idir,q,q1,q2,ptcut,p(k3,1:4),p(k1,1:4),w)
					ELSE
						CALL phegas_twobody_had_pt(idir,q,q1,q2,PTFirst,p(k3,1:4),p(k1,1:4),w)
					ENDIF
				ENDIF
				IF(w*w2.EQ.0)THEN
					iweight=-1
					IF(idir.EQ.0)iwarning(25)=iwarning(25)+1
					IF(idir.EQ.1)iwonders(25)=iwonders(25)+1
					RETURN						
				ENDIF
				IF(idir.EQ.0)THEN
					pq=q
					CALL boostl(pq,pboo,p(k3,1:4))
					CALL boostl(pq,pboo,p(k1,1:4))
					kk1=Quarkonium(i-1,2)
					kk2=Quarkonium(i-1,3)
					p(2**(kk1-1),1:4)=parmas(ifl(kk1))/(parmas(ifl(kk1))+parmas(ifl(kk2)))&
					*p(k3,1:4)
					p(2**(kk2-1),1:4)=parmas(ifl(kk2))/(parmas(ifl(kk1))+parmas(ifl(kk2)))&
					*p(k3,1:4)
				ENDIF
				weight=weight*w*w2
			ELSE
				w2=1
				q2=q2min
				IF(i.GT.2.OR.flagreturn.EQ.0.OR.Coll_Type.EQ.3)THEN
					CALL phegas_twobody(idir,q,q1,q2,p(k3,1:4),p(k1,1:4),w)
				ELSE
					IF(.NOT.dSigmadPTQ)THEN
						ptcut=ptc(Quarkonium(1,2))+ptc(Quarkonium(1,3))
						CALL phegas_twobody_had(idir,q,q1,q2,ptcut,p(k3,1:4),p(k1,1:4),w)
					ELSE
						CALL phegas_twobody_had_pt(idir,q,q1,q2,PTFirst,p(k3,1:4),p(k1,1:4),w)
					ENDIF
				ENDIF
				IF(w*w2.EQ.0)THEN
					iweight=-1
					IF(idir.EQ.0)iwarning(25)=iwarning(25)+1
					IF(idir.EQ.1)iwonders(25)=iwonders(25)+1
					RETURN						
				ENDIF
				IF(idir.EQ.0)THEN
					pq=q
					CALL boostl(pq,pboo,p(k3,1:4))
					CALL boostl(pq,pboo,p(k1,1:4))
					kk1=Quarkonium(i-1,2)
					kk2=Quarkonium(i-1,3)
					p(2**(kk1-1),1:4)=parmas(ifl(kk1))/(parmas(ifl(kk1))+parmas(ifl(kk2)))&
					*p(k3,1:4)
					p(2**(kk2-1),1:4)=parmas(ifl(kk2))/(parmas(ifl(kk1))+parmas(ifl(kk2)))&
					*p(k3,1:4)
					kk1=Quarkonium(i,2)
					kk2=Quarkonium(i,3)
					p(2**(kk1-1),1:4)=parmas(ifl(kk1))/(parmas(ifl(kk1))+parmas(ifl(kk2)))&
					*p(k1,1:4)
					p(2**(kk2-1),1:4)=parmas(ifl(kk2))/(parmas(ifl(kk1))+parmas(ifl(kk2)))&
					*p(k1,1:4)
				ENDIF
				weight=weight*w*w2
				EXIT
			ENDIF
			q=q2
			pboo(1:4)=p(k1,1:4)
		ENDDO
	ENDIF
ENDIF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IF(flagreturn.EQ.1)THEN
	IF(weight.EQ.0)THEN
		iweight=0
		IF(idir.EQ.0)iwarning(22)=iwarning(22)+1
		IF(idir.EQ.1)iwonders(22)=iwonders(22)+1
		RETURN
	ENDIF
	DO i=3,n
		j=2**(i-1)
		pt=DSQRT(p(j,1)**2+p(j,2)**2)
		pq=DSQRT(pt**2+p(j,3)**2)
		zq(i,1)= DCMPLX(p(j,4)+p(j,3),p(j,3))
		zq(i,2)= DCMPLX(p(j,4)-p(j,3),pt)
		zq(i,3)= DCMPLX(p(j,1), p(j,2))
		zq(i,4)= DCMPLX(p(j,1),-p(j,2))
		zq(i,5)= DCMPLX( parmas(ifl(i)) , pq )
		Phegas_pmom(i,1:4)=p(j,1:4)
		Phegas_pmom(i,5)=parmas(ifl(i))
	ENDDO
	IF(idir.EQ.0)CALL Phegas_testmom(iweight,xp1,xp2)
        IF((imode.EQ.1.OR.NDecayChains.NE.0).AND.(istruc.EQ.1.OR..NOT.labeqcoll).AND.emep_ISR_shower.EQ.0)THEN
                CALL cmstolab() ! in QEDPS, the boost is done in QEDPSEVNT
		CALL trans_p_to_zq(n,Phegas_pmom(1:n,1:5),zq(1:n,1:5))
	ENDIF
	RETURN
ENDIF
!  ---- when idir=1 
IF(idir.EQ.1)THEN
    i=n-1
    i0=0
    DO WHILE(i0.EQ.0)
       i=i-1
       i0=is(ig,i,1)
    ENDDO
    k0=i
ENDIF
IF(idir.EQ.0)k=1
IF(idir.EQ.1)k=k0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 10    continue
loop1:DO
  DO
    IF(k.EQ.0)RETURN
 
    i0=is(ig,k,1)
    m0=is(ig,k,2)
    i1=is(ig,k,3)
    m1=is(ig,k,4)
    i2=is(ig,k,5)
    m2=is(ig,k,6)
    i3=is(ig,k,7)
    m3=is(ig,k,8)
    IF(i0.EQ.0)THEN  !goto1000
	   EXIT
    ELSE      
       IF(i3.EQ.0.AND.ForwardQ(i1))THEN
           IF(idir.EQ.0)THEN
             k=k+1
			 IF(i2.GT.2)THEN
				p(i2,1:4)=p(i0,1:4)
				k1=i1
				IF(MOD(i1/2,2).NE.0)THEN
					p(i2,1:4)=p(2,1:4)+p(i2,1:4)
					k1=k1-2
				ENDIF
				IF(k1.NE.0)THEN
					CALL Helac_bin(n,k1,yykk)
					DO k2=3,n
						IF(yykk(n+1-k2).NE.0)p(i2,1:4)=p(i2,1:4)-p(2**(k2-1),1:4)
					ENDDO
				ENDIF
			 ENDIF
			 IF(i1.GT.2)p(i1,1:4)=p(i0,1:4)-p(i2,1:4)
           ENDIF
           IF(idir.EQ.1)THEN
             k=k-1
             p(i0,1:4)=p(i2,1:4)
			 k1=i1
			 IF(MOD(i1/2,2).NE.0)THEN
				p(i0,1:4)=p(i0,1:4)-p(2,1:4)
				k1=k1-2
			 ENDIF
			 IF(k1.NE.0)THEN
				CALL Helac_bin(n,k1,yykk)
				DO k2=3,n
					IF(yykk(n+1-k2).NE.0)p(i0,1:4)=p(i0,1:4)+p(2**(k2-1),1:4)
				ENDDO
			 ENDIF
           ENDIF
!       goto10
       ELSEIF(i3.EQ.0.AND.ForwardQ(i2))THEN
           IF(idir.EQ.0)THEN
             k=k+1
			 IF(i1.GT.0)THEN
				p(i1,1:4)=p(i0,1:4)
				k1=i2
				IF(MOD(i2/2,2).NE.0)THEN
					p(i1,1:4)=p(i1,1:4)+p(2,1:4)
					k1=k1-2
				ENDIF
				IF(k1.NE.0)THEN
					CALL Helac_bin(n,k1,yykk)
					DO k2=3,n
						IF(yykk(n+1-k2).NE.0)p(i1,1:4)=p(i1,1:4)-p(2**(k2-1),1:4)
					ENDDO
				ENDIF
			 ENDIF
			 IF(i2.GT.2)p(i2,1:4)=p(i0,1:4)-p(i1,1:4)
           ENDIF
           IF(idir.EQ.1)THEN
             k=k-1
             p(i0,1:4)=p(i1,1:4)
			 k1=i2
			 IF(MOD(i2/2,2).NE.0)THEN
				p(i0,1:4)=p(i0,1:4)-p(2,1:4)
				k1=k1-2
			 ENDIF
			 IF(k1.NE.0)THEN
				CALL Helac_bin(n,k1,yykk)
				DO k2=3,n
					IF(yykk(n+1-k2).NE.0)p(i0,1:4)=p(i0,1:4)+p(2**(k2-1),1:4)
				ENDDO
			 ENDIF
           ENDIF
!       goto10
       ELSE
	       EXIT
       ENDIF
    ENDIF
  ENDDO

!  IF(i0.EQ.0)EXIT

! 11    continue
  loop2:DO
    flagif1=0
    if1:IF(i0.NE.0)THEN       ! i0=0 we end the loop2 and loop1
      CALL Helac_bin(n,i1,y1)
	  IF(i2.EQ.0)EXIT
      CALL Helac_bin(n,i2,y2)
	  y3(1:n)=0
      IF(i3.NE.0)CALL Helac_bin(n,i3,y3)

      j1=0
      j2=0
      j3=0
      mm10=0
	  mm20=0
	  mm30=0
       
! LOOK IF IN THE SUBSEQUENT DECAY THERE IS A 
! SINGLE PARTICLE AND TAKE ITS FLAVOUR
!    --------
!       |
!       |
! 2  ---O----
!
      IF(y1(n-1).EQ.1.OR.ContainHad(y1(1:n)))THEN
         j1=1
         CALL Phegas_look(is(ig,1:8,1:8),k,i1,mm10)
      ENDIF
      IF(y2(n-1).EQ.1.OR.ContainHad(y2(1:n)))THEN
         j2=1
         CALL Phegas_look(is(ig,1:8,1:8),k,i2,mm20)
      ENDIF
      IF(i3.NE.0)THEN
         IF(y3(n-1).EQ.1.OR.ContainHad(y3(1:n)))THEN
            j3=1
            CALL Phegas_look(is(ig,1:8,1:8),k,i3,mm30)
         ENDIF
      ENDIF
      if2:IF(j1.EQ.0.AND.j2.EQ.0.AND.j3.EQ.0)THEN
        
!  ---- call s1-s2 phase space 

	     IF(idir.EQ.1)THEN
            p(i0,1:4)=p(i1,1:4)+p(i2,1:4)
            IF(i3.NE.0)p(i0,1:4)=p(i0,1:4)+p(i3,1:4)
         ENDIF

         q=scalar_product(p(i0,1:4),p(i0,1:4))

         q1min=cmas(i1)
         q2min=cmas(i2)
         q3min=cmas(i3)
       
         IF(q-(q1min+q2min+q3min)**2.LE.0)THEN
! ---------------------------------------------------------------
          IF(idebug.EQ.1)THEN
	           WRITE(nunit2,*)'----------------------------------------------------'
               WRITE(nunit2,*)'warning: 1',q,q1min,q2min,q3min,idir,ig,k,i0,i1,i2,i3
          ENDIF
! ---------------------------------------------------------------
          IF(idir.EQ.0)iwarning(1)=iwarning(1)+1
          IF(idir.EQ.1)iwonders(1)=iwonders(1)+1
          iweight=0
          RETURN
        ENDIF
! ---------------------------------------------------------------
        q=DSQRT(q)

! SYMMETRIZATION

        lsym=(q.LT.parmas(m1)+parmas(m2)+parmas(m3)&
        .AND.parwid(m1).GT.0.AND.parwid(m2)+parwid(m3).GT.0)
!      lsym=.false.
        IF(lsym)THEN
          weight_s=0
          cha=2
          IF(i3.GT.0)cha=6
          idir0=idir
          icase0=-1
          CALL init_choice(i1,m1,i2,m2,i3,m3)
          IF(idir.EQ.0)THEN
              icase=INT(Helac_rnmy(0)*cha+1)
          ELSE
              icase=1
          ENDIF
          CALL symm_choice(icase,i1,m1,i2,m2,i3,m3)
        ENDIF
! 12    continue
        DO     
          q1min=cmas(i1)
          q2min=cmas(i2)
          q3min=cmas(i3)
       
! ---- generate q1
       
          q1max=q-q2min-q3min
          IF(Helac_level(n,i1).EQ.1)q1max=q1min
          IF(q1max.GT.q1min)THEN
            IF(idir.EQ.1)THEN
               q1=scalar_product(p(i1,1:4),p(i1,1:4))
            ENDIF
            CALL Phegas_gen(idir,q1min,q1max,parmas(m1),parwid(m1),q1,w1)
            w1=w1/(2*pi)
            q1=DSQRT(q1)
          ELSE
            q1=q1min
            w1=1
          ENDIF
       
! ---- generate q2

          q2max=q-q1-q3min
          IF(Helac_level(n,i2).EQ.1)q2max=q2min
          IF(q2max.GT.q2min)THEN
             IF(idir.EQ.1)THEN
                q2=scalar_product(p(i2,1:4),p(i2,1:4)) 
             ENDIF
			 CALL Phegas_gen(idir,q2min,q2max,parmas(m2),parwid(m2),q2,w2)
			 w2=w2/(2*pi)
			 q2=DSQRT(q2)
          ELSE
             q2=q2min
             w2=1
          ENDIF
       

! ---- generate q3

          IF(i3.NE.0)THEN
             q3max=q-q2-q1
             IF(Helac_level(n,i3).EQ.1)q3max=q3min
             IF(q3max.GT.q3min)THEN
                IF(idir.EQ.1)THEN
                   q3=scalar_product(p(i3,1:4),p(i3,1:4)) 
                ENDIF
!				IF(SoniumQ2(i3+i1).AND.Helac_level(n,i1).EQ.2)THEN
!					w3=1
!					IF(idir.EQ.0)THEN
!						CALL Phegas_genQ2(idir,q3,q1,i1,i3)
!					ENDIF
!					q3=-1  ! DSQRT(q3)
!				ELSEIF(SoniumQ2(i3+i2).AND.Helac_level(n,i2).EQ.2)THEN
!					w3=1
!					IF(idir.EQ.0)THEN
!						CALL Phegas_genQ2(idir,q3,q2,i2,i3)
!					ENDIF
!					q3=-1 !DSQRT(q3)
!				ELSEIF((SoniumQ3(i1+i2+i3).AND.Helac_level(n,i3).EQ.2).OR.SoniumQ4(i1+i2,i3))THEN
!					w3=1
!					IF(idir.EQ.0)THEN
!						CALL Phegas_genQ3(idir,q3,q1,q2,i1,i2,i3)
!					ENDIF
!					q3=-1 !DSQRT(q3)
!				ELSE
!					IF(SoniumQ4(i1+i2,i3))THEN
!					ELSE
				CALL Phegas_gen(idir,q3min,q3max,parmas(m3),parwid(m3),q3,w3)
				w3=w3/(2*pi)
				q3=DSQRT(q3)
!					ENDIF
!				ENDIF
             ELSE
                q3=q3min
                w3=1
             ENDIF
          ELSE
             w3=1
          ENDIF
       
! ---- generate p(i1),p(i2),p(i3)

          IF(i3.EQ.0)THEN
! ---- generate p(i1),p(i2) for Quarkonium
!			 IF(SoniumQ(i0))THEN
!				IF(idir.EQ.0)THEN
!					p(i1,1:4)=parmas(ifl(i1))/(parmas(ifl(i1))+parmas(ifl(i2)))*p(i0,1:4)
!					p(i2,1:4)=parmas(ifl(i2))/(parmas(ifl(i1))+parmas(ifl(i2)))*p(i0,1:4)
!				ENDIF
!				w=1
!			 ELSE
			 CALL Phegas_twobody(idir,q,q1,q2,p(i1,1:4),p(i2,1:4),w)
!			 ENDIF
          ELSEIF(i3.NE.0)THEN
!			 IF(SoniumQ(i1+i3))THEN
!			 IF(idir.EQ.1)THEN
!				p(i1+i3,1:4)=p(i1,1:4)+p(i3,1:4)
!				q23= scalar_product(p(i1+i3,1:4),p(i1+i3,1:4)) 
!				IF(q23-(q1+q3)**2.LT.0)THEN
! ---------------------------------------------------------------
!					IF(idebug.EQ.1)THEN
!						WRITE(nunit2,*)'----------------------------------------------------'
!						WRITE(nunit2,*)'warning: 2',q23,q1,q3,k,i0,i2,i1,i3
!					ENDIF
! ---------------------------------------------------------------
!					IF(idir.EQ.0)iwarning(2)=iwarning(2)+1
!					IF(idir.EQ.1)iwonders(2)=iwonders(2)+1
!					iweight=0
!					RETURN
!				ENDIF
! ---------------------------------------------------------------
!			ENDIF
!			CALL Phegas_threebody(idir,q,q2,q1,q3,q23,p(i2,1:4),p(i1,1:4),p(i3,1:4),w)
!			ELSEIF(SoniumQ(i1+i2))THEN
!				IF(idir.EQ.1)THEN
!					p(i1+i2,1:4)=p(i1,1:4)+p(i2,1:4)
!					q23= scalar_product(p(i1+i2,1:4),p(i1+i2,1:4)) 
!					IF(q23-(q1+q2)**2.LT.0)THEN
! ---------------------------------------------------------------
!						IF(idebug.EQ.1)THEN
!							WRITE(nunit2,*)'----------------------------------------------------'
!							WRITE(nunit2,*)'warning: 2',q23,q1,q2,k,i0,i3,i1,i2
!						ENDIF
! ---------------------------------------------------------------
!						IF(idir.EQ.0)iwarning(2)=iwarning(2)+1
!						IF(idir.EQ.1)iwonders(2)=iwonders(2)+1
!						iweight=0
!						RETURN
!					ENDIF
! ---------------------------------------------------------------
!				ENDIF
!				CALL Phegas_threebody(idir,q,q3,q1,q2,q23,p(i3,1:4),p(i1,1:4),p(i2,1:4),w)
!			ELSE
			 IF(idir.EQ.1)THEN
				p(i2+i3,1:4)=p(i2,1:4)+p(i3,1:4)
				q23= scalar_product(p(i2+i3,1:4),p(i2+i3,1:4)) 
				IF(q23-(q2+q3)**2.LT.0)THEN
! ---------------------------------------------------------------
					IF(idebug.EQ.1)THEN
						WRITE(nunit2,*)'----------------------------------------------------'
						WRITE(nunit2,*)'warning: 2',q23,q2,q3,k,i0,i1,i2,i3
					ENDIF
! ---------------------------------------------------------------
					IF(idir.EQ.0)iwarning(2)=iwarning(2)+1
					IF(idir.EQ.1)iwonders(2)=iwonders(2)+1
					iweight=0
					RETURN
				ENDIF
! ---------------------------------------------------------------
			ENDIF
			CALL Phegas_threebody(idir,q,q1,q2,q3,q23,p(i1,1:4),p(i2,1:4),p(i3,1:4),w)
!			ENDIF
		 ENDIF
       
! ---------------------------------------------------------------

          IF(w*w1*w2*w3.EQ.0)THEN
             iweight=0
             IF(idir.EQ.0)iwarning(23)=iwarning(23)+1
             IF(idir.EQ.1)iwonders(23)=iwonders(23)+1
             RETURN
          ENDIF
       
! ---- boost to the over all frame

          IF(idir.EQ.0)THEN
             pboo(1:4)=p(i0,1:4)
             pq=DSQRT( scalar_product(pboo,pboo) )
             CALL boostl(pq,pboo,p(i1,1:4))
             CALL boostl(pq,pboo,p(i2,1:4))
             IF(i3.NE.0)CALL boostl(pq,pboo,p(i3,1:4))
          ENDIF

! SYMMETRIZATION

          IF(lsym)THEN
             weight_s=weight_s+(w1*w2*w3*w)**(-1)/cha
             IF(idir.EQ.0)THEN
               icase0=icase
               icase=0
             ENDIF
             icase=icase+1
             IF(icase.EQ.icase0)icase=icase+1
		     flag12=0
             IF(icase.LE.cha)THEN  ! goto14
                CALL symm_choice(icase,i1,m1,i2,m2,i3,m3)      
                idir=1
!        goto12
             ELSE
		        flag12=1
	         ENDIF
		     IF(flag12.EQ.1)THEN
                idir=idir0    ! 14
                CALL fini_choice(i1,m1,i2,m2,i3,m3)
                weight_s=weight_s**(-1)
		        EXIT
		     ENDIF
          ELSE
             weight_s=w1*w2*w3*w
	         EXIT
          ENDIF
       ENDDO
       weight=weight*weight_s
       EXIT   ! exit loop2
! forward !      
!       IF(idir.EQ.1)THEN
!         k=k-1
!     goto10
!     ELSE
!        k=k+1
!     ENDIF
      
!       k=k+1
!       goto10      
	         
!------------------------------------------------------------------------------
!------------------------------------------------------------------------------       
     ELSE if2   ! there is at least one non-vanishing element of j1,j2,j3 
!------------------------------------------------------------------------------       
!------------------------------------------------------------------------------
! ---- call t1-s2 phase space  
       w_costh=1
       w_phi=1
       w_q1=1
       w_q2=1
      
! ----------------------------------------------------------
! CHECK IF THE II1 PARTICLE IS JUST 2 THEN NO GENERATION
! IS NEDDED
!            |
!        i0  | /<---- ii2 --*
!            |/             |
!    2 ------*----    <-----*
!
	   IF(j1.EQ.1)THEN
			IF(ForwardQ(i1))THEN
				ii2=i2+i3
				IF(idir.EQ.0)THEN
					IF(ii2.GT.2)THEN
						p(ii2,1:4)=p(i0,1:4)
						IF(y1(n-1).NE.0)THEN
							p(ii2,1:4)=p(ii2,1:4)+p(2,1:4)
						ENDIF
						IF(i1.GT.2)THEN
							DO k2=3,n
								IF(y1(n+1-k2).NE.0)p(ii2,1:4)=p(ii2,1:4)-p(2**(k2-1),1:4)
							ENDDO
						ENDIF
					ENDIF
					IF(i1.GT.2)p(i1,1:4)=p(i0,1:4)-p(ii2,1:4)
				ENDIF
				IF(idir.EQ.1)THEN
					IF(i2.EQ.2)THEN
						p(ii2,1:4)=p(i3,1:4)-p(i2,1:4)
					ELSEIF(i3.EQ.2)THEN
						p(ii2,1:4)=p(i2,1:4)-p(i3,1:4)
					ELSE
						p(ii2,1:4)=p(i2,1:4)+p(i3,1:4)
					ENDIF
					p(i0,1:4)=p(ii2,1:4)
					IF(y1(n-1).NE.0)THEN
						p(i0,1:4)=p(i0,1:4)-p(2,1:4)
					ENDIF
					IF(i1.GT.2)THEN
						DO k2=3,n
							IF(y1(n+1-k2).NE.0)p(i0,1:4)=p(i0,1:4)+p(2**(k2-1),1:4)
						ENDDO
					ENDIF
				ENDIF
				i1=i2
				m1=m2
				i2=i3
				m2=m3
				i0=ii2
				i3=0
				m3=0
				CYCLE  ! cycle loop2  ! goto 11    ! to generate the momentum of the rest two
			ENDIF
	   ENDIF

	   IF(j2.EQ.1)THEN
			IF(ForwardQ(i2))THEN
				ii2=i1+i3
				IF(idir.EQ.0)THEN
					IF(ii2.GT.2)THEN
						p(ii2,1:4)=p(i0,1:4)
						IF(y2(n-1).NE.0)THEN
							p(ii2,1:4)=p(ii2,1:4)+p(2,1:4)
						ENDIF
						IF(i2.GT.2)THEN
							DO k2=3,n
								IF(y2(n+1-k2).NE.0)p(ii2,1:4)=p(ii2,1:4)-p(2**(k2-1),1:4)
							ENDDO
						ENDIF
					ENDIF
					IF(i2.GT.2)p(i2,1:4)=p(i0,1:4)-p(ii2,1:4)
				ENDIF
				IF(idir.EQ.1)THEN
					IF(i1.EQ.2)THEN
						p(ii2,1:4)=p(i3,1:4)-p(i1,1:4)
					ELSEIF(i3.EQ.2)THEN
						p(ii2,1:4)=p(i1,1:4)-p(i3,1:4)
					ELSE
						p(ii2,1:4)=p(i1,1:4)+p(i3,1:4)
					ENDIF
					p(i0,1:4)=p(ii2,1:4)
					IF(y2(n-1).NE.0)THEN
						p(i0,1:4)=p(i0,1:4)-p(2,1:4)
					ENDIF
					IF(i2.NE.2)THEN
						DO k2=3,n
							IF(y2(n+1-k2).NE.0)p(i0,1:4)=p(i0,1:4)+p(2**(k2-1),1:4)
						ENDDO
					ENDIF
				ENDIF
				i1=i1
				m1=m1
				i2=i3
				m2=m3
				i0=ii2
				i3=0
				m3=0
				CYCLE  ! cycle loop2  ! goto 11    ! to generate the momentum of the rest two
			ENDIF
	   ENDIF

	   IF(j3.EQ.1)THEN
			IF(ForwardQ(i3))THEN
				ii2=i2+i1
				IF(idir.EQ.0)THEN
					IF(ii2.GT.2)THEN
						p(ii2,1:4)=p(i0,1:4)
						IF(y3(n-1).NE.0)THEN
							p(ii2,1:4)=p(ii2,1:4)+p(2,1:4)
						ENDIF
						IF(i3.NE.2)THEN
							DO k2=3,n
								IF(y3(n+1-k2).NE.0)p(ii2,1:4)=p(ii2,1:4)-p(2**(k2-1),1:4)
							ENDDO
						ENDIF
					ENDIF
					IF(i3.GT.2)p(i3,1:4)=p(i0,1:4)-p(ii2,1:4)
				ENDIF
				IF(idir.EQ.1)THEN
					IF(i1.EQ.2)THEN
						p(ii2,1:4)=p(i2,1:4)-p(i1,1:4)
					ELSEIF(i2.EQ.2)THEN
						p(ii2,1:4)=p(i1,1:4)-p(i2,1:4)
					ELSE
						p(ii2,1:4)=p(i2,1:4)+p(i1,1:4)
					ENDIF
					p(i0,1:4)=p(ii2,1:4)
					IF(y3(n-1).NE.0)THEN
						p(i0,1:4)=p(i0,1:4)-p(2,1:4)
					ENDIF
					IF(i3.NE.2)THEN
						DO k2=3,n
							IF(y3(n+1-k2).NE.0)p(i0,1:4)=p(i0,1:4)+p(2**(k2-1),1:4)
						ENDDO
					ENDIF
				ENDIF
				i1=i1
				m1=m1
				i2=i2
				m2=m2
				i0=ii2
				i3=0
				m3=0
				CYCLE  ! cycle loop2  ! goto 11    ! to generate the momentum of the rest two
			ENDIF
	   ENDIF
! ----------------------------------------------------------
! CHECK IF THERE ARE >3 PARTICLES IN THE VERTEX
! AND COMBINE THEM
       itogo=0
       iq2max=0
       IF(j1.EQ.1)THEN
          ii1=i1
          mm1=m1
		  yykk(1:n)=y1(1:n)
          ii2=i2+i3
          mm2=0
		  kk1=SubtractHad(ii1)
		  kk2=SubtractHad(ii2)
          IF(i3.EQ.0)mm2=mm20
		  mm10=mm10
		  y2(1:n)=y2(1:n)+y3(1:n)
		  y1(1:n)=y1(1:n)+y2(1:n)
          IF(idir.EQ.1)THEN
			p(ii2,1:4)=p(i2,1:4)+p(i3,1:4)
			IF(kk2.NE.ii2.AND.i3.NE.0)THEN
				p(kk2,1:4)=p(ii2,1:4)
				IF(y2(n-1).NE.0)p(kk2,1:4)=p(kk2,1:4)+p(2,1:4)
				DO k1=3,n
					IF(y2(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)p(kk2,1:4)=p(kk2,1:4)&
					-p(2**(k1-1),1:4)
				ENDDO
			ENDIF
		  ENDIF
          IF(Helac_level(n,kk1).EQ.1)itogo=1
          IF(Helac_level(n,kk2).EQ.1)iq2max=1
       ELSEIF(j2.EQ.1)THEN
          ii1=i2
          mm1=m2
		  yykk(1:n)=y2(1:n)
          ii2=i1+i3
          mm2=0
		  kk1=SubtractHad(ii1)
		  kk2=SubtractHad(ii2)
          IF(i3.EQ.0)mm2=mm10
		  mm10=mm20
		  y1(1:n)=y1(1:n)+y2(1:n)+y3(1:n)
		  y2(1:n)=y1(1:n)-y2(1:n)
          IF(idir.EQ.1)THEN
			p(ii2,1:4)=p(i1,1:4)+p(i3,1:4)
			IF(kk2.NE.ii2.AND.i3.NE.0)THEN
				p(kk2,1:4)=p(ii2,1:4)
				IF(y2(n-1).NE.0)p(kk2,1:4)=p(kk2,1:4)+p(2,1:4)
				DO k1=3,n
					IF(y2(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)p(kk2,1:4)=p(kk2,1:4)&
					-p(2**(k1-1),1:4)
				ENDDO
			ENDIF
		  ENDIF
          IF(Helac_level(n,kk1).EQ.1)itogo=1
          IF(Helac_level(n,kk2).EQ.1)iq2max=1
       ELSEIF(j3.EQ.1)THEN
          ii1=i3
          mm1=m3
		  yykk(1:n)=y3(1:n)
          ii2=i1+i2
          mm2=0
		  mm10=mm30
		  kk1=SubtractHad(ii1)
		  kk2=SubtractHad(ii2)
		  y2(1:n)=y1(1:n)+y2(1:n)
		  y1(1:n)=y2(1:n)+y3(1:n)
          IF(idir.EQ.1)THEN
			p(ii2,1:4)=p(i1,1:4)+p(i2,1:4)
			IF(kk2.NE.ii2)THEN
				p(kk2,1:4)=p(ii2,1:4)
				IF(y2(n-1).NE.0)p(kk2,1:4)=p(kk2,1:4)+p(2,1:4)
				DO k1=3,n
					IF(y2(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)p(kk2,1:4)=p(kk2,1:4)&
					-p(2**(k1-1),1:4)
				ENDDO
			ENDIF
		  ENDIF
          IF(Helac_level(n,kk1).EQ.1)itogo=1
          IF(Helac_level(n,kk2).EQ.1)iq2max=1
       ENDIF
       
       q1min=cmas(kk1)   ! there is no cut to the particle 2
       q2min=cmas(kk2)
       
       IF(idir.EQ.1)THEN
          p(kk1+kk2,1:4)=p(kk1,1:4)+p(kk2,1:4)
          p(i0,1:4)=p(ii1,1:4)+p(ii2,1:4)
       ENDIF
       IF(idir.EQ.0)THEN
			p(kk1+kk2,1:4)=p(i0,1:4)
			IF(y1(n-1).NE.0)p(kk1+kk2,1:4)=p(kk1+kk2,1:4)+p(2,1:4)
			DO k1=3,n
				IF(y1(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)p(kk1+kk2,1:4)=p(kk1+kk2,1:4)&
				-p(2**(k1-1),1:4)
			ENDDO
	   ENDIF    
       q=scalar_product(p(kk1+kk2,1:4),p(kk1+kk2,1:4))
       IF(q-(q1min+q2min)**2.LE.0)THEN
! ---------------------------------------------------------------
          IF(idebug.EQ.1)THEN
	          WRITE(nunit2,*)'--------------------------------------------'
              WRITE(nunit2,*)'warning: 3'
              WRITE(nunit2,*)'q,q1min,q2min,k,i0,i1,m1,i2,m2,i3,m3,ig,idir'
              WRITE(nunit2,*) q,q1min,q2min,k,i0,i1,m1,i2,m2,i3,m3,ig,idir
          ENDIF
! ---------------------------------------------------------------
          IF(idir.EQ.0)iwarning(3)=iwarning(3)+1
          IF(idir.EQ.1)iwonders(3)=iwonders(3)+1
          iweight=0
          RETURN
       ENDIF
! ---------------------------------------------------------------
       q=DSQRT(q)
       
! ---- generate q2
       
       q2max=q-q1min
       
! print*,'q2max,q2min'
! print*, q2max,q2min 

       IF(iq2max.EQ.1)q2max=q2min  ! when the q2 is the external particles
       IF(q2max.GT.q2min)THEN
           IF(idir.EQ.1)THEN
               q2=scalar_product(p(kk2,1:4),p(kk2,1:4)) 
           ENDIF
           IF(i3.EQ.0)THEN
               IF(q-parmas(mm2)-parmas(mm10).LT.0.AND.&
			   q-parmas(mm10).GT.q2min.AND.q-parmas(mm10).LT.q2max)THEN
!if(q-parmas(mm10).gt.q2min.and.q-parmas(mm10).lt.q2max)then
                  CALL Phegas_gen1(idir,q2min,q2max,parmas(mm2),parwid(mm2),&
                         parmas(mm10),parwid(mm10),q,q1min,q2,w_q2)
               ELSE
                  CALL Phegas_gen(idir,q2min,q2max,parmas(mm2),parwid(mm2),q2,w_q2)
               ENDIF
           ELSEIF(i3.NE.0)THEN
               CALL Phegas_gen(idir,q2min,q2max,r_0,-r_1,q2,w_q2)
           ENDIF
           w_q2=w_q2/(2*pi)
           q2=DSQRT(q2)
        ELSE
           q2=q2min
           w_q2=1
        ENDIF
       
        IF(idir.EQ.1)THEN
           q1=scalar_product(p(kk1,1:4),p(kk1,1:4))
!---------------------------------------------------------------------
           IF(q1.LT.q1min**2.AND.Helac_level(n,kk1).GT.1)THEN
!---------------------------------------------------------------------
                IF(idebug.EQ.1)THEN
	                WRITE(nunit2,*)'---------------------------------------------'
                    WRITE(nunit2,*)'warning 18,idir=',idir
                    WRITE(nunit2,*)'q1,q1min**2,ii1,ig,k,p(kk1,1:4)'
                    WRITE(nunit2,*) q1,q1min**2,ii1,ig,k,p(kk1,1:4)
                ENDIF
!---------------------------------------------------------------------
                IF(idir.EQ.0)iwarning(18)=iwarning(18)+1
                IF(idir.EQ.1)iwonders(18)=iwonders(18)+1
                iweight=0
	            RETURN
            ENDIF
!---------------------------------------------------------------------
        ENDIF
! ---- generate q1
        gm1=parmas(mm10)
        gg1=parwid(mm10)
       
        q1max=q-q2
        IF(itogo.EQ.1)q1max=q1min
        IF(q1max.GT.q1min)THEN
            CALL Phegas_gen(idir,q1min,q1max,gm1,gg1,q1,w_q1)
            w_q1=w_q1/(2*pi)
            q1=DSQRT(q1)
        ELSE
            q1=q1min
            w_q1=1
        ENDIF

!if(idir.eq.0)print*,ig,mm10,gg1,gm1,q1min,q1max,q,q1,q2

! ---------------------------------------------------------------
        IF(w_q1*w_q2.EQ.0)THEN
            iweight=0
            IF(idir.EQ.0)iwarning(24)=iwarning(24)+1
            IF(idir.EQ.1)iwonders(24)=iwonders(24)+1
            RETURN
        ENDIF       
! ---- boost to p(kk1+kk2,1:4) rest frame

        pboo(4)  =p(kk1+kk2,4)
        pboo(1:3)=-p(kk1+kk2,1:3)
        pq=DSQRT( scalar_product(pboo,pboo) )
		pp1(1:4)=0
		IF(yykk(n-1).NE.0)pp1(1:4)=pp1(1:4)+p(2,1:4)
		DO k1=3,n
			IF(yykk(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)pp1(1:4)=pp1(1:4)-p(2**(k1-1),1:4)
		ENDDO
        CALL boost(pq,pboo,pp1(1:4),p2)
        IF(idir.EQ.1)CALL boost(pq,pboo,p(kk1,1:4),pdummy)
        em=scalar_product(pp1(1:4),pp1(1:4)) ! parmas(ifl(2))
        e2=p2(4)     
! ---- find t1max t1min

!        t1max=tfunp(q1**2)
!        t1min=tfunm(q1**2)
! Eq.(5) in hep-ph/0007335
        t1max=Phegas_tfunp(q1)   ! The maxmal invariant mass of i1 (if i1 contains 2)
        t1min=Phegas_tfunm(q1)   ! The minimal invariant mass of i1 (if i1 contains 2)
        a1=( t1min+t1max)/2
        a2=(-t1min+t1max)/2
        IF(tm(ii1).NE.0)THEN
           IF(tm(ii1).LT.t1min)THEN
! ---------------------------------------------------------------
                IF(idebug.EQ.1)THEN
                      WRITE(nunit2,*)'----------------------------------------------'
                      WRITE(nunit2,*)'warning 4: for t1min < tm(',ii1,') in 1a'
                      WRITE(nunit2,*)i0,ii1,ii2,kk1,kk2,idir,ig,ig0,k
                      WRITE(nunit2,*)tm(ii1),t1min,t1l,t1u
                      WRITE(nunit2,*)q,q1,q2,e2,em,q1min,q1max
                ENDIF
! ---------------------------------------------------------------
                IF(idir.EQ.0)iwarning(4)=iwarning(4)+1
                IF(idir.EQ.1)iwonders(4)=iwonders(4)+1
                iweight=-1
                RETURN
             ENDIF
! ---------------------------------------------------------------
             t1max=MIN(tm(ii1),t1max)
          ENDIF

          t1l=-t1max+parmas(mm1)**2
          t1u=-t1min+parmas(mm1)**2
          IF(t1l.GT.t1u)THEN
! ---------------------------------------------------------------
             IF(idebug.EQ.1)THEN
                WRITE(nunit2,*)'----------------------------------------------'
                WRITE(nunit2,*)'warning 7: for t1l > t1u'
                WRITE(nunit2,*)i0,ii1,ii2,idir,ig,ig0,k
                WRITE(nunit2,*)t1l,t1u,t1max,t1min,q1,mm1
                WRITE(nunit2,*)q,q1,q2,e2,em
             ENDIF
! ---------------------------------------------------------------
             IF(idir.EQ.0)iwarning(7)=iwarning(7)+1
             IF(idir.EQ.1)iwonders(7)=iwonders(7)+1
             iweight=0
             RETURN
          ENDIF
		  a1=a1-parmas(mm1)**2
		  IF(t1u.LT.0)THEN
			t1l=-t1l
			t1u=-t1u
			a1=-a1
			a2=-a2
		  ENDIF

          IF(t1l/t1u.LT.-10d0*EPSILON(t1l))THEN
! ---------------------------------------------------------------
             IF(idebug.EQ.1)THEN
                  WRITE(nunit2,*)'----------------------------------------------'
                  WRITE(nunit2,*)'warning 9: for t1l < -10*epsilon(t1l)'
                  WRITE(nunit2,*)i0,ii1,ii2,idir,ig,ig0,k
                  WRITE(nunit2,*)t1l,t1u,t1max,t1min,q1,mm1,EPSILON(t1l)
                  WRITE(nunit2,*)q,q1,q2,e2,em
             ENDIF
! ---------------------------------------------------------------
             IF(idir.EQ.0)iwarning(9)=iwarning(9)+1
             IF(idir.EQ.1)iwonders(9)=iwonders(9)+1
             iweight=0
             RETURN
           ENDIF
		   IF(t1l.LT.0.AND.t1l/t1u.GE.-10d0*EPSILON(t1l))t1l=0d0
!           IF(t1l.LE.0)THEN
! ---------------------------------------------------------------
!              IF(idebug.EQ.1)THEN
!                 WRITE(nunit2,*)'----------------------------------------------'
!                 WRITE(nunit2,*)'warning 8: for t1l < 0'
!                 WRITE(nunit2,*)i0,ii1,ii2,idir,ig,ig0,k
!                 WRITE(nunit2,*)t1l,t1u,t1max,t1min,q1,mm1
!                 WRITE(nunit2,*)q,q1,q2,e2,em
!              ENDIF
! ---------------------------------------------------------------
!              IF(idir.EQ.0)iwarning(8)=iwarning(8)+1
!              IF(idir.EQ.1)iwonders(8)=iwonders(8)+1
!              iweight=0
!              RETURN
!           ENDIF
       
! ---- generate costh
     
           IF(idir.EQ.1)THEN
              cth=cosij(p2(1:3),pdummy(1:3))  ! in Kinetic_Func
           ENDIF

!      call givemea1a2(q1,a1,a2)
           CALL Phegas_gen_x(idir,-a1,-a2,t1l,t1u,cth,w_costh)    
		   GLOBALINIT_gen_x=1  
! -------------------------------------------------
           IF(w_costh.LE.0.OR.(1-cth*cth).LT.0)THEN
!       print*,idir,ig,ig0,i0,m0,i1,m1,i2,m2,i3,m3
!       print*,a1,a2,q,q1,q2
               IF(idir.EQ.0)iwarning(5)=iwarning(5)+1
               IF(idir.EQ.1)iwonders(5)=iwonders(5)+1
               iweight=0
               RETURN
           ENDIF
           t1=a1+a2*cth
      
! ---- generate phi
!		   flagphi=0
!		   IF(i3.EQ.0.AND.kk2.LT.ii2)THEN
!			flagphi=1
!			CALL boost(pq,pboo,-p(i0,1:4)+p(kk1+kk2,1:4)-pp1(1:4),pffp3(1:4))
!			CALL Phegas_tfunpm2(q,q1,q2,p2(1:4),pffp3(1:4),cth,t2max,t2min)
!			IF(tm(ii2).NE.0)THEN
!				IF(tm(ii2).LT.t2min)THEN
! ---------------------------------------------------------------
!					IF(idebug.EQ.1)THEN
!						WRITE(nunit2,*)'----------------------------------------------'
!						WRITE(nunit2,*)'warning 4: for t2min < tm(',ii2,') in 1a'
!						WRITE(nunit2,*)i0,ii1,ii2,kk1,kk2,idir,ig,ig0,k
!						WRITE(nunit2,*)tm(ii2),t2min
!					ENDIF
! ---------------------------------------------------------------
!					IF(idir.EQ.0)iwarning(4)=iwarning(4)+1
!					IF(idir.EQ.1)iwonders(4)=iwonders(4)+1
!					iweight=-1
!					RETURN
!				ENDIF
! ---------------------------------------------------------------
!				t2max=MIN(tm(ii2),t2max)
!			ENDIF
!			IF(t2min.GT.t2max)THEN
! ---------------------------------------------------------------
!				IF(idebug.EQ.1)THEN
!					WRITE(nunit2,*)'----------------------------------------------'
!					WRITE(nunit2,*)'warning 7: for t2min > t2max'
!					WRITE(nunit2,*)i0,ii1,ii2,idir,ig,ig0,k
!					WRITE(nunit2,*)t2min,t2max
!					WRITE(nunit2,*)q,q1,q2
!				ENDIF
! ---------------------------------------------------------------
!				IF(idir.EQ.0)iwarning(7)=iwarning(7)+1
!				IF(idir.EQ.1)iwonders(7)=iwonders(7)+1
!				iweight=0
!				RETURN
!			ENDIF
!			IF(t2max*t2min.LT.0)THEN
!				flagphi=0
!				IF(idir.EQ.0)THEN
!					r=Helac_rnmy(0)
!					phi=2*pi*r
!				ENDIF
!				w_phi=2*pi
!			ELSE
!				IF(idir.EQ.1)THEN
!					cth0=p2(3)/DSQRT(p2(1)**2+p2(2)**2+p2(3)**2)
!					IF(p2(1)**2+p2(2)**2.EQ.0)THEN
!						IF(p2(3).GT.0)cth0=1
!						IF(p2(3).LT.0)cth0=-1
!					ENDIF
!					phi0=ph4(p2(1),p2(2),p2(3))
!				! rotate p2 to z -axis
!					CALL rotate_inv(cth0,DCOS(phi0),DSIN(phi0),pffp3(1:3))
!					CALL rotate_inv(cth0,DCOS(phi0),DSIN(phi0),p(kk1,1:3))
				! rotate p3 to x,z-plane
!					phi0=ph4(pffp3(1),pffp3(2),pffp3(3))
!					CALL rotate_inv(1d0,DCOS(phi0),DSIN(phi0),p(kk1,1:3))
!					phi=ph4(p(kk1,1),p(kk1,2),p(kk1,3))
!				ENDIF
!				CALL Phegas_gen_x2(idir,parmas(m2),t2min,t2max,phi,w_phi)
!			IF(w_phi.LE.0)THEN
!				IF(idir.EQ.0)iwarning(5)=iwarning(5)+1
!				IF(idir.EQ.1)iwonders(5)=iwonders(5)+1
!				iweight=0
!				RETURN
!			ENDIF
!			ENDIF
!		   ELSE
			IF(idir.EQ.0)THEN
               r=Helac_rnmy(0)
               phi=2*pi*r
			ENDIF
			w_phi=2*pi
!		   ENDIF
 ! generate the momentums of kk1 and kk2 in   p(kk1+kk2,1:4) rest frame     
           IF(idir.EQ.0)THEN
       
               sth=DSQRT(1-cth*cth)
               p(kk1,4)=q/2+q1**2/q/2-q2**2/q/2
               p(kk2,4)= q-p(kk1,4)
!      pm=p(ii1-2,4)*dsqrt(1-(q1/p(ii1-2,4))**2)
!rewind(100)
!write(100,*)p(ii1-2,4),pm,q1,p(ii1-2,4)**2-pm**2-q1**2
!      if(q2.gt.q1)pm=p(ii2,4)*dsqrt(1-(q2/p(ii2,4))**2)
!write(100,*)p(ii2,4),pm,q2,p(ii2,4)**2-pm**2-q2**2
               pm=Phegas_clambda(q,q1,q2)/q/2
       
               p(kk1,3)=pm*cth
               p(kk1,2)=pm*sth*DSIN(phi)
               p(kk1,1)=pm*sth*DCOS(phi)
               p(kk2,1:3)=-p(kk1,1:3)

! ---- rotate  			   

               cth0=p2(3)/DSQRT(p2(1)**2+p2(2)**2+p2(3)**2)
               IF(p2(1)**2+p2(2)**2.EQ.0)THEN
                   IF(p2(3).GT.0)cth0=1
                   IF(p2(3).LT.0)cth0=-1
               ENDIF
               phi0=ph4(p2(1),p2(2),p2(3))

!			   IF(i3.EQ.0.AND.kk2.LT.ii2.AND.flagphi.EQ.1)THEN
!					CALL rotate_inv(cth0,DCOS(phi0),DSIN(phi0),pffp3(1:3))
!					phi02=ph4(pffp3(1),pffp3(2),pffp3(3))
!					CALL Phegas_rotate(1d0,DCOS(phi02),DSIN(phi02),p(kk1,1:3))
!					CALL Phegas_rotate(1d0,DCOS(phi02),DSIN(phi02),p(kk2,1:3))
!			   ENDIF
       
               CALL Phegas_rotate(cth0,DCOS(phi0),DSIN(phi0),p(kk1,1:3))
               CALL Phegas_rotate(cth0,DCOS(phi0),DSIN(phi0),p(kk2,1:3))
!call rotate(cth0,dcos(phi0),dsin(phi0),p2(1:3))
       

! ---- boost to the over all frame

               pboo(1:4)=p(kk1+kk2,1:4)
               pq=DSQRT( scalar_product(pboo,pboo) )
               CALL boostl(pq,pboo,p(kk1,1:4))
               CALL boostl(pq,pboo,p(kk2,1:4))

       
               t1new=scalar_product(p(kk1,1:4),p(kk1,1:4))&
                   +scalar_product(pp1(1:4),pp1(1:4))&
                   -2*scalar_product(p(kk1,1:4),pp1(1:4))
!       data it1new/0/
               IF(ABS((t1new-t1)/t1).GT.1.AND.it1new.LT.10)THEN
                    it1new=it1new+1
                    IF(idebug.EQ.1)THEN
                        WRITE(nunit2,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
                        WRITE(nunit2,*)'abs((t1new-t1)/t1),t1,t1new'
                        WRITE(nunit2,*) ABS((t1new-t1)/t1),t1,t1new
                        WRITE(nunit2,*)'t1min,t1max,q2max,q2min,q1max,q1min'
                        WRITE(nunit2,*) t1min,t1max,q2max,q2min,q1max,q1min 
                        WRITE(nunit2,*)'q,q1,q2'
                        WRITE(nunit2,*) q,q1,q2 
                        WRITE(nunit2,*)'ig,k,i0,m0,i1,m1,i2,m2,i3,m3'
                        WRITE(nunit2,*) ig,k,i0,m0,i1,m1,i2,m2,i3,m3 
                        WRITE(nunit2,*)'try an other set of cuts'
                        WRITE(nunit2,*)'OR change to QP mode'
                        WRITE(nunit2,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
                     ENDIF
!       stop
               ENDIF

           ENDIF
       
! ---- calculate total weight
       
!      w_lam=clambda(q**2,q1**2,q2**2)
           w_lam=Phegas_clambda(q,q1,q2)
           weight=weight*w_q1*w_q2*w_costh*w_phi*w_lam/(32*pi**2)/q**2
     
!     if(weight.eq.0)then
!     print*,'weight,w_q1,w_q2,w_costh,w_phi,q,k'
!     print*, weight,w_q1,w_q2,w_costh,w_phi,q,k
!     endif
           IF(idir.EQ.0)THEN
				 p(ii1,1:4)=p(kk1,1:4)
				 IF(yykk(n-1).NE.0)p(ii1,1:4)=p(ii1,1:4)-p(2,1:4)
				 DO k1=3,n
					IF(yykk(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)p(ii1,1:4)=p(ii1,1:4)+&
					p(2**(k1-1),1:4)
				 ENDDO
!                 IF(j1.EQ.1)THEN
!                     p(i1,1:4)=p(ii1-2,1:4)-p(2,1:4)
!                 ELSEIF(j2.EQ.1)THEN
!                     p(i2,1:4)=p(ii1-2,1:4)-p(2,1:4)
!                 ELSEIF(j3.EQ.1)THEN
!                     p(i3,1:4)=p(ii1-2,1:4)-p(2,1:4)
!                 ENDIF
           ENDIF


		   IF(ii2.NE.kk2)THEN
				IF(idir.EQ.0)THEN
					y1(1:n)=y1(1:n)-yykk(1:n)
					p(ii2,1:4)=p(kk2,1:4)
					IF(y1(n-1).NE.0)p(ii2,1:4)=p(ii2,1:4)-p(2,1:4)
					DO k1=3,n
						IF(y1(n+1-k1).NE.0.AND.Quarkonium3(k1).NE.0)p(ii2,1:4)=p(ii2,1:4)&
						+p(2**(k1-1),1:4)
					ENDDO
				ENDIF
		   ENDIF
       
! ---------------------------------------------------------------
! Go back (11) to generate the cluster of two if i3=/=0
! VVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVVV
           IF(i3.NE.0)THEN
                 i0=ii2
                 IF(j1.EQ.1)THEN
                      i1=i2
                      m1=m2
                      i2=i3
                      m2=m3
                 ELSEIF(j2.EQ.1)THEN
                      i1=i1
                      i2=i3
                      m2=m3
                 ELSEIF(j3.EQ.1)THEN
                      i1=i1
                      i2=i2
                 ENDIF
                 i3=0
                 m3=0
                 CYCLE      !goto11
           ENDIF
! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
! Go back (11) to generate the cluster of two if i3=/=0
! ---------------------------------------------------------------
           EXIT
 !          IF(idir.EQ.1)THEN
 !                k=k-1
 !                goto10
!		   ELSE
!		         k=k+1
 !          ENDIF
!       k=k+1
 !      goto10      
    ENDIF  if2
  
  ELSE if1
   flagif1=1
   EXIT 
  ENDIF if1
  ENDDO loop2
! forward
  IF(flagif1.EQ.1)THEN
    EXIT
  ELSE
    IF(idir.EQ.1)THEN
       k=k-1
!     goto10
    ELSE
       k=k+1
    ENDIF
   ENDIF      
ENDDO loop1
! 1000  continue

IF(weight.EQ.0)THEN
   iweight=0
   IF(idir.EQ.0)iwarning(22)=iwarning(22)+1
   IF(idir.EQ.1)iwonders(22)=iwonders(22)+1
   RETURN
ENDIF

! ---- define zq's 
! in calculation of matrix elements, one only need the external zq.
DO i=3,n
   j=2**(i-1)
   pt=DSQRT(p(j,1)**2+p(j,2)**2)
   pq=DSQRT(pt**2+p(j,3)**2)
   zq(i,1)= DCMPLX(p(j,4)+p(j,3),p(j,3))
   zq(i,2)= DCMPLX(p(j,4)-p(j,3),pt)
   zq(i,3)= DCMPLX(p(j,1), p(j,2))
   zq(i,4)= DCMPLX(p(j,1),-p(j,2))
   zq(i,5)= DCMPLX( parmas(ifl(i)) , pq )
   Phegas_pmom(i,1:4)=p(j,1:4)
   Phegas_pmom(i,5)=parmas(ifl(i))
!   PRINT *,Phegas_pmom(i,1:4)
!   PRINT *,scalar_product(Phegas_pmom(i,1:4),Phegas_pmom(i,1:4))
ENDDO

IF(idir.EQ.0)CALL Phegas_testmom(iweight,xp1,xp2)

IF((imode.EQ.1.OR.NDecayChains.NE.0).AND.(istruc.EQ.1.OR..NOT.labeqcoll).AND.emep_ISR_shower.EQ.0)THEN
	CALL cmstolab() ! in QEDPS, the boost is done in QEDPSEVNT
	CALL trans_p_to_zq(n,Phegas_pmom(1:n,1:5),zq(1:n,1:5))
ENDIF

END SUBROUTINE Phegas

!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
!&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
       
FUNCTION Phegas_tfunm(xss)
!       include 'declare.h'
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::xss
!       common/phegas1/q,q2,e2,em
REAL(KIND=DBL)::Phegas_tfunm
REAL(KIND=DBL)::xx1gg,xx2gg,p2gg

p2gg=ABS(e2*DSQRT(1-em/e2**2))

xx1gg=(q**2+xss**2-q2**2)/q/2
xx2gg=Phegas_clambda(q,xss,q2)/q/2
!      sqrtarg=1-(2*q/xx1)**2*x
!      if(sqrtarg.lt.0d0)sqrtarg=0d0
!      xx2=xx1*dsqrt(sqrtarg)
Phegas_tfunm=em+xss**2-2*e2*xx1gg-2*p2gg*xx2gg
END FUNCTION Phegas_tfunm
       
! This rotates a vector an angle theta around Y and ph around Z
! so that a (0,0,1) -> (sin(th)cos(ph),sin(th)sin(ph),cos(th))

SUBROUTINE Phegas_rotate(cth,cph,sph,a)
!      include 'declare.h'
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::cth,cph,sph
REAL(KIND=DBL),DIMENSION(3),INTENT(INOUT)::a
REAL(KIND=DBL),DIMENSION(3)::b
REAL(KIND=DBL)::sth
sth=DSQRT(1-cth*cth)
b(1)= a(1)*cth*cph-a(2)*sph+a(3)*sth*cph
b(2)= a(1)*cth*sph+a(2)*cph+a(3)*sth*sph
b(3)=-a(1)*sth             +a(3)*cth
a(1:3)=b(1:3)
END SUBROUTINE Phegas_rotate

! This rotates a vector an angle -ph around Z and -th around Y
! so that a (sin(th)cos(ph),sin(th)sin(ph),cos(th)) -> (0,0,1)

SUBROUTINE rotate_inv(cth,cph,sph,a)
!       include 'declare.h'
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(3),INTENT(INOUT)::a
REAL(KIND=DBL),DIMENSION(3)::b
REAL(KIND=DBL),INTENT(IN)::cth,cph,sph
REAL(KIND=DBL)::sth
sth=DSQRT(1-cth*cth)
b(1)= a(1)*cth*cph+a(2)*cth*sph-a(3)*sth
b(2)=-a(1)*sph    +a(2)*cph    
b(3)= a(1)*sth*cph+a(2)*sth*sph+a(3)*cth
a(1:3)=b(1:3)
END SUBROUTINE rotate_inv
       
       
SUBROUTINE Phegas_twobody(idir,q,q1,q2,p1,p2,w)
!       include 'declare.h'
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::q,q1,q2
INTEGER::init=0
REAL(KIND=DBL)::pi,r_0,r_1,r1,r2,e1,e2,p,cth,sth,ph
SAVE pi
REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::p1,p2
REAL(KIND=DBL),INTENT(OUT)::w
r_0=0
r_1=1
      
IF(init.EQ.0)THEN
   CALL Helac_mypi(pi)
   init=1
ENDIF
       
w=0
IF(q.LE.q1+q2)RETURN
       
e1=(q**2+q1**2-q2**2)/q/2
e2=q-e1
p=Phegas_clambda(q,q1,q2)/q/2

w=p/q/(4*pi)
       
IF(idir.EQ.1)RETURN
       
p1(4)=e1
p2(4)=e2
       
r1=Helac_rnmy(0)
r2=Helac_rnmy(0)

cth=2*r1-1
sth=DSQRT(1-cth*cth)
ph=2*pi*r2
       
p1(3)=p*cth
p1(2)=p*sth*DSIN(ph)
p1(1)=p*sth*DCOS(ph)
       
p2(1:3)=-p1(1:3)

END SUBROUTINE Phegas_twobody

SUBROUTINE Phegas_twobody_had(idir,q,q1,q2,ptcut,p1,p2,w)
!       include 'declare.h'
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::q,q1,q2,ptcut
INTEGER::init=0
REAL(KIND=DBL)::pi,r_0,r_1,r1,r2,e1,e2,p,cth,sth,ph,up
SAVE pi
REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::p1,p2
REAL(KIND=DBL),INTENT(OUT)::w
r_0=0
r_1=1
      
IF(init.EQ.0)THEN
   CALL Helac_mypi(pi)
   init=1
ENDIF
       
w=0
IF(q.LE.q1+q2)RETURN
       
e1=(q**2+q1**2-q2**2)/q/2
e2=q-e1
p=Phegas_clambda(q,q1,q2)/q/2
up=1-ptcut**2/p**2
IF(up.LT.0)RETURN
up=DSQRT(up)
w=p/q/(4*pi)*up
       
IF(idir.EQ.1)RETURN
       
p1(4)=e1
p2(4)=e2
     
r1=Helac_rnmy(0)
r2=Helac_rnmy(0)

cth=(2*r1-1)*up
sth=DSQRT(1-cth*cth)
ph=2*pi*r2
       
p1(3)=p*cth
p1(2)=p*sth*DSIN(ph)
p1(1)=p*sth*DCOS(ph)
       
p2(1:3)=-p1(1:3)

END SUBROUTINE Phegas_twobody_had

SUBROUTINE Phegas_twobody_had_pt(idir,q,q1,q2,pthad,p1,p2,w)
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::q,q1,q2,pthad
INTEGER::init=0
REAL(KIND=DBL)::pi,r1,r2,e1,e2,p,cth,sth,ph,up
SAVE pi,init
REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::p1,p2
REAL(KIND=DBL),INTENT(OUT)::w
IF(init.EQ.0)THEN
   CALL Helac_mypi(pi)
   init=1
ENDIF       
w=0
IF(q.LE.q1+q2)RETURN
       
e1=(q**2+q1**2-q2**2)/q/2
e2=q-e1
p=Phegas_clambda(q,q1,q2)/q/2
IF(p.LT.pthad)RETURN
w=1/q/(4*pi)*DSQRT(pthad**2/(p**2-pthad**2))
       
IF(idir.EQ.1)RETURN
       
p1(4)=e1
p2(4)=e2
up=1d0
r1=Helac_rnmy(0)
IF(r1.GT.0.5d0)up=-1d0
r2=Helac_rnmy(0)

sth=pthad/p
cth=up*DSQRT(1-sth*sth)
ph=2*pi*r2
       
p1(3)=p*cth
p1(2)=p*sth*DSIN(ph)
p1(1)=p*sth*DCOS(ph)
       
p2(1:3)=-p1(1:3)

END SUBROUTINE Phegas_twobody_had_pt
       
SUBROUTINE Phegas_tgen(idir,xl,xu,gm,xss,w)
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::xl,xu,gm
REAL(KIND=DBL)::ex,r,r_0,r_1
REAL(KIND=DBL),INTENT(INOUT)::xss
REAL(KIND=DBL),INTENT(OUT)::w
INTEGER::init=0
SAVE ex,init
       
r_0=0
r_1=1

ex=dnou(1)/100d0
IF(xl.EQ.0.OR.xu.EQ.0)ex=1d0
IF(init.EQ.0)THEN
   WRITE(nunit2,*)'Phegas_tgen: ex=',ex,1-ex
   init=1
ENDIF
       
IF(idir.EQ.0)THEN
   r=Helac_rnmy(0)
   xss=(r * xu**ex + (r_1-r) * xl**ex)**(r_1/ex)
ENDIF
       
w=(xu**ex-xl**ex) * xss**(1-ex) /ex

END SUBROUTINE Phegas_tgen

SUBROUTINE Phegas_gen_had(idir,x0,x1,xss,w)
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::x0,x1
REAL(KIND=DBL),INTENT(INOUT)::xss
REAL(KIND=DBL),INTENT(OUT)::w
REAL(KIND=DBL)::xu,xl,r_0,r_1,ex,r
r_0=0
r_1=1
ex=dnou(1)
xu=x1**2
xl=x0**2
IF(idir.EQ.0)THEN
	r=Helac_rnmy(0)
	xss=( xu**ex *  r + (r_1-r) * xl**ex )**(r_1/ex)
ENDIF
w=xss**(r_1-ex)*(xu**ex-xl**ex)/ex
END SUBROUTINE Phegas_gen_had
       
SUBROUTINE Phegas_gen(idir,x0,x1,xm,gm,xss,w)
!       include 'declare.h'
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::xm,gm,x0,x1
REAL(KIND=DBL),INTENT(INOUT)::xss
REAL(KIND=DBL),INTENT(OUT)::w
REAL(KIND=DBL)::r_0,r_1,r,ex,xu,xl
r_0=0
r_1=1
       
IF(xm.GT.0.AND.gm.GT.0.)THEN
   xu=DATAN( (x1**2-xm**2)/(xm*gm) )
   xl=DATAN( (x0**2-xm**2)/(xm*gm) )
   IF(idir.EQ.0)THEN
      r=Helac_rnmy(0)
      xss=xm**2+xm*gm*TAN(r*xu+(r_1-r)*xl)
   ENDIF

   w=(xss-xm**2)**2+xm**2*gm**2
   w=w/xm/gm*(xu-xl)

ELSE
      
   ex=dnou(1)/100d0
   IF(gm.LT.0)ex=dnou(95)/100d0
   IF(x1.EQ.0.OR.x0.EQ.0)ex=dnou(1)
   xu=x1**2
   xl=x0**2
   IF(idir.EQ.0)THEN
      r=Helac_rnmy(0)
      xss=( xu**ex *  r + (r_1-r) * xl**ex )**(r_1/ex)
   ENDIF
   w=xss**(r_1-ex)*(xu**ex-xl**ex)/ex
       
ENDIF

END SUBROUTINE Phegas_gen

SUBROUTINE Phegas_gen1(idir,x0,x1,xm,gm,xn,gn,q,y0,xss,w)
!      include 'declare.h'
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::x0,x1,xm,gm,xn,gn,y0
REAL(KIND=DBL),INTENT(IN)::q
REAL(KIND=DBL),INTENT(INOUT)::xss
REAL(KIND=DBL),INTENT(OUT)::w
REAL(KIND=DBL)::r_1,r_0,xu1,xl1,xu2,xl2,yu1,yl1,yu2,yl2,y1,vol1,vol2,alpha1,alpha2,r,ex,xu,xl
r_0=0
r_1=1
       
IF(xm.GT.0.AND.gm.GT.0.AND.x1/xm.GT.r_1/100)THEN
   xu1=DATAN( (x1**2-xm**2)/(xm*gm) )
   xl1=DATAN( ((q-xn)**2-xm**2)/(xm*gm) )
   y1=MAX(q-xm,(xn+y0)/2,y0)
   yu1=DATAN( (y1**2-xn**2)/(xn*gn) )
   yl1=DATAN( (y0**2-xn**2)/(xn*gn) )
   vol1=(xu1-xl1)/xm/gm
   alpha1=(yu1-yl1)*vol1
  
   xu2=DATAN( ((q-xn)**2-xm**2)/(xm*gm) )
   xl2=DATAN( (x0**2-xm**2)/(xm*gm) )
   y1=q-x0
   yu2=DATAN( (y1**2-xn**2)/(xn*gn) )
   yl2=DATAN( (y0**2-xn**2)/(xn*gn) )
   vol2=(xu2-xl2)/xm/gm
   alpha2=(yu2-yl2)*vol2

   alpha1=alpha1/(alpha1+alpha2)
   alpha2=1-alpha1
       
   IF(idir.EQ.0)THEN
       r=Helac_rnmy(0)
       IF(r.LT.alpha1)THEN
         r=Helac_rnmy(0)
         xss=xm**2+xm*gm*TAN(r*xu1+(r_1-r)*xl1)
         w=(xss-xm**2)**2+xm**2*gm**2
         w=w*vol1/alpha1
       ELSE
         r=Helac_rnmy(0)
         xss=xm**2+xm*gm*TAN(r*xu2+(r_1-r)*xl2)
         w=(xss-xm**2)**2+xm**2*gm**2
         w=w*vol2/alpha2
       ENDIF
    ELSE
       w=(xss-xm**2)**2+xm**2*gm**2
       IF(xss.GT.(q-xn)**2)THEN
         w=w*vol1/alpha1
       ELSE
         w=w*vol2/alpha2
       ENDIF
    ENDIF


ELSE
      
    ex=dnou(1)/100d0
    IF(gm.LT.0)ex=dnou(95)/100d0
	IF(x1.EQ.0.OR.x0.EQ.0)ex=dnou(1)
    xu=x1**2
    xl=x0**2
    IF(idir.EQ.0)THEN
       r=Helac_rnmy(0)
       xss=( xu**ex *  r + (r_1-r) * xl**ex )**(r_1/ex)
    ENDIF
    w=xss**(r_1-ex)*(xu**ex-xl**ex)/ex
       

ENDIF
     
END SUBROUTINE Phegas_gen1

      
SUBROUTINE Phegas_threebody(idir,q,q1,q2,q3,q23,p1,p2,p3,w)
!       include 'declare.h'
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::q,q1,q2,q3
REAL(KIND=DBL),INTENT(INOUT)::q23
REAL(KIND=DBL),DIMENSION(4),INTENT(INOUT)::p1,p2,p3
REAL(KIND=DBL),INTENT(OUT)::w
REAL(KIND=DBL),DIMENSION(4)::p23
INTEGER::init=0
REAL(KIND=DBL)::pi,r_0,r_1,q23min,q23max,w_23,w_2,w_1
SAVE pi,init
       
r_0=0
r_1=1
       
IF(init.EQ.0)THEN
   CALL Helac_mypi(pi)
   init=1
ENDIF
       
q23min=q2+q3
q23max=q-q1
IF(q23max.GT.q23min)THEN
!q23=0
	CALL Phegas_gen(idir,q23min,q23max,r_0,-r_1,q23,w_23)
	w_23=w_23/(2*pi)
	q23=DSQRT(q23)
	CALL Phegas_twobody(idir,q23,q2,q3,p2,p3,w_2)
ELSE
	q23=q23min
	w_23=1
	IF(idir.EQ.0)THEN
		p2(4)=q2
		p3(4)=q3
		p2(1:3)=0
		p3(1:3)=0
	ENDIF
	w_2=1
ENDIF
CALL Phegas_twobody(idir,q,q1,q23,p1,p23,w_1)

       
w=w_23*w_1*w_2
      
IF(idir.EQ.1)RETURN
CALL boostl(q23,p23,p2(1:4))
CALL boostl(q23,p23,p3(1:4))
       
END SUBROUTINE Phegas_threebody
       
FUNCTION Phegas_clambda(x1,x2,x3)
!      include 'declare.h'
!       include 'common_print.h'
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::x1,x2,x3
REAL(KIND=DBL)::Phegas_clambda,fgf      
!      if(x2.lt.eps(0).and.x3.lt.eps(0))then
!       clambda=x1
!       return
!      endif
!      if(x2.lt.eps(0))then
!       clambda=x1-x3
!       return
!      endif
!      if(x3.lt.eps(0))then
!       clambda=x1-x2
!       return
!      endif
fgf=(x1+x2+x3)*(x1-x2+x3)*(x1+x2-x3)*(x1-x2-x3)
IF(fgf.LT.0)THEN
   WRITE(nunit2,*)'warning in Phegas_clambda: fgf,x1,x2,x3',fgf,x1,x2,x3
   STOP
ENDIF
Phegas_clambda=DSQRT(fgf)
       
END FUNCTION Phegas_clambda
       

SUBROUTINE Phegas_setmin(cmas)
IMPLICIT NONE
INTEGER,DIMENSION(20)::y0,y1,y2
REAL(KIND=DBL),DIMENSION(*),INTENT(INOUT)::cmas
INTEGER::init=0,l
INTEGER::ix1,ix2,k,i0,j0,i1,i2,j1,j2
REAL(KIND=DBL)::e1gg,e2gg,p1,p2,xgg,rm1,rm2
REAL(KIND=DBL)::tranmas,massres
REAL(KIND=DBL),DIMENSION(3:20)::ec_old,c1_old,c2_old,ptc_old
REAL(KIND=DBL),DIMENSION(3:20,3:20)::cc_old
SAVE init,ec_old,c1_old,c2_old,ptc_old,cc_old
IF(GLOBALINIT_setmin.EQ.0)THEN
	init=0
ENDIF
IF(init.EQ.0)THEN
     y0(1:20)=0
     y1(1:20)=0
     y2(1:20)=0
     ec_old(3:20)=ec(3:20)
     c1_old(3:20)=c1(3:20)
     c2_old(3:20)=c2(3:20)
     ptc_old(3:20)=ptc(3:20)
     cc_old(3:20,3:20)=cc(3:20,3:20)
     init=1
ENDIF
IF(emep_ISR_shower.EQ.1)THEN
   ec(3:20)=0d0 ! a nonzero small value to avoid IR divergence
   c1(3:20)=1d0
   c2(3:20)=1d0
   ptc(3:20)=0d0 ! a nonzero small value to avoid IR divergence
   cc(3:20,3:20)=1d0
ENDIF
DO l=1,n-1
     ix1=0
     DO WHILE(ix1.EQ.0)
        CALL Helac_id(n,l,ix1,i0)
        IF(ix1.EQ.1)CYCLE

        IF(MOD(i0,4).NE.0)CYCLE   ! excluding the initial states
        IF(l.EQ.1)THEN
            CALL Helac_bin(n,i0,y0)
            DO k=1,n
                 IF(y0(k).NE.0)j0=n-k+1
            ENDDO
            cmas(i0)=parmas(ifl(j0))
        ENDIF
        IF(l.LT.2)CYCLE
        ix2=0
        DO WHILE(ix2.EQ.0)
            CALL Helac_id2(n,l,i0,ix2,i1,i2)
            IF(ix2.EQ.1)CYCLE
       
            IF(Helac_level(n,i1).EQ.1.AND.Helac_level(n,i2).EQ.1)THEN
                CALL Helac_bin(n,i1,y1)
                CALL Helac_bin(n,i2,y2)
                DO k=1,n
                   IF(y1(k).NE.0)j1=n-k+1
                   IF(y2(k).NE.0)j2=n-k+1
                ENDDO
                rm1=parmas(ifl(j1))
                rm2=parmas(ifl(j2))
                IF(rm1.EQ.0.AND.rm2.EQ.0)THEN
                   e1gg=ec(j1)
                   e2gg=ec(j2)
                   p1=0
                   p2=0
                   IF(e1gg.GT.0)p1=e1gg*DSQRT(DMAX1(1-(rm1/e1gg)**2,0d0))
                   IF(e2gg.GT.0)p2=e2gg*DSQRT(DMAX1(1-(rm2/e2gg)**2,0d0))
                   xgg=rm1**2+rm2**2+2*e1gg*e2gg-2*p1*p2*cc(j1,j2)
                   cmas(i0)=MAX(DSQRT(xgg),gmas(j1,j2),cmas(i1)+cmas(i2))
                ELSE
                   cmas(i0)=MAX(gmas(j1,j2),cmas(i1)+cmas(i2))
                ENDIF
           ELSE
                cmas(i0)=MAX(cmas(i0),cmas(i1)+cmas(i2))
           ENDIF
       ENDDO
    ENDDO
ENDDO
tranmas=0
IF(.NOT.dSigmadPTQ)THEN
	DO l=3,n
		tranmas=tranmas+DSQRT(parmas(ifl(l))**2+ptc(l)**2)
	ENDDO
ELSE
	IF(nhad.EQ.4.AND.n.EQ.6)THEN
		tranmas=DSQRT((parmas(ifl(3))+parmas(ifl(4)))**2+PTFirst**2)&
		+DSQRT((parmas(ifl(5))+parmas(ifl(6)))**2+PTFirst**2)
	ELSEIF(n-nhad.NE.0)THEN
                massres=0d0
                DO l=5,n
                   massres=massres+parmas(ifl(l))
                ENDDO
		tranmas=DSQRT((parmas(ifl(3))+parmas(ifl(4)))**2+PTFirst**2)&
		+DSQRT(massres**2+PTFirst**2)
!		tranmas=DSQRT((parmas(ifl(3))+parmas(ifl(4)))**2+PTFirst**2)*2
	ELSE
                massres=0d0
                DO l=5,n
                   massres=massres+parmas(ifl(l))
                ENDDO
		tranmas=DSQRT(parmas(ifl(3))**2+PTFirst**2)+DSQRT(massres**2+PTFirst**2)
	ENDIF
ENDIF
cmas(2**n-4)=MAX(cmas(2**n-4),tranmas)
IF(emep_ISR_shower.EQ.1)THEN
   ec(3:20)=ec_old(3:20)
   c1(3:20)=c1_old(3:20)
   c2(3:20)=c2_old(3:20)
   ptc(3:20)=ptc_old(3:20)
   cc(3:20,3:20)=cc_old(3:20,3:20)
ENDIF
END SUBROUTINE Phegas_setmin
       
SUBROUTINE Phegas_setmax(tm,cmas,xp1,xp2)
REAL(KIND=DBL),DIMENSION(*),INTENT(IN)::cmas
REAL(KIND=DBL),DIMENSION(*),INTENT(INOUT)::tm
REAL(KIND=DBL),INTENT(IN)::xp1,xp2
REAL(KIND=DBL)::exp1,exp2
INTEGER::init=0,i,i0,l
REAL(KIND=DBL)::w,cc1,cc2,rm1,rm2,q2,e1,e2,p1,p2,dtp,alpha
REAL(KIND=DBL),DIMENSION(3:20)::ec_old,c1_old,c2_old
SAVE init,ec_old,c1_old,c2_old
IF(GLOBALINIT_setmax.EQ.0)THEN
	init=0
ENDIF
IF(init.EQ.0)THEN
   ec_old(3:20)=ec(3:20)
   c1_old(3:20)=c1(3:20)
   c2_old(3:20)=c2(3:20)
ENDIF
w=2*scalar_product(phegas_pmom(1,1:4),phegas_pmom(2,1:4))&
+scalar_product(phegas_pmom(1,1:4),phegas_pmom(1,1:4))&
+  scalar_product(phegas_pmom(2,1:4),phegas_pmom(2,1:4))
w=DSQRT(w)
IF(emep_ISR_shower.EQ.1)THEN
   ec(3:20)=0d0 ! a small nonzero value to avoid IR divergence
   c1(3:20)=1d0
   c2(3:20)=1d0
ENDIF
DO i=3,n

   l=2**(i-1)
   i0=2**n-2-l  ! excluding the 1st and ith particle

   cc1=c1(i)
   cc2=c2(i)

   rm1=parmas(ifl(1)) ! we always negelect the mass of initial partonic particle 
   rm2=parmas(ifl(i))
   exp1=ebeam(1)*xp1
   exp2=ebeam(2)*xp2
   IF(rm2.EQ.0)cc1=((exp2-exp1)+c1(i)*(exp2+exp1))/((exp2+exp1)+c1(i)*(exp2-exp1))
   q2=cmas(i0-2)**2
   e1=phegas_pmom(1,4)
   e2=(w**2+rm2**2-q2)/2/w

   IF(init.EQ.0)THEN
       WRITE(nunit2,*)'i,i0,2**n-2-l-2,q2,w,e1,e2'
       WRITE(nunit2,*) i,i0,2**n-2-l-2,q2,w,e1,e2
   ENDIF
   p1=0
   p2=0
   IF(e1.GT.0)p1=e1*DSQRT(1-(rm1/e1)**2)   ! the momentum(3-dim) of particle 1
   IF(e2.GT.0)p2=e2*DSQRT(1-(rm2/e2)**2) 
   dtp=e1/w-p1/p2*(w**2+rm2**2-q2)/(2*w**2)*cc1      
   IF(init.EQ.0)WRITE(nunit2,*)'dtp',dtp
   IF(dtp.GT.0)THEN
        alpha=1-p1**2*cc1**2/e1**2
        IF(alpha.GT.0)THEN
           q2=DMAX1(w**2+rm2**2-2*w*rm2/DSQRT(alpha),q2)
           IF(q2.GT.w**2+rm2**2-2*w*ec(i))q2=w**2+rm2**2-2*w*ec(i)
        ELSE
           q2=w**2+rm2**2-2*w*ec(i)
        ENDIF
        e2=(w**2+rm2**2-q2)/2/w
        p2=0
        IF(e2.GT.0)p2=e2*DSQRT(DMAX1(1-(rm2/e2)**2,0d0))
   ENDIF
   tm(i0)=rm1**2+rm2**2-2*e1*e2+2*p1*p2*cc1
   IF(e1*e2.GT.0.AND.rm1.GT.0.AND.rm2.GT.0.AND.&
    rm1/e1.LT.1/dnou(1000).AND.rm2/e2.LT.1/dnou(1000))THEN
        tm(i0)=rm1**2*(1-e2/e1)+rm2**2*(1-e1/e2)-rm1**4*e2/4/e1**3 &
              -rm2**4*e1/4/e2**3+rm1**2*rm2**2/2/e1/e2
   ENDIF
   IF(init.EQ.0)THEN
        WRITE(nunit2,*)'i0,rm1,rm2,e1,e2,p1,p2,cc1,tm(i0)'
        WRITE(nunit2,*) i0,rm1,rm2,e1,e2,p1,p2,cc1,tm(i0)
   ENDIF

   i0=l+2   ! including particle 2 and i

   rm1=parmas(ifl(2))
   rm2=parmas(ifl(i))
   IF(rm2.EQ.0)cc2=((exp1-exp2)+c2(i)*(exp1+exp2))/((exp1+exp2)+c2(i)*(exp1-exp2))
   q2=cmas(2**n-2-i0)**2
   e1=phegas_pmom(2,4)
   e2=(w**2+rm2**2-q2)/2/w
   IF(init.EQ.0)THEN
       WRITE(nunit2,*)'i,i0,2**n-2-l-2,q2,w,e1,e2'
       WRITE(nunit2,*) i,i0,2**n-2-l-2,q2,w,e1,e2
   ENDIF
   p1=0
   p2=0
   IF(e1.GT.0)p1=e1*DSQRT(DMAX1(1-(rm1/e1)**2,0d0))
   IF(e2.GT.0)p2=e2*DSQRT(DMAX1(1-(rm2/e2)**2,0d0))
       
   dtp=e1/w-p1/p2*(w**2+rm2**2-q2)/(2*w**2)*cc2       
   IF(dtp.GT.0)THEN
        alpha=1-p1**2*cc2**2/e1**2
        IF(alpha.GT.0)THEN
            q2=DMAX1(w**2+rm2**2-2*w*rm2/DSQRT(alpha),q2)
            IF(q2.GT.w**2+rm2**2-2*w*ec(i))q2=w**2+rm2**2-2*w*ec(i)
        ELSE
            q2=w**2+rm2**2-2*w*ec(i)
        ENDIF
        e2=(w**2+rm2**2-q2)/2/w
        p2=0
        IF(e2.GT.0)p2=e2*DSQRT(DMAX1(1-(rm2/e2)**2,0d0))
   ENDIF

       
   tm(i0)=rm1**2+rm2**2-2*e1*e2+2*p1*p2*cc2
   IF(e1*e2.GT.0.AND.rm1.GT.0.AND.rm2.GT.0.AND.&
   rm1/e1.LT.1/dnou(1000).AND.rm2/e2.LT.1/dnou(1000))THEN
        tm(i0)=rm1**2*(1-e2/e1)+rm2**2*(1-e1/e2)-rm1**4*e2/4/e1**3 &
          -rm2**4*e1/4/e2**3+rm1**2*rm2**2/2/e1/e2
   ENDIF
       
   IF(init.EQ.0)THEN
        WRITE(nunit2,*)'i0,rm1,rm2,e1,e2,p1,p2,cc2,tm(i0)'
        WRITE(nunit2,*) i0,rm1,rm2,e1,e2,p1,p2,cc2,tm(i0)
   ENDIF
ENDDO
IF(emep_ISR_shower.EQ.1)THEN
   ec(3:20)=ec_old(3:20)
   c1(3:20)=c1_old(3:20)
   c2(3:20)=c2_old(3:20)
ENDIF
init=1
END SUBROUTINE Phegas_setmax
       
SUBROUTINE Phegas_look(isis,k,i1,mm)
INTEGER,DIMENSION(8,8),INTENT(IN)::isis
INTEGER,INTENT(IN)::k,i1
INTEGER,INTENT(INOUT)::mm
INTEGER::i
mm=0
DO i=k,8
   IF(isis(i,1).EQ.i1)THEN
      IF(isis(i,7).EQ.0)THEN
		IF(ForwardQ(isis(i,3))) mm=isis(i,6) 
		IF(ForwardQ(isis(i,5))) mm=isis(i,4)
	  ELSE
		IF(ForwardQ(isis(i,3)).AND.ForwardQ(isis(i,5)))mm=isis(i,8)
		IF(ForwardQ(isis(i,5)).AND.ForwardQ(isis(i,7)))mm=isis(i,4)
		IF(ForwardQ(isis(i,3)).AND.ForwardQ(isis(i,7)))mm=isis(i,6)
	  ENDIF 
   ENDIF
ENDDO
END SUBROUTINE Phegas_look
!---------------------------------------------------------------


SUBROUTINE Phegas_givemea1a2(rq1,ra1,ra2)
!     implicit real(kind(1.q0)) (a-h,o-q,s-y)
!      implicit real(kind(1.d0)) (a-h,o-q,s-y)
!      implicit real(kind(1.d0)) (r)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::rq1
REAL(KIND=DBL),INTENT(OUT)::ra1,ra2
!common/phegas1/rq,rq2,re2,rem
REAL(KIND=DBL)::eps=1d-5,q,q2gg,e2gg,em,q1gg,a1,a2,p,pe
SAVE eps
q=rq
q2gg=rq2
e2gg=re2 
em=rem
q1gg=rq1
a1=em**2+q1gg**2-2*e2gg*(q**2+q1gg**2-q2gg**2)/2/q

!         if(q1/q.lt.eps.and.q2/q.lt.eps)then
!      p=q/2
!    .  -(q1**2+q2**2)/q/2
!    .  -2*q1**2*q2**2/q**3/2
!    .  -2*q1**2*q2**2*(q1**2+q2**2)/q**5/2
!     elseif(q1/q.lt.eps)then
!      p=(q**2-q2**2)/q/2
!    .  -(q**2+q2**2)/(q**2-q2**2)/q * q1**2/2
!    .  -2*q2**2/(q**2-q2**2)**3*q   * q1**4/2
!     elseif(q2/q.lt.eps)then
!      p=(q**2-q1**2)/q/2
!    .  -(q**2+q1**2)/(q**2-q1**2)/q * q2**2/2
!    .  -2*q1**2/(q**2-q1**2)**3*q   * q2**4/2
!     else
!      p=sqrt(q**4+q1**4+q2**4-2*q**2*q1**2
!    .          -2*q1**2*q2**2-2*q2**2*q**2)/2/q
!     endif
!     print*,'---------------------'
!     p1=p
!     print*,p1
p=Phegas_clambda(q,q1gg,q2gg)/2/q
!     print*,p,p-p1

! print*,p,sqrt(q**4+q1**4+q2**4-2*q**2*q1**2
!.          -2*q1**2*q2**2-2*q2**2*q**2)/2/q

IF(em/e2gg.LT.eps)THEN
    pe=e2gg*(1-(em/e2gg)**2/2-(em/e2gg)**4/8-(em/e2gg)**6/16)
ELSE
    pe=SQRT(e2gg**2-em**2)
ENDIF

a2=2*pe*p

ra1=a1
ra2=a2

!     if(abs(ra1+ra2).lt.1d-15)then
!      print*,a1+a2,ra1+ra2
!      ra12=a1+a2
!      ra2=ra12-ra1
!      print*,a1+a2,ra1+ra2
!      print*,q,q1,q2
!      print*,e2,em
!     endif
END SUBROUTINE Phegas_givemea1a2

SUBROUTINE Phegas_gen_x(idir,a1,a2,t1l,t1u,xgg,w)       
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::a1,a2,t1l,t1u
REAL(KIND=DBL),INTENT(INOUT)::xgg
REAL(KIND=DBL),INTENT(OUT)::w
REAL(KIND=DBL)::ex,r_0,r_1,r,f
INTEGER::init=0
SAVE ex,init
r_0=0
r_1=1
ex=dnou(1)/100d0
IF(t1l.EQ.0.OR.t1u.EQ.0)ex=dnou(1)
IF(GLOBALINIT_gen_x.EQ.0)THEN
	init=0
ENDIF
IF(init.EQ.0)THEN
  WRITE(nunit2,*)'Phegas_gen_x: ex=',ex,1-ex
  init=1
ENDIF  
IF(idir.EQ.0)THEN
  r=Helac_rnmy(0)
  f=(r * t1u**ex + (r_1-r) * t1l**ex)**(r_1/ex)
  xgg=(f-a1)/a2
ENDIF
       
w=((t1l)**ex-(t1u)**ex)/ex*(a1+a2*xgg)**(r_1-ex)/a2
END SUBROUTINE Phegas_gen_x

SUBROUTINE Phegas_gen_x2(idir,mass2,t2min,t2max,phi,w_phi)
IMPLICIT NONE
INTEGER,INTENT(IN)::idir
REAL(KIND=DBL),INTENT(IN)::mass2,t2min,t2max
REAL(KIND=DBL),INTENT(INOUT)::phi
REAL(KIND=DBL),INTENT(OUT)::w_phi
REAL(KIND=DBL)::f,r,pi,init=0,b,tem
SAVE init,pi
IF(init.EQ.0)THEN
	pi=DATAN(1d0)*4
	init=1
ENDIF
b=(t2max-t2min)/(t2min-mass2**2)
IF(idir.EQ.0)THEN
	r=Helac_rnmy(0)
	tem=DTAN(pi*(r-0.5d0))
	f=tem/DSQRT(1+b+tem**2)
	phi=2d0*DACOS(f)
	f=f**2
ELSE
	f=(DCOS(phi)+1d0)/2d0
ENDIF
w_phi=ABS(2*pi*(b*f+1)/DSQRT(1+b))
END SUBROUTINE Phegas_gen_x2

SUBROUTINE Phegas_testmom(iweight,xp1,xp2)       
IMPLICIT NONE
INTEGER,INTENT(OUT)::iweight
REAL(KIND=DBL),INTENT(IN)::xp1,xp2
REAL(KIND=DBL)::exp1,exp2
!REAL(KIND=DBL),DIMENSION(1:4):: phegas_sum
REAL(KIND=DBL),DIMENSION(20)::rm
COMPLEX(KIND=DBL),DIMENSION(20,4)::zy
COMPLEX(KIND=DBL),DIMENSION(4)::zv
INTEGER::k,i,j
REAL(KIND=DBL)::cc1,cc2,cc1c,cc2c
REAL(KIND=DBL)::xgg
REAL(KIND=DBL),DIMENSION(3:20)::c1_old,c2_old
INTEGER::init=0
SAVE init,c1_old,c2_old
IF(init.EQ.0)THEN
   c1_old(3:20)=c1(3:20)
   c2_old(3:20)=c2(3:20)
   init=1
ENDIF
IF(emep_ISR_shower.EQ.1)THEN
   c1(3:20)=1d0
   c2(3:20)=1d0
ENDIF
iweight=1
exp1=ebeam(1)*xp1
exp2=ebeam(2)*xp2
DO k=3,n
   zy(k,1)=DCMPLX(DREAL(zq(k,1)),dnou(0))*io(k)
   zy(k,2)=DCMPLX(DREAL(zq(k,2)),dnou(0))*io(k)
   zy(k,3)=zq(k,3)*io(k)
   zy(k,4)=zq(k,4)*io(k)
   cc1=phegas_pmom(k,3)/phegas_pmom(k,4)
   cc2=-phegas_pmom(k,3)/phegas_pmom(k,4)
   IF(parmas(ifl(k)).EQ.0)THEN
      cc1c=((exp2-exp1)+c1(k)*(exp2+exp1))/((exp2+exp1)+c1(k)*(exp2-exp1))
      IF(cc1.GT.cc1c)THEN
          iweight=-1
          RETURN
      ENDIF
      cc2c=((exp1-exp2)+c2(k)*(exp1+exp2))/((exp1+exp2)+c2(k)*(exp1-exp2))
      IF(cc2.GT.cc2c)THEN
          iweight=-1
          RETURN
      ENDIF
   ENDIF
ENDDO

DO i=3,n-1
   DO j=i+1,n
      zv(1:4)=zy(i,1:4)+zy(j,1:4)
      xgg=Helac_zprod(zv,zv)
      IF(xgg.LT.gmas(i,j)**2)THEN
         iweight=-1
         RETURN
      ENDIF
   ENDDO
ENDDO
IF(emep_ISR_shower.EQ.1)THEN
   c1(3:20)=c1_old(3:20)
   c2(3:20)=c2_old(3:20)
ENDIF
END SUBROUTINE Phegas_testmom
      
SUBROUTINE Phegas_bininv(ww)
INTEGER,DIMENSION(n)::y
REAL(KIND=DBL),DIMENSION(4)::p
REAL(KIND=DBL)::s
REAL(KIND=DBL),INTENT(IN)::ww
INTEGER::i,j,k
j=0
DO i=2,2**n-4,2
   IF(Helac_level(n,i).EQ.1)CYCLE
   CALL Helac_bin(n,i,y)
   p(1:4)=0
   DO k=1,n
      IF(y(k).GT.0)p(1:4)=p(1:4)+phegas_pmom(n-k+1,1:4)*io(n-k+1)
   ENDDO
   j=j+1
   s=scalar_product(p,p)
   IF(s.LT.0)THEN
      s=DSQRT(-s)
      CALL Helac_histo1(j,100,0d0,1000d0,s,ww)
   ELSE
      s=DSQRT(s)
      CALL Helac_histo1(j,100,0d0,1000d0,s,ww)
   ENDIF
ENDDO
END SUBROUTINE Phegas_bininv

SUBROUTINE Phegas_printinv
INTEGER,DIMENSION(n)::y
REAL(KIND=DBL),DIMENSION(4)::p,sum
REAL(KIND=DBL)::s
INTEGER::i,k
DO i=2,2**n-2,2
   IF(Helac_level(n,i).EQ.1)CYCLE
   CALL Helac_bin(n,i,y)
   p(1:4)=0
   DO k=1,n
      IF(y(k).GT.0)p(1:4)=p(1:4)+phegas_pmom(n-k+1,1:4)*io(n-k+1)
   ENDDO
   s=scalar_product(p,p)
   IF(s.LT.0)THEN
        WRITE(*,*)i,'-',SQRT(-s)
   ELSE
        WRITE(*,*)i,'+',SQRT(s)
   ENDIF
ENDDO
sum(1:4)=0
DO i=1,n
   sum(1:4)=sum(1:4)+phegas_pmom(i,1:4)*io(i)
ENDDO
WRITE(*,*)sum
END SUBROUTINE Phegas_printinv

SUBROUTINE Phegas_findlimits
INTEGER,DIMENSION(n)::y
REAL(KIND=DBL),DIMENSION(4)::p
INTEGER::init,i,k
REAL(KIND=DBL)::s
SAVE init
!IF(GLOBALINIT_findlimits.EQ.0)THEN
!	init=0
!ENDIF
IF(init.EQ.0)THEN
  x(1:2**n)=1d24
  init=1
ENDIF
DO i=2,2**n-2,2
   IF(Helac_level(n,i).EQ.1)CYCLE
   CALL Helac_bin(n,i,y)
   IF(y(n-1).GT.0)CYCLE
   p(1:4)=0
   DO k=1,n
      IF(y(k).GT.0)p(1:4)=p(1:4)+phegas_pmom(n-k+1,1:4)
   ENDDO
   s=scalar_product(p,p)
   IF(x(i).GT.s)THEN
       x(i)=s
   ENDIF
ENDDO
END SUBROUTINE Phegas_findlimits

SUBROUTINE Phegas_printlimits
INTEGER,DIMENSION(n)::y
INTEGER::i
DO i=2,2**n-2,2
   IF(Helac_level(n,i).EQ.1)CYCLE
   CALL Helac_bin(n,i,y)
   IF(y(n-1).GT.0)CYCLE
   WRITE(*,*)'limit',i,x(i)
ENDDO
END SUBROUTINE Phegas_printlimits

FUNCTION Phegas_tfunp(rxss)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::rxss
REAL(KIND=DBL)::Phegas_tfunp
REAL(KIND=DBL)::xgg,tplus,qgg,q2gg,e2gg,emgg,p2gg,xx1gg,xx2gg
!       common/phegas1/rq,rq2,re2,rem
qgg=q   !rq
q2gg=q2 !rq2
e2gg=e2 !re2
emgg=em  !rem
xgg=rxss
p2gg=ABS(e2gg*DSQRT(1-emgg/e2gg**2))   ! the momentum of particle 2 in q1+q2 c.m. rest frame
xx1gg=qgg+xgg**2/qgg-q2gg**2/qgg
xx2gg=Phegas_clambda(qgg,xgg,q2gg)/qgg
tplus=emgg+xgg**2-e2gg*xx1gg+p2gg*xx2gg
Phegas_tfunp=tplus
END FUNCTION Phegas_tfunp

SUBROUTINE Phegas_tfunpm2(Q,Q1,Q2,p2,p3,costh,t2max,t2min)
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::Q,Q1,Q2,costh
REAL(KIND=DBL),DIMENSION(1:4),INTENT(IN)::p2,p3
REAL(KIND=DBL),INTENT(OUT)::t2max,t2min
REAL(KIND=DBL)::xx1gg,q32,q3gg,q2gg,E2gg,E3gg,q2q3,sinth,costhp,sinthp,part1,part2
q2q3=scalar_product(p2(1:4),p3(1:4))
q32=scalar_product(p3(1:4),p3(1:4))
q2gg=DSQRT(p2(1)**2+p2(2)**2+p2(3)**2)
E3gg=p3(4)
q3gg=ABS(E3gg*DSQRT(1-q32/E3gg**2))
E2gg=p2(4)
IF(1-costh**2.LT.0)THEN
	sinth=0
ELSE
	sinth=DSQRT(1-costh**2)
ENDIF
costhp=(E2gg*E3gg-q2q3)/q2gg/q3gg
IF(1-costhp**2.LT.0)THEN
	IF(costhp.LT.0)THEN
		costhp=-1d0
	ELSE
		costhp=1d0
	ENDIF
	sinthp=0
ELSE
	sinthp=DSQRT(1-costhp**2)
ENDIF
xx1gg=Phegas_clambda(Q,Q1,Q2)/Q
part1=Q2**2+q32-(E3gg/Q)*(Q**2-Q1**2+Q2**2)-q3gg*xx1gg*costh*costhp
part2=q3gg*sinthp*xx1gg*sinth
t2max=part1+part2
t2min=part1-part2
RETURN
END SUBROUTINE Phegas_tfunpm2

!-----------------------------------------------

FUNCTION Phegas_costh(rq1,rt1)
!      implicit real(kind(1.q0)) (a-h,o-q,s-y)
!       implicit real(kind(1.d0)) (a-h,o-q,s-y)
!       implicit real(kind(1.d0)) (r)
!       common/phegas1/rq,rq2,re2,rem
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::rq1,rt1
REAL(KIND=DBL)::Phegas_costh
REAL(KIND=DBL)::qgg,q2gg,e2gg,emgg,t1gg,q1gg,agg,bgg
qgg=rq
q2gg=rq2
e2gg=re2
emgg=rem
t1gg=rt1
q1gg=rq1
agg=emgg**2+q1gg**2-2*e2gg*(qgg**2+q1gg**2-q2gg**2)/2/qgg
bgg=2*DSQRT(e2gg**2-emgg**2)*DSQRT(qgg**4+q1gg**4+q2gg**4-2*qgg**2*q1gg**2 &
     -2*q1gg**2*q2gg**2-2*q2gg**2*qgg**2)/2/qgg
Phegas_costh=(t1gg-agg)/bgg
END FUNCTION Phegas_costh

!-----------------------------------------------

SUBROUTINE Phegas_solve2(rt1,ry_1,ry_2,ra,istop)
!      implicit real(kind(1.q0)) (a-h,o-q,s-y)
!       implicit real(kind(1.d0)) (a-h,o-q,s-y)
!       implicit real(kind(1.d0)) (r)
!       common/phegas1/rq,rq2,re2,rem
IMPLICIT NONE
REAL(KIND=DBL),INTENT(IN)::rt1
REAL(KIND=DBL),INTENT(OUT)::ry_1,ry_2,ra
INTEGER,INTENT(OUT)::istop
REAL(KIND=DBL)::qgg,q2gg,e2gg,emgg,t1gg,a3gg,a2gg,a11gg,a1gg,a12gg,agg,bgg,cgg,dgg,&
                y_1,y_2,x_1,x_2,d1gg
INTEGER::flag
qgg=rq
q2gg=rq2
e2gg=re2
emgg=rem
t1gg=rt1

istop=0
a3gg=1-e2gg/qgg
a2gg=(e2gg/qgg)**2-(emgg/qgg)**2
a11gg=t1gg-emgg*emgg
a12gg=e2gg/qgg*(qgg**2-q2gg**2)
a1gg=a11gg+a12gg
       
agg=1-2*e2gg/qgg+(emgg/qgg)**2
IF(ABS(agg).LE.EPSILON(agg))THEN 
	agg=0
ENDIF
bgg=-2*( a3gg*a1gg-a2gg*(qgg**2+q2gg**2) )
cgg=a11gg**2+2*a11gg*a12gg+(emgg/qgg)**2*(qgg**2-q2gg**2)**2

flag=0       
IF(agg.EQ.0)THEN
	IF(bgg.EQ.0)THEN
		istop=1
		flag=1
	ELSE
		y_1=-cgg/bgg
		y_2=bgg
		flag=1
    ENDIF
ENDIF

IF(flag.EQ.0)THEN       
	d1gg=bgg**2-4*agg*cgg
	dgg=4*(e2gg**2-emgg**2)/qgg**2 &
		*((t1gg-agg*qgg**2-q2gg**2)**2-4*agg*qgg**2*q2gg**2)
       
	IF(dgg.LT.0)THEN 
		WRITE(*,*)'d<0 in solve2',dgg,d1gg,agg,bgg,cgg
		WRITE(*,*)'q,q2,t1,e2,em',qgg,q2gg,t1gg,e2gg,emgg
		istop=2
		flag=1
	ENDIF
    
	IF(flag.EQ.0)THEN          
		IF(bgg.GT.0)THEN
			x_1=(-bgg-SQRT(dgg))/(2*agg)
			x_2=cgg/agg/x_1
		ELSE
			x_1=(-bgg+SQRT(dgg))/(2*agg)
			x_2=cgg/agg/x_1
		ENDIF

		y_1=MIN(x_1,x_2)
		y_2=MAX(x_1,x_2)
	ENDIF
ENDIF
! 1      continue
ry_1=y_1
ry_2=y_2
ra=agg

END SUBROUTINE Phegas_solve2

FUNCTION SoniumQ(ii)
IMPLICIT NONE
INTEGER,INTENT(IN)::ii
LOGICAL::SoniumQ
INTEGER,DIMENSION(n)::yy
INTEGER::k
SoniumQ=.FALSE.
IF(Helac_level(n,ii).NE.2)RETURN
CALL Helac_bin(n,ii,yy)
DO k=3,nhad
	IF(Quarkonium(k,1).EQ.0)CYCLE
	IF(yy(n+1-Quarkonium(k,2)).NE.0.AND.yy(n+1-Quarkonium(k,3)).NE.0)THEN
		SoniumQ=.TRUE.
		RETURN
	ENDIF
ENDDO
END FUNCTION SoniumQ

FUNCTION SoniumQ2(ii)
IMPLICIT NONE
INTEGER,INTENT(IN)::ii
LOGICAL::SoniumQ2
INTEGER,DIMENSION(n)::yy
INTEGER::k1,k2
SoniumQ2=.FALSE.
IF(n-nhad.LE.1)RETURN
IF(Helac_level(n,ii).NE.4)RETURN
CALL Helac_bin(n,ii,yy)
DO k1=3,nhad-1
	IF(Quarkonium(k1,1).EQ.0)CYCLE
	IF(yy(n+1-Quarkonium(k1,2)).EQ.0.OR.yy(n+1-Quarkonium(k1,3)).EQ.0)CYCLE
	DO k2=k1+1,nhad
		IF(Quarkonium(k2,1).EQ.0)CYCLE
		IF(yy(n+1-Quarkonium(k2,2)).NE.0.AND.yy(n+1-Quarkonium(k2,3)).NE.0)THEN
			SoniumQ2=.TRUE.
			RETURN
		ENDIF
	ENDDO
ENDDO
END FUNCTION SoniumQ2

FUNCTION SoniumQ3(ii)
IMPLICIT NONE
INTEGER,INTENT(IN)::ii
LOGICAL::SoniumQ3
INTEGER,DIMENSION(n)::yy
INTEGER::k1,k2,k3
SoniumQ3=.FALSE.
IF(n-nhad.LE.2)RETURN
IF(Helac_level(n,ii).NE.6)RETURN
CALL Helac_bin(n,ii,yy)
DO k1=3,nhad-2
	IF(Quarkonium(k1,1).EQ.0)CYCLE
	IF(yy(n+1-Quarkonium(k1,2)).EQ.0.OR.yy(n+1-Quarkonium(k1,3)).EQ.0)CYCLE
	DO k2=k1+1,nhad-1
		IF(Quarkonium(k2,1).EQ.0)CYCLE
		IF(yy(n+1-Quarkonium(k2,2)).EQ.0.OR.yy(n+1-Quarkonium(k2,3)).EQ.0)CYCLE
		DO k3=k2+1,nhad
			IF(Quarkonium(k3,1).EQ.0)CYCLE
			IF(yy(n+1-Quarkonium(k3,2)).NE.0.AND.yy(n+1-Quarkonium(k3,3)).NE.0)THEN
				SoniumQ3=.TRUE.
				RETURN
			ENDIF
		ENDDO
	ENDDO
ENDDO
END FUNCTION SoniumQ3

FUNCTION SoniumQ4(i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::i1,i2
INTEGER,DIMENSION(n)::yy1   !,yy2
LOGICAL::SoniumQ4
INTEGER::k,hadnum
SoniumQ4=.FALSE.
IF(Helac_level(n,i1).NE.3)RETURN
CALL Helac_bin(n,i1,yy1)
!CALL Helac_bin(n,i2,yy2)
hadnum=0
DO k=3,nhad
	IF(Quarkonium(k,1).EQ.0)CYCLE
	IF(yy1(n+1-Quarkonium(k,2)).NE.0.OR.yy1(n+1-Quarkonium(k,3)).NE.0)hadnum=hadnum+1
ENDDO
IF(hadnum.EQ.3)SoniumQ4=.TRUE.
RETURN
END FUNCTION SoniumQ4

SUBROUTINE Phegas_genQ2(idir,q2,q1,i1,i2)
IMPLICIT NONE
INTEGER,INTENT(IN)::idir,i1,i2
REAL(KIND=DBL),INTENT(IN)::q1
REAL(KIND=DBL),INTENT(OUT)::q2
INTEGER,DIMENSION(n)::y1   !,y2
REAL(KIND=DBL)::m11,m12,m21,m22
INTEGER::k,k2
CALL Helac_bin(n,i1,y1)
!CALL Helac_bin(n,i2,y2)
m11=-1
m12=-1
m21=-1
m22=-1
k2=0
DO k=3,nhad-1
	IF(Quarkonium(k,1).EQ.0)CYCLE
	IF(y1(n+1-Quarkonium(k,2)).NE.0)THEN
		m11=parmas(ifl(Quarkonium(k,2)))
		m12=parmas(ifl(Quarkonium(k,3)))
		k2=k
		EXIT
	ENDIF
	IF(y1(n+1-Quarkonium(k,3)).NE.0)THEN
		m11=parmas(ifl(Quarkonium(k,3)))
		m12=parmas(ifl(Quarkonium(k,2)))
		k2=k
		EXIT
	ENDIF
ENDDO
IF(k2.EQ.0)THEN
	PRINT *,"Wrong (1) in Phegas_genQ2 ! STOP !"
	STOP
ENDIF
DO k=k2+1,nhad
	IF(Quarkonium(k,1).EQ.0)CYCLE
	IF(y1(n+1-Quarkonium(k,2)).NE.0)THEN
		m21=parmas(ifl(Quarkonium(k,2)))
		m22=parmas(ifl(Quarkonium(k,3)))
		EXIT
	ENDIF
	IF(y1(n+1-Quarkonium(k,3)).NE.0)THEN
		m21=parmas(ifl(Quarkonium(k,3)))
		m22=parmas(ifl(Quarkonium(k,2)))
		EXIT
	ENDIF
ENDDO
IF(m21.EQ.-1)THEN
	PRINT *,"Wrong (2) in Phegas_genQ2 ! STOP !"
	STOP
ENDIF
q2=m12**2+m22**2+(q1**2-m11**2-m21**2)*m12*m22/m11/m21
RETURN
END SUBROUTINE Phegas_genQ2

SUBROUTINE Phegas_genQ3(idir,q3,q1,q2,i1,i2,i3)
IMPLICIT NONE
INTEGER,INTENT(IN)::idir,i1,i2,i3
REAL(KIND=DBL),INTENT(IN)::q1,q2
REAL(KIND=DBL),INTENT(OUT)::q3
REAL(KIND=DBL)::qq1,qq2
INTEGER::ii
q3=0
RETURN
END SUBROUTINE Phegas_genQ3

FUNCTION ForwardQ(ii)
IMPLICIT NONE
INTEGER,INTENT(IN)::ii
LOGICAL::ForwardQ
INTEGER::jj,k
INTEGER,DIMENSION(n)::yy
ForwardQ=.FALSE.
IF(ii.EQ.2)THEN
	ForwardQ=.TRUE.
	RETURN
ENDIF
jj=ii
IF(MOD(ii/2,2).NE.0)jj=ii-2
CALL Helac_bin(n,jj,yy)
DO k=3,n
	IF(yy(n+1-k).NE.0.AND.Quarkonium3(k).EQ.0)RETURN
ENDDO
ForwardQ=.TRUE.
RETURN
END FUNCTION ForwardQ

FUNCTION ContainHad(ii)
IMPLICIT NONE
INTEGER,DIMENSION(n),INTENT(IN)::ii
LOGICAL::ContainHad
INTEGER::k
ContainHad=.TRUE.
DO k=3,n
	IF(ii(n+1-k).NE.0.AND.Quarkonium3(k).NE.0)RETURN
ENDDO
ContainHad=.FALSE.
RETURN
END FUNCTION ContainHad

FUNCTION SubtractHad(ii)
IMPLICIT NONE
INTEGER,INTENT(IN)::ii
INTEGER::SubtractHad
INTEGER,DIMENSION(n)::yy
INTEGER::k
CALL Helac_bin(n,ii,yy)
SubtractHad=ii
IF(yy(n-1).NE.0)SubtractHad=SubtractHad-2
DO k=3,n
	IF(yy(n+1-k).NE.0.AND.Quarkonium3(k).NE.0)SubtractHad=SubtractHad-2**(k-1)
ENDDO
RETURN
END FUNCTION SubtractHad
END MODULE Phegas_mod
