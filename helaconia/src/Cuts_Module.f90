MODULE Cuts_Module
USE Helac_Global
USE Kinetic_Func
USE  Constants
IMPLICIT NONE
CONTAINS
SUBROUTINE readcuts_HC
IMPLICIT NONE
CHARACTER(len=24)::file
REAL(KIND(1d0))::cutoff,pi
INTEGER::i,j,i1,j1,l,k
REAL(KIND(1d0))::r
REAL(KIND(1d0))::ptl,ptq,ptbottom,pttop,ptg,etal,etaq,etabottom,etatop,etag,&
                drll,drlq,drqq,drgx,gqq,gll,glq,ggx,drqb,gqb,drbb,gbb,ptmiss,ptcharm,&
                etacharm,etabonia,etaconia,etaBconia,ptconia,ptbonia,ptBconia
REAL(KIND(1d0))::maxptl,maxptq,maxptbottom,maxpttop,maxptcharm,maxptg,maxptconia,maxptbonia,maxptBconia
REAL(KIND(1d0))::ycl,ycq,ycbottom,yctop,ycg,yccharm,ycbonia,ycconia,ycBconia,&
                 ycllow,ycqlow,ycbottomlow,yctoplow,ycglow,yccharmlow,ycbonialow,ycconialow,ycBconialow
REAL(KIND(1d0))::xFcl,xFcq,xFcbottom,xFctop,xFcg,xFccharm,xFcbonia,xFcconia,xFcBconia,&
     xFcllow,xFcqlow,xFcbottomlow,xFctoplow,xFcglow,xFccharmlow,xFcbonialow,xFcconialow,xFcBconialow
REAL(KIND(1d0))::gbeamq
INTEGER::iounit,flag=0
LOGICAL::lexist
! open default input file
INQUIRE(FILE=TRIM(input_dir)//"default.inp",EXIST=lexist)
IF(.NOT.lexist)THEN
	PRINT *,"Warning: the file default.inp does not exist ! STOP !"
	STOP
ENDIF
INQUIRE(FILE=TRIM(input_dir)//"default.inp",OPENED=lexist)
IF(lexist)THEN
    INQUIRE(FILE=TRIM(input_dir)//"default.inp",NUMBER=iounit)
	IF(iounit.NE.udefault)THEN
		PRINT *,"WARNING: the default.inp has been linked with another unit ! Close and reopen !"
		CLOSE(UNIT=iounit)
		OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
	ENDIF
ELSE
	OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
ENDIF
! open user's input file
IF(TRIM(Input_File)/="default.inp")THEN
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
	IF(.NOT.lexist)THEN
		PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
		STOP
	ENDIF
	INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
	IF(lexist)THEN
		INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
		IF(iounit.NE.uinput)THEN
			PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
			CLOSE(UNIT=iounit)
			OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
		ENDIF
	ELSE
		OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
	ENDIF
ELSE
    flag=1
ENDIF

cutoff=readvalue_r("cutoffp",flag)

PRINT *,"WARNING CUTOFF SET:",cutoff
ptc(3:20)=cutoff
maxptc(3:20)=-1d0
drc(3:20,3:20)=0
etac(3:20)=20
yycut(3:20)=1d9
IF(absrap)THEN
   yycutlow(3:20)=0d0
ELSE
   yycutlow(3:20)=-1d9
ENDIF
xFcut(3:20)=1d0
xFcutlow(3:20)=-1d0
xFcutflag=.FALSE.
ec(3:20)=cutoff
c1(3:20)=1d0
c2(3:20)=1d0
cc(3:20,3:20)=1d0
gmas(3:20,3:20)=cutoff
gbeammass(1:2,3:20)=0d0
pi=DACOS(-1d0)
y1cup=30d0
y1clow=0d0
! maxrapidity for the first final particle
y1cup=readvalue_r("maxy1c",flag)
! minrapidity for the first final particle
y1clow=readvalue_r("miny1c",flag)
! minimum lepton pt
ptl=readvalue_r("minptl",flag)
! maximum lepton pt
maxptl=readvalue_r("maxptl",flag)
!PRINT *,'ptl',ptl
! minimum quark pt
ptq=readvalue_r("minptq",flag)
! maximum quark pt
maxptq=readvalue_r("maxptq",flag)
! minimum charm pt
ptcharm=readvalue_r("minptc",flag)
! maximum charm pt
maxptcharm=readvalue_r("maxptc",flag)
!PRINT *,'ptq',ptq
! minimum bottom pt
ptbottom=readvalue_r("minptb",flag)
! maximum bottom pt
maxptbottom=readvalue_r("maxptb",flag)
!PRINT *,'ptbottom',ptbottom
! minimum top pt
pttop=readvalue_r("minptt",flag)
! maximum top pt
maxpttop=readvalue_r("maxptt",flag)
!PRINT *,'pttop',pttop
! minimum photon pt
ptg=readvalue_r("minptp",flag)
! maximum photon pt
maxptg=readvalue_r("maxptp",flag)
! minimum charmonia pt
ptconia=readvalue_r("minptconia",flag)
! maximum charmonia pt
maxptconia=readvalue_r("maxptconia",flag)
! minimum bottomonia pt
ptbonia=readvalue_r("minptbonia",flag)
! maximum bottomonia pt
maxptbonia=readvalue_r("maxptbonia",flag)
! minimum Bc pt
ptBconia=readvalue_r("minptBconia",flag)
! maximum Bc pt
maxptBconia=readvalue_r("maxptBconia",flag)
!PRINT *,'ptg',ptg
! maximum rapidity lepton
etal=readvalue_r("maxrapl",flag)
! maximum y rapidity lepton
ycl=readvalue_r("maxyrapl",flag)
! minimum y rapidity lepton
ycllow=readvalue_r("minyrapl",flag)
! maximum Feynman parameter xF
xFcl=readvalue_r("maxxFl",flag)
! minimum Feynman parameter xF
xFcllow=readvalue_r("minxFl",flag)
!PRINT *,'etal',etal
! maximum rapidity quark 
etaq=readvalue_r("maxrapq",flag)
! maximum y rapidity quark
ycq=readvalue_r("maxyrapq",flag)
! minimum y rapidity quark
ycqlow=readvalue_r("minyrapq",flag)
! maximum Feynman parameter xF
xFcq=readvalue_r("maxxFq",flag)
! minimum Feynman parameter xF
xFcqlow=readvalue_r("minxFq",flag)
! maximum rapidity charm
etacharm=readvalue_r("maxrapc",flag)
! maximum y rapidity charm
yccharm=readvalue_r("maxyrapc",flag)
! minimum y rapidity charm
yccharmlow=readvalue_r("minyrapc",flag)
! maximum Feynman parameter xF
xFccharm=readvalue_r("maxxFc",flag)
! minimum Feynman parameter xF
xFccharmlow=readvalue_r("minxFc",flag)
!PRINT *,'etaq',etaq
! maximum rapidity bottom 
etabottom=readvalue_r("maxrapb",flag)
! maximum y rapidity bottom
ycbottom=readvalue_r("maxyrapb",flag)
! minimum y rapidity bottom
ycbottomlow=readvalue_r("minyrapb",flag)
! maximum Feynman parameter xF
xFcbottom=readvalue_r("maxxFb",flag)
! minimum Feynman parameter xF
xFcbottomlow=readvalue_r("minxFb",flag)
!PRINT *,'etabottom',etabottom
! maximum rapidity top 
etatop=readvalue_r("maxrapt",flag)
! maximum y rapidity top
yctop=readvalue_r("maxyrapt",flag)
! minimum y rapidity top
yctoplow=readvalue_r("minyrapt",flag)
! maximum Feynman parameter xF
xFctop=readvalue_r("maxxFt",flag)
! minimum Feynman parameter xF
xFctoplow=readvalue_r("minxFt",flag)
!PRINT *,'etatop',etatop
! maximum rapidity photon
etag=readvalue_r("maxrapp",flag)
! maximum y rapidity photon
ycg=readvalue_r("maxyrapp",flag)
! minimum y rapidity photon
ycglow=readvalue_r("minyrapp",flag)
! maximum Feynman parameter xF
xFcg=readvalue_r("maxxFp",flag)
! minimum Feynman parameter xF
xFcglow=readvalue_r("minxFp",flag)
! maximum rapidity charmonium
etaconia=readvalue_r("maxrapconia",flag)
! maximum y rapidity charmonium
ycconia=readvalue_r("maxyrapconia",flag)
! minimum y rapidity charmonium
ycconialow=readvalue_r("minyrapconia",flag)
! maximum Feynman parameter xF
xFcconia=readvalue_r("maxxFconia",flag)
! minimum Feynman parameter xF
xFcconialow=readvalue_r("minxFconia",flag)
! maximum rapidity bottomnium
etabonia=readvalue_r("maxrapbonia",flag)
! maximum y rapidity bottomnium
ycbonia=readvalue_r("maxyrapbonia",flag)
! minimum y rapidity bottomnium
ycbonialow=readvalue_r("minyrapbonia",flag)
! maximum Feynman parameter xF 
xFcbonia=readvalue_r("maxxFbonia",flag)
! minimum Feynman parameter xF 
xFcbonialow=readvalue_r("minxFbonia",flag)
! maximum rapidity Bc
etaBconia=readvalue_r("maxrapBconia",flag)
! maximum y rapidity Bc
ycBconia=readvalue_r("maxyrapBconia",flag)
! minimum y rapidity Bc
ycBconialow=readvalue_r("minyrapBconia",flag)
! maximum Feynman parameter xF
xFcBconia=readvalue_r("maxxFBconia",flag)
! minimum Feynman parameter xF
xFcBconialow=readvalue_r("minxFBconia",flag)
!PRINT *,'etag',etag
! minimum DR lepton with lepton
drll=readvalue_r("mindrll",flag)
! minimum mass lepton with lepton
gll=readvalue_r('minmll',flag)
!PRINT *,'drll',drll
! minimum DR lepton with quark
drlq=readvalue_r("mindrlq",flag)
! minimum mass lepton with quark
glq=readvalue_r('minmlq',flag)
!PRINT *,'drlq',drlq
! minimum DR quark with quark
drqq=readvalue_r("mindrqq",flag)
!PRINT *,'drqq',drqq
! minimum DR photon with fermion
drgx=readvalue_r("mindrpf",flag)
! minimum mass photon with fermion
ggx=readvalue_r("minmpf",flag)
!PRINT *,'drgx',drgx
! minimum mass quark with quark
gqq=readvalue_r("minmqqp",flag)
!PRINT *,'gqq',gqq
! minimum DR quark with b-quark
drqb=readvalue_r("mindrqb",flag)
!PRINT *,'drqb',drqb
! minimum mass quark with b-quark
gqb=readvalue_r("minmqb",flag)
!PRINT *,'gqb',gqb
! minimum DR b-quark with b-quark
drbb=readvalue_r("mindrbb",flag)
!PRINT *,'drbb',drbb
! minimum mass b-quark with b-quark
gbb=readvalue_r("minmbb",flag)
! minimum mass u,d,s quarks and gluon  with partonic beam
gbeamq=readvalue_r("minmqbeam",flag)
!PRINT *,'gbb',gbb
! minimum missing transverse momentum
ptmiss=readvalue_r("minptmiss",flag)
!PRINT *,'ptmiss',ptmiss
!END SUBROUTINE readcuts_HC
CLOSE(UNIT=udefault)
CLOSE(UNIT=uinput)

!SUBROUTINE Auto_HC_Cuts
!IMPLICIT NONE
!INTEGER::i,j,i1,j1,l,k
!REAL(KIND(1d0))::r
DO i=3,n
   IF(ifl(i).LT.13)THEN
       i1=IABS(MOD(ifl(i),4))
   ELSE
       i1=ifl(i)
   ENDIF
   IF(i1.EQ.35)i1=0
!- charm quarks
   IF(IABS(ifl(i)).EQ.7)i1=7
!- bottom quarks
   IF(IABS(ifl(i)).EQ.12)i1=12
!- top quarks
   IF(IABS(ifl(i)).EQ.11)i1=11
!       goto 111
! -------
! special for mc4lhc 2003
![
!       if(ifl(i).eq.11)i1=-1
!       if(ifl(i).eq.-11)i1=-1
!       if(ifl(i).eq.12)i1=-1
!       if(ifl(i).eq.-12)i1=-1
!]
! 111    continue
!   photons
    IF(i1.EQ.31)THEN
        ptc(i)=ptg
        IF(maxptg.GE.0d0)maxptc(i)=maxptg
        etac(i)=etag
        yycut(i)=ycg
        yycutlow(i)=ycglow
        xFcut(i)=xFcg
        xFcutlow(i)=xFcglow
    ENDIF
!   neutrinos
    IF(i1.EQ.1)THEN
        ptmissc=ptmiss
    ENDIF
!   ch. leptons
    IF(i1.EQ.2)THEN
!      if(i1.le.2)then
        ptc(i)=ptl
        IF(maxptl.GE.0d0)maxptc(i)=maxptl
        etac(i)=etal
        yycut(i)=ycl
        yycutlow(i)=ycllow
        xFcut(i)=xFcl
        xFcutlow(i)=xFcllow
    ENDIF
!   quarks
    IF(i1.EQ.3.OR.i1.EQ.0)THEN
        ptc(i)=ptq
        IF(maxptq.GE.0d0)maxptc(i)=maxptq
        etac(i)=etaq
        yycut(i)=ycq
        yycutlow(i)=ycqlow
        xFcut(i)=xFcq
        xFcutlow(i)=xFcqlow
        gbeammass(1,i)=gbeamq
        gbeammass(2,i)=gbeamq
    ENDIF
!  charm
    IF(i1.EQ.7)THEN
!	ptc(i)=ptcharm
        IF(parton2hadrontype(i).EQ.0)THEN
           ! open charm
           ptc(i)=ptcharm
           IF(maxptcharm.GE.0d0)maxptc(i)=maxptcharm
           etac(i)=etacharm
           yycut(i)=yccharm
           yycutlow(i)=yccharmlow
           xFcut(i)=xFccharm
           xFcutlow(i)=xFccharmlow
        ELSEIF(parton2hadrontype(i).EQ.1)THEN
           ! charmonium
           ptc(i)=ptconia/2d0
           IF(maxptconia.GE.0d0)maxptc(i)=maxptconia/2d0
           etac(i)=etaconia
           yycut(i)=ycconia
           yycutlow(i)=ycconialow
           xFcut(i)=xFcconia
           xFcutlow(i)=xFcconialow
        ELSEIF(parton2hadrontype(i).EQ.3)THEN
           ! Bc
           ptc(i)=ptBconia*rcm/(rcm+rbm)
           IF(maxptBconia.GE.0d0)maxptc(i)=maxptBconia*rcm/(rcm+rbm)
           etac(i)=etaBconia
           yycut(i)=ycBconia
           yycutlow(i)=ycBconialow
           xFcut(i)=xFcBconia
           xFcutlow(i)=xFcBconialow
        ELSE
           PRINT *,"Wrong with the hadron type for parton ",i
           STOP
        ENDIF
    ENDIF
!  bottom 
    IF(i1.EQ.12)THEN
!        ptc(i)=ptbottom
        IF(parton2hadrontype(i).EQ.0)THEN
           ! open bottom
           ptc(i)=ptbottom
           IF(maxptbottom.GE.0d0)maxptc(i)=maxptbottom
           etac(i)=etabottom
           yycut(i)=ycbottom
           yycutlow(i)=ycbottomlow
           xFcut(i)=xFcbottom
           xFcutlow(i)=xFcbottomlow
        ELSEIF(parton2hadrontype(i).EQ.2)THEN
           ! bottonium
           ptc(i)=ptbonia/2d0
           IF(maxptbonia.GE.0d0)maxptc(i)=maxptbonia/2d0
           etac(i)=etabonia
           yycut(i)=ycbonia
           yycutlow(i)=ycbonialow
           xFcut(i)=xFcbonia
           xFcutlow(i)=xFcbonialow
        ELSEIF(parton2hadrontype(i).EQ.3)THEN
           ! Bc
           ptc(i)=ptBconia*rbm/(rcm+rbm)
           IF(maxptBconia.GE.0d0)maxptc(i)=maxptBconia*rbm/(rcm+rbm)
           etac(i)=etaBconia
           yycut(i)=ycBconia
           yycutlow(i)=ycBconialow
           xFcut(i)=xFcBconia
           xFcutlow(i)=xFcBconialow
        ELSE
           PRINT *, "Wrong with the hadron type for parton ",i
           STOP
        ENDIF
    ENDIF
!   top
    IF(i1.EQ.11)THEN
        ptc(i)=pttop
        IF(maxpttop.GE.0d0)maxptc(i)=maxpttop
        etac(i)=etatop
        yycut(i)=yctop
        yycutlow(i)=yctoplow
        xFcut(i)=xFctop
        xFcutlow(i)=xFctoplow
    ENDIF
ENDDO
    
DO i=3,n-1
    IF(ifl(i).LT.13)THEN
        i1=IABS(MOD(ifl(i),4))
    ELSE
        i1=ifl(i)
    ENDIF
    IF(i1.EQ.35)i1=0
!- bottom quarks                                           
    IF(IABS(ifl(i)).EQ.12)i1=12
!- top quarks                                                                                                               
    IF(IABS(ifl(i)).EQ.11)i1=11
!- quarkonium
    IF(parton2hadrontype(i).NE.0)i1=1
!        goto 222
! -------
! special for mc4lhc 2003
![
!       if(ifl(i).eq.11)i1=-1
!       if(ifl(i).eq.-11)i1=-1
!       if(ifl(i).eq.12)i1=-1
!       if(ifl(i).eq.-12)i1=-1
!]

! 222    continue
    DO j=i+1,n
       IF(ifl(j).LT.13)THEN
          j1=IABS(MOD(ifl(j),4))
       ELSE
          j1=ifl(j)
       ENDIF
       IF(j1.EQ.35)j1=0
!- bottom quarks                                                 
                                           
       IF(IABS(ifl(j)).EQ.12)j1=12
!- top quarks                                                         
                                            
       IF(IABS(ifl(j)).EQ.11)j1=11
!- quarkonium
       IF(parton2hadrontype(j).NE.0)j1=1
!         goto 333

! -------
! special for mc4lhc 2003
![
!       if(ifl(j).eq.11)j1=-1
!       if(ifl(j).eq.-11)j1=-1
!       if(ifl(j).eq.12)j1=-1
!       if(ifl(j).eq.-12)j1=-1
!]
! 333    continue
!  n-n
        IF(i1.EQ.1.OR.j1.EQ.1)THEN
!  l-q
        ELSEIF(i1.EQ.2.AND.(j1.EQ.3.OR.j1.EQ.0.OR.j1.EQ.12))THEN
             drc(i,j)=drlq
             gmas(i,j)=MAX(glq,gmas(i,j))
!  l-l
        ELSEIF(i1.EQ.2.AND.j1.EQ.2)THEN
             drc(i,j)=drll
             gmas(i,j)=MAX(gll,gmas(i,j))
!  q-l
        ELSEIF(j1.EQ.2.AND.(i1.EQ.3.OR.i1.EQ.0.OR.i1.EQ.12))THEN
             drc(i,j)=drlq
             gmas(i,j)=MAX(glq,gmas(i,j))
!  g-x  photon-charged
        ELSEIF(i1.EQ.31.AND.j1.NE.1)THEN
             drc(i,j)=drgx
             gmas(i,j)=MAX(ggx,gmas(i,j))
!  x-g  charged-photon
        ELSEIF(j1.EQ.31.AND.i1.NE.1)THEN
             drc(i,j)=drgx
             gmas(i,j)=MAX(ggx,gmas(i,j))
!  q-q
        ELSEIF((i1.EQ.3.OR.i1.EQ.0).AND.(j1.EQ.3.OR.j1.EQ.0))THEN
             drc(i,j)=drqq
             gmas(i,j)=MAX(gqq,gmas(i,j))
!  q-b
        ELSEIF((i1.EQ.3.OR.i1.EQ.0).AND.(j1.EQ.12))THEN
             drc(i,j)=drqb
             gmas(i,j)=MAX(gqb,gmas(i,j))
!  b-q
        ELSEIF((j1.EQ.3.OR.j1.EQ.0).AND.(i1.EQ.12))THEN
             drc(i,j)=drqb
             gmas(i,j)=MAX(gqb,gmas(i,j))
!  b-b
        ELSEIF((i1.EQ.12).AND.(j1.EQ.12))THEN
             drc(i,j)=drbb
             gmas(i,j)=MAX(gbb,gmas(i,j))
        ENDIF
    ENDDO
ENDDO
 
DO i=3,n
    IF(ptc(i).GT.0)THEN
         ec(i)=MAX(ec(i),ptc(i))
    ELSE
         ec(i)=MAX(ec(i),parmas(ifl(i)))
    ENDIF
    r=EXP(2*etac(i))
    c1(i)=(r-1)/(r+1)
    c2(i)=c1(i)
ENDDO
   
!      do i=3,n-1
!       do j=i+1,n
!        cc(i,j)=cos(drc(i,j))
!       enddo
!      enddo

DO i=3,n-1
    DO j=i+1,n
        gmas(i,j)=MAX(gmas(i,j),DSQRT( 2*ptc(i)*ptc(j)*(1-COS(drc(i,j))) ) )
    ENDDO
ENDDO

DO l=4,20
    DO k=3,l-1
        drc(l,k)= drc(k,l)
        cc(l,k)=  cc(k,l)
        gmas(l,k)=  gmas(k,l)
    ENDDO
ENDDO
  
WRITE(*,*)'---------------------------------------------------'
WRITE(*,*)'        the cuts        '
DO i=3,n
    WRITE(*,*)'pt     of  ',i,'   particle   ',ptc(i)
    IF(maxptc(i).GE.0d0)THEN
       WRITE(*,*)'max pt of  ',i,'   particle   ',maxptc(i)
       IF(maxptc(i).LE.ptc(i))THEN
          WRITE(*,*)"ERROR: One final state ",i," was cut off by pt cut"
          STOP
       ENDIF
    ELSE
       WRITE(*,*)'no max pt cut of  ',i,'   particle   '
    ENDIF
    WRITE(*,*)'energy of  ',i,'   particle   ',ec(i)
ENDDO
DO i=3,n
    WRITE(*,*)'rapidity of  ',i,'   particle   ',etac(i)
ENDDO
DO i=3,n
   WRITE(*,*)'max y rapidity of ',i,'   particle   ',yycut(i)
   WRITE(*,*)'min y rapidity of ',i,'   particle   ',yycutlow(i)
   IF(yycut(i).LE.yycutlow(i))THEN
      WRITE(*,*)"ERROR:One final state ",i," was cut off by y rapidity"
      STOP
   ENDIF
   WRITE(*,*)'max Feynman parameter xF of ',i,' particle ',xFcut(i)
   IF(xFcut(i).LT.1d0)xFcutflag=.TRUE.
   WRITE(*,*)'min Feynman parameter xF of ',i,' particle ',xFcutlow(i)
   IF(xFcutlow(i).GT.-1d0)xFcutflag=.TRUE.
   IF(xFcut(i).LE.xFcutlow(i))THEN
      WRITE(*,*)"ERROR:One final state ",i," was cut off by xF "
      STOP
   ENDIF
ENDDO
WRITE(*,*)'The maxrapidity of the first particle',y1cup
WRITE(*,*)'The minrapidity of the first particle',y1clow
DO i=3,n
    WRITE(*,*)'cos-beam1 of ',i,'   particle   ',c1(i)
ENDDO
DO i=3,n
    WRITE(*,*)'cos-beam2 of ',i,'   particle   ',c2(i)
ENDDO

DO i=3,n-1
    DO j=i+1,n
       IF(i.NE.j)THEN
          WRITE(*,*)'DR     ',i,'  with  ',j,drc(i,j)
          WRITE(*,*)'cos of ',i,'  with  ',j,cc(i,j)
       ENDIF
    ENDDO
ENDDO

DO i=3,n-1
    DO j=i+1,n
       IF(i.NE.j)WRITE(*,*)'mass of ',i,'  with  ',j,gmas(i,j)
    ENDDO
ENDDO
WRITE(*,*)'minimum missing transverse momentum',ptmissc
DO i=1,2
   DO j=3,n
      WRITE(*,*)'mass of partonic beam ',i,' with final ',j,gbeammass(i,j)
   ENDDO
ENDDO
WRITE(*,*)'---------------------------------------------------'
!END SUBROUTINE Auto_HC_Cuts
END SUBROUTINE readcuts_HC     


SUBROUTINE cuts_HC(icut)
IMPLICIT NONE
REAL(KIND=DBL),DIMENSION(20,4)::p
REAL(KIND=DBL),DIMENSION(4)::pboo,ponia,pboo2
REAL(KIND=DBL)::pi,q,e,pt,eta,d1,d2,dr,pmiss_E,pmiss_x,pmiss_y,pmiss_z,pt_miss,s,t
REAL(KIND=DBL)::exp1,exp2,aaa,bbb,ptcut
INTEGER,INTENT(OUT)::icut
INTEGER::i,l,flag,l1,l2,i_nu,i1
!    -----
!     special for mc4lhc
!    -----
pi=DACOS(-1d0)
DO i=3,n
   p(i,1:4)=phegas_pmom(i,1:4)
ENDDO
icut=0
! invariant mass cuts
DO l1=3,n-1
   DO l2=l1+1,n
      IF(l1.EQ.l2.OR.Quarkonium3(l1).EQ.Quarkonium3(l2).AND.Quarkonium3(l1).NE.0)CYCLE
      s=2*scalar_product(phegas_pmom(l1,1:4),phegas_pmom(l2,1:4))&
         +scalar_product(phegas_pmom(l1,1:4),phegas_pmom(l1,1:4))&
         +scalar_product(phegas_pmom(l2,1:4),phegas_pmom(l2,1:4))
      IF(s.LT.gmas(l1,l2)**2)RETURN
   ENDDO
ENDDO
! special for the beam cuts for NLO* or NNLO*
DO l1=1,2
   DO l2=3,n
      t=scalar_product(phegas_pmom(l1,1:4),phegas_pmom(l1,1:4))&
           +scalar_product(phegas_pmom(l2,1:4),phegas_pmom(l2,1:4))&
           -2*scalar_product(phegas_pmom(l1,1:4),phegas_pmom(l2,1:4))
      IF(DABS(t).LT.gbeammass(l1,l2)**2)RETURN
   ENDDO
ENDDO

IF((istruc.EQ.1.OR..NOT.labeqcoll).AND.imode.EQ.0.AND.NDecayChains.EQ.0)THEN ! why imode.EQ.0 ? because it has been boosted into lab frame when imode.NE.0 or with decay chains (spin correlations)
	q=ehat
	e=q/SQRT(xp1*xp2)
        exp1=ebeam(1)*xp1
        exp2=ebeam(2)*xp2
        IF(fixtarget)THEN
           IF(.NOT.fixtargetrev)THEN
              exp1=exp1+FT_M*xp1/2d0
              exp2=exp2/2d0
           ELSE
              exp1=exp1/2d0
              exp2=exp2+FT_M*xp2/2d0
           ENDIF
        ENDIF
	pboo(4)=(exp1+exp2)
	pboo(3)=(exp1-exp2)
	pboo(1:2)=0
! boost to the lab frame
	DO i=3,n
		CALL boostl(q,pboo,p(i,1:4))
	ENDDO
ENDIF
flag=0
DO l=3,n
   pt=SQRT(p(l,1)**2+p(l,2)**2)
   IF(pt.LT.ptc(l).AND.nhad.GT.3)THEN
	  flag=1
	  EXIT
   ENDIF
   IF(maxptc(l).GE.0d0.AND.nhad.GT.3)THEN
      IF(pt.GT.maxptc(l))THEN
         flag=1
         EXIT
      ENDIF
   ENDIF
   eta=prapidity(p(l,1:4))
   IF(ABS(eta).GT.etac(l).AND.nhad.GT.3)THEN 
      flag=1
      EXIT
   ENDIF
   eta=rapidity(p(l,1:4))
   IF(absrap)eta=ABS(eta)
   IF(eta.GT.yycut(l).OR.eta.LT.yycutlow(l))THEN
      flag=1
      EXIT
   ENDIF
   ! special for the first particle , which can be used to calculate the y distribution
   IF(l.EQ.3)THEN
!		eta=rapidity(p(l,1:4))
      IF(eta.GT.y1cup.OR.eta.LT.y1clow)THEN
	flag=1
	EXIT
      ENDIF
   ENDIF
ENDDO

IF(flag.EQ.0)THEN
   DO l1=3,n
      DO l2=3,n
         IF(l1.EQ.l2.OR.Quarkonium3(l1).EQ.Quarkonium3(l2).AND.Quarkonium3(l1).NE.0)CYCLE
         d1=prapidity(p(l1,1:4))-prapidity(p(l2,1:4))
         d2=ph4(p(l1,1),p(l1,2),p(l1,3))-ph4(p(l2,1),p(l2,2),p(l2,3))
         d2=MIN(DABS(d2),2*pi-DABS(d2))
         IF(d2/pi.GT.1.d0)WRITE(*,*)d2/pi
         dr=SQRT(d1**2+d2**2)
         IF(dr.LT.drc(l1,l2))THEN 
			 flag=1
			 EXIT
         ENDIF
      ENDDO
	  IF(flag.EQ.1)EXIT
   ENDDO
   IF(flag.EQ.0)THEN
      i_nu=0
      pmiss_E=0.d0
      pmiss_x=0.d0
      pmiss_y=0.d0
      pmiss_z=0.d0
      DO l=3,n
         IF(ifl(l).LT.13)THEN
            i1=IABS(MOD(ifl(l),4))
         ENDIF
         IF(i1.EQ.1)THEN        !neutrino found
            i_nu=1
            pmiss_E = pmiss_E + p(l,4)
            pmiss_x = pmiss_x + p(l,1)
            pmiss_y = pmiss_y + p(l,2)
            pmiss_z = pmiss_z + p(l,3)
         ENDIF
      ENDDO
      pt_miss=SQRT((pmiss_x)**2+(pmiss_y)**2)
      IF(i_nu.NE.1 .OR. pt_miss.GE.ptmissc.OR.nhad.EQ.3)THEN
         icut=1
	  ENDIF
	ENDIF
ENDIF
! special for the beam cuts for NLO* or NNLO*
! special cutoffs in the literature
IF(icut.EQ.1.AND.flag.EQ.0)THEN
   IF(literature_cutoffs.EQ.14060484)THEN
      l1=3
      DO l=3,nhad
         IF(iflh(l).GT.440001.AND.iflh(l).LT.449999)THEN
            ! charmonia
            ponia(1:4)=p(l1,1:4)+p(l1+1,1:4)
            eta=rapidity(ponia(1:4))
            IF(ABS(eta).LT.1.2d0)THEN
               pt=SQRT(ponia(1)**2+ponia(2)**2)
               IF(pt.LE.6.5d0)THEN
                  icut=0
                  EXIT
               ENDIF
            ELSEIF(ABS(eta).GT.1.2d0.AND.ABS(eta).LT.1.43d0)THEN
               pt=SQRT(ponia(1)**2+ponia(2)**2)
               aaa=-8.695652173913045d0
               bbb=16.934782608695656d0
               ptcut=aaa*ABS(eta)+bbb
               IF(pt.LE.ptcut)THEN
                  icut=0
                  EXIT
               ENDIF
            ELSEIF(ABS(eta).GT.1.43d0.AND.ABS(eta).LT.2.2d0)THEN
               pt=SQRT(ponia(1)**2+ponia(2)**2)
               IF(pt.LE.4.5d0)THEN
                  icut=0
                  EXIT
               ENDIF
            ELSE
               icut=0
               EXIT
            ENDIF
            l1=l1+2
         ELSEIF(ABS(iflh(l)).GT.100)THEN
            l1=l1+2
         ELSE
            l1=l1+1
         ENDIF
      ENDDO
   ENDIF
ENDIF
! 1000 continue

!     ---
!     boost from lab to the collision frame for xF cut
!     ---
IF(xFcutflag.AND.icut.EQ.1.AND.flag.EQ.0)THEN
   q=ehat
   e=q/SQRT(xp1*xp2)
   IF(.NOT.labeqcoll)THEN
      IF(.NOT.fixtarget)THEN
         pboo2(4)=(ebeam(1)+ebeam(2))
         pboo2(3)=-(ebeam(1)-ebeam(2))
      ELSE
         IF(.NOT.fixtargetrev)THEN
            pboo2(4)=(ebeam(1)+ebeam(2))
            pboo2(3)=-ebeam(1)
         ELSE
            pboo2(4)=(ebeam(1)+ebeam(2))
            pboo2(3)=ebeam(2)
         ENDIF
      ENDIF
      pboo2(1:2)=0
      ! boost from the lab frame to the collision frame
      DO i=3,n
         CALL boostl(e,pboo2,p(i,1:4))
      ENDDO
   ENDIF
   DO l=3,n
      eta=xFeynman(p(l,1:4),e)
      IF(eta.GT.xFcut(l).OR.eta.LT.xFcutlow(l))THEN
         icut=0
         EXIT
      ENDIF
   ENDDO
   IF(.NOT.labeqcoll)THEN
      pboo2(3)=-pboo2(3)
      ! boost back to the lab frame
      DO i=3,n
         CALL boostl(e,pboo2,p(i,1:4))
      ENDDO
   ENDIF
ENDIF

!     ---
!     boost to the cms
!     ---
IF((istruc.EQ.1.OR..NOT.labeqcoll).AND.imode.EQ.0.AND.NDecayChains.EQ.0)THEN ! why imode=0? because it has been boosted into lab frame when imode.NE.0 or withdecay chains (spin correlations)
    pboo(3)=-pboo(3)
    DO i=3,n
       CALL boostl(q,pboo,p(i,1:4))
    ENDDO
ENDIF
END SUBROUTINE cuts_HC

SUBROUTINE readcuts_epem
IMPLICIT NONE
CHARACTER(len=24)::file
REAL(KIND(1d0))::cutoff,pi
INTEGER::i,j,i1,j1,l,k
REAL(KIND(1d0))::r
REAL(KIND(1d0))::el,eq,eg,cl,cq,cg,cll,clq,cqq,cgx,gqq,ptmiss
INTEGER,PARAMETER::udefault=45138,uinput=45139
INTEGER::iounit,flag=0
LOGICAL::lexist
! open default input file
INQUIRE(FILE=TRIM(input_dir)//"default.inp",EXIST=lexist)
IF(.NOT.lexist)THEN
	PRINT *,"Warning: the file default.inp does not exist ! STOP !"
	STOP
ENDIF
INQUIRE(FILE=TRIM(input_dir)//"default.inp",OPENED=lexist)
IF(lexist)THEN
    INQUIRE(FILE=TRIM(input_dir)//"default.inp",NUMBER=iounit)
	IF(iounit.NE.udefault)THEN
		PRINT *,"WARNING: the default.inp has been linked with another unit ! Close and reopen !"
		CLOSE(UNIT=iounit)
		OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
	ENDIF
ELSE
	OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
ENDIF
! open user's input file
IF(TRIM(Input_File)/="default.inp")THEN
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
	IF(.NOT.lexist)THEN
		PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
		STOP
	ENDIF
	INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
	IF(lexist)THEN
		INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
		IF(iounit.NE.uinput)THEN
			PRINT *,"WARNING: the "//TRIM(Input_File)//" has been linked with another unit ! Close and reopen !"
			CLOSE(UNIT=iounit)
			OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
		ENDIF
	ELSE
		OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
	ENDIF
ELSE
    flag=1
ENDIF

cutoff=readvalue_r("cutoffe",flag)

PRINT *,'WARNING: CUTOFF SET AT',cutoff
ec(3:20)=0
c1(3:20)=1
c2(3:20)=1
cc(3:20,3:20)=1
gmas(3:20,3:20)=cutoff
pi=DACOS(-1.d0)

! minimum lepton energy
el=readvalue_r("minenl",flag)
! minimum quark  energy
eq=readvalue_r("minenq",flag)
! minimum photon energy
eg=readvalue_r("minenp",flag)
! maximum cos lepton with beam
cl=COS(readvalue_r("minanglb",flag)/180.d0*pi)
! maximum cos quark with beam
cq=COS(readvalue_r("minangqb",flag)/180.d0*pi)
! maximum cos photon with beam
cg=COS(readvalue_r("minangpb",flag)/180.d0*pi)
! maximum cos lepton with lepton
cll=COS(readvalue_r("minangll",flag)/180.d0*pi)
! maximum cos lepton with quark
clq=COS(readvalue_r("minanglq",flag)/180.d0*pi)
! maximum cos quark with quark
cqq=COS(readvalue_r("minangqq",flag)/180.d0*pi)
! maximum cos photon with fermions
cgx=COS(readvalue_r("minangpf",flag)/180.d0*pi)
! minimum mass quark with quark
gqq=readvalue_r("minmqqe",flag)
! minimum missing transverse momentum
ptmiss=readvalue_r("minptmiss",flag)

!END SUBROUTINE readcuts_epem
CLOSE(UNIT=udefault)
CLOSE(UNIT=uinput)

!SUBROUTINE Auto_epem_cuts
!IMPLICIT NONE
!INTEGER::i,i1,j,j1,l,k
DO i=3,n
   IF(ifl(i).LT.13)THEN
      i1=IABS(MOD(ifl(i),4))
   ELSE
      i1=ifl(i)
   ENDIF
   IF(i1.EQ.35)i1=0
   ec(i)=parmas(ifl(i))
!   photons
   IF(i1.EQ.31)THEN
        ec(i)=eg
        c1(i)=cg
        c2(i)=cg
    ENDIF
!   neutrinos
    IF(i1.EQ.1)THEN
        ptmissc=ptmiss
    ENDIF
!   ch. leptons
    IF(i1.EQ.2)THEN
        ec(i)=MAX(el,ec(i))
        c1(i)=cl
        c2(i)=cl
    ENDIF
!   quarks
    IF(i1.EQ.3.OR.i1.EQ.0)THEN
        ec(i)=MAX(eq,ec(i))
        c1(i)=cq
        c2(i)=cq
    ENDIF
ENDDO
    
DO i=3,n-1
    IF(ifl(i).LT.13)THEN
        i1=IABS(MOD(ifl(i),4))
    ELSE
        i1=ifl(i)
    ENDIF
    IF(i1.EQ.35)i1=0
    DO j=i+1,n
       IF(ifl(j).LT.13)THEN
           j1=IABS(MOD(ifl(j),4))
       ELSE
           j1=ifl(j)
       ENDIF
       IF(j1.EQ.35)j1=0
!  n-n
       IF(i1.EQ.1.OR.j1.EQ.1)THEN
!  l-q
       ELSEIF(i1.EQ.2.AND.(j1.EQ.3.OR.j1.EQ.0))THEN
           cc(i,j)=clq
!  l-l
       ELSEIF(i1.EQ.2.AND.j1.EQ.2)THEN
           cc(i,j)=cll
!  q-l
       ELSEIF(j1.EQ.2.AND.(i1.EQ.3.OR.i1.EQ.0))THEN
           cc(i,j)=clq
!  g-x  photon-charged
       ELSEIF(i1.EQ.31.AND.j1.NE.1)THEN
           cc(i,j)=cgx
!  x-g  charged-photon
       ELSEIF(j1.EQ.31.AND.i1.NE.1)THEN
           cc(i,j)=cgx
!  q-q
       ELSEIF((i1.EQ.3.OR.i1.EQ.0).AND.(j1.EQ.3.OR.j1.EQ.0))THEN
           cc(i,j)=cqq
           gmas(i,j)=gqq
!
       ENDIF
   ENDDO
ENDDO

DO l=4,20
   DO k=3,l-1
       gmas(l,k)=gmas(k,l)
       cc(l,k)=  cc(k,l)
   ENDDO
ENDDO

WRITE(*,*)'---------------------------------------------------'
WRITE(*,*)'        the cuts        '
DO i=3,n
   WRITE(*,*)'energy of  ',i,ifl(i),'   particle   ',ec(i)
ENDDO
DO i=3,n
   WRITE(*,*)'cos-beam1 of ',i,ifl(i),'   particle   ',c1(i)
ENDDO
DO i=3,n
   WRITE(*,*)'cos-beam2 of ',i,ifl(i),'   particle   ',c2(i)
ENDDO

DO i=3,n-1
   DO j=i+1,n
      IF(i.NE.j)WRITE(*,*)'cos of ',i,ifl(i),'  with  ',j,ifl(j),'  ',cc(i,j)
   ENDDO
ENDDO

DO i=3,n-1
   DO j=i+1,n
      IF(i.NE.j)WRITE(*,*)'mass of ',i,ifl(i),'  with  ',j,ifl(j),'  ',gmas(i,j)
   ENDDO
ENDDO
WRITE(*,*)'minimum missing transverse momentum',ptmissc
WRITE(*,*)'---------------------------------------------------'

!END SUBROUTINE Auto_epem_cuts
END SUBROUTINE readcuts_epem
      
SUBROUTINE cuts_epem(icut)
IMPLICIT NONE
INTEGER,INTENT(OUT)::icut
REAL(KIND=DBL),DIMENSION(20,4)::p
INTEGER::l,l1,l2,i_nu,i1
REAL(KIND=DBL)::c,s,pmiss_E,pmiss_x,pmiss_y,pmiss_z,pt_miss
REAL(KIND=DBL),DIMENSION(4)::pboo
REAL(KIND=DBL)::exp1,exp2,q
INTEGER::i
DO i=3,n
   p(i,1:4)=phegas_pmom(i,1:4)
ENDDO
! boost to lab frame when the collision frame is not equal to it 
IF((.NOT.labeqcoll).AND.imode.EQ.0.AND.NDecayChains.EQ.0)THEN ! why imode.EQ.0 ? because it has been boosted into lab frame when imode.NE.0 or with decay chains (spin correlations)
        q=ehat
        exp1=ebeam(1)
        exp2=ebeam(2)
        pboo(4)=(exp1+exp2)
        pboo(3)=(exp1-exp2)
        pboo(1:2)=0
! boost to the lab frame
        DO i=3,n
           CALL boostl(q,pboo,p(i,1:4))
        ENDDO
ENDIF  
icut=0
DO l=3,n 
   IF(p(l,4).LT.ec(l))RETURN
   c=p(l,3)/p(l,4)
   IF( c.GT.c1(l))RETURN
   IF(-c.GT.c2(l))RETURN
ENDDO
       
DO l1=3,n
   DO l2=3,n
      IF(l1.EQ.l2.OR.Quarkonium3(l1).EQ.Quarkonium3(l2).AND.Quarkonium3(l1).NE.0)CYCLE
      IF(cosij(p(l1,1:3),p(l2,1:3)).GT.cc(l1,l2))RETURN
      s=2*scalar_product(p(l1,1:4),p(l2,1:4))&
         +scalar_product(p(l1,1:4),p(l1,1:4))&
         +scalar_product(p(l2,1:4),p(l2,1:4))
      IF(s.LT.gmas(l1,l2)**2)RETURN
   ENDDO
ENDDO

i_nu=0
pmiss_E=0.d0
pmiss_x=0.d0
pmiss_y=0.d0
pmiss_z=0.d0
DO l=3,n
   IF(ifl(l).LT.13)THEN
      i1=IABS(MOD(ifl(l),4))
   ENDIF
   IF(i1.EQ.1)THEN       !neutrino found
      i_nu=1
      pmiss_E = pmiss_E + p(l,4)
      pmiss_x = pmiss_x + p(l,1)
      pmiss_y = pmiss_y + p(l,2)
      pmiss_z = pmiss_z + p(l,3)
   ENDIF
ENDDO
pt_miss=DSQRT((pmiss_x)**2+(pmiss_y)**2)
IF(i_nu.EQ.1 .AND. pt_miss.LT.ptmissc.OR.nhad.EQ.3)RETURN
icut=1
IF((.NOT.labeqcoll).AND.imode.EQ.0.AND.NDecayChains.EQ.0)THEN ! why imode=0? because it has been boosted into lab frame when imode.NE.0 or with decay chains (spin correlations)
    pboo(3)=-pboo(3)
    DO i=3,n
       CALL boostl(q,pboo,p(i,1:4))
    ENDDO
ENDIF
END SUBROUTINE cuts_epem
END MODULE Cuts_Module
