MODULE pp_aajj_DPS_cuts
  USE Helac_Global
  USE pp_aajj_dps_global
  USE Constants
  USE Kinetic_Func
  USE MC_VEGAS
  IMPLICIT NONE
CONTAINS
  SUBROUTINE readcuts_pp_aajj_dps
    IMPLICIT NONE
    CHARACTER(len=24)::file
    LOGICAL::lexist
    INTEGER::iounit,flag=0,i,i1,j,j1
    REAL(KIND(1d0))::ptq,ptg,etaq,ycq,ycqlow,etag,ycg,ycglow
    REAL(KIND(1d0))::cutoff,xFcq,xFcqlow,xFcg,xFcglow
    REAL(KIND(1d0))::gbeamq,drqq,gqq
    INTEGER::nnhad
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
          PRINT *,"WARNING: the default.inp has been linked with another unit! Close and reopen !"
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
      drc(3:20,3:20)=0d0
      etac(3:20)=20
      yycut(3:20)=1d9
      IF(absrap)THEN
         yycutlow(3:20)=0d0
      ELSE
         yycutlow(3:20)=-1d9
      ENDIF
      xFcut(3:20)=1d0
      xFcutlow(3:20)=1d0
      xFcutflag=.FALSE.
      ec(3:20)=cutoff
      c1(3:20)=1d0
      c2(3:20)=1d0
      cc(3:20,3:20)=1d0
      gmas(3:20,3:20)=cutoff
      gbeammass(1:2,3:20)=0d0
      y1cup=30d0
      y1clow=0d0
      y1cup=readvalue_r("maxy1c",flag)
      y1clow=readvalue_r("miny1c",flag)
      ! minimum quark pt (5-flavor quarks and gluon)
      ptq=readvalue_r("minptq",flag)
      ! minimum photon pt
      ptg=readvalue_r("minptp",flag)
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
      ! minimum mass quark with quark
      gqq=readvalue_r("minmqqp",flag)
      ! minimum mass u,d,s quarks and gluon  with partonic beam
      gbeamq=readvalue_r("minmqbeam",flag)
      ! minimum dr quark with quark
      drqq=readvalue_r("mindrqq",flag)
      CLOSE(UNIT=udefault)
      CLOSE(UNIT=uinput)
      IF(dpsorsps.EQ.1)THEN
         nnhad=nhad
      ELSE
         nnhad=2*nhad
      ENDIF
      DO i=5,nnhad
         ! the first two final states are photon
         ! while the last two final states are jets
         IF(i.GT.6)THEN
            ! jets
            ptc(i)=ptq
            etac(i)=etaq
            yycut(i)=ycq
            yycutlow(i)=ycqlow
            xFcut(i)=xFcq
            xFcutlow(i)=xFcqlow
         ELSE
            ! photon
            ptc(i)=ptg
            etac(i)=etag
            yycut(i)=ycg
            yycutlow(i)=ycglow
            xFcut(i)=xFcg
            xFcutlow(i)=xFcglow
         ENDIF
      ENDDO
      DO i=5,nnhad
         DO j=i+1,nnhad
            IF(i.EQ.7.AND.j.EQ.8)THEN
               ! dijets
               drc(i,j)=drqq
               gmas(i,j)=MAX(gqq,gmas(i,j))
            ENDIF
         ENDDO
      ENDDO
      DO i=5,nnhad-1
         DO j=i+1,nnhad
            gmas(i,j)=MAX(gmas(i,j),DSQRT(2*ptc(i)*ptc(j)*(1-COS(drc(i,j)))))
         ENDDO
      ENDDO
      DO i=6,nnhad
         DO j=5,i-1
            drc(i,j)=drc(j,i)
            gmas(i,j)=gmas(j,i)
         ENDDO
      ENDDO
      WRITE(*,*)'---------------------------------------------------'
      WRITE(*,*)'    the cuts for p p > diphoton+dijet + X '
      WRITE(*,*)'    with double parton scattering (DPS) '
      DO i=5,nnhad
         WRITE(*,*)'pt     of  ',i,'   particle   ',ptc(i)
         WRITE(*,*)'energy of  ',i,'   particle   ',ec(i)
         WRITE(*,*)'rapidity of  ',i,'   particle   ',etac(i)
         WRITE(*,*)'max y rapidity of ',i,'   particle   ',yycut(i)
         WRITE(*,*)'min y rapidity of ',i,'   particle   ',yycutlow(i)
         WRITE(*,*)'max Feynman parameter xF of ',i,' particle ',xFcut(i)
         IF(xFcut(i).LT.1d0)xFcutflag=.TRUE.
         WRITE(*,*)'min Feynman parameter xF of ',i,' particle ',xFcutlow(i)
         IF(xFcutlow(i).GT.-1d0)xFcutflag=.TRUE.
      ENDDO
      WRITE(*,*)'The maxrapidity of the first particle',y1cup
      WRITE(*,*)'The minrapidity of the first particle',y1clow
      DO i=5,nnhad-1
         DO j=i+1,nnhad
            WRITE(*,*)'DR     ',i,'  with  ',j,drc(i,j)
            WRITE(*,*)'mass of ',i,'  with  ',j,gmas(i,j)
         ENDDO
      ENDDO
      WRITE(*,*)'---------------------------------------------------'
    END SUBROUTINE readcuts_pp_aajj_dps

    SUBROUTINE Cuts_pp_aajj_dps(icut)
      IMPLICIT NONE
      INTEGER,INTENT(OUT)::icut
      INTEGER::l,l1,l2,flag
      REAL(KIND(1d0))::s,d1,d2,dr,pt,eta,aaa,bbb,ptcut
      REAL(KIND(1d0)),DIMENSION(4)::ponia,pboo2
      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
      REAL(KIND(1d0))::e,q
      !INCLUDE "pp_psipsi_dps.inc"
      icut=0
      ! invariant mass cuts
      DO l1=5,nhad-1
         DO l2=l1+1,nhad
            s=2*scalar_product(dps_hadron_pmom(l1,1:4),dps_hadron_pmom(l2,1:4))&
                 +scalar_product(dps_hadron_pmom(l1,1:4),dps_hadron_pmom(l1,1:4))&
                 +scalar_product(dps_hadron_pmom(l2,1:4),dps_hadron_pmom(l2,1:4))
            IF(s.LT.gmas(l1,l2)**2)RETURN
         ENDDO
      ENDDO
      flag=0
      DO l=5,nhad
         pt=DSQRT(dps_hadron_pmom(l,1)**2+dps_hadron_pmom(l,2)**2)
         IF(pt.LT.ptc(l))THEN
            flag=1
            EXIT
         ENDIF
         eta=prapidity(dps_hadron_pmom(l,1:4))
         IF(ABS(eta).GT.etac(l))THEN
            flag=1
            EXIT
         ENDIF
         eta=rapidity(dps_hadron_pmom(l,1:4))
         IF(absrap)eta=ABS(eta)
         IF(eta.GT.yycut(l).OR.eta.LT.yycutlow(l))THEN
            flag=1
            EXIT
         ENDIF
      ENDDO
      IF(flag.EQ.0)THEN
         DO l1=5,nhad-1
            DO l2=l1+1,nhad
               d1=prapidity(dps_hadron_pmom(l1,1:4))-prapidity(dps_hadron_pmom(l2,1:4))
               d2=ph4(dps_hadron_pmom(l1,1),dps_hadron_pmom(l1,2),dps_hadron_pmom(l1,3))&
                    -ph4(dps_hadron_pmom(l2,1),dps_hadron_pmom(l2,2),dps_hadron_pmom(l2,3))
               d2=MIN(DABS(d2),2*pi-DABS(d2))
               IF(d2/pi.GT.1d0)WRITE(*,*)d2/pi
               dr=SQRT(d1**2+d2**2)
               IF(dr.LT.drc(l1,l2))THEN
                  flag=1
                  EXIT
               ENDIF
            ENDDO
            IF(flag.EQ.1)EXIT
         ENDDO
         IF(flag.EQ.0)icut=1
      ENDIF
      ! special cutoffs in the literature
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
            DO l=5,nhad
               CALL boostl(e,pboo2,dps_hadron_pmom(l,1:4))
            ENDDO
         ENDIF
         DO l=5,nhad
            eta=xFeynman(dps_hadron_pmom(l,1:4),e)
            IF(eta.GT.xFcut(l).OR.eta.LT.xFcutlow(l))THEN
               icut=0
               EXIT
            ENDIF
         ENDDO
         IF(.NOT.labeqcoll)THEN
            pboo2(3)=-pboo2(3)
            ! boost back to the lab frame
            DO l=5,nhad
               CALL boostl(e,pboo2,dps_hadron_pmom(l,1:4))
            ENDDO
         ENDIF
      ENDIF
    END SUBROUTINE Cuts_pp_aajj_dps

    SUBROUTINE Cuts_pp_jj_aj_aa(icut)
      IMPLICIT NONE
      INTEGER,INTENT(OUT)::icut
      INTEGER::l,l1,l2,flag,iil,iil1,iil2
      REAL(KIND(1d0))::s,d1,d2,dr,pt,eta,aaa,bbb,ptcut
      REAL(KIND(1d0)),DIMENSION(4)::ponia,pboo2
      REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
      REAL(KIND(1d0))::e,q
      icut=0
      ! invariant mass cuts
      DO l1=3,nhad-1
         DO l2=l1+1,nhad
            s=2*scalar_product(dps_pmom(1,l1,1:4),dps_pmom(1,l2,1:4))&
                 +scalar_product(dps_pmom(1,l1,1:4),dps_pmom(1,l1,1:4))&
                 +scalar_product(dps_pmom(1,l2,1:4),dps_pmom(1,l2,1:4))
            IF((ifldps(1,l1).EQ.22.AND.ifldps(1,l2).NE.22).OR.&
                 (ifldps(1,l1).NE.22.AND.ifldps(1,l2).EQ.22))THEN
               ! a j
               IF(s.LT.gmas(5,7)**2)RETURN
            ELSEIF(ifldps(1,l1).EQ.22.AND.ifldps(1,l2).EQ.22)THEN
               ! a a
               IF(s.LT.gmas(5,6)**2)RETURN
            ELSE
               ! j j
               IF(s.LT.gmas(7,8)**2)RETURN
            ENDIF
         ENDDO
      ENDDO
      flag=0
      DO l=3,nhad
         IF(ifldps(1,l).EQ.22)THEN
            ! a
            iil=5
         ELSE
            ! j
            iil=7
         ENDIF
         pt=SQRT(dps_pmom(1,l,1)**2+dps_pmom(1,l,2)**2)
         IF(pt.LT.ptc(iil))THEN
            flag=1
            EXIT
         ENDIF
         eta=prapidity(dps_pmom(1,l,1:4))
         IF(ABS(eta).GT.etac(iil))THEN
            flag=1
            EXIT
         ENDIF
         eta=rapidity(dps_pmom(1,l,1:4))
         IF(absrap)eta=ABS(eta)
         IF(eta.GT.yycut(iil).OR.eta.LT.yycutlow(iil))THEN
            flag=1
            EXIT
         ENDIF
      ENDDO
      IF(flag.EQ.0)THEN
         DO l1=3,nhad-1
            IF(ifldps(1,l1).EQ.22)THEN
               ! a
               iil1=5
            ELSE
               ! j
               iil1=7
            ENDIF
            DO l2=l1+1,nhad
               IF(ifldps(1,l2).EQ.22)THEN
                  ! a2
                  iil2=5
               ELSE
                  ! j
                  iil2=7
               ENDIF
               d1=prapidity(dps_pmom(1,l1,1:4))-prapidity(dps_pmom(1,l2,1:4))
               d2=ph4(dps_pmom(1,l1,1),dps_pmom(1,l1,2),dps_pmom(1,l1,3))&
                    -ph4(dps_pmom(1,l2,1),dps_pmom(1,l2,2),dps_pmom(1,l2,3))
               d2=MIN(DABS(d2),2*pi-DABS(d2))
               IF(d2/pi.GT.1d0)WRITE(*,*)d2/pi
               dr=SQRT(d1**2+d2**2)
               IF(dr.LT.drc(iil1,iil2))THEN
                  flag=1
                  EXIT
               ENDIF
            ENDDO
            IF(flag.EQ.1)EXIT
         ENDDO
         IF(flag.EQ.0)icut=1
      ENDIF
      ! special cutoffs in the literature  
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
            DO l=3,nhad
               CALL boostl(e,pboo2,dps_pmom(1,l,1:4))
            ENDDO
         ENDIF
         DO l=3,nhad
            IF(ifldps(1,l).EQ.22)THEN
               ! a 
               iil=5
            ELSE
               ! j
               iil=7
            ENDIF
            eta=xFeynman(dps_pmom(1,l,1:4),e)
            IF(eta.GT.xFcut(iil).OR.eta.LT.xFcutlow(iil))THEN
               icut=0
               EXIT
            ENDIF
         ENDDO
         IF(.NOT.labeqcoll)THEN
            pboo2(3)=-pboo2(3)
            ! boost back to the lab frame
            DO l=3,nhad
               CALL boostl(e,pboo2,dps_pmom(1,l,1:4))
            ENDDO
         ENDIF
      ENDIF
    END SUBROUTINE Cuts_pp_jj_aj_aa

    SUBROUTINE readcuts_pp_aajj_dps_cuts
      IMPLICIT NONE
      LOGICAL::lexist
      CHARACTER(len=24)::file
      INTEGER::iounit, flag=0
      ! open pp_aajj_DPS_cuts.inp
      INQUIRE(FILE=TRIM(input_dir)//"pp_aajj_DPS_cuts.inp",EXIST=lexist)
      IF(.NOT.lexist)THEN
         PRINT *,"Warning: the file pp_aajj_DPS_cuts.inp does not exist ! STOP !"
         STOP
      ENDIF
      INQUIRE(FILE=TRIM(input_dir)//"pp_aajj_DPS_cuts.inp",OPENED=lexist)
      IF(lexist)THEN
         INQUIRE(FILE=TRIM(input_dir)//"pp_aajj_DPS_cuts.inp",NUMBER=iounit)
         IF(iounit.NE.udefault)THEN
            PRINT *,"WARNING: the pp_aajj_DPS_cuts.inp has been linked with another unit ! Close and reopen !"
            CLOSE(UNIT=iounit)
            OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"pp_aajj_DPS_cuts.inp")
         ENDIF
      ELSE
         OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"pp_aajj_DPS_cuts.inp")
      ENDIF
      flag=1
      aajj_pta1u=-1d0
      aajj_pta1l=0d0
      aajj_etaa1=30d0
      aajj_pta2u=-1d0
      aajj_pta2l=0d0
      aajj_etaa2=30d0
      aajj_ptj1u=-1d0
      aajj_ptj1l=0d0
      aajj_etaj1=30d0
      aajj_ptj2u=-1d0
      aajj_ptj2l=0d0
      aajj_etaj2=0d0
      aajj_pta1u=readvalue_r("pta1u",flag)
      aajj_pta1l=readvalue_r("pta1l",flag)
      IF(aajj_pta1u.GE.0d0.AND.aajj_pta1u.LE.aajj_pta1l)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the leading photon cut"
         WRITE(*,*)"ERROR:Because pta1u<pta1l when pta1u>=0"
         STOP
      ENDIF
      aajj_etaa1=readvalue_r("etaa1",flag)
      IF(aajj_etaa1.LE.0d0)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the leading photon eta cut"
         WRITE(*,*)"ERROR:Because etaa1 <= 0"
         STOP
      ENDIF
      aajj_pta2u=readvalue_r("pta2u",flag)
      aajj_pta2l=readvalue_r("pta2l",flag)
      IF(aajj_pta2u.GE.0d0.AND.aajj_pta2u.LE.aajj_pta2l)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the subleading photon cut"
         WRITE(*,*)"ERROR:Because pta2u<pta2l when pta2u>=0"
         STOP
      ENDIF
      aajj_etaa2=readvalue_r("etaa2",flag)
      IF(aajj_etaa2.LE.0d0)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the subleading photon eta cut"
         WRITE(*,*)"ERROR:Because etaa2 <= 0"
         STOP
      ENDIF
      aajj_ptj1u=readvalue_r("ptj1u",flag)
      aajj_ptj1l=readvalue_r("ptj1l",flag)
      IF(aajj_ptj1u.GE.0d0.AND.aajj_ptj1u.LE.aajj_ptj1l)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the leading jet cut"
         WRITE(*,*)"ERROR:Because ptj1u<ptj1l when ptj1u>=0"
         STOP
      ENDIF
      aajj_etaj1=readvalue_r("etaj1",flag)
      IF(aajj_etaj1.LE.0d0)THEN
          WRITE(*,*)"ERROR:All of events failed to pass the leading jet eta cut"
          WRITE(*,*)"ERROR:Because etaj1 <= 0"
          STOP
      ENDIF
      aajj_ptj2u=readvalue_r("ptj2u",flag)
      aajj_ptj2l=readvalue_r("ptj2l",flag)
      IF(aajj_ptj2u.GE.0d0.AND.aajj_ptj2u.LE.aajj_ptj2l)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the subleading jet cut"
         WRITE(*,*)"ERROR:Because ptj2u<ptj2l when ptj2u>=0"
         STOP
      ENDIF
      aajj_etaj2=readvalue_r("etaj2",flag)
      IF(aajj_etaj2.LE.0d0)THEN
         WRITE(*,*)"ERROR:All of events failed to pass the subleading jet eta cut"
         WRITE(*,*)"ERROR:Because etaj2 <= 0"
         STOP
      ENDIF
      CLOSE(UNIT=udefault)
      WRITE(*,*)'---------------------------------------------------'
      WRITE(*,*)'    the special cuts for p p > diphoton+dijet + X '
      WRITE(*,*)'    with double parton scattering (DPS) '
      WRITE(*,*)'    the lower pt cut of leading photon    = ',aajj_pta1l
      WRITE(*,*)'    the upper pt cut of leading photon    = ',aajj_pta1u
      WRITE(*,*)'    the eta cut of leading photon         = ',aajj_etaa1
      WRITE(*,*)'    the lower pt cut of subleading photon = ',aajj_pta2l
      WRITE(*,*)'    the upper pt cut of subleading photon = ',aajj_pta2u
      WRITE(*,*)'    the eta cut of subleading photon      = ',aajj_etaa2
      WRITE(*,*)'    the lower pt cut of leading jet       = ',aajj_ptj1l
      WRITE(*,*)'    the upper pt cut of leading jet       = ',aajj_ptj1u
      WRITE(*,*)'    the eta cut of leading jet            = ',aajj_etaj1
      WRITE(*,*)'    the lower pt cut of subleading jet    = ',aajj_ptj2l
      WRITE(*,*)'    the upper pt cut of subleading jet    = ',aajj_ptj2u
      WRITE(*,*)'    the eta cut of subleading jet         = ',aajj_etaj2
      WRITE(*,*)'---------------------------------------------------'
      RETURN
    END SUBROUTINE readcuts_pp_aajj_dps_cuts

    SUBROUTINE Cuts_pp_aajj_dps_cuts(icut)
      IMPLICIT NONE
      INTEGER,INTENT(OUT)::icut
      REAL(KIND(1d0)),DIMENSION(4)::pmom_a1,pmom_a2,pmom_j1,pmom_j2
      REAL(KIND(1d0))::pt1,pt2,eta
      icut=0
      pt1=dps_hadron_pmom(5,1)**2+dps_hadron_pmom(5,2)**2
      pt2=dps_hadron_pmom(6,1)**2+dps_hadron_pmom(6,2)**2
      IF(pt1.GE.pt2)THEN
         pmom_a1(1:4)=dps_hadron_pmom(5,1:4)
         pmom_a2(1:4)=dps_hadron_pmom(6,1:4)
      ELSE
         pmom_a1(1:4)=dps_hadron_pmom(6,1:4)
         pmom_a2(1:4)=dps_hadron_pmom(5,1:4)
      ENDIF
      pt1=dps_hadron_pmom(7,1)**2+dps_hadron_pmom(7,2)**2
      pt2=dps_hadron_pmom(8,1)**2+dps_hadron_pmom(8,2)**2
      IF(pt1.GE.pt2)THEN
         pmom_j1(1:4)=dps_hadron_pmom(7,1:4)
         pmom_j2(1:4)=dps_hadron_pmom(8,1:4)
      ELSE
         pmom_j1(1:4)=dps_hadron_pmom(8,1:4)
         pmom_j2(1:4)=dps_hadron_pmom(7,1:4)
      ENDIF
      pt1=DSQRT(pmom_a1(1)**2+pmom_a1(2)**2)
      IF(pt1.LT.aajj_pta1l)RETURN
      IF(aajj_pta1u.GE.0d0.AND.pt1.GT.aajj_pta1u)RETURN
      eta=prapidity(pmom_a1(1:4))
      IF(ABS(eta).GT.aajj_etaa1)RETURN
      pt1=DSQRT(pmom_a2(1)**2+pmom_a2(2)**2)
      IF(pt1.LT.aajj_pta2l)RETURN
      IF(aajj_pta2u.GE.0d0.AND.pt1.GT.aajj_pta2u)RETURN
      eta=prapidity(pmom_a2(1:4))
      IF(ABS(eta).GT.aajj_etaa2)RETURN
      pt1=DSQRT(pmom_j1(1)**2+pmom_j1(2)**2)
      IF(pt1.LT.aajj_ptj1l)RETURN
      IF(aajj_ptj1u.GE.0d0.AND.pt1.GT.aajj_ptj1u)RETURN
      eta=prapidity(pmom_j1(1:4))
      IF(ABS(eta).GT.aajj_etaj1)RETURN
      pt1=DSQRT(pmom_j2(1)**2+pmom_j2(2)**2)
      IF(pt1.LT.aajj_ptj2l)RETURN
      IF(aajj_ptj2u.GE.0d0.AND.pt1.GT.aajj_ptj2u)RETURN
      eta=prapidity(pmom_j2(1:4))
      IF(ABS(eta).GT.aajj_etaj2)RETURN
      icut=1
      RETURN
    END SUBROUTINE Cuts_pp_aajj_dps_cuts
END MODULE pp_aajj_DPS_cuts
