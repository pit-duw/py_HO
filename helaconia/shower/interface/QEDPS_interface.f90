MODULE QEDPS_interface
!**************************************************************
!*              HELAC-Onia - QEDPS interface                  *
!*------------------------------------------------------------*
!*  Written by H.-S.Shao, erdissshaw@gmail.com, 05 Jan 2014   *
!*============================================================*
!*                                                            *
!*        NPTCL     :  Number of particles                    *
!*        PLPTN(1, ):  The lightcone fraction                 *
!*        PLPTN(2, ):  virtual mass squared                   *
!*        PLPTN(3, ):  Pt**2                                  *
!*        PLPTN(4, ):  x-component of the four momentum       *
!*        PLPTN(5, ):  y-component of the four momentum       *
!*        PLPTN(6, ):  z-component of the four momentum       *
!*        PLPTN(7, ):  E-component of the four momentum       *
!*        PLPTN(8, ):  + component of the lightcone momentum  *
!*        PLPTN(9, ):  - component of the lightcone momentum  *
!*                                                            *
!*        NLPTN(1, ): particle ID                             *
!*                     electron:11, positron:-11              *
!*                     photon :22                             *
!*                    (PDG standard)                          *
!*        NLPTN(2, ): relative address of the parent          *
!*        NLPTN(3, ): number of the children                  *
!*        NLPTN(4, ): relative address of the first child     *
!*        NLPTN(5, ): status of the process                   *
!*        NLPTN(6, ): spacelike(-1) or timelike(+1)           *
!*                                                            *
!*============================================================*
!*                                                            *
!* WARNING:                                                   *
!*     - The first incident particle is electron,while the    *
!*       second one is positron                               *
!*     - Only photon radiation from ISR                       *
!*     - Momenta (both initial and final) are in electron     *
!*       positron (before shower) center of mass frame        *
!*                                                            *
!**************************************************************
  USE Helac_Global
  IMPLICIT NONE
!  REAL(KIND(1d0)),PUBLIC::Q2MAX_QEDPS
  REAL(KIND(1d0)),DIMENSION(2),PRIVATE::ebeam_old
!  REAL(KIND(1d0)),PUBLIC::Q2OUT_QEDPS
  INTEGER::NVPH_QEDPS
  REAL(KIND(1d0)),PARAMETER::amel=0.511d-3
  LOGICAL,PRIVATE::labeqcoll_QEDPS=.TRUE.
  SAVE ebeam_old,NVPH_QEDPS,labeqcoll_QEDPS
CONTAINS
  SUBROUTINE QEDPSINIT
    IMPLICIT NONE
    ! QEDPS parameters
    !REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    !INTEGER,DIMENSION(10,1000)::NLPTN
    !INTEGER::NPTCL
    !COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    !REAL(KIND(1d0))::q2max
    !REAL(KIND(1d0)),PARAMETER::amel=0.511d-3
    REAL(KIND(1d0))::ainmas
    INTEGER::iseed
    IF(iflh(1).NE.2.OR.iflh(2).NE.-2)THEN
       STOP "ERROR:Wrong beams, please make sure e- e+ (not e+ e-)"
    ENDIF
    ! Mass of incident particles
    ainmas=amel
    ! Note that in HELAC-Onia, the initial particles' masses are always zero
    ! but in QEDPS we should keep it
    ! Q-square maximum
    ! ainmas -> larger,q2max -> smaller
    q2max_QEDPS=2d0*ebeam(1)*ebeam(2)+2d0*DSQRT(ebeam(1)**2-ainmas**2)&
         *DSQRT(ebeam(2)**2-ainmas**2)+ainmas**2*2d0
    ! Random number seed (if negative, use defualt value) for QEDPS
    iseed=-999
    
    CALL QPINIT(q2max_QEDPS,ainmas,iseed)
    ebeam_old(1:2)=ebeam(1:2)
    !EBMUP1=ebeam(1)
    !EBMUP2=ebeam(2)
    IF(ABS(ebeam(1)-ebeam(2))/MAX(ebeam(1)+ebeam(2),1d-6).LT.1d-7)THEN
       labeqcoll_QEDPS=.TRUE.
    ELSE
       labeqcoll_QEDPS=.FALSE.
    ENDIF
    !OPEN(UNIT=20111,FILE=TRIM(output_dir)//"/QEDPS.out",STATUS="unknown")
  END SUBROUTINE QEDPSINIT

  SUBROUTINE QEDPS_ISR_SHOWER
    IMPLICIT NONE
    ! QEDPS parameters
    !REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    !INTEGER,DIMENSION(10,1000)::NLPTN
    !INTEGER::NPTCL
    !COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    REAL(KIND(1d0))::q2in
    ! Q2IN: (CM energy)**2 before radiative photon emission
    ! Q2OUT: (CM energy)**2 after radiative photon emission
    q2in=Q2MAX_QEDPS
    CALL QPGEN(q2in,q2out_QEDPS)
    RETURN
  END SUBROUTINE QEDPS_ISR_SHOWER

  ! after generating the momenta of final states in partonic c.m. frame
  ! NVPH related to the total momenta of electron-positron 
  ! after photon emission in collision frame
  
  SUBROUTINE QEDPSEVNT(icut)
    USE Helac_Func_1
    USE Kinetic_Func
    IMPLICIT NONE
    INTEGER,INTENT(OUT)::icut
    !REAL(KIND(1d0)),INTENT(OUT)::wgt
    ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    INTEGER::PDG_f
    REAL(KIND(1d0))::ef,px_f,py_f,pz_f
    INTEGER::IFRAG
    INTEGER::i,kk
    REAL(KIND(1d0)),DIMENSION(20,4)::p,phad
    INTEGER::nsave=0
    SAVE nsave
    NVPH_QEDPS=NPTCL
    DO i=1,n
       p(i,1:4)=phegas_pmom(i,1:4)
    ENDDO
    ! now everything is generated in partonic c.m frame
    ! when imode=1 or have decay chains, it is generated in lab frame
    ! so, one should make sure ebeam(1)=ebeam(2) and istruct=0
    ! hence,make sure there is no cut in generating phase space
    ! and pTQ=.FALSE. when there is ISR
    kk=1
    DO i=1,nhad
       IF(ABS(iflh(i)).GT.100)THEN
          phad(i,1:4)=p(kk,1:4)+p(kk+1,1:4)
          kk=kk+2
       ELSE
          phad(i,1:4)=p(kk,1:4)
          kk=kk+1
       ENDIF
    ENDDO
    ! make sure nhad > 3
    DO i=3,nhad
       PDG_f=pdgt(iflh(i))
       ef=phad(i,4)
       px_f=phad(i,1)
       py_f=phad(i,2)
       pz_f=phad(i,3)
       IF(i.EQ.3)THEN
          ! first call of QPSET
          ! set necessary information in 'COMMON/QPLIST'
          IFRAG=0
       ELSEIF(i.EQ.nhad)THEN
          ! last call of QPSET
          ! set necessary information in 'COMMON/QPLIST'
          ! and boost them to collision frame
          IFRAG=1
       ELSE
          ! set necessary information in 'COMMON/QPLIST'
          IFRAG=2
       ENDIF
       CALL QPSET(PDG_f,ef,px_f,py_f,pz_f,IFRAG)
    ENDDO
    !ebeam(1:2)=ebeam_old(1:2) ! it is dangerous
    ! cuts in lab
    CALL cuts_emep_QEDPS(icut)
    IF(icut.EQ.1)THEN
       nsave=nsave+1
       !IF(nsave.LE.10)THEN
       !   CALL QPDUMP
       !ENDIF
       IF(imode.EQ.1.OR.NDecayChains.NE.0)THEN
          ! boost from cms frame to lab frame
          CALL cmstolab_QEDPS()
          CALL trans_p_to_zq(n,Phegas_pmom(1:n,1:5),zq(1:n,1:5))
       ENDIF
    ENDIF
  END SUBROUTINE QEDPSEVNT

  SUBROUTINE QEDPSEND
    IMPLICIT NONE
    !CLOSE(20111)
    RETURN
  END SUBROUTINE QEDPSEND

  SUBROUTINE cmstolab_QEDPS()
    USE Kinetic_Func
    ! keep the mass of initial electron-positron
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(4)::pboo
    INTEGER::i
    REAL(KIND(1d0))::Q
    ! first cms to collision frame
    CALL cmstocollision_QEDPS()
    IF(labeqcoll_QEDPS)RETURN
    ! second collision frame to lab frame
    DO i=1,n
       p(i,1:4)=phegas_pmom(i,1:4)
    ENDDO
    pboo(4)=ebeam_old(1)+ebeam_old(2)
    pboo(3)=DSQRT(ebeam_old(1)**2-amel**2)&
         -DSQRT(ebeam_old(2)**2-amel**2)
    pboo(1:2)=0d0
    Q=DSQRT(Q2MAX_QEDPS)
    DO i=1,n
       CALL BOOSTL(Q,pboo,p(i,1:4))
    ENDDO
    phegas_pmom(1:n,1:4)=p(1:n,1:4)
  END SUBROUTINE cmstolab_QEDPS

  SUBROUTINE labtocms_QEDPS()
    USE Kinetic_Func
    ! keep the mass of initial electron-positron
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(4)::pboo
    INTEGER::i
    REAL(KIND(1d0))::Q
    IF(.NOT.labeqcoll_QEDPS)THEN
       ! first lab to collision
       DO i=1,n
          p(i,1:4)=phegas_pmom(i,1:4)
       ENDDO
       pboo(4)=ebeam_old(1)+ebeam_old(2)
       pboo(3)=-DSQRT(ebeam_old(1)**2-amel**2)+&
            DSQRT(ebeam_old(2)**2-amel**2)
       pboo(1:2)=0d0
       Q=DSQRT(Q2MAX_QEDPS)
       DO i=1,n
          CALL BOOSTL(Q,pboo,p(i,1:4))
       ENDDO
       phegas_pmom(1:n,1:4)=p(1:n,1:4)
    ENDIF
    ! second collision to cms
    CALL collisiontocms_QEDPS()
    RETURN
  END SUBROUTINE labtocms_QEDPS

  SUBROUTINE cmstocollision_QEDPS()
    USE Kinetic_Func
    ! keep the mass of initial electron-positron
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(4)::pboo
    INTEGER::i
    REAL(KIND(1d0))::Q
    ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    DO i=1,n
       p(i,1:4)=phegas_pmom(i,1:4)
    ENDDO
    pboo(1)=PLPTN(4,NVPH_QEDPS)
    pboo(2)=PLPTN(5,NVPH_QEDPS)
    pboo(3)=PLPTN(6,NVPH_QEDPS)
    pboo(4)=PLPTN(7,NVPH_QEDPS)
    Q=DSQRT(Q2OUT_QEDPS)
    DO i=1,n
       CALL BOOSTL(Q,pboo,p(i,1:4))
    ENDDO
    phegas_pmom(1:n,1:4)=p(1:n,1:4)
    RETURN
  END SUBROUTINE cmstocollision_QEDPS

  SUBROUTINE collisiontocms_QEDPS()
    USE Kinetic_Func
    ! keep the mass of initial electron-positron
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(4)::pboo
    INTEGER::i
    REAL(KIND(1d0))::Q
    ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    DO i=1,n
       p(i,1:4)=phegas_pmom(i,1:4)
    ENDDO
    pboo(1)=PLPTN(4,NVPH_QEDPS)
    pboo(2)=PLPTN(5,NVPH_QEDPS)
    pboo(3)=-PLPTN(6,NVPH_QEDPS)
    pboo(4)=PLPTN(7,NVPH_QEDPS)
    Q=DSQRT(Q2OUT_QEDPS)
    DO i=1,n
       CALL BOOSTL(Q,pboo,p(i,1:4))
    ENDDO
    phegas_pmom(1:n,1:4)=p(1:n,1:4)
    RETURN
  END SUBROUTINE collisiontocms_QEDPS

  SUBROUTINE cuts_emep_QEDPS(icut)
    USE Kinetic_Func
    IMPLICIT NONE
    INTEGER,INTENT(OUT)::icut
    INTEGER::l,l1,l2,i_nu,i1
    REAL(KIND=DBL)::c,s,pmiss_E,pmiss_x,pmiss_y,pmiss_z,pt_miss
    REAL(KIND=DBL),DIMENSION(4)::pboo
    REAL(KIND=DBL)::exp1,exp2,q
    INTEGER::i,kk,kk1,kk2,addkk
        ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    icut=0
    ! boost from collision frame to lab frame
    CALL collision2lab_QEDPS
    kk=3
    DO i=1,NPTCL
       IF(NLPTN(6,i).EQ.-1)CYCLE
       ! don't cut on the PS photon
       IF(NLPTN(1,i).EQ.22.AND.i.LE.NVPH_QEDPS)CYCLE
       IF(i-NVPH_QEDPS+2.LT.1.OR.i-NVPH_QEDPS+2.GT.nhad)CYCLE
       IF(ABS(iflh(i-NVPH_QEDPS+2)).GT.100)THEN
          IF(PLPTN(7,i).LT.ec(kk)+ec(kk+1))RETURN
          addkk=2
       ELSE
          IF(PLPTN(7,i).LT.ec(kk))RETURN
          addkk=1
       ENDIF
       c=PLPTN(6,i)/PLPTN(7,i)
       IF(c.GT.c1(kk))RETURN
       IF(-c.GT.c2(kk))RETURN
       kk=addkk+kk
    ENDDO
    kk1=3
    DO l1=1,NPTCL
       IF(NLPTN(6,l1).EQ.-1)CYCLE
       IF(NLPTN(1,l1).EQ.22.AND.l1.LE.NVPH_QEDPS)CYCLE
       kk2=3
       DO l2=1,NPTCL
          IF(NLPTN(6,l1).EQ.-1)CYCLE
          IF(NLPTN(1,l1).EQ.22.AND.l1.LE.NVPH_QEDPS)CYCLE
          IF(l1.EQ.l2)CYCLE
          IF(cosij(PLPTN(4:6,l1),PLPTN(4:6,l2)).GT.cc(kk1,kk2))RETURN
          IF(ABS(iflh(l2)).GT.100)THEN
             kk2=kk2+2
          ELSE
             kk2=kk2+1
          ENDIF
       ENDDO
       IF(ABS(iflh(l1)).GT.100)THEN
          kk1=kk1+2
       ELSE
          kk1=kk1+1
       ENDIF
    ENDDO
    DO l1=3,n
       DO l2=3,n
          IF(l1.EQ.l2.OR.Quarkonium3(l1).EQ.Quarkonium3(l2).AND.Quarkonium3(l1).NE.0)CYCLE
          s=2*scalar_product(phegas_pmom(l1,1:4),phegas_pmom(l2,1:4))&
               +scalar_product(phegas_pmom(l1,1:4),phegas_pmom(l1,1:4))&
               +scalar_product(phegas_pmom(l2,1:4),phegas_pmom(l2,1:4))
          IF(s.LT.gmas(l1,l2)**2)RETURN
       ENDDO
    ENDDO
    i_nu=0
    pmiss_E=0.d0
    pmiss_x=0.d0
    pmiss_y=0.d0
    pmiss_z=0.d0
    kk=3
    DO l=1,NPTCL
       IF(NLPTN(6,l).EQ.-1)CYCLE
       IF(NLPTN(1,l).EQ.22.AND.l.LE.NVPH_QEDPS)CYCLE
       kk=kk+1
       IF(iflh(kk).LT.13)THEN
          i1=IABS(MOD(iflh(kk),4))
       ENDIF
       IF(i1.EQ.1)THEN ! neutrino found
          i_nu=1
          pmiss_E = pmiss_E + PLPTN(7,l)
          pmiss_x = pmiss_x + PLPTN(4,l)
          pmiss_y = pmiss_y + PLPTN(5,l)
          pmiss_z = pmiss_z + PLPTN(6,l)
       ENDIF
    ENDDO
    pt_miss=DSQRT(pmiss_x**2+pmiss_y**2)
    IF(i_nu.EQ.1.AND.pt_miss.LT.ptmissc.OR.nhad.EQ.3)RETURN
    icut=1
    RETURN
  END SUBROUTINE cuts_emep_QEDPS

  SUBROUTINE collision2lab_QEDPS
    IMPLICIT NONE
    ! QEDPS parameters                                                      
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    REAL(KIND(1d0)),DIMENSION(4)::PREF,pp
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    INTEGER::i,j
    !REAL(KIND(1d0)),PARAMETER::amel=0.511d-3 ! electron's mass in GeV
    LOGICAL::init=.TRUE.
    REAL(KIND(1d0)),DIMENSION(4,100)::PPPP
    REAL(KIND(1d0))::QQ
    SAVE init,PREF,QQ
    IF(labeqcoll_QEDPS)RETURN
    IF(init)THEN
       PREF(1:2)=0d0
       PREF(4)=ebeam_old(1)+ebeam_old(2)
       PREF(3)=DSQRT(ebeam_old(1)**2-amel**2)&
            -DSQRT(ebeam_old(2)**2-amel**2)
       init=.FALSE.
       QQ=DSQRT(Q2MAX_QEDPS)
    ENDIF
    DO i=1,NPTCL
       DO j=1,4
          PPPP(J,I)=PLPTN(J+3,I)
       ENDDO
    ENDDO
    CALL QPBOST(PPPP,PREF,1,NPTCL)
    DO I=1,NPTCL
       DO J=1,4
          PLPTN(J+3,I)=PPPP(J,I)
       ENDDO
    ENDDO
  END SUBROUTINE collision2lab_QEDPS

  SUBROUTINE unwei_writer_QEDPS
    USE setscale
    IMPLICIT NONE
    ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    REAL(KIND(1d0))::scale1
    INTEGER::i,kk,idup,idprup,istup,imothup1,imothup2
    INTEGER,DIMENSION(2,1000)::children
    INTEGER,DIMENSION(1000)::imoth,mother
    INTEGER::i1max,i2max,icol1,icol2
    REAL(KIND(1d0))::px,py,pz,p0,pmass,scalup,vtime,vspin,xwgtup
    idprup = 81
    xwgtup=1
    CALL qcdscale(scale1)
    scalup=scale1
    WRITE(nunit3)NPTCL-1,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    ! the cuts already boost PLPTN into lab frame
    !IF(imode.EQ.0)CALL cmstolab_QEDPS()
    imoth(1:1000)=0
    mother(1:1000)=0
    children(1:2,1:1000)=0
    imoth(1)=1
    imoth(2)=2
    i1max=1
    i2max=2
    DO i=1,NVPH_QEDPS-1 ! only from parton shower
       IF(NLPTN(6,i).EQ.-1)THEN
          IF(NLPTN(3,i).EQ.2.AND.NLPTN(4,i).NE.0)THEN
             children(1,i)=NLPTN(4,i)+i
             children(2,i)=children(1,i)+1
          ENDIF
          IF(NLPTN(2,i).NE.0)THEN
             mother(i)=i-NLPTN(2,i)
             mother(i+1)=mother(i)
             imoth(i)=imoth(mother(i))
             imoth(i+1)=imoth(mother(i))
             IF(imoth(i).EQ.1.AND.i.GT.i1max)i1max=i
             IF(imoth(i).EQ.2.AND.i.GT.i2max)i2max=i
          ENDIF
       ENDIF
    ENDDO
    kk=3
    DO i=1,NPTCL
       IF(i.EQ.NVPH_QEDPS)CYCLE
       idup=NLPTN(1,i)
       istup=NLPTN(6,i)
       IF(i.LT.NVPH_QEDPS)THEN
          IF(imoth(i).EQ.1)THEN
             imothup1=mother(i)
             imothup2=0
          ELSEIF(imoth(i).EQ.2)THEN
             imothup1=0
             imothup2=mother(i)
          ELSE
             imothup1=0
             imothup2=0
          ENDIF
       ELSE
          imothup1=i1max
          imothup2=i2max
       ENDIF
       IF(NLPTN(6,i).EQ.-1)THEN
          icol1=0
          icol2=0
       ELSE
          IF(i.LT.NVPH_QEDPS)THEN
             icol1=0
             icol2=0
          ELSE
             icol1=icol_un(i-NVPH_QEDPS+2,2)+100
             IF(icol1.EQ.100)icol1=0
             icol2=icol_un(i-NVPH_QEDPS+2,1)+100
             IF(icol2.EQ.100)icol2=0
          ENDIF
       ENDIF
       px=PLPTN(4,i)
       py=PLPTN(5,i)
       pz=PLPTN(6,i)
       p0=PLPTN(7,i)
       IF(ABS(NLPTN(1,i)).EQ.11.AND.i.LT.NVPH_QEDPS)THEN
          pmass=amel
       ELSEIF(NLPTN(1,i).EQ.22.OR.i.LT.NVPH_QEDPS)THEN
          pmass=0d0
       ELSE
          pmass=parmas(iflh(i-NVPH_QEDPS+2))
       ENDIF
       vtime=0
       vspin=9
       WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2,&
            px,py,pz,p0,pmass,vtime,vspin
    ENDDO
  END SUBROUTINE unwei_writer_QEDPS
  
END MODULE QEDPS_INTERFACE
