MODULE plot_user
!============================================================
!
! PLOT FILE FOR USER
! p p > psi + j j j in LHCb 2.0<y<4.5
!
!============================================================
  USE Helac_Global
  USE Kinetic_Func ! provides kinematic functions
  IMPLICIT NONE
  INTEGER::nwgt_analysis
  INTEGER::num_plots=6
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE plot_begin(nwgt,max_weight,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt,max_weight
    !INTEGER,PARAMETER::num_plots=1
    CHARACTER(len=*),DIMENSION(max_weight),INTENT(IN)::weights_info
    INTEGER::i,kk,l
!    character*5 cc(2)
!      data cc/'     ',' Born'/
    INCLUDE '../hbook/dbookf90.inc'
    CALL INIHIST
    nwgt_analysis=nwgt
    IF(emep_ISR_shower.EQ.1)THEN
       WRITE(*,*)"Please change a plot_user.f90 !!"
       STOP
    ENDIF
    IF(nwgt_analysis*16.GT.nplots/4)THEN
       WRITE(*,*) 'error in plot_begin: '//&
            'too many histograms, increase NPLOTS to',&
            nwgt_analysis*num_plots*4 ! 4 is internal in dbook,each plot has 4 histogram to estimate error
       STOP
    ENDIF
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       ! bookup is defined in analysis/hbook/dbook.f
       ! NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE
       !
       ! enlarge num_plots at the begining
       !
       ! NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE
       CALL bookup(l+ 1,'total rate    '//weights_info(kk),&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'pt of psi in 2.0<y<2.5 '//weights_info(kk),&
            1.0d0,5d0,105d0)
       CALL bookup(l+ 3,'pt of psi in 2.5<y<3.0 '//weights_info(kk),&
            1.0d0,5d0,105d0)
       CALL bookup(l+ 4,'pt of psi in 3.0<y<3.5 '//weights_info(kk),&
            1.0d0,5d0,105d0)
       CALL bookup(l+ 5,'pt of psi in 3.5<y<4.0 '//weights_info(kk),&
            1.0d0,5d0,105d0)
       CALL bookup(l+ 6,'pt of psi in 4.0<y<4.5 '//weights_info(kk),&
            1.0d0,5d0,105d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_begin

  SUBROUTINE plot_end(xnorm)
    USE Helac_Global
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l
    INCLUDE '../hbook/dbookf90.inc'
    CHARACTER(len=400)::psfile,rootfile
    !IF(topdrawer_output)CALL open_topdrawer_file
    !IF(gnuplot_output)CALL open_gnuplot_file
    !IF(root_output)CALL open_root_file
    CALL mclear
    DO i=1,NPLOTS
       CALL mopera(i,'+',i,i,xnorm,0.d0)
       CALL mfinal(i)
    ENDDO
    ytit='sigma per bin '
    IF(topdrawer_output)THEN
       CALL open_topdrawer_file
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL multitop(l+ 1,n_width,n_height,'total rate   ',ytit,'LIN')
          CALL multitop(l+ 2,n_width,n_height,'pt of psi in 2.0<y<2.5  ',ytit,'LOG')
          CALL multitop(l+ 3,n_width,n_height,'pt of psi in 2.5<y<3.0  ',ytit,'LOG')
          CALL multitop(l+ 4,n_width,n_height,'pt of psi in 3.0<y<3.5  ',ytit,'LOG')
          CALL multitop(l+ 5,n_width,n_height,'pt of psi in 3.5<y<4.0  ',ytit,'LOG')
          CALL multitop(l+ 6,n_width,n_height,'pt of psi in 4.0<y<4.5  ',ytit,'LOG')
       END DO
       CALL close_topdrawer_file
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file
       psfile=TRIM(process(2:10*nhad+1))//".ps"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 2,'pt of psi in 2.0<y<2.5 ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 3,'pt of psi in 2.5<y<3.0 ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 4,'pt of psi in 3.0<y<3.5 ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 5,'pt of psi in 3.5<y<4.0 ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 6,'pt of psi in 4.0<y<4.5 ',ytit,'LOG',psfile)
       ENDDO
       CALL close_gnuplot_file
    ENDIF
    IF(root_output)THEN
       CALL open_root_file
       rootfile=TRIM(process(2:10*nhad+1))//".root"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 2,'pt of psi in 2.0<y<2.5  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 3,'pt of psi in 2.5<y<3.0  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 4,'pt of psi in 3.0<y<3.5  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 5,'pt of psi in 3.5<y<4.0  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 6,'pt of psi in 4.0<y<4.5  ',ytit,'LOG',rootfile)
       ENDDO
       ! wirte the end of the .C file for root
       WRITE (96, 100)
100    FORMAT(/1X,&
            ' hohisto -> cd();',/1X,&
            ' if (histos -> GetEntries() > 0 ) then {',/1X,&
            '  histos->Write();',/1X,&
            '  hohisto -> Close();',/1X,&
            ' }',/1X,'}')
       CALL close_root_file
    ENDIF
    RETURN                
  END SUBROUTINE plot_end

  SUBROUTINE plot_fill(wgts)
    USE Helac_Global
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(20,4)::phad
    REAL(KIND(1d0)),DIMENSION(4)::pboo,pini
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,q,e
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,kk,l
    REAL(KIND(1d0))::www,yy,pt
    LOGICAL::y2025,y2530,y3035,y3540,y4045
    ! WARNING*WARNING*WARNING*WARNING*WARNING*WARNING*WARNING
    !
    ! 1)The generated momenta phegas_pmom is in the partonic c.m. frame when imode.EQ.0 without decay
    ! boost it to lab frame first
    ! 2)phegas_pmom is the partonic momenta,i.e. if there is a final QQ~
    ! bound state, the hadronic momenta is pQ+pQ~
    ! 
    ! WARNING*WARNING*WARNING*WARNING*WARNING*WARNING*WARNING
    DO i=1,n
       p(i,1:4)=phegas_pmom(i,1:4)
    ENDDO
    IF((istruc.EQ.1.OR..NOT.labeqcoll).AND.imode.EQ.0.AND.NDecayChains.EQ.0)THEN
       q=ehat
       e=q/SQRT(xp1*xp2)
       exp1=xp1*ebeam(1)
       exp2=xp2*ebeam(2)
       pboo(4)=(exp1+exp2)
       pboo(3)=(exp1-exp2)
       pboo(1:2)=0
       ! boost to the lab frame
       DO i=1,n
          CALL boostl(q,pboo,p(i,1:4))
       ENDDO
    ENDIF
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
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    var=1.0d0 ! total cross section
    pt=DSQRT(phad(3,1)**2+phad(3,2)**2)
    yy=rapidity(phad(3,1:4))
    y2025=yy.GE.2.0d0.AND.yy.LT.2.5d0
    y2530=yy.GE.2.5d0.AND.yy.LT.3.0d0
    y3035=yy.GE.3.0d0.AND.yy.LT.3.5d0
    y3540=yy.GE.3.5d0.AND.yy.LT.4.0d0
    y4045=yy.GE.4.0d0.AND.yy.LE.4.5d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       IF(y2025)CALL mfill(l+2,pt,www)
       IF(y2530)CALL mfill(l+3,pt,www)
       IF(y3035)CALL mfill(l+4,pt,www)
       IF(y3540)CALL mfill(l+5,pt,www)
       IF(y4045)CALL mfill(l+6,pt,www)
    ENDDO
    ! boost back to c.m. frame
    IF((istruc.EQ.1.OR..NOT.labeqcoll).AND.imode.EQ.0.AND.NDecayChains.EQ.0)THEN
       pboo(3)=-(exp1-exp2)
       DO i=1,n
          CALL boostl(q,pboo,p(i,1:4))
       ENDDO
    ENDIF
    RETURN      
  END SUBROUTINE plot_fill

  SUBROUTINE open_topdrawer_file
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    !common/cuseitmax/useitmax                                                  
    OPEN(unit=99,FILE=TRIM(output_dir)//TRIM(process(2:10*nhad+1))//'.top',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_topdrawer_file

  SUBROUTINE close_topdrawer_file
    IMPLICIT NONE
    CLOSE(99)
    RETURN
  END SUBROUTINE close_topdrawer_file

  SUBROUTINE open_gnuplot_file
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=97,FILE=TRIM(output_dir)//TRIM(process(2:10*nhad+1))//'.gnu',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_gnuplot_file

  SUBROUTINE close_gnuplot_file
    IMPLICIT NONE
    CLOSE(97)
    RETURN
  END SUBROUTINE close_gnuplot_file

  SUBROUTINE open_root_file
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=96,FILE=TRIM(output_dir)//TRIM(process(2:10*nhad+1))//'.C',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_root_file

  SUBROUTINE close_root_file
    IMPLICIT NONE
    CLOSE(96)
    RETURN
  END SUBROUTINE close_root_file

  SUBROUTINE plot_fill_Decay(wgts)
    USE Decay_interface
    USE QEDPS_interface
    USE Helac_Func_1
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    INCLUDE "../../decay/stdhep.inc"
    REAL(KIND(1d0))::var,costhl,Q,www
    INTEGER::i,m,l,kk
    REAL(KIND(1d0)),DIMENSION(4)::Pl,PV,PBOO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    var=1.0d0 ! total cross section
    !nphoton=0
    !DO i=1,NHEP
    !   IF(ISHEP(i).EQ.1.AND.ISTHEP(i).EQ.0.AND.JMOHEP(1,i).GT.0.AND.JMOHEP(1,i).LE.NHEP)THEN
    !      IF(IDHEP(i).EQ.22.AND.ISHEP(JMOHEP(1,i)).EQ.-1)THEN
    !         nphoton=nphoton+1
    !      ENDIF
    !   ENDIF
    !ENDDO
    ! select the first decayed stable lepton
    DO i=1,NHEP
       IF(ISHEP(i).EQ.-1)CYCLE
       IF(ISTHEP(i).EQ.1)CYCLE
       IF(JMOHEP(1,i).LE.0.OR.JMOHEP(1,i).GT.NHEP)CYCLE
       IF(ISTHEP(JMOHEP(1,i)).NE.1.OR.ISHEP(JMOHEP(1,i)).EQ.-1)CYCLE
       IF(.NOT.lepton_pdg(IDHEP(i)))CYCLE
       Pl(1:4)=PHEP(1:4,i)
       m=JMOHEP(1,i)
       PV(1:4)=PHEP(1:4,m)
       EXIT   
    ENDDO
    Q=scalar_product(PV,PV)
    IF(Q.LT.0d0)THEN
       PRINT *,"Q<0 in plot_fill_Decay !"
       STOP
    ENDIF
    Q=DSQRT(Q)
    PBOO(1:3)=-PV(1:3)
    PBOO(4)=PV(4)
    CALL BOOSTL(Q,PBOO,Pl) ! boost to the rest frame of V
    costhl=cosij(PV(1:3),Pl(1:3))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,costhl,www)
       ! dummy ones
       ! please replace the correct one if you want to include the ISR from QEDPS
       CALL mfill(l+3,1d0,www)
       CALL mfill(l+4,0.95d0,www)
    ENDDO
    RETURN
  END SUBROUTINE plot_fill_Decay

  SUBROUTINE plot_fill_QEDPS(wgts)
    USE QEDPS_interface
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    ! QEDPS parameters                                                      
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    REAL(KIND(1d0)),DIMENSION(4)::pini
    REAL(KIND(1d0))::wgt,var,frache
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,kk,l,nphoton
    REAL(KIND(1d0))::www,yftot
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! 
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    var=1.0d0 ! total cross section
    pini(1:4)=PLPTN(4:7,1)+PLPTN(4:7,2)
    yftot=rapidity(pini)
    nphoton=0
    DO i=1,NPTCL
       IF(NLPTN(1,i).EQ.22)THEN
          nphoton=nphoton+1
       ENDIF
    ENDDO
    frache=MIN(Q2OUT_QEDPS/Q2MAX_QEDPS,1d0)
    frache=DSQRT(frache)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,yftot,www)
       CALL mfill(l+3,DBLE(nphoton),www)
       CALL mfill(l+4,frache,www)
    ENDDO
    RETURN
  END SUBROUTINE plot_fill_QEDPS
END MODULE plot_user
