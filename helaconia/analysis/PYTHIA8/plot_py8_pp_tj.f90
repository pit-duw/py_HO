MODULE plot_py8_user
!============================================================
!
! PLOT FILE WITH PYTHIA8 FOR USER
! This is an example file for single top p p > t j
!
!============================================================
  IMPLICIT NONE
  INTEGER::nwgt_analysis
  INTEGER::num_plots=20
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE plot_py8_begin(nwgt,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt
    INTEGER,PARAMETER::max_weight=1
    CHARACTER(len=15),DIMENSION(max_weight),INTENT(IN)::weights_info
    CHARACTER(len=15)::weights_info2
    INTEGER::i,kk,l
!    character*5 cc(2)
!      data cc/'     ',' Born'/
    !INCLUDE '../hbook/dbookf90.inc'
    INCLUDE 'dbookf90.inc' ! it will be copied to shower directory
    CALL INIHIST
    CALL INIHIST_3D
    nwgt_analysis=nwgt
    IF(nwgt_analysis*16.GT.nplots/4)THEN
       WRITE(*,*) 'error in plot_begin: '//&
            'too many histograms, increase NPLOTS to',&
            nwgt_analysis*num_plots*4 ! 4 is internal in dbook,each plot has 4 histogram to estimate error
       STOP
    ENDIF
    DO kk=1,nwgt_analysis
       weights_info2=weights_info(kk)
       l=(kk-1)*num_plots
       ! bookup is defined in analysis/hbook/dbook.f
       ! NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE
       !
       ! enlarge num_plots at the begining
       !
       ! NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE
       CALL bookup(l+ 1,'total rate    '//TRIM(weights_info2),&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'t pt  '//TRIM(weights_info2),&
            2d0,0d0,200d0)
       CALL bookup(l+ 3,'t y '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+ 4,'t eta '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+ 5,'j1 pt '//TRIM(weights_info2),&
            2d0,0d0,200d0)
       CALL bookup(l+ 6,'j1 y '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+ 7,'j1 eta '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+ 8,'j2 pt '//TRIM(weights_info2),&
            2d0,0d0,200d0)
       CALL bookup(l+ 9,'j2 y '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+10,'j2 eta '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+11,'bj1 pt '//TRIM(weights_info2),&
            2d0,0d0,200d0)
       CALL bookup(l+12,'bj1 y '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+13,'bj1 eta '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+14,'bj2 pt '//TRIM(weights_info2),&
            2d0,0d0,200d0)
       CALL bookup(l+15,'bj2 y '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+16,'bj2 eta '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+17,'syst pt '//TRIM(weights_info2),&
            2d0,0d0,200d0)
       CALL bookup(l+18,'syst y '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+19,'syst eta '//TRIM(weights_info2),&
            0.1d0,-5d0,5d0)
       CALL bookup(l+20,'syst mass '//TRIM(weights_info2),&
            5d0,0d0,500d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_py8_begin

  SUBROUTINE plot_py8_end(xnorm)
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l
    !INCLUDE '../hbook/dbookf90.inc'
    INCLUDE 'dbookf90.inc' ! it will be copied to shower directory
    CHARACTER(len=400)::psfile,rootfile
    CALL mclear
    CALL mclear_3d
    DO i=1,NPLOTS
       CALL mopera(i,'F',i,i,xnorm,0.d0)
       CALL mfinal(i)
       CALL mopera_3d(i,'F',i,i,xnorm,0.d0)
       CALL mfinal_3d(i)
    ENDDO
    ytit='sigma per bin '
    CALL open_topdrawer_file
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL multitop(l+ 1,n_width,n_height,'total rate   ',ytit,'LIN')
       CALL multitop(l+ 2,n_width,n_height,'t pt  ',ytit,'LOG')
       CALL multitop(l+ 3,n_width,n_height,'t y ',ytit,'LOG')
       CALL multitop(l+ 4,n_width,n_height,'t eta ',ytit,'LOG')
       CALL multitop(l+ 5,n_width,n_height,'j1 pt ',ytit,'LOG')
       CALL multitop(l+ 6,n_width,n_height,'j1 y ',ytit,'LOG')
       CALL multitop(l+ 7,n_width,n_height,'j1 eta ',ytit,'LOG')
       CALL multitop(l+ 8,n_width,n_height,'j2 pt ',ytit,'LOG')
       CALL multitop(l+ 9,n_width,n_height,'j2 y ',ytit,'LOG')
       CALL multitop(l+10,n_width,n_height,'j2 eta ',ytit,'LOG')
       CALL multitop(l+11,n_width,n_height,'bj1 pt ',ytit,'LOG')
       CALL multitop(l+12,n_width,n_height,'bj1 y ',ytit,'LOG')
       CALL multitop(l+13,n_width,n_height,'bj1 eta ',ytit,'LOG')
       CALL multitop(l+14,n_width,n_height,'bj2 pt ',ytit,'LOG')
       CALL multitop(l+15,n_width,n_height,'bj2 y ',ytit,'LOG')
       CALL multitop(l+16,n_width,n_height,'bj2 eta ',ytit,'LOG')
       CALL multitop(l+17,n_width,n_height,'syst pt ',ytit,'LOG')
       CALL multitop(l+18,n_width,n_height,'syst y ',ytit,'LOG')
       CALL multitop(l+19,n_width,n_height,'syst eta ',ytit,'LOG')
       CALL multitop(l+20,n_width,n_height,'syst mass ',ytit,'LOG')
    END DO
    CALL close_topdrawer_file
    CALL open_gnuplot_file
    psfile="helaconia_py8.ps"
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
       CALL MGNUPLOT(l+ 2,'t pt  ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 3,'t y ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 4,'t eta ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 5,'j1 pt ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 6,'j1 y ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 7,'j1 eta ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 8,'j2 pt ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 9,'j2 y ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+10,'j2 eta ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+11,'bj1 pt ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+12,'bj1 y ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+13,'bj1 eta ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+14,'bj2 pt ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+15,'bj2 y ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+16,'bj2 eta ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+17,'syst pt ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+18,'syst y ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+19,'syst eta ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+20,'syst mass ',ytit,'LOG',psfile)
    ENDDO
    CALL close_gnuplot_file
    CALL open_root_file
    rootfile="helaconia_py8.root"
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
       CALL MROOTPLOT(l+ 2,'t pt  ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 3,'t y ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 4,'t eta ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 5,'j1 pt ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 6,'j1 y ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 7,'j1 eta ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 8,'j2 pt ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 9,'j2 y ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+10,'j2 eta ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+11,'bj1 pt ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+12,'bj1 y ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+13,'bj1 eta ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+14,'bj2 pt ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+15,'bj2 y ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+16,'bj2 eta ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+17,'syst pt ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+18,'syst y ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+19,'syst eta ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+20,'syst mass ',ytit,'LOG',rootfile)
    ENDDO
    ! wirte the end of the .C file for root
    WRITE (96, 100)
100 FORMAT(/1X,&
         ' hohisto -> cd();',/1X,&
         ' if (histos -> GetEntries() > 0 ) then {',/1X,&
         '  histos->Write();',/1X,&
         '  hohisto -> Close();',/1X,&
         ' }',/1X,'}')
    CALL close_root_file
    RETURN                
  END SUBROUTINE plot_py8_end

  SUBROUTINE plot_py8_fill(wgts)
    IMPLICIT NONE
    !INCLUDE '../include/HEPMC90.INC'
    INCLUDE 'HEPMC90.INC' ! it will be copied to shower directory
    REAL(KIND(1d0)),DIMENSION(4)::psyst
    REAL(KIND(1d0)),DIMENSION(250),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,q,e
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,j,l,kk,IHEP,NB,NT,NBJET,NTRACKS,NJET,IST,ID,ID1
    INTEGER::count_j,count_bj
    REAL(KIND(1d0))::www,JET_KTRADIUS,JET_KTPTMIN,PALG
    INTEGER,PARAMETER::MAXTRACK=2048,MAXJET=2048,MAXNUM=30
    REAL(KIND(1d0)),DIMENSION(4,MAXNUM)::P_TOP,P_B,P_BJET
    REAL(KIND(1d0)),DIMENSION(4,MAXTRACK)::PJET,PTRACK
    INTEGER,DIMENSION(MAXNUM)::IB,BTRACK
    LOGICAL,DIMENSION(MAXNUM)::IS_B_JET
    INTEGER,DIMENSION(MAXTRACK)::JETVEC
    REAL(KIND(1d0))::pttop,ytop,etatop,ptj1,yj1,etaj1
    REAL(KIND(1d0))::ptj2,yj2,etaj2,ptbj1,ybj1,etabj1
    REAL(KIND(1d0))::ptbj2,ybj2,etabj2,ptsyst,ysyst,etasyst,msyst
    ! INCOMING PARTONS MAY TRAVEL IN THE SAME DIRECTION: IT''S 
    ! A POWER-SUPPRESSED EFFECT, SO THROW THE EVENT AWAY
    IF(SIGN(1.D0,PHEP(3,1)).EQ.SIGN(1.D0,PHEP(3,2)))THEN
       WRITE(*,*)"WARNING 1 in plot_py8_fill"
       GOTO 999
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! IMPORT INFORMATION AND DEFINE JETS
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NT=0 ! number of top
    NB=0 ! number of b
    NBJET=0 ! number of b jet
    NTRACKS=0 ! number of tracks
    NJET=0 ! number of non-b jet
    
    DO IHEP=1,NHEP
       IST=ISTHEP(IHEP) ! status
       ID=IDHEP(IHEP) ! pdg
       ID1=IHADR(ID) ! the massive quark's PDG in hadron, defined later

       ! top quark
       IF(ABS(ID).EQ.6)THEN
          DO i=1,4
             P_TOP(i,1)=PHEP(i,IHEP)
          ENDDO
       ENDIF

       ! define particles that to into jet.
       IF(IST.EQ.1.AND.ABS(ID).GE.100)THEN ! final + hadron
          NTRACKS=NTRACKS+1
          IF(ABS(ID1).EQ.5)THEN
             ! FOUND a stable B hadron
             NB=NB+1
             IB(NB)=IHEP ! the index for b hadron
             DO i=1,4
                P_B(i,NB)=PHEP(i,IHEP)
             ENDDO
             BTRACK(NB)=NTRACKS ! the index of tracks (with b hadron)
          ENDIF
          DO i=1,4
             PTRACK(i,NTRACKS)=PHEP(i,IHEP)
          ENDDO
          IF(NTRACKS.EQ.MAXTRACK)THEN
             WRITE(*,*)"ERROR: Too many particles in plot_py8_fill, please increase MATRACK !"
             STOP
          ENDIF
       ENDIF
    ENDDO

    IF(NTRACKS.EQ.0)THEN
       WRITE(*,*)'WARNING:NO TRACKS FOUND, DROP THIS EVENT'
       GOTO 999
    ENDIF

    ! KT Algorithm with fastjet
    NJET=0
    JET_KTRADIUS = 0.7d0 ! R0
    JET_KTPTMIN = 5d0  ! pt cut on the jet
    PALG=1d0 ! 1d0: kt algorithm, 0d0: C/A, -1d0: anti-kt algorithm
    CALL fastjetppgenkt(PTRACK,NTRACKS,JET_KTRADIUS,JET_KTPTMIN,PALG,&
         PJET,NJET,JETVEC) ! fastjet internal function: the returned 
                           ! PJET (momentum),NJET(number of jets)
                           ! JETVEC(itrack) is the jet index for itrack TRACK

    ! check that jets are orderd in pt
    DO i=1,NJET-1
       IF(GETPT(PJET(1:4,i)).LT.GETPT(PJET(1:4,i+1)))THEN
          WRITE(*,*)"ERROR: Jets are not ordered in plot_py8_fill"
          STOP
       ENDIF
    ENDDO

    ! b-jet
    DO i=1,NJET
       is_b_jet(i)=.FALSE.
       DO j=1,NB
          IF(JETVEC(BTRACK(j)).EQ.i)THEN
             is_b_jet(i)=.TRUE.
             EXIT
          ENDIF
       ENDDO
    ENDDO

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    var=1.0d0 ! total cross section
    pttop=GETPT(P_TOP(1:4,1))
    ytop=rapidity(P_TOP(1:4,1))
    etatop=prapidity(P_TOP(1:4,1))

    count_j=0 ! number of non-b jets
    count_bj=0 ! number of b jets
    DO i=1,NJET
       IF(.NOT.IS_B_JET(i))THEN
          count_j=count_j+1
          IF(count_j.EQ.1)THEN
             ! leading non-b jet
             ptj1=GETPT(PJET(1:4,i))
             yj1=rapidity(PJET(1:4,i))
             etaj1=prapidity(PJET(1:4,i))
             DO j=1,4
                psyst(j)=p_top(j,1)+pjet(j,i)
             ENDDO
             ptsyst = GETPT(psyst)
             ysyst = rapidity(psyst)
             etasyst = prapidity(psyst)
             msyst = MINV(psyst)
          ELSEIF(count_j.EQ.2)THEN
             ! subleading non-b jet (can come from shower)
             ptj2 = GETPT(PJET(1:4,i))
             yj2 = rapidity(PJET(1:4,i))
             etaj2 = prapidity(PJET(1:4,i))
          ENDIF
       ELSEIF(IS_B_JET(i))THEN
          count_bj=count_bj+1
          IF(count_bj.EQ.1)THEN
             ! leading b jet
             ptbj1=GETPT(PJET(1:4,i))
             ybj1=rapidity(PJET(1:4,i))
             etabj1=prapidity(PJET(1:4,i))
          ELSEIF(count_bj.EQ.2)THEN
             ! subleading b jet
             ptbj2=GETPT(PJET(1:4,i))
             ybj2=rapidity(PJET(1:4,i))
             etabj2=prapidity(PJET(1:4,i))
          ENDIF
       ENDIF
    ENDDO
    NBJET=count_bj
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,pttop,www)
       CALL mfill(l+3,ytop,www)
       CALL mfill(l+4,etatop,www)
       IF(count_j.GE.1)THEN
          CALL mfill(l+5,ptj1,www)
          CALL mfill(l+6,yj1,www)
          CALL mfill(l+7,etaj1,www)
          CALL mfill(l+17,ptsyst,www)
          CALL mfill(l+18,ysyst,www)
          CALL mfill(l+19,etasyst,www)
          CALL mfill(l+20,msyst,www)
       ENDIF
       IF(count_j.GE.2)THEN
          CALL mfill(l+8,ptj2,www)
          CALL mfill(l+9,yj2,www)
          CALL mfill(l+10,etaj2,www)
       ENDIF
       IF(NBJET.GE.1)THEN
          CALL mfill(l+11,ptbj1,www)
          CALL mfill(l+12,ybj1,www)
          CALL mfill(l+13,etabj1,www)
       ENDIF
       IF(NBJET.GE.2)THEN
          CALL mfill(l+14,ptbj2,www)
          CALL mfill(l+15,ybj2,www)
          CALL mfill(l+16,etabj2,www)
       ENDIF
    ENDDO
999 RETURN      
  END SUBROUTINE plot_py8_fill

  SUBROUTINE open_topdrawer_file
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE='helaconia_py8.top',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_topdrawer_file

  SUBROUTINE close_topdrawer_file
    IMPLICIT NONE
    CLOSE(99)
    RETURN
  END SUBROUTINE close_topdrawer_file

  SUBROUTINE open_gnuplot_file
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=97,FILE='helaconia_py8.gnu',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_gnuplot_file

  SUBROUTINE close_gnuplot_file
    IMPLICIT NONE
    CLOSE(97)
    RETURN
  END SUBROUTINE close_gnuplot_file

  SUBROUTINE open_root_file
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=96,FILE='helaconia_py8.C',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_root_file

  SUBROUTINE close_root_file
    IMPLICIT NONE
    CLOSE(96)
    RETURN
  END SUBROUTINE close_root_file

! Kinematic functions
  FUNCTION rapidity(p)
    ! rapidity
    REAL(KIND(1d0)),DIMENSION(4)::p
    REAL(KIND(1d0))::rapidity,c
    IF(p(4).EQ.0d0)THEN
       rapidity = 0d0
       RETURN
    ENDIF
    c=p(3)/ABS(p(4))
    IF(ABS(c).GE.1d0)THEN
       rapidity =0d0
    ELSE
       rapidity=0.5d0*LOG((1+c)/(1-c))
    ENDIF
  END FUNCTION rapidity

  FUNCTION prapidity(p)
    ! psudo-rapidity
    REAL(KIND(1d0)),DIMENSION(4)::p
    REAL(KIND(1d0))::prapidity,c
    IF(p(1)**2+p(2)**2+p(3)**2.LE.0d0)THEN
       prapidity = 0d0
       RETURN
    ENDIF
    c=p(3)/SQRT(p(1)**2+p(2)**2+p(3)**2)
    IF(ABS(c).EQ.1d0)THEN
       prapidity = 0d0
    ELSE
       prapidity=0.5d0*LOG((1+c)/(1-c))
    ENDIF
  END FUNCTION prapidity

  FUNCTION GETPT(p)
    ! pt
    REAL(KIND(1d0)),DIMENSION(4)::p
    REAL(KIND(1d0))::GETPT
    GETPT=DSQRT(p(1)**2+p(2)**2)
    RETURN
  END FUNCTION GETPT

  FUNCTION MINV(p)
    ! invariant mass (NOT mass squared)
    ! timelike : positive
    ! spacelike : negative
    REAL(KIND(1d0)),DIMENSION(4)::p
    REAL(KIND(1d0))::MINV
    MINV=p(4)**2-p(1)**2-p(2)**2-p(3)**2
    IF(MINV.GE.0d0)THEN
       MINV=DSQRT(MINV)
    ELSE
       MINV=-DSQRT(ABS(MINV))
    ENDIF
    RETURN
  END FUNCTION MINV

  ! returns the PDG code of the heavier quark in the hadron of PDG code ID
  FUNCTION IHADR(ID)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::ID
    INTEGER::IHADR,ID1
    IF(ID.NE.0)THEN
       ID1=ABS(ID)
       IF(ID1.GT.10000)ID1=ID1-1000*INT(ID1/1000)
       IHADR=ID1/(10**INT(LOG10(DFLOAT(ID1))))
    ELSE
       IHADR=0
    ENDIF
    RETURN
  END FUNCTION IHADR

END MODULE plot_py8_user
