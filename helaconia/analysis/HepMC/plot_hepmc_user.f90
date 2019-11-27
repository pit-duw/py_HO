MODULE plot_hepmc_user
!============================================================
!
! PLOT FILE WITH HepMC FOR USER
!
!============================================================
  IMPLICIT NONE
  INTEGER::nwgt_analysis
  INTEGER::num_plots=2
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE plot_hepmc_begin(nwgt,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt
    INTEGER,PARAMETER::max_weight=1
    CHARACTER(len=15),DIMENSION(max_weight),INTENT(IN)::weights_info
    INTEGER::i,kk,l
!    character*5 cc(2)
!      data cc/'     ',' Born'/
    INCLUDE '../hbook/dbookf90.inc'
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
       l=(kk-1)*num_plots
       ! bookup is defined in analysis/hbook/dbook.f
       ! NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE
       !
       ! enlarge num_plots at the begining
       !
       ! NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE**NOTE
       CALL bookup(l+ 1,'total rate    '//weights_info(kk),&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'yini  '//weights_info(kk),&
            0.1d0,-5d0,5d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_hepmc_begin

  SUBROUTINE plot_hepmc_end(xnorm)
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l
    INCLUDE '../hbook/dbookf90.inc'
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
       CALL multitop(l+ 2,n_width,n_height,'yini  ',ytit,'LOG')
    END DO
    CALL close_topdrawer_file
    CALL open_gnuplot_file
    psfile="helaconia_py8.ps"
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
       CALL MGNUPLOT(l+ 2,'yini ',ytit,'LOG',psfile)
    ENDDO
    CALL close_gnuplot_file
    CALL open_root_file
    rootfile="helaconia_py8.root"
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
       CALL MROOTPLOT(l+ 2,'yini ',ytit,'LOG',rootfile)
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
  END SUBROUTINE plot_hepmc_end

  SUBROUTINE plot_hepmc_fill(wgts)
    IMPLICIT NONE
    INCLUDE '../include/HEPMC90.INC'
    REAL(KIND(1d0)),DIMENSION(4)::pboo,pini
    REAL(KIND(1d0)),DIMENSION(250),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,q,e
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,kk,l
    REAL(KIND(1d0))::www,yftot
    ! INCOMING PARTONS MAY TRAVEL IN THE SAME DIRECTION: IT''S 
    ! A POWER-SUPPRESSED EFFECT, SO THROW THE EVENT AWAY
    IF(SIGN(1.D0,PHEP(3,1)).EQ.SIGN(1.D0,PHEP(3,2)))THEN
       WRITE(*,*)"WARNING 1 in plot_hepmc_fill"
       GOTO 999
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    var=1.0d0 ! total cross section
    pini(1:4)=PHEP(1:4,1)+PHEP(1:4,2) ! (px,py,pz,E)
    yftot=rapidity(pini)
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
    ENDDO
999 RETURN      
  END SUBROUTINE plot_hepmc_fill

  SUBROUTINE open_topdrawer_file
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE='helaconia_hepmc.top',&
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
    OPEN(unit=97,FILE='helaconia_hepmc.gnu',&
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
    OPEN(unit=96,FILE='helaconia_hepmc.C',&
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

END MODULE plot_hepmc_user
