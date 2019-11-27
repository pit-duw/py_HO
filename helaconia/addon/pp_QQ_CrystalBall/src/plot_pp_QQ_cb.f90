MODULE plot_pp_QQ_cb
  USE HELAC_Global
  USE pp_QQ_cb_global
  USE Kinetic_Func 
  IMPLICIT NONE
  INTEGER,PARAMETER::max_weight=1
  INTEGER::nwgt_analysis
  INTEGER::num_plots=3
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE initplot_pp_QQ_cb
    IMPLICIT NONE
    INTEGER::nwgt
    CHARACTER(len=15),DIMENSION(max_weight)::weights_info
    INTEGER::i
    nwgt=1
    weights_info(nwgt)=" CB"
    ! output plot files
    CALL plot_begin_pp_QQ_cb(nwgt,weights_info)
    RETURN
  END SUBROUTINE initplot_pp_QQ_cb

  SUBROUTINE plotout_pp_QQ_cb
    USE MC_VEGAS
    IMPLICIT NONE
    LOGICAL::usexinteg=.FALSE.,mint=.FALSE. ! I have not implemented MINT 
    LOGICAL::useitmax=.TRUE.
    REAL(KIND(1d0))::xnorm
    xnorm=1.d0
    IF(useitmax.AND.gener.EQ.3)xnorm=xnorm/float(itmx) ! use vegas with itmax
    ! output plot files  
    CALL plot_end_pp_QQ_cb(xnorm)
    RETURN
  END SUBROUTINE plotout_pp_QQ_cb

  SUBROUTINE outfun_pp_QQ_cb(www)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::www
    INTEGER::i,j
    INTEGER::nwgt
    REAL(KIND(1d0)),DIMENSION(max_weight)::wgts,wgtden
    nwgt=1
    wgts(1)=www
    ! output plot file 
    CALL plot_fill_pp_QQ_cb(wgts)
    RETURN
  END SUBROUTINE outfun_pp_QQ_cb

  SUBROUTINE plot_begin_pp_QQ_cb(nwgt,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt
    CHARACTER(len=*),DIMENSION(max_weight),INTENT(IN)::weights_info
    INTEGER::i,kk,l
    INCLUDE '../../../analysis/hbook/dbookf90.inc'
    CALL INIHIST
    !CALL INIHIST_3D
    nwgt_analysis=nwgt
    IF(nwgt_analysis*16.GT.nplots/4)THEN
       WRITE(*,*) 'error in plot_begin_pp_QQ_cb: ',&
            nwgt_analysis*num_plots*4 ! 4 is internal in dbook
       STOP
    ENDIF
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL bookup(l+ 1,'total rate    '//weights_info(kk),&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'Q Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 3,'Q y '//weights_info(kk),&
            0.1d0,0d0,10d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_begin_pp_QQ_cb

  SUBROUTINE plot_end_pp_QQ_cb(xnorm)
    USE Helac_Global
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l
    INCLUDE '../../../analysis/hbook/dbookf90.inc'
    CHARACTER(len=400)::psfile,rootfile
    CALL mclear
    !CALL mclear_3d
    DO i=1,NPLOTS
       CALL mopera(i,'+',i,i,xnorm,0.d0)
       CALL mfinal(i)
       !CALL mopera_3d(i,'+',i,i,xnorm,0.d0)
       !CALL mfinal_3d(i)
    ENDDO
    ytit='sigma per bin '
    IF(topdrawer_output)THEN
       CALL open_topdrawer_file_pp_QQ_cb
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL multitop(l+ 1,n_width,n_height,'total rate   ',ytit,'LIN')
          CALL multitop(l+ 2,n_width,n_height,'Q Pt  ',ytit,'LOG')
          CALL multitop(l+ 3,n_width,n_height,'Q y ',ytit,'LOG')
       END DO
       CALL close_topdrawer_file_pp_QQ_cb
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file_pp_QQ_cb
       psfile="pp_QQ_crystalball.ps"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 2,'Q Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 3,'Q y ',ytit,'LOG',psfile)
       ENDDO
       CALL close_gnuplot_file_pp_QQ_cb
    ENDIF
    IF(root_output)THEN
       CALL open_root_file_pp_QQ_cb
       rootfile="pp_QQ_crystalball.root"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 2,'Q Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 3,'Q y ',ytit,'LOG',rootfile)
       ENDDO
       ! wirte the end of the .C file for root
       WRITE (96, 100)
100    FORMAT(/1X,&
            ' hohisto -> cd();',/1X,&
            ' if (histos -> GetEntries() > 0 ) then {',/1X,&
            '  histos->Write();',/1X,&
            '  hohisto -> Close();',/1X,&
            ' }',/1X,'}')
       CALL close_root_file_pp_QQ_cb
    ENDIF
    RETURN
  END SUBROUTINE plot_end_pp_QQ_cb

  SUBROUTINE plot_fill_pp_QQ_cb(wgts)
    USE Helac_Global
    IMPLICIT NONE
    !INCLUDE "pp_psipsi_dps.inc"
    REAL(KIND(1d0)),DIMENSION(4)::pboo,pini
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,q,e
    REAL(KIND(1d0))::exp1,exp2,pt
    INTEGER::i,kk,l
    REAL(KIND(1d0))::www,yftot
    !IF(NDecayChains.GT.0)THEN
    !   CALL plot_fill_Decay_pp_psiX_cb(wgts)
    !   RETURN
    !ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                                   
    ! DEFINE OBSERVABLES                                                                                
    !                                                                                                   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    var=1.0d0 ! total cross section
    pini(1:4)=cb_hadron_pmom(1,1:4)+cb_hadron_pmom(2,1:4) ! (px,py,pz,E)
    yftot=rapidity(cb_hadron_pmom(3,1:4))
    pt=DSQRT(cb_hadron_pmom(3,1)**2+cb_hadron_pmom(3,2)**2)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                               
    !                                                                                                   
    ! START TO FILL HISTOGRAM                                                                           
    !                                                                                                   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,pt,www)
       CALL mfill(l+3,yftot,www)
    ENDDO
    ! boost back to c.m. frame
    RETURN
  END SUBROUTINE plot_fill_pp_QQ_cb

  SUBROUTINE open_topdrawer_file_pp_QQ_cb
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE=TRIM(output_dir)//'pp_QQ_crystalball.top',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_topdrawer_file_pp_QQ_cb

  SUBROUTINE close_topdrawer_file_pp_QQ_cb
    IMPLICIT NONE
    CLOSE(99)
    RETURN
  END SUBROUTINE close_topdrawer_file_pp_QQ_cb

  SUBROUTINE open_gnuplot_file_pp_QQ_cb
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=97,FILE=TRIM(output_dir)//'pp_QQ_crystalball.gnu',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_gnuplot_file_pp_QQ_cb

  SUBROUTINE close_gnuplot_file_pp_QQ_cb
    IMPLICIT NONE
    CLOSE(97)
    RETURN
  END SUBROUTINE close_gnuplot_file_pp_QQ_cb

  SUBROUTINE open_root_file_pp_QQ_cb
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=96,FILE=TRIM(output_dir)//'pp_QQ_crystalball.C',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_root_file_pp_QQ_cb

  SUBROUTINE close_root_file_pp_QQ_cb
    IMPLICIT NONE
    CLOSE(96)
    RETURN
  END SUBROUTINE close_root_file_pp_QQ_cb

END MODULE plot_pp_QQ_cb
