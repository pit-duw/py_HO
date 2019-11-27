MODULE plot_pp_psiX_cb
  USE HELAC_Global
  USE pp_psiX_cb_global
  USE Kinetic_Func 
  IMPLICIT NONE
  INTEGER,PARAMETER::max_weight=1
  INTEGER::nwgt_analysis
  INTEGER::num_plots=3
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE initplot_pp_psiX_cb
    IMPLICIT NONE
    INTEGER::nwgt
    CHARACTER(len=15),DIMENSION(max_weight)::weights_info
    INTEGER::i
    nwgt=1
    weights_info(nwgt)=" CB"
    ! output plot files
    CALL plot_begin_pp_psiX_cb(nwgt,weights_info)
    RETURN
  END SUBROUTINE initplot_pp_psiX_cb

  SUBROUTINE plotout_pp_psiX_cb
    USE MC_VEGAS
    IMPLICIT NONE
    LOGICAL::usexinteg=.FALSE.,mint=.FALSE. ! I have not implemented MINT 
    LOGICAL::useitmax=.TRUE.
    REAL(KIND(1d0))::xnorm
    xnorm=1.d0
    IF(useitmax.AND.gener.EQ.3)xnorm=xnorm/float(itmx) ! use vegas with itmax
    ! output plot files  
    CALL plot_end_pp_psiX_cb(xnorm)
    RETURN
  END SUBROUTINE plotout_pp_psiX_cb

  SUBROUTINE outfun_pp_psiX_cb(www)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::www
    INTEGER::i,j
    INTEGER::nwgt
    REAL(KIND(1d0)),DIMENSION(max_weight)::wgts,wgtden
    nwgt=1
    wgts(1)=www
    ! output plot file 
    CALL plot_fill_pp_psiX_cb(wgts)
    RETURN
  END SUBROUTINE outfun_pp_psiX_cb

  SUBROUTINE plot_begin_pp_psiX_cb(nwgt,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt
    CHARACTER(len=*),DIMENSION(max_weight),INTENT(IN)::weights_info
    INTEGER::i,kk,l
    INCLUDE '../../../analysis/hbook/dbookf90.inc'
    CALL INIHIST
    !CALL INIHIST_3D
    nwgt_analysis=nwgt
    IF(nwgt_analysis*16.GT.nplots/4)THEN
       WRITE(*,*) 'error in plot_begin_pp_psiX_cb: ',&
            nwgt_analysis*num_plots*4 ! 4 is internal in dbook
       STOP
    ENDIF
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL bookup(l+ 1,'total rate    '//weights_info(kk),&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'psi Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 3,'psi y '//weights_info(kk),&
            0.1d0,0d0,10d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_begin_pp_psiX_cb

  SUBROUTINE plot_end_pp_psiX_cb(xnorm)
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
       CALL open_topdrawer_file_pp_psiX_cb
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL multitop(l+ 1,n_width,n_height,'total rate   ',ytit,'LIN')
          CALL multitop(l+ 2,n_width,n_height,'psi Pt  ',ytit,'LOG')
          CALL multitop(l+ 3,n_width,n_height,'psi y ',ytit,'LOG')
       END DO
       CALL close_topdrawer_file_pp_psiX_cb
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file_pp_psiX_cb
       psfile="pp_psiX_crystalball.ps"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 2,'psi Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 3,'psi y ',ytit,'LOG',psfile)
       ENDDO
       CALL close_gnuplot_file_pp_psiX_cb
    ENDIF
    IF(root_output)THEN
       CALL open_root_file_pp_psiX_cb
       rootfile="pp_psiX_crystalball.root"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 2,'psi Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 3,'psi y ',ytit,'LOG',rootfile)
       ENDDO
       ! wirte the end of the .C file for root
       WRITE (96, 100)
100    FORMAT(/1X,&
            ' hohisto -> cd();',/1X,&
            ' if (histos -> GetEntries() > 0 ) then {',/1X,&
            '  histos->Write();',/1X,&
            '  hohisto -> Close();',/1X,&
            ' }',/1X,'}')
       CALL close_root_file_pp_psiX_cb
    ENDIF
    RETURN
  END SUBROUTINE plot_end_pp_psiX_cb

  SUBROUTINE plot_fill_pp_psiX_cb(wgts)
    USE Helac_Global
    IMPLICIT NONE
    !INCLUDE "pp_psipsi_dps.inc"
    REAL(KIND(1d0)),DIMENSION(4)::pboo,pini
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,q,e
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,kk,l
    REAL(KIND(1d0))::www,yftot
    IF(NDecayChains.GT.0)THEN
       CALL plot_fill_Decay_pp_psiX_cb(wgts)
       RETURN
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !                                                                                                   
    ! DEFINE OBSERVABLES                                                                                
    !                                                                                                   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    var=1.0d0 ! total cross section
    pini(1:4)=cb_hadron_pmom(1,1:4)+cb_hadron_pmom(2,1:4) ! (px,py,pz,E)
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
    ! boost back to c.m. frame
    RETURN
  END SUBROUTINE plot_fill_pp_psiX_cb

  SUBROUTINE open_topdrawer_file_pp_psiX_cb
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE=TRIM(output_dir)//'pp_psiX_crystalball.top',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_topdrawer_file_pp_psiX_cb

  SUBROUTINE close_topdrawer_file_pp_psiX_cb
    IMPLICIT NONE
    CLOSE(99)
    RETURN
  END SUBROUTINE close_topdrawer_file_pp_psiX_cb

  SUBROUTINE open_gnuplot_file_pp_psiX_cb
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=97,FILE=TRIM(output_dir)//'pp_psiX_crystalball.gnu',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_gnuplot_file_pp_psiX_cb

  SUBROUTINE close_gnuplot_file_pp_psiX_cb
    IMPLICIT NONE
    CLOSE(97)
    RETURN
  END SUBROUTINE close_gnuplot_file_pp_psiX_cb

  SUBROUTINE open_root_file_pp_psiX_cb
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=96,FILE=TRIM(output_dir)//'pp_psiX_crystalball.C',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_root_file_pp_psiX_cb

  SUBROUTINE close_root_file_pp_psiX_cb
    IMPLICIT NONE
    CLOSE(96)
    RETURN
  END SUBROUTINE close_root_file_pp_psiX_cb

  SUBROUTINE plot_fill_Decay_pp_psiX_cb(wgts)
    USE Helac_Func_1
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    INCLUDE "stdhep_pp_psiX_cb.inc"
    REAL(KIND(1d0))::var,Q,www,y_psi,pt_psi
    INTEGER::i,m,l,kk
    REAL(KIND(1d0)),DIMENSION(1,4)::PJpsi
    REAL(KIND(1d0)),DIMENSION(2,4)::Plep
    REAL(KIND(1d0)),DIMENSION(4)::plep1,Plep2
    INTEGER,DIMENSION(2)::lep2jpsi
    INTEGER,DIMENSION(1)::newmoth
    INTEGER::i_jpsi,i_lep,kk1,kk2
    INTEGER::definemoth
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                               
    !                                                                                                   
    ! DEFINE OBSERVABLES                                                                                
    !                                                                                                   
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    var=1.0d0 ! total cross section
    ! select the first four decayed stable lepton 
    i_jpsi=0
    i_lep=0
    lep2jpsi(1:2)=0
    newmoth(1:1)=0
    DO i=1,NHEP
       IF(ISHEP(i).EQ.-1)CYCLE
       IF(ISTHEP(i).EQ.1)CYCLE
       IF(JMOHEP(1,i).LE.0.OR.JMOHEP(1,i).GT.NHEP)CYCLE
       IF(ISTHEP(JMOHEP(1,i)).NE.1.OR.ISHEP(JMOHEP(1,i)).EQ.-1)CYCLE
       IF(.NOT.lepton_pdg(IDHEP(i)))CYCLE
       i_lep=i_lep+1
       Plep(i_lep,1:4)=PHEP(1:4,i)
       m=JMOHEP(1,i)
       definemoth=0
       DO kk=1,i_jpsi
          IF(newmoth(kk).EQ.m)THEN
             definemoth=kk
             EXIT
          END IF
       ENDDO
       IF(definemoth.NE.0)THEN
          ! the mother has been defined                                                                 
          lep2jpsi(i_lep)=definemoth
       ELSE
          ! create the new mother
          i_jpsi=i_jpsi+1
          PJpsi(i_jpsi,1:4)=PHEP(1:4,m)
          newmoth(i_jpsi)=m
          lep2jpsi(i_lep)=i_jpsi
          IF(i_lep.EQ.2)EXIT
       ENDIF
    ENDDO
    y_psi=DABS(rapidity(PJpsi(1,1:4)))
    pt_psi=PJpsi(1,1)**2+PJpsi(1,2)**2
    pt_psi=DSQRT(DABS(pt_psi))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                                          
    !                                                                                              
    ! START TO FILL HISTOGRAM                                                                      
    !                                                                                              
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,pt_psi,www)
       CALL mfill(l+3,y_psi,www)
    ENDDO
    RETURN
  END SUBROUTINE plot_fill_Decay_pp_psiX_cb
END MODULE plot_pp_psiX_cb
