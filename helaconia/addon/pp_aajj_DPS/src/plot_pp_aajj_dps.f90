MODULE plot_pp_aajj_dps
  USE HELAC_Global
  USE pp_aajj_dps_global
  USE Kinetic_Func 
  IMPLICIT NONE
  INTEGER,PARAMETER::max_weight=1
  INTEGER::nwgt_analysis
  INTEGER::num_plots=21
  INTEGER,PARAMETER::n_width=1,n_height=1
CONTAINS
  SUBROUTINE initplot_pp_aajj_dps
    IMPLICIT NONE
    INTEGER::nwgt
    CHARACTER(len=15),DIMENSION(max_weight)::weights_info
    INTEGER::i
    nwgt=1
    weights_info(nwgt)=" DPS"
    ! output plot files
    CALL plot_begin_pp_aajj_dps(nwgt,weights_info)
    RETURN
  END SUBROUTINE initplot_pp_aajj_dps

  SUBROUTINE plotout_pp_aajj_dps
    USE MC_VEGAS
    IMPLICIT NONE
    LOGICAL::usexinteg=.FALSE.,mint=.FALSE. ! I have not implemented MINT 
    LOGICAL::useitmax=.TRUE.
    REAL(KIND(1d0))::xnorm
    xnorm=1.d0
    IF(useitmax.AND.gener.EQ.3)xnorm=xnorm/float(itmx) ! use vegas with itmax
    ! output plot files  
    CALL plot_end_pp_aajj_dps(xnorm)
    RETURN
  END SUBROUTINE plotout_pp_aajj_dps

  SUBROUTINE outfun_pp_aajj_dps(www)
    IMPLICIT NONE
    REAL(KIND(1d0)),INTENT(IN)::www
    INTEGER::i,j
    INTEGER::nwgt
    REAL(KIND(1d0)),DIMENSION(max_weight)::wgts,wgtden
    nwgt=1
    wgts(1)=www
    ! output plot file 
    CALL plot_fill_pp_aajj_dps(wgts)
    RETURN
  END SUBROUTINE outfun_pp_aajj_dps

  SUBROUTINE plot_begin_pp_aajj_dps(nwgt,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt
    CHARACTER(len=*),DIMENSION(max_weight),INTENT(IN)::weights_info
    INTEGER::i,kk,l
    INCLUDE '../../../analysis/hbook/dbookf90.inc'
    CALL INIHIST
    CALL INIHIST_3D
    nwgt_analysis=nwgt
    IF(nwgt_analysis*16.GT.nplots/4)THEN
       WRITE(*,*) 'error in plot_begin_pp_aajj_dps: ',&
            nwgt_analysis*num_plots*4 ! 4 is internal in dbook
       STOP
    ENDIF
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL bookup(l+ 1,'total rate    '//weights_info(kk),&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'|Delta eta(a,a)|  '//weights_info(kk),&
            0.05d0,0d0,5d0)
       CALL bookup(l+ 3,'|Delta phi(a,a)| '//weights_info(kk),&
            3.1417d0/100d0,-0.0001d0,3.1416d0)
       CALL bookup(l+ 4,'a a inv. mass '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 5,'a a trans mom '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 6,'max |Delta eta(a,j)|'//weights_info(kk),&
            0.05d0,0d0,5d0)
       CALL bookup(l+ 7,'max |Delta phi(a,j)| '//weights_info(kk),&
            3.1417d0/100d0,-0.0001d0,3.1416d0)
       CALL bookup(l+ 8,'max a j inv. mass '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 9,'max a j trans mom '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 10,'|Delta eta(j,j)|'//weights_info(kk),&
            0.05d0,0d0,5d0)
       CALL bookup(l+ 11,'|Delta phi(j,j)| '//weights_info(kk),&
            3.1417d0/100d0,-0.0001d0,3.1416d0)
       CALL bookup(l+ 12,'j j inv. mass '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 13,'j j trans mom '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 14,'a1 Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 15,'a1 Abs(eta) '//weights_info(kk),&
            1d0,5d0,0.1d0)
       CALL bookup(l+ 16,'a2 Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 17,'a2 Abs(eta) '//weights_info(kk),&
            1d0,5d0,0.1d0)
       CALL bookup(l+ 18,'j1 Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 19,'j1 Abs(eta) '//weights_info(kk),&
            1d0,5d0,0.1d0)
       CALL bookup(l+ 20,'j2 Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 21,'j2 Abs(eta) '//weights_info(kk),&
            1d0,5d0,0.1d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_begin_pp_aajj_dps

  SUBROUTINE plot_end_pp_aajj_dps(xnorm)
    USE Helac_Global
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l
    INCLUDE '../../../analysis/hbook/dbookf90.inc'
    CHARACTER(len=400)::psfile,rootfile
    CALL mclear
    CALL mclear_3d
    DO i=1,NPLOTS
       CALL mopera(i,'+',i,i,xnorm,0.d0)
       CALL mfinal(i)
       CALL mopera_3d(i,'+',i,i,xnorm,0.d0)
       CALL mfinal_3d(i)
    ENDDO
    ytit='sigma per bin '
    IF(topdrawer_output)THEN
       CALL open_topdrawer_file_pp_aajj_dps
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL multitop(l+ 1,n_width,n_height,'total rate   ',ytit,'LIN')
          CALL multitop(l+ 2,n_width,n_height,'|Delta eta(a,a)|  ',ytit,'LOG')
          CALL multitop(l+ 3,n_width,n_height,'|Delta phi(a,a)| ',ytit,'LOG')
          CALL multitop(l+ 4,n_width,n_height,'a a inv. mass ',ytit,'LOG')
          CALL multitop(l+ 5,n_width,n_height,'a a trans mom ',ytit,'LOG')
          CALL multitop(l+ 6,n_width,n_height,'max |Delta eta(a,j)| ',ytit,'LOG')
          CALL multitop(l+ 7,n_width,n_height,'max |Delta phi(a,j)| ',ytit,'LOG')
          CALL multitop(l+ 8,n_width,n_height,'max a j inv. mass ',ytit,'LOG')
          CALL multitop(l+ 9,n_width,n_height,'max a j trans mom ',ytit,'LOG')
          CALL multitop(l+ 10,n_width,n_height,'|Delta eta(j,j)| ',ytit,'LOG')
          CALL multitop(l+ 11,n_width,n_height,'|Delta phi(j,j)| ',ytit,'LOG')
          CALL multitop(l+ 12,n_width,n_height,'j j inv. mass ',ytit,'LOG')
          CALL multitop(l+ 13,n_width,n_height,'j j trans mom ',ytit,'LOG')
          CALL multitop(l+ 14,n_width,n_height,'a1 Pt ',ytit,'LOG')
          CALL multitop(l+ 15,n_width,n_height,'a1 Abs(eta) ',ytit,'LOG')
          CALL multitop(l+ 16,n_width,n_height,'a2 Pt ',ytit,'LOG')
          CALL multitop(l+ 17,n_width,n_height,'a2 Abs(eta) ',ytit,'LOG')
          CALL multitop(l+ 18,n_width,n_height,'j1 Pt ',ytit,'LOG')
          CALL multitop(l+ 19,n_width,n_height,'j1 Abs(eta) ',ytit,'LOG')
          CALL multitop(l+ 20,n_width,n_height,'j2 Pt ',ytit,'LOG')
          CALL multitop(l+ 21,n_width,n_height,'j2 Abs(eta) ',ytit,'LOG')
       END DO
       CALL close_topdrawer_file_pp_aajj_dps
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file_pp_aajj_dps
       psfile="pp_aajj_dps.ps"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 2,'|Delta eta(a,a)| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 3,'|Delta phi(a,a)| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 4,'a a inv. mass ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 5,'a a trans mom ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 6,'max |Delta eta(a,j)| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 7,'max |Delta phi(a,j)| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 8,'max a j inv. mass ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 9,'max a j trans mom ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 10,'|Delta eta(j,j)| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 11,'|Delta phi(j,j)| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 12,'j j inv. mass ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 13,'j j trans mom ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 14,'a1 Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 15,'a1 Abs(eta)',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 16,'a2 Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 17,'a2 Abs(eta)',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 18,'j1 Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 19,'j1 Abs(eta)',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 20,'j2 Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 21,'j2 Abs(eta)',ytit,'LOG',psfile)
       ENDDO
       CALL close_gnuplot_file_pp_aajj_dps
    ENDIF
    IF(root_output)THEN
       CALL open_root_file_pp_aajj_dps
       rootfile="pp_aajj_dps.root"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 2,'|Delta eta(a,a)| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 3,'|Delta phi(a,a)| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 4,'a a inv. mass ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 5,'a a trans mom ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 6,'max |Delta eta(a,j)| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 7,'max |Delta phi(a,j)| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 8,'max a j inv. mass ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 9,'max a j trans mom ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 10,'|Delta eta(j,j)| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 11,'|Delta phi(j,j)| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 12,'j j inv. mass ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 13,'j j trans mom ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 14,'a1 Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 15,'a1 Abs(eta) ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 16,'a2 Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 17,'a2 Abs(eta) ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 18,'j1 Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 19,'j1 Abs(eta) ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 20,'j2 Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 21,'j2 Abs(eta) ',ytit,'LOG',rootfile)
       ENDDO
       ! wirte the end of the .C file for root
       WRITE (96, 100)
100    FORMAT(/1X,&
            ' hohisto -> cd();',/1X,&
            ' if (histos -> GetEntries() > 0 ) then {',/1X,&
            '  histos->Write();',/1X,&
            '  hohisto -> Close();',/1X,&
            ' }',/1X,'}')
       CALL close_root_file_pp_aajj_dps
    ENDIF
    RETURN
  END SUBROUTINE plot_end_pp_aajj_dps

  SUBROUTINE plot_fill_pp_aajj_dps(wgts)
    USE Helac_Global
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(2,4)::pmom_a,pmom_j
    REAL(KIND(1d0)),DIMENSION(4)::pa1,pa2,pj1,pj2,psum
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,q,e
    REAL(KIND(1d0))::exp1,exp2,pt1,pt2,pta1,pta2,ptj1,ptj2
    REAL(KIND(1d0))::etaa1,etaa2,etaj1,etaj2,maa,mjj,majmax
    REAL(KIND(1d0))::Ptaa,Ptjj,Ptajmax,phia1,phia2,phij1,phij2
    REAL(KIND(1d0))::Deltaetaaa,Deltaetajj,Deltaetaajmax
    REAL(KIND(1d0))::Deltaphiaa,Deltaphijj,Deltaphiajmax
    INTEGER::i,kk,l
    REAL(KIND(1d0))::www,yftot
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    var=1.0d0 ! total cross section
    ! the first two final states are photons
    pmom_a(1,1:4)=dps_hadron_pmom(5,1:4)
    pmom_a(2,1:4)=dps_hadron_pmom(6,1:4)
    ! the last two final states are jets
    pmom_j(1,1:4)=dps_hadron_pmom(7,1:4)
    pmom_j(2,1:4)=dps_hadron_pmom(8,1:4)
    ! sort photon into their pt
    pt1=pmom_a(1,1)**2+pmom_a(1,2)**2
    pt2=pmom_a(2,1)**2+pmom_a(2,2)**2
    IF(pt1.GE.pt2)THEN
       pa1(1:4)=pmom_a(1,1:4)
       pa2(1:4)=pmom_a(2,1:4)
    ELSE
       pa1(1:4)=pmom_a(2,1:4)
       pa2(1:4)=pmom_a(1,1:4)
    ENDIF
    psum(1:4)=pa1(1:4)+pa2(1:4)
    pta1=DSQRT(pa1(1)**2+pa1(2)**2)
    pta2=DSQRT(pa2(1)**2+pa2(2)**2)
    etaa1=prapidity(pa1)
    etaa2=prapidity(pa2)
    maa=psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2
    IF(maa.LE.0d0)THEN
       maa=0d0
    ELSE
       maa=DSQRT(maa)
    ENDIF
    phia1=ph4(pa1(1),pa1(2),pa1(3))
    phia2=ph4(pa2(1),pa2(2),pa2(3))
    Deltaetaaa=DABS(etaa1-etaa2)
    Ptaa=DSQRT(psum(1)**2+psum(2)**2)
    Deltaphiaa=DABS(phia1-phia2)
    ! sort jet into their pt
    pt1=pmom_j(1,1)**2+pmom_j(1,2)**2
    pt2=pmom_j(2,1)**2+pmom_j(2,2)**2
    IF(pt1.GE.pt2)THEN
       pj1(1:4)=pmom_j(1,1:4)
       pj2(1:4)=pmom_j(2,1:4)
    ELSE
       pj1(1:4)=pmom_j(2,1:4)
       pj2(1:4)=pmom_j(1,1:4)
    ENDIF
    psum(1:4)=pj1(1:4)+pj2(1:4)
    ptj1=DSQRT(pj1(1)**2+pj1(2)**2)
    ptj2=DSQRT(pj2(1)**2+pj2(2)**2)
    etaj1=prapidity(pj1)
    etaj2=prapidity(pj2)
    mjj=psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2
    IF(mjj.LE.0d0)THEN
       mjj=0d0
    ELSE
       mjj=DSQRT(mjj)
    ENDIF
    phij1=ph4(pj1(1),pj1(2),pj1(3))
    phij2=ph4(pj2(1),pj2(2),pj2(3))
    Deltaetajj=DABS(etaj1-etaj2)
    Ptjj=DSQRT(psum(1)**2+psum(2)**2)
    Deltaphijj=DABS(phij1-phij2)
    ! for a+j
    ! slect the max ones
    Ptajmax=0d0
    Deltaetaajmax=0d0
    Deltaphiajmax=0d0
    majmax=0d0
    DO i=1,2
       DO kk=1,2
          psum(1:4)=pmom_a(i,1:4)+pmom_j(kk,1:4)
          pt1=DSQRT(psum(1)**2+psum(2)**2)
          IF(pt1.GT.Ptajmax)Ptajmax=pt1
          pt1=prapidity(pmom_a(i,1:4))
          pt2=prapidity(pmom_j(kk,1:4))
          pt1=DABS(pt1-pt2)
          IF(pt1.GT.Deltaetaajmax)Deltaetaajmax=pt1
          pt1=ph4(pmom_a(i,1),pmom_a(i,2),pmom_a(i,3))
          pt2=ph4(pmom_j(kk,1),pmom_j(kk,2),pmom_j(kk,3))
          pt1=DABS(pt1-pt2)
          IF(pt1.GT.Deltaphiajmax)Deltaphiajmax=pt1
          pt1=psum(4)**2-psum(1)**2-psum(2)**2-psum(3)**2
          IF(pt1.LE.0)THEN
             pt1=0d0
          ELSE
             pt1=DSQRT(pt1)
          ENDIF
          IF(pt1.GT.majmax)majmax=pt1
       ENDDO
    ENDDO
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,Deltaetaaa,www)
       CALL mfill(l+3,Deltaphiaa,www)
       CALL mfill(l+4,maa,www)
       CALL mfill(l+5,Ptaa,www)
       CALL mfill(l+6,Deltaetaajmax,www)
       CALL mfill(l+7,Deltaphiajmax,www)
       CALL mfill(l+8,majmax,www)
       CALL mfill(l+9,Ptajmax,www)
       CALL mfill(l+10,Deltaetajj,www)
       CALL mfill(l+11,Deltaphijj,www)
       CALL mfill(l+12,mjj,www)
       CALL mfill(l+13,Ptjj,www)
       CALL mfill(l+14,pta1,www)
       CALL mfill(l+15,DABS(etaa1),www)
       CALL mfill(l+16,pta2,www)
       CALL mfill(l+17,DABS(etaa2),www)
       CALL mfill(l+18,ptj1,www)
       CALL mfill(l+19,DABS(etaj1),www)
       CALL mfill(l+20,ptj2,www)
       CALL mfill(l+21,DABS(etaj2),www)
    ENDDO
    ! boost back to c.m. frame
    RETURN
  END SUBROUTINE plot_fill_pp_aajj_dps

  SUBROUTINE open_topdrawer_file_pp_aajj_dps
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE=TRIM(output_dir)//'pp_aajj_dps.top',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_topdrawer_file_pp_aajj_dps

  SUBROUTINE close_topdrawer_file_pp_aajj_dps
    IMPLICIT NONE
    CLOSE(99)
    RETURN
  END SUBROUTINE close_topdrawer_file_pp_aajj_dps

  SUBROUTINE open_gnuplot_file_pp_aajj_dps
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=97,FILE=TRIM(output_dir)//'pp_aajj_dps.gnu',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_gnuplot_file_pp_aajj_dps

  SUBROUTINE close_gnuplot_file_pp_aajj_dps
    IMPLICIT NONE
    CLOSE(97)
    RETURN
  END SUBROUTINE close_gnuplot_file_pp_aajj_dps

  SUBROUTINE open_root_file_pp_aajj_dps
    USE Helac_Global
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=96,FILE=TRIM(output_dir)//'pp_aajj_dps.C',&
         status='unknown')
    useitmax=.FALSE.
  END SUBROUTINE open_root_file_pp_aajj_dps

  SUBROUTINE close_root_file_pp_aajj_dps
    IMPLICIT NONE
    CLOSE(96)
    RETURN
  END SUBROUTINE close_root_file_pp_aajj_dps

END MODULE plot_pp_aajj_dps
