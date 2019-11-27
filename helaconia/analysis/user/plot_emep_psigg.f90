MODULE plot_user
!============================================================
!
! PLOT FILE FOR USER
! SPECIAL EXAMPLE FOR e- e+ > psi+gg
!
!============================================================
  USE Helac_Global
  USE Kinetic_Func ! provides kinematic functions
  IMPLICIT NONE
  INTEGER::nwgt_analysis
  INTEGER::num_plots=3
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
    IF(emep_ISR_shower.EQ.1)num_plots=num_plots+2
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
       CALL bookup(l+ 2,'total partonic y rap  '//weights_info(kk),&
            0.2d0,-5d0,5d0)
       CALL bookup(l+ 3,'momentum of psi ',0.05d0,0d0,5d0)
       IF(emep_ISR_shower.EQ.1)THEN
          CALL bookup(l+ 4,'number of photon '//weights_info(kk),&
               1d0,-0.5d0,11.5d0)
          CALL bookup(l+ 5,'fraction of hard energy '//weights_info(kk),&
               0.02d0,-0.01d0,1.01d0)
       ENDIF
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
          CALL multitop(l+ 2,n_width,n_height,'total partonic y rap  ',ytit,'LOG')
          CALL multitop(l+ 3,n_width,n_height,'momentum of psi ',ytit,'LIN')
          IF(emep_ISR_shower.EQ.1)THEN
             CALL multitop(l+ 4,n_width,n_height,'number of photon ',ytit,'LOG')
             CALL multitop(l+ 5,n_width,n_height,'fraction of hard energy ',ytit,'LOG')
          ENDIF
       END DO
       CALL close_topdrawer_file
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file
       psfile=TRIM(process(2:10*nhad+1))//".ps"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 2,'total partonic y rap  ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 3,'momentum of psi ',ytit,'LIN',psfile)
          IF(emep_ISR_shower.EQ.1)THEN
             CALL MGNUPLOT(l+ 4,'number of photon ',ytit,'LOG',psfile)
             CALL MGNUPLOT(l+ 5,'fraction of hard energy ',ytit,'LOG',psfile)
          ENDIF
       ENDDO
       CALL close_gnuplot_file
    ENDIF
    IF(root_output)THEN
       CALL open_root_file
       rootfile=TRIM(process(2:10*nhad+1))//".root"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 2,'total partonic y rap  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 3,'momentum of psi ',ytit,'LIN',rootfile)
          IF(emep_ISR_shower.EQ.1)THEN
             CALL MROOTPLOT(l+ 4,'number of photon ',ytit,'LOG',rootfile)
             CALL MROOTPLOT(l+ 5,'fraction of hard energy ',ytit,'LOG',rootfile)
          ENDIF
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
    REAL(KIND(1d0))::ppsi
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,kk,l
    REAL(KIND(1d0))::www,yftot
    IF(emep_ISR_shower.EQ.1)THEN
       CALL plot_fill_QEDPS(wgts)
       RETURN
    ENDIF
    ! WARNING*WARNING*WARNING*WARNING*WARNING*WARNING*WARNING
    !
    ! 1)The generated momenta phegas_pmom is in the partonic c.m. frame when imode.EQ.0
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
    pini(1:4)=phad(1,1:4)+phad(2,1:4) ! (px,py,pz,E)
    yftot=rapidity(pini)
    ppsi=phad(3,1)**2+phad(3,2)**2+phad(3,3)**2
    ppsi=DSQRT(ppsi)
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
       CALL mfill(l+3,ppsi,www)
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
    REAL(KIND(1d0))::exp1,exp2,ppsi
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
    ppsi=0d0
    DO i=1,NPTCL
       IF(NLPTN(1,i).EQ.22)THEN
          nphoton=nphoton+1
       ELSEIF(NLPTN(1,i).EQ.443)THEN
          ppsi=PLPTN(4,i)**2+PLPTN(5,i)**2+PLPTN(6,i)**2
          ppsi=DSQRT(ppsi)
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
       CALL mfill(l+3,ppsi,www)
       CALL mfill(l+4,DBLE(nphoton),www)
       CALL mfill(l+5,frache,www)
    ENDDO
    RETURN
  END SUBROUTINE plot_fill_QEDPS
END MODULE plot_user
