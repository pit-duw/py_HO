MODULE plot_user
!============================================================
!
! PLOT FILE FOR USER
!
!============================================================
  USE Helac_Global
  USE Kinetic_Func ! provides kinematic functions
  IMPLICIT NONE
  INTEGER::nwgt_analysis
  INTEGER::num_plots=14
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
       CALL bookup(l+ 2,'|Delta y|  '//weights_info(kk),&
            0.05d0,0d0,5d0)
       CALL bookup(l+ 3,'|Delta phi| '//weights_info(kk),&
            3.1417d0/100d0,-0.0001d0,3.1416d0)
       CALL bookup(l+ 4,'psi psi inv. mass '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 5,'psi psi inv. mass zoom '//weights_info(kk),&
            0.05d0,5.5d0,10.5d0)
       CALL bookup(l+ 6,'psi psi trans mom '//weights_info(kk),&
            1d0,0d0,100d0)
       ! q12^2 = (Ejpsi1+Ejpsi2)^2-(pvecjpsi1+pvecjpsi2)^2-(2mjpsi)^2
       CALL bookup(l+ 7,'q12 '//weights_info(kk),&
            0.5d0,0d0,50d0)
       CALL bookup(l+ 8,'q12 zoom '//weights_info(kk),&
            0.04d0,0d0,4d0)
       !CALL bookup(l+ 9,'Cos(alpha) '//weights_info(kk),&
       !     0.01d0,0d0,1d0)
       CALL bookup(l+ 9,'|Delta eta| '//weights_info(kk),&
            0.04d0,0d0,4d0)
       CALL bookup(l+ 10,'random Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 11,'leading Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 12,'subleading Pt '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 13,'Pt diff '//weights_info(kk),&
            1d0,0d0,100d0)
       CALL bookup(l+ 14,'|Pt1-Pt2| '//weights_info(kk),&
            1d0,0d0,100d0)
       !CALL bookup_3d(l+ 16,'Cos(th1)-Cos(th2) '//weights_info(kk),&
       !     0.02d0,-1d0,1d0,0.02d0,-1d0,1d0)
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
    CALL mclear_3d
    DO i=1,NPLOTS
       CALL mopera(i,'+',i,i,xnorm,0.d0)
       CALL mfinal(i)
       CALL mopera_3d(i,'+',i,i,xnorm,0.d0)
       CALL mfinal_3d(i)
    ENDDO
    ytit='sigma per bin '
    IF(topdrawer_output)THEN
       CALL open_topdrawer_file
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL multitop(l+ 1,n_width,n_height,'total rate   ',ytit,'LIN')
          CALL multitop(l+ 2,n_width,n_height,'|Delta y|  ',ytit,'LOG')
          CALL multitop(l+ 3,n_width,n_height,'|Delta phi| ',ytit,'LOG')
          CALL multitop(l+ 4,n_width,n_height,'psi psi inv. mass ',ytit,'LOG')
          CALL multitop(l+ 5,n_width,n_height,'psi psi inv. mass zoom ',ytit,'LOG')
          CALL multitop(l+ 6,n_width,n_height,'psi psi trans mom ',ytit,'LOG')
          CALL multitop(l+ 7,n_width,n_height,'q12 ',ytit,'LOG')
          CALL multitop(l+ 8,n_width,n_height,'q12 zoom ',ytit,'LOG')
          !CALL multitop(l+ 9,n_width,n_height,'Cos(alpha) ',ytit,'LIN')
          CALL multitop(l+ 9,n_width,n_height,'|Delta eta| ',ytit,'LOG')
          CALL multitop(l+ 10,n_width,n_height,'random Pt ',ytit,'LOG')
          CALL multitop(l+ 11,n_width,n_height,'leading Pt ',ytit,'LOG')
          CALL multitop(l+ 12,n_width,n_height,'subleading Pt ',ytit,'LOG')
          CALL multitop(l+ 13,n_width,n_height,'Pt diff ',ytit,'LOG')
          CALL multitop(l+ 14,n_width,n_height,'|Pt1-Pt2| ',ytit,'LOG')
       END DO
       CALL close_topdrawer_file
    ENDIF
    IF(gnuplot_output)THEN
       CALL open_gnuplot_file
       psfile=TRIM(process(2:10*nhad+1))//".ps"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 2,'|Delta y| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 3,'|Delta phi| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 4,'psi psi inv. mass ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 5,'psi psi inv. mass zoom ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 6,'psi psi trans mom ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 7,'q12 ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 8,'q12 zoom ',ytit,'LOG',psfile)
          !CALL MGNUPLOT(l+ 9,'Cos(alpha) ',ytit,'LIN',psfile)
          CALL MGNUPLOT(l+ 9,'|Delta eta| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 10,'random Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 11,'leading Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 12,'subleading Pt ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 13,'Pt diff ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+ 14,'|Pt1-Pt2| ',ytit,'LOG',psfile)
          !CALL MGNUPLOT_3D(l+ 16,'Cos(th1)','Cos(th2)',ytit,'LIN',psfile)
       ENDDO
       CALL close_gnuplot_file
    ENDIF
    IF(root_output)THEN
       CALL open_root_file
       rootfile=TRIM(process(2:10*nhad+1))//".root"
       DO kk=1,nwgt_analysis
          l=(kk-1)*num_plots
          CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 2,'|Delta y| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 3,'|Delta phi| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 4,'psi psi inv. mass ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 5,'psi psi inv. mass zoom ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 6,'psi psi trans mom ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 7,'q12 ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 8,'q12 zoom ',ytit,'LOG',rootfile)
          !CALL MROOTPLOT(l+ 9,'Cos(alpha) ',ytit,'LIN',rootfile)
          CALL MROOTPLOT(l+ 9,'|Delta eta| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 10,'random Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 11,'leading Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 12,'subleading Pt ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 13,'Pt diff ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+ 14,'|Pt1-Pt2| ',ytit,'LOG',rootfile)
          !CALL MROOTPLOT_3D(l+ 16,'Cos(th1)','Cos(th2)',ytit,'LIN',rootfile)
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
    USE Helac_Func_1
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(20,4)::phad
    REAL(KIND(1d0)),DIMENSION(:),INTENT(IN)::wgts
    REAL(KIND(1d0))::var,Q,www,dy_psipsi,dphi_psipsi,M_psipsi,Pt_psipsi,deta_psipsi,Pt1Pt2_psipsi
    REAL(KIND(1d0))::q12_psipsi,ptpsi1,ptpsi2,ptpsimax,ptpsimin,ptpsidiff,e
    REAL(KIND(1d0))::exp1,exp2
    INTEGER::i,m,l,kk
    REAL(KIND(1d0)),DIMENSION(4)::PBOO
    REAL(KIND(1d0)),DIMENSION(2,4)::PJpsi,PJpsiCMS
    !REAL(KIND(1d0)),DIMENSION(4,4)::Plep
    !REAL(KIND(1d0)),DIMENSION(4)::plep1,Plep2
    !REAL(KIND(1d0)),DIMENSION(3)::Pplain1,Pplain2
    !REAL(KIND(1d0)),DIMENSION(3)::Plepvec1,Plepvec2
    !REAL(KIND(1d0)),DIMENSION(2,3)::PJpsi1l,PJpsi2l
    !REAL(KIND(1d0))::PP1,PP2
    !REAL(KIND(1d0))::costh1,costh2
    !INTEGER,DIMENSION(4)::lep2jpsi
    INTEGER,DIMENSION(2)::newmoth
    INTEGER::i_jpsi,i_lep,kk1,kk2
    INTEGER::definemoth
    REAL(KIND(1d0)),PARAMETER::pi=3.1415926535897932384626433832795d0
    !IF(NDecayChains.GT.0)THEN
    !   CALL plot_fill_Decay(wgts)
    !   RETURN
    !ENDIF
    !IF(emep_ISR_shower.EQ.1)THEN
    !   CALL plot_fill_QEDPS(wgts)
    !   RETURN
    !ENDIF
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
    PJpsi(1,1:4)=phad(3,1:4) ! (px,py,pz,E)
    PJpsi(2,1:4)=phad(4,1:4)
    IF(literature_cutoffs.EQ.14080000)THEN
       ! LHCb Kinematics Region
       ptpsi1=transverse(PJpsi(1,1:4))
       IF(ptpsi1.GT.10d0)THEN
          RETURN
       ENDIF
       ptpsi2=transverse(PJpsi(2,1:4))
       IF(ptpsi2.GT.10d0)THEN
          RETURN
       ENDIF
    ENDIF
    dy_psipsi=DABS(rapidity(PJpsi(1,1:4))-rapidity(PJpsi(2,1:4)))
    deta_psipsi=DABS(prapidity(PJpsi(1,1:4))-prapidity(PJpsi(2,1:4)))
    dphi_psipsi=ph4(PJpsi(1,1),PJpsi(1,2),PJpsi(1,3))
    dphi_psipsi=dphi_psipsi-ph4(PJpsi(2,1),PJpsi(2,2),PJpsi(2,3))
    dphi_psipsi=DABS(dphi_psipsi)
    IF(dphi_psipsi.GT.pi)dphi_psipsi=2*pi-dphi_psipsi
    M_psipsi=scalar_product(PJpsi(1,1:4)+PJpsi(2,1:4),PJpsi(1,1:4)+PJpsi(2,1:4))
    M_psipsi=DSQRT(M_psipsi)
    Q=scalar_product(PJpsi(1,1:4),PJpsi(1,1:4))
    Q=DSQRT(Q)
    ptpsi1=transverse(PJpsi(1,1:4))
    ptpsi2=transverse(PJpsi(2,1:4))
    ptpsimax=MAX(ptpsi1,ptpsi2)
    ptpsimin=MIN(ptpsi1,ptpsi2)
    ptpsidiff=ptpsimax-ptpsimin
    Pt_psipsi=transverse(PJpsi(1,1:4)+PJpsi(2,1:4))
    Pt1Pt2_psipsi=transverse(PJpsi(1,1:4)-PJpsi(2,1:4))
    q12_psipsi=M_psipsi**2-4*Q**2
    q12_psipsi=DSQRT(DABS(q12_psipsi))
    !PBOO(1:3)=-PJpsi(1,1:3)-PJpsi(2,1:3)
    !PBOO(4)=PJpsi(1,4)+PJpsi(2,4)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! START TO FILL HISTOGRAM
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,dy_psipsi,www)
       CALL mfill(l+3,dphi_psipsi,www)
       CALL mfill(l+4,M_psipsi,www)
       CALL mfill(l+5,M_psipsi,www)
       CALL mfill(l+6,Pt_psipsi,www)
       CALL mfill(l+7,q12_psipsi,www)
       CALL mfill(l+8,q12_psipsi,www)
       !CALL mfill(l+9,cosalpha,www)
       CALL mfill(l+9,deta_psipsi,www)
       CALL mfill(l+10,ptpsi1,www)
       CALL mfill(l+11,ptpsimax,www)
       CALL mfill(l+12,ptpsimin,www)
       CALL mfill(l+13,ptpsidiff,www)
       CALL mfill(l+14,Pt1Pt2_psipsi,www)
       !CALL mfill_3d(l+16,costh1,costh2,www)
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
    INCLUDE "../../decay/stdhep.inc"! decay variables are defined in stdhep.inc
    REAL(KIND(1d0))::var,Q,www,dy_psipsi,dphi_psipsi,M_psipsi,Pt_psipsi,deta_psipsi,Pt1Pt2_psipsi
    REAL(KIND(1d0))::q12_psipsi,cosalpha,ptpsi1,ptpsi2,ptpsimax,ptpsimin,ptpsidiff
    INTEGER::i,m,l,kk
    REAL(KIND(1d0)),DIMENSION(4)::PBOO
    REAL(KIND(1d0)),DIMENSION(2,4)::PJpsi,PJpsiCMS
    REAL(KIND(1d0)),DIMENSION(4,4)::Plep
    REAL(KIND(1d0)),DIMENSION(4)::plep1,Plep2
    REAL(KIND(1d0)),DIMENSION(3)::Pplain1,Pplain2
    REAL(KIND(1d0)),DIMENSION(3)::Plepvec1,Plepvec2
    REAL(KIND(1d0)),DIMENSION(2,3)::PJpsi1l,PJpsi2l
    REAL(KIND(1d0))::PP1,PP2
    REAL(KIND(1d0))::costh1,costh2
    INTEGER,DIMENSION(4)::lep2jpsi
    INTEGER,DIMENSION(2)::newmoth
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
    lep2jpsi(1:4)=0
    newmoth(1:2)=0
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
          IF(i_lep.EQ.4)EXIT
       ENDIF
    ENDDO
    IF(literature_cutoffs.EQ.14080000)THEN
       ! LHCb Kinematics Region
       ptpsi1=transverse(PJpsi(1,1:4))
       IF(ptpsi1.GT.10d0)THEN
          RETURN
       ENDIF
       ptpsi2=transverse(PJpsi(2,1:4))
       IF(ptpsi2.GT.10d0)THEN
          RETURN
       ENDIF
    ENDIF
    dy_psipsi=DABS(rapidity(PJpsi(1,1:4))-rapidity(PJpsi(2,1:4)))
    deta_psipsi=DABS(prapidity(PJpsi(1,1:4))-prapidity(PJpsi(2,1:4)))
    dphi_psipsi=ph4(PJpsi(1,1),PJpsi(1,2),PJpsi(1,3))
    dphi_psipsi=dphi_psipsi-ph4(PJpsi(2,1),PJpsi(2,2),PJpsi(2,3))
    dphi_psipsi=DABS(dphi_psipsi)
    IF(dphi_psipsi.GT.pi)dphi_psipsi=2*pi-dphi_psipsi
    M_psipsi=scalar_product(PJpsi(1,1:4)+PJpsi(2,1:4),PJpsi(1,1:4)+PJpsi(2,1:4))
    M_psipsi=DSQRT(M_psipsi)
    Q=scalar_product(PJpsi(1,1:4),PJpsi(1,1:4))
    Q=DSQRT(Q)
    ptpsi1=transverse(PJpsi(1,1:4))
    ptpsi2=transverse(PJpsi(2,1:4))
    ptpsimax=MAX(ptpsi1,ptpsi2)
    ptpsimin=MIN(ptpsi1,ptpsi2)
    ptpsidiff=ptpsimax-ptpsimin
    Pt_psipsi=transverse(PJpsi(1,1:4)+PJpsi(2,1:4))
    Pt1Pt2_psipsi=transverse(PJpsi(1,1:4)-PJpsi(2,1:4))
    q12_psipsi=M_psipsi**2-4*Q**2
    q12_psipsi=DSQRT(DABS(q12_psipsi))
    PBOO(1:3)=-PJpsi(1,1:3)-PJpsi(2,1:3)
    PBOO(4)=PJpsi(1,4)+PJpsi(2,4)
    DO i=1,4
       CALL BOOSTL(M_psipsi,PBOO,Plep(i,1:4))
    ENDDO
    DO i=1,2
       PJpsiCMS(i,1:4)=PJpsi(i,1:4)
       CALL BOOSTL(M_psipsi,PBOO,PJpsiCMS(i,1:4))
    ENDDO
    kk1=0
    kk2=0
    DO i=1,4
       IF(lep2jpsi(i).EQ.1)THEN
          kk1=kk1+1
          PJpsi1l(kk1,1:3)=Plep(i,1:3)
          IF(kk1.EQ.1)THEN
             Plep1(1:4)=Plep(i,1:4)
          ENDIF
       ELSEIF(lep2jpsi(i).EQ.2)THEN
          kk2=kk2+1
          PJpsi2l(kk2,1:3)=Plep(i,1:3)
          IF(kk2.EQ.1)THEN
             Plep2(1:4)=Plep(i,1:4)
          ENDIF
       ENDIF
    ENDDO
    Pplain1(1:3)=vec_crossprod(PJpsi1l(1,1:3),PJpsi1l(2,1:3))
    Pplain2(1:3)=vec_crossprod(PJpsi2l(1,1:3),PJpsi2l(2,1:3))
    PP1=Pplain1(1)**2+Pplain1(2)**2+Pplain1(3)**2
    PP1=DSQRT(PP1)
    PP2=Pplain2(1)**2+Pplain2(2)**2+Pplain2(3)**2
    PP2=DSQRT(PP2)
    IF(PP1*PP2.EQ.0d0)THEN
       cosalpha=0d0
    ELSE
       cosalpha=Pplain1(1)*Pplain2(1)+Pplain1(2)*Pplain2(2) &
       +Pplain1(3)*Pplain2(3)
       cosalpha=DABS(cosalpha/(PP1*PP2))
    ENDIF
    ! costh1
    Q=scalar_product(PJpsiCMS(1,1:4),PJpsiCMS(1,1:4))
    Q=DSQRT(Q)
    PBOO(1:3)=-PJpsiCMS(1,1:3)
    PBOO(4)=PJpsiCMS(1,4)
    CALL BOOSTL(Q,PBOO,Plep1)
    costh1=cosij(Plep1(1:3),PJpsiCMS(1,1:3))
    ! costh2
    Q=scalar_product(PJpsiCMS(2,1:4),PJpsiCMS(2,1:4))
    Q=DSQRT(Q)
    PBOO(1:3)=-PJpsiCMS(2,1:3)
    PBOO(4)=PJpsiCMS(2,4)
    CALL BOOSTL(Q,PBOO,Plep2)
    costh2=cosij(Plep2(1:3),PJpsiCMS(2,1:3))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    !
    ! START TO FILL HISTOGRAM 
    !  
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    DO kk=1,nwgt_analysis
       www=wgts(kk)
       l=(kk-1)*num_plots
       CALL mfill(l+1,var,www)
       CALL mfill(l+2,dy_psipsi,www)
       CALL mfill(l+3,dphi_psipsi,www)
       CALL mfill(l+4,M_psipsi,www)
       CALL mfill(l+5,M_psipsi,www)
       CALL mfill(l+6,Pt_psipsi,www)
       CALL mfill(l+7,q12_psipsi,www)
       CALL mfill(l+8,q12_psipsi,www)
       CALL mfill(l+9,cosalpha,www)
       CALL mfill(l+10,deta_psipsi,www)
       CALL mfill(l+11,ptpsi1,www)
       CALL mfill(l+12,ptpsimax,www)
       CALL mfill(l+13,ptpsimin,www)
       CALL mfill(l+14,ptpsidiff,www)
       CALL mfill(l+15,Pt1Pt2_psipsi,www)
       CALL mfill_3d(l+16,costh1,costh2,www)
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
