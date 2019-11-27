MODULE plot_MLM_py8_user
!============================================================
!
! PLOT FILE WITH PYTHIA8 and MLM merging FOR USER
! This is an example file for z+njet p p > z+n-jet
!
!============================================================
  USE MLM_merge
  IMPLICIT NONE
  INTEGER::nwgt_analysis
  INTEGER::num_plots=61
  INTEGER,PARAMETER::N_JET_ALG=3
  INTEGER,PARAMETER::n_width=1,n_height=1
! set max_multiplicity=.true. if the current sample is that of largest
! multiplicity, and equal to .false. otherwise
  LOGICAL,PARAMETER::max_multiplicity=.FALSE.
! matching_scale: this is the \mu_Q of the paper
  REAL(KIND(1d0)),PARAMETER::matching_scale=50d0
! n_match identified the i-parton sample of the paper: that's the number
! of partons AT THE BORN LEVEL in a given multiplicity
  INTEGER,PARAMETER::n_match=1
CONTAINS
  SUBROUTINE plot_MLM_py8_begin(nwgt,weights_info)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nwgt
    INTEGER,PARAMETER::max_weight=1
    CHARACTER(len=15),DIMENSION(max_weight),INTENT(IN)::weights_info
    CHARACTER(len=15)::weights_info2
    CHARACTER(len=7),DIMENSION(N_JET_ALG)::cc=(/'kt 0.5 ','akt 0.5','akt 1.0'/)
    INTEGER::i,kk,l,j,k
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
       CALL bookup(l+ 1,'total rate    ',&
            1.0d0,0.5d0,5.5d0)
       CALL bookup(l+ 2,'z pt  ',&
            10d0,0d0,1000d0)
       CALL bookup(l+ 3,'z |y| ',&
            0.2d0,0d0,2.2d0)
       CALL bookup(l+ 4,'z y ',&
            0.2d0,-4d0,4d0)
       DO j=1,N_JET_ALG
          k=(j-1)*19
          CALL bookup(l+k+5,'# jets '//cc(j),&
               1d0,-0.5d0,10d0)
          CALL bookup(l+k+6,'pt[j1] '//cc(j),&
               1d0,0d0,100d0)
          CALL bookup(l+k+7,'pt[j2] '//cc(j),&
               1d0,0d0,100d0)
          CALL bookup(l+k+8,'pt[j3] '//cc(j),&
               1d0,0d0,100d0)
          CALL bookup(l+k+9,'pt[j4] '//cc(j),&
               1d0,0d0,100d0)
          CALL bookup(l+k+10,'|y[j1]| '//cc(j),&
               0.2d0,0d0,2.4d0)
          CALL bookup(l+k+11,'|y[j2]| '//cc(j),&
               0.2d0,0d0,2.4d0)
          CALL bookup(l+k+12,'|y[j3]| '//cc(j),&
               0.2d0,0d0,2.4d0)
          CALL bookup(l+k+13,'|y[j4]| '//cc(j),&
               0.2d0,0d0,2.4d0)
! the smoothness of d0 in z+0jet and z+1jet is a signal of goodness of the merging
          CALL bookup(l+k+14,'d0      '//cc(j),&
               5d0,0.d0,500.d0)
          CALL bookup(l+k+15,'log10[d0] '//cc(j),&
               0.05d0,-0.5d0,3.5d0)
          CALL bookup(l+k+16,'d1       '//cc(j),&
               5d0,0.d0,500.d0)
          CALL bookup(l+k+17,'log10[d1] '//cc(j),&
               0.05d0,-0.5d0,3.5d0)
          CALL bookup(l+k+18,'d2       '//cc(j),&
               5d0,0.d0,500.d0)
          CALL bookup(l+k+19,'log10[d2] '//cc(j),&
               0.05d0,-0.5d0,3.5d0)
          CALL bookup(l+k+20,'d3       '//cc(j),&
               5d0,0.d0,500.d0)
          CALL bookup(l+k+21,'log10[d3] '//cc(j),&
               0.05d0,-0.5d0,3.5d0)
! see arXiv:1310.3082
! ydiff=|y[ee]-y[j1]|/2
          CALL bookup(l+k+22,'ydiff, 30 '//cc(j),&
               0.2d0,0.0d0,1.8d0)
! ysum=|y[ee]+y[j1]|/2
          CALL bookup(l+k,23,'ysum, 30 '//cc(j),&
               0.2d0,0.0d0,2.2d0)
    ENDDO
    RETURN
  END SUBROUTINE plot_MLM_py8_begin

  SUBROUTINE plot_MLM_py8_end(xnorm)
    IMPLICIT NONE
    CHARACTER(len=14)::ytit
    REAL(KIND(1d0)),INTENT(IN)::xnorm
    INTEGER::i
    INTEGER::kk,l,k,j
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
       CALL multitop(l+ 2,n_width,n_height,'z pt  ',ytit,'LOG')
       CALL multitop(l+ 3,n_width,n_height,'z |y| ',ytit,'LOG')
       CALL multitop(l+ 4,n_width,n_height,'z y   ',ytit,'LOG')
       DO j=1,N_JET_ALG
          k=(j-1)*19
          CALL multitop(l+k+5,n_width,n_height,'# jets  ',ytit,'LOG')
          CALL multitop(l+k+6,n_width,n_height,'pt[j1]  ',ytit,'LOG')
          CALL multitop(l+k+7,n_width,n_height,'pt[j2]  ',ytit,'LOG')
          CALL multitop(l+k+8,n_width,n_height,'pt[j3]  ',ytit,'LOG')
          CALL multitop(l+k+9,n_width,n_height,'pt[j4]  ',ytit,'LOG')
          CALL multitop(l+k+10,n_width,n_height,'|y[j1]|',ytit,'LOG')
          CALL multitop(l+k+11,n_width,n_height,'|y[j2]|',ytit,'LOG')
          CALL multitop(l+k+12,n_width,n_height,'|y[j3]|',ytit,'LOG')
          CALL multitop(l+k+13,n_width,n_height,'|y[j4]|',ytit,'LOG')
          CALL multitop(l+k+14,n_width,n_height,'d0     ',ytit,'LOG')
          CALL multitop(l+k+15,n_width,n_height,'log10[d0] ',ytit,'LOG')
          CALL multitop(l+k+16,n_width,n_height,'d1     ',ytit,'LOG')
          CALL multitop(l+k+17,n_width,n_height,'log10[d1] ',ytit,'LOG')
          CALL multitop(l+k+18,n_width,n_height,'d2     ',ytit,'LOG')
          CALL multitop(l+k+19,n_width,n_height,'log10[d2] ',ytit,'LOG')
          CALL multitop(l+k+20,n_width,n_height,'d3     ',ytit,'LOG')
          CALL multitop(l+k+21,n_width,n_height,'log10[d3] ',ytit,'LOG')
          CALL multitop(l+k+22,n_width,n_height,'ydiff, 30 ',ytit,'LOG')
          CALL multitop(l+k+23,n_width,n_height,'ysum, 30 ',ytit,'LOG')
       ENDDO
    END DO
    CALL close_topdrawer_file
    CALL open_gnuplot_file
    psfile="helaconia_MLM_py8.ps"
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL MGNUPLOT(l+ 1,'total rate   ',ytit,'LIN',psfile)
       CALL MGNUPLOT(l+ 2,'z pt  ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 3,'z |y| ',ytit,'LOG',psfile)
       CALL MGNUPLOT(l+ 4,'z  y  ',ytit,'LOG',psfile)
       DO j=1,N_JET_ALG
          k=(j-1)*19
          CALL MGNUPLOT(l+k+5,'# jets ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+6,'pt[j1] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+7,'pt[j2] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+8,'pt[j3] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+9,'pt[j4] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+10,'|y[j1]| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+11,'|y[j2]| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+12,'|y[j3]| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+13,'|y[j4]| ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+14,'d0      ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+15,'log10[d0] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+16,'d1      ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+17,'log10[d1] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+18,'d2      ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+19,'log10[d2] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+20,'d3      ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+21,'log10[d3] ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+22,'ydiff, 30 ',ytit,'LOG',psfile)
          CALL MGNUPLOT(l+k+23,'ysum, 30 ',ytit,'LOG',psfile)
       ENDDO
    ENDDO
    CALL close_gnuplot_file
    CALL open_root_file
    rootfile="helaconia_MLM_py8.root"
    DO kk=1,nwgt_analysis
       l=(kk-1)*num_plots
       CALL MROOTPLOT(l+ 1,'total rate  ',ytit,'LIN',rootfile)
       CALL MROOTPLOT(l+ 2,'z pt  ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 3,'z |y| ',ytit,'LOG',rootfile)
       CALL MROOTPLOT(l+ 4,'z  y  ',ytit,'LOG',rootfile)
       DO j=1, N_JET_ALG
          k=(j-1)*19
          CALL MROOTPLOT(l+k+5,'# jets ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+6,'pt[j1]  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+7,'pt[j2]  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+8,'pt[j3]  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+9,'pt[j4]  ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+10,'|y[j1]| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+11,'|y[j2]| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+12,'|y[j3]| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+13,'|y[j4]| ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+14,'d0      ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+15,'log10[d0] ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+16,'d1      ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+17,'log10[d1] ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+18,'d2      ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+19,'log10[d2] ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+20,'d3      ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+21,'log10[d3] ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+22,'ydiff, 30 ',ytit,'LOG',rootfile)
          CALL MROOTPLOT(l+k+23,'ysum, 30  ',ytit,'LOG',rootfile)
       ENDDO
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
  END SUBROUTINE plot_MLM_py8_end

  SUBROUTINE plot_MLM_py8_fill(wgts)
    IMPLICIT NONE
    !INCLUDE '../include/HEPMC90.INC'
    INCLUDE 'HEPMC90.INC' ! it will be copied to shower directory
    REAL(KIND(1d0)),DIMENSION(4)::psyst
    REAL(KIND(1d0)),DIMENSION(250),INTENT(IN)::wgts
    REAL(KIND(1d0))::wgt,var,yz,ptz,ptj,etaj,yj
    REAL(KIND(1d0))::d01,d12,d23,d34,ptj1,yj1,ptj2,yj2,ptj3,yj3,ptj4,yj4
    INTEGER::i,j,l,k,kk,IHEP,NPart,NTRACK,IFZ,NJET30
    REAL(KIND(1d0))::www,ydiff,ysum
    INTEGER,PARAMETER::MAXTRACK=2048,MAXJET=2048,MAXNUM=30
    REAL(KIND(1d0)),DIMENSION(4)::P_Z
    REAL(KIND(1d0)),DIMENSION(4,MAXTRACK)::PJET,PTRACK,P_PART
    INTEGER,DIMENSION(MAXTRACK)::JETVEC
    REAL(KIND(1d0)),EXTERNAL::fastjetdmergemax
    ! INCOMING PARTONS MAY TRAVEL IN THE SAME DIRECTION: IT''S 
    ! A POWER-SUPPRESSED EFFECT, SO THROW THE EVENT AWAY
    IF(SIGN(1.D0,PHEP(3,1)).EQ.SIGN(1.D0,PHEP(3,2)))THEN
       WRITE(*,*)"WARNING 1 in plot_MLM_py8_fill"
       GOTO 999
    ENDIF
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! IMPORT INFORMATION AND DEFINE JETS
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    NTRACK=0
    NPART=0
    IFZ=0
    DO IHEP=1,NHEP
       IST=ISTHEP(IHEP) ! status
       ID=IDHEP(IHEP) ! pdg

       ! z boson
       IF(IST.EQ.22.AND.ID.EQ.23)THEN
          IFZ=IFZ+1
          DO i=1,4
             P_Z(i,1)=PHEP(i,IHEP)
          ENDDO
       ENDIF

       ! define particles that to into jet.
       IF(IST.EQ.1.AND.ABS(ID).GE.100)THEN ! final + hadron
          NTRACKS=NTRACKS+1
          DO i=1,4
             PTRACK(i,NTRACKS)=PHEP(i,IHEP)
          ENDDO
          IF(NTRACKS.EQ.MAXTRACK)THEN
             WRITE(*,*)"ERROR: Too many particles in plot_MLM_py8_fill, please increase MATRACK !"
             STOP
          ENDIF
       ENDIF
       ! Keep QCD partons from primary shower
       ! see PYTHIA8 status code
       ! http://home.thep.lu.se/~torbjorn/pythia81html/Welcome.html
       ! Study Output -> Particle Properties
       ! 11-19 beam particles
       ! 21-29 particles of the hardest subprocess
       !   21 incoming
       !   22 intermediate
       !   23 outgoing
       !   24 outgoing, nonperturbatively kicked out in diffraction
       ! 31-39 particles of subsequent subprocesses
       ! 41-49 particles produced by initial-state-showers
       ! 51-59 particles produced by final-state-showers
       ! 61-69 particles produced by beam-remnant treatment
       ! 71-79 partons in preparation of hadronization process
       ! 81-89 primary hadrons produced by hadronization process
       ! 91-99 particles produced in decay process, or by Bose-Einstein effects
       ! 101-109 particles in the handling of R-hadron production and decay
       IF(IST.GE.71..AND.IST.LE.79.AND.((ABS(ID).GE.1.AND.ABS(ID).LE.5).OR.&
            ID.EQ.21))THEN ! five light flavor
          NPART=NPART+1
          DO I=1,4
             PPART(I,NPART)=PHEP(I,IHEP)
          ENDDO
          IF(NPART.EQ.MAXTRACK)THEN
             WRITE(*,*)"ERROR: Too many parton particles (after showering before hadronization) in plot_MLM_py8_fill, please increase MATRACK !"
             STOP
          ENDIF
       ENDIF
    ENDDO

    IF(IFZ.NE.1)THEN
       WRITE(*,*)"ERROR: The number of Z boson is not 1"
       STOP
    ENDIF
    ! rejecting events via MLM
    CALL MLM_matching(PPART,NPART,matching_scale,n_match,&
         max_multiplicity,matched)
    IF(.NOT.matched)RETURN
    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !
    ! DEFINE OBSERVABLES
    !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    var=1.0d0 ! total cross section
    ptz=GETPT(P_Z(1:4))
    yz=rapidity(P_Z(1:4))

    kk=1
    www=wgts(kk)
    l=(kk-1)*num_plots
    CALL mfill(l+1,var,www)
    CALL mfill(l+2,ptz,www)
    CALL mfill(l+3,DABS(yz),www)
    CALL mfill(l+4,yz,www)
    DO j=1,N_JET_ALG
       k=(j-1)*19
       IF(j.EQ.1)THEN
          ! kt with R=0.5
          palg=1d0
          JET_KTRADIUS=0.5d0
          JET_KTPTMIN=0d0
       ELSEIF(j.EQ.2)THEN
          ! anti-kt with R=0.5
          palg=1d0
          JET_KTRADIUS=0.5d0
          JET_KTPTMIN=0d0
       ELSEIF(j.EQ.3)THEN
          ! anti-kt with R=1.0
          palg=-1d0
          JET_KTRADIUS=1.0d0
          JET_KTPTMIN=0d0
       ENDIF
       CALL fastjetppgenkt(PTRACK,NTRACKS,JET_KTRADIUS,JET_KTPTMIN,PALG,PJET,NJET,JETVEC) ! fastjet internal function: the returned
           ! PJET (momentum),NJET(number of jets)
           ! JETVEC(itrack) is the jet index for itrack TRACK
       ! check that jets are ordered in pt
       DO i=1,NJET-1
          IF(GETPT(PJET(1:4,i)).LT.GETPT(PJET(1:4,i+1)))THEN
             WRITE(*,*)"ERROR: Jets are not ordered in plot_MLM_py8_fill"
             STOP
       ENDDO
       IF(NTRACKS.GE.1)THEN
! it is defined in fastjetfortran.cc or fjcorefortran.cc
! return the maximum of the dmin encountered during all recombinations
! up to the one that led to an n-jet final state; identical to
! exclusive_dmerge, except in cases where the dmin do not increase
! monotonically.
          d01=DSQRT(fastjetdmergemax(0))
       ELSE
          d01=-1d8
       ENDIF
       IF(NTRACKS.GE.2)THEN
          d12=DSQRT(fastjetdmergemax(1))
       ELSE
          d12=-1d8
       ENDIF
       IF(NTRACKS.GE.3)THEN
          d23=DSQRT(fastjetdmergemax(2))
       ELSE
          d23=-1d8
       ENDIF
       NJET30=0
       ptj1=-1d0
       yj1=-1d8
       ptj2=-1d0
       yj2=-1d8
       ptj3=-1d0
       yj3=-1d8
       ptj4=-1d0
       yj4=-1d8
       !  njet30 is pT > 30 and |eta|<2.4
       DO I=1,NJET
          ptj=GETPT(PJET(1:4,I))
          etaj=prapidity(PJET(1:4,I))
          yj=rapidity(PJET(1:4,I))
          IF(ptj.GE.30.d0.AND.DABS(etaj).LT.2.4d0)THEN
             njet30=njet30+1
             IF(njet30.EQ.1)THEN
                ptj1=ptj
                yj1=yj
             ELSEIF(njet30.EQ.2)THEN
                ptj2=ptj
                yj2=yj
             ELSEIF(njet30.EQ.3)THEN
                ptj3=ptj
                yj3=yj
             ELSEIF(njet30.EQ.4)THEN
                ptj4=ptj
                yj4=yj
            ENDIF
          ENDIF
       ENDDO
       CALL mfill(l+k+5,DBLE(njet30),www)
       CALL mfill(l+k+6,ptj1,www)
       CALL mfill(l+k+7,ptj2,www)
       CALL mfill(l+k+8,ptj3,www)
       CALL mfill(l+k+9,ptj4,www)
       CALL mfill(l+k+10,DABS(yj1),www)
       CALL mfill(l+k+11,DABS(yj2),www)
       CALL mfill(l+k+12,DABS(yj3),www)
       CALL mfill(l+k+13,DABS(yj4),www)
       CALL mfill(l+k+14,d01,www)
       IF(d01.GT.0d0)CALL mfill(l+k+15,LOG10(d01),www)
       CALL mfill(l+k+16,d12,www)
       IF(d12.GT.0d0)CALL mfill(l+k+17,LOG10(d12),www)
       CALL mfill(l+k+18,d23,www)
       IF(d23.GT.0d0)CALL mfill(l+k+19,LOG10(d23),www)
       CALL mfill(l+k+20,d34,www)
       IF(d34.GT.0d0)CALL mfill(l+k+21,LOG10(d34),www)
       IF(njet30.EQ.1)THEN
          ydiff=DABS(yz-yj1)/2d0
          ysum=DABS(yz+yj1)/2d0
          CALL mfill(l+k+22,ydiff,www)
          CALL mfill(l+k+23,ysum,www)
       ENDIF
    ENDDO

999 RETURN      
  END SUBROUTINE plot_MLM_py8_fill

  SUBROUTINE open_topdrawer_file
    IMPLICIT NONE
    LOGICAL::useitmax
    OPEN(unit=99,FILE='helaconia_MLM_py8.top',&
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
    OPEN(unit=97,FILE='helaconia_MLM_py8.gnu',&
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
    OPEN(unit=96,FILE='helaconia_MLM_py8.C',&
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

END MODULE plot_MLM_py8_user
