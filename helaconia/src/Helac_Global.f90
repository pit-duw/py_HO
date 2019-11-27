MODULE Helac_Global
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! This file includes the global variables in the helac-onia
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
IMPLICIT NONE
SAVE
INTEGER,PARAMETER::SGL=SELECTED_REAL_KIND(p=6)
INTEGER,PARAMETER::DBL=SELECTED_REAL_KIND(p=13)
!REAL(KIND=DBL),PARAMETER::pi=3.1415926535897932385d0
! vertex  common_coup.h
COMPLEX(KIND=DBL),DIMENSION(31:35,31:35,31:35),PUBLIC::zgv3
COMPLEX(KIND=DBL),DIMENSION(31:36,31:36,31:36),PUBLIC::zgvvx
COMPLEX(KIND=DBL),DIMENSION(31:35,31:35,31:35,31:35),PUBLIC::zgv4
COMPLEX(KIND=DBL),DIMENSION(31:35,-12:-1,1:12),PUBLIC::zgvffl
COMPLEX(KIND=DBL),DIMENSION(31:35,-12:-1,1:12),PUBLIC::zgvffr
COMPLEX(KIND=DBL),DIMENSION(31:35,31:35,41:44),PUBLIC::zgvvs
COMPLEX(KIND=DBL),DIMENSION(41:44,31:35,31:35),PUBLIC::zgsvv
COMPLEX(KIND=DBL),DIMENSION(31:35,41:44,41:44),PUBLIC::zgvss
COMPLEX(KIND=DBL),DIMENSION(41:44,41:44,31:35),PUBLIC::zgssv
COMPLEX(KIND=DBL),DIMENSION(31:35,31:35,41:44,41:44),PUBLIC::zgvvss
COMPLEX(KIND=DBL),DIMENSION(41:44,41:44,31:35,31:35),PUBLIC::zgssvv
COMPLEX(KIND=DBL),DIMENSION(41:44,-12:-1,1:12),PUBLIC::zgsffl
COMPLEX(KIND=DBL),DIMENSION(41:44,-12:-1,1:12),PUBLIC::zgsffr
COMPLEX(KIND=DBL),DIMENSION(41:45,41:45,41:45),PUBLIC::zgs3
COMPLEX(KIND=DBL),DIMENSION(41:44,41:44,41:44,41:44),PUBLIC::zgs4
! masses & widths   common_masses.h
REAL(KIND=DBL),DIMENSION(-12:44),PUBLIC::parmas,parwid
! flags: iflag=0 sum over all helicity configurations
!              1 specific helicity configuration
!        iunitary=0  Feynman gauge
!                 1  unitary gauge
!        ihiggs=0    not inclusion the higgs particle as an intermediate state
!               1    inclusion the higgs particle as an intermediate state
!        iwidth=0    fixed scheme for the introduction of the width of W and Z
!               1    complex scheme for the introduction of the width of W and Z
! common_flags.h
INTEGER,PUBLIC::iflag,iunitary,ihiggs,iwidth
! common_print.h
INTEGER,PUBLIC::nunit1=6,nunit2,nunit3
! common_strf.h
REAL(KIND=DBL),PUBLIC::wjac,xp1,xp2,scale,ehat
INTEGER,PUBLIC::istruc,iPDFGUP1,iPDFGUP2,iPDFSUP1,iPDFSUP2  ! istruc=1 need pdfs
                                                         ! istruc=0 don't need pdfs, such as e+e-
REAL(KIND=DBL),DIMENSION(2),PUBLIC::ebeam           ! beam energy, so that the lab frame don't have to in collision c.m.
LOGICAL,PUBLIC::absrap=.TRUE.,labeqcoll=.TRUE.
!common_lha.h
REAL(KIND=DBL),PUBLIC::aqedup,aqcdup
!common_qcdrun.h
INTEGER::irun
LOGICAL::onlyqcd,withqcd
LOGICAL::useMCFMrun=.TRUE. ! use the MCFM alpha_s running (new)
LOGICAL::uselhapdf=.FALSE. ! use the LHAPDF (new)
LOGICAL::lhapdfwrapinit=.FALSE. ! a tag for using lhapdfwarp, which forbid it before calling LHAPDF but after alphasPDF
CHARACTER(len=150)::LHAPath="" ! the path of the LHAPDF
LOGICAL::reweight_pdf=.FALSE. ! use reweighting to estimate pdf uncertainty
INTEGER::ipdf_set_min=21101   ! the first of the error pdf sets
INTEGER::ipdf_set_max=21140   ! the last of the error pdf sets
LOGICAL::reweight_scale=.FALSE. ! use reweighting to estimate scale uncertainty
REAL(KIND=DBL)::rw_Rscale_down=0.5d0   ! lower bound for ren scale variations
REAL(KIND=DBL)::rw_Rscale_up=2.0d0     ! upper bound for ren scale variations
REAL(KIND=DBL)::rw_Fscale_down=0.5d0   ! lower bound for fac scale variations
REAL(KIND=DBL)::rw_Fscale_up=2.0d0     ! upper bound for fac scale variations 
LOGICAL::reweight_scale_phase=.FALSE.
INTEGER::alphas_power=0,ho_nscale=0,ho_npdf=0
REAL(KIND=DBL),DIMENSION(3,3)::wgtxsecmu
INTEGER::reweight_Fscale_phase=0
REAL(KIND=DBL),DIMENSION(:),ALLOCATABLE::wgtxsecpdf
INTEGER,DIMENSION(:),ALLOCATABLE::idpdf
! constants.h
REAL(KIND=DBL),PUBLIC::gfermi=1.16639d-5 ! Fermi coupling constant in GeV^(-2)
REAL(KIND=DBL),PUBLIC::rzm=91.188d0      ! the mass of Z in GeV
REAL(KIND=DBL),PUBLIC::rwm=80.419d0      ! the mass of W in GeV
REAL(KIND=DBL),PUBLIC::rzw=2.446d0       ! the width of Z in GeV
REAL(KIND=DBL),PUBLIC::rww=2.048d0       ! the width of W in GeV
REAL(KIND=DBL),PUBLIC::alphaQCD2=0.118d0  ! alphas(mZ)=0.118
REAL(KIND=DBL),PUBLIC::rhm=130.0d0       ! the mass of Higgs in GeV
REAL(KIND=DBL),PUBLIC::rhw=4.291d-3      ! the width of Higgs in GeV
REAL(KIND=DBL),PUBLIC::rcm=1.5d0           ! the mass of Charm quark in GeV
REAL(KIND=DBL),PUBLIC::rbm=4.5d0           ! the mass of Beauty quark in GeV
REAL(KIND=DBL),PUBLIC::rtm=174.3d0       ! the mass of Top quark in GeV
REAL(KIND=DBL),PUBLIC::rtw=1.6d0         ! the width of Top quark in GeV
REAL(KIND=DBL),PUBLIC::alphaem=1d0/137d0     ! alpha_em
INTEGER,PUBLIC::nloop=1           ! running of alphaQCD with order (new added one)
CHARACTER(len=7),PUBLIC::pdlabel="cteq6l1" ! pdf label (new added one) 
! common_int.h
INTEGER,PUBLIC::n,nhad                            ! number of parton/hadron particles involved
INTEGER,DIMENSION(20),PUBLIC::io,ioh,&                ! =1 incoming =-1 outgoing
                              ifl,iflh            ! flavor for partons and hadrons
! common_feyn.h
!REAL(KIND=DBL),DIMENSION(20,5),PUBLIC::pmom
INTEGER,PUBLIC,PARAMETER::ISMAX=20000
INTEGER,PUBLIC::ng=0
INTEGER,DIMENSION(ISMAX,11,8),PUBLIC::is
! helicity and colors common_helc.h
INTEGER,DIMENSION(20),PUBLIC::ipol                !ipol(k)=2:rigth-handed 
                                                  !ipol(k)=1:left-handed 
                                                  !ipol(k)=3:long
INTEGER,DIMENSION(20,2),PUBLIC::icol
! common_new.h
REAL(KIND=DBL),PUBLIC::avhel,&                    ! average over helicity
                       avcol,&                    ! average over color
					   symet                      ! symmetry factors         
INTEGER,PUBLIC::ncc                               ! the number of color configurations
INTEGER,PUBLIC::nhc                               ! the number of helicity configurations
INTEGER,PUBLIC::nphyhc                            ! the number of physical helicity configurations
INTEGER,PUBLIC::imc,iranhel                       ! iranhel=0,we use the specific helicity
                                                  ! iranhel=1, use random helicity in Helac_pan1
												  ! iranhel=2, use random helicity for eps(L_z) too
												  ! iranhel=3, use random helicity for Lz and Sz
												  ! similar for imc in Helac_master_f (color selection)
REAL(KIND=DBL),PUBLIC::repeat                     ! repeat=0 run both phase of Helac_Phegas
                                                  ! repeat=1 run only the first phase
												  ! repeat=2 run only the second phase
! common_unweight.h
LOGICAL,PUBLIC::lwri   
! common_unw.h         
INTEGER,DIMENSION(20,2),PUBLIC::icol_un
! common_onep.h
LOGICAL,PUBLIC::onep                         ! onep=.TRUE. the user supply the phase space point
CHARACTER(100)::momfile='Mom'      ! if onep=.TRUE. the name of the input file of the momenta
CHARACTER(100)::momout='Momout'               ! name of output file if onep=.TRUE.
! parameter_NCOL.h
INTEGER,PARAMETER,PUBLIC::NCOL=3                  ! SUNN=3
! nguess.h
INTEGER,PARAMETER,PUBLIC::ngues=1500
! common_mom.h
!zq(1:n,1)=(p0+pz,pz);zq(1:n,2)=(p0-pz,pt);zq(1:n,3)=(px,py);
!zq(1:n,4)=(px,-py);zq(1:n,5)=(m,p)
! zQq is the momenta for hadron array
COMPLEX(KIND=DBL),DIMENSION(20,5),PUBLIC::zq   !,zQq
!common_phegas_mom.h
REAL(KIND=DBL),DIMENSION(20,5),PUBLIC::phegas_pmom
REAL(KIND=DBL),DIMENSION(20,5),PUBLIC::hadron_pmom
! common_norm.h
REAL(KIND=DBL),PUBLIC::wrest
! common_warn.h
INTEGER,DIMENSION(25),PUBLIC::iwarning,iwonders
! common_psp.h
INTEGER,PUBLIC::i_psp
! common_debug.h
INTEGER,PUBLIC::idebug=0  ! idebug=1 input some information on the screen
                        ! in Subroutine Helac_optimize in Main_Program.f90
! common_cuts.h various cut parameters
REAL(KIND=DBL),DIMENSION(3:20),PUBLIC::ec,c1,c2
REAL(KIND=DBL),DIMENSION(3:20,3:20),PUBLIC::cc,gmas
REAL(KIND=DBL),DIMENSION(3:20),PUBLIC::ptc,etac,yycut,yycutlow,xFcut,xFcutlow,maxptc
LOGICAL::xFcutflag=.FALSE.
REAL(KIND=DBL),DIMENSION(3:20,3:20),PUBLIC::drc
REAL(KIND=DBL),PUBLIC::ptmissc
! special for NLO* and NNLO* in pp(bar) collisions
REAL(KIND=DBL),DIMENSION(1:2,3:20),PUBLIC::gbeammass
! Input filename
CHARACTER(len=20),PARAMETER::Input_File="user.inp"
CHARACTER(len=24),PARAMETER::input_dir="./input/"
CHARACTER(len=24),PARAMETER::tmp_dir="./tmp/"
CHARACTER(len=24),PARAMETER::output_dir="./output/"
CHARACTER(len=24),PARAMETER::decay_dir="./decay/"
! some keywords
LOGICAL::ktrw=.FALSE.
INTEGER::COLL_TYPE=1  ! colpar   1=pp  2=ppbar  3=epem
INTEGER::CONSTNUM=0   ! 0 use the default constants value, 1 use the user defined
INTEGER::MaxAdapt=15   ! the value for maxadap_adapt in Adapt.f90
CHARACTER(len=20)::histo_file="hi_file",error_file='err_file'
LOGICAL::SUMMATIONQ=.FALSE.
INTEGER::GLOBALINIT_Phegas=0,GLOBALINIT_Bookin=0,GLOBALINIT_Sethel=0,&
         GLOBALINIT_master=0,GLOBALINIT_optimize=0,GLOBALINIT_feyn=0,GLOBALINIT_average=0,&
		 GLOBALINIT_checkgluons=0,GLOBALINIT_id=0,GLOBALINIT_id2=0,GLOBALINIT_id3=0,&
		 GLOBALINIT_unwei=0,GLOBALINIT_redo=0,GLOBALINIT_setmin=0,GLOBALINIT_setmax=0,&
		 GLOBALINIT_gen_x=0,GLOBALINIT_adapt_gen=0,GLOBALINIT_physicalpol=0,GLOBALINIT_mader=0,&
		 GLOBALINIT_ls2j=0,GLOBALINIT_physicalpol2=0,GLOBALINIT_idoctet=0
! appendix A of Phys.Rev.D  57 (1998)  4258
! 1: Recoil(Helicity) frame  2:Collins-Soper frame  3: Gottfried-Jackson frame  4: Target frame
INTEGER::PolarFrame=1       ! Polarization frames
INTEGER,DIMENSION(10,3)::Quarkonium
INTEGER,DIMENSION(20)::Quarkonium2,Quarkonium3   ! hadron and parton
REAL(KIND=DBL),DIMENSION(20)::Qhelran          ! the random number for Quarkonium when iranhel=2
REAL(KIND=DBL),DIMENSION(20)::Decayran
                                                 ! i.e.r(eps(L_z))
REAL(KIND=DBL),DIMENSION(20)::QSzran           ! the random number for Quarkonium when iranhel=3
                                                ! i.e. r(eps(S_z))
! the Pt for the first hadron
REAL(KIND=DBL)::PTFirst=3d0
LOGICAL::dSigmadPTQ=.FALSE.
INTEGER::Pwavenum=0                               ! the number of P wave
INTEGER,DIMENSION(10)::Pwave                     ! returns the number of hadron in p wave array
INTEGER,DIMENSION(20)::PwaveQ                    ! returns the number of p wave in hadron array
! used in helicity selection of Helac_iniqq in Helac_pan1.f90
!INTEGER::rannum
! Color Singlet Long Distance Matrix Elements
! Long Distance Matrix Element <O(3S1[1])>=|R(0)|^2/4/Pi
! in JHEP 02 (2008) 102  <O(3S1[1])>=(2J+1)*2Nc*|R(0)|^2/4/Pi,
! i.e. LDME****1=<O(2S+1)LJ[1]>/2Nc/(2J+1)
! For p-wave <O(3P0[1])>=<O(3P1[1])>=<O(3P2[1])>=3*|R'(0)|^2/4/Pi
! Charmonium System
REAL(KIND=DBL),PUBLIC::LDME441001=0.0644578d0,LDME443011=0.0644578d0,&
LDME441111=0.0179049d0,LDME443101=0.0179049d0,LDME443111=0.0179049d0,&
LDME443121=0.0179049d0
! Bottonium System
REAL(KIND=DBL),PUBLIC::LDME553011=0.515556d0,LDME551001=0.515556d0,&
LDME551111=0.515556d0,LDME553101=0.515556d0,LDME553111=0.515556d0,&
LDME553121=0.515556d0
! Bc System
REAL(KIND=DBL),PUBLIC::LDME453011=0.122667d0,LDME451001=0.122667d0,&
LDME451111=0.0478333d0,LDME453101=0.0478333d0,LDME453111=0.0478333d0,&
LDME453121=0.0478333d0
! Color Octet(3*3bar) Long Distance Matrix Elements
! in JHEP 02 (2008) 102 LDME****8=<O((2S+1)LJ[8])>/(Nc^2-1)/(2J+1)
! Charmonium System
REAL(KIND=DBL),PUBLIC::LDME441008=0.0095d0,LDME443018=0.00265d0,&
LDME441118=0.134287d-2,LDME443108=0.134287d-2,LDME443118=0.134287d-2,&
LDME443128=0.134287d-2
! Bottonium System
REAL(KIND=DBL),PUBLIC::LDME551008=0.386667d-2,LDME553018=0.386667d-2,&
LDME551118=0.386667d-2,LDME553108=0.386667d-2,LDME553118=0.386667d-2,&
LDME553128=0.386667d-2
! Bc System
REAL(KIND=DBL),PUBLIC::LDME451008=0.920004d-3,LDME453018=0.920004d-3,&
LDME451118=0.35875d-3,LDME453108=0.35875d-3,LDME453118=0.35875d-3,&
LDME453128=0.35875d-3
REAL(KIND=DBL),PUBLIC::LDMEwt
! the rapidity (not pesudo) cut for the first final particle
REAL(KIND=DBL),PUBLIC::y1cup,y1clow
! the rapidity (not pesudo) cut for the second final particle
!REAL(KIND=DBL),PUBLIC::y2cup,y2clow
! The scale scheme
INTEGER,PUBLIC::nscheme=0
! The scale value in fixed scheme (nscheme=0)
REAL(KIND=DBL),PUBLIC::fschemevalue=12.8d0,scalefactor=1d0
LOGICAL,PUBLIC::exp3pjQ
INTEGER,PUBLIC::imode,iLSJ,ihel1,ihel2
INTEGER,PUBLIC::Num3pj
INTEGER,DIMENSION(10),PUBLIC::Pwave3pj     ! returns the number of hadron in 3pj wave array
INTEGER,DIMENSION(20),PUBLIC::Pwave3pjQ    ! returns the number of 3pj wave in hadron array
INTEGER,DIMENSION(10),PUBLIC::JNum3pj      ! returns J number in 3pj wave array
INTEGER,PUBLIC::PhyPolNum=0,SDPart=0,iqnum
LOGICAL,PUBLIC::FirstHelicityQ=.TRUE.,qsumQ=.FALSE.
INTEGER,PUBLIC::gener
INTEGER,DIMENSION(20),PUBLIC::hadron2parton  ! return the first parton index in hadron array
INTEGER,DIMENSION(20),PUBLIC::parton2hadron  ! return the  hadron index in parton array
INTEGER,DIMENSION(20),PUBLIC::parton2hadrontype ! return whether the parton i is charmonium (1), bottonium (2), Bc(3), no (0)
INTEGER,PUBLIC::octetmsinglet=0 ! the integer for choosing whether the nonet component or singlet component in octet
LOGICAL,PUBLIC::octetQ=.FALSE.
INTEGER,PUBLIC::octetnum
INTEGER,DIMENSION(10),PUBLIC::octetlist
! kt-measure
INTEGER,PUBLIC::ktmeasure=1 ! 1: ordinary ktmeasure 2: ktmeasure in SHERPA or MG
REAL(KIND=DBL),PUBLIC::RR=1d0 ! parameter R or D
! a flag for whether only measure speed or not
LOGICAL,PUBLIC::measurespeed
! number of lhe events
INTEGER,PUBLIC::Nevents=0
INTEGER,PUBLIC::ntotps=0
! str for the process
CHARACTER(len=400)::process
! plot flags
LOGICAL::plot_output=.FALSE.
LOGICAL::topdrawer_output=.FALSE.
LOGICAL::gnuplot_output=.FALSE.
LOGICAL::root_output=.FALSE.
! shower variables
! emep_ISR_shower=0 or 1
! 0 is no shower
! 1 is QEDPS
INTEGER::emep_ISR_shower=0
REAL(KIND(1d0))::Q2OUT_QEDPS,Q2MAX_QEDPS
! Decay informations
INTEGER,PARAMETER::MAX_DecayChain=99
INTEGER,DIMENSION(MAX_DecayChain,-1:5)::DecayChains ! decay chains
INTEGER::NDecayChains ! number of decay chains in the process
INTEGER::AllNDecayChains ! total number of decay chains
REAL(KIND(1d0)),DIMENSION(MAX_DecayChain)::DecayBR ! branching ratios
INTEGER,DIMENSION(20,0:MAX_DecayChain)::iflh2DecayChains ! iflh -> DecayChains
INTEGER,DIMENSION(100,0:MAX_DecayChain)::JMO2DecayChains  ! for the second and later decays
INTEGER::NDecayIflh
LOGICAL::MCoHelicity=.FALSE.
REAL(KIND(1d0))::RMCoH,weight_br=1d0
! special cuts
INTEGER::literature_cutoffs
! fixed-target or not
LOGICAL::fixtarget=.FALSE.,fixtargetrev=.FALSE.
REAL(KIND(1d0))::FT_M,FT_E1
END MODULE Helac_Global
