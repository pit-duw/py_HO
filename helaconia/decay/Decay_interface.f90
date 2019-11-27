MODULE Decay_interface
  USE Helac_Global
  USE HOVll
  USE Helac_Func_1
  USE Kinetic_Func
  USE QEDPS_interface
  USE Constants
  USE Helac_ranmar_mod
  USE HO_t2bw
  USE HO_chi2psia
  IMPLICIT NONE
  ! cut on the decayed leptons
  REAL(KIND(1d0))::ho_dptcl,ho_dycllow,ho_dycl,ho_detacl,ho_decl,ho_dcl
  ! decay params
  REAL(KIND(1d0))::ho_dmass_jpsi,ho_dmass_chic0,ho_dmass_chic1,ho_dmass_chic2,&
       ho_dmass_hc,ho_dmass_etac,ho_dmass_psi2S,ho_dmass_Y1S,ho_dmass_chib11P,&
       ho_dmass_chib01P,ho_dmass_chib21P,ho_dmass_hb1P,ho_dmass_chib02P,&
       ho_dmass_chib12P,ho_dmass_chib22P,ho_dmass_Y2S,ho_dmass_Y3S,&
       ho_dmass_chib3P
  INTEGER::ho_colmax=0
  SAVE
CONTAINS
  SUBROUTINE HO_Decay(wgtbr)
    IMPLICIT NONE
    INCLUDE "stdhep.inc"
    INTEGER::i,j,k,kk,m,ndecay
    REAL(KIND(1d0)),INTENT(OUT)::wgtbr
    REAL(KIND(1d0))::br,r
    REAL(KIND(1d0)),DIMENSION(MAX_DecayChain)::braccum
    ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL
    INTEGER,DIMENSION(2,1000)::children
    INTEGER,DIMENSION(1000)::imoth,mother
    INTEGER::i1max,i2max
    ho_colmax=0
    wgtbr=1d0
    ndecay=0
    IF(emep_ISR_shower.EQ.1)THEN
       NEVHEP=1
       NHEP=NPTCL-1
       imoth(1:1000)=0
       mother(1:1000)=0
       children(1:2,1:1000)=0
       imoth(1)=1
       imoth(2)=2
       i1max=1
       i2max=2
       DO i=1,NVPH_QEDPS-1 ! only from parton shower
          IF(NLPTN(6,i).EQ.-1)THEN
             IF(NLPTN(3,i).EQ.2.AND.NLPTN(4,i).NE.0)THEN
                children(1,i)=NLPTN(4,i)+i
                children(2,i)=children(1,i)+1
             ENDIF
             IF(NLPTN(2,i).NE.0)THEN
                mother(i)=i-NLPTN(2,i)
                mother(i+1)=mother(i)
                imoth(i)=imoth(mother(i))
                imoth(i+1)=imoth(mother(i))
                IF(imoth(i).EQ.1.AND.i.GT.i1max)i1max=i
                IF(imoth(i).EQ.2.AND.i.GT.i2max)i2max=i
             ENDIF
          ENDIF
       ENDDO
       DO i=1,NPTCL
          IF(i.EQ.NVPH_QEDPS)CYCLE
          IF(i.LT.NVPH_QEDPS)THEN
             ISTHEP(i)=0
             IDHEP(i)=NLPTN(1,i)
             ISHEP(i)=NLPTN(6,i)
             IF(imoth(i).EQ.1)THEN
                JMOHEP(1,i)=mother(i)
                JMOHEP(2,i)=0
             ELSEIF(imoth(i).EQ.2)THEN
                JMOHEP(1,i)=0
                JMOHEP(2,i)=mother(i)
             ELSE
                JMOHEP(1,i)=0
                JMOHEP(2,i)=0
             ENDIF
             JDAHEP(1,i)=children(1,i)
             JDAHEP(2,i)=children(2,i)
             PHEP(1:4,i)=PLPTN(4:7,i)
             IF(ABS(NLPTN(1,i)).EQ.11)THEN
                PHEP(5,i)=amel
             ELSE
                PHEP(5,i)=0d0
             ENDIF
             VHEP(1:4,i)=0d0
             IPOLHEP(i)=99
             ICOLHEP(1:2,i)=0
          ELSE
             ISTHEP(i-1)=0
             IDHEP(i-1)=NLPTN(1,i)
             ISHEP(i-1)=NLPTN(6,i)
             JMOHEP(1,i-1)=i1max
             JMOHEP(2,i-1)=i2max
             JDAHEP(1:2,i-1)=0
             PHEP(1:4,i-1)=PLPTN(4:7,i)
             PHEP(5,i-1)=parmas(iflh(i-NVPH_QEDPS+2))
             VHEP(1:4,i-1)=0d0
             IPOLHEP(i-1)=99
             IF(ho_colmax.LT.icol_un(i-NVPH_QEDPS+2,2))ho_colmax=icol_un(i-NVPH_QEDPS+2,2)
             ICOLHEP(1,i-1)=icol_un(i-NVPH_QEDPS+2,2)+100
             IF(ICOLHEP(1,i-1).EQ.100)ICOLHEP(1,i-1)=0
             IF(ho_colmax.LT.icol_un(i-NVPH_QEDPS+2,1))ho_colmax=icol_un(i-NVPH_QEDPS+2,1)
             ICOLHEP(2,i-1)=icol_un(i-NVPH_QEDPS+2,1)+100
             IF(ICOLHEP(2,i-1).EQ.100)ICOLHEP(2,i-1)=0
          ENDIF
       ENDDO
       JDAHEP(1,i1max)=NVPH_QEDPS
       JDAHEP(2,i1max)=NPTCL-1
       JDAHEP(1,i2max)=NVPH_QEDPS
       JDAHEP(2,i2max)=NPTCL-1
    ELSE
       CALL partonmom2hadronmom
       NEVHEP=1
       NHEP=nhad
       DO i=1,nhad
          ISTHEP(i)=0
          IDHEP(i)=pdgt(iflh(i))
          IF(i.GT.2.AND.i.LE.nhad)THEN
             ISHEP(i)=1
             JMOHEP(1,i)=1
             JMOHEP(2,i)=2
             JDAHEP(1,i)=0
             JDAHEP(2,i)=0
             IF(ho_colmax.LT.icol_un(i,2))ho_colmax=icol_un(i,2)
             ICOLHEP(1,i)=icol_un(i,2)+100
             IF(ICOLHEP(1,i).EQ.100)ICOLHEP(1,i)=0
             IF(ho_colmax.LT.icol_un(i,1))ho_colmax=icol_un(i,1)
             ICOLHEP(2,i)=icol_un(i,1)+100
             IF(ICOLHEP(2,i).EQ.100)ICOLHEP(2,i)=0
          ELSE
             ISHEP(i)=-1
             JMOHEP(1,i)=0
             JMOHEP(2,i)=0
             JDAHEP(1,i)=3
             JDAHEP(2,i)=nhad
             IF(ho_colmax.LT.icol_un(i,1))ho_colmax=icol_un(i,1)
             ICOLHEP(1,i)=icol_un(i,1)+100
             IF(ICOLHEP(1,i).EQ.100)ICOLHEP(1,i)=0
             IF(ho_colmax.LT.icol_un(i,2))ho_colmax=icol_un(i,2)
             ICOLHEP(2,i)=icol_un(i,2)+100
             IF(ICOLHEP(2,i).EQ.100)ICOLHEP(2,i)=0
          ENDIF
          PHEP(1:5,i)=hadron_pmom(i,1:5)
          VHEP(1:4,i)=0d0
          IPOLHEP(i)=99
       ENDDO
    ENDIF
    ! first decay
    DO i=1,nhad
       j=iflh2DecayChains(i,0)
       IF(j.GT.0)THEN
          br=0d0
          braccum(1:MAX_DecayChain)=0d0
          DO k=1,j
             m=iflh2DecayChains(i,k)
             br=br+DecayBR(m)
             braccum(k)=br
          ENDDO
          ! randomly select the decay process
          IF(j.EQ.1)THEN
             kk=1
          ELSE
             ndecay=ndecay+1
             r=Decayran(ndecay)*br
             DO k=1,j
                IF(braccum(k).GT.r)THEN
                   kk=k
                   EXIT
                ELSEIF(k.EQ.j)THEN
                   kk=j
                   EXIT
                ENDIF
             ENDDO
          ENDIF
          m=iflh2DecayChains(i,kk)
          CALL HO_One_Decay(i,m)
          wgtbr=wgtbr*br
       ENDIF
    ENDDO
    ! the next decays (after the first one)
    CALL HO_Decay_next(wgtbr)
    RETURN
  END SUBROUTINE HO_Decay

  RECURSIVE SUBROUTINE HO_Decay_next(wgtbr)
    REAL(KIND(1d0)),INTENT(INOUT)::wgtbr
    INTEGER::i,j,k,m,lll,kk,ndecay
    REAL(KIND(1d0))::br,r
    LOGICAL::decayq
    REAL(KIND(1d0)),DIMENSION(MAX_DecayChain)::braccum
    INCLUDE "stdhep.inc"
    ! first to assign the decay chains
    decayq=.FALSE.
    DO i=1,NHEP
       JMO2DecayChains(i,0)=0
       IF(ISHEP(i).EQ.-1)CYCLE
       IF(ISTHEP(i).EQ.1)CYCLE
       lll=0
       DO j=1,AllNDecayChains
          IF(IDHEP(i).EQ.pdgt(DecayChains(j,0)))THEN
             JMO2DecayChains(i,0)=JMO2DecayChains(i,0)+1
             JMO2DecayChains(i,lll+1)=j
             lll=lll+1
             decayq=.TRUE.
          ENDIF
       ENDDO
    ENDDO
    IF(.NOT.decayq)RETURN
    ! the next decays
    DO i=1,NHEP
       j=JMO2DecayChains(i,0)
       IF(j.GT.0)THEN
          br=0d0
          braccum(1:MAX_DecayChain)=0d0
          DO k=1,j
             m=JMO2DecayChains(i,k)
             br=br+DecayBR(m)
             braccum(k)=br
          ENDDO
          ! randomly select the decay process
          IF(j.EQ.1)THEN
             kk=1
          ELSE
             ndecay=ndecay+1
             r=Helac_rnmy(0)*br
             DO k=1,j
                IF(braccum(k).GT.r)THEN
                   kk=k
                   EXIT
                ELSEIF(k.EQ.j)THEN
                   kk=j
                   EXIT
                ENDIF
             ENDDO
          ENDIF
          m=JMO2DecayChains(i,kk)
          CALL HO_One_Decay_next(i,m)
          wgtbr=wgtbr*br
       ENDIF
    ENDDO
    CALL HO_Decay_next(wgtbr)
    RETURN
  END SUBROUTINE HO_Decay_next

  SUBROUTINE HO_One_Decay(imoth,idecay)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::imoth,idecay
    INTEGER,PARAMETER::available_ndecay=5
    CHARACTER(len=20),DIMENSION(available_ndecay)::available_decay
    INTEGER::i,isign
    INTEGER::lavail
    INTEGER,DIMENSION(3)::index2
    REAL(KIND(1d0)),DIMENSION(0:3)::PM,PD1,PD2
    REAL(KIND(1d0)),DIMENSION(3)::svec
    REAL(KIND(1d0))::YL,YR,mchi,mpsi,lambdath,sigma1100,rr
    INTEGER::iheld1,iheld2,anti
    INCLUDE "stdhep.inc"
    ! QEDPS parameters
    REAL(KIND(1d0)),DIMENSION(10,1000)::PLPTN
    INTEGER,DIMENSION(10,1000)::NLPTN
    INTEGER::NPTCL
    COMMON/QPLIST/PLPTN,NLPTN,NPTCL

    CALL HO_judge_avail_decay(imoth,idecay,lavail)
    IF(lavail.GT.0)THEN
       SELECT CASE(lavail)
       CASE(1)
          ! 3S11 > l+ l-
          CALL Qeveryhel(index2,ipol(imoth),imoth)
          IF(emep_ISR_shower.EQ.1)THEN
             ! isr in e-e+ collisions
             PM(0)=PLPTN(7,imoth-2+NVPH_QEDPS)
             PM(1:3)=PLPTN(4:6,imoth-2+NVPH_QEDPS)
          ELSE
             PM(0)=hadron_pmom(imoth,4)
             PM(1:3)=hadron_pmom(imoth,1:3)
          END IF
          svec(1:3)=PM(1:3)
          IF(emep_ISR_shower.EQ.1)THEN
             JMOHEP(1,NHEP+1:NHEP+2)=imoth+NVPH_QEDPS-3
          ELSE
             JMOHEP(1,NHEP+1:NHEP+2)=imoth
          ENDIF
          JMOHEP(2,NHEP+1:NHEP+2)=0
          ! for testing
          !lambdath=-0.5d0
          !sigma1100=(1d0+lambdath)/(1d0-lambdath)
          !sigma1100=1d0/(1d0+2d0*sigma1100)
          !rr=Helac_rnmy(0)
          !IF(rr.LE.sigma1100)THEN
          !   ! longitudial polarized in HX
          !   CALL HO_Vll00(PM,svec,PD1,PD2)
          !ELSE
          !   ! transverse polarized in HX
          !   CALL HO_Vll11(PM,svec,PD1,PD2)
          !ENDIF
          ! end for testing
          IF(index2(3).LE.2)THEN
             ! transverse polarized in HX
             CALL HO_Vll11(PM,svec,PD1,PD2)
             IF(index2(3).EQ.1)THEN
                IPOLHEP(JMOHEP(1,NHEP+1))=1
             ELSE
                IPOLHEP(JMOHEP(1,NHEP+1))=-1
             ENDIF
          ELSE
             ! longitudial polarized in HX
             CALL HO_Vll00(PM,svec,PD1,PD2)
             IPOLHEP(JMOHEP(1,NHEP+1))=0
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          PHEP(5,NHEP+1:NHEP+2)=0d0
          PHEP(1:3,NHEP+1)=PD1(1:3)
          PHEP(1:3,NHEP+2)=PD2(1:3)
          PHEP(4,NHEP+1)=PD1(0)
          PHEP(4,NHEP+2)=PD2(0)
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          IPOLHEP(NHEP+1)=99 ! decay products will not decay anymore
          IPOLHEP(NHEP+2)=99
          ICOLHEP(1:2,NHEP+1:NHEP+2)=0
          NHEP=NHEP+2
       CASE(2)
          ! W+- > f1 f2
          IF(emep_ISR_shower.EQ.1)THEN
             ! isr in e-e+ collisions
             PM(0)=PLPTN(7,imoth-2+NVPH_QEDPS)
             PM(1:3)=PLPTN(4:6,imoth-2+NVPH_QEDPS)
          ELSE
             PM(0)=hadron_pmom(imoth,4)
             PM(1:3)=hadron_pmom(imoth,1:3)
          END IF
          svec(1:3)=PM(1:3)
          IF(emep_ISR_shower.EQ.1)THEN
             JMOHEP(1,NHEP+1:NHEP+2)=imoth+NVPH_QEDPS-3
          ELSE
             JMOHEP(1,NHEP+1:NHEP+2)=imoth
          ENDIF
          JMOHEP(2,NHEP+1:NHEP+2)=0
          IF(ipol(imoth).LE.2)THEN
             ! transverse
             isign=2*ipol(imoth)-3
             IF(DecayChains(idecay,1).GT.0)THEN
                isign=-isign
             ENDIF
             CALL HO_VllW11(PM,svec,isign,PD1,PD2)
             IF(ipol(imoth).EQ.1)THEN
                IPOLHEP(JMOHEP(1,NHEP+1))=1
             ELSE
                IPOLHEP(JMOHEP(1,NHEP+1))=-1
             ENDIF
          ELSE
             ! longitudial
             CALL HO_Vll00(PM,svec,PD1,PD2)
             IPOLHEP(JMOHEP(1,NHEP+1))=0
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JMOHEP(2,NHEP+1:NHEP+2)=0
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          PHEP(5,NHEP+1:NHEP+2)=0d0
          PHEP(1:3,NHEP+1)=PD1(1:3)
          PHEP(1:3,NHEP+2)=PD2(1:3)
          PHEP(4,NHEP+1)=PD1(0)
          PHEP(4,NHEP+2)=PD2(0)
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          IPOLHEP(NHEP+1)=99 ! decay product will not decay anymore
          IPOLHEP(NHEP+2)=99
          IF(MOD(ABS(DecayChains(idecay,1)),4).EQ.0.OR.&
               MOD(ABS(DecayChains(idecay,1)),4).EQ.3)THEN
             IF(DecayChains(idecay,1).GT.0)THEN
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+1)=0
                ICOLHEP(1,NHEP+2)=ho_colmax+101
                ICOLHEP(2,NHEP+2)=0
                ho_colmax=ho_colmax+1
             ELSE
                ICOLHEP(2,NHEP+1)=0
                ICOLHEP(1,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+2)=0
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ho_colmax=ho_colmax+1
             ENDIF
          ELSE
             ICOLHEP(1:2,NHEP+1:NHEP+2)=0
          ENDIF
          NHEP=NHEP+2
       CASE(3)
          ! Z > f f~
          IF(emep_ISR_shower.EQ.1)THEN
             ! isr in e-e+ collisions
             PM(0)=PLPTN(7,imoth-2+NVPH_QEDPS)
             PM(1:3)=PLPTN(4:6,imoth-2+NVPH_QEDPS)
          ELSE
             PM(0)=hadron_pmom(imoth,4)
             PM(1:3)=hadron_pmom(imoth,1:3)
          END IF
          svec(1:3)=PM(1:3)
          IF(emep_ISR_shower.EQ.1)THEN
             JMOHEP(1,NHEP+1:NHEP+2)=imoth+NVPH_QEDPS-3
          ELSE
             JMOHEP(1,NHEP+1:NHEP+2)=imoth
          ENDIF
          JMOHEP(2,NHEP+1:NHEP+2)=0
          IF(ipol(imoth).LE.2)THEN
             ! transverse
             isign=2*ipol(imoth)-3
             IF(DecayChains(idecay,1).GT.0)THEN
                isign=-isign
             ENDIF
             CALL Coup_Zff(DecayChains(idecay,1),YL,YR)
             CALL HO_Zff11(PM,svec,isign,YL,YR,PD1,PD2)
             IF(ipol(imoth).EQ.1)THEN
                IPOLHEP(JMOHEP(1,NHEP+1))=1
             ELSE
                IPOLHEP(JMOHEP(1,NHEP+1))=-1
             ENDIF
          ELSE
             ! longitudial
             CALL HO_Vll00(PM,svec,PD1,PD2)
             IPOLHEP(JMOHEP(1,NHEP+1))=0
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          PHEP(5,NHEP+1:NHEP+2)=0d0
          PHEP(1:3,NHEP+1)=PD1(1:3)
          PHEP(1:3,NHEP+2)=PD2(1:3)
          PHEP(4,NHEP+1)=PD1(0)
          PHEP(4,NHEP+2)=PD2(0)
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          IPOLHEP(NHEP+1)=99 ! decay product will not decay anymore
          IPOLHEP(NHEP+2)=99
          IF(MOD(ABS(DecayChains(idecay,1)),4).EQ.0.OR.&
               MOD(ABS(DecayChains(idecay,1)),4).EQ.3)THEN
             IF(DecayChains(idecay,1).GT.0)THEN
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+1)=0
                ICOLHEP(1,NHEP+2)=ho_colmax+101
                ICOLHEP(2,NHEP+2)=0
                ho_colmax=ho_colmax+1
             ELSE
                ICOLHEP(2,NHEP+1)=0
                ICOLHEP(1,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+2)=0
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ho_colmax=ho_colmax+1
             ENDIF
          ELSE
             ICOLHEP(1:2,NHEP+1:NHEP+2)=0
          ENDIF
          NHEP=NHEP+2
       CASE(4)
          ! t > b w+ or t~ > b~ w-
          IF(emep_ISR_shower.EQ.1)THEN
             ! isr in e-e+ collisions
             PM(0)=PLPTN(7,imoth-2+NVPH_QEDPS)
             PM(1:3)=PLPTN(4:6,imoth-2+NVPH_QEDPS)
          ELSE
             PM(0)=hadron_pmom(imoth,4)
             PM(1:3)=hadron_pmom(imoth,1:3)
          END IF
          svec(1:3)=PM(1:3)
          isign=2*ipol(imoth)-3
          isign=-isign
          IF(iflh(imoth).GT.0)THEN
             anti=1
          ELSE
             anti=-1
          ENDIF
          IF(emep_ISR_shower.EQ.1)THEN
             JMOHEP(1,NHEP+1:NHEP+2)=imoth+NVPH_QEDPS-3
          ELSE
             JMOHEP(1,NHEP+1:NHEP+2)=imoth
          ENDIF
          JMOHEP(2,NHEP+1:NHEP+2)=0
          CALL HO_t2bw_hel(PM,PD1,PD2,svec,isign,iheld1,iheld2,anti)
          IF(ipol(imoth).EQ.1)THEN
             IPOLHEP(JMOHEP(1,NHEP+1))=1
          ELSE
             IPOLHEP(JMOHEP(1,NHEP+1))=-1
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          IF(ABS(DecayChains(idecay,1)).GT.ABS(DecayChains(idecay,2)))THEN
             ! w b
             PHEP(5,NHEP+1)=parmas(33)
             PHEP(5,NHEP+2)=0d0
             PHEP(1:3,NHEP+1)=PD2(1:3)
             PHEP(1:3,NHEP+2)=PD1(1:3)
             PHEP(4,NHEP+1)=PD2(0)
             PHEP(4,NHEP+2)=PD1(0)
             IPOLHEP(NHEP+1)=iheld2 ! decay product will not decay anymore
             IPOLHEP(NHEP+2)=iheld1
             ICOLHEP(1:2,NHEP+1)=0
             ICOLHEP(1:2,NHEP+2)=ICOLHEP(1:2,JMOHEP(1,NHEP+1))
          ELSE
             PHEP(5,NHEP+2)=parmas(33)
             PHEP(5,NHEP+1)=0d0
             PHEP(1:3,NHEP+1)=PD1(1:3)
             PHEP(1:3,NHEP+2)=PD2(1:3)
             PHEP(4,NHEP+1)=PD1(0)
             PHEP(4,NHEP+2)=PD2(0)
             IPOLHEP(NHEP+1)=iheld1 ! decay product will not decay anymore
             IPOLHEP(NHEP+2)=iheld2
             ICOLHEP(1:2,NHEP+1)=ICOLHEP(1:2,JMOHEP(1,NHEP+1))
             ICOLHEP(1:2,NHEP+2)=0
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          NHEP=NHEP+2
       CASE(5,6)
          ! chi1,chi2 > psi + gamma
          CALL Qeveryhel(index2,ipol(imoth),imoth)
          IF(emep_ISR_shower.EQ.1)THEN
             ! isr in e-e+ collisions
             PM(0)=PLPTN(7,imoth-2+NVPH_QEDPS)
             PM(1:3)=PLPTN(4:6,imoth-2+NVPH_QEDPS)
          ELSE
             PM(0)=hadron_pmom(imoth,4)
             PM(1:3)=hadron_pmom(imoth,1:3)
          END IF
          svec(1:3)=PM(1:3)
          isign=2*index2(3)-3
          isign=-isign
          isign=isign-(2*index2(2)-3)
          IF(iflh(imoth).GT.440000.AND.iflh(imoth).LT.449999)THEN
             ! charmonia decay
             IF(iflh(imoth).EQ.443111)THEN
                mchi=ho_dmass_chic1
             ELSE
                mchi=ho_dmass_chic2
             ENDIF
             mpsi=ho_dmass_jpsi
          ELSE
             ! bottonia decay,only restricted to 1P -> 1S
             IF(iflh(imoth).EQ.553111)THEN
                mchi=ho_dmass_chib11P
             ELSE
                mchi=ho_dmass_chib21P
             ENDIF
             mpsi=ho_dmass_Y1S
          ENDIF
          IF(emep_ISR_shower.EQ.1)THEN
             JMOHEP(1,NHEP+1:NHEP+2)=imoth+NVPH_QEDPS-3
          ELSE
             JMOHEP(1,NHEP+1:NHEP+2)=imoth
          ENDIF
          JMOHEP(2,NHEP+1:NHEP+2)=0
          ! test decay angular distributions
          !rr=Helac_rnmy(0)
          !IF(rr.LE.2d0/6d0)THEN
          !   isign=2
          !ELSEIF(rr.GT.2d0/6d0.AND.rr.LE.4d0/6d0)THEN
          !   isign=1
          !ELSE
          !   isign=0
          !ENDIF
          !isign=0
          !lavail=6
          ! end of test
          IF(lavail.EQ.5)THEN
             CALL HO_chi12psia_hel(PM,PD1,PD2,svec,isign,iheld1,iheld2,mchi,mpsi)
          ELSE
             CALL HO_chi22psia_hel(PM,PD1,PD2,svec,isign,iheld1,iheld2,mchi,mpsi)
          ENDIF
          IPOLHEP(JMOHEP(1,NHEP+1))=isign
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          IF(ABS(DecayChains(idecay,1)).GT.100)THEN
             ! psi gamma
             PHEP(5,NHEP+2)=0d0
             PHEP(1:3,NHEP+1)=PD1(1:3)
             PHEP(1:3,NHEP+2)=PD2(1:3)
             PHEP(4,NHEP+1)=PD1(0)
             PHEP(4,NHEP+2)=PD2(0)
             PHEP(5,NHEP+1)=DSQRT(PD1(0)**2-PD1(1)**2-PD1(2)**2-PD1(3)**2)
             IPOLHEP(NHEP+1)=iheld1 ! decay product will decay
             IPOLHEP(NHEP+2)=iheld2 ! decay product will decay
             ICOLHEP(1:2,NHEP+1)=0
             ICOLHEP(1:2,NHEP+2)=0
          ELSE
             ! gamma psi
             PHEP(5,NHEP+1)=0d0
             PHEP(1:3,NHEP+1)=PD2(1:3)
             PHEP(1:3,NHEP+2)=PD1(1:3)
             PHEP(4,NHEP+1)=PD2(0)
             PHEP(4,NHEP+2)=PD1(0)
             PHEP(5,NHEP+2)=DSQRT(PD1(0)**2-PD1(1)**2-PD1(2)**2-PD1(3)**2)
             IPOLHEP(NHEP+1)=iheld2 ! decay product will decay
             IPOLHEP(NHEP+2)=iheld1
             ICOLHEP(1:2,NHEP+1)=0
             ICOLHEP(1:2,NHEP+2)=0
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          NHEP=NHEP+2
       END SELECT
    ELSE
       CALL HO_avail_decays(available_ndecay,available_decay)
       PRINT *,"Only the following decay processes are available in HELAC-Onia"
       DO i=1,available_ndecay
          PRINT *,available_decay(i)
       ENDDO
       STOP
    ENDIF
    RETURN
  END SUBROUTINE HO_One_Decay

  SUBROUTINE HO_One_Decay_next(imoth,idecay)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::imoth,idecay
    INTEGER,PARAMETER::available_ndecay=5
    CHARACTER(len=20),DIMENSION(available_ndecay)::available_decay
    INTEGER::i,isign,iheld1,iheld2,anti
    INTEGER::lavail,imm
    REAL(KIND(1d0)),DIMENSION(0:3)::PM,PD1,PD2
    REAL(KIND(1d0)),DIMENSION(4)::PMM,PM2,PBOO
    REAL(KIND(1d0)),DIMENSION(3)::svec
    REAL(KIND(1d0))::YL,YR,QMM,mchi,mpsi
    INCLUDE "stdhep.inc"
    CALL HO_judge_avail_decay_next(imoth,idecay,lavail)
    IF(lavail.GT.0)THEN
       SELECT CASE(lavail)
       CASE(1)
          ! 3S11 > l+ l-
          ! mother 3S11's mother
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! 3S11 is coming from decay
             imm=JMOHEP(1,imoth)
             PMM(1:4)=PHEP(1:4,imm)
             PBOO(1:3)=-PMM(1:3)
             PBOO(4)=PMM(4)
             QMM=DSQRT(ABS(PMM(4)**2-PMM(1)**2-PMM(2)**2-PMM(3)**2))
             PM2(1:4)=PHEP(1:4,imoth)
             CALL BOOSTL(QMM,PBOO,PM2)
             ! boost to the rest frame of 3S11's mother particle in the decaychain
          ELSE
             PM2(1:4)=PHEP(1:4,imoth)
          ENDIF
          PM(0)=PM2(4)
          PM(1:3)=PM2(1:3)
          svec(1:3)=PM2(1:3)
          ! KEEP IN MIND:
          ! THE FOLLOWING MOMENTA PM,PD1,PD2 are in PMM rest frame
          ! until the boostl
          IF(IPOLHEP(imoth).NE.0)THEN
             IF(IPOLHEP(imoth).EQ.99)THEN
                WRITE(*,*)"ERROR:Cannot decay helicity-summed 3S11 onium"
                STOP
             ENDIF
             ! transverse polarized in HX
             CALL HO_Vll11(PM,svec,PD1,PD2)
          ELSE
             ! longitudial polarized in HX
             CALL HO_Vll00(PM,svec,PD1,PD2)
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JMOHEP(1,NHEP+1:NHEP+2)=imoth
          JMOHEP(2,NHEP+1:NHEP+2)=0
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          PHEP(5,NHEP+1:NHEP+2)=0d0
          PHEP(1:3,NHEP+1)=PD1(1:3)
          PHEP(1:3,NHEP+2)=PD2(1:3)
          PHEP(4,NHEP+1)=PD1(0)
          PHEP(4,NHEP+2)=PD2(0)
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! boost back to lab frame
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+1))
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+2))
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          IPOLHEP(NHEP+1)=99 ! decay products will not decay anymore
          IPOLHEP(NHEP+2)=99
          ICOLHEP(1:2,NHEP+1:NHEP+2)=0
          NHEP=NHEP+2
       CASE(2)
          ! W+- > f1 f2
          ! mother W's mother
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! W+- is coming from decay
             imm=JMOHEP(1,imoth)
             PMM(1:4)=PHEP(1:4,imm)
             PBOO(1:3)=-PMM(1:3)
             PBOO(4)=PMM(4)
             QMM=DSQRT(ABS(PMM(4)**2-PMM(1)**2-PMM(2)**2-PMM(3)**2))
             PM2(1:4)=PHEP(1:4,imoth)
             CALL BOOSTL(QMM,PBOO,PM2)
             ! boost to the rest frame of W's mother particle in the decaychain
          ELSE
             PM2(1:4)=PHEP(1:4,imoth)
          ENDIF
          PM(0)=PM2(4)
          PM(1:3)=PM2(1:3)
          svec(1:3)=PM(1:3)
          ! KEEP IN MIND:
          ! THE FOLLOWING MOMENTA PM,PD1,PD2 are in PMM rest frame
          ! until the boostl 
          IF(IPOLHEP(imoth).NE.0)THEN
             IF(IPOLHEP(imoth).EQ.99)THEN
                WRITE(*,*)"ERROR:cannot decay helicity-summed W boson"
                STOP
             ENDIF
             ! transverse
             isign=-IPOLHEP(imoth)
             IF(DecayChains(idecay,1).GT.0)THEN
                isign=-isign
             ENDIF
             CALL HO_VllW11(PM,svec,isign,PD1,PD2)
          ELSE
             ! longitudial
             CALL HO_Vll00(PM,svec,PD1,PD2)
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JMOHEP(1,NHEP+1:NHEP+2)=imoth
          JMOHEP(2,NHEP+1:NHEP+2)=0
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          PHEP(5,NHEP+1:NHEP+2)=0d0
          PHEP(1:3,NHEP+1)=PD1(1:3)
          PHEP(1:3,NHEP+2)=PD2(1:3)
          PHEP(4,NHEP+1)=PD1(0)
          PHEP(4,NHEP+2)=PD2(0)
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! boost back to lab frame
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+1))
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+2))
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          IPOLHEP(NHEP+1)=99 ! decay product will not decay anymore
          IPOLHEP(NHEP+2)=99
          IF(MOD(ABS(DecayChains(idecay,1)),4).EQ.0.OR.&
               MOD(ABS(DecayChains(idecay,1)),4).EQ.3)THEN
             IF(DecayChains(idecay,1).GT.0)THEN
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+1)=0
                ICOLHEP(1,NHEP+2)=ho_colmax+101
                ICOLHEP(2,NHEP+2)=0
                ho_colmax=ho_colmax+1
             ELSE
                ICOLHEP(2,NHEP+1)=0
                ICOLHEP(1,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+2)=0
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ho_colmax=ho_colmax+1
             ENDIF
          ELSE
             ICOLHEP(1:2,NHEP+1:NHEP+2)=0
          ENDIF
          NHEP=NHEP+2
       CASE(3)
          ! Z > f f~
          ! mother Z's mother
          IF(JMOHEP(2,imoth).EQ.0)THEN
             imm=JMOHEP(1,imoth)
             PMM(1:4)=PHEP(1:4,imm)
             PBOO(1:3)=-PMM(1:3)
             PBOO(4)=PMM(4)
             QMM=DSQRT(ABS(PMM(4)**2-PMM(1)**2-PMM(2)**2-PMM(3)**2))
             PM2(1:4)=PHEP(1:4,imoth)
             CALL BOOSTL(QMM,PBOO,PM2)
             ! boost to the rest frame of Z's mother particle in the decaychain
          ELSE
             PM2(1:4)=PHEP(1:4,imoth)
          ENDIF
          PM(0)=PM2(4)
          PM(1:3)=PM2(1:3)
          svec(1:3)=PM(1:3)
          ! KEEP IN MIND:
          ! THE FOLLOWING MOMENTA PM,PD1,PD2 are in PMM rest frame
          ! until the boostl
          PM(0)=PHEP(4,imoth)
          PM(1:3)=PHEP(1:3,imoth)
          svec(1:3)=PM(1:3)
          IF(IPOLHEP(imoth).NE.0)THEN
             IF(IPOLHEP(imoth).EQ.99)THEN
                WRITE(*,*)"ERROR:Cannot decay helicity-summed Z boson"
                STOP
             ENDIF
             ! transverse
             isign=-IPOLHEP(imoth)
             IF(DecayChains(idecay,1).GT.0)THEN
                isign=-isign
             ENDIF
             CALL Coup_Zff(DecayChains(idecay,1),YL,YR)
             CALL HO_Zff11(PM,svec,isign,YL,YR,PD1,PD2)
          ELSE
             ! longitudial
             CALL HO_Vll00(PM,svec,PD1,PD2)
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JMOHEP(1,NHEP+1:NHEP+2)=imoth
          JMOHEP(2,NHEP+1:NHEP+2)=0
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          PHEP(5,NHEP+1:NHEP+2)=0d0
          PHEP(1:3,NHEP+1)=PD1(1:3)
          PHEP(1:3,NHEP+2)=PD2(1:3)
          PHEP(4,NHEP+1)=PD1(0)
          PHEP(4,NHEP+2)=PD2(0)
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! boost back to lab frame
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+1))
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+2))
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          IPOLHEP(NHEP+1)=99 ! decay product will not decay anymore
          IPOLHEP(NHEP+2)=99
          IF(MOD(ABS(DecayChains(idecay,1)),4).EQ.0.OR.&
               MOD(ABS(DecayChains(idecay,1)),4).EQ.3)THEN
             IF(DecayChains(idecay,1).GT.0)THEN
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+1)=0
                ICOLHEP(1,NHEP+2)=ho_colmax+101
                ICOLHEP(2,NHEP+2)=0
                ho_colmax=ho_colmax+1
             ELSE
                ICOLHEP(2,NHEP+1)=0
                ICOLHEP(1,NHEP+1)=ho_colmax+101
                ICOLHEP(1,NHEP+2)=0
                ICOLHEP(2,NHEP+1)=ho_colmax+101
                ho_colmax=ho_colmax+1
             ENDIF
          ELSE
             ICOLHEP(1:2,NHEP+1:NHEP+2)=0
          ENDIF
          NHEP=NHEP+2
       CASE(4)
          ! t > b w+ or t~ > b~ w-
          ! mother t's mother
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! t is coming from another decay chain
             imm=JMOHEP(1,imoth)
             PMM(1:4)=PHEP(1:4,imm)
             PBOO(1:3)=-PMM(1:3)
             PBOO(4)=PMM(4)
             QMM=DSQRT(ABS(PMM(4)**2-PMM(1)**2-PMM(2)**2-PMM(3)**2))
             PM2(1:4)=PHEP(1:4,imoth)
             CALL BOOSTL(QMM,PBOO,PM2)
             ! boost to the rest frame of t's mother particle in the decay chain
          ELSE
             PM2(1:4)=PHEP(1:4,imoth)
          ENDIF
          PM(0)=PM2(4)
          PM(1:3)=PM2(1:3)
          svec(1:3)=PM(1:3)
          ! KEEP IN MIND:
          ! THE FOLLOWING MOMENTA PM,PD1,PD2 are in PMM rest frame
          ! until the boostl
          IF(IPOLHEP(imoth).EQ.99)THEN
             WRITE(*,*)"ERROR: cannot decay helicity-summed top quark"
             STOP
          ENDIF
          isign=IPOLHEP(imoth)
          IF(IDHEP(imoth).GT.0)THEN
             anti=1
          ELSE
             anti=-1
          ENDIF
          CALL HO_t2bw_hel(PM,PD1,PD2,svec,isign,iheld1,iheld2,anti)
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JMOHEP(1,NHEP+1:NHEP+2)=imoth
          JMOHEP(2,NHEP+1:NHEP+2)=0
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          IF(ABS(DecayChains(idecay,1)).GT.ABS(DecayChains(idecay,2)))THEN
             ! w b
             PHEP(5,NHEP+1)=parmas(33)
             PHEP(5,NHEP+2)=0d0
             PHEP(1:3,NHEP+1)=PD2(1:3)
             PHEP(1:3,NHEP+2)=PD1(1:3)
             PHEP(4,NHEP+1)=PD2(0)
             PHEP(4,NHEP+2)=PD1(0)
             IPOLHEP(NHEP+1)=iheld2 ! decay product will not decay anymore
             IPOLHEP(NHEP+2)=iheld1
             ICOLHEP(1:2,NHEP+1)=0
             ICOLHEP(1:2,NHEP+2)=ICOLHEP(1:2,imoth)
          ELSE
             PHEP(5,NHEP+2)=parmas(33)
             PHEP(5,NHEP+1)=0d0
             PHEP(1:3,NHEP+1)=PD1(1:3)
             PHEP(1:3,NHEP+2)=PD2(1:3)
             PHEP(4,NHEP+1)=PD1(0)
             PHEP(4,NHEP+2)=PD2(0)
             IPOLHEP(NHEP+1)=iheld1 ! decay product will not decay anymore
             IPOLHEP(NHEP+2)=iheld2
             ICOLHEP(1:2,NHEP+1)=ICOLHEP(1:2,imoth)
             ICOLHEP(1:2,NHEP+2)=0
          ENDIF
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! boost back to lab frame
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+1))
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+2))
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          NHEP=NHEP+2
       CASE(5,6)
          ! chi1,chi2 > psi gamma
          ! mother chi1's mother
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! chi12 is coming from another decay
             imm=JMOHEP(1,imoth)
             PMM(1:4)=PHEP(1:4,imm)
             PBOO(1:3)=-PMM(1:3)
             PBOO(4)=PMM(4)
             QMM=DSQRT(ABS(PMM(4)**2-PMM(1)**2-PMM(2)**2-PMM(3)**2))
             PM2(1:4)=PHEP(1:4,imoth)
             CALL BOOSTL(QMM,PBOO,PM2)
             ! boost to the rest frame of chi1's mother particle in the decay chain
          ELSE
             ! chi12 is coming from hard generation
             PM2(1:4)=PHEP(1:4,imoth)
          ENDIF
          PM(0)=PM2(4)
          PM(1:3)=PM2(1:3)
          svec(1:3)=PM(1:3)
          ! KEEP IN MIND:
          ! THE FOLLOWING MOMENTA PM,PD1,PD2 are in PMM rest frame
          ! until the boostl
          IF(IPOLHEP(imoth).EQ.99)THEN
             WRITE(*,*)"ERROR: cannot decay helicity-summed chi1 meson"
             STOP
          ENDIF
          isign=IPOLHEP(imoth)
          IF(IDHEP(imoth).EQ.20443)THEN
             mchi=ho_dmass_chic1
             mpsi=ho_dmass_jpsi
          ELSEIF(IDHEP(imoth).EQ.445)THEN
             mchi=ho_dmass_chic2
             mpsi=ho_dmass_jpsi
          ELSEIF(IDHEP(imoth).EQ.555)THEN
             mchi=ho_dmass_chib21P
             mpsi=ho_dmass_Y1S
          ELSE
             mchi=ho_dmass_chib11P
             mpsi=ho_dmass_Y1S
          ENDIF
          IF(lavail.EQ.5)THEN
             CALL HO_chi12psia_hel(PM,PD1,PD2,svec,isign,iheld1,iheld2,mchi,mpsi)
          ELSE
             CALL HO_chi22psia_hel(PM,PD1,PD2,svec,isign,iheld1,iheld2,mchi,mpsi)
          ENDIF
          ISTHEP(imoth)=1
          ISHEP(NHEP+1)=1
          ISHEP(NHEP+2)=1
          ISTHEP(NHEP+1)=0
          ISTHEP(NHEP+2)=0
          IDHEP(NHEP+1)=pdgt(DecayChains(idecay,1))
          IDHEP(NHEP+2)=pdgt(DecayChains(idecay,2))
          JMOHEP(1,NHEP+1:NHEP+2)=imoth
          JMOHEP(2,NHEP+1:NHEP+2)=0
          JDAHEP(1,imoth)=NHEP+1
          JDAHEP(2,imoth)=NHEP+2
          IF(ABS(DecayChains(idecay,1)).GT.100)THEN
             ! psi gamma
             PHEP(5,NHEP+2)=0d0
             PHEP(1:3,NHEP+1)=PD1(1:3)
             PHEP(1:3,NHEP+2)=PD2(1:3)
             PHEP(4,NHEP+1)=PD1(0)
             PHEP(4,NHEP+2)=PD2(0)
             PHEP(5,NHEP+1)=DSQRT(PD1(0)**2-PD1(1)**2-PD1(2)**2-PD1(3)**2)
             IPOLHEP(NHEP+1)=iheld1 ! decay product will decay                   
             IPOLHEP(NHEP+2)=iheld2
             ICOLHEP(1:2,NHEP+1)=0
             ICOLHEP(1:2,NHEP+2)=0
          ELSE
             ! gamma psi
             PHEP(5,NHEP+1)=0d0
             PHEP(1:3,NHEP+1)=PD2(1:3)
             PHEP(1:3,NHEP+2)=PD1(1:3)
             PHEP(4,NHEP+1)=PD2(0)
             PHEP(4,NHEP+2)=PD1(0)
             PHEP(5,NHEP+2)=DSQRT(PD1(0)**2-PD1(1)**2-PD1(2)**2-PD1(3)**2)
             IPOLHEP(NHEP+1)=iheld2 ! decay product will decay
             IPOLHEP(NHEP+2)=iheld1
             ICOLHEP(1:2,NHEP+1)=0
             ICOLHEP(1:2,NHEP+2)=0
          ENDIF
          IF(JMOHEP(2,imoth).EQ.0)THEN
             ! chi12 is coming from another decay
             ! boost back to lab frame
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+1))
             CALL BOOSTL(QMM,PMM,PHEP(1:4,NHEP+2))
          ENDIF
          VHEP(1:4,NHEP+1:NHEP+2)=0d0
          NHEP=NHEP+2
       END SELECT
    ELSE
       CALL HO_avail_decays(available_ndecay,available_decay)
       PRINT *,"Only the following decay processes are available in HELAC-Onia"
       DO i=1,available_ndecay
          PRINT *,available_decay(i)
       ENDDO
       STOP
    ENDIF
    RETURN
  END SUBROUTINE HO_One_Decay_next

  SUBROUTINE HO_avail_decays(nnn,char_decay)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::nnn
    CHARACTER(len=20),DIMENSION(nnn),INTENT(OUT)::char_decay
    IF(nnn.LE.0)RETURN
    char_decay(1)="3S11 > l+ l-"
    IF(nnn.LE.1)RETURN
    char_decay(2)="W+- > f1 f2"
    IF(nnn.LE.2)RETURN
    char_decay(3)="Z > f f~"
    IF(nnn.LE.3)RETURN
    char_decay(4)="t > b w+ | t~ > b~ w-"
    IF(nnn.LE.4)RETURN
    char_decay(5)="chic1 > psi a | chib1 > Y a"
    IF(nnn.LE.5)RETURN
    char_decay(6)="chic2 > psi a | chib2 > Y a"
    RETURN
  END SUBROUTINE HO_avail_decays

  SUBROUTINE HO_judge_avail_decay(imoth,idecay,lavail)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::imoth,idecay
    INTEGER,INTENT(OUT)::lavail
    INTEGER::ll,sum,i,f1,f2
    lavail=0
    IF(iflh(imoth).EQ.443011.OR.iflh(imoth).EQ.553011)THEN
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       sum=0
       DO i=1,ll
          sum=sum+DecayChains(idecay,i)
       ENDDO
       IF(sum.NE.0)RETURN
       IF(ABS(DecayChains(idecay,1)).NE.2.AND.ABS(DecayChains(idecay,1)).NE.6)RETURN
       lavail=1
       RETURN
    ELSEIF(iflh(imoth).EQ.33.OR.iflh(imoth).EQ.34)THEN
       ! W+- > f1 f2
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       ! only include light fermions
       IF(ABS(f1).GT.12.OR.ABS(f1).EQ.11)RETURN
       f2=DecayChains(idecay,2)
       IF(ABS(f2).GT.12.OR.ABS(f2).EQ.11)RETURN
       IF(f1*f2.GT.0)RETURN
       i=MOD(ABS(f1),4)
       IF(MOD(ABS(f2),4)+i.NE.3)RETURN
       IF(i.EQ.1.OR.i.EQ.2)THEN
          ! leptonic decay
          IF(ABS(ABS(f1)-ABS(f2)).NE.1)RETURN ! no generation mixing
          ! charge conservation
          IF(f1.GT.0)THEN
             IF(i.EQ.1.AND.iflh(imoth).NE.33)RETURN
             IF(i.EQ.2.AND.iflh(imoth).NE.34)RETURN
          ELSE
             IF(i.EQ.1.AND.iflh(imoth).NE.34)RETURN
             IF(i.EQ.2.AND.iflh(imoth).NE.33)RETURN
          ENDIF
       ELSE
          ! hadronic decay
          ! charge conservation
          IF(f1.GT.0)THEN
             IF(i.EQ.3.AND.iflh(imoth).NE.33)RETURN
             IF(i.EQ.0.AND.iflh(imoth).NE.34)RETURN
          ELSE
             IF(i.EQ.3.AND.iflh(imoth).NE.34)RETURN
             IF(i.EQ.0.AND.iflh(imoth).NE.33)RETURN
          ENDIF
       ENDIF
       lavail=2
       RETURN
    ELSEIF(iflh(imoth).EQ.32)THEN
       ! Z > f f~
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       ! only include light fermions
       IF(ABS(f1).GT.12.OR.ABS(f1).EQ.11)RETURN
       f2=DecayChains(idecay,2)
       IF(ABS(f2).GT.12.OR.ABS(f2).EQ.11)RETURN
       IF(f1+f2.NE.0)RETURN
       lavail=3
       RETURN
    ELSEIF(ABS(iflh(imoth)).EQ.11)THEN
       ! t > b W+ or t~ > b~ W-
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       IF(ABS(f1).NE.12.OR.f1.NE.33.OR.f1.NE.34)RETURN
       f2=DecayChains(idecay,2)
       IF(ABS(f2).NE.12.OR.f2.NE.33.OR.f2.NE.34)RETURN
       IF(ABS(f1).EQ.12.AND.(f2.NE.33.OR.f2.NE.34))RETURN
       IF(ABS(f2).EQ.12.AND.(f1.NE.33.OR.f1.NE.34))RETURN
       IF(f1*f2*iflh(imoth).LT.0)RETURN
       IF(iflh(imoth).GT.0)THEN
          IF(f1.NE.33.AND.f2.NE.33)RETURN
          IF(f1.NE.12.AND.f2.NE.12)RETURN
       ELSE
          IF(f1.NE.34.AND.f2.NE.34)RETURN
          IF(f1.NE.-12.AND.f2.NE.-12)RETURN
       ENDIF
       lavail=4
       RETURN
    ELSEIF(iflh(imoth).EQ.443111.OR.iflh(imoth).EQ.553111)THEN
       ! chi1 > psi a
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       f2=DecayChains(idecay,2)
       IF(f1.NE.31.AND.f2.NE.31)RETURN
       IF(iflh(imoth).EQ.443111.AND.f1.NE.443011.AND.f2.NE.443011)RETURN
       IF(iflh(imoth).EQ.553111.AND.f1.NE.553011.AND.f2.NE.553011)RETURN
       lavail=5
       RETURN
    ELSEIF(iflh(imoth).EQ.443121.OR.iflh(imoth).EQ.553121)THEN
       ! chi2 > psi a
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       f2=DecayChains(idecay,2)
       IF(f1.NE.31.AND.f2.NE.31)RETURN
       IF(iflh(imoth).EQ.443121.AND.f1.NE.443011.AND.f2.NE.443011)RETURN
       IF(iflh(imoth).EQ.553121.AND.f1.NE.553011.AND.f2.NE.553011)RETURN
       lavail=6
       RETURN
    ENDIF
    RETURN
  END SUBROUTINE HO_judge_avail_decay

  SUBROUTINE HO_judge_avail_decay_next(imoth,idecay,lavail)
    IMPLICIT NONE
    INTEGER,INTENT(IN)::imoth,idecay
    INTEGER,INTENT(OUT)::lavail
    INTEGER::ll,sum,i,f1,f2
    INCLUDE "stdhep.inc"
    lavail=0
    IF(IDHEP(imoth).EQ.pdgt(443011).OR.IDHEP(imoth).EQ.pdgt(553011))THEN
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       sum=0
       DO i=1,ll
          sum=sum+DecayChains(idecay,i)
       ENDDO
       IF(sum.NE.0)RETURN
       IF(ABS(DecayChains(idecay,1)).NE.2.AND.ABS(DecayChains(idecay,1)).NE.6)RETURN
       lavail=1
       RETURN
    ELSEIF(IDHEP(imoth).EQ.pdgt(33).OR.IDHEP(imoth).EQ.pdgt(34))THEN
       ! W+- > f1 f2
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       ! only include light fermions
       IF(ABS(f1).GT.12.OR.ABS(f1).EQ.11)RETURN
       f2=DecayChains(idecay,2)
       IF(ABS(f2).GT.12.OR.ABS(f2).EQ.11)RETURN
       IF(f1*f2.GT.0)RETURN
       i=MOD(ABS(f1),4)
       IF(MOD(ABS(f2),4)+i.NE.3)RETURN
       IF(i.EQ.1.OR.i.EQ.2)THEN
          ! leptonic decay
          IF(ABS(ABS(f1)-ABS(f2)).NE.1)RETURN ! no generation mixing
          ! charge conservation
          IF(f1.GT.0)THEN
             IF(i.EQ.1.AND.iflh(imoth).NE.33)RETURN
             IF(i.EQ.2.AND.iflh(imoth).NE.34)RETURN
          ELSE
             IF(i.EQ.1.AND.iflh(imoth).NE.34)RETURN
             IF(i.EQ.2.AND.iflh(imoth).NE.33)RETURN
          ENDIF
       ELSE
          ! hadronic decay
          ! charge conservation
          IF(f1.GT.0)THEN
             IF(i.EQ.3.AND.iflh(imoth).NE.33)RETURN
             IF(i.EQ.0.AND.iflh(imoth).NE.34)RETURN
          ELSE
             IF(i.EQ.3.AND.iflh(imoth).NE.34)RETURN
             IF(i.EQ.0.AND.iflh(imoth).NE.33)RETURN
          ENDIF
       ENDIF
       lavail=2
       RETURN
   ELSEIF(IDHEP(imoth).EQ.pdgt(32))THEN
       ! Z > f f~
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       ! only include light fermions
       IF(ABS(f1).GT.12.OR.ABS(f1).EQ.11)RETURN
       f2=DecayChains(idecay,2)
       IF(ABS(f2).GT.12.OR.ABS(f2).EQ.11)RETURN
       IF(f1+f2.NE.0)RETURN
       lavail=3
       RETURN
    ELSEIF(ABS(IDHEP(imoth)).EQ.pdgt(11))THEN
       ! t > b W+ or t~ > b~ W-
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       IF(ABS(f1).NE.12.OR.f1.NE.33.OR.f1.NE.34)RETURN
       f2=DecayChains(idecay,2)
       IF(ABS(f2).NE.12.OR.f2.NE.33.OR.f2.NE.34)RETURN
       IF(ABS(f1).EQ.12.AND.(f2.NE.33.OR.f2.NE.34))RETURN
       IF(ABS(f2).EQ.12.AND.(f1.NE.33.OR.f1.NE.34))RETURN
       IF(f1*f2*IDHEP(imoth).LT.0)RETURN
       IF(IDHEP(imoth).GT.0)THEN
          IF(f1.NE.33.AND.f2.NE.33)RETURN
          IF(f1.NE.12.AND.f2.NE.12)RETURN
       ELSE
          IF(f1.NE.34.AND.f2.NE.34)RETURN
          IF(f1.NE.-12.AND.f2.NE.-12)RETURN
       ENDIF
       lavail=4
       RETURN
    ELSEIF(IDHEP(imoth).EQ.pdgt(443111).OR.IDHEP(imoth).EQ.pdgt(553111))THEN
       ! chi1 > psi a
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       f2=DecayChains(idecay,2)
       IF(f1.NE.31.AND.f2.NE.31)RETURN
       IF(IDHEP(imoth).EQ.pdgt(443111).AND.f1.NE.443011.AND.f2.NE.443011)RETURN
       IF(IDHEP(imoth).EQ.pdgt(553111).AND.f1.NE.553011.AND.f2.NE.553011)RETURN
       lavail=5
       RETURN
    ELSEIF(IDHEP(imoth).EQ.pdgt(443121).OR.IDHEP(imoth).EQ.pdgt(553121))THEN
       ! chi2 > psi a
       ll=DecayChains(idecay,-1)
       IF(ll.NE.2)RETURN
       f1=DecayChains(idecay,1)
       f2=DecayChains(idecay,2)
       IF(f1.NE.31.AND.f2.NE.31)RETURN
       IF(IDHEP(imoth).EQ.pdgt(443121).AND.f1.NE.443011.AND.f2.NE.443011)RETURN
       IF(IDHEP(imoth).EQ.pdgt(553121).AND.f1.NE.553011.AND.f2.NE.553011)RETURN
       lavail=6
       RETURN
    ENDIF
    RETURN
  END SUBROUTINE HO_judge_avail_decay_next

  SUBROUTINE partonmom2hadronmom
    IMPLICIT NONE
    INTEGER::i,j
    REAL(KIND(1d0)),DIMENSION(20,4)::p
    REAL(KIND(1d0)),DIMENSION(4)::pboo
    REAL(KIND(1d0))::exp1,exp2
    REAL(KIND(1d0))::m1,m2,q,e
    DO i=1,n
       p(i,1:4)=phegas_pmom(i,1:4)
    ENDDO
    ! we have boosted it to the lab frame with decay chains
    !IF((istruc.EQ.1.OR..NOT.labeqcoll).AND.imode.EQ.0)THEN
    !   q=ehat
    !   e=q/SQRT(xp1*xp2)
    !   exp1=ebeam(1)*xp1
    !   exp2=ebeam(2)*xp2
    !   pboo(4)=(exp1+exp2)
    !   pboo(3)=(exp1-exp2)
    !   pboo(1:2)=0
       ! boost to the lab frame                                                        
    !   DO i=1,n
    !      CALL boostl(q,pboo,p(i,1:4))
    !   ENDDO
    !ENDIF
    j=1
    DO i=1,nhad
       IF(ABS(iflh(i)).GT.100)THEN
          m1=phegas_pmom(j,5)
          m2=phegas_pmom(j+1,5)
          hadron_pmom(i,1:4)=p(j,1:4)+p(j+1,1:4)
          hadron_pmom(i,5)=m1+m2
          j=j+2
       ELSE
          hadron_pmom(i,1:4)=p(j,1:4)
          hadron_pmom(i,5)=phegas_pmom(j,5)
          j=j+1
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE partonmom2hadronmom

  SUBROUTINE unwei_writer_Decay
    USE setscale
    IMPLICIT NONE
    INCLUDE "stdhep.inc"
    REAL(KIND(1d0))::scale1
    INTEGER::icol1,icol2,idup,idprup,istup,imothup1,imothup2
    REAL(KIND(1d0))::px,py,pz,p0,pmass,scalup,vtime,vspin,xwgtup
    INTEGER::i
    idprup = 81
    xwgtup=1
    CALL qcdscale(scale1)
    scalup=scale1
    WRITE(nunit3)NHEP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP
    DO i=1,NHEP
       idup=IDHEP(i)
       istup=ISHEP(i)
       imothup1=JMOHEP(1,i)
       imothup2=JMOHEP(2,i)
       IF(emep_ISR_shower.EQ.1)THEN
          IF(i.LT.NVPH_QEDPS)THEN
             icol1=0
             icol2=0
          ELSE
             IF(i-NVPH_QEDPS+3.GT.nhad)THEN
                icol1=ICOLHEP(1,i)
                icol2=ICOLHEP(2,i)
             ELSE
                icol1=icol_un(i-NVPH_QEDPS+3,2)+100
                IF(icol1.EQ.100)icol1=0
                icol2=icol_un(i-NVPH_QEDPS+3,1)+100
                IF(icol2.EQ.100)icol2=0
             ENDIF
          ENDIF
       ELSE
          IF(i.GT.nhad)THEN
             icol1=ICOLHEP(1,i)
             icol2=ICOLHEP(2,i)
          ELSE
             IF(istup.LT.0)THEN
                icol1=icol_un(i,1)+100
                IF(icol1.EQ.100)icol1=0
                icol2=icol_un(i,2)+100
                IF(icol2.EQ.100)icol2=0
             ELSE
                icol1=icol_un(i,2)+100
                IF(icol1.EQ.100)icol1=0
                icol2=icol_un(i,1)+100
                IF(icol2.EQ.100)icol2=0
             ENDIF
          ENDIF
       ENDIF
       px=PHEP(1,i)
       py=PHEP(2,i)
       pz=PHEP(3,i)
       p0=PHEP(4,i)
       pmass=PHEP(5,i)
       vtime=0
       vspin=9
       WRITE(nunit3)idup,istup,imothup1,imothup2,icol1,icol2,&
            px,py,pz,p0,pmass,vtime,vspin
    ENDDO
  END SUBROUTINE unwei_writer_Decay

  SUBROUTINE readcuts_Decay
    IMPLICIT NONE
    CHARACTER(len=24)::file
    LOGICAL::lexist
    INTEGER::iounit,flag=0
    INQUIRE(FILE=TRIM(input_dir)//"default.inp",EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file default.inp does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//"default.inp",OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//"default.inp",NUMBER=iounit)
       CLOSE(UNIT=iounit)
    ENDIF
    OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"default.inp")
    ! open user's input file 
    IF(TRIM(Input_File)/="default.inp")THEN
       INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
       IF(.NOT.lexist)THEN
          PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
          STOP
       ENDIF
       INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),OPENED=lexist)
       IF(lexist)THEN
          INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),NUMBER=iounit)
          CLOSE(UNIT=iounit)
       ENDIF
       OPEN(UNIT=uinput,FILE=TRIM(input_dir)//TRIM(Input_File))
    ELSE
       flag=1
    ENDIF
    ! lepton cutoffs
    ho_dptcl=0d0
    ho_dycl=1d9
    IF(absrap)THEN
       ho_dycllow=0d0
    ELSE
       ho_dycllow=-1d9
    ENDIF
    ho_detacl=1d9
    ho_decl=0d0
    ho_dcl=1d0
    ho_dptcl=readvalue_r("decay_minptl",flag)
    ho_dycl=readvalue_r("decay_maxycl",flag)
    ho_dycllow=readvalue_r("decay_minycl",flag)
    ho_detacl=readvalue_r("decay_maxrapl",flag)
    ho_decl=readvalue_r("decay_minel",flag)
    ho_dcl=readvalue_r("decay_maxcl",flag)
    CLOSE(UNIT=udefault)
    CLOSE(UNIT=uinput)
    WRITE(*,*)'---------------------------------------------------'
    WRITE(*,*)'        the decay cuts        '
    WRITE(*,*)"pt  of  lepton   ",ho_dptcl
    WRITE(*,*)"energy of lepton  ",ho_decl
    WRITE(*,*)"max y rapidity of lepton ",ho_dycl
    WRITE(*,*)"min y rapidity of lepton ",ho_dycllow
    WRITE(*,*)"max rapidity of lepton ",ho_detacl
    WRITE(*,*)"cos-beam of lepton ",ho_dcl
    WRITE(*,*)'---------------------------------------------------'
    IF(literature_cutoffs.NE.0)THEN
       WRITE(*,*)"impose the special cutoffs in literature arXiv:",literature_cutoffs
    ENDIF
    ! read the decay params
    INQUIRE(FILE=TRIM(input_dir)//"decay_param_default.inp",EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file decay_param_default_default.inp does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//"decay_param_default.inp",OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//"decay_param_default.inp",NUMBER=iounit)
       CLOSE(UNIT=iounit)
    ENDIF
    OPEN(UNIT=udefault,FILE=TRIM(input_dir)//"decay_param_default.inp")
    ! open user's decay_param_user.inp file
    INQUIRE(FILE=TRIM(input_dir)//TRIM(Input_File),EXIST=lexist)
    IF(.NOT.lexist)THEN
       PRINT *,"Warning: the file "//TRIM(Input_File)//" does not exist ! STOP !"
       STOP
    ENDIF
    INQUIRE(FILE=TRIM(input_dir)//"decay_param_user.inp",OPENED=lexist)
    IF(lexist)THEN
       INQUIRE(FILE=TRIM(input_dir)//"decay_param_user.inp",NUMBER=iounit)
       CLOSE(UNIT=iounit)
    ENDIF
    OPEN(UNIT=uinput,FILE=TRIM(input_dir)//"decay_param_user.inp")
    ho_dmass_jpsi=readvalue_r("decay_mass_jpsi",flag)
    ho_dmass_chic0=readvalue_r("decay_mass_chic0",flag)
    ho_dmass_chic1=readvalue_r("decay_mass_chic1",flag)
    ho_dmass_chic2=readvalue_r("decay_mass_chic2",flag)
    ho_dmass_hc=readvalue_r("decay_mass_hc",flag)
    ho_dmass_etac=readvalue_r("decay_mass_etac",flag)
    ho_dmass_psi2S=readvalue_r("decay_mass_psi2S",flag)
    ho_dmass_Y1S=readvalue_r("decay_mass_Y1S",flag)
    ho_dmass_chib11P=readvalue_r("decay_mass_chib11P",flag)
    ho_dmass_chib21P=readvalue_r("decay_mass_chib21P",flag)
    ho_dmass_chib01P=readvalue_r("decay_mass_chib01P",flag)
    ho_dmass_hb1P=readvalue_r("decay_mass_hb1P",flag)
    ho_dmass_chib02P=readvalue_r("decay_mass_chib02P",flag)
    ho_dmass_chib12P=readvalue_r("decay_mass_chib12P",flag)
    ho_dmass_chib22P=readvalue_r("decay_mass_chib22P",flag)
    ho_dmass_Y2S=readvalue_r("decay_mass_Y2S",flag)
    ho_dmass_Y3S=readvalue_r("decay_mass_Y3S",flag)
    ho_dmass_chib3P=readvalue_r("decay_mass_chib3P",flag)
    CLOSE(UNIT=udefault)
    CLOSE(UNIT=uinput)
    RETURN
  END SUBROUTINE readcuts_Decay


  SUBROUTINE cuts_Decay(icut)
    ! only cut on the decay particles
    IMPLICIT NONE
    INTEGER,INTENT(OUT)::icut
    INTEGER::i,kk
    REAL(KIND(1d0))::eta,pt,c
    REAL(KIND(1d0))::ptminmuon
    iNTEGER::flag
    INTEGER::muon4psiATLAS
    INCLUDE "stdhep.inc"
    icut=0
    flag=0
    muon4psiATLAS=0
    DO i=1,NHEP
       IF(ISHEP(i).EQ.-1)CYCLE ! only cut on the final state
       IF(ISTHEP(i).EQ.1)CYCLE ! exclude the mother particles, which have been cutted before
       IF(JMOHEP(1,i).LE.0.OR.JMOHEP(1,i).GT.NHEP)CYCLE
       IF(ISTHEP(JMOHEP(1,i)).NE.1.OR.ISHEP(JMOHEP(1,i)).EQ.-1)CYCLE
       ! cut is only applied on the decay particles
       IF(lepton_pdg(IDHEP(i)))THEN
          ! cut on the leptons
          pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
          IF(pt.LT.ho_dptcl)THEN
             flag=1
             EXIT
          ENDIF
          eta=prapidity(PHEP(1:4,i))
          IF(ABS(eta).GT.ho_detacl)THEN
             flag=2
             EXIT
          ENDIF
          eta=rapidity(PHEP(1:4,i))
          IF(absrap)eta=ABS(eta)
          IF(eta.GT.ho_dycl.OR.eta.LT.ho_dycllow)THEN
             flag=3
             EXIT
          ENDIF
          IF(PHEP(4,i).LT.ho_decl)THEN
             flag=4
             EXIT
          ENDIF
          c=PHEP(3,i)/PHEP(4,i)
          IF(ABS(c).GT.ho_dcl)THEN
             flag=5
             EXIT
          ENDIF
          ! special cutoffs
          IF(literature_cutoffs.EQ.14062380)THEN
             ! DZero (arXiv:1406.2380) special cutoffs
             eta=prapidity(PHEP(1:4,i))
             IF(ABS(eta).LT.1.35d0)THEN
                pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                IF(pt.LE.2d0)THEN
                   flag=6
                   EXIT
                ENDIF
             ELSEIF(ABS(eta).LT.2d0.AND.ABS(eta).GT.1.35d0)THEN
                IF(PHEP(4,i).LE.4d0)THEN
                   flag=7
                   EXIT
                ENDIF
             ELSE
                flag=8
                EXIT
             ENDIF
          ELSEIF(literature_cutoffs.EQ.14060000)THEN
             ! ATLAS condition, speical cutoffs
             ! there is no arXiv number now
             ! At least one of the Jpsi must have two muons with pT>4 GeV each
             ! muon pT>2.5 GeV and muon |eta|<2.3
             eta=prapidity(PHEP(1:4,i))
             IF(ABS(eta).GE.2.3d0)THEN
                flag=9
                EXIT
             ENDIF
             pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
             IF(pt.LE.2.5d0)THEN
                flag=10
                EXIT
             ENDIF
             IF(IDHEP(JMOHEP(1,i)).EQ.443.AND.muon4psiATLAS.LT.2)THEN
                IF(pt.GT.4d0)THEN
                   muon4psiATLAS=muon4psiATLAS+1
                ELSE
                   muon4psiATLAS=muon4psiATLAS-1
                ENDIF
             ENDIF
             IF(muon4psiATLAS.EQ.-2)muon4psiATLAS=0
          ELSEIF(literature_cutoffs.EQ.14070000)THEN
             ! CMS condition (from Junquan Tao for psi+psi)
             eta=prapidity(PHEP(1:4,i))
             IF(ABS(eta).LT.1.2d0)THEN
                pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                IF(pt.LE.3.5d0)THEN
                   flag=12
                   EXIT
                ENDIF
             ELSEIF(ABS(eta).GT.1.2d0.AND.ABS(eta).LT.1.6d0)THEN
                ptminmuon=3.5d0-(ABS(eta)-1.2d0)*1.5d0/0.4d0
                pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                IF(pt.LE.ptminmuon)THEN
                   flag=13
                   EXIT
                ENDIF
             ELSEIF(ABS(eta).GT.1.6d0.AND.ABS(eta).LT.2.4d0)THEN
                pt=SQRT(PHEP(1,i)**2+PHEP(2,i)**2)
                IF(pt.LE.2.0d0)THEN
                   flag=14
                   EXIT
                ENDIF
             ELSE
                flag=15
                EXIT
             ENDIF
          ENDIF
       ENDIF
    ENDDO
    IF(flag.EQ.0.AND.literature_cutoffs.EQ.14060000)THEN
       ! ATLAS condition, speical cutoffs
       ! there is no arXiv number now
       ! At least one of the Jpsi must have two muons with pT>4 GeV each
       IF(muon4psiATLAS.LT.2)THEN
          flag=11
       ENDIF
    ENDIF
    IF(flag.EQ.0)THEN
       icut=1
    ENDIF
    RETURN
  END SUBROUTINE cuts_Decay
END MODULE Decay_interface
