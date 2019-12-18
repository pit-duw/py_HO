MODULE MLM_merge
  IMPLICIT NONE
CONTAINS
  SUBROUTINE MLM_matching(p_shower,n_shower,matching_scale,n_match,&
       max_multiplicity, matched)
    IMPLICIT NONE
    INCLUDE "HEPEUPf90.inc" ! hard events information
    INTEGER::N_HARD,i,j
    REAL(KIND(1d0)),DIMENSION(4,MAXNUP)::P_HARD
    INTEGER::Njet_hard
    INTEGER,DIMENSION(MAXNUP)::jet_hard
    INTEGER,INTENT(IN)::n_match
    INTEGER,INTENT(IN)::n_shower
    INTEGER::njet_shower
    INTEGER,DIMENSION(n_shower)::jet_shower
    REAL(KIND(1d0))::palg,rfj,sycut,etacut
    REAL(KIND(1d0)),DIMENSION(4,MAXNUP)::pjet_hard
    REAL(KIND(1d0)),INTENT(IN)::matching_scale !,matching_eta
    REAL(KIND(1d0)),DIMENSION(4,n_shower),INTENT(IN)::p_shower
    REAL(KIND(1d0)),DIMENSION(4,n_shower)::pjet_shower
    LOGICAL,INTENT(OUT)::matched
    LOGICAL,INTENT(IN)::max_multiplicity
    LOGICAL::jet_shower_matched
    matched=.TRUE.
    palg=1d0
    sycut=matching_scale
    !etacut=matching_eta
    rfj=1d0
    CALL fastjetppgenkt_etamax(p_shower,n_shower,rfj,sycut,etacut,palg,pjet_shower,njet_shower,jet_shower)
    ! First check that the number of jets after showering is okay if the number
    ! of jets required in the matching
    IF(max_multiplicity)THEN
       IF(njet_shower.LT.n_match)THEN
          matched=.FALSE.
          RETURN
       ENDIF
    ELSE
       IF(njet_shower.NE.n_match)THEN
          matched=.FALSE.
          RETURN
       ENDIF
    ENDIF
    n_hard=0
    DO i=1,NUP
       IF(ISTUP(i).EQ.1.AND.(ABS(IDUP(i)).LT.6.OR.IDUP(i).EQ.21))THEN
          !it is the final QCD parton
          n_hard=n_hard+1
          DO j=1,4
             p_hard(j,n_hard)=PUP(j,i)
          ENDDO
       ENDIF
    ENDDO
    sycut=0d0
    CALL fastjetppgenkt_etamax(p_hard,n_hard,rfj,sycut,etacut,palg,pjet_hard,njet_hard,jet_hard)
    ! check that the hardest jet after showering match with the ardest jets before showering
    DO i=1,n_match
       jet_shower_matched=.FALSE.
       DO j=1,njet_hard
          IF(getdrv(pjet_shower(1,i),pjet_hard(1,j)).LT.1.5d0*rfj)THEN
             jet_shower_matched=.TRUE.
             EXIT
          ENDIF
       ENDDO
       IF(.NOT.jet_shower_matched)THEN
          matched=.FALSE.
          RETURN
          EXIT
       ENDIF
    ENDDO
    RETURN
  END SUBROUTINE MLM_matching

  FUNCTION getdrv(p1,p2)
    IMPLICIT NONE
    REAL(KIND(1d0)),DIMENSION(5),INTENT(IN)::p1,p2
    REAL(KIND(1d0))::getdrv
    getdrv=getdr(p1(4),p1(1),p1(2),p1(3),&
         p2(4),p2(1),p2(2),p2(3))
    RETURN
  END FUNCTION getdrv

  FUNCTION getdr(en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2)
    IMPLICIT NONE
    REAL(KIND(1d0))::getdr
    REAL(KIND(1d0)),INTENT(IN)::en1,ptx1,pty1,pl1,en2,ptx2,pty2,pl2
    REAL(KIND(1d0))::deta,dphi
    deta=getpseudorap(en1,ptx1,pty1,pl1)-getpseudorap(en2,ptx2,pty2,pl2)
    dphi=getdelphi(ptx1,pty1,ptx2,pty2)
    getdr=DSQRT(dphi**2+deta**2)
    RETURN
  END FUNCTION getdr

  FUNCTION getdelphi(ptx1,pty1,ptx2,pty2)
    IMPLICIT NONE
    REAL(KIND(1d0))::getdelphi
    REAL(KIND(1d0)),INTENT(IN)::ptx1,pty1,ptx2,pty2
    REAL(KIND(1d0)),PARAMETER::tiny=1d-5
    REAL(KIND(1d0))::pt1,pt2,tmp
    pt1=DSQRT(ptx1**2+pty1**2)
    pt2=DSQRT(ptx2**2+pty2**2)
    IF(pt1.NE.0d0.AND.pt2.NE.0d0)THEN
       tmp=ptx1*ptx2+pty1*pty2
       tmp=tmp/(pt1*pt2)
       IF(ABS(tmp).GT.1d0+tiny)THEN
          WRITE(*,*)'ERROR:Cos(phi) is larger than 1'
          STOP
       ELSEIF(ABS(tmp).GE.1d0)THEN
          tmp=SIGN(1d0,tmp)
       ENDIF
       tmp=DACOS(tmp)
    ELSE
       tmp=1d8
    ENDIF
    getdelphi=tmp
    RETURN
  END FUNCTION getdelphi

  FUNCTION getpseudorap(en,ptx,pty,pl)
    IMPLICIT NONE
    REAL(KIND(1d0))::getpseudorap
    REAL(KIND(1d0)),INTENT(IN)::en,ptx,pty,pl
    REAL(KIND(1d0)),PARAMETER::tiny=1d-5
    REAL(KIND(1d0))::pt,eta,th
    pt=DSQRT(ptx**2+pty**2)
    IF(pt.LT.tiny.AND.ABS(pl).LT.tiny)THEN
       eta=SIGN(1d0,pl)*1d8
    ELSE
       th=ATAN2(pt,pl)
       eta=-LOG(TAN(th/2d0))
    ENDIF
    getpseudorap=eta
    RETURN
  END FUNCTION getpseudorap

  FUNCTION getrapidity(en,pl)
    IMPLICIT NONE
    REAL(KIND(1d0))::getrapidity
    REAL(KIND(1d0)),INTENT(IN)::en,pl
    REAL(KIND(1d0)),PARAMETER::tiny=1d-8
    REAL(KIND(1d0))::xplus,xminus,y
    xplus=en+pl
    xminus=en-pl
    IF(xplus.GT.tiny.AND.xminus.GT.tiny)THEN
       IF((xplus/xminus).GT.tiny.AND.(xminus/xplus).GT.tiny)THEN
          y=0.5d0*LOG(xplus/xminus)
       ELSE
          y=SIGN(1d0,pl)*1d8
       ENDIF
    ELSE
       y=SIGN(1d0,pl)*1d8
    ENDIF
    getrapidity=y
    RETURN
  END FUNCTION getrapidity

  FUNCTION PTCALC(P)
    IMPLICIT NONE
    REAL(KIND(1d0))::PTCALC
    REAL(KIND(1d0)),DIMENSION(4),INTENT(IN)::P
    REAL(KIND(1d0))::PTSQ
    PTSQ=P(1)**2+P(2)**2
    IF(PTSQ.EQ.0d0)THEN
       PTCALC=0d0
    ELSE
       PTCALC=DSQRT(PTSQ)
    ENDIF
    RETURN
  END FUNCTION PTCALC
END MODULE MLM_merge
