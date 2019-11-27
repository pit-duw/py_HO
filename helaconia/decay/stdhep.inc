!*******************************************************
!  ---------------------------------------------------------------------
!     The stdhep event common
!  ---------------------------------------------------------------------
! This is the "official" event common of the HEP community for MC events.
! It is the same with EvtGen
!  NEVHEP      = event number
!  NHEP        = number of particles ("lines")
!  ISTHEP(N)   = stability code; 0 if stable, 1 if decays.
!  IDHEP(I)    = particle ID code; see PDG '92, III.73
!  JMOHEP(1,I) = line number of mother particle
!  JMOHEP(2,I) = line number of last mother particle, if formed by string
!  JDAHEP(1,I) = line number of 1st daughter particle
!  JDAHEP(2,I) = line number of last daughter particle - MUST BE CONSECUTIVE!
!  ISHEP(I)    = -1 if initial state, 1 if final state
!  PHEP(5,I)   = px,py,pz,E,m of particle
!  VHEP(4,I)   = x,y,z,L of production vertex
!  IPOLHEP(I)  = polarization of particle i, 99 means polarized summed
!  ICOLHEP(1,I)   = color of particle
!  ICOLHEP(2,I)   = color of particle

INTEGER,PARAMETER::nmxhep=4000  ! maximum number of particles
REAL(KIND(1d0)),DIMENSION(5,nmxhep)::phep ! four-momentum, mass [GeV]
REAL(KIND(1d0)),DIMENSION(4,nmxhep)::vhep ! vertex [mm]
INTEGER::nevhep ! serial number
INTEGER::nhep ! number of particles
INTEGER,DIMENSION(nmxhep)::isthep ! status code
INTEGER,DIMENSION(nmxhep)::idhep ! particle ident KF
INTEGER,DIMENSION(2,nmxhep)::jmohep ! parent particles
INTEGER,DIMENSION(2,nmxhep)::jdahep ! children particles
INTEGER,DIMENSION(nmxhep)::ISHEP  ! initial or final states
INTEGER,DIMENSION(nmxhep)::IPOLHEP ! polarization, if 99 means polarized summed
INTEGER,DIMENSION(2,nmxhep)::ICOLHEP ! color
COMMON/hepevt/nevhep,nhep,isthep,idhep,jmohep,jdahep,ishep,phep,vhep,ipolhep,ICOLHEP

!$$$      INTEGER NMXHEP,NEVHEP,NHEP,ISTHEP,IDHEP,JMOHEP,JDAHEP
!$$$      REAL    PHEP,VHEP
!$$$      PARAMETER (NMXHEP=2000)
!$$$      COMMON/STDHEP/NEVHEP,NHEP,ISTHEP(NMXHEP),IDHEP(NMXHEP),
!$$$     *JMOHEP(2,NMXHEP),JDAHEP(2,NMXHEP),PHEP(5,NMXHEP),VHEP(4,NMXHEP)
!$$$      SAVE /STDHEP/
