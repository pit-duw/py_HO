      ! this is stored for hard events
      ! i.e. those from Les Houches Event file
      INTEGER,PARAMETER::MAXNUP=500
      INTEGER::NUP,IDPRUP
      INTEGER,DIMENSION(MAXNUP)::IDUP,ISTUP
      INTEGER,DIMENSION(2,MAXNUP)::MOTHUP,ICOLUP
      REAL(KIND(1d0))::XWGTUP,SCALUP,AQEDUP,AQCDUP
      REAL(KIND(1d0)),DIMENSION(5,MAXNUP)::PUP
      REAL(KIND(1d0)),DIMENSION(MAXNUP)::VTIMUP,SPINUP
      COMMON/HEPEUP/NUP,IDPRUP,XWGTUP,SCALUP,AQEDUP,AQCDUP,&
	             IDUP,ISTUP,MOTHUP,ICOLUP,PUP,&
		     VTIMUP,SPINUP
     