    Pl1(0)=mass/2d0
    Pl1(1)=mass/2d0*sinth*DCOS(phi)
    Pl1(2)=mass/2d0*sinth*DSIN(phi)
    Pl1(3)=mass/2d0*costh
    Pl2(0)=mass/2d0
    Pl2(1:3)=-Pl1(1:3)
    PBOO(1:3)=PV(1:3)
    PBOO(4)=PV(0)
    PL(1:3)=Pl1(1:3)
    PL(4)=Pl1(0)
    ! Rotate and boost Pl1 
    CALL ROTATEL(squ,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pl1(0)=PL(4)
    Pl1(1:3)=PL(1:3)
    ! Rotate and boost Pl2
    PL(1:3)=Pl2(1:3)
    PL(4)=Pl2(0)
    CALL ROTATEL(squ,PL)
    CALL BOOSTL(mass,PBOO,PL)
    Pl2(0)=PL(4)
    Pl2(1:3)=PL(1:3)