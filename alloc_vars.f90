SUBROUTINE ALLOC_VARS
  ! Allocate and Initialize Most Arrays
  USE ALL_VARS
  IMPLICIT NONE
  INTEGER :: NCT

  NCT = N*3

  !  ALLOCATE
  ! Grid Metrics
  ALLOCATE(XC(0:N))            ;XC   = ZERO   !!X-COORD AT FACE CENTER 
  ALLOCATE(YC(0:N))            ;YC   = ZERO   !!Y-COORD AT FACE CENTER
  ALLOCATE(VX(0:M))            ;VX   = ZERO   !!X-COORD AT GRID POINT
  ALLOCATE(VY(0:M))            ;VY   = ZERO   !!Y-COORD AT GRID POINT

  ! Node, Boundary Condition, and Control Volume
  ALLOCATE(NV(0:N,4))           ;NV       = 0  !!NODE NUMBERING FOR ELEMENTS
  ALLOCATE(NBE(0:N,3))          ;NBE      = 0  !!INDICES OF ELEMENT NEIGHBORS
  ALLOCATE(NTVE(0:M))           ;NTVE     = 0 
  ALLOCATE(ISONB(0:M))          ;ISONB    = 0  !!NODE MARKER = 0,1,2
  ALLOCATE(ISBCE(0:N))          ;ISBCE    = 0 

  ! 1-d arrays for the sigma coordinate
  ALLOCATE(Z(KB))               ; Z      = ZERO    !!SIGMA COORDINATE VALUE 
  ALLOCATE(ZZ(KB))              ; ZZ     = ZERO    !!INTRA LEVEL SIGMA VALUE
  ALLOCATE(DZ(KB))              ; DZ     = ZERO    !!DELTA-SIGMA VALUE
  ALLOCATE(DZZ(KB))             ; DZZ    = ZERO    !!DELTA OF INTRA LEVEL SIGMA 

  ! 2-d flow variable arrays at nodes
  ALLOCATE(H(0:M))       ;H    = ZERO       !!BATHYMETRIC DEPTH   
  ALLOCATE(D(0:M))       ;D    = ZERO       !!DEPTH   
  ALLOCATE(EL(0:M))      ;EL   = ZERO       !!SURFACE ELEVATION
  ALLOCATE(ET(0:M))      ;ET  = ZERO       !!SURFACE ELEVATION PREVIOUS TIMESTEP

  ! internal mode arrays-(element based)
  ALLOCATE(U(0:N,KB))       ;U     = ZERO   !!X-VELOCITY
  ALLOCATE(V(0:N,KB))       ;V     = ZERO   !!Y-VELOCITY
  ALLOCATE(W(0:N,KB))       ;W     = ZERO   !!VERTICAL VELOCITY IN SIGMA SYSTEM
  ALLOCATE(WW(0:N,KB))      ;WW    = ZERO   !!Z-VELOCITY
  ALLOCATE(UT(0:N,KB))      ;UT    = ZERO   !!X-VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(VT(0:N,KB))      ;VT    = ZERO   !!Y-VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(WT(0:N,KB))      ;WT    = ZERO   !!VERTICAL VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(WWT(0:N,KB))     ;WWT   = ZERO   !!Z-VELOCITY FROM PREVIOUS TIMESTEP
  ALLOCATE(KH(0:N,KB))     ;KH    = ZERO   !!TURBULENT QUANTITY

  ! 3d variable arrays-(node based)
  ALLOCATE(T1(0:M,KB))       ;T1     = ZERO  !!TEMPERATURE AT NODES
  ALLOCATE(S1(0:M,KB))       ;S1     = ZERO  !!SALINITY AT NODES               
  ALLOCATE(R1(0:M,KB))       ;R1   = ZERO  !!DENSITY AT NODES
  ALLOCATE(TT1(0:M,KB))      ;TT1    = ZERO  !!TEMPERATURE FROM PREVIOUS TIME
  ALLOCATE(ST1(0:M,KB))      ;ST1    = ZERO  !!SALINITY FROM PREVIOUS TIME 
  ALLOCATE(RT1(0:M,KB))      ;RT1 = ZERO
  ALLOCATE(WTS(0:M,KB))      ;WTS    = ZERO  !!VERTICAL VELOCITY IN SIGMA SYSTEM

  ! Shape coefficient arrays and control volume metrics
  ALLOCATE(A1U(0:N,4))         ;A1U   = ZERO
  ALLOCATE(A2U(0:N,4))         ;A2U   = ZERO 
  ALLOCATE(AWX(0:N,3))         ;AWX   = ZERO 
  ALLOCATE(AWY(0:N,3))         ;AWY   = ZERO 
  ALLOCATE(AW0(0:N,3))         ;AW0   = ZERO 

  RETURN
END SUBROUTINE ALLOC_VARS
