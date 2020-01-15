subroutine ALLOC_VARS
  ! allocate and Initialize Most Arrays
  use ALL_VARS
  implicit none
  INTEGER :: NCT

  NCT = N*3

  ! Grid Metrics
  allocate(XC(0:N))            ;XC   = zero   !!X-COORD AT FACE CENTER 
  allocate(YC(0:N))            ;YC   = zero   !!Y-COORD AT FACE CENTER
  allocate(VX(0:M))            ;VX   = zero   !!X-COORD AT GRID POINT
  allocate(VY(0:M))            ;VY   = zero   !!Y-COORD AT GRID POINT

  ! Node, Boundary Condition, and Control Volume
  allocate(NV(0:N,4))           ;NV       = 0  !!NODE NUMBERING FOR ELEMENTS
  allocate(NBE(0:N,3))          ;NBE      = 0  !!INDICES OF ELEMENT NEIGHBORS
  allocate(NTVE(0:M))           ;NTVE     = 0 
  allocate(ISONB(0:M))          ;ISONB    = 0  !!NODE MARKER = 0,1,2
  allocate(ISBCE(0:N))          ;ISBCE    = 0 

  ! 1-d arrays for the sigma coordinate
  allocate(Z(KB))               ; Z      = zero    !!SIGMA COORDINATE VALUE 
  allocate(ZZ(KB))              ; ZZ     = zero    !!INTRA LEVEL SIGMA VALUE
  allocate(DZ(KB))              ; DZ     = zero    !!DELTA-SIGMA VALUE
  allocate(DZZ(KB))             ; DZZ    = zero    !!DELTA OF INTRA LEVEL SIGMA 

  ! 2-d flow variable arrays at nodes
  allocate(H(0:M))       ;H    = zero       !!BATHYMETRIC DEPTH   
  allocate(D(0:M))       ;D    = zero       !!DEPTH   
  allocate(EL(0:M))      ;EL   = zero       !!SURFACE ELEVATION
  allocate(ET(0:M))      ;ET  = zero       !!SURFACE ELEVATION PREVIOUS TIMESTEP

  ! internal mode arrays-(element based)
  allocate(U(0:N, KB))       ;U     = zero   !!X-VELOCITY
  allocate(V(0:N, KB))       ;V     = zero   !!Y-VELOCITY
  allocate(W(0:N, KB))       ;W     = zero   !!VERTICAL VELOCITY IN SIGMA SYSTEM
  allocate(WW(0:N, KB))      ;WW    = zero   !!Z-VELOCITY
  allocate(UT(0:N, KB))      ;UT    = zero   !!X-VELOCITY FROM PREVIOUS TIMESTEP
  allocate(VT(0:N, KB))      ;VT    = zero   !!Y-VELOCITY FROM PREVIOUS TIMESTEP
  allocate(WT(0:N, KB))      ;WT    = zero   !!VERTICAL VELOCITY FROM PREVIOUS TIMESTEP
  allocate(WWT(0:N, KB))     ;WWT   = zero   !!Z-VELOCITY FROM PREVIOUS TIMESTEP
  allocate(KH(0:N, KB))     ;KH    = zero   !!TURBULENT QUANTITY

  ! 3d variable arrays-(node based)
  allocate(T1(0:M, KB))       ;T1     = zero  !!TEMPERATURE AT NODES
  allocate(S1(0:M, KB))       ;S1     = zero  !!SALINITY AT NODES
  allocate(R1(0:M, KB))       ;R1   = zero  !!DENSITY AT NODES
  allocate(TT1(0:M, KB))      ;TT1    = zero  !!TEMPERATURE FROM PREVIOUS TIME
  allocate(ST1(0:M, KB))      ;ST1    = zero  !!SALINITY FROM PREVIOUS TIME
  allocate(RT1(0:M, KB))      ;RT1 = zero
  allocate(WTS(0:M, KB))      ;WTS    = zero  !!VERTICAL VELOCITY IN SIGMA SYSTEM

  ! Shape coefficient arrays and control volume metrics
  allocate(A1U(0:N, 4))         ;A1U   = zero
  allocate(A2U(0:N, 4))         ;A2U   = zero
  allocate(AWX(0:N, 3))         ;AWX   = zero
  allocate(AWY(0:N, 3))         ;AWY   = zero
  allocate(AW0(0:N, 3))         ;AW0   = zero

  return
end subroutine
