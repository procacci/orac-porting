module sasamod
  use parst
  implicit none
  !-------------------- sasa stuff --------------------------------------------
  character*7  sasa_betb(400)                  ! max 400 SASA types
  INTEGER sasa_nbtype(m1)
  REAL*8 sasa(m1),sasaener,r2cut,rsolv 
  real*8       rvdw(400),psasa(400),ssasa(400)  ! max 400 SASA types
  !-------------------- end of sasa stuff --------------------------------------
end module sasamod
