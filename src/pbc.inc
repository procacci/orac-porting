#include "config.H"

#ifdef MAGIC
      REAL*8 PBC,x,magic
      DATA magic/6755399441055744.D0/
      PBC(x)=((0.5D0*x+magic)-magic)
#endif
#ifdef ANINT 
      REAL*8 PBC,x
      PBC(x)=DNINT(0.5D0*x)
#endif
