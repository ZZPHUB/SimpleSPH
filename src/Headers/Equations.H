#ifndef __EQUATIONS_H__
#define __EQUATIONS_H__

/* Include Headers Here */
#include "SPH.H"

/* Extern Functions Here */
extern void fluid_ptc_pressure(SPH_PARTICLE *);
extern void ptc_dif_density(SPH_PARTICLE *,SPH_PAIR *,SPH_KERNEL *);
extern void ptc_acc(SPH_PARTICLE *,SPH_PAIR *,SPH_KERNEL *);

#endif