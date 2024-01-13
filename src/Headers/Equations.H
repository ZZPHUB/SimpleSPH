#ifndef __EQUATIONS_H__
#define __EQUATIONS_H__

/* Include Headers Here */
#include "SPH.H"

/* Extern Functions Here */
extern void fluid_ptc_pressure(SPH *);
extern void ptc_dif_density(SPH *);
extern void ptc_acc(SPH *);
extern void ptc_viscous(SPH *);

#endif