#ifndef __EQUATIONS_H__
#define __EQUATIONS_H__

/* Include Headers Here */
#include "SPH.H"

/* Extern Functions Here */
extern void ptc_fluid_pressure(SPH *);
extern void ptc_dif_density(SPH *);
extern void ptc_acc(SPH *);
extern void ptc_viscous(SPH *);
extern void sph_governing(SPH *);

#endif