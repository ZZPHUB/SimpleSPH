#ifndef __LIB_H__
#define __LIB_H__

/* Include Headers Here*/
#include "SPH.H"
#include "sys/resource.h"
#include "cmath"
#include <iostream>

/* Extern Function Here*/
extern void get_max_men(void);
extern void set_max_men(void);
extern void nnps_direct(SPH_PARTICLE *,SPH_PAIR *);
extern void ptc_kernel(SPH_PARTICLE *,SPH_PAIR *,SPH_KERNEL *);

/* Define Symbols Here*/
#define PI 3.14159265358
#define ALPHA 15.0/(7*PI*pow(PTC_SML,2))


#endif