#ifndef __POSTPROCESS_H__
#define __POSTPROCESS_H__

/* Include Headers Here */
#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>

/* Extern Functions Here */
extern void sph_save_single(SPH *,char *);
extern void sph_save_last(SPH *);
#endif