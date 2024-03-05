#ifndef __POSTPROCESS_H__
#define __POSTPROCESS_H__

/* Include Headers Here */
#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <vector>

/* Extern Functions Here */
extern void sph_save_single(SPH *);
extern void sph_save_last(SPH *);
extern void sph_save_rigid(SPH *);
#endif