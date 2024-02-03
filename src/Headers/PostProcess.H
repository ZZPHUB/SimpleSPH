#ifndef __POSTPROCESS_H__
#define __POSTPROCESS_H__

/* Include Headers Here */
#include "SPH.H"
#include <fstream>
#include <iomanip>
#include <stdlib.h>

/* Extern Functions Here */
extern void ptc_vtk_direct(SPH *,char *);
extern void ptc_pressure_conv(SPH *);
#endif