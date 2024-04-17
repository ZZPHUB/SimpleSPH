#ifndef __IO_H__
#define __IO_H__

/* Headers Include Here*/
#include "SPH.cuh"

/* Extern Function Here*/
extern void ptc_rigid_generate(SPH *);
extern unsigned int ptc_rigid_num(void);
extern void ptc_rigid_init(SPH *);
extern void ptc_read_vtk(SPH *);
extern void ptc_generate(SPH *);
extern void ptc_info_init(SPH *);
extern void ptc_init(SPH *);

extern void sph_save_single(SPH *);
extern void sph_save_last(SPH *);

#endif