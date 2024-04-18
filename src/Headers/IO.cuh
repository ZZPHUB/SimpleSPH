#ifndef __IO_H__
#define __IO_H__

/* Headers Include Here*/
#include "SPH.cuh"
#include <string>
#include <fstream>
#include <iostream>
#include <iomanip>

/* Extern Function Here*/
extern void sph_read_vtk(SPH *);
extern void sph_save_single(SPH *);
extern void sph_save_last(SPH *);
extern void sph_read_info(SPH *);
extern void sph_write_info(SPH *);

#endif