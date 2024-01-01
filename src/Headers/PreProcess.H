#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__


/* Headers Include Here*/
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include "SPH.H"

/* Extern Function Here*/
extern void solid_ptc_generate(SPH_PARTICLE *);
extern unsigned int solid_ptc_num(void);
extern void ptc_generate(SPH_PARTICLE *);
extern void ptc_init(SPH_PARTICLE *,RIGID *,RIGID *);
#endif