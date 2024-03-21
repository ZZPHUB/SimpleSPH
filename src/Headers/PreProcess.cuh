#ifndef __PREPROCESS_H__
#define __PREPROCESS_H__


/* Headers Include Here*/
#include <vtkUnstructuredGridReader.h>
#include <vtkSmartPointer.h>
#include <vtkType.h>
#include <vtkUnstructuredGrid.h>
#include <vtkPointSet.h>
#include <vtkDataSetReader.h>
#include <vtkPointData.h>
#include <vtkDataArray.h>
#include <vtkDoubleArray.h>

#include "SPH.cuh"

/* Extern Function Here*/
extern void ptc_rigid_generate(SPH *);
extern unsigned int ptc_rigid_num(void);
extern void ptc_rigid_init(SPH *);
extern void ptc_read_vtk(SPH *);
extern void ptc_generate(SPH *);
extern void ptc_info_init(SPH *);
extern void ptc_init(SPH *);
#endif