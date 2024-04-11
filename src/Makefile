VPATH = ./Equation:./Lib:./Post_Processing:./Pre_Processing:./SPH:
# 正则表达式表示目录下所有.c文件
SRC = $(notdir $(wildcard ./*/*.cu))
# 头文件路径
INCLUDE = -I./Headers
# 生成文SRC
BUILD_DIR = ./build
OBJ = $(addprefix $(BUILD_DIR)/,$(patsubst %.cu,%.o,$(SRC))) 
# vtk链接库路径以及需要链接的库
LD = -L /usr/lib/x86_64-linux-gnu -lvtkChartsCore-9.1 -lvtkCommonColor-9.1 -lvtkCommonComputationalGeometry-9.1 -lvtkCommonCore-9.1 -lvtkCommonDataModel-9.1 -lvtkCommonExecutionModel-9.1 -lvtkCommonMath-9.1 -lvtkCommonMisc-9.1 -lvtkCommonPython-9.1 -lvtkCommonSystem-9.1 -lvtkCommonTransforms-9.1 -lvtkDICOMParser-9.1 -lvtkDomainsChemistry-9.1 -lvtkDomainsChemistryOpenGL2-9.1 -lvtkDomainsParallelChemistry-9.1 -lvtkFiltersAMR-9.1 -lvtkFiltersCore-9.1 -lvtkFiltersExtraction-9.1 -lvtkFiltersFlowPaths-9.1 -lvtkFiltersGeneral-9.1 -lvtkFiltersGeneric-9.1 -lvtkFiltersGeometry-9.1 -lvtkFiltersHybrid-9.1 -lvtkFiltersHyperTree-9.1 -lvtkFiltersImaging-9.1 -lvtkFiltersModeling-9.1 -lvtkFiltersParallel-9.1 -lvtkFiltersParallelGeometry-9.1 -lvtkFiltersParallelImaging-9.1 -lvtkFiltersParallelMPI-9.1 -lvtkFiltersParallelVerdict-9.1 -lvtkFiltersPoints-9.1 -lvtkFiltersProgrammable-9.1 -lvtkFiltersPython-9.1 -lvtkFiltersSMP-9.1 -lvtkFiltersSelection-9.1 -lvtkFiltersSources-9.1 -lvtkFiltersStatistics-9.1 -lvtkFiltersTexture-9.1 -lvtkFiltersTopology-9.1 -lvtkFiltersVerdict-9.1 -lvtkGeovisCore-9.1 -lvtkIOAMR-9.1 -lvtkIOAsynchronous-9.1 -lvtkIOCGNSReader-9.1 -lvtkIOCONVERGECFD-9.1 -lvtkIOChemistry-9.1 -lvtkIOCityGML-9.1 -lvtkIOCore-9.1 -lvtkIOEnSight-9.1 -lvtkIOExodus-9.1 -lvtkIOExport-9.1 -lvtkIOExportGL2PS-9.1 -lvtkIOExportPDF-9.1 -lvtkIOGeometry-9.1 -lvtkIOHDF-9.1 -lvtkIOIOSS-9.1 -lvtkIOImage-9.1 -lvtkIOImport-9.1 -lvtkIOInfovis-9.1 -lvtkIOLSDyna-9.1 -lvtkIOLegacy-9.1 -lvtkIOMINC-9.1 -lvtkIOMPIImage-9.1 -lvtkIOMPIParallel-9.1 -lvtkIOMotionFX-9.1 -lvtkIOMovie-9.1 -lvtkIONetCDF-9.1 -lvtkIOOggTheora-9.1 -lvtkIOPLY-9.1 -lvtkIOParallel-9.1 -lvtkIOParallelNetCDF-9.1 -lvtkIOParallelXML-9.1 -lvtkIOSQL-9.1 -lvtkIOSegY-9.1 -lvtkIOTecplotTable-9.1 -lvtkIOVeraOut-9.1 -lvtkIOVideo-9.1 -lvtkIOXML-9.1 -lvtkIOXMLParser-9.1 -lvtkImagingColor-9.1 -lvtkImagingCore-9.1 -lvtkImagingFourier-9.1 -lvtkImagingGeneral-9.1 -lvtkImagingHybrid-9.1 -lvtkImagingMath-9.1 -lvtkImagingMorphological-9.1 -lvtkImagingSources-9.1 -lvtkImagingStatistics-9.1 -lvtkImagingStencil-9.1 -lvtkInfovisCore-9.1 -lvtkInfovisLayout-9.1 -lvtkInteractionImage-9.1 -lvtkInteractionStyle-9.1 -lvtkInteractionWidgets-9.1 -lvtkJava-9.1 -lvtkParallelCore-9.1 -lvtkParallelDIY-9.1 -lvtkParallelMPI-9.1 -lvtkParallelMPI4Py-9.1 -lvtkPythonContext2D-9.1 -lvtkPythonInterpreter-9.1 -lvtkRenderingAnnotation-9.1 -lvtkRenderingContext2D-9.1 -lvtkRenderingContextOpenGL2-9.1 -lvtkRenderingCore-9.1 -lvtkRenderingFreeType-9.1 -lvtkRenderingGL2PSOpenGL2-9.1 -lvtkRenderingImage-9.1 -lvtkRenderingLOD-9.1 -lvtkRenderingLabel-9.1 -lvtkRenderingOpenGL2-9.1 -lvtkRenderingSceneGraph-9.1 -lvtkRenderingTk-9.1 -lvtkRenderingUI-9.1 -lvtkRenderingVolume-9.1 -lvtkRenderingVolumeOpenGL2-9.1 -lvtkRenderingVtkJS-9.1 -lvtkTestingRendering-9.1 -lvtkViewsContext2D-9.1 -lvtkViewsCore-9.1 -lvtkViewsInfovis-9.1 -lvtkWrappingPythonCore3.10-9.1 -lvtkWrappingTools-9.1 -lvtkcgns-9.1 -lvtkexodusII-9.1 -lvtkfmt-9.1 -lvtkioss-9.1 -lvtkkissfft-9.1 -lvtkloguru-9.1 -lvtkmetaio-9.1 -lvtkpugixml-9.1 -lvtksys-9.1 -lvtkverdict-9.1

ifeq ($(ver), debug)
CPPFLAGS = -Xcompiler -fopenmp -dc -g -Ddebug -G -arch=compute_75 -code=sm_75 
else
CPPFLAGS = -Xcompiler -fopenmp -dc -arch=compute_75 -code=sm_75
endif

TARGET = run

CC = nvcc

.PHONY:all clean

all: $(TARGET)

# 第一行依赖关系：冒号后面为依赖的文件，相当于Hello: main.o a.o b.o
# 第二行规则：$@表示目标文件，$^表示所有依赖文件，$<表示第一个依赖文件
$(TARGET): $(OBJ)
	$(CC) $(CPPFLAGS) $^ -o $@  $(LD)

# 上一句目标文件依赖一大堆.o文件，这句表示所有.o都由相应名字的.cpp文件自动生成
$(BUILD_DIR)/%.o:%.cu
	$(CC) $(CPPFLAGS) $(INCLUDE) -I /usr/include/vtk-9.1 -o $@ $^ 

# make clean删除所有.o和目标文件
clean:
	rm -f $(OBJ) $(TARGET) ../data/postprocess/csv/* ../data/postprocess/vtk/* 

