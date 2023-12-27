import vtk

fd = open("./wedge_1.vtk","w+")

reader = vtk.vtkUnstructuredGridReader()
reader.SetFileName("./wedge.vtk")
reader.ReadAllScalarsOn()
reader.ReadAllVectorsOn()
reader.Update()

vtkdata = reader.GetOutput()

n=0
for i in range(vtkdata.GetNumberOfPoints()):
    if vtkdata.GetPoint(i)[2] == 0 and (vtkdata.GetPoint(i)[0] != 0 or vtkdata.GetPoint(i)[1]!=0):
        n = n+1

fd.write("# vtk DataFile Version 5.1\nvtk output\nASCII\nDATASET UNSTRUCTURED_GRID\nPOINTS {} float\n".format(n))

print(vtkdata.GetNumberOfPoints())
for i in range(vtkdata.GetNumberOfPoints()):
    if vtkdata.GetPoint(i)[2] == 0 and (vtkdata.GetPoint(i)[0] != 0 or vtkdata.GetPoint(i)[1]!=0):
        '''
        fd.write("{.4f}".format(vtkdata.GetPoint(i)[0])+" {.4f}".format(vtkdata.GetPoint(i)[1])+ \
                " {.4f}".format(vtkdata.GetPoint(i)[2]))
        '''
        fd.write(str(vtkdata.GetPoint(i)[0]/1000)+" "+str(vtkdata.GetPoint(i)[1]/1000)+ \
                " "+str(vtkdata.GetPoint(i)[2]/1000)+"\n")
        print(vtkdata.GetPoint(i))
    
    #print(vtkdata.GetPoint(i))
fd.close()
