#!/usr/bin/env python
import sys
# noinspection PyUnresolvedReferences
import vtkmodules.vtkInteractionStyle
# noinspection PyUnresolvedReferences
import vtkmodules.vtkRenderingOpenGL2

import glob
from vtkmodules.vtkCommonColor import vtkNamedColors
from vtkmodules.vtkIOGeometry import vtkSTLReader
from vtkmodules.vtkRenderingCore import (
    vtkActor,
    vtkPolyDataMapper,
    vtkRenderWindow,
    vtkRenderWindowInteractor,
    vtkRenderer
)


def main():

    folder = r"C:\emdlab-win64\tmp\cads3d"   # << change this
    stl_files = glob.glob(folder + "/*.stl")

    colors = vtkNamedColors()

    # Renderer and window
    ren = vtkRenderer()
    renWin = vtkRenderWindow()
    renWin.AddRenderer(ren)
    renWin.SetWindowName("EMDLAB: mesh zones")

    iren = vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    # Load each STL file and add to renderer
    for file_path in stl_files:

        reader = vtkSTLReader()
        reader.SetFileName(file_path)
        reader.Update()

        mapper = vtkPolyDataMapper()
        mapper.SetInputConnection(reader.GetOutputPort())

        actor = vtkActor()
        actor.SetMapper(mapper)

        # Basic appearance
        actor.GetProperty().SetDiffuse(0.8)
        actor.GetProperty().SetSpecular(0.3)
        actor.GetProperty().SetSpecularPower(60.0)

        # Optional: different color for each STL
        actor.GetProperty().SetColor(colors.GetColor3d("LightSteelBlue"))

        ren.AddActor(actor)

    ren.SetBackground(colors.GetColor3d("White"))

    iren.Initialize()
    renWin.Render()
    iren.Start()


if __name__ == "__main__":
    main()
