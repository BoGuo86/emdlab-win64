# run_geo.py
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# Load .geo file
gmsh.open("C:\\emdlab-win64\\tmp\\emdlab_gmsh_geoFile.geo") 

# Write step file
gmsh.write("C:\\emdlab-win64\\tmp\\emdlab_gmsh_stpFile.step")

# Write stl file
gmsh.write("C:\\emdlab-win64\\tmp\\emdlab_gmsh_stlFile.stl")

gmsh.finalize()
