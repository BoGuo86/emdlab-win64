# run_geo.py
import gmsh
import sys

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.option.setNumber("Mesh.Smoothing", 10)

# Load .geo file
gmsh.open("C:\\emdlab-win64\\tmp\\emdlab_gmsh_geoFile.geo") 

# Generate 2D mesh
gmsh.model.mesh.generate(2)

# Write mesh to file
gmsh.write("C:\\emdlab-win64\\tmp\\emdlab_gmsh_mshFile.m")

gmsh.finalize()
