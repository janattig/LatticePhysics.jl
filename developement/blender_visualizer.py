# imports
import bpy
import math
import mathutils



# DEFINE A FUNCTION TO ADD A SPHERE
def addSphere(x,y,z, radius):
    x = bpy.context.scene.cursor_location[0] + x
    y = bpy.context.scene.cursor_location[1] + y
    z = bpy.context.scene.cursor_location[2] + z
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, size=radius, location=(x,y,z))



# file
filename = "/home/jattig/PackageDevelopement/LatticePhysics.jl/test_2.txt"

# read in file
with open(filename) as f:
    lines = f.readlines()

# you may also want to remove whitespace characters like `\n` at the end of each line
lines = [x.strip() for x in lines]

for l in lines:
    if l.startswith("site:\t"):
        sphere_data = [float(x) for x in l.split("\t")[1].split(", ")]
        addSphere(sphere_data[0],sphere_data[1],sphere_data[2],sphere_data[3])

for obj in bpy.data.objects:
    if "subd" not in obj.modifiers:
        obj.modifiers.new("subd", "SUBSURF")

for obj in bpy.data.objects:
    if "subd" in obj.modifiers:
        obj.modifiers["subd"].levels = 2
        obj.modifiers["subd"].render_levels = 5
        print(obj)