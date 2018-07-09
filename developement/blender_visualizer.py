# imports
import bpy
import math
import mathutils



# DEFINE A FUNCTION TO ADD A SPHERE
def addSphere(site_index, x,y,z, radius, color):
    # set the correct position
    x = bpy.context.scene.cursor_location[0] + x
    y = bpy.context.scene.cursor_location[1] + y
    z = bpy.context.scene.cursor_location[2] + z
    # add an ico sphere at that point
    bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=2, size=radius, location=(x,y,z))
    # get the current object
    obj = bpy.context.active_object
    # set the name
    obj.name = "site_{0}".format(site_index)
    # add subsurf modifier to make surface look smoother
    obj.modifiers.new("subd", "SUBSURF")
    obj.modifiers["subd"].levels = 2
    obj.modifiers["subd"].render_levels = 3
    # set smooth shading
    bpy.ops.object.shade_smooth()

# DEFINE A FUNCTION TO ADD A TUBE
def addTube(bond_index, x_from,y_from,z_from, x_to,y_to,z_to, radius, color):
    # set the correct positions
    x_from = bpy.context.scene.cursor_location[0] + x_from
    y_from = bpy.context.scene.cursor_location[1] + y_from
    z_from = bpy.context.scene.cursor_location[2] + z_from
    x_to   = bpy.context.scene.cursor_location[0] + x_to
    y_to   = bpy.context.scene.cursor_location[1] + y_to
    z_to   = bpy.context.scene.cursor_location[2] + z_to
    # difference in coordinates
    dx = x_to - x_from
    dy = y_to - y_from
    dz = z_to - z_from
    # angles
    phi   = math.atan2(dy, dx)
    theta = math.pi/2 - math.atan2(dz,math.sqrt(dx**2+dy**2))
    # find the length of the bond
    bond_length = math.sqrt(dx**2 + dy**2 + dz**2)
    # add a cylinder of the desired data
    bpy.ops.mesh.primitive_cylinder_add(vertices=16, radius=radius, depth=bond_length, location=((x_to+x_from)/2, (y_to+y_from)/2, (z_to+z_from)/2))
    bpy.ops.transform.rotate(value=theta, axis=(0,1,0))
    bpy.ops.transform.rotate(value=phi, axis=(0,0,1))
    # get the current object
    obj = bpy.context.active_object
    # set the name
    obj.name = "bond_{0}".format(bond_index)
    # set smooth shading
    bpy.ops.object.shade_smooth()




# file
filename = "/home/jattig/PackageDevelopement/LatticePhysics.jl/test_2.txt"

# read in file
with open(filename) as f:
    lines = f.readlines()

# remove whitespace characters like `\n` at the end of each line
lines = [x.strip() for x in lines]


# go through all lines and check if a site or bond has to be added
for l in lines:
    # line specifying a site
    if l.startswith("site:\t"):
        # get the data from the line
        sphere_data = [float(x) for x in l.split("\t")[1].split(", ")]
        # add the site
        addSphere(int(sphere_data[0]),sphere_data[1],sphere_data[2],sphere_data[3],sphere_data[4],[int(sphere_data[5]),int(sphere_data[6]),int(sphere_data[7])])
    # line specifying a bond
    if l.startswith("bond:\t"):
        # get the data from the line
        bond_data = [float(x) for x in l.split("\t")[1].split(", ")]
        # add the bond
        addTube(int(bond_data[0]),bond_data[1],bond_data[2],bond_data[3],bond_data[4],bond_data[5],bond_data[6],bond_data[7],[int(bond_data[8]),int(bond_data[9]),int(bond_data[10])])
