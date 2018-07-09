# imports
import bpy
import os

from bpy.props import StringProperty, BoolProperty
from bpy_extras.io_utils import ImportHelper
from bpy.types import Operator

import math
import mathutils


# BLENDER ADDON INFORMATION
bl_info = {
    "name": "LatticePhysics Support",
    "author": "Jan Attig",
    "version": (1, 0, 0),
    "blender": (2, 76, 0),
    "location": "Import-Export",
    "description": "Loading LatticePhysics dumps of lattices",
    "category": "Import-Export",
}






##########---------------------------------
# STEP 1 #  FUNCTIONS FOR ADDING PARTS
##########---------------------------------


# DEFINE A FUNCTION TO ADD A MATERIAL
def addMaterial(name, color):
    # get the material
    material = bpy.data.materials.get(name)
    if material is None:
        # create material
        material = bpy.data.materials.new(name=name)
    # set the color
    material.diffuse_color = color
    # use nodes for rendering
    material.use_nodes = True
    # clear the default nodes
    material.node_tree.nodes.clear()
    # add the necessary nodes
    diffuse = material.node_tree.nodes.new(type = 'ShaderNodeBsdfDiffuse')
    glossy  = material.node_tree.nodes.new(type = 'ShaderNodeBsdfGlossy')
    output  = material.node_tree.nodes.new(type = 'ShaderNodeOutputMaterial')
    mixer   = material.node_tree.nodes.new(type = 'ShaderNodeMixShader')
    # link the nodes
    material.node_tree.links.new(diffuse.outputs['BSDF'],   mixer.inputs[1])
    material.node_tree.links.new( glossy.outputs['BSDF'],   mixer.inputs[2])
    material.node_tree.links.new(  mixer.outputs['Shader'], output.inputs['Surface'])
    # set some default values
    mixer.inputs[0].default_value             = 0.15  # mixing between diffuse and glossy
    diffuse.inputs[0].default_value           = (color[0],color[1],color[2],1.0)  # color of the material
    diffuse.inputs['Roughness'].default_value = 0.015 # Roughness of the diffuse part of material
    glossy.inputs['Roughness'].default_value  = 0.005 # Roughness of the glossy part of material
    glossy.inputs[0].default_value            = (1.0,1.0,1.0, 1.0)  # color of the gloss
    # fix location of nodes
    output.location = (200, 0)
    diffuse.location = (-200, 100)
    glossy.location = (-200, -100)


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
    # set the material (color)
    material = bpy.data.materials.get(color)
    if obj.data.materials:
        obj.data.materials[0] = material
    else:
        obj.data.materials.append(material)
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
    # set the material (color)
    material = bpy.data.materials.get(color)
    if obj.data.materials:
        obj.data.materials[0] = material
    else:
        obj.data.materials.append(material)
    # set smooth shading
    bpy.ops.object.shade_smooth()



##########---------------------------------
# STEP 2 #  FUNCTION FOR LOADING LATTICE
##########---------------------------------
def loadLatticeDump(filename):

    # read in file
    with open(filename) as f:
        lines = f.readlines()

    # remove whitespace characters like `\n` at the end of each line
    lines = [x.strip() for x in lines]

    # set the render engine
    bpy.context.scene.render.engine = "CYCLES"

    # go through all lines and check if a site or bond has to be added
    for l in lines:
        # line specifying a material
        if l.startswith("material:\t"):
            # Get material properties
            name=l.split("\t")[1]
            color = [float(x)/255.0 for x in l.split("\t")[2].split(", ")]
            # add material
            addMaterial(name, color)
        # line specifying a site
        if l.startswith("site:\t"):
            # get the data from the line
            sphere_data = [float(x) for x in l.split("\t")[1].split(", ")]
            # add the site
            addSphere(int(sphere_data[0]),sphere_data[1],sphere_data[2],sphere_data[3],sphere_data[4],l.split("\t")[2])
        # line specifying a bond
        if l.startswith("bond:\t"):
            # get the data from the line
            bond_data = [float(x) for x in l.split("\t")[1].split(", ")]
            # add the bond
            addTube(int(bond_data[0]),bond_data[1],bond_data[2],bond_data[3],bond_data[4],bond_data[5],bond_data[6],bond_data[7],l.split("\t")[2])





##########---------------------------------
# STEP 3 #  BLENDER ADDON DEFINITION
##########---------------------------------
class LatticePhysicsBlenderAddon(Operator, ImportHelper):

    # id and label
    bl_idname = "latticephysics.open_filebrowser"
    bl_label = "Import LatticePhysics dump file"

    # global filter options
    filter_glob = StringProperty(
        default='*.lpbd;*.txt',
        options={'HIDDEN'}
    )

    # EXECUTE THE ADDON
    def execute(self, context):
        """Import selected LatticePhysics lattice file into Blender"""

        # get filename and extension
        filename, extension = os.path.splitext(self.filepath)

        # print the selected results
        print('Selected file:', self.filepath)
        print('File name:', filename)
        print('File extension:', extension)

        # add lattice into blender
        addLatticeDump(self.filepath)

        # return
        return {'FINISHED'}


# REGISTER FUNCTION (to register addon)
def register():
    bpy.utils.register_class(LatticePhysicsBlenderAddon)


# UNREGISTER FUNCTION (to unregister addon)
def unregister():
    bpy.utils.unregister_class(LatticePhysicsBlenderAddon)


# Call register function if the file is executed
if __name__ == "__main__":
    register()
