# imports
import bpy
import os

from bpy.props import StringProperty, BoolProperty, IntProperty
from bpy_extras.io_utils import ImportHelper
from bpy.types import Operator

import math
import mathutils


# BLENDER ADDON INFORMATION
bl_info = {
    "name": "LatticePhysics Support",
    "author": "Jan Attig",
    "version": (2, 0, 0),
    "blender": (2, 80, 0),
    "location": "Import-Export",
    "description": "Loading LatticePhysics XML-dumps of lattices",
    "category": "Import-Export",
}









##########---------------------------------
# STEP 1 #  BLENDER ADDON DEFINITION
##########---------------------------------
class LatticePhysicsBlenderAddon(Operator, ImportHelper):

    # id and label
    bl_idname = "latticephysics.open_filebrowser"
    bl_label = "Import LatticePhysics dump file"

    # global filter options
    filter_glob = StringProperty(
        default='*.lpbd;*.txt;*.xml',
        options={'HIDDEN'}
    )


    # Options for importing

    # resolution of sites (passed to icosphere)
    resolution_sites = IntProperty(
        name = "Site resolution",
        default = 2,
        min = 1,
        max = 16
    )

    # number of subdivisions of sites (passed to the subsurf modifier)
    subdivision_sites = IntProperty(
        name = "Site subdivision",
        default = 2,
        min = 0,
        max = 8
    )
    # number of subdivisions of sites (passed to the subsurf modifier)
    subdivision_bonds = IntProperty(
        name = "Bond subdivision",
        default = 2,
        min = 0,
        max = 8
    )





    ##########---------------------------------
    # STEP 2 #  FUNCTIONS FOR ADDING PARTS
    ##########---------------------------------


    # DEFINE A FUNCTION TO ADD A MATERIAL
    def addMaterial(self, name, color):
        # get the material
        material = bpy.data.materials.get(name)
        if material is None:
            # create material
            material = bpy.data.materials.new(name=name)
        # set the color
        material.diffuse_color = (color[0],color[1],color[2],1.0)
        # use nodes for rendering
        material.use_nodes = True
        # clear the default nodes
        material.node_tree.nodes.clear()
        # add the necessary nodes
        diffuse = material.node_tree.nodes.new(type = 'ShaderNodeBsdfDiffuse')
        glossy  = material.node_tree.nodes.new(type = 'ShaderNodeBsdfGlossy')
        output  = material.node_tree.nodes.new(type = 'ShaderNodeOutputMaterial')
        mixer   = material.node_tree.nodes.new(type = 'ShaderNodeMixShader')
        fresnel = material.node_tree.nodes.new(type = 'ShaderNodeFresnel')
        # link the nodes
        material.node_tree.links.new(fresnel.outputs['Fac'],    mixer.inputs['Fac'])
        material.node_tree.links.new(diffuse.outputs['BSDF'],   mixer.inputs[1])
        material.node_tree.links.new( glossy.outputs['BSDF'],   mixer.inputs[2])
        material.node_tree.links.new(  mixer.outputs['Shader'], output.inputs['Surface'])
        material.node_tree.links.new(  mixer.outputs['Shader'], output.inputs['Surface'])
        # set some default values
        diffuse.inputs[0].default_value           = (color[0],color[1],color[2],1.0)  # color of the material
        diffuse.inputs['Roughness'].default_value = 0.75  # Roughness of the diffuse part of material
        glossy.inputs['Roughness'].default_value  = 0.05  # Roughness of the glossy part of material
        glossy.inputs[0].default_value            = (1.0,1.0,1.0, 1.0)  # color of the gloss
        glossy.distribution                       = "MULTI_GGX"
        # fix location of nodes
        output.location  = ( 200, 0)
        fresnel.location = (-200, 150)
        diffuse.location = (-200, 0)
        glossy.location  = (-200, -150)


    # DEFINE A FUNCTION TO ADD A SPHERE
    def addSphere(self, site_index, x,y,z, radius, color):
        # set the correct position
        x = bpy.context.scene.cursor.location[0] + x
        y = bpy.context.scene.cursor.location[1] + y
        z = bpy.context.scene.cursor.location[2] + z
        # add an ico sphere at that point
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=self.resolution_sites, radius=radius, location=(x,y,z))
        # get the current object
        obj = bpy.context.active_object
        # set the name
        obj.name = "site_{0}".format(site_index)
        # add subsurf modifier to make surface look smoother
        if self.subdivide_sites:
            obj.modifiers.new("site_subsurf", "SUBSURF")
            obj.modifiers["site_subsurf"].levels = 1
            obj.modifiers["site_subsurf"].render_levels = self.subdivision_sites
        # set the material (color)
        material = bpy.data.materials.get(color)
        if obj.data.materials:
            obj.data.materials[0] = material
        else:
            obj.data.materials.append(material)
        # set smooth shading
        bpy.ops.object.shade_smooth()


    # DEFINE A FUNCTION TO ADD A TUBE
    def addTube(self, bond_index, x_from,y_from,z_from, x_to,y_to,z_to, radius, color):
        # set the correct positions
        x_from = bpy.context.scene.cursor.location[0] + x_from
        y_from = bpy.context.scene.cursor.location[1] + y_from
        z_from = bpy.context.scene.cursor.location[2] + z_from
        x_to   = bpy.context.scene.cursor.location[0] + x_to
        y_to   = bpy.context.scene.cursor.location[1] + y_to
        z_to   = bpy.context.scene.cursor.location[2] + z_to
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
        bpy.ops.mesh.primitive_cylinder_add(vertices=self.resolution_bonds, radius=radius, depth=bond_length, location=((x_to+x_from)/2, (y_to+y_from)/2, (z_to+z_from)/2))
        bpy.ops.transform.rotate(value=-theta, orient_axis='Y')
        bpy.ops.transform.rotate(value=-phi, orient_axis='Z')
        # get the current object
        obj = bpy.context.active_object
        # set the name
        obj.name = "bond_{0}".format(bond_index)
        # add subsurf modifier to make surface look smoother
        if self.subdivide_bonds:
            # first edge split (to seperate the faces of the tubes from subsurfing)
            obj.modifiers.new("bond_edge_split", "EDGE_SPLIT")
            obj.modifiers["bond_edge_split"].split_angle = 1.0
            # then subsurf
            obj.modifiers.new("bond_subsurf", "SUBSURF")
            obj.modifiers["bond_subsurf"].levels = 1
            obj.modifiers["bond_subsurf"].render_levels = self.subdivision_bonds
        # set the material (color)
        material = bpy.data.materials.get(color)
        if obj.data.materials:
            obj.data.materials[0] = material
        else:
            obj.data.materials.append(material)
        # set smooth shading
        bpy.ops.object.shade_smooth()


    ##########---------------------------------
    # STEP 3 #  FUNCTION FOR LOADING LATTICE
    ##########---------------------------------
    def loadLatticeDump(self, filename):

        # read in file
        with open(filename) as f:
            lines = f.readlines()

        # remove whitespace characters like `\n` at the end of each line
        lines = [x.strip() for x in lines]


        # set the render engine
        bpy.context.scene.render.engine = "CYCLES"

        # set the units
        bpy.context.scene.unit_settings.system = "METRIC"
        bpy.context.scene.unit_settings.system_rotation = "DEGREES"

        # try to set to filmic color space
        try:
            bpy.context.scene.view_settings.view_transform = "Filmic"
        except:
            print("Filmic color space not available")


        # try to obtain lists for sites and bonds
        if not lines[0].startswith("<LATTICEGRAPH"):
            print("ERROR, not a Lattice Graph file")

        # make lists for sites and bonds
        N_sites = int(lines[0].split(" ")[1].split("=")[1][1:-1])
        N_bonds = int(lines[0].split(" ")[2].split("=")[1][1:-1])
        print(N_sites, " sites found")
        print(N_bonds, " bonds found")
        sites_x = [0.0 for s in range(N_sites)]
        sites_y = [0.0 for s in range(N_sites)]
        sites_z = [0.0 for s in range(N_sites)]
        sites_l = ["l" for s in range(N_sites)]
        bonds_f = [  0 for b in range(N_bonds)]
        bonds_t = [  0 for b in range(N_bonds)]
        bonds_l = ["l" for b in range(N_bonds)]

        # go through all lines and check if a site or bond has to be added
        for l in lines:
            # line specifying a site
            if l.startswith("<SITE"):
                # sample line:
                # <SITE index="27" label="4" X="2.4748737341529163" Y="-1.0606601717798216" Z="1.7677669529663687" >
                # get the data from the line
                site_i = int(l.split(" ")[1].split("=")[1][1:-1]) - 1
                site_l = l.split(" ")[2].split("=")[1][1:-1]
                site_x = float(l.split(" ")[3].split("=")[1][1:-1])
                site_y = float(l.split(" ")[3].split("=")[1][1:-1])
                site_z = float(l.split(" ")[3].split("=")[1][1:-1])
                # add the site information to the lists
                sites_l[site_i] = site_l
                sites_x[site_i] = site_x
                sites_y[site_i] = site_y
                sites_z[site_i] = site_z
                #print("site found, index=", site_i, ", label is \"", sites_l[site_i],"\"")
            # line specifying a bond
            if l.startswith("<BOND"):
                # sample line:
                # <BOND index="1" label="1" from="1" to="2" >
                # get the data from the line
                bond_i = int(l.split(" ")[1].split("=")[1][1:-1]) - 1
                bond_l = l.split(" ")[2].split("=")[1][1:-1]
                bond_f = float(l.split(" ")[3].split("=")[1][1:-1])
                bond_t = float(l.split(" ")[4].split("=")[1][1:-1])
                # add the site information to the lists
                bonds_l[bond_i] = bond_l
                bonds_f[bond_i] = bond_f
                bonds_t[bond_i] = bond_t
                #print("bond found, index=", bond_i, ", label is \"", bonds_l[bond_i],"\"")
                # get the data from the line
                # add the bond
                #self.addTube(int(bond_data[0]),bond_data[1],bond_data[2],bond_data[3],bond_data[4],bond_data[5],bond_data[6],bond_data[7],l.split("\t")[2])





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
        self.loadLatticeDump(self.filepath)

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
