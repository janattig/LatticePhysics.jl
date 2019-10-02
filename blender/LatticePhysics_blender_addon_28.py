# imports
import bpy
import os

from bpy.props import StringProperty, BoolProperty, IntProperty, FloatProperty
from bpy_extras.io_utils import ImportHelper
from bpy.types import Operator

import math
import mathutils
from mathutils import *
from math import *


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



# HELPER FUNCTIONS
def find_collection(context, item):
    collections = item.users_collection
    if len(collections) > 0:
        return collections[0]
    return context.scene.collection

def make_collection(collection_name, parent_collection):
    if collection_name in bpy.data.collections: # Does the collection already exist?
        return bpy.data.collections[collection_name]
    else:
        new_collection = bpy.data.collections.new(collection_name)
        parent_collection.children.link(new_collection) # Add the new collection under a parent
        return new_collection





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
    radius_sites = FloatProperty(
        name = "Site radius (relative)",
        default = 0.2,
        min = 0.0,
        max = 10.0
    )
    # number of subdivisions of sites (passed to the subsurf modifier)
    subdivision_sites = IntProperty(
        name = "Subdivision of sites",
        default = 2,
        min = 0,
        max = 8
    )
    radius_bonds = FloatProperty(
        name = "Bond radius (relative)",
        default = 0.08,
        min = 0.0,
        max = 1.0
    )
    # number of subdivisions of sites (passed to the subsurf modifier)
    subdivision_bonds = IntProperty(
        name = "Subdivision of bonds",
        default = 2,
        min = 0,
        max = 8
    )

    # if there is a subsurf modifier for sites
    hook_bonds_to_sites = BoolProperty(
        name = "Hook bonds to sites",
        default = True
    )




    ##########---------------------------------
    # STEP 2 #  FUNCTIONS FOR ADDING PARTS
    ##########---------------------------------

    # DEFINE A FUNCTION TO ADD A SPHERE
    def addSite(self, site_index, p, radius, lbl):
        # set the correct position
        x = bpy.context.scene.cursor.location[0] + p[0]
        y = bpy.context.scene.cursor.location[1] + p[1]
        z = bpy.context.scene.cursor.location[2] + p[2]
        # add an ico sphere at that point
        bpy.ops.mesh.primitive_ico_sphere_add(subdivisions=self.resolution_sites, radius=radius, location=(x,y,z))
        # get the current object
        obj = bpy.context.active_object
        # set the name
        obj.name = "site_{0}".format(site_index)
        # add subsurf modifier to make surface look smoother
        if self.subdivision_sites > 0:
            obj.modifiers.new("site_subsurf", "SUBSURF")
            obj.modifiers["site_subsurf"].levels = 1
            obj.modifiers["site_subsurf"].render_levels = self.subdivision_sites
        # set the material (color)
        material_name = "site_material_"+lbl
        # get the material
        material = bpy.data.materials.get(material_name)
        if material is None:
            # create material
            material = bpy.data.materials.new(name=material_name)
            # set the color
            material.diffuse_color = (1.0,1.0,1.0,1.0)
            # use nodes for rendering
            material.use_nodes = True
        # set the material
        if obj.data.materials:
            obj.data.materials[0] = material
        else:
            obj.data.materials.append(material)
        # set smooth shading
        bpy.ops.object.shade_smooth()
        # deselect the object
        obj.select_set(False)
        # return the object
        return obj


    # DEFINE A FUNCTION TO ADD A TUBE
    def addBond(self, bond_index, p_from, p_to, radius, lbl, hook_from=-1, hook_to=-1):
        # set the correct positions
        x_from = bpy.context.scene.cursor.location[0] + p_from[0]
        y_from = bpy.context.scene.cursor.location[1] + p_from[1]
        z_from = bpy.context.scene.cursor.location[2] + p_from[2]
        x_to   = bpy.context.scene.cursor.location[0] + p_to[0]
        y_to   = bpy.context.scene.cursor.location[1] + p_to[1]
        z_to   = bpy.context.scene.cursor.location[2] + p_to[2]
        # adding curve from
        # https://blender.stackexchange.com/questions/120074/how-to-make-a-curve-path-from-scratch-given-a-list-of-x-y-z-points
        # make a new curve
        crv = bpy.data.curves.new("crv_{0}".format(bond_index), 'CURVE')
        crv.dimensions = '3D'
        # make a new spline in that curve
        spline = crv.splines.new(type='NURBS')
        # a spline point for each point
        spline.points.add(1) # theres already one point by default
        # assign the point coordinates to the spline points
        spline.points[0].co = (x_from,y_from,z_from,1)
        spline.points[1].co = (  x_to,  y_to,  z_to,1)
        # make a new object with the curve
        obj = bpy.data.objects.new("bond_{0}".format(bond_index), crv)
        # TODO render bevel (Probably change in the future)
        obj.data.splines[0].use_bezier_u = True
        obj.data.splines[0].use_bezier_u = False
        # set the bevel so that the bonds are actually tubes
        obj.data.bevel_depth = radius
        # link it to the scene
        bpy.context.scene.collection.objects.link(obj)
        # maybe hook the bond to the sites
        if self.hook_bonds_to_sites and hook_from>=0 and hook_to>=0:
            # add the hook modifiers
            obj.modifiers.new("follow_site_{0}".format(hook_from), "HOOK")
            obj.modifiers["follow_site_{0}".format(hook_from)].object = bpy.data.objects["site_{0}".format(hook_from)]
            obj.modifiers["follow_site_{0}".format(hook_from)].vertex_indices_set([0,])
            obj.modifiers.new("follow_site_{0}".format(hook_to), "HOOK")
            obj.modifiers["follow_site_{0}".format(hook_to)].object   = bpy.data.objects["site_{0}".format(hook_to)]
            obj.modifiers["follow_site_{0}".format(hook_to)].vertex_indices_set([1,])
        # add subsurf modifier to make surface look smoother
        if self.subdivision_bonds > 0:
            # then subsurf
            obj.modifiers.new("bond_subsurf", "SUBSURF")
            obj.modifiers["bond_subsurf"].levels = 1
            obj.modifiers["bond_subsurf"].render_levels = self.subdivision_bonds
        # set the material (color)
        material_name = "bond_material_"+lbl
        # get the material
        material = bpy.data.materials.get(material_name)
        if material is None:
            # create material
            material = bpy.data.materials.new(name=material_name)
            # set the color
            material.diffuse_color = (0.5,0.5,0.5,1.0)
            # use nodes for rendering
            material.use_nodes = True
        # set the material
        if obj.data.materials:
            obj.data.materials[0] = material
        else:
            obj.data.materials.append(material)
        # set smooth shading
        bpy.ops.object.shade_smooth()
        # select the object
        #obj.select_set(True)
        # set origin to geometry
        #bpy.ops.object.origin_set(type='ORIGIN_GEOMETRY', center='MEDIAN')
        # deselect the object
        #obj.select_set(False)
        # return the object
        return obj





    ##########---------------------------------
    # STEP 3 #  FUNCTION FOR LOADING LATTICE
    ##########---------------------------------
    def loadLatticeDump(self, filename):

        # read in file
        with open(filename) as f:
            lines = f.readlines()

        # refactor filename
        fn = filename.split("/")[-1]
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
                site_y = float(l.split(" ")[4].split("=")[1][1:-1])
                site_z = float(l.split(" ")[5].split("=")[1][1:-1])
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
                bond_f = int(l.split(" ")[3].split("=")[1][1:-1]) - 1
                bond_t = int(l.split(" ")[4].split("=")[1][1:-1]) - 1
                # add the site information to the lists
                bonds_l[bond_i] = bond_l
                bonds_f[bond_i] = bond_f
                bonds_t[bond_i] = bond_t
                #print("bond found, index=", bond_i, ", label is \"", bonds_l[bond_i],"\"")
                # get the data from the line
                # add the bond
                #self.addTube(int(bond_data[0]),bond_data[1],bond_data[2],bond_data[3],bond_data[4],bond_data[5],bond_data[6],bond_data[7],l.split("\t")[2])


        # center all sites
        cm_x = sum(sites_x) / len(sites_x)
        cm_y = sum(sites_y) / len(sites_y)
        cm_z = sum(sites_z) / len(sites_z)
        sites_x = [s - cm_x for s in sites_x]
        sites_y = [s - cm_y for s in sites_y]
        sites_z = [s - cm_z for s in sites_z]

        # determine lenght of all bonds
        bond_lengths = [sqrt((sites_x[bonds_f[i]]-sites_x[bonds_t[i]])**2 + (sites_y[bonds_f[i]]-sites_y[bonds_t[i]])**2) for i in range(N_bonds)]
        min_bond_length = min(bond_lengths)

        # deselect all objects
        bpy.ops.object.select_all(action='DESELECT')

        # check if lattice already added
        if "Lattice {0}".format(fn) in bpy.data.collections: # Does the top level collection already exist?
            print("ERROR: Lattice already loaded")
            return

        # add all sites and collect objects into list
        site_objects = [ self.addSite(s+1, [sites_x[s],sites_y[s],sites_z[s]], self.radius_sites * min_bond_length, sites_l[s]) for s in range(N_sites) ]

        # add all bonds and collect objects into list
        bond_objects = [ self.addBond(b+1, [sites_x[bonds_f[b]],sites_y[bonds_f[b]],sites_z[bonds_f[b]]], [sites_x[bonds_t[b]],sites_y[bonds_t[b]],sites_z[bonds_t[b]]], self.radius_bonds * min_bond_length, bonds_l[b], bonds_f[b]+1, bonds_t[b]+1) for b in range(N_bonds) ]

        # make the collections right
        general_collection = bpy.data.collections.new("Lattice {0}".format(fn))
        bpy.context.scene.collection.children.link(general_collection)
        # General site and bond collections
        site_collection = make_collection("Sites ({0})".format(fn), general_collection)
        bond_collection = make_collection("Bonds ({0})".format(fn), general_collection)

        site_collections = dict()
        bond_collections = dict()
        # make new collections for the respective labels
        for s in sites_l:
            site_collections[s] = make_collection("Sites (label={0}) ({1})".format(s, fn), site_collection)
        for b in bonds_l:
            bond_collections[b] = make_collection("Bonds (label={0}) ({1})".format(b, fn), bond_collection)
        # push all sites into theses lists
        for s in range(N_sites):
            old_site_collection = find_collection(bpy.context, site_objects[s])
            site_collections[sites_l[s]].objects.link(site_objects[s])
            old_site_collection.objects.unlink(site_objects[s])
        # push all bonds into theses lists
        for b in range(N_bonds):
            old_bond_collection = find_collection(bpy.context, bond_objects[b])
            bond_collections[bonds_l[b]].objects.link(bond_objects[b])
            old_bond_collection.objects.unlink(bond_objects[b])


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
