# class for rendering a sphere (Earth) with blender
#
import sys
import os
import time

# blender
import bpy

print("")
print("blender: version ",bpy.app.version_string)
print("")

###############################################################################################

## needs additional python paths to import scipy, matplotlib

# see also: adding option --python-use-system-env to use PYTHONPATH
paths = sys.path.copy()
#print("paths: ",paths)
envs = os.environ.copy()
#print("envs: ",envs)

# BLENDER_SYSTEM_PATH
# check if BLENDER_SYSTEM_PATH has been set
# for example:
#  export BLENDER_SYSTEM_PYTHON='/opt/local/Library/Frameworks/Python.framework/Versions/3.10/'
if 'BLENDER_SYSTEM_PYTHON' in envs.keys():
    print("setting: using BLENDER_SYSTEM_PYTHON")
else:
    print("setting: adding local python path")
    sys.path.clear()
    sys.path.append('/opt/local/Library/Frameworks/Python.framework/Versions/3.10/lib/python3.10/site-packages')
    for path in paths:
        #print("add path ",path)
        sys.path.append(path)

#print("path:")
#for path in sys.path: print("  ",path)

###############################################################################################

# Constants
PI = 3.141592653589793
DEGREE_TO_RAD = PI / 180.0

###############################################################################################


# class to avoid long stdout output by renderer
# see: https://stackoverflow.com/questions/24277488/in-python-how-to-capture-the-stdout-from-a-c-shared-library-to-a-variable/29834357
class SuppressStream(object):
    def __init__(self, stream=sys.stderr,suppress=False):
        # turns on/off suppressing stdout of renderer process
        self.SUPPRESS_STDOUT = suppress

        if self.SUPPRESS_STDOUT:
            self.orig_stream_fileno = stream.fileno()

    def __enter__(self):
        if self.SUPPRESS_STDOUT:
            self.orig_stream_dup = os.dup(self.orig_stream_fileno)
            self.devnull = open(os.devnull, 'w')
            os.dup2(self.devnull.fileno(), self.orig_stream_fileno)

    def __exit__(self, type, value, traceback):
        if self.SUPPRESS_STDOUT:
            os.close(self.orig_stream_fileno)
            os.dup2(self.orig_stream_dup, self.orig_stream_fileno)
            os.close(self.orig_stream_dup)
            self.devnull.close()


###############################################################################################

# blender class for sphere rendering

class blender_sphere(object):
    """
    class for rendering (earth) sphere with blender
    """
    def __init__(self,
                 verbose=False,
                 texture_globe="",texture_clouds="",texture_night="",texture_topo="",
                 img_size_X=100,
                 img_size_Y=100,
                 animation=False,
                 animation_rotation_degree=0.1):

        ## initializations
        # verbosity
        self.verbose = verbose

        # render image size
        self.img_size_X = img_size_X
        self.img_size_Y = img_size_Y

        # textures
        self.texture_globe  = texture_globe
        self.texture_clouds = texture_clouds
        self.texture_night  = texture_night
        self.texture_topo   = texture_topo

        # animation info
        self.animation = animation
        self.animation_number_of_keyframes = 10
        self.animation_keyframe_interval   = 10
        self.animation_rotation_increment  = animation_rotation_degree * DEGREE_TO_RAD

        # rendering settings
        self.blender_engine = 'BLENDER_EEVEE'  # 'BLENDER_EEVEE','CYCLES'
        self.blender_device = 'CPU'            # 'CPU','GPU'

        self.blender_img_resolution_X = 600
        self.blender_img_resolution_Y = 600

        # cleanup
        self.clean_scene()

    def __str__(self):
        info = "helper class for rendering globe with blender"

        return info


    def clean_scene(self):
        """
        delete all existing items in scene
        """
        if self.verbose:
            print("scene: clean default objects")

        bpy.ops.object.select_all(action='SELECT')
        bpy.ops.object.delete(use_global=False)


    def add_sphere(self):
        """
        creates sphere
        """
        position = (0, 0, 0)

        if self.verbose:
            print("sphere: ",position)

        bpy.ops.mesh.primitive_uv_sphere_add(enter_editmode=False, align='WORLD',
                                             location=position, scale=(1.0, 1.0, 1.0))

        # current object
        sphere = bpy.context.object
        sphere.name = "Sphere"

        # set quality to highest and apply shade smooth.
        bpy.ops.object.modifier_add(type='SUBSURF')
        sphere.modifiers["Subdivision"].quality = 4
        sphere.modifiers["Subdivision"].levels = 4
        sphere.modifiers["Subdivision"].render_levels = 4
        bpy.ops.object.modifier_apply(modifier="Subdivision")
        bpy.ops.object.shade_smooth()

        # adds material
        self.add_sphere_material(sphere)


    def add_sphere_material(self,sphere):
        """
        adds material to sphere
        """
        # simple example material
        #mat = bpy.data.materials.new(name="Material")
        #mat.diffuse_color = [0.5,0.8,1.0,1.0]
        #mat.specular_intensity = 0.8

        # create material
        mat = bpy.data.materials.new(name="Material Globe")

        # enable node-graph edition mode
        mat.use_nodes = True

        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        # clears out starter nodes
        if self.verbose:
            print("material:")
            print("  nodes: ",nodes.keys())

        # sphere image
        self.add_sphere_image(nodes,links)

        # topography
        self.add_sphere_topo(nodes,links)

        # night image
        self.add_sphere_night(nodes,links)

        ## specular emission
        if False:
            bsdf = nodes["Principled BSDF"]
            #links.new(bsdf.inputs['Metallic'], texTopo.outputs['Color'])
            # sets metallic surface to 1
            #bpy.data.materials["Material"].node_tree.nodes["Principled BSDF"].inputs[4].default_value = 1
            mat.node_tree.links.new(bsdf.inputs['Specular'], texImage.outputs['Color'])

        # adds material to sphere
        sphere.data.materials.append(mat)
        # gets active sphere object
        #obj = bpy.context.view_layer.objects.active
        # adds material to sphere
        #obj.data.materials.append(mat)


    def add_sphere_image(self,nodes,links):
        """
        adds image texture
        """
        # checks if anything to do
        if not self.texture_globe:
            print("no globe image")
            return

        # user output
        if self.verbose:
            print("image:")
            print("  loading texture: ",self.texture_globe)

        # checks if file exists
        if not os.path.isfile(self.texture_globe):
            print("Please check if texture file exists: ",self.texture_globe)
            sys.exit(1)

        texImage = nodes.new('ShaderNodeTexImage')
        texImage.image = bpy.data.images.load(self.texture_globe)

        # adds Principled BSDF shader
        bsdf = nodes["Principled BSDF"]
        # texture as coloring
        links.new(bsdf.inputs['Base Color'], texImage.outputs['Color'])
        links.new(bsdf.inputs['Emission'], texImage.outputs['Color'])
        bsdf.inputs['Emission Strength'].default_value = 0.1

        bsdf.inputs['Metallic'].default_value = 0.1
        bsdf.inputs['Roughness'].default_value = 0.8
        bsdf.inputs['Specular'].default_value = 0.2


    def add_sphere_topo(self,nodes,links):
        """
        adds topography texture
        """
        # checks if anything to do
        if not self.texture_topo:
            print("no topo")
            return

        # user output
        if self.verbose:
            print("topo:")
            print("  loading texture: ",self.texture_topo)

        # checks if file exists
        if not os.path.isfile(self.texture_topo):
            print("Please check if texture file exists: ",self.texture_topo)
            sys.exit(1)

        texTopo = nodes.new('ShaderNodeTexImage')
        texTopo.image = bpy.data.images.load(self.texture_topo)

        # sets colorspace to Non-Color
        texTopo.image.colorspace_settings.name = 'Non-Color'
        #texTopo.projection = 'SPHERE'
        texTopo.interpolation = 'Cubic' # 'Cubic', 'Linear'
        #texTopo.extension = 'EXTEND'

        # adds shader
        #shader_mat = nodes["Material Output"]
        bump_shader = nodes.new('ShaderNodeBump')
        bump_shader.inputs['Strength'].default_value = 0.1
        links.new(bump_shader.inputs['Height'], texTopo.outputs['Color'])
        #bump_shader.inputs['Height'].default_value = 1.0
        #bump_shader.inputs['Distance'].default_value = 1000.0
        #seperate_rgb = nodes.new('ShaderNodeSeparateRGB')
        #links.new(seperate_rgb.inputs[0], texTopo.outputs['Color'])
        #links.new(bump_shader.inputs['Strength'], seperate_rgb.outputs['G'])
        #mix_shader = nodes.new(type='ShaderNodeMixShader')
        #links.new(bump_shader.inputs['Normal'], texTopo.outputs['Color'])
        #links.new(mix_shader.inputs[0], bump_shader.outputs['Normal'] )
        #links.new(mix_shader.inputs[1], bump_shader.outputs['Normal'])
        #links.new(shader_mat.inputs['Displacement'], mix_shader.outputs[0])
        # displacement
        #links.new(bump_shader.inputs['Normal'], texTopo.outputs['Color'])
        #links.new(shader_mat.inputs['Displacement'], bump_shader.outputs['Normal'])
        #links.new(shader_mat.inputs['Displacement'], bump_shader.outputs[0])

        bsdf = nodes["Principled BSDF"]
        links.new(bsdf.inputs['Normal'], bump_shader.outputs['Normal'])


    def add_sphere_night(self,nodes,links):
        """
        adds night lights texture
        """
        # checks if anything to do
        if not self.texture_night:
            print("no night")
            return

        # user output
        if self.verbose:
            print("night:")
            print("  loading texture: ",self.texture_night)

        # checks if file exists
        if not os.path.isfile(self.texture_night):
            print("Please check if texture file exists: ",self.texture_night)
            sys.exit(1)

        texNight = nodes.new('ShaderNodeTexImage')
        texNight.image = bpy.data.images.load(self.texture_night)

        # main shader
        bsdf = nodes["Principled BSDF"]

        # ramp for creating alpha to only have main lights
        ramp = nodes.new('ShaderNodeValToRGB')
        ramp.color_ramp.elements[0].position = 0.02
        ramp.color_ramp.elements[0].color = [0,0,0,0]
        ramp.color_ramp.elements[1].position = 1.0
        ramp.color_ramp.elements[1].color = [0.95,1.0,0.46,1.0]
        links.new(ramp.inputs['Fac'],texNight.outputs['Color'])

        #night_multiply = nodes.new('ShaderNodeMath')
        #night_multiply.operation = 'MULTIPLY'
        #night_multiply.inputs[1].default_value = 0.8 #0.65

        #gamma = nodes.new('ShaderNodeGamma')
        #gamma.inputs['Gamma'].default_value = 1.5
        #links.new(gamma.inputs['Color'],texNight.outputs['Color'])

        #colormix = nodes.new('ShaderNodeMix')
        #colormix.data_type = 'RGBA'
        #colormix.blend_type = 'OVERLAY'
        #links.new(colormix.inputs['A'],gamma.outputs['Color'])
        #links.new(colormix.inputs['B'],texImage.outputs['Color'])

        #colormix = nodes.new('ShaderNodeMixRGB')
        #colormix.blend_type = 'OVERLAY'
        #links.new(colormix.inputs['Color1'],gamma.outputs['Color'])
        #links.new(colormix.inputs['Color2'],texImage.outputs['Color'])

        # emission node
        emission = nodes.new('ShaderNodeEmission')
        links.new(emission.inputs['Color'],ramp.outputs['Color'])
        emission.inputs['Strength'].default_value = 5.0

        # emission bsdf node (with subsurface lighting for glow effect)
        #emission_bsdf = nodes.new('ShaderNodeBsdfPrincipled')
        #links.new(emission_bsdf.inputs['Base Color'],ramp.outputs['Color'])
        #links.new(emission_bsdf.inputs['Emission'],ramp.outputs['Color'])
        #emission_bsdf.inputs['Emission Strength'].default_value = 20.0
        #links.new(emission_bsdf.inputs['Subsurface Color'],ramp.outputs['Color'])
        #emission_bsdf.inputs['Subsurface'].default_value = 100.0

        #tex_coord = nodes.new(type = 'ShaderNodeTexCoord')

        # mix light emission and main image
        mix = nodes.new('ShaderNodeMixShader')
        mix.name = "Mix Shader"
        links.new(mix.inputs[1], bsdf.outputs[0])
        #links.new(mix.inputs[2], emission_bsdf.outputs[0])  # bsdf glow
        links.new(mix.inputs[2], emission.outputs[0])

        # adds night lights only on dark earth side
        # see: https://github.com/ptabriz/geodesign_with_blender/wiki/advanced_shading#iv-illuminating-the-earth-at-night
        cam = nodes.new('ShaderNodeCameraData')

        map = nodes.new('ShaderNodeMapping')
        map.inputs['Rotation'].default_value = [0 * DEGREE_TO_RAD, 95.0 * DEGREE_TO_RAD, 0 * DEGREE_TO_RAD]
        map.inputs['Scale'].default_value = [10.0, 0.5, 1.0]

        ramp = nodes.new('ShaderNodeValToRGB')
        ramp.color_ramp.interpolation = 'B_SPLINE'
        ramp.color_ramp.elements[0].position = 0.15
        ramp.color_ramp.elements[1].position = 0.32

        links.new(map.inputs['Vector'], cam.outputs['View Vector'])
        links.new(ramp.inputs[0], map.outputs['Vector'])
        links.new(mix.inputs['Fac'], ramp.outputs['Color'])

        # output to material shader
        shader_mat = nodes["Material Output"]
        links.new(shader_mat.inputs['Surface'], mix.outputs[0])


    def create_moviedata_textureImage(self,data_gridded,movdata,counter,power_scaling=0.0,colormap_max=None):
        """
        creates moviedata as texture image
        """
        # adds data as image texture
        # see: https://bertvandenbroucke.netlify.app/assets/code/rotating_sphere.py
        import numpy as np
        #import scipy.interpolate as interpol
        import matplotlib as mpl
        import matplotlib.cm as cm

        print("moviedata: data   shape     = ",np.shape(data_gridded))
        print("moviedata: data   min/max   = ",data_gridded.min(),data_gridded.max())

        # color power scaling
        if power_scaling > 0.0:
            print("moviedata: color power scaling  = ",power_scaling)
            data_gridded = np.where(data_gridded >= 0.0,
                                    data_gridded**power_scaling,
                                    -abs(data_gridded)**power_scaling)

        # limit size
        total_max = abs(data_gridded).max()
        total_max = 0.9 * total_max  # sets limit at 90% of the maximum

        if total_max != 0.0:
            total_max = 1.0 * 10**(int(np.log10(total_max))-1)  # example: 2.73e-11 limits to 1.e-11
            #total_max = 1.0 * 10**(int(np.log10(total_max))-2)  # example: 2.73e-11 limits to 1.e-12
        else:
            total_max = 0.0
            print("moviedata: color data min/max   = ",data_gridded.min(),data_gridded.max())
            print("moviedata: zero color data - nothing to show")
            # nothing left to do
            return

        # checks if fixing maximum value
        if colormap_max:
            total_max = colormap_max

        print("moviedata: color data min/max   = ",data_gridded.min(),data_gridded.max())
        if colormap_max:
            print("moviedata: color data total max = ",total_max," (fixed)")
        else:
            print("moviedata: color data total max = ",total_max)

        # limits range [-total_max,total_max]
        #data_gridded = np.where(data_gridded < -total_max,-total_max,data_gridded)
        #data_gridded = np.where(data_gridded > total_max,total_max,data_gridded)

        # convert the mapped Z values to RGBA colours using a matplotlib colormap
        if movdata.use_component == 4:
            # data norm component range [0,max]
            norm = mpl.colors.Normalize(vmin=0.0, vmax=total_max)
        else:
            # single component: 1==Z,2==N,3==E range [-max,max]
            norm = mpl.colors.Normalize(vmin=-total_max, vmax=total_max)

        mapper = cm.ScalarMappable(norm=norm)
        pixels = mapper.to_rgba(data_gridded).flatten()

        print("moviedata: pixels min/max   = ",pixels.min(),pixels.max())

        # material from globe sphere
        mat = bpy.data.materials["Material Globe"]
        # node-graph
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        # create a texture image
        h = np.shape(data_gridded)[0]  # height
        w = np.shape(data_gridded)[1]  # width

        name = "out.textureimage.{:04d}".format(counter)
        texImg = bpy.data.images.new(name, width=w, height=h)

        # overwrite pixel data
        texImg.pixels[:] = pixels
        # update texture
        texImg.update()

        # save texture as image

        # filepath
        # use absolute path to be able to reload texture as an image sequence
        # w/ relative path
        #filename = name + ".png"
        # w/ absolute path
        dir = os.getcwd()
        filename = dir + "/" + name + ".png"

        texImg.filepath_raw = filename
        texImg.file_format = 'PNG'
        texImg.save()
        print("moviedata: texture image ",filename)
        print("")

        # Pack the image into .blend so it gets saved with it
        texImg.pack()


    def add_moviedata(self,movdata):
        """
        adds moviedata texture
        """
        #print("moviedata: texture images ",bpy.data.images.keys())
        #print("moviedata: nodes ",nodes.keys())

        # material from globe sphere
        mat = bpy.data.materials["Material Globe"]
        # node-graph
        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        # creates texture node with the image
        texData = nodes.new('ShaderNodeTexImage')
        texData.name = "Moviedata Image Texture"

        texData.interpolation = 'Cubic'
        texData.extension = 'EXTEND'

        # note: here, we setup the shader nodes and the corresponding shader links.
        #       for this, we can assign only the first moviedata image moviedata_textureimage.0001
        #       as input image to the shader as we will reload the textures as an image sequence in the output routine
        #       to change the moviedata textures.
        #
        # gets corresponding texture image
        name = "out.textureimage.{:04d}".format(1)
        texImg = bpy.data.images[name]
        texImg.update()

        # sets image
        texData.image = texImg
        texData.update()

        # Connect the Texture Coordinate node to the texture.
        # This uses the active UV map of the object.
        #tex_coord = nodes.new('ShaderNodeTexCoord')
        #tex_coord.from_instancer = True
        #tex_coord.object = bpy.data.objects["Sphere"]
        #tex_coord.texture_coords = 'ORCO'
        #tex_coord.mapping = 'SPHERE'
        #links.new(texData.inputs['Vector'], tex_coord.outputs['UV'])

        # ramp for creating alpha to only have main lights
        ramp = nodes.new('ShaderNodeValToRGB')
        ramp.name = "Moviedata ColorRamp"

        ramp.color_ramp.interpolation = 'B_SPLINE'

        if movdata.use_component == 4:
            # norm: range [0,1]
            ramp.color_ramp.elements[0].position = 0.3
            ramp.color_ramp.elements[0].color = [0,0,0,0.0]
            ramp.color_ramp.elements[1].position = 0.6
            ramp.color_ramp.elements[1].color = [1.0,0.6,0.0,1.0]  # yellow
            # add additional ramp point
            ramp.color_ramp.elements.new(1.0)
            ramp.color_ramp.elements[2].position = 1.0
            ramp.color_ramp.elements[2].color = [1.0,1.0,1.0,1.0]  # white

        else:
            # single component: 1==Z,2==N,3==E
            ramp.color_ramp.elements[0].position = 0.0
            ramp.color_ramp.elements[0].color = [1.0,0.8,0.5,1.0]  # yellow
            ramp.color_ramp.elements[1].position = 0.3
            ramp.color_ramp.elements[1].color = [0.0,0.0,0.0,0.0]
            # add two additional ramp points
            ramp.color_ramp.elements.new(0.7)
            ramp.color_ramp.elements[2].position = 0.7
            ramp.color_ramp.elements[2].color = [0.0,0.0,0.0,0.0]  # white/alpha
            ramp.color_ramp.elements.new(1.0)
            ramp.color_ramp.elements[3].position = 1.0
            ramp.color_ramp.elements[3].color = [1.0,0.8,0.5,1.0]  # yellow

        links.new(ramp.inputs['Fac'],texData.outputs['Color'])

        # emission node
        emission = nodes.new('ShaderNodeEmission')
        emission.name = "Moviedata Emission"

        links.new(emission.inputs['Color'],ramp.outputs['Color'])
        emission.inputs['Strength'].default_value = 5.0

        # emission bsdf node (with subsurface lighting for glow effect)
        #emission_bsdf = nodes.new('ShaderNodeBsdfPrincipled')
        #links.new(emission_bsdf.inputs['Base Color'],ramp.outputs['Color'])
        #links.new(emission_bsdf.inputs['Alpha'],ramp.outputs['Alpha'])
        #links.new(emission_bsdf.inputs['Emission'],ramp.outputs['Color'])
        #emission_bsdf.inputs['Emission Strength'].default_value = 2.0
        #links.new(emission_bsdf.inputs['Subsurface Color'],ramp.outputs['Color'])
        #emission_bsdf.inputs['Subsurface'].default_value = 1.0

        # mix light emission and main image
        mix_2 = nodes.new('ShaderNodeMixShader')
        mix_2.name = "Moviedata Mix Shader"

        links.new(mix_2.inputs['Fac'], ramp.outputs['Alpha'])

        # output from main shader
        if self.texture_night:
            # takes output from night mix shader
            mix_1 = nodes["Mix Shader"]
        else:
            # takes output from main BSDF node
            mix_1 = nodes["Principled BSDF"]
        links.new(mix_2.inputs[1], mix_1.outputs[0])

        # output from moviedata emission node
        links.new(mix_2.inputs[2], emission.outputs[0])
        # or emission bsdf node (glow effect)
        #links.new(mix_2.inputs[2], emission_bsdf.outputs[0])

        # output to material shader
        shader_mat = nodes["Material Output"]
        links.new(shader_mat.inputs['Surface'], mix_2.outputs[0])


    def add_clouds(self):
        """
        adds clouds on second sphere
        """
        # checks if anything to do
        if not self.texture_clouds:
            print("no clouds")
            return

        if self.verbose:
            print("clouds:")
            print("  loading texture: ",self.texture_clouds)

        # checks if file exists
        if not os.path.isfile(self.texture_clouds):
            print("Please check if texture file exists: ",self.texture_clouds)
            sys.exit(1)

        # adds clouds sphere
        # raised a factor 1.001 ~ 7km above.
        bpy.ops.mesh.primitive_uv_sphere_add(enter_editmode=False, align='WORLD',
                                             location=(0, 0, 0), scale=(1.001, 1.001, 1.001))

        # current object
        sphere = bpy.context.object
        sphere.name = "Sphere Cloud"

        # set quality to highest and apply shade smooth.
        bpy.ops.object.modifier_add(type='SUBSURF')
        sphere.modifiers["Subdivision"].quality = 4
        sphere.modifiers["Subdivision"].levels = 4
        sphere.modifiers["Subdivision"].render_levels = 4
        bpy.ops.object.modifier_apply(modifier="Subdivision")
        bpy.ops.object.shade_smooth()

        # simple example material
        #mat = bpy.data.materials.new(name="Material Cloud")
        #mat.diffuse_color = [0.5,0.8,1.0,1.0]
        #mat.specular_intensity = 0.8
        #mat.alpha = 0.0

        # create material
        mat = bpy.data.materials.new(name="Material Cloud")

        # enable transparency for eevee
        if self.blender_engine == 'BLENDER_EEVEE':
            mat.blend_method = 'BLEND'   # 'OPAQUE', 'BLEND', ..

        # enable node-graph edition mode
        mat.use_nodes = True

        nodes = mat.node_tree.nodes
        links = mat.node_tree.links

        texClouds = nodes.new('ShaderNodeTexImage')
        texClouds.image = bpy.data.images.load(self.texture_clouds)

        # adds texture as coloring
        bsdf = nodes["Principled BSDF"]
        bsdf.inputs['Alpha'].default_value = 0.0
        #bsdf.inputs['Transmission'].default_value = 1.0
        #bsdf.inputs['Emission'].default_value = [1,1,1,0.1]
        links.new(bsdf.inputs['Base Color'], texClouds.outputs['Color'])
        links.new(bsdf.inputs['Emission'], texClouds.outputs['Color'])
        bsdf.inputs['Emission Strength'].default_value = 0.1
        # for eevee engine
        if self.blender_engine == 'BLENDER_EEVEE':
            bsdf.inputs['Subsurface'].default_value = 0.2

        # adds a bump effect to make clouds more "volumetric"
        bump_shader = nodes.new('ShaderNodeBump')
        # displacement
        if self.blender_engine == 'BLENDER_EEVEE':
            bump_shader.inputs['Strength'].default_value = 0.02
        else:
            bump_shader.inputs['Strength'].default_value = 1.0
        links.new(bump_shader.inputs['Height'], texClouds.outputs['Color'])

        diff = nodes.new('ShaderNodeBsdfDiffuse')
        links.new(diff.inputs['Color'], texClouds.outputs['Color'])
        links.new(diff.inputs['Normal'], bump_shader.outputs['Normal'])

        # for transparency
        texClouds_gray = nodes.new('ShaderNodeTexImage')
        texClouds_gray.image = bpy.data.images.load(self.texture_clouds)
        # sets colorspace to Non-Color
        texClouds_gray.image.colorspace_settings.name = 'Non-Color'
        #inv = nodes.new('ShaderNodeInvert')
        #inv.inputs[0].default_value = 0
        #links.new(inv.inputs[1], texClouds.outputs['Color'])

        ramp = nodes.new('ShaderNodeValToRGB')
        ramp.color_ramp.interpolation = 'B_SPLINE'
        ramp.color_ramp.elements[0].position = 0.08
        ramp.color_ramp.elements[0].color = [1,1,1,0.0]
        ramp.color_ramp.elements[1].position = 1.0
        ramp.color_ramp.elements[1].color = [1,1,1,1.0]
        links.new(ramp.inputs[0], texClouds_gray.outputs['Color'])

        links.new(bsdf.inputs['Alpha'], ramp.outputs['Alpha'])
        #links.new(bsdf.inputs['Emission'], texClouds.outputs['Color'])

        # mixer
        mix = nodes.new('ShaderNodeMixShader')
        # adds texture as coloring
        links.new(mix.inputs[1], bsdf.outputs[0])
        links.new(mix.inputs[2], diff.outputs[0])

        shader_mat = nodes["Material Output"]
        #links.new(shader_mat.inputs[0], inv.outputs[0])
        links.new(shader_mat.inputs['Surface'], mix.outputs[0])

        # cycles rendering setting
        #mat.cycles.displacement_method = 'BOTH'  # BUMP, DISPLACEMENT, BOTH
        #mat.cycles.use_sss_translucency = True

        # adds material
        sphere.data.materials.append(mat)
        # gets active sphere object
        #obj = bpy.context.view_layer.objects.active
        # adds material to sphere
        #obj.data.materials.append(mat)

    def add_sun(self):
        """
        adds sun light
        """
        # initializes position
        position = (0, 3, 0.5)

        if self.verbose:
            print("sun:",position)

        # adds sun position
        bpy.ops.object.light_add(type='SUN', radius=1, align='WORLD',
                                 location=position,
                                 rotation=(60.0 * DEGREE_TO_RAD, 90.0 * DEGREE_TO_RAD, 150.0 * DEGREE_TO_RAD),
                                 scale=(1, 1, 1))
        # current object
        light = bpy.context.object
        light.name = "Sun"

        # intensity
        light.data.energy = 8
        # angle
        light.data.angle = 0


    def add_moon(self):
        """
        adds moon light
        """
        # initializes position
        position = (1, -1, 1.5)

        if self.verbose:
            print("moon:",position)

        # adds sun position
        bpy.ops.object.light_add(type='AREA', radius=0.53, align='WORLD',
                                 location=position,
                                 rotation=(60.0 * DEGREE_TO_RAD, 20.0 * DEGREE_TO_RAD, 30.0 * DEGREE_TO_RAD),
                                 scale=(1, 1, 1))
        # current object
        light = bpy.context.object
        light.name = "Moon"

        # intensity
        light.data.energy = 5
        # shape
        light.data.shape = 'DISK'


    def add_camera(self):
        """
        adds camera object
        """
        # initializes position
        position = (3, 0, 0.02)

        if self.verbose:
            print("camera:",position)

        # adds camera position
        bpy.ops.object.camera_add(enter_editmode=False, align='VIEW',
                                  location=position,
                                  rotation=(89.6 * DEGREE_TO_RAD, 0.0 * DEGREE_TO_RAD, 90.0 * DEGREE_TO_RAD),
                                  scale=(1, 1, 1))
        # current object
        cam = bpy.context.object
        cam.name = "Camera"


    def add_scene_effects(self):
        """
        adds scene defaults and composite effects
        """
        if self.verbose:
            print("scene:")
            print("  render engine: ",self.blender_engine)
            print("  render device: ",self.blender_device)
            print("  resolution   : (X,Y) = ({},{})".format(self.blender_img_resolution_X,self.blender_img_resolution_Y))
            print("")

        # gets scene
        #scene = bpy.data.scenes["Scene"]
        scene = bpy.context.scene

        # using cycles renderer engine (to have clouds sphere with transparency effect)
        scene.render.engine = self.blender_engine
        scene.cycles.device = self.blender_device

        # turns on bloom
        if scene.render.engine == 'BLENDER_EEVEE':
            scene.eevee.use_bloom = True

        # render resolution
        scene.render.resolution_x = self.blender_img_resolution_X
        scene.render.resolution_y = self.blender_img_resolution_Y

        # debug: print all scene names in a list
        #print("scenes : ",bpy.data.scenes.keys())
        #print("objects: ",bpy.data.objects.keys())
        #for obj in bpy.data.objects:
        #    print("  object: ",obj.name)
        #print("")

        # sets black background color
        bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (0, 0, 0, 1)

        # adds glare effect
        scene.use_nodes = True
        nodes = scene.node_tree.nodes
        links = scene.node_tree.links

        if self.verbose:
            print("  adding glare effect")
            print("  composite nodes: ",nodes.keys())
            print("")

        comp = nodes['Composite']
        render = nodes['Render Layers']

        glare = nodes.new('CompositorNodeGlare')
        glare.glare_type = 'FOG_GLOW'
        glare.quality = 'HIGH'
        glare.threshold = 0.5
        glare.size = 6
        links.new(glare.inputs['Image'], render.outputs['Image'])
        links.new(comp.inputs['Image'], glare.outputs['Image'])


    def setup_image_frame(self,counter,num_datafiles):
        """
        sets corresponding moviedata texture image frame for output rendering
        """
        # checks if anything to do
        if num_datafiles == 0: return

        # material from globe sphere
        mat = bpy.data.materials["Material Globe"]
        # node-graph
        nodes = mat.node_tree.nodes

        # sets corresponding moviedata texture image
        # note: we can either just use the stored images in the bpy.data.images
        #       or reload the textures as an image sequence.
        #       since we will use image sequences for the animation part,
        #       we also do it here for the still image rendering

        #debug
        #print("debug: bpy.data.images ",bpy.data.images.keys())

        # gets existing texture node
        texData = nodes["Moviedata Image Texture"]

        # gets corresponding texture image
        #name = "moviedata_textureimage.{:04d}".format(counter)
        #texImg = bpy.data.images[name]

        # sets image
        #texData.image = texImg
        #texData.update()

        # to animate texture as an image sequence
        if counter == 1:
            # Point to the first Element
            # gets corresponding texture image
            name = "out.textureimage.{:04d}".format(1)
            # reload
            texImg = bpy.data.images[name]
            texImg.reload()
            # same as:
            #dir = os.getcwd()
            #filename = dir + "/" + name + ".png"
            #texImg = bpy.data.images.load(filename)

            # set to image sequence
            texImg.source = "SEQUENCE"

            # don't run texImg.update(), it will say:
            #   RuntimeError: Error: Image 'moviedata_textureimage.0001.png' does not have any image data

            # sets image
            texData.image = texImg

            # sets image sequence
            texData.image_user.use_auto_refresh = True
            texData.image_user.frame_duration = num_datafiles
            texData.image_user.frame_start = 1

            # current offset
            # filepath naming will be using: **.0001.png for start==1 and offset==0
            # filepath naming will be using: **.0002.png for start==1 and offset==1
            # filepath naming will be using: **.0003.png for start==1 and offset==2
            # ..
            texData.image_user.frame_offset = 0

        else:
            # updates frame offset
            # filepath naming will be using: **.0001.png for start==1 and offset==0
            # filepath naming will be using: **.0002.png for start==1 and offset==1
            # filepath naming will be using: **.0003.png for start==1 and offset==2
            # ..
            texData.image_user.frame_offset = (counter-1)

        #print("debug: frame ",texData.image_user.frame_current)
        #print("debug: filepath ",texData.image.filepath_from_user(image_user=texData.image_user))


    def output_image(self,dir,fov=50.0,appendix="",counter=None,num_datafiles=None,suppress=False):
        """
        renders a jpeg image
        """
        if self.verbose:
            print("image:")
            print("  render engine: ",self.blender_engine)
            print("  render device: ",self.blender_device)
            print("")

        # timing
        tic = time.perf_counter()

        # gets scene
        scene = bpy.context.scene

        # link camera
        cam = bpy.data.objects["Camera"]
        scene.camera = cam
        # Set camera fov in degrees
        scene.camera.data.angle = float(fov * DEGREE_TO_RAD)

        # Set camera rotation in euler angles
        #scene.camera.rotation_mode = 'XYZ'
        #scene.camera.rotation_euler[0] = 0.0 * DEGREE_TO_RAD
        #scene.camera.rotation_euler[1] = 0.0 * DEGREE_TO_RAD
        #scene.camera.rotation_euler[2] = -30.0 * DEGREE_TO_RAD

        # Set camera translation
        #scene.camera.location.x = 0.0
        #scene.camera.location.y = 0.0
        #scene.camera.location.z = 8.0

        # output image
        scene.render.image_settings.file_format = 'JPEG'   # 'PNG'
        name = './out'
        if appendix: name = name + '.' + appendix
        name = name + '.jpg'
        scene.render.filepath = name
        scene.render.use_file_extension = False

        # redirect stdout to null
        print("rendering image: {} ...\n".format(name))
        # to avoid long stdout output by the renderer:
        #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Sun
        #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Camera
        #  ..
        with SuppressStream(sys.stdout,suppress):
            # Render Scene and store the scene
            bpy.ops.render.render(write_still=True)

        # save blend file
        name = 'out.blend'
        filename = dir + "/" + name
        bpy.ops.wm.save_as_mainfile(filepath=filename)
        if self.verbose:
            #print("blend file written to: ",filename)
            print("")

        # timing
        toc = time.perf_counter()
        print("elapsed time for image render is {:0.4f} seconds\n".format(toc - tic))


    def setup_animation_frames(self,num_datafiles):
        """
        setup animation settings
        """
        if self.verbose:
            print("animation:")

        # gets scene
        scene = bpy.context.scene

        # moviedata files
        if num_datafiles > 0:
            # sets number of keyframes for animation to the number of moviedata files (if possible)
            self.animation_number_of_keyframes = num_datafiles

        # animation
        # animation variables
        total_frames = self.animation_number_of_keyframes * self.animation_keyframe_interval

        # checks if anything to do
        if total_frames == 0: return

        if self.verbose:
            print("  number of keyframes    = ",self.animation_number_of_keyframes)
            print("  keyframe interval      = ",self.animation_keyframe_interval)
            print("  total number of frames = ",total_frames)
            print("")

        # define a one hundred frame timeline
        scene.frame_end = total_frames
        scene.frame_start = 0

        # moviedata files
        if num_datafiles > 0:
            # material from globe sphere
            mat = bpy.data.materials["Material Globe"]
            # node-graph
            nodes = mat.node_tree.nodes

            # sets corresponding moviedata texture image
            # note: we can either just use the stored images in the bpy.data.images
            #       or reload the textures as an image sequence.
            #       since we will use image sequences for the animation part,
            #       we also do it here for the still image rendering

            #debug
            #print("debug: bpy.data.images ",bpy.data.images.keys())

            # gets existing texture node
            texData = nodes["Moviedata Image Texture"]

            # to animate texture as an image sequence
            # Point to the first Element
            # gets corresponding texture image
            name = "out.textureimage.{:04d}".format(1)
            # reload
            texImg = bpy.data.images[name]
            texImg.reload()
            # same as:
            #dir = os.getcwd()
            #filename = dir + "/" + name + ".png"
            #texImg = bpy.data.images.load(filename)

            # set to image sequence
            texImg.source = "SEQUENCE"

            # don't run texImg.update(), it will say:
            #   RuntimeError: Error: Image 'moviedata_textureimage.0001.png' does not have any image data

            # sets image
            texData.image = texImg

            # sets image sequence
            texData.image_user.use_auto_refresh = True
            #texData.image_user.frame_duration = self.animation_keyframe_interval
            #texData.image_user.frame_start = 0

            # current offset
            # filepath naming will be using: **.0001.png for start==1 and offset==0
            # filepath naming will be using: **.0002.png for start==1 and offset==1
            # filepath naming will be using: **.0003.png for start==1 and offset==2
            # ..
            #texData.image_user.frame_offset = 1

        # add keyframes
        # rotates sphere object, keeps light & camera fixed
        counter = 0
        for frame in range(0, total_frames + 1, self.animation_keyframe_interval):
            scene.frame_set(frame)
            # rotates globe sphere
            sphere = bpy.data.objects['Sphere']
            sphere.rotation_euler[2] = frame * self.animation_rotation_increment
            sphere.keyframe_insert(data_path='rotation_euler')

            # rotates clouds sphere
            if self.texture_clouds:
                clouds = bpy.data.objects['Sphere Cloud']
                clouds.rotation_euler[2] = frame * self.animation_rotation_increment
                clouds.keyframe_insert(data_path='rotation_euler')

            # changes texture
            if num_datafiles > 0:
                # texture counter
                counter += 1
                # current offset
                # filepath naming will be using: **.0001.png for start==1 and offset==0
                # filepath naming will be using: **.0002.png for start==1 and offset==1
                # filepath naming will be using: **.0003.png for start==1 and offset==2
                # ..
                texData.image_user.frame_duration = 1
                texData.image_user.frame_start = 0
                # image number offset
                texData.image_user.frame_offset = counter - 1
                # last keyframe uses last texture again
                if frame == total_frames:
                    texData.image_user.frame_offset = num_datafiles - 1

                texData.image_user.keyframe_insert("frame_offset", frame=frame)

                # sets constant interpolation (instead of linear default)
                if frame == total_frames:
                    #print("debug: actions ",bpy.data.actions.keys())
                    action = bpy.data.actions['Shader NodetreeAction']
                    # shader has only a single f-curve object
                    fcurve = action.fcurves[0]
                    #debug
                    #for fc in action.fcurves:
                    #    print("debug: actions f-curve data_path: ",fc.data_path)
                    # sets interpolation type on keyframe points
                    for point in fcurve.keyframe_points:
                        point.interpolation = 'CONSTANT'

                # or
                # adds an additional frame before next keyframe;
                # avoids having texture image offset increased by intermediate evaluations using the default, linear interpolation
                #if counter < num_datafiles:
                #    texData.image_user.keyframe_insert("frame_offset", frame=frame+(self.animation_keyframe_interval-1))


    def output_animation(self,dir,fov=50.0,appendix="",suppress=False):
        """
        renders animation
        """
        # timing
        tic = time.perf_counter()

        # gets scene
        scene = bpy.context.scene

        # link camera
        cam = bpy.data.objects["Camera"]
        scene.camera = cam
        # Set camera fov in degrees
        scene.camera.data.angle = float(fov * DEGREE_TO_RAD)

        if self.verbose:
            print("animation:")
            print("  field of view          = ",fov)
            print("  total number of frames = ",scene.frame_end - scene.frame_start)
            print("")

        # render settings
        # see: https://docs.blender.org/api/current/bpy.types.FFmpegSettings.html#bpy.types.FFmpegSettings
        # ffmpeg
        scene.render.image_settings.file_format = 'FFMPEG'   # 'AVI_JPEG', 'FFMPEG'
        # codec
        scene.render.ffmpeg.codec = 'H264'                  # compression: 'MPEG4', 'H264'
        #scene.render.ffmpeg.constant_rate_factor = 'HIGH'
        scene.render.ffmpeg.format = 'MPEG4'                # output file container format
        # frames per second (blender default is 24)
        scene.render.fps = 20

        # output movie
        name = './out.anim'
        if appendix: name = name + '.' + appendix
        name = name + '.mp4'
        scene.render.filepath = name
        scene.render.use_file_extension = False

        if self.verbose:
            print("  movie fps                    = ",scene.render.fps)
            print("  movie format                 = ",scene.render.ffmpeg.format )
            print("")

        # redirect stdout to null
        print("rendering animation: {} ...\n".format(name))
        # to avoid long stdout output by the renderer:
        #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Sun
        #    Fra:1 Mem:189.14M (Peak 190.26M) | Time:00:00.68 | Syncing Camera
        #  ..
        with SuppressStream(sys.stdout,suppress):
            # render animation
            bpy.ops.render.render(animation=True)

        # save blend file
        name = 'out.blend'
        filename = dir + "/" + name
        bpy.ops.wm.save_as_mainfile(filepath=filename)
        if self.verbose:
            #print("blend file written to: ",filename)
            print("")

        # timing
        toc = time.perf_counter()
        if toc - tic < 100.0:
            print("elapsed time for animation render is {:0.4f} seconds\n".format(toc - tic))
        else:
            min = int((toc-tic)/60.0)
            sec = (toc - tic) - min * 60
            print("elapsed time for animation render is {} min {:0.4f} sec\n".format(min,sec))

