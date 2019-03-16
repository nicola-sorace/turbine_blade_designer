import bpy
import bmesh
import mathutils
import math

import csv


def make_loop(profile, p, c, t):
    edges = []
    verts = []
    v_old = None
    
    for i in range(res): # Iterate through points in profile
        v = bm.verts.new( (profile[i][0], profile[i][1], p) )
        verts.append(v)
        
        if v_old != None:
            edges.append(bm.edges.new( (v,v_old) ))
        else:
            v_first = v
        
        v_old = v
    
    edges.append(bm.edges.new( (v_first,v) ))  # Close the loop
    bmesh.ops.scale(bm, vec=[c,-c,1], verts=verts)
    bmesh.ops.rotate(bm, matrix=mathutils.Matrix.Rotation(math.radians(-t+90),3,'Z'), cent=(0,0,0), verts=verts)
    
    return( (edges, verts) )

def get_profile(name):
    profile = []
    with open(bpy.path.abspath("//profiles/"+name+".clean"), 'r') as lines:
        for line in lines:
            x,y = line.split(' ')
            profile.append( (float(x),float(y)) )
    
    # Offset points so max chord length is center:
    #ys = [y for (x,y) in profile]
    #center_offset = profile[ys.index(max(ys))][0]
    center_offset = 0.30 #0.225
    profile = [(x-center_offset,y) for (x,y) in profile]
    return profile

def lerp(a, b, f):
    return a*(1-f) + b*f

def interp_foil(FS, f):
    Fa = FS[ int(math.floor(f)) ]
    Fb = FS[ int(math.ceil(f)) ]
    f = f%1

    return [ (lerp(x1,x2,f),lerp(y1,y2,f)) for ((x1,y1),(x2,y2)) in zip(Fa, Fb) ]


### Setup mesh editing:
obj = bpy.data.objects['blade']

scene = bpy.context.scene
scene.objects.active = obj
obj.select = True

mesh = bpy.context.object.data
bm = bmesh.new()


### Load airfoil profiles from data files:
res = len(get_profile('naca0018_hires')) # Number of points in profile (resolution) ####

r = 0.5
circle = [(r*math.cos(2*math.pi*t/res),r*math.sin(2*math.pi*t/res)) for t in range(res)]

FOILS = [circle, get_profile('naca0018_hires'), get_profile('naca4418_hires'), get_profile('naca4418_hires')]

### Load blade paramters:
ps = []
cs = []
ts = []
fs = []

with open(bpy.path.abspath("//matlab/blade_blend.csv")) as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        ps.append(row[0])
        cs.append(row[1])
        ts.append(row[2])
        fs.append(row[3])

### Generate blade mesh:
edges_old = None
for n in range(0, len(ps)): # Iterate through sections
    
    profile = interp_foil(FOILS, float(fs[n])-1)
    
    (edges, verts) = make_loop(profile, float(ps[n]), float(cs[n]), float(ts[n]))
    
    
    if edges_old != None:
        bmesh.ops.bridge_loops(bm, edges=edges_old+edges)
    else:
        base_verts = verts
    edges_old = edges
    
    if n==len(ps)-1:
        bm.faces.new(verts) # Close top

### Clean up:
#bmesh.ops.remove_doubles(bm, verts=bm.verts, dist=0.000001)
bmesh.ops.recalc_face_normals(bm, faces=bm.faces)
#bmesh.ops.normals_make_consistent(bm, inside=False) #TODO make work


### Finalize mesh editing:
bm.to_mesh(mesh)
bm.free()