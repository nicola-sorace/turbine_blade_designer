import csv
import Part,PartGui
import math

from FreeCAD import Base

DIR = "/home/nicola/Documents/uni/windsOfChange/"

def get_polyshape(points):
	poly = Part.makePolygon(points)
	shape = doc.addObject("Part::Feature")
	shape.Shape = poly
	return shape

def get_profile(name):
	profile = []
	with open(DIR+"profiles/"+name+".clean", 'r') as lines:
	    for line in lines:
	        x,y = line.split(' ')
	        profile.append( (float(x),float(y)) )

	center_offset = 0.30
	profile = [(x-center_offset,y) for (x,y) in profile]
	return profile

def lerp(a, b, f):
	return a*(1-f) + b*f

def interp_foil(FS, f):
	Fa = FS[ int(math.floor(f)) ]
	Fb = FS[ int(math.ceil(f)) ]
	f = f%1
	return [ (lerp(x1,x2,f),lerp(y1,y2,f)) for ((x1,y1),(x2,y2)) in zip(Fa, Fb) ]

res = len(get_profile('naca0018_hires'))
r = 0.5                                                                                                      
circle = [(r*math.cos(2*math.pi*t/res),r*math.sin(2*math.pi*t/res)) for t in range(res)]
FOILS = [circle, get_profile('naca0018_hires'), get_profile('naca4418_hires'), get_profile('naca4418_hires')]

### Load blade paramters:
ps = []
cs = []
ts = []
fs = []

with open(DIR+"matlab/blade_blend.csv", 'r') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        ps.append( float(row[0]) )
        cs.append( float(row[1]) )
        ts.append( float(row[2]) )
        fs.append( float(row[3]) )

### Generate blade:
doc = App.newDocument()

loft = doc.addObject("Part::Loft", "Blade")

sections = []
for n in range(len(ps)):
	if n in range(9,17):
		continue
	profile = interp_foil(FOILS, fs[n]-1)
	profile.append(profile[0]) # Close the loop
	
	profile = [ (x*cs[n], y*cs[n]) for (x,y) in profile ] # Scale
	a = math.radians( -ts[n] + 90 )
	profile = [ (x*math.cos(a)-y*math.sin(a), x*math.sin(a)+y*math.cos(a) ) for (x,y) in profile ]
	
	scale = 1000 # Meters to millimeters
	points = [ (x*scale, y*scale, ps[n]*scale) for (x,y) in profile ]
	shape = get_polyshape(points)
	
	sections.append(shape)

loft.Sections = sections
loft.Solid = True
loft.Ruled = True

doc.recompute()