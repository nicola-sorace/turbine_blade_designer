#!/usr/bin/python
import csv

####

length = 1  # meters
segments = 20

def chord(x):
    return 0.2-x*0.15

def twist(x):
    return x*10 + 5

####

import subprocess as sp

def typ(s):
    sp.call(['xdotool', 'type', str(s)])

def tab(n = 1):
    for i in range(n):
        sp.call(['xdotool', 'key', 'Tab'])

def down(n = 1):
    for i in range(n):
        sp.call(['xdotool', 'key', 'Down'])

#ps = [0, 0.1, 0.2, 0.3, 0.4]
#cs = [0.2, 0.1, 0.05, 0.01, 0.01]
#ts = [0, 1, 2, 3, 4]

""" Test values
ps = [x*(length/segments) for x in range(segments)] #list(range(0, length, length/segments))
cs = list(map(chord, ps))
ts = list(map(twist, ps))
"""

ps = [] #Positions
cs = [] #Chord lengths
ts = [] #Twist angles
fs = [] #Airfoil profile index

with open('blade_qblade.csv') as csv_file:
    csv_reader = csv.reader(csv_file)
    for row in csv_reader:
        ps.append(row[0])
        cs.append(row[1])
        ts.append(row[2])
        fs.append(row[3])


print('Click on \'Insert after section 1\' button (or press escape)')
brand_new = sp.call(['slop'], stdout=sp.PIPE)
if brand_new == 0:
    print('working...')
    for i in range(len(ps)-2):
        sp.call(['xdotool', 'click', '1'])
else:
    print('Updating old blade. Skipping section creation and airfoil selection steps.')

print('Great! Now click on first position.')
sp.check_call(['slop'], stdout=sp.PIPE)
print('working...')
sp.call(['xdotool', 'click', '1 --repeat 2'])

for i in range(len(ps)):
    if i != 0:
        typ(ps[i])
    tab()
    typ(cs[i])
    tab()
    typ(ts[i])
    tab()
    down(int(fs[i]))
    tab(2)

print('All done!')
