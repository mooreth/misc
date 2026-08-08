from samson import *
import math



#some variables
b = 1.567 #this is the C-C bond length in lonsdaleite in angstrom
a = 2.714 #this is center to center spacing of two rings
R = a / math.sqrt(3) #this is the
s = 30 #number of rings in a side edge
d = 59 #number of rings in the diameter

theta0 = 30

atoms = {}

if d != s*2 - 1:
    print("this is not a regular hexagon")
    
for q in range(-(s-1), s):
    rmin = max(-(s-1), -q-(s-1))
    rmax = min( (s-1), -q+(s-1))

    for r in range(rmin, rmax+1):
        
        x = a * (q + r/2)
        y = a * (math.sqrt(3)/2) * r      
        print(q,r, x, y)
        
        cx = a * (q + r/2)
        cy = a * (math.sqrt(3)/2) * r        
        for i in range(6):
        
            theta = math.radians(theta0 + i * 60)
        
            x = cx + R * math.cos(theta)
            y = cy + R * math.sin(theta)
        
            if i % 2 == 0:
                z = 0.273
            else:
                z = -0.273   
            key = (
                round(x, 5),
                round(y, 5),
                round(z, 5)
            )
            
            atoms[key] = (x, y, z)
            


# --------------------------------
# Calculate centroid
# --------------------------------

n = len(atoms)

centroid_x = sum(p[0] for p in atoms.values()) / n
centroid_y = sum(p[1] for p in atoms.values()) / n
centroid_z = sum(p[2] for p in atoms.values()) / n

print("Centroid:")
print(centroid_x, centroid_y, centroid_z)


#create SAMSON strucutral model
sam_atoms = []

# --------------------------------
# Create structural model
# --------------------------------

sm = SBStructuralModel()


# --------------------------------
# Layer 1 - A
# --------------------------------

for x, y, z in atoms.values():

    atom = SBAtom(SBElement.Carbon)

    atom.setPosition(
        SBPosition3(
            SBQuantity.angstrom(x),
            SBQuantity.angstrom(y),
            SBQuantity.angstrom(z)
        )
    )

    sm.addChild(atom)


# --------------------------------
# Layer 2 - B
# Rotate 60 degrees and shift -Z
# --------------------------------

layer2 = {}

angle = math.radians(60)

cos_a = math.cos(angle)
sin_a = math.sin(angle)

for x, y, z in atoms.values():

    # Move centroid to origin
    x0 = x - centroid_x
    y0 = y - centroid_y

    # Rotate 60 degrees
    x_rot = x0 * cos_a - y0 * sin_a
    y_rot = x0 * sin_a + y0 * cos_a

    # Move centroid back
    x_new = x_rot + centroid_x
    y_new = y_rot + centroid_y

    # Shift down one bond length
    z_new = z - b

    layer2[
        (round(x_new, 5),
         round(y_new, 5),
         round(z_new, 5))
    ] = (x_new, y_new, z_new)


for x, y, z in layer2.values():

    atom = SBAtom(SBElement.Carbon)

    atom.setPosition(
        SBPosition3(
            SBQuantity.angstrom(x),
            SBQuantity.angstrom(y),
            SBQuantity.angstrom(z)
        )
    )

    sm.addChild(atom)


# --------------------------------
# Layer 3 - A
# Same XY as layer 1
# Shift down 2*b
# --------------------------------

layer3 = {}

for x, y, z in atoms.values():

    x_new = x
    y_new = y
    z_new = z - 2 * b

    layer3[
        (round(x_new, 5),
         round(y_new, 5),
         round(z_new, 5))
    ] = (x_new, y_new, z_new)


for x, y, z in layer3.values():

    atom = SBAtom(SBElement.Carbon)

    atom.setPosition(
        SBPosition3(
            SBQuantity.angstrom(x),
            SBQuantity.angstrom(y),
            SBQuantity.angstrom(z)
        )
    )

    sm.addChild(atom)


# --------------------------------
# Add complete 3-layer model
# --------------------------------

with SAMSON.holding('Add lonsdaleite hexagon'):

    SAMSON.hold(sm)

    sm.create()

SAMSON.getActiveDocument().addChild(sm)


print("Layer 1 atoms =", len(atoms))
print("Layer 2 atoms =", len(layer2))
print("Layer 3 atoms =", len(layer3))
print("Total atoms   =", len(atoms) + len(layer2) + len(layer3))