#
# Reverse Engineering 3D Printer Files to CAD models.
#
# This takes a set file created by ply2set.py and
# outputs an OpenSCAD CAD model from it.
#
# RepRap Ltd
# https://reprapltd.com
#
# Written by: Adrian Bowyer
#
# 3 September 2021
#
# Licence: GPL
#

import math as maths
import numpy as np


# Fudge factor
small = 0.000001

# This takes a halfspace and makes an OpenSCAD module() for it.
# OpenSCAD doesn't deal in halfSpaces, so we use a giant cube with one of
# its faces in the right place.  The module name is sNNN where NNN is the
# index of the half space used in the set-theoretic expression that describes the object.
def OpenSCADHalfSpace(halfSpaceAndIndex):
 global negativeCorner, positiveCorner, diagonal, halfSpaceList, set
 hsIndex = halfSpaceAndIndex[0]
 halfSpace = halfSpaceAndIndex[1]
 hsString = "module s"+str(round(hsIndex)) + "()\n{\n"

 # The final translation. At this stage we have our cube with its relevant face going through
 # the origin and the normal of that face pointing in the right direction. All we have to do
 # is translate it along that normal vector by it's distance from the origin.
 translate = np.multiply(halfSpace[0], -halfSpace[1])
 translate = " translate([" + str(translate[0]) + ", " + str(translate[1]) + ", " + str(translate[2]) +"])\n"
 hsString += translate;

 # Our cube will have its top face in the XY plane, with the normal of that face
 # pointing up the Z axis. The outer product of that normal and the halfspace's
 # normal will be the axis vector we need to rotate the cube around, and the magnitude
 # of that vector will be the sine of the angle we need to rotate by.
 cross = np.cross(halfSpace[0], [0, 0, 1])
 sine = maths.sqrt(np.dot(cross, cross))
 # Stop small rounding errors messing with the asin() function.
 if sine > 1:
  sine = 1
 if sine < -1:
  sine = -1

 # Create the cube and shift it down so it's top face is in the XY plane.
 # Test if we actually need to rotate it or not.
 d = -5.0*diagonal
 if abs(sine) < small:
  cosine = np.dot(halfSpace[0], [0, 0, 1])
  if cosine > 0:
   rotateAndTranslate = " translate([" + str(d) + ", " + str(d) + ", " + str(2.0*d) +"])\n"
  else:
   rotateAndTranslate = " translate([" + str(d) + ", " + str(d) + ", 0.0])\n"
 else:
  sine = maths.asin(sine)
  rotateAndTranslate = " rotate(a = " + str(-180.0*sine/maths.pi) + ", v = [" + str(cross[0]) + ", " + str(cross[1]) + ", " + str(cross[2]) + "])\n"
  rotateAndTranslate += " translate([" + str(d) + ", " + str(d) + ", " + str(2.0*d) +"])\n"

 # Finally make the cube, which all of the above is applied to.
 d = 10.0*diagonal
 hsString += rotateAndTranslate + " cube([" + str(d) + ", " + str(d) + ", " + str(d) +"]);\n}\n\n"
 return hsString


# Take a set theoretic expression like 7 & 12 (where 7 and 12 are half-space
# indices) and create the OpenSCAD expression that does that with module s7()
# and so on.
def SetExpression(o, a, b):
 if o == '|':
  result = "union()"
 elif o == '&':
  result = "intersection()"
 else:
  result = "difference()"
 result += "{"
 if a[0].isdigit():
  result += "s" + str(a) + "();"
 else:
  result += a + ";"
 if b[0].isdigit():
  result += "s" + str(b) + "();"
 else:
  result += b + ";"
 result += "}"
 return result

# Is o a set theory operator?
def IsOperator(o):
 if o == '~':
  return True
 if o == '|':
  return True
 if o == '&':
  return True
 return False

# Parse the reverse polish set-theoretic expression creating its OpenSCAD equivalent.
def ParseSet():
 global negativeCorner, positiveCorner, diagonal, halfSpaceList, set
 a=[]
 for s in set:
  if IsOperator(s):
   a.append(SetExpression(s, a.pop(), a.pop()))
  else:
   a.append(s)
 set = a[0]

# Small function to turn things like "[0.0, -0.28, 0.77, -15.9]" from a file
# and turn them into a Python list. There's probably a library function to do this
# if I could be bothered to look, but it's quicker to write it myself...
def GetListFromLine(file):
 lst = file.readline()
 lst = lst.replace("[", "")
 lst = lst.replace("]", "")
 lst = " ".join( lst.split() )
 lst = lst.split(" ")
 result = []
 for item in lst:
  result.append(float(item))
 return result

# Read the set file from the ply2set.py program and from it
# create a list of OpenSCAD halfspace module()s and the OpenSCAD
# set theory to tie them together.
def ReadSetFile(fileName):
 global negativeCorner, positiveCorner, diagonal, halfSpaceList, set
 file = open(fileName, "r")
 negativeCorner = GetListFromLine(file)
 positiveCorner = GetListFromLine(file)
 diagonal = np.subtract(positiveCorner, negativeCorner)
 diagonal = maths.sqrt(np.dot(diagonal, diagonal))
 halfSpaceCount = round(GetListFromLine(file)[0])
 halfSpaceList = []
 for hsIndex in range(halfSpaceCount):
  hs = GetListFromLine(file)
  normal = [hs[1], hs[2], hs[3]]
  d = hs[4]
  hs = (hs[0], [normal, d])
  halfSpaceList.append(hs)
 set = file.readline()
 set = set.replace("\n", "")
 set = set.split(" ")
 newList = []
 for hsIndex in range(len(halfSpaceList)):
  newList.append(OpenSCADHalfSpace(halfSpaceList[hsIndex]))
 halfSpaceList = newList
 ParseSet()



ReadSetFile("STL2CSG-test-objects-woo-2.set")
for hs in halfSpaceList:
 print(hs)
print(set)


