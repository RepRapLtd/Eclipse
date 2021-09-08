#
# Reverse Engineering 3D Printer Files to CAD models.
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

def OpenSCADHalfSpace(hsIndex, halfSpace):
 hsString = "module s"+str(hsIndex) + "()\n{\n"

 translate = np.multiply(halfSpace[0], halfSpace[1])
 translate = " translate([" + str(translate[0]) + ", " + str(translate[1]) + ", " + str(translate[2]) +"])\n"
 hsString += translate;

 cross = np.cross(halfSpace[0], [0, 0, 1])
 sine = maths.sqrt(np.dot(cross, cross))

 d = -5.0*diagonal
 if abs(sine) < small:
  cosine = np.dot(halfSpace[0], [0, 0, 1])
  if cosine > 0:
   rotateAndTranslate = " translate([" + str(d) + ", " + str(d) + ", " + str(2.0*d) +"])\n"
  else:
   rotateAndTranslate = " translate([" + str(d) + ", " + str(d) + ", 0.0])\n"
 else:
  sine = maths.sin(sine)
  rotateAndTranslate = " rotate(a = " + str(180.0*sine/maths.pi) + ", v = [" + str(cross[0]) + ", " + str(cross[1]) + ", " + str(cross[2]) + "])\n"
  rotateAndTranslate += " translate([" + str(d) + ", " + str(d) + ", " + str(2.0*d) +"])\n"

 d = 10.0*diagonal
 hsString += rotateAndTranslate + " cube([" + str(d) + ", " + str(d) + ", " + str(d) +"]);\n}\n\n"
 return hsString

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

def IsOperator(symbol):
 if symbol == 'u':
  return True
 if symbol == '^':
  return True
 if symbol == '-':
  return True
 return False

def ParseSet():
 global negativeCorner, positiveCorner, diagonal, halfSpaceList, set
 a=[]
 b={'-': lambda x,y: "difference() {s" + str(x) + "(); s" + str(y)+ "();}",
    'u': lambda x,y: "union() {s" + str(x) + "(); s" + str(y)+ "();}",
    '^': lambda x,y: "intersection() {s" + str(x) + "(); s" + str(y)+ "();}"}
 for c in set:
  if c in b:
   a.append(b[c](a.pop(),a.pop()))
  else:
   a.append(c)
 s = a[0]
 s = s.replace("sintersection", "intersection")
 s = s.replace("sunion", "union")
 s = s.replace("sdifference", "difference")
 s = s.replace("}()", "}")
 return s



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
  normal = [hs[0], hs[1], hs[2]]
  d = hs[3]
  hs = [normal, d]
  halfSpaceList.append(hs)
 set = file.readline()
 set = set.replace("\n", "")
 set = set.split(" ")
 newList = []
 for hsIndex in range(len(halfSpaceList)):
  newList.append(OpenSCADHalfSpace(hsIndex, halfSpaceList[hsIndex]))
 halfSpaceList = newList
 set = ParseSet()



ReadSetFile("woo2.set")
for hs in halfSpaceList:
 print(hs)
print(set)


