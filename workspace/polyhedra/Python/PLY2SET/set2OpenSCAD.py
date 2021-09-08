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
import numbers

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


def SetExpression(op, a, b):
 result = op + "{"
 if isinstance(a, (int)):
  result += "s" + str(a) + "();"
 else:
  result += a + ";"
 if isinstance(b, (int)):
  result += "s" + str(b) + "();"
 else:
  result += b + ";"
 result += "}"
 return result

def IsOperator(o):
 if o == '-':
  return True
 if o == 'u':
  return True
 if o == '^':
  return True
 return False

def ParseSet():
 global negativeCorner, positiveCorner, diagonal, halfSpaceList, set
 newSet = []
 for s in set:
  if IsOperator(s):
   newSet.append(s)
  else:
   newSet.append(int(s))
 set = newSet
 a=[]
 b={'-': lambda x,y: SetExpression("difference()", x, y),
    'u': lambda x,y: SetExpression("union()", x, y),
    '^': lambda x,y: SetExpression("intersection()", x, y)}
 for c in set:
  if c in b:
   a.append(b[c](a.pop(),a.pop()))
  else:
   a.append(c)
 set = a[0]



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
 ParseSet()



ReadSetFile("woo2.set")
for hs in halfSpaceList:
 print(hs)
print(set)


