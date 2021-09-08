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

def OpenSCADHalfSpace(hsIndex, halfSpace):
 hsString = "s"+str(hsIndex) + " = "
 d = str(-5.0*diagonal)
 hsString += "translate([" + d + ", " + d + ", 0])\n"
 d = str(10.0*diagonal)
 hsString += " cube([" + d + ", " + d + ", " + d +"]);\n"
 print(hsString)
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
 newList = []
 for hsIndex in range(len(halfSpaceList)):
  newList.append(OpenSCADHalfSpace(hsIndex, halfSpaceList[hsIndex]))
 halfSpaceList = newList



ReadSetFile("woo2.set")


