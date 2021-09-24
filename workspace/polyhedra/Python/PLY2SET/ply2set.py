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
# ply2set.py is a program to convert PLY triangulation files to set-theoretic or CSG solid
# models. There is a companion program - set2OpenSCAD.py that converts the results
# to openSCAD models, which can be used in OpenSCAD (clearly) and also FreeCAD.
#
# PLY files are easily created from STL files, which are universally used as input
# to 3D printers. For example the Java program Triangulation.java in this
# repository will do that conversion. Thus this program allows the reverse engineering
# of a 3D Printer file back to a CAD file.
#
# The program uses a modification of Tony Woo's alternating sum of volumes algorithm, which is described
# in his two papers in the Documentation directory of this repository. At the moment it is
# a partial implementation; there are some pathological cases with which it won't deal.
#
# Here are some useful and explanatory web pages:
#
# Ply files: https://en.wikipedia.org/wiki/PLY_(file_format)
# Ply files in Python: https://github.com/dranjan/python-plyfile
# Convex hulls: https://en.wikipedia.org/wiki/Convex_hull
# Convex hulls in Python: https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.ConvexHull.html
# Set-theoretic or CSG solid models: https://en.wikipedia.org/wiki/Constructive_solid_geometry
# Boolean logic in Python: https://booleanpy.readthedocs.io/en/latest/index.html
# OpenSCAD: https://openscad.org/
# FreeCAD: https://www.freecadweb.org/
# STL to PLY conversion free online: https://products.aspose.app/3d/conversion/stl-to-ply
# Infix to postfix/reverse polish conversion: http://csis.pace.edu/~wolf/CS122/infix-postfix.htm
# Planar half space geometric models: http://adrianbowyer.com/Publications/Better-faster-pictures-Bowyer-Woodwark.pdf
#

import copy
import sys
import random
import re
import pyglet
from pyglet.gl import *
from pyglet import window
import math as maths
import numpy as np
import boolean
from plyfile import PlyData
from scipy.spatial import ConvexHull

# Fudge factor
small = 0.000001

# What (if anything) to plot
#
# 0 - no graphics
# 1 - input triangles at each stage
# 2 - input triangles at each stage plus the CH
# 3 - input and output triangles at each stage plus the convex hulls
graphics = 0

# True to print debugging information
debug = False

# The bounding box and a point near the middle
positiveCorner = [-sys.float_info.max, -sys.float_info.max, -sys.float_info.max]
negativeCorner = [sys.float_info.max, sys.float_info.max, sys.float_info.max]
diagonal = None
middle = [0, 0, 0]

# We want there to be just one global instance of these
halfSpaces = None
finalSet = None

# How deep to run the recursion before bailing out because
# we may have a pathological infinite-recursion case
maximumLevel = 10

# The Set class uses this to represent the set theoretic expression
booleanAlgebra = boolean.BooleanAlgebra()
TRUE, FALSE, NOT, AND, OR, symbol = booleanAlgebra.definition()


# *******************************************************************************************************

# Small holding class to load the PLY file and provide simple access to its data.

class TriangleFileData:
    def __init__(self, fileName):
        self.plyFileData = PlyData.read(fileName)
        # Ply files can hold polygons other than triangles. Check...
        for face in self.plyFileData['face']:
            if len(face[0]) != 3:
                error = "\nA polygon in ply file " + fileName + " has " + str(
                    len(face[0])) + " vertices (i.e. it's not a triangle).\n"
                error += "ply2set.py can only work with triangles.\n"
                sys.exit(error)

    def VertexCount(self):
        return len(self.plyFileData['vertex'])

    def Vertex(self, index):
        return list(self.plyFileData['vertex'][index])

    def TriangleCount(self):
        return len(self.plyFileData['face'])

    def Triangle(self, index):
        return self.plyFileData['face'][index][0]


# ************************************************************************************************************

# Class to hold a single triangle. The three vertices are indices into the list of vertices
# that have been previously loaded from a PLY file into ply.

class Triangle:
    def __init__(self, vertices, colour, ply):
        self.vertices = vertices
        self.neighbours = [-1, -1, -1] # Set up later
        self.ply = ply
        self.points = self.Points()
        self.normal = np.cross(np.subtract(self.points[1], self.points[0]), np.subtract(self.points[2], self.points[0]))
        s2 = maths.sqrt(np.dot(self.normal, self.normal))
        if s2 < small:
            print("Triangle: vertices are collinear!")
        self.normal = np.multiply(self.normal, 1.0 / s2)
        self.colour = RandomShade(colour)
        self.centroid = np.multiply(np.add(np.add(self.points[0], self.points[1]), self.points[2]), 1.0 / 3.0)
        # This will be an index into halfSpaces, not the halfspace itself
        self.halfSpace = -1
        self.Reset()

    # The three actual space (x, y, z) coordinates of the vertices
    def Points(self):
        return [self.ply.Vertex(self.vertices[0]), self.ply.Vertex(self.vertices[1]), self.ply.Vertex(self.vertices[2])]

    # The list of vertex indices.
    def Vertices(self):
        return self.vertices

    def SetColour(self, colour):
        self.colour = RandomShade(colour)

    def Reset(self):
        self.flag = False

    # Change the sense of the triangle. The normal should point from solid to air.
    def Invert(self, hspaces):
        temp = self.vertices[0]
        self.vertices[0] = self.vertices[1]
        self.vertices[1] = temp
        temp = self.neighbours[0]
        self.neighbours[0] = self.neighbours[1]
        self.neighbours[1] = temp
        self.halfSpace = hspaces.Get(self.halfSpace).complement
        self.normal = [-self.normal[0], -self.normal[1], -self.normal[2]]
        self.points = self.Points()

    # Is a point inside the triangle?
    # Note - this is doing an essentially 2D calculation in 3D space.
    # To see how this works write it out as vector algebra.
    def ContainsPoint(self, point):
        for p0 in range(3):
            p1 = (p0 + 1) % 3
            p2 = (p0 + 2) % 3
            p10 = np.subtract(self.points[p1], self.points[p0])
            p20 = np.subtract(self.points[p2], self.points[p0])
            d12 = np.cross(p10, p20)
            pp0 = np.subtract(point, self.points[p0])
            dp0 = np.cross(p10, pp0)
            if np.dot(d12, dp0) < 0:
                return False
        return True

    # Is this triangle contiguous with another?
    def MakeNeighbourIfSharesEdge(self, triangle):
        count = 0
        for v in range(len(self.vertices)):
            myV0 = self.vertices[v]
            myV1 = self.vertices[(v + 1) % 3]
            for tv in range(len(triangle.vertices)):
                tV0 = triangle.vertices[tv]
                tV1 = triangle.vertices[(tv + 1) % 3]
                if ((myV0 == tV0) and (myV1 == tV1)) or ((myV1 == tV0) and (myV0 == tV1)):
                    self.neighbours[(v + 2) % 3] = triangle
                    triangle.neighbours[(tv + 2) % 3] = self
                    return

    def __str__(self):
        corners = self.Points()
        result = "<\n"
        for c in corners:
            result += " " + str(c) + "\n"
        result += ">\n"
        return result


def ConstructNeighbourLists(triangles):
    for tr0 in range(len(triangles)):
        t0 = triangles[tr0]
        for tr1 in range(tr0 + 1, len(triangles)):
            t1 = triangles[tr1]
            if t0 is not t1:
                t0.MakeNeighbourIfSharesEdge(t1)


# Take a base RGB colour and perturb it a bit. Used to allow coplanar triangles
# to be distinguished.
def RandomShade(baseColour):
    colour = [random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1), random.uniform(-0.1, 0.1)]
    for c in range(3):
        colour[c] += baseColour[c]
        if colour[c] < 0:
            colour[c] = 0
        if colour[c] > 1:
            colour[c] = 1
    return colour


# *************************************************************************************************************

# A planar half space of the form Ax + By + Cz + D <= 0, where (A, B, C) is the
# plane's normal vector and -D is the distance from it to the origin.
#
# It is the plane through a triangle.

class HalfSpace:
    def __init__(self, triangle=None):
        self.complement = -1
        self.flag = False
        if triangle is None:
            return
        self.normal = triangle.normal
        self.d = -np.dot(self.normal, triangle.centroid)

    # Note: two coincident half spaces of opposite sense are not considered equal.
    def __eq__(self, halfSpace):
        if 1.0 - np.dot(halfSpace.normal, self.normal) > small:
            return False
        if abs(self.d - halfSpace.d) > small:
            return False
        return True

    # Is a halfspace my complement?
    def Opposite(self, halfSpace):
        if 1.0 + np.dot(halfSpace.normal, self.normal) > small:
            return False
        if abs(self.d + halfSpace.d) > small:
            return False
        return True

    # Work out the complement
    def Complement(self):
        result = HalfSpace()
        result.normal = [-self.normal[0], -self.normal[1], -self.normal[2]]
        result.d = -self.d
        result.complement = -1
        return result

    # Used for output so that only halfSpaces that are actually referred to
    # in a set-theoretic expression are saved.
    def SetFlag(self, f):
        self.flag = f

    # Get the coeficients as an array.
    def Coeficients(self):
        return np.array([self.normal[0], self.normal[1], self.normal[2], self.d])

    # The potential value of a point [x, y, z]. -ve means in solid; +ve in air.
    # Thus this is a membership test.
    def Value(self, point):
        return np.dot(self.normal, point) + self.d

    def __str__(self):
        result = "{" + str(self.normal[0]) + "x + " + str(self.normal[1]) + "y + " + \
                 str(self.normal[2]) + "z + " + str(self.d) + " <= 0}"
        result = result.replace("+ -", "- ")  # Sigh!
        return result

    # This takes a halfspace and makes an OpenSCAD module() for it.
    # OpenSCAD doesn't deal in halfSpaces, so we use a giant cube with one of
    # its faces in the right place.  The module name is sNNN where NNN is the
    # index of the half space used in the set-theoretic expression that describes the object.
    def ToOpenSCAD(self, hsIndex):
        hsString = "module s"+str(round(hsIndex)) + "()\n{\n"

    # The final translation. At this stage we have our cube with its relevant face going through
    # the origin and the normal of that face pointing in the right direction. All we have to do
    # is translate it along that normal vector by it's distance from the origin.
        translate = np.multiply(self.normal, -self.d)
        translate = " translate([" + str(translate[0]) + ", " + str(translate[1]) + ", " + str(translate[2]) +"])\n"
        hsString += translate;

    # Our cube will have its top face in the XY plane, with the normal of that face
    # pointing up the Z axis. The outer product of that normal and the halfspace's
    # normal will be the axis vector we need to rotate the cube around, and the magnitude
    # of that vector will be the sine of the angle we need to rotate by.
        cross = np.cross(self.normal, [0, 0, 1])
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
            cosine = np.dot(self.normal, [0, 0, 1])
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

def MakeHalfSpaceFromCoefficients(coefficients):
    result = HalfSpace()
    result.normal = [coefficients[0], coefficients[1], coefficients[2]]
    result.d = coefficients[3]
    return result

# ********************************************************************************************************

# A list of half spaces

class HalfSpaceList:
    def __init__(self):
        self.halfSpaceList = []

    # Is halfSpace in the list? If so return an index to it. If not, return -1
    def LookUp(self, halfSpace):
        for hs in range(len(self.halfSpaceList)):
            listHalfSpace = self.halfSpaceList[hs]
            if halfSpace == listHalfSpace:
                return hs
        return -1

    # Add the halfspace to the list, unless that half space is
    # already in the list. Return the index of the half space in the list.
    # Also add its complement, so we can flip between them.
    def AddSpace(self, halfSpace):
        last = len(self.halfSpaceList)
        inList = self.LookUp(halfSpace)
        if inList < 0:
            self.halfSpaceList.append(halfSpace)
            self.halfSpaceList.append(halfSpace.Complement())
            halfSpace.complement = last + 1
            self.halfSpaceList[last + 1].complement = last
            return last
        else:
            return inList

    # Add the halfspace corresponding to a triangle using the function above.
    def AddTriangle(self, triangle):
        halfSpace = HalfSpace(triangle)
        index = self.AddSpace(halfSpace)
        triangle.halfSpace = index
        return index

    # Set all the flags false.
    def ResetFlags(self):
        for hs in self.halfSpaceList:
            hs.SetFlag(False)

    # Get an actual half space from its index.
    def Get(self, index):
        return self.halfSpaceList[index]

    def __len__(self):
        return len(self.halfSpaceList)


# ********************************************************************************************************

# Take a set theoretic expression like 7 12 & (where 7 and 12 are half-space
# indices) and create the OpenSCAD expression that does that with module s7()
# and so on.
def OpenSCADSetExpression(a, b, o):
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

# Unions and intersections of half-spaces.
# Set() gives the empty set
# Set(-n) gives the universal set

class Set:

    # Create a set, optionally corresponding to a halfspace.
    # If we convert to reverse polish (self.rpExpression) that
    # will be remembered for future lazy evaluation.
    def __init__(self, halfSpaceIndex=None):
        self.rpExpression = None
        if halfSpaceIndex is None:
            self.expression = FALSE
            self.complexity = 0
            return
        if halfSpaceIndex < 0:
            self.expression = TRUE
            self.complexity = 0
            return
        # We either use the halfspace or its complement, whichever is
        # lexically earlier.
        hsComplement = halfSpaces.Get(halfSpaceIndex).complement
        if halfSpaceIndex < hsComplement:
            self.expression = booleanAlgebra.Symbol(halfSpaceIndex)
        else:
            self.expression = NOT(booleanAlgebra.Symbol(hsComplement))
        self.complexity = 1

    # I hope the following are self-explanatory.

    def Unite(self, set2):
        result = Set()
        result.expression = OR(self.expression, set2.expression)
        result.complexity = self.complexity + set2.complexity
        return result

    def Intersect(self, set2):
        result = Set()
        result.expression = AND(self.expression, set2.expression)
        result.complexity = self.complexity + set2.complexity
        return result

    def Subtract(self, set2):
        result = self.Intersect(set2.Complement())
        return result

    def Complement(self):
        self.Simplify()
        self.ToRP()
        stack = []
        for r in self.rpExpression:
            if r[0].isdigit():
                hsIndex = int(r)
                c = halfSpaces.Get(hsIndex).complement
                if c < hsIndex:
                    stack.append(booleanAlgebra.Symbol(c))
                else:
                    stack.append(NOT(booleanAlgebra.Symbol(hsIndex)))
            else:
                b = stack.pop()
                a = stack.pop()
                if r == '&':
                    stack.append(OR(a, b))
                elif r == '|':
                    stack.append(AND(a, b))
                else:
                    print("Set.Complement(): illegal operator: " + r)
        result = Set()
        result.expression = stack[0]
        result.complexity = self.complexity
        return result

    def Simplify(self):
        self.expression = self.expression.simplify()
        self.ToRP()
        self.complexity = 0
        for r in self.rpExpression:
            if r[0].isdigit():
                self.complexity += 1

    # Membership test. If the value of the point is negative, it is in the solid
    # region; if it is positive it is in air.
    def Value(self, point):
        self.Simplify()
        self.ToRP()
        stack = []
        for r in self.rpExpression:
            if r[0].isdigit():
                v = halfSpaces.Get(int(r)).Value(point)
                stack.append(v)
            else:
                b = stack.pop()
                a = stack.pop()
                if r == '&':
                    stack.append(max(a, b))
                elif r == '|':
                    stack.append(min(a, b))
                else:
                    print("Set.Value(): illegal operator: " + r)
        return stack[0]

    def __str__(self):
        return str(self.expression)

    # Convert the set's expression to reverse polish. It remembers the result for future
    # Lazy evaluation. Note this depends on the fact that boolean expressions are
    # constructed in such a way that the NOT/~ operator only applies to single half spaces, not
    # compound expressions.
    def ToRP(self):
        if self.rpExpression is not None:
            return self.rpExpression

        # Start by splitting the expression into unique tags, which are either brackets, operators, or symbols.
        # The symbols are the half-space indices (i.e. positive integers).
        s = str(self.expression)
        s = re.split(r'(\d+)', s)

        # This splits chains like ')&~' into their individual elements
        expression = []
        for e in s:
            if len(e) <= 0:
                pass
            elif len(e) == 1:
                expression.append(e)
            elif e[0].isdigit():
                expression.append(e)
            else:
                for d in e:
                    expression.append(d)

        # Now remove complement operators and replace the halfspaces that they complement
        # with the index of that halfspace's actual complement. Note this is the bit that
        # depends on the fact that boolean expressions are constructed in such a way that
        # the NOT/~ operator only applies to single half spaces, not to compound expressions.
        s = expression
        expression = []
        e = 0
        while e < len(s):
            tag = s[e]
            if tag == '~':
                e += 1
                tag = s[e]
                hs = halfSpaces.Get(int(tag)).complement
                expression.append(str(hs))
            else:
                expression.append(tag)
            e += 1

        # Now run through the expression converting it to RP
        rpExpression = []
        stack = []
        e = 0
        while e < len(expression):
            tag = expression[e]
            if tag[0].isdigit():
                rpExpression.append(tag)
            else:
                if len(stack) == 0 or stack[-1] == '(':
                    stack.append(tag)
                elif tag == '(':
                    stack.append(tag)
                elif tag == ')':
                    s = stack.pop()
                    while s != '(':
                        rpExpression.append(s)
                        s = stack.pop()
                else:
                    stack.append(tag)
            e += 1
        while len(stack) > 0:
            rpExpression.append(stack.pop())
        self.rpExpression = rpExpression
        return rpExpression



    # Parse the reverse polish set-theoretic expression creating its OpenSCAD equivalent.
    def ToOpenSCAD(self):
        self.Simplify()
        self.ToRP()
        stack=[]
        for s in self.rpExpression:
            if s[0].isdigit():
                stack.append(s)
            else:
                b = stack.pop()
                a = stack.pop()
                stack.append(OpenSCADSetExpression(a, b, s))
        return stack[0]
#********************************************************************************************************

# DORA style set-theoretic geometric model
# See http://adrianbowyer.com/Publications/Better-faster-pictures-Bowyer-Woodwark.pdf

class Box:
    def __init__(self, corner0, corner1):
        self.xs = [corner0[0], corner1[0]]
        self.ys = [corner0[1], corner1[1]]
        self.zs = [corner0[2], corner1[2]]

    def Corner(self, index):
        x = index & 1
        y = (index >> 1) & 1
        z = (index >> 2) & 1
        return [self.xs[x], self.ys[y], self.zs[z]]

    def XLength(self):
        return self.xs[1] - self.xs[0]
    def YLength(self):
        return self.ys[1] - self.ys[0]
    def ZLength(self):
        return self.zs[1] - self.zs[0]

    def Diagonal2(self):
        return self.XLength()**2 + self.YLength()**2 + self.ZLength()**2

    def Divide(self):
        a0 = [self.xs[0], self.ys[0], self.zs[0]]
        b1 = [self.xs[1], self.ys[1], self.zs[1]]
        if self.XLength() > self.YLength():
            if self.XLength() > self.ZLength():
                xHalf = self.xs[0] + self.XLength() / 2.0
                a1 = [xHalf, self.ys[1], self.zs[1]]
                b0 = [xHalf, self.ys[0], self.zs[0]]
            else:
                zHalf = self.zs[0] + self.ZLength() / 2.0
                a1 = [self.xs[1], self.ys[1], zHalf]
                b0 = [self.xs[0], self.ys[0], zHalf]
        else:
            if self.YLength() > self.ZLength():
                yHalf = self.ys[0] + self.YLength() / 2.0
                a1 = [self.xs[1], yHalf, self.zs[1]]
                b0 = [self.xs[0], yHalf, self.zs[0]]
            else:
                zHalf = self.zs[0] + self.ZLength() / 2.0
                a1 = [self.xs[1], self.ys[1], zHalf]
                b0 = [self.xs[0], self.ys[0], zHalf]
        return (Box(a0, a1), Box(b0, b1))

    def IsSolid(self, hsIndex):
        halfSpace = halfSpaces.Get(hsIndex)
        v0 = halfSpace.Value(self.Corner(7))
        if abs(v0) < small:
            return 0
        for c in range(7):
            v = halfSpace.Value(self.Corner(c))
            if abs(v) < small:
                return 0
            if v*v0 < 0:
                return 0
        if v0 < 0:
            return -1
        return 1

smallestDiagonal2Ratio = 0.001
simplestSet = 3

class Model:
    def __init__(self, set, box, parent = None):
        self.set = set
        self.box = box
        self.child0 = None
        self.child1 = None
        self.parent = parent

    def SimplifySet(self):
        self.set.Simplify()
        rpExpression = self.set.ToRP()
        stack = []
        for r in rpExpression:
            if r[0].isdigit():
                ri = int(r)
                s = self.box.IsSolid(ri)
                if s == 0:
                    stack.append(Set(ri))
                elif s < 0:
                    stack.append(Set(-1))
                else:
                    stack.append(Set())
            else:
                b = stack.pop()
                a = stack.pop()
                if r == '|':
                    stack.append(a.Unite(b))
                elif r == '&':
                    stack.append(a.Intersect(b))
                else:
                    print("Model.SimplifySet(): illegal RP operator: " + r)
        self.set = stack[0]
        self.set.Simplify()

    def Root(self):
        if self.parent is None:
            return self
        return self.parent.Root()

    def Divide(self):
        if self.child0 is not None or self.child1 is not None:
            print("Model.Divide(): dividing an already divided model.")
        self.SimplifySet()
        rootDiagonal2 = self.Root().box.Diagonal2()
        if self.box.Diagonal2()/rootDiagonal2 < smallestDiagonal2Ratio:
            return
        if self.set.complexity <= simplestSet:
            return
        boxes = self.box.Divide()
        self.child0 = Model(self.set, boxes[0], parent = self)
        self.child1 = Model(self.set, boxes[1], parent = self)
        self.child0.Divide()
        self.child1.Divide()



# *******************************************************************************************************

# Graphics
# ----------

# A scene consists of a list of triangles for display. The pattern can be rotated on the screen
# and zoomed using the mouse so it can be easily seen.

class Scene:
    def __init__(self, triangles):
        self.triangles = triangles
        self.xangle = 0
        self.yangle = 0
        diag = np.subtract(positiveCorner[0], negativeCorner[0])
        self.distance = -3.0 * (maths.sqrt(np.dot(diag, diag)))

    def MouseDrag(self, dx, dy):
        self.xangle += dx
        self.xangle %= 360
        self.yangle += dy
        self.yangle %= 360

    def MouseScroll(self, dy):
        self.distance += dy

    def Draw(self):
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        gluPerspective(60, 1, 0.1, 1000)
        glMatrixMode(GL_MODELVIEW)

        glEnable(GL_DEPTH_TEST)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        glEnable(GL_LIGHTING)
        glEnable(GL_LIGHT0)
        glEnable(GL_COLOR_MATERIAL)
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE)

        glLoadIdentity()
        glTranslatef(0, 0, self.distance)
        glRotatef(self.xangle, 0, 1, 0)
        glRotatef(self.yangle, -1, 0, 0)

        for triangle in self.triangles:
            glBegin(GL_POLYGON)
            for point in triangle.Points():
                glColor3f(triangle.colour[0], triangle.colour[1], triangle.colour[2])
                glNormal3f(triangle.normal[0], triangle.normal[1], triangle.normal[2])
                glVertex3f(point[0], point[1], point[2])
            glEnd()

        glFlush()


# The World is the container for models (see above) for display.

class World:
    def __init__(self):
        self.element = []

    def AddModel(self, model):
        self.element.append(model)

    def Draw(self):
        for obj in self.element:
            obj.Draw()

    def Setup(self):
        glEnable(GL_DEPTH_TEST)

    def MouseDrag(self, dx, dy):
        for obj in self.element:
            obj.MouseDrag(dx, dy)

    def MouseScroll(self, scroll_y):
        for obj in self.element:
            obj.MouseScroll(scroll_y)


def PutWorldInWindow(world, title):
    win = window.Window(fullscreen=False, vsync=True, resizable=False, height=600, width=600, caption=title)
    range = np.multiply(np.subtract(positiveCorner, negativeCorner), [2, 2, 5])
    negP = np.subtract(middle, range)
    posP = np.add(middle, range)

    @win.event
    def on_resize(width, height):
        glViewport(0, 0, width, height)
        glMatrixMode(GL_PROJECTION)
        glLoadIdentity()
        glOrtho(negP[0], posP[0], negP[1], posP[1], negP[2], posP[2])
        glMatrixMode(GL_MODELVIEW)
        return pyglet.event.EVENT_HANDLED

    @win.event
    def on_mouse_drag(x, y, dx, dy, buttons, modifiers):
        world.MouseDrag(dx, dy)

    @win.event
    def on_mouse_scroll(x, y, scroll_x, scroll_y):
        world.MouseScroll(scroll_y)

    @win.event
    def on_draw():
        glClearColor(0.2, 0.2, 0.2, 0.8)
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT)
        world.Draw()

    world.Setup()


# Set up a window for a list of triangles
def MakeTriangleWindow(triangles, title):
    if graphics == 0:
        return;
    # They may be given a different colour after this function is called.
    triangles = copy.deepcopy(triangles)
    world = World()
    m = Scene(triangles)
    world.AddModel(m)
    PutWorldInWindow(world, title)


# *************************************************************************************************

# The array triangles is assumed to represent the complete surface of a polyhedron.
# This uses the Jordan Curve Theorem to count intersections between the polyhedron
# and a straight line to the centroid of each triangle to decide if the triangle's
# normal is pointing from inside to outside. If it isn't the triangle is inverted.
# The classification is done by counting the number of times the line crosses the polyhedron's
# surface (not including each triangle of interest). If that's odd the line points
# roughly from inside to outside at the triangle. If then the inner product of the line's
# direction and the triangle's normal is negative, the triangle needs to be inverted.
# Similarly, if the count is even and the inner product is positive, the triangle also needs
# to be inverted. Otherwise it is left alone.

def AdjustAttitudes(triangles, hspaces):
    farPoint = np.multiply(positiveCorner, [random.uniform(20, 30), random.uniform(20, 30), random.uniform(20, 30)])
    for triangle in triangles:
        gradient = np.subtract(triangle.centroid, farPoint)
        crossCount = 0
        for t in triangles:
            if t is not triangle:
                hs = hspaces.Get(t.halfSpace)
                parameter = -(np.dot(farPoint, hs.normal) + hs.d) / np.dot(gradient, hs.normal)
                if parameter > 0 and parameter < 1:
                    point = np.add(farPoint, np.multiply(gradient, parameter))
                    if t.ContainsPoint(point):
                        crossCount += 1
        if crossCount % 2 == 1:
            if np.dot(triangle.normal, gradient) < 0:
                triangle.Invert(hspaces)
        else:
            if np.dot(triangle.normal, gradient) > 0:
                triangle.Invert(hspaces)


# *************************************************************************************************

# Walk the graph of triangles making a list of contiguous ones

def GraphWalk(triangle, triangleList):
    if not triangle.flag:
        return
    triangleList.append(triangle)
    triangle.flag = False
    for neighbour in triangle.neighbours:
        GraphWalk(neighbour, triangleList)


# *************************************************************************************************

# Recursive procedure for Woo's alternating sum of volumes algorithm.

def WooStep(trianglesAndPly):
    global finalSet, halfSpaces

    # The triangles for which we want the next convex hull
    # Anything to do?
    triangles = trianglesAndPly[0]
    if len(triangles) <= 0:
        return

    # The original data with coordinates in etc.
    triangleFileData = trianglesAndPly[1]

    # The level of recursion
    level = trianglesAndPly[2]
    if level > maximumLevel:
        exit("WooStep(): maximum recursion exceeded.")

    title = "Level: " + str(level)
    if graphics > 0:
        MakeTriangleWindow(triangles, title + " original triangles")

    # The vertices of those triangles with no duplicates. These are indices of the
    # points in the original ply file.
    pointIndices = []

    # Add a halfspace for each triangle as long as it is not already in the halfspace list
    # Also gather all the vertices
    for triangle in triangles:
        halfSpaces.AddTriangle(triangle)
        for v in triangle.Vertices():
            pointIndices.append(v)

    # Remove duplicate point indices
    pointIndices = list(dict.fromkeys(pointIndices))

    # The actual (x, y, z) coordinates of the points of which to compute the convex hull
    points = []
    for p in pointIndices:
        points.append(triangleFileData.Vertex(p))

    # Find the convex hull
    hull = ConvexHull(points)

    # Make a list of the hull triangles
    hullTriangles = []
    for face in hull.simplices:
        vertices = []
        for f in face:
            # We need to record using the global indexing system
            vertices.append(pointIndices[f])
        triangle = Triangle(vertices, [0.5, 0.5, 1.0], triangleFileData)
        hullTriangles.append(triangle)
        halfSpaces.AddTriangle(triangle)

    # hullTriangles form a complete polyhedron.
    AdjustAttitudes(hullTriangles, halfSpaces)

    # Now make a list of the corresponding halfspaces
    hullHalfSpaces = []
    for triangle in hullTriangles:
        hullHalfSpaces.append(triangle.halfSpace)
    # Remove duplicates
    hullHalfSpaces = list(dict.fromkeys(hullHalfSpaces))


    if graphics > 1:
        MakeTriangleWindow(hullTriangles, title + " CH")

    hullSet = Set(-1)
    for h in hullHalfSpaces:
        hullSet = hullSet.Intersect(Set(h))

    if level == 0:
        finalSet = hullSet
    elif level % 2 == 1:
        finalSet = finalSet.Subtract(hullSet)
    else:
        finalSet = finalSet.Unite(hullSet)

    # Simplfying as we go along increases execution time. But it also allows the simplification to
    # achieve greater reduction.
    finalSet.Simplify()

    # Put all the complements in hullHalfspaces as for the next bit
    # orientation doesn't matter.
    newHullHalfSpaces = []
    for h in hullHalfSpaces:
        newHullHalfSpaces.append(h)
        newHullHalfSpaces.append(halfSpaces.Get(h).complement)

    # I don't THINK this is needed...
    hullHalfSpaces = list(dict.fromkeys(newHullHalfSpaces))

    # Classify each triangle we came in with as being on
    # the hull or not. On means it coincides with a hull halfspace; not, not.
    newTriangles = []
    for triangle in triangles:
        if triangle.halfSpace not in hullHalfSpaces:
            triangle.SetColour([1, 0.5, 0.5])
            newTriangles.append(triangle)

    if graphics > 2:
        MakeTriangleWindow(newTriangles, title + " Inside triangles")

    if len(newTriangles) > 0:
        for triangle in newTriangles:
            triangle.flag = True

        contiguousTriangleLists = []
        for triangle in newTriangles:
            if triangle.flag:
                contiguousTriangles = []
                GraphWalk(triangle, contiguousTriangles)
                contiguousTriangleLists.append(contiguousTriangles)

        # Carry on down...
        for contiguousTriangles in contiguousTriangleLists:
            trianglesAndPly = (contiguousTriangles, triangleFileData, level + 1)
            WooStep(trianglesAndPly)


# **************************************************************************************************

# Create the .set file

def ToSetFile(fileName, halfSpaces, set):
    file = open(fileName, "w")

    # The first things in the file are the bounding box
    file.write(str(negativeCorner) + "\n")
    file.write(str(positiveCorner) + "\n")

    # Get the RP expression to write out at the end
    set.Simplify()
    rp = set.ToRP()

    # Count the number of unique halfSpaces in the RP expression and set their flags
    count = 0
    for item in rp:
        if item[0].isdigit():
            hs = int(item)
            halfSpace = halfSpaces.Get(hs)
            if not halfSpace.flag:
                halfSpace.SetFlag(True)
                count += 1

    # Write out the count of the halfSpaces followed by the coefficients of each one
    file.write(str(count) + "\n")
    for c in range(len(halfSpaces)):
        hs = halfSpaces.Get(c)
        if hs.flag:
            file.write(str(c) + " " + str(hs.Coeficients()) + "\n")

    # Write out the reverse polish expression
    for item in rp:
        file.write(item + " ")
    file.write("\n")

    # Tidy up
    file.close()
    halfSpaces.ResetFlags()


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

# Read the set file
def ReadSetFile(fileName):
    global negativeCorner, positiveCorner, diagonal, middle, halfSpaces, finalSet
    file = open(fileName, "r")
    negativeCorner = GetListFromLine(file)
    positiveCorner = GetListFromLine(file)
    diagonal = np.subtract(positiveCorner, negativeCorner)
    middle = np.add(negativeCorner, np.multiply(diagonal, 0.5))
    diagonal = maths.sqrt(np.dot(diagonal, diagonal))
    halfSpaceCount = round(GetListFromLine(file)[0])
    hsIndices = []
    hsList = []
    halfSpaces = HalfSpaceList()
    for hsIndex in range(halfSpaceCount):
        hs = GetListFromLine(file)
        hsIndices.append(hs[0])
        del hs[0]
        hs = MakeHalfSpaceFromCoefficients(hs)
        hsList.append(hs)
        halfSpaces.AddSpace(hs)
    set = file.readline()
    set = set.replace("\n", "")
    set = set.split(" ")
    stack = []
    for s in set:
        if not s == '':
            if s[0].isdigit():
                hs = hsList[hsIndices.index(int(s))]
                stack.append(Set( halfSpaces.LookUp(hs) ))
            else:
                b = stack.pop()
                a = stack.pop()
                if s == '&':
                    stack.append(a.Intersect(b))
                elif s == '|':
                    stack.append(a.Unite(b))
                else:
                    print("ReadSetFile(): illegal set operator: " + s)
    finalSet = stack[0]


# Create the .set file

def ToOpenSCAD(fileName, halfSpaces, set):
    file = open(fileName, "w")

    # Get the RP expression to write out at the end
    set.Simplify()
    rp = set.ToRP()

    # Count the number of unique halfSpaces in the RP expression and set their flags
    count = 0
    for item in rp:
        if item[0].isdigit():
            hs = int(item)
            halfSpace = halfSpaces.Get(hs)
            if not halfSpace.flag:
                halfSpace.SetFlag(True)
                count += 1

    for c in range(len(halfSpaces)):
        hs = halfSpaces.Get(c)
        if hs.flag:
            file.write(hs.ToOpenSCAD(c))


    file.write(set.ToOpenSCAD())
    file.write("\n")

    # Tidy up
    file.close()
    halfSpaces.ResetFlags()


# ***********************************************************************************************************

# Run the conversion
# -------------------

#ReadSetFile("STL2CSG-test-objects-woo-2.set")
#exit(0)

halfSpaces = HalfSpaceList()

name = "STL2CSG-test-objects-woo-2"
#name = "projections-bottom-end"
place = "../../"
# fileName = '../../cube.ply'
# fileName = '../../two-disjoint-cubes.ply'
# fileName = '../../two-overlapping-cubes.ply'
# fileName = '../../hole-enclosed-in-cylinder.ply'
# fileName = '../../two-nonmanifold-cubes.ply'
# fileName = '../../two-nasty-nonmanifold-cubes.ply'
# fileName = '../../554.2-extruder-drive-pneumatic.ply'
# fileName = '../../cube-1-cylinder-1.ply'
# fileName = '../../STL2CSG-test-objects-woo-1.ply'
# fileName = '../../STL2CSG-test-objects-woo-2.ply'
# fileName = '../../STL2CSG-test-objects-cube-cylinder.ply'
# fileName = '../../STL2CSG-test-objects-cubePlusCylinder.ply'

# Load up the .ply file of triangles
triangleFileData = TriangleFileData(place + name + ".ply")

# Find the bounding box
for v in range(triangleFileData.VertexCount()):
    vertex = triangleFileData.Vertex(v)
    middle = np.add(middle, vertex)
    positiveCorner = np.maximum(positiveCorner, vertex)
    negativeCorner = np.minimum(negativeCorner, vertex)
middle = np.multiply(middle, 1.0 / triangleFileData.VertexCount())
diagonal = np.subtract(positiveCorner, negativeCorner)
diagonal = maths.sqrt(np.dot(diagonal, diagonal))

# Make a list of all the original triangles
originalTriangles = []
for t in range(triangleFileData.TriangleCount()):
    triangle = Triangle(triangleFileData.Triangle(t), [0.5, 1.0, 0.5], triangleFileData)
    originalTriangles.append(triangle)

ConstructNeighbourLists(originalTriangles)

# The recursive alternating sum of volumes function needs
# the list of triangles to work with, the .ply file data
# to get vertex coordinates from, and the current level
# of recursion (here at the start, 0).
trianglesAndPly = (originalTriangles, triangleFileData, 0)

# Run the algorithm
WooStep(trianglesAndPly)

# Save and maybe plot the results
finalSet.Simplify()
'''
points = [
    [-1, -1, -1],
    [5, 5, 5],
    [5, 5, 11],
    [5, 5, 16]
]
for point in points:
    print(str(point) + " has value " + str(finalSet.Value(point)))
'''

model = Model(finalSet, Box(negativeCorner, positiveCorner))
model.Divide()

ToOpenSCAD(name + ".scad", halfSpaces, finalSet)
if graphics > 0:
    pyglet.app.run()
