/**
 * This program reads a list of triangles from an STL file; STLs have no imposed structure (Boo!).
 * 
 * It attempts to stitch them together to form shells of triangles where each triangle has
 * three corners in space and opposite them three sides.  Each side connects to another triangle 
 * in the shell.
 * 
 * A certain amount of checking for integrity and solidity is done.
 * 
 * There are a number of TODOs.  These mainly relate to imposing a logarithmic spatial structure
 * such as a BSP tree on the whole problem, which would obviously make spatial searching much more efficient.
 * 
 */

/**
 * @author Adrian Bowyer
 * 
 * RepRap Ltd
 * https://reprapltd.com
 * 
 * Date: about 24 October 2018
 * 
 * Licence: GPL
 *
 */

package com.reprapltd.triangulator;

import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import org.j3d.renderer.java3d.loaders.STLLoader;
import com.sun.j3d.loaders.Scene;
import javax.vecmath.Point3d;
import java.util.Hashtable;
import javax.media.j3d.BoundingBox;
import javax.media.j3d.GeometryArray;
import javax.media.j3d.Shape3D;
import javax.vecmath.Matrix3d;
import javax.vecmath.Vector3d;


public class Triangulator 
{
	/**
	 * This is the enclosing class that constructs and maintains the entire triangulation.
	 * 
	 * The triangulation consists of one or more shells.  These are accessed by a single triangle
	 * that is a part of each shell stored in the list meshes.  Given one triangle, every triangle
	 * in a shell can be visited by a recursive neighbour graph walk.
	 */
	
	private boolean debugging = true;
	
	private CornerList cornerList = null;
	private double smallD2;
	private ArrayList<Triangle> shells = null;
	private int cornerCount;  // Number of points read in >> cornerList.corners.size()
	private int triangleCount;
	private BoundingBox boundingBox = null;
	private Point3d cloudCentroid = null;
	//private Vector3d principalComponents = new Vector3d(0, 0, 0);
	private double shortestEdge2;
	private double longestEdge2;
	private String STLLocation = null;
	
	private void Error(String s, boolean serious)
	{
		if(!serious && !debugging)
			return;
		System.err.println(s);
	}
	
	
	private class CornerList
	{
		/**
		 * This class stores the list of vertices/corners of all the triangles as points in space.  Each
		 * point will have a ring of triangles around it.  The points all have one arbitrary triangle
		 * stored alongside them (from which the ring can be derived by a local search).
		 */
		
		private ArrayList<Point3d> corners = null;
		private ArrayList<Triangle> triangles = null;
		
		private CornerList()
		{
			corners = new ArrayList<Point3d>();
			triangles = new ArrayList<Triangle>();
		}
		
		/**
		 * Add a triangle corner to the list, together with any
		 * triangle that has that corner.  It returns the index of
		 * the corner.
		 * 
		 * @param corner
		 * @param triangle
		 * @return
		 */
		private int AddCorner(Point3d corner, Triangle triangle)
		{
			corners.add(corner);
			triangles.add(triangle);
			return corners.size() - 1;
		}
		
		private Point3d GetCorner(int i)
		{
			return corners.get(i);
		}
		
		private Triangle GetTriangle(int i)
		{
			return triangles.get(i);
		}
		
		private int CornerCount()
		{
			return corners.size();
		}
		
		/**
		 * Sanity check.
		 * @return
		 */
		private boolean Sane()
		{
			boolean result = (corners.size() == triangles.size());
			if(!result)
			{
				Error("Sane(): These should be equal - corners: " + String.valueOf(corners.size()) +
						", triangles: " + String.valueOf(triangles.size()) , false);
			}
			return result;
		}
		
		private void SavePoints(PrintStream out)
		{
			out.println(String.valueOf(corners.size()));
			for(int i = 0; i < corners.size(); i++)
			{
				Point3d p = corners.get(i);
				out.print(String.valueOf(p.x) + " ");
				out.print(String.valueOf(p.y) + " ");
				out.print(String.valueOf(p.z) + " ");
				out.println(String.valueOf(triangles.get(i).Index()));
			}
		}
		
		/**
		 * This finds the point in the list corresponding to target and
		 * returns its index.  Any point with squared distance to target
		 * less than or equal to smallD2 is considered a match.  If no such 
		 * point is found -1 is returned.
		 * 
		 * TODO: Obviously this should not do an exhaustive search; it
		 * should use a spatial structure (BSP tree?).
		 * 
		 * @param target
		 * @return
		 */
		private int Find(Point3d target)
		{
			for(int i = 0; i < corners.size(); i++)
			{
				if(corners.get(i).distanceSquared(target) <= smallD2)
					return i;
			}
			return -1;
		}
		
		/**
		 * This finds the triangle in the lists that shares the edge corner1-corner2 with notThisOne.
		 * 
		 * TODO: Again obviously this should not do an exhaustive search, but should use spatial structure.
		 * 
		 * @param corner0
		 * @param corner1
		 * @param extraTriangles
		 * @param notThisOne
		 * @return
		 */
		private Triangle FindTriangleWithEdge(int corner0, int corner1, ArrayList<Triangle> extraTriangles, Triangle notThisOne)
		{	
			for(int i = 0; i < triangles.size(); i++)
			{
				Triangle triangle = triangles.get(i);
				if(triangle != notThisOne)
					if(triangle.HasCorners(corner0, corner1))
						return triangle;
			}
			
			for(int i = 0; i < extraTriangles.size(); i++)
			{
				Triangle triangle = extraTriangles.get(i);
				if(triangle != notThisOne)
					if(triangle.HasCorners(corner0, corner1))
						return triangle;
			}
			
			return null;
		}
		
		/**
		 * Run round the three edges of a triangle finding its neighbours that share that edge
		 * and record the neighbour relationship.
		 * 
		 * TODO: AGAIN obviously this should set the neighbour relationship for both triangles at once.
		 * 
		 * @param triangle
		 * @param extraTriangles
		 * @return
		 */
		
		private boolean StitchUpThisTriangle(Triangle triangle, ArrayList<Triangle> extraTriangles)
		{
			boolean result = true;
			
			for(int j = 0; j < 3; j++)
			{
				int j1 = (j+1)%3;
				int corner0 = triangle.GetCorner(j);
				int corner1 = triangle.GetCorner(j1);
				Triangle neighbour = FindTriangleWithEdge(corner0, corner1, extraTriangles, triangle);
				if(neighbour == null)
				{
					Error("StitchUpThisTriangle() - triangle with missing neighbour found.", false);
					result = false;
				}
				triangle.SetEdge(neighbour, 3 - (j + j1));
			}
			
			return result;
		}
		
		/**
		 * This sets all the neighbours of all the triangles in all the shells.
		 * 
		 * @param extraTriangles
		 * @return
		 */
		private boolean StitchUp(ArrayList<Triangle> extraTriangles)
		{
			boolean result = true;
			
			for(int i = 0; i < triangles.size(); i++)
			{
				Triangle triangle = triangles.get(i);
				if(!StitchUpThisTriangle(triangle, extraTriangles))
					result = false;	
			}
			
			for(int i = 0; i < extraTriangles.size(); i++)
			{
				Triangle triangle = extraTriangles.get(i);
				if(!StitchUpThisTriangle(triangle, extraTriangles))
					result = false;	
			}
			
			return result;
		}
		
		/**
		 * Search the corners' triangles to find an unvisited triangle.
		 * If none is found, return null.
		 * 
		 * This cannot avoid a global search.  But at least it should be quick.
		 * 
		 * @return
		 */
		private Triangle FindAnUnvisitedTriangle()
		{
			for(int i = 0; i < triangles.size(); i++)
				if(!triangles.get(i).Visited())
					return triangles.get(i);
			return null;
		}
	}
	
	private class Triangle 
	{

		/**
		 * This holds the representation of a single triangle.  A triangle
		 * has corners and, opposite them, neighbouring triangles.  Each neighbour is the adjacent triangle
		 * in the mesh.  For a complete mesh no neighbours should be null.  The corners
		 * are stored as indices into cornerList, where the actual point coordinates
		 * in space are kept.  
		 */
		
		private int corners[];
		private Triangle neighbours[];
		
		/**
		 * The index is just used for input and output; every index
		 * should be unique.
		 */
		private int index;
		
		/**
		 * visited is a flag to facilitate recursively walking over the mesh
		 */
		private boolean visited;
		
		private Triangle()
		{
			visited = false;
			index = -1;
			corners = new int[3];
			neighbours = new Triangle[3];
			for(int i = 0; i < 3; i++)
			{
				corners[i] = -1;
				neighbours[i] = null;
			}
		}
		
		/**
		 * Constructor to make a triangle from three points in space.
		 * 
		 * If visited is true on return, the new triangle is also in
		 * the corner list and visited should be reset to false.  If it is false
		 * it is not in the corner list (because all three of its corners
		 * were put in the list by previous triangles) and it may need to
		 * be recorded elsewhere.
		 * 
		 * @param aC
		 * @param bC
		 * @param cC
		 * @param ix
		 */
		
		private Triangle(Point3d aC, Point3d bC, Point3d cC, int ix) 
		{
			this();
			
			index = ix;
			
			int c = cornerList.Find(aC);
			if(c < 0)
			{
				c = cornerList.AddCorner(aC,  this);
				visited = true;
			}
			corners[0] = c;
			
			c = cornerList.Find(bC);
			if(c < 0)
			{
				c = cornerList.AddCorner(bC,  this);
				visited = true;
			}
			corners[1] = c;
			
			c = cornerList.Find(cC);
			if(c < 0)
			{
				c = cornerList.AddCorner(cC,  this);
				visited = true;
			}
			corners[2] = c;
		}
		
		private void ResetVisited()
		{
			visited = false;
		}
		
		private void SetVisited()
		{
			visited = true;
		}
		
		private boolean Visited()
		{
			return visited;
		}
		
		private int GetCorner(int i)
		{
			return corners[i];
		}
		
		private Triangle GetNeighbour(int i)
		{
			return neighbours[i];
		}
		
		private int Index()
		{
			return index;
		}		
	
		/**
		 * Set all the triangles in one shell to visited, starting with this one.
		 * Note this works with partial triangulations with null edges as long as they
		 * are not disjoint, though this will produce a debug warning.
		 * It counts the triangles as a side-effect and returns the total.
		 */
		private int Set()
		{
			int count = 1;
			SetVisited();
			for(int i = 0; i < 3; i++)
			{
				if(neighbours[i] != null)
				{
					if(!neighbours[i].Visited())
						count += neighbours[i].Set();
				} else
					Error("Set(): triangle with null neighbour found.", false);
			}
			return count;
		}
		
		/**
		 * Reset all the triangles in one shell to unvisited, starting with this one.
		 * Note this works with partial triangulations with null edges as long as they
		 * are not disjoint.
		 * It counts the triangles as a side-effect and returns the total.
		 */
		private int Reset()
		{
			int count = 1;
			ResetVisited();
			for(int i = 0; i < 3; i++)
			{
				if(neighbours[i] != null)
				{
					if(neighbours[i].Visited())
						count += neighbours[i].Reset();
				} else
					Error("Reset(): triangle with null neighbour found.", false);
			}
			return count;
		}
		
		/**
		 * Do a set and reset to count the triangles (twice).
		 * 
		 * TODO: the reset could record the count in each triangle for
		 * instant return in the future; lazy evaluation.
		 * 
		 * @return
		 */
		private int Count()
		{
			int count0 = Set();
			int count1 = Reset();
			
			if(count1 != count0)
				Error("Count(): set count (" + String.valueOf(count0) + ") not equal to reset count (" + String.valueOf(count1) + ").", true);
			
			return count0;
		}
		
		
		private void SaveTrianglesR(PrintStream out)
		{
			out.print(String.valueOf(index) + " ");
			for(int i = 0; i < 3; i++)
				out.print(String.valueOf(corners[i]) + " ");
			for(int i = 0; i < 2; i++)
				out.print(String.valueOf(neighbours[i].Index()) + " ");
			out.println(String.valueOf(neighbours[2].Index()));
			
			SetVisited();
			for(int i = 0; i < 3; i++)
			{
				if(neighbours[i] != null)
				{
					if(!neighbours[i].Visited())
						neighbours[i].SaveTrianglesR(out);
				} else
					Error("SaveTrianglesR(): triangle with null neighbour found.", true);
			}			
		}
		
		
		private void SaveTriangles(PrintStream out)
		{
			int count = Count();
			out.println(String.valueOf(count) + " ");
			SaveTrianglesR(out);
			Reset();
		}
		
		/**
		 * The vector product of corners ([1] - [0])x([2] - [0]) faces
		 * outward from solid to air.  This inverts a triangle that has 
		 * the wrong sense. 
		 */
		
		private void Invert()
		{
			int cTemp = corners[2];
			corners[2] = corners[1];
			corners[1] = cTemp;
			Triangle nTemp = neighbours[2];
			neighbours[2] = neighbours[1];
			neighbours[1] = nTemp;
		}
		
		private void SetEdge(Triangle neighbour, int i)
		{
			neighbours[i] = neighbour;
		}
		
		private boolean HasCorner(int corner)
		{
			for(int i = 0; i < 3; i++)
			{
				if(corners[i] == corner)
					return true;
			}		
			return false;
		}		
		
		private boolean HasCorners(int corner1, int corner2)
		{
			return HasCorner(corner1) && HasCorner(corner2);
		}
	
	}


	/**
	 * 
	 */
	public Triangulator(String location) 
	{
		cornerList = new CornerList();
		smallD2 = 0.001;
		shells = new ArrayList<Triangle>();
		cornerCount = 0;  // Number of points read in >> cornerList.corners.size()
		triangleCount = 0;
		boundingBox = new BoundingBox();
		cloudCentroid = new Point3d(0, 0, 0);
		//principalComponents = new Vector3d(0, 0, 0);
		shortestEdge2 = Double.MAX_VALUE;
		longestEdge2 = 0;
		STLLocation = "";
		
		ArrayList<Triangle> extraTriangles = new ArrayList<Triangle>();
		
		STLLoader loader = new STLLoader();
		Scene scene = null;
		
		STLLocation = location;
		
		try 
        {
            scene = loader.load(location);
        } catch (Exception e) 
	    {
	        Error("Triangulator(): Exception loading STL file from: " + location, false);
	        e.printStackTrace();
	        System.exit(1);
	    }
		
		if(scene == null)
		{
			Error("Triangulator(): no triangulation to process", false);
			System.exit(2);
		}
        
        /*
         * First go through all the triangles gathering statistics.
         */
        
        Hashtable<?,?> namedObjects = scene.getNamedObjects();
        java.util.Enumeration<?> enumValues = namedObjects.elements();
        
        if(enumValues != null) 
        {
            while(enumValues.hasMoreElements()) 
            {
            	Shape3D stlTriangles = (Shape3D)enumValues.nextElement();
                UpdateStatistics((GeometryArray)stlTriangles.getGeometry());
            }
        }
        
        // It seems reasonable to assume that any points closer together than a quarter of the
        // length of the shortest edge are coincident. 16 is because we use squared distances throughout.
        
        smallD2 = shortestEdge2/16;
        
        // After the stats are gathered the centroid is the sum of all the points
        
        cloudCentroid.scale(1.0/(double)cornerCount);
        
        /*
         * Now go through all the triangles again adding all the triangles to the structure.
         */
        
        triangleCount = 0; // Defensive...
        namedObjects = scene.getNamedObjects();
        enumValues = namedObjects.elements();
        
        if(enumValues != null) 
        {
            while(enumValues.hasMoreElements()) 
            {
            	Shape3D stlTriangles = (Shape3D)enumValues.nextElement();
                AddSomeTriangles((GeometryArray)stlTriangles.getGeometry(), extraTriangles);
            }
        }
        
        if(!cornerList.Sane())
        {
        	Error("Triangulator(): Sane() false after adding all the triangles.", false);
        	System.exit(3);
        }
        
        /*
         * Stitch the triangulation together.
         */
        
        cornerList.StitchUp(extraTriangles);
        
        /*
         * Discover and record the shells. 
         */
        
        Triangle t;
        do
        {
        	t = cornerList.FindAnUnvisitedTriangle();
        	if(t != null)
        	{
                shells.add(t);
                t.Set();
        	}
        } while (t != null);
        
        for(int i = 0; i < shells.size(); i++)
        	shells.get(i).Reset();
	}
	
	private void UpdatePointStatistics(Point3d corner)
	{
		boundingBox.combine(corner);
		cloudCentroid.add(corner);
	}
	
	private void UpdateEdgeStatistics(Point3d corner0, Point3d corner1)
	{
		double d2 = corner0.distanceSquared(corner1);
		if(d2 > longestEdge2)
			longestEdge2 = d2;
		if(d2 < shortestEdge2)
			shortestEdge2 = d2;
	}
	
	private void UpdateStatistics(GeometryArray g)
	{
		if(g == null)
			return;

		Point3d corner0 = new Point3d();
		Point3d corner1 = new Point3d();
		Point3d corner2 = new Point3d();

		for(int i = 0; i < g.getVertexCount(); i += 3) 
		{
			g.getCoordinate(i, corner0);
			g.getCoordinate(i+1, corner1);
			g.getCoordinate(i+2, corner2);
			
			UpdatePointStatistics(corner0);
			UpdatePointStatistics(corner1);
			UpdatePointStatistics(corner2);
			
			UpdateEdgeStatistics(corner0, corner1);
			UpdateEdgeStatistics(corner1, corner2);
			UpdateEdgeStatistics(corner2, corner0);
			
			cornerCount += 3;
		}
	}
	
	private void AddATriangle(Point3d corner0, Point3d corner1, Point3d corner2, ArrayList<Triangle> extraTriangles)
	{
		Triangle triangle = new Triangle(corner0, corner1, corner2, triangleCount);
		if(triangle.Visited())
			triangle.ResetVisited();
		else
			extraTriangles.add(triangle);
		triangleCount++;
	}
	
	private void AddSomeTriangles(GeometryArray g, ArrayList<Triangle> extraTriangles)
	{
		if(g == null)
			return;

		for(int i = 0; i < g.getVertexCount(); i += 3) 
		{
			Point3d corner0 = new Point3d();
			Point3d corner1 = new Point3d();
			Point3d corner2 = new Point3d();
			
			g.getCoordinate(i, corner0);
			g.getCoordinate(i+1, corner1);
			g.getCoordinate(i+2, corner2);
			AddATriangle(corner0, corner1, corner2, extraTriangles);
		}
	}
	
	/**
	 * TODO: Maybe this should use XML (like LandXML, which is overkill).  But there
	 * doesn't seem to be a simple defined XML triangle mesh format, and there's not much point in
	 * making one up as only this would use it.
	 * 
	 * @param location
	 */
	private void Save(String location)
	{
		PrintStream out = null;
		try 
		{
			out = new PrintStream(location);
			cornerList.SavePoints(out);
			out.println(String.valueOf(shells.size()));
			for(int i = 0; i < shells.size(); i++)
				shells.get(i).SaveTriangles(out);
			out.close();
		} catch (IOException e) 
		{
			Error("Save(): Error writing to " + location, true);
			e.printStackTrace();
		}
	}
	
	private void PrintStatistics()
	{
		System.out.println("Statistics for STL file " + STLLocation +":");
		System.out.println(" Bounding box: " + boundingBox.toString());
		System.out.println(" Point cloud Centroid: " + cloudCentroid.toString());
		System.out.println(" Shortest edge: " + String.valueOf(Math.sqrt(shortestEdge2)) );
		System.out.println(" Longest edge: " + String.valueOf(Math.sqrt(longestEdge2)) );
		System.out.println(" Points coincide if closer than: " + String.valueOf(Math.sqrt(smallD2)) );
		System.out.println(" Number of unique corner/vertex points: " + String.valueOf(cornerList.CornerCount()) );
		System.out.println(" Number of corner points in the STL file (including duplicates): " + String.valueOf(cornerCount) );
		System.out.println(" Number of triangles read: " + String.valueOf(triangleCount) );
		System.out.println(" Number of shells: " + String.valueOf(shells.size()) );
		for(int i = 0; i < shells.size(); i++)
			System.out.println("  Shell " + String.valueOf(i) + " has " + String.valueOf(shells.get(i).Count()) + " triangles.");
	}
	
    private double tetVolume(Point3d a, Point3d b, Point3d c, Point3d d)
    {
    	Matrix3d m = new Matrix3d(b.x - a.x, c.x - a.x, d.x - a.x, b.y - a.y, c.y - a.y, d.y - a.y, b.z - a.z, c.z - a.z, d.z - a.z);
    	return m.determinant()/6.0;
    }
	
    private double prismVolume(Point3d a, Point3d b, Point3d c)
    {
    	Point3d d = new Point3d(a.x, a.y, 0); 
    	Point3d e = new Point3d(b.x, b.y, 0);  
    	Point3d f = new Point3d(c.x, c.y, 0);
    	return tetVolume(a, b, c, e) +
    		tetVolume(a, e, c, d) +
    		tetVolume(e, f, c, d);
    }
		
	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		Triangulator t = new Triangulator("file:///home/ensab/Desktop/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/triangulator/test-cube.stl");
		//Triangulator t = new Triangulator("file:///home/ensab/Desktop/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/triangulator/two-disjoint-cubes.stl");
		//Triangulator t = new Triangulator("file:///home/ensab/Desktop/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/triangulator/two-overlapping-cubes.stl");
		//Triangulator t = new Triangulator("file:///home/ensab/Desktop/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/triangulator/hole-enclosed-in-cylinder.stl");
		//Triangulator t = new Triangulator("file:///home/ensab/Desktop/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/triangulator/two-nonmanifold-cubes.stl");
		//Triangulator t = new Triangulator("file:///home/ensab/Desktop/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/triangulator/two-nasty-nonmanifold-cubes.stl");
		t.PrintStatistics();
		
		t.Save("triangulation.tri");
	}

}
