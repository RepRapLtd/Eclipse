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
	
	private CornerList cornerList = new CornerList();
	private double smallD2 = 0.001;
	private ArrayList<Triangle> shells = new ArrayList<Triangle>();
	private int triangleCount = 0;
	private int cornerCount = 0;
	private boolean debug = true;
	private BoundingBox boundingBox = new BoundingBox();
	private Point3d cloudCentroid = new Point3d(0, 0, 0);
	private Vector3d principalComponents = new Vector3d(0, 0, 0);
	private double shortestEdge2 = Double.MAX_VALUE;
	private double longestEdge2 = 0;
	
	private void Error(String s)
	{
		System.err.println(s);
	}
	
	
	public class CornerList
	{
		/**
		 * This class stores the list of vertices/corners of all the triangles as points in space.  Each
		 * point will have a ring of triangles around it.  The points all have one arbitrary triangle
		 * stored alongside them (from which the ring can be derived by a local search).
		 */
		
		private ArrayList<Point3d> corners = new ArrayList<Point3d>();
		private ArrayList<Triangle> triangles = new ArrayList<Triangle>();
		
		/**
		 * Add a triangle corner to the list, together with any
		 * triangle that has that corner.  It returns the index of
		 * the corner.
		 * 
		 * @param corner
		 * @param triangle
		 * @return
		 */
		public int AddCorner(Point3d corner, Triangle triangle)
		{
			corners.add(corner);
			triangles.add(triangle);
			return corners.size() - 1;
		}
		
		public Point3d GetCorner(int i)
		{
			return corners.get(i);
		}
		
		public Triangle GetTriangle(int i)
		{
			return triangles.get(i);
		}
		
		/**
		 * Sanity check.
		 * @return
		 */
		public boolean Sane()
		{
			boolean result = (corners.size() == triangles.size());
			if(!result)
			{
				Error("Sane(): These should be equal - corners: " + String.valueOf(corners.size()) +
						", triangles: " + String.valueOf(triangles.size())	);
			}
			return result;
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
		public int Find(Point3d target)
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
		
		public boolean StitchUpThisTriangle(Triangle triangle, ArrayList<Triangle> extraTriangles)
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
					Error("StitchUpThisTriangle() - triangle with missing neighbour found.");
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
		public boolean StitchUp(ArrayList<Triangle> extraTriangles)
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
	}
	
	public class Triangle 
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
		 * visited is a flag to facilitate recursively walking over the mesh
		 */
		
		private boolean visited;
		
		private Triangle()
		{
			visited = false;
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
		 */
		
		public Triangle(Point3d aC, Point3d bC, Point3d cC) 
		{
			this();
			
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
		
		public void ResetVisited()
		{
			visited = false;
		}
		
		public void MarkVisited()
		{
			visited = true;
		}
		
		public boolean Visited()
		{
			return visited;
		}
		
		public int GetCorner(int i)
		{
			return corners[i];
		}
		
		public Triangle GetNeighbour(int i)
		{
			return neighbours[i];
		}
		
		
		/**
		 * Reset all the triangles in one shell to unvisited, starting with this one.
		 * Note this works with partial triangulations with null edges as long as they
		 * are not disjoint.
		 */
		
		public void Reset()
		{
			ResetVisited();
			for(int i = 0; i < 3; i++)
			{
				if(neighbours[i] != null)
					if(neighbours[i].Visited())
						neighbours[i].Reset();
			}
		}
		
		/**
		 * The vector product of corners ([1] - [0])x([2] - [0]) faces
		 * outward from solid to air.  This inverts a triangle that has 
		 * the wrong sense. 
		 */
		
		public void Invert()
		{
			int cTemp = corners[2];
			corners[2] = corners[1];
			corners[1] = cTemp;
			Triangle nTemp = neighbours[2];
			neighbours[2] = neighbours[1];
			neighbours[1] = nTemp;
		}
		
		public void SetEdge(Triangle neighbour, int i)
		{
			neighbours[i] = neighbour;
		}
		
		public boolean HasCorner(int corner)
		{
			for(int i = 0; i < 3; i++)
			{
				if(corners[i] == corner)
					return true;
			}		
			return false;
		}		
		
		public boolean HasCorners(int corner1, int corner2)
		{
			return HasCorner(corner1) && HasCorner(corner2);
		}
	
	}


	/**
	 * 
	 */
	public Triangulator(String location) 
	{
		ArrayList<Triangle> extraTriangles = new ArrayList<Triangle>();
		
		STLLoader loader = new STLLoader();
		Scene scene = null;
		
		try 
        {
            scene = loader.load(location);
        } catch (Exception e) 
	    {
	        Error("Triangulator(): Exception loading STL file from: " + location);
	        e.printStackTrace();
	    }
		
		if(scene == null)
		{
			Error("Triangulator(): no triangulation to process");
			System.exit(1);
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
        
        // After the stats are gathered the centroid is the sum of the points
        
        cloudCentroid.scale(1.0/(double)cornerCount);
        
        /*
         * Now go through all the triangles again adding all the triangles to the structure.
         */
        
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
        	Error("Triangulator(): Sane() false after adding all the triangles.");
        	System.exit(2);
        }
        
        /*
         * Finally stitch the triangulation together.
         */
        
        cornerList.StitchUp(extraTriangles);
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
			
			triangleCount++;
			cornerCount += 3;
		}
	}
	
	private void AddATriangle(Point3d corner0, Point3d corner1, Point3d corner2, ArrayList<Triangle> extraTriangles)
	{
		Triangle triangle = new Triangle(corner0, corner1, corner2);
		if(triangle.Visited())
			triangle.ResetVisited();
		else
			extraTriangles.add(triangle);
	}
	
	private void AddSomeTriangles(GeometryArray g, ArrayList<Triangle> extraTriangles)
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
			AddATriangle(corner0, corner1, corner2, extraTriangles);
		}
	}
	
	public void PrintStatistics()
	{
		System.out.println("Bounding box: " + boundingBox.toString());
		System.out.println("Cloud Centroid: " + cloudCentroid.toString());
		System.out.println("Shortest edge: " + String.valueOf(Math.sqrt(shortestEdge2)) );
		System.out.println("Longest edge: " + String.valueOf(Math.sqrt(longestEdge2)) );
		System.out.println("Points coincide if closer than: " + String.valueOf(Math.sqrt(smallD2)) );
		System.out.println("Number of triangles: " + String.valueOf(triangleCount) );
		System.out.println("Number of corners: " + String.valueOf(cornerCount) );
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
		t.PrintStatistics();
	}

}
