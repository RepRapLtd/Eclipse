/**
 * This program reads a list of triangles from an STL file, which have no imposed structure.
 * 
 * It attempts to stitch them together to form meshes of triangles where each triangle has
 * three corners in space and opposite them three sides.  Each side connects to another triangle 
 * in the mesh.  The result is then written to a new file with all that structure intact.
 * 
 * A certain amount of checking for integrity and solidity is done.
 * 
 */
package com.reprapltd.triangulator;

import org.j3d.renderer.java3d.loaders.STLLoader;
import javax.vecmath.Point3d;

/**
 * @author Adrian Bowyer
 * 
 * RepRap Ltd
 * https://reprapltd.com
 * 
 * Licence: GPL
 *
 */
public class Triangulator 
{
	
	public class CornerList
	{
		/**
		 * This stores the list of vertices/corners of all the triangles as points in space.  Each
		 * point will have a ring of triangles around it.  The points all have one arbitrary one 
		 * of those triangles stored alongside them (from which the ring can be derived by a local search).
		 */
	}
	
	public class Triangle 
	{

		/**
		 * This holds the representation of a single triangle.  A triangle
		 * has corners and, opposite them, edges.  Each edge is the adjacent triangle
		 * in the mesh.  For a complete mesh no edges should be null.  The corners
		 * are stored as indices into the CornerList, where the actual point coordinates
		 * in space are kept.
		 */
		
		private int cornerIndices[];
		private Triangle edges[];
		
		/**
		 * visited is a flag to facilitate recursively walking over the mesh
		 */
		
		private boolean visited;  
		
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
		
		/**
		 * Reset all the triangles in one mesh to unvisited, starting with this one.
		 * Note this works with partial triangulations with null edges as long as they
		 * are not disjoint.
		 */
		
		public void Reset()
		{
			ResetVisited();
			for(int i = 0; i < 3; i++)
			{
				if(edges[i] != null)
					if(edges[i].Visited())
						edges[i].Reset();
			}
		}
		
		/**
		 * The vector product of corners ([1] - [0])x([2] - [0]) faces
		 * outward from solid to air.  This inverts a triangle that has 
		 * the wrong sense. 
		 */
		
		public void Invert()
		{
			int cTemp = cornerIndices[2];
			cornerIndices[2] = cornerIndices[1];
			cornerIndices[1] = cTemp;
			Triangle eTemp = edges[2];
			edges[2] = edges[1];
			edges[1] = eTemp;
		}
		
		public void SetEdge(Triangle e, int i)
		{
			edges[i] = e;
		}
		
		public Triangle()
		{
			cornerIndices = new int[3];
			edges = new Triangle[3];
			for(int i = 0; i < 3; i++)
			{
				cornerIndices[i] = -1;
				edges[i] = null;
			}
		}
		
		public Triangle(int aC, int bC, int cC) 
		{
			this();
			cornerIndices[0] = aC;
			cornerIndices[1] = bC;
			cornerIndices[2] = cC;
		}

	}
	
	private Triangle meshes[];


	/**
	 * 
	 */
	public Triangulator() 
	{
		
	}

	/**
	 * @param args
	 */
	public static void main(String[] args) 
	{
		

	}

}
