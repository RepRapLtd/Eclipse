package com.reprapltd.polyhedra;

import javax.media.j3d.*;
import javax.swing.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point3d;

import com.reprapltd.polyhedra.Triangulation.Triangle;
import com.sun.j3d.utils.behaviors.mouse.*;
import com.sun.j3d.utils.universe.SimpleUniverse;

import java.awt.*;
import java.awt.geom.*;
import java.util.ArrayList;
    
/**
 * Plot a bunch of stuff in 3D.
 * 
 *     Thanks to: https://stackoverflow.com/questions/12313917/how-to-draw-simple-3d-pointsx-y-z-in-java-using-java3d-api
 *     which saved a lot of document reading...
 *     
 * @author ensab
 *
 */
public final class Graphics extends JPanel
{
	/**
	 * Recursive function to walk over all the triangles in one shell and record them
	 * in triangle strips for display.
	 * 
	 * @param out
	 */
	private Point3d[] stripCorners = null;
	private int[] triangleCounts = null;
	private final Integer three = 3;
	
	/**
	 * Small class to define a triangle strip.
	 * 
	 * A strip is completely defined by its start triangle,
	 * the triangle visited next after that, its end triangle,
	 * and the first and second corners at the edge of the start triangle
	 * from which the next triangle was visited.
	 *  
	 * @author ensab
	 *
	 */
	private class StripDefine
	{
		private Triangle start, next, end;
		private int firstCorner;
		private int secondCorner;
		
		private StripDefine(Triangle s)
		{
			start = null;
			next = null;
			end = null;
			firstCorner = -1;
			secondCorner = -1;
			
			// If this triangle has already been visited, it can't contribute to the strips again.
			
			if(s.Visited())
			{
				return;
			}
			
			// No - s hasn't been visited; it can start a strip.
			
			start = s;
			
			// Find an unvisited neighbour
			
			int nIndex = 0;
			int free = -1;
			while(nIndex < 3 && free < 0)
			{
				Triangle neighbour = start.GetNeighbour(nIndex);
				if(!neighbour.Visited())
					free = nIndex;
				nIndex++;
			}
			
			// If there are no unvisited neighbours, this triangle is a singleton.
			
			if(free < 0)
				return;
			
			// Record which triangle we will visit next
			
			next = start.GetNeighbour(free);
			
			// Draw a little diagram to show how this works...
			
			firstCorner = (free + 1)%3;
			secondCorner = (free + 2)%3;
		}
		
		private void SetEnd(Triangle e)
		{
			end = e;
		}
	}
	
	/**
	 * This starts at Triangle start and walks across the triangulation adding a strip of triangles to the list for
	 * display.
	 * 
	 * @param start
	 * @param corners
	 * @param stripCornersList
	 * @param triangleCountsList
	 * @return
	 */
	private StripDefine WalkStrip(Triangle start, Triangulation.CornerList corners, ArrayList<Point3d> stripCornersList, ArrayList<Integer> triangleCountsList)
	{
		StripDefine strip = new StripDefine(start);
		
		// Anything to do?
		
		if(strip.start == null)
			return strip;
		
		// Yes - add the first triangle's corners to the strip
		// Note we have to do this so the last corner added is 
		// strip.secondCorner if there is a next triangle to visit.
		
		int corner = 0;
		if(strip.next != null)
			corner = (strip.secondCorner + 1)%3;
		
		for(int count = 0; count < 3; count++)
		{
			Point3d newPoint = corners.GetCorner(strip.start.GetCorner(corner));
			stripCornersList.add(newPoint);
			corner++;
			if(corner >= 3)
				corner = 0;
		}
		triangleCountsList.add(three);
		strip.start.SetVisited();
		
		Triangle next = strip.next;
		
		// Anything more to do?
		
		if(next == null)
		{
			strip.SetEnd(strip.start);
			return strip;
		}
		
		// We need the actual corners in the corner list so we can find them
		// in Triangle next's list of 3 corners.
		
		int c1 = strip.start.GetCorner(strip.firstCorner);
		int c2 = strip.start.GetCorner(strip.secondCorner);
		
		while(!next.Visited())
		{
			// Find the edge c1-c2 in next's edges
			
			c1 = next.FindCorner(c1);
			if(c1 < 0)
			{
				Debug.Error("Graphics.WalkStrip(): previous triangle in a strip did not share a corner 1 with next triangle.", true);
				return strip;
			}
			c2 = next.FindCorner(c2);
			if(c2 < 0)
			{
				Debug.Error("Graphics.WalkStrip(): previous triangle in a strip did not share a corner 2 with next triangle.", true);
				return strip;
			}
			
			// We need to add the point opposite that edge
			
			int nextPointIndex = next.Opposite(c1, c2);
			Point3d newPoint = corners.GetCorner(next.GetCorner(nextPointIndex));
			stripCornersList.add(newPoint);
			Integer newCount = triangleCountsList.get(triangleCountsList.size() - 1) + 1;
			triangleCountsList.set(triangleCountsList.size() - 1, newCount);
			next.SetVisited();
			
			// Update the edge to be the one between c2 and the newly added point from triangle next.
			
			c1 = c2;
			c2 = nextPointIndex;
			
			// The next triangle to visit is the neighbour of next that is opposite that edge.
			
			Triangle n = next.GetNeighbour(next.Opposite(c1, c2));
			
			// Find the actual corners in the corner list for the edge finder at the start of this loop.
			
			c1 = next.GetCorner(c1);
			c2 = next.GetCorner(c2);
			
			// Continuously update the end of the strip
			
			strip.SetEnd(next);
			next = n;
		}
		
		return strip;
	}
	
	
	/**
	 * This walks along a strip that has been added finding unvisited neighbours of the triangles in it, and
	 * walks along each of them in turn.
	 * 
	 * The walking logic here is (and MUST be) identical to that in WalkStrip().
	 * 
	 * @param strip
	 * @param corners
	 * @param stripCornersList
	 * @param triangleCountsList
	 */
	private void RecurseOnStrip(StripDefine strip, Triangulation.CornerList corners, ArrayList<Point3d> stripCornersList, ArrayList<Integer> triangleCountsList)
	{
		// Anything to do?
		
		if(strip.start == null)
			return;
		
		// Yes - run WalkStrip on any unvisited neighbours of the start triangle.
		
		StripDefine newStrip;
		
		for(int tri = 0; tri < 3; tri++)
		{
			Triangle neighbour = strip.start.GetNeighbour(tri);
			newStrip = WalkStrip(neighbour, corners, stripCornersList, triangleCountsList);
			RecurseOnStrip(newStrip, corners, stripCornersList, triangleCountsList);
		}
		
		Triangle next = strip.next;
		
		// Anything more to do?
		
		if(next == null)
			return;
		
		// We need the actual corners in the corner list so we can find them
		// in Triangle next's list of 3 corners.
		
		int c1 = strip.start.GetCorner(strip.firstCorner);
		int c2 = strip.start.GetCorner(strip.secondCorner);
		
		while(true) // Hmmm....
		{
			// Find the edge c1-c2 in next's edges
			
			c1 = next.FindCorner(c1);
			if(c1 < 0)
			{
				Debug.Error("Graphics.RecurseOnStrip(): previous triangle in a strip did not share a corner 1 with next triangle.", true);
				return;
			}
			c2 = next.FindCorner(c2);
			if(c2 < 0)
			{
				Debug.Error("Graphics.RecurseOnStrip(): previous triangle in a strip did not share a corner 2 with next triangle.", true);
				return;
			}
			
			// We need the point opposite that edge to find the next triangle
			
			int nextPointIndex = next.Opposite(c1, c2);
			
			// Visit any of next's neighbours that haven't been visited
			
			for(int tri = 0; tri < 3; tri++)
			{
				Triangle neighbour = next.GetNeighbour(tri);
				newStrip = WalkStrip(neighbour, corners, stripCornersList, triangleCountsList);
				RecurseOnStrip(newStrip, corners, stripCornersList, triangleCountsList);
			}
			
			// If we are at the end, there is no more to do.
			
			if(next == strip.end)
				return;
			
			// Update the edge to be the one between c2 and the newly added point from triangle next.
			
			c1 = c2;
			c2 = nextPointIndex;
			
			// The next triangle to visit is the neighbour of next that is opposite that edge.
			
			Triangle n = next.GetNeighbour(next.Opposite(c1, c2));
			
			// Find the actual corners in the corner list for the edge finder at the start of this loop.
			
			c1 = next.GetCorner(c1);
			c2 = next.GetCorner(c2);
			next = n;
		}
	}
	
	private void VisitTrianglesR(Triangle start, Triangulation.CornerList corners, ArrayList<Point3d> stripCornersList,
			ArrayList<Integer> triangleCountsList)
	{

			
	}
	
	/**
	 * Non-recursive function to plot all the shells by calling the recursive
	 * function above.
	 * 
	 * @param out
	 */
	private void VisitTriangles(Triangulation triangulation)
	{
		ArrayList<Point3d> stripCornersList = new ArrayList<Point3d>();
		ArrayList<Integer> triangleCountsList = new ArrayList<Integer>();
		
		for(int shell = 0; shell < triangulation.Shells().size(); shell++)
		{
			triangulation.Shells().get(shell).Reset(); // Shouldn't be needed
			VisitTrianglesR(null, triangulation.Shells().get(shell), triangulation.Corners(), stripCornersList, triangleCountsList, false);
			triangulation.Shells().get(shell).Reset();
		}
		
		stripCorners = new Point3d[stripCornersList.size()];
		triangleCounts = new int[triangleCountsList.size()];
		
		for(int corner = 0; corner < stripCornersList.size(); corner++)
		{
			stripCorners[corner] = stripCornersList.get(corner);
		}
		
		for(int countIndex = 0; countIndex < triangleCountsList.size(); countIndex++)
		{
			triangleCounts[countIndex] = triangleCountsList.get(countIndex);
		}
		
	}

	public Graphics() 
	{
		setLayout(new BorderLayout());
		GraphicsConfiguration gc=SimpleUniverse.getPreferredConfiguration();
		Canvas3D canvas3D = new Canvas3D(gc);//See the added gc? this is a preferred config
		add("Center", canvas3D);

		BranchGroup scene = createSceneGraph();
		scene.compile();

		// SimpleUniverse is a Convenience Utility class
		SimpleUniverse simpleU = new SimpleUniverse(canvas3D);


		// This moves the ViewPlatform back a bit so the
		// objects in the scene can be viewed.
		simpleU.getViewingPlatform().setNominalViewingTransform();

		simpleU.addBranchGraph(scene);
	}
	
	public BranchGroup createSceneGraph() 
	{
	       BranchGroup lineGroup = new BranchGroup();
	        Appearance app = new Appearance();
	        ColoringAttributes ca = new ColoringAttributes(new Color3f(204.0f, 204.0f,204.0f), ColoringAttributes.SHADE_FLAT);
	        app.setColoringAttributes(ca);
	        
	        Point3d[] plaPts = new Point3d[8];
	        plaPts[0] = new Point3d(0.2,0,0);
	        plaPts[1] = new Point3d(0,0,0);
	        plaPts[2] = new Point3d(0.2,0.1,0);
	        plaPts[3] = new Point3d(0,0.2,0);
	        plaPts[4] = new Point3d(0.3,0.3,0);
	        plaPts[5] = new Point3d(0.3,0,0);
	        plaPts[6] = new Point3d(0.4,0,0);
	        plaPts[7] = new Point3d(0.4,0.1,0);	        
	        int[]intArr=new int[2];
	        intArr[0]=5;
	        intArr[1]=3;

	        TriangleStripArray pla = new TriangleStripArray(8, GeometryArray.COLOR_3|GeometryArray.NORMALS|GeometryArray.COORDINATES, intArr);
	        pla.setCoordinates(0, plaPts);
	        
	        PointAttributes a_point_just_bigger=new PointAttributes();
	        a_point_just_bigger.setPointSize(10.0f);//10 pixel-wide point
	        a_point_just_bigger.setPointAntialiasingEnable(true);//now points are sphere-like(not a cube)
	        app.setPointAttributes(a_point_just_bigger);
	        PolygonAttributes la=new PolygonAttributes();
	        la.setPolygonMode(PolygonAttributes.POLYGON_FILL);
	        la.setCullFace(PolygonAttributes.CULL_NONE);
	        app.setPolygonAttributes(la);
	        
	        Material matt=new Material();
	        matt.setAmbientColor(new Color3f(1,1,1));
	        matt.setDiffuseColor(new Color3f(0.5f,0.5f,0.7f));
	        matt.setEmissiveColor(new Color3f(0.2f,0.2f,0.3f));
	        matt.setShininess(0.5f);
	        matt.setSpecularColor(new Color3f(0.4f,0.6f,0.9f));
	        matt.setLightingEnable(true);
	        app.setMaterial(matt);
	        
	        RenderingAttributes ra=new RenderingAttributes();
	        ra.setIgnoreVertexColors(true);
	        app.setRenderingAttributes(ra);
	        
	        Shape3D plShape = new Shape3D(pla, app);

	        TransformGroup objRotate = new TransformGroup();

	        MouseRotate mr=new MouseRotate();
	        mr.setTransformGroup(objRotate);
	        mr.setSchedulingBounds(new BoundingSphere());
	        lineGroup.addChild(mr);
	        MouseZoom mz=new MouseZoom();
	        mz.setTransformGroup(objRotate);
	        mz.setSchedulingBounds(new BoundingSphere());
	        lineGroup.addChild(mz);
	        MouseTranslate msl=new MouseTranslate();
	        msl.setTransformGroup(objRotate);
	        msl.setSchedulingBounds(new BoundingSphere());
	        lineGroup.addChild(msl);


	        objRotate.setCapability(TransformGroup.ALLOW_TRANSFORM_WRITE);
	        objRotate.addChild(plShape);
	        AmbientLight al=new AmbientLight();
	      //  al.addScope(objRotate);
	        al.setBounds(new BoundingSphere());
	        al.setEnable(true);
	        al.setColor(new Color3f(1,1,1));

	        lineGroup.addChild(objRotate);
	        lineGroup.addChild(al);
	        return lineGroup;
	}

	public static void main(String[] args) 
	{
		JFrame frame = new JFrame();
		frame.add(new JScrollPane(new Graphics()));
		frame.setSize(800, 800);
		frame.setVisible(true);
		frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	}

}
    



