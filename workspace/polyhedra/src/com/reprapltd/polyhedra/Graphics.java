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
	
	private void VisitTrianglesR(Triangle lastTriangle, Triangle tri, Triangulation.CornerList corners, ArrayList<Point3d> stripCornersList,
			ArrayList<Integer> triangleCountsList, boolean sameStrip)
	{
		/*
		 * TODO - this is not quite right.  When running along a strip, the LAST two points added must be the base of the next triangle.
		 */
		
		if(sameStrip)
		{
			if(lastTriangle == null)
			{
				Debug.Error("Graphics.VisitTrianglesR(): previous triangle in a strip was null.", true);
				return;
			}
			boolean notGotOne = true;
			for(int corner = 0; corner < 3; corner++)
			{
				int nextCorner = (corner + 1)%3;
				if(lastTriangle.HasCorners(tri.GetCorner(corner), tri.GetCorner(nextCorner)))
				{
					int newCorner = tri.OppositeCorner(corner, nextCorner);
					Point3d newPoint = corners.GetCorner(newCorner);
					stripCornersList.add(newPoint);
					Integer newCount = triangleCountsList.get(triangleCountsList.size() - 1) + 1;
					triangleCountsList.set(triangleCountsList.size() - 1, newCount);
					notGotOne = false;
					break;
				}
			}
			if(notGotOne)
			{
				Debug.Error("Graphics.VisitTrianglesR(): previous triangle in a strip did not have current triangles edge.", true);
				return;				
			}
		} else
		{
			for(int corner = 0; corner < 3; corner++)
			{
				int cornerListIndex = tri.GetCorner(corner);
				Point3d newPoint = corners.GetCorner(cornerListIndex);
				stripCornersList.add(newPoint);
			}
			triangleCountsList.add(three);
		}

		// Now set me as having been dealt with then recursively visit my neighbours.

		tri.SetVisited();

		if(tri.GetNeighbour(0) != null)
		{
			if(!tri.GetNeighbour(0).Visited())
				VisitTrianglesR(tri, tri.GetNeighbour(0), corners, stripCornersList, triangleCountsList, true);
		} else
		{
			Debug.Error("Graphics.VisitTrianglesR() - 0: triangle with null neighbour found.", true);
		}
		
		if(tri.GetNeighbour(1) != null)
		{
			if(!tri.GetNeighbour(1).Visited())
				VisitTrianglesR(tri, tri.GetNeighbour(1), corners, stripCornersList, triangleCountsList, false);
		} else
		{
			Debug.Error("Graphics.VisitTrianglesR() - 1: triangle with null neighbour found.", true);
		}
		
		if(tri.GetNeighbour(2) != null)
		{
			if(!tri.GetNeighbour(2).Visited())
				VisitTrianglesR(tri, tri.GetNeighbour(2), corners, stripCornersList, triangleCountsList, false);
		} else
		{
			Debug.Error("Graphics.VisitTrianglesR() - 2: triangle with null neighbour found.", true);
		}
			
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
	        plaPts[1] = new Point3d(0,0,0);
	        plaPts[0] = new Point3d(0.2,0,0);
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
    



