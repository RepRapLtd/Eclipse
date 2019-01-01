package com.reprapltd.polyhedra;

import javax.media.j3d.AmbientLight;
import javax.media.j3d.Appearance;
import javax.media.j3d.BoundingSphere;
import javax.media.j3d.BranchGroup;
import javax.media.j3d.Canvas3D;
import javax.media.j3d.ColoringAttributes;
import javax.media.j3d.GeometryArray;
import javax.media.j3d.Material;
import javax.media.j3d.PointArray;
import javax.media.j3d.PointAttributes;
import javax.media.j3d.PolygonAttributes;
import javax.media.j3d.RenderingAttributes;
import javax.media.j3d.Shape3D;
import javax.media.j3d.TransformGroup;
import javax.media.j3d.TriangleStripArray;
import javax.swing.*;
import javax.vecmath.Color3f;
import javax.vecmath.Point3f;

import com.reprapltd.polyhedra.Triangulation.Triangle;
import com.sun.j3d.utils.behaviors.mouse.MouseRotate;
import com.sun.j3d.utils.behaviors.mouse.MouseTranslate;
import com.sun.j3d.utils.behaviors.mouse.MouseZoom;
import com.sun.j3d.utils.universe.SimpleUniverse;

import java.awt.*;
import java.awt.geom.*;
    
    
public final class Graphics extends JPanel
{


	int s = 0, count = 0;
	
	/**
	 * Recursive function to walk over all the triangles in one shell and output
	 * them to a .ply file.
	 * 
	 * @param out
	 */
	private static void VisitTrianglesR(Triangle tri, Triangulation triangulation, Graphics2D g2)
	{
		Path2D path = new Path2D.Double();

		int i0 = tri.GetCorner(0);
		int i1 = tri.GetCorner(1);
		int i2 = tri.GetCorner(2);
		Point3D p0 = triangulation.Corners().GetCorner(i0);
		Point3D p1 = triangulation.Corners().GetCorner(i1);
		Point3D p2 = triangulation.Corners().GetCorner(i2);

		path.moveTo(p0.x(), p0.y());
		path.lineTo(p1.x(), p1.y());
		path.lineTo(p2.x(), p2.y());
		path.closePath();
		g2.draw(path);

		// Now recursively visit my neighbours.

		tri.SetVisited();
		for(int i = 0; i < 3; i++)
		{
			if(tri.GetNeighbour(i) != null)
			{
				if(!tri.GetNeighbour(i).Visited())
					VisitTrianglesR(tri.GetNeighbour(i), triangulation, g2);
			} else
				Debug.Error("Graphics.VisitTrianglesR(): triangle with null neighbour found.", true);
		}			
	}
	
	/**
	 * Brief non-recursive function to plot all the shells by calling the recursive
	 * function above.
	 * 
	 * @param out
	 */
	private static void VisitTriangles(Triangulation triangulation, Graphics2D g2)
	{
		for(int i = 0; i < triangulation.Shells().size(); i++)
		{
			triangulation.Shells().get(i).Reset(); // Shouldn't be needed
			VisitTrianglesR(triangulation.Shells().get(i), triangulation, g2);
			triangulation.Shells().get(i).Reset();
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
	        Point3f[] plaPts = new Point3f[5];
	        for (int i = 0; i < 2; i++) 
	        {
	            for (int j = 0; j <2; j++) 
	            {
	                plaPts[count] = new Point3f(i/10.0f,j/10.0f,0);
	                count++;
	            }
	        }
	        plaPts[count] = new Point3f(3.0f/10.0f,2.0f/10.0f,0);
	        int[]intArr=new int[5];
	        intArr[0]=3;intArr[1]=4;intArr[2]=4;intArr[3]=3;intArr[4]=3;

	        TriangleStripArray pla =new TriangleStripArray(20, GeometryArray.COLOR_3|GeometryArray.NORMALS|GeometryArray.COORDINATES,intArr);
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
    



