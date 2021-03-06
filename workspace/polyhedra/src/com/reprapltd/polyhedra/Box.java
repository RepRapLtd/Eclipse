package com.reprapltd.polyhedra;

import com.reprapltd.polyhedra.Interval;
import com.reprapltd.polyhedra.Point3D;

public class Box 
{
	/**
	 * Compass directions
	 */
	public static final byte rr_N = 1;
	public static final byte rr_E = 2;
	public static final byte rr_S = 4;
	public static final byte rr_W = 8;
	public static final byte rr_U = 16;
	public static final byte rr_D = 32;
	/**
	 * X interval
	 */
	private Interval x = null;
	
	/**
	 * Y interval
	 */
	private Interval y = null;
	
	/**
	 * Z interval
	 */
	private Interval z = null;
	
	/**
	 * Anyone home?
	 */
	private boolean empty;
	

	
	/**
	 * Default is empty
	 */ 
	public Box()
	{
		empty = true;
	}
	
	
	/**
	 * Copy constructor
	 */
	public Box(Box b)
	{
		x = new Interval(b.x);
		y = new Interval(b.y);
		z = new Interval(b.z);
		empty = b.empty;
	}
	
	/**
	 * Make from any diagonal corners
	 * @param sw
	 * @param ne
	 */
	public Box(Point3D a, Point3D b)
	{
		x = new Interval(Math.min(a.x(), b.x()), Math.max(a.x(), b.x()));
		y = new Interval(Math.min(a.y(), b.y()), Math.max(a.y(), b.y()));
		z = new Interval(Math.min(a.z(), b.z()), Math.max(a.z(), b.z()));
		empty = x.empty() || y.empty() || z.empty();
	}
	
	/**
	 * Make from X and Y intervals
	 * @param sw
	 * @param ne
	 */
	public Box(Interval xi, Interval yi, Interval zi)
	{
		x = new Interval(xi);
		y = new Interval(yi);
		z = new Interval(zi);
		empty = x.empty() || y.empty() || z.empty();
	}
	
	
	/**
	 * @return Return the x interval
	 */
	public Interval x() { return x; }
	
	
	/**
	 * @return Return the y interval
	 * 
	 */
	public Interval y() { return y; }
	
	/**
	 * @return Return the z interval
	 * 
	 */
	public Interval z() { return z; }
	
	/**
	 * 
	 * @return
	 */
	public boolean empty() { return empty; }
	
	
	/**
	 * Expand the box to incorporate another box or a point
	 * @param a
	 */
	public void expand(Box a)
	{
		if(a.empty)
			return;
		if(empty)
		{
			empty = false;
			x = new Interval(a.x);
			y = new Interval(a.y);
			z = new Interval(a.z);
		} else
		{
			x.expand(a.x);
			y.expand(a.y);
			z.expand(a.z);
		}
	}
	
	/**
	 * Shrink or grow by a given distance
	 * @param dist
	 * @return
	 */
	public Box offset(double dist)
	{
		return new Box(new Interval(x.low() - dist, x.high() + dist), 
				new Interval(y.low() - dist, y.high() + dist),
				new Interval(z.low() - dist, z.high() + dist));
	}
	
	/**
	 * Move somewhere else
	 * @param p
	 * @return
	 */
	public Box translate(Point3D p)
	{
		return new Box(new Interval(x.low() + p.x(), x.high() + p.x()), 
				new Interval(y.low() + p.y(), y.high() + p.y()),
				new Interval(z.low() + p.z(), z.high() + p.z()));
	}
	
	/**
	 * @param a
	 */
	public void expand(Point3D a)
	{
		if(empty)
		{
			empty = false;
			x = new Interval(a.x(), a.x());
			y = new Interval(a.y(), a.y());
			z = new Interval(a.z(), a.z());
		} else
		{
			x.expand(a.x());
			y.expand(a.y());
			z.expand(a.z());
		}
	}
	
	/**
	 * Corner points and center
	 * @return NE cornerpoint
	 */
	public Point3D nel()
	{
		return new Point3D(x.high(), y.high(), z.low());
	}
	public Point3D neh()
	{
		return new Point3D(x.high(), y.high(), z.high());
	}
	
	/**
	 * @return SW cornerpoint
	 */
	public Point3D swl()
	{
		return new Point3D(x.low(), y.low(), z.low());
	}
	public Point3D swh()
	{
		return new Point3D(x.low(), y.low(), z.high());
	}
	
	/**
	 * @return SE cornerpoint
	 */
	public Point3D sel()
	{
		return new Point3D(x.high(), y.low(), z.low());
	}
	public Point3D seh()
	{
		return new Point3D(x.high(), y.low(), z.high());
	}
	
	/**
	 * @return NW cornerpoint
	 */
	public Point3D nwl()
	{
		return new Point3D(x.low(), y.high(), z.low());
	}   
	public Point3D nwh()
	{
		return new Point3D(x.low(), y.high(), z.high());
	}   
	
	/**
	 * @return Centre point
	 */
	public Point3D centre()
	{
		return new Point3D(x.cen(), y.cen(), z.cen());
	}	
	/**
	 * Scale the box by a factor about its center
	 * @param f
	 * @return scaled box object
	 */
	public Box scale(double f)
	{
		Box r = new Box();
		if(empty)
			return r;
		f = 0.5*f;
		Point3D p = new Point3D(x.length()*f, y.length()*f, z.length()*f);
		Point3D c = centre();
		r.expand(Point3D.add(c, p));
		r.expand(Point3D.sub(c, p));
		return r;
	}
	
	/**
	 * Convert to a string
	 */
	public String toString()
	{
		if(empty)
			return "<empty>";
		return "<BOX x:" + x.toString() + ", y:" + y.toString() + ", z:" + z.toString() + ">";
	}
	
	/**
	 * Squared diagonal
	 * @return squared diagonal of the box
	 */
	public double dSquared()
	{
		if(empty)
			return 0;
		return Point3D.dSquared(swl(), neh());
	}

//	/**
//	 * Squared distance to a point
//	 * @return minimal squared distance to a point from one of the corners of the box
//	 */
//	public double dSquared(Point3D p)
//	{
//		if(empty)
//			return Double.POSITIVE_INFINITY;
//		byte b = pointRelative(p);
//		double d1 = Double.POSITIVE_INFINITY, d2 = Double.POSITIVE_INFINITY;
//		switch(b)
//		{
//		case 0:
//			return 0;
//			
//		case 1:
//			d1 = Point2D.dSquared(p, nw());
//			d2 = Point2D.dSquared(p, ne());
//			break;
//			
//		case 2:
//			d1 = Point2D.dSquared(p, ne());
//			d2 = Point2D.dSquared(p, se());
//			break;
//			
//		case 3:
//			return Point2D.dSquared(p, ne());
//			
//		case 4:
//			d1 = Point2D.dSquared(p, sw());
//			d2 = Point2D.dSquared(p, se());
//			break;
//			
//		case 6:
//			return Point2D.dSquared(p, se());
//			
//		case 8:
//			d1 = Point2D.dSquared(p, sw());
//			d2 = Point2D.dSquared(p, nw());
//			break;
//			
//		case 9:
//			return Point2D.dSquared(p, nw());
//			
//		case 12:
//			return Point2D.dSquared(p, sw());
//			
//		default:
//			Debug.e("RrRectangle.dSquared(): dud value from point_relative()!");	
//		}
//		
//		return Math.min(d1, d2);
//
//	}
	
	/**
	 * Where is a point relative to a box?
	 * @param p 
	 * @return relative position of a point p (E, W, N or S)
	 */
	public byte pointRelative(Point3D p)
	{
		byte result = 0;
		if(p.x() >= x.high())
			result |= rr_E;
		if(p.x() < x.low())
			result |= rr_W;
		if(p.y() >= y.high())
			result |= rr_N;
		if(p.y() < y.low())
			result |= rr_S;
		if(p.z() >= z.high())
			result |= rr_U;
		if(p.z() < z.low())
			result |= rr_D;
		return result;
	}
	
	/**
	 * Intersection
	 * @param a
	 * @param b
	 * @return
	 */
	public static Box intersection(Box a, Box b)
	{
		if(a.empty)
			return a;
		if(b.empty)
			return b;
		return new Box(Interval.intersection(a.x, b.x), Interval.intersection(a.y, b.y), Interval.intersection(a.z, b.z));	
	}
	
	/**
	 * Union
	 * @param a
	 * @param b
	 * @return
	 */
	public static Box union(Box a, Box b)
	{
		if(a.empty)
			return b;
		if(b.empty)
			return a;
		return new Box(Interval.union(a.x, b.x), Interval.union(a.y, b.y), Interval.union(a.z, b.z));	
	}
}

