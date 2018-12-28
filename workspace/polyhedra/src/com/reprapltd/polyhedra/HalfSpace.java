package com.reprapltd.polyhedra;

import javax.vecmath.Matrix4d;

import com.reprapltd.polyhedra.Interval;

public class HalfSpace 
{
	
	/**
	 * The half-plane is normal*(x, y) + offset <= 0 
	 */
	private Point3D normal = null; 
	private double offset;
	

	
	/**
	 * Make one from three points in it
	 * @param a
	 * @param b
	 */
	public HalfSpace(Point3D a, Point3D b, Point3D c)
	{
		normal = Point3D.op(Point3D.sub(b, a), Point3D.sub(c, a)).norm();
		offset = -Point3D.mul(normal, a);
	}   
	
	/**
	 * Make one from a normal and one point in it
	 * @param a
	 * @param b
	 */
	public HalfSpace(Point3D n, Point3D a)
	{
		normal = n.norm();
		offset = -Point3D.mul(normal, a);
	}  
	
	/**
	 * Make one from a normal and a distance from the origin
	 * @param a
	 * @param b
	 */
	public HalfSpace(Point3D n, double a)
	{
		normal = n.norm();
		offset = -a;
	}  
	
	/**
	 * Deep copy
	 * @param a
	 */
	public HalfSpace(HalfSpace a)
	{
		normal = new Point3D(a.normal);
		offset = a.offset;
	}
	
	
	/**
	 * Return the plane as a string
	 * @return string representation
	 */
	public String toString()
	{
		return "|" + normal.toString() + ", " + Double.toString(offset) + "|";
	} 
	
	
	/**
	 * Get the components
	 * @return components?
	 */
	public Point3D normal() { return normal; }
	public double offset() { return offset; }
	
	/**
	 * Is another plane the same within a tolerance?
	 * @param a
	 * @param b
	 * @param tolerance
	 * @return 0 if the distance between halfplane a and b is less then the tolerance, -1 if one
	 * is the complement of the other within the tolerance, otherwise 1
	 */
	public static int same(HalfSpace a, HalfSpace b, double tolerance)
	{
		if(a == b)
			return 0;
		
		int result = 0;
		if(Math.abs(a.normal.x() - b.normal.x()) > tolerance)
		{
			if(Math.abs(a.normal.x() + b.normal.x()) > tolerance)
				return 1;
			result = -1;
		}
		if(Math.abs(a.normal.y() - b.normal.y()) > tolerance)
		{
			if(Math.abs(a.normal.y() + b.normal.y()) > tolerance || result != -1)
				return 1;
		}
		double rms = Math.sqrt((a.offset*a.offset + b.offset*b.offset)*0.5);
		if(Math.abs(a.offset - b.offset) > tolerance*rms)
		{
			if(Math.abs(a.offset + b.offset) > tolerance*rms || result != -1)
				return 1;
		}
		
		return result;
	}

	
	/**
	 * Change the sense
	 * @return complent of half plane
	 */
	public HalfSpace complement()
	{
		HalfSpace r = new HalfSpace(this);
		r.normal = r.normal.neg();
		r.offset = -r.offset;
		return r;
	}
	
	/**
	 * Move somewhere else.  NOTE THIS EXPECTS THE INVERSE OF THE TRANSFORM
	 * @param iM
	 * @return
	 */
	public HalfSpace transform(Matrix4d iM)
	{
		Point3D n = new Point3D(iM.m00*normal.x() + iM.m10*normal.y() + iM.m20*normal.z() + iM.m30*offset,
				iM.m01*normal.x() + iM.m11*normal.y() + iM.m21*normal.z() + iM.m31*offset,
				iM.m02*normal.x() + iM.m12*normal.y() + iM.m22*normal.z() + iM.m32*offset);
		double o = iM.m03*normal.x() + iM.m13*normal.y() + iM.m23*normal.z() + iM.m33*offset;
		double d = n.mod();
		HalfSpace result = new HalfSpace(this);
		result.normal = Point3D.div(n, d);
		result.offset = o/d;
		return result;
	}
	
	/**
	 * Move
	 * @param d
	 * @return offset halfplane
	 */
	public HalfSpace offset(double d)
	{
		HalfSpace r = new HalfSpace(this);
		r.offset = r.offset - d;
		return r;
	}
	
	
	/**
	 * Find the potential value of a point
	 * @param p
	 * @return potential value of point p
	 */
	public double value(Point3D p)
	{
		return offset + Point3D.mul(normal, p);
	}
	
	
	/**
	 * Find the potential interval of a box
	 * @param b
	 * @return potential interval of box b
	 */
	public Interval value(Box b)
	{
		Interval x = Interval.mul(b.x(), normal.x());
		Interval y = Interval.mul(b.y(), normal.y());
		Interval z = Interval.mul(b.z(), normal.z());
		return Interval.add(Interval.add(Interval.add(x, y), z), offset);
	}

}
