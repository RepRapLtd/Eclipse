package com.reprapltd.polyhedra;

import javax.vecmath.Matrix4d;
import javax.vecmath.Point3d;

public class Point3D 
{
	/**
	 * 
	 */
	private double x, y, z;
	
	/**
	 * Default to the origin
	 */
	public Point3D()
	{
		x = 0;
		y = 0;
		z = 0;
	}
	
	/**
	 * Usual constructor
	 * @param a
	 * @param b
	 */
	public Point3D(double a, double b, double c)
	{
		x = a;
		y = b;
		z = c;
	}
	
	/**
	 * Make one from the javax.vecmath equivalent
	 * @param p
	 */
	public Point3D(Point3d p)
	{
		x = p.x;
		y = p.y;
		z = p.z;
	}	
	
	/**
	 * Copy
	 * @param r Rr2Point to copy from
	 */
	public Point3D(Point3D r)
	{
		x = r.x;
		y = r.y;
		z = r.z;
	}
	
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString()
	{
		return Double.toString(x) + " " + Double.toString(y) + " " + Double.toString(z);
	}
	
	/**
	 * Coordinates
	 */
	public double x() { return x; }
	public double y() { return y; }
	public double z() { return z; }
	
	/**
	 * Arithmetic
	 * @return neg of point
	 */
	public Point3D neg()
	{
		return new Point3D(-x, -y, -z);
	}
	
	/**
	 * Go somewhere else
	 * @param m
	 * @return
	 */
	public Point3D transform(Matrix4d m)
	{
		Point3D result = new Point3D(m.m00*x + m.m10*y + m.m20*z + m.m30,
				m.m01*x + m.m11*y + m.m21*z + m.m31,
				m.m02*x + m.m12*y + m.m22*z + m.m32);
		double d = m.m03*x + m.m13*y + m.m23*z + m.m33;
		return Point3D.div(result, d);
	}
	
	
	/**
	 * @param a
	 * @param b
	 * @return a new point based on a vector addition of points a and b
	 */
	public static Point3D add(Point3D a, Point3D b)
	{
		Point3D r = new Point3D(a);
		r.x += b.x;
		r.y += b.y;
		r.z += b.z;
		return r;
	}
	
	/**
	 * @param a
	 * @param b
	 * @return a new point based on a vector subtraction of a - b
	 */
	public static Point3D sub(Point3D a, Point3D b)
	{
		return add(a, b.neg());
	}
	
	
	/**
	 * Scale a point
	 * @param b An R2rPoint
	 * @param factor A scale factor
	 * @return The point Rr2Point scaled by a factor of factor
	 */
	public static Point3D mul(Point3D b, double factor)
	{
		return new Point3D(b.x*factor, b.y*factor, b.z*factor);
	}
	
	/**
	 * @param a
	 * @param b
	 * @return the point Rr2Point scaled by a factor of a
	 */
	public static Point3D mul(double a, Point3D b)
	{
		return mul(b, a);
	}
	
	/**
	 * Downscale a point
	 * @param b An R2rPoint
	 * @param factor A scale factor
	 * @return The point Rr2Point divided by a factor of a
	 */
	public static Point3D div(Point3D b, double factor)
	{
		return mul(b, 1/factor);
	}
	
	/**
	 * Inner product
	 * @param a
	 * @param b
	 * @return The scalar product of the points
	 */
	public static double mul(Point3D a, Point3D b)
	{
		return a.x*b.x + a.y*b.y + a.z*b.z;
	}
	

	/**
	 * Outer product
	 * @param a
	 * @param b
	 * @return oute product
	 */
	public static Point3D op(Point3D a, Point3D b)
	{
		return new Point3D(a.y*b.z - a.z*b.y, a.z*b.x - a.x*b.z,  a.x*b.y - a.y*b.x);
	}
	
	
	/**
	 * Modulus
	 * @return modulus
	 */
	public double mod()
	{
		return Math.sqrt(mul(this, this));
	}
	
	
	/**
	 * Unit length normalization
	 * @return normalized unit lenght 
	 */
	public Point3D norm()
	{
		return div(this, mod());
	}

	
	/**
	 * Squared distance
	 * @param a
	 * @param b
	 * @return squared distance
	 */
	public static double dSquared(Point3D a, Point3D b)
	{
		Point3D c = sub(a, b);
		return mul(c, c);
	}
	
	/**
	 * distance
	 * @param a
	 * @param b
	 * @return distance
	 */
	public static double d(Point3D a, Point3D b)
	{
		return Math.sqrt(dSquared(a, b));
	}
	
	/**
	 * The same, withing tolerance?
	 * @param a
	 * @param b
	 * @param tol_2
	 * @return true if the squared distance between points a and b 
	 * is within tolerance tol_2, otherwise false
	 */
	public static boolean same(Point3D a, Point3D b, double tol_2)
	{
		return dSquared(a, b) < tol_2;
	}

}
