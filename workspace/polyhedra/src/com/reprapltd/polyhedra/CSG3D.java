package com.reprapltd.polyhedra;

import java.util.ArrayList;
import javax.vecmath.Matrix4d;

import com.reprapltd.polyhedra.CSGOp;
import com.reprapltd.polyhedra.Debug;
import com.reprapltd.polyhedra.Interval;

public class CSG3D 
{
	
	/**
	 * Universal set 
	 */
	private static final CSG3D u = new CSG3D(true);  
	
	/**
	 * Null set  
	 */
	private static final CSG3D n = new CSG3D(false); 
	
	/**
	 * Leaf half plane 
	 */
	private HalfSpace hp = null;
	
	/**
	 * Type of set
	 */
	private CSGOp op;
	
	/**
	 * Non-leaf child operands 
	 */
	private CSG3D c1 = null;
	private CSG3D c2 = null; 
	
	/**
	 * The complement (if there is one) 
	 */
	private CSG3D comp = null;        
	
	/**
	 * How much is in here (leaf count)?
	 */
	private int complexity;
		
	/**
	 * Make a leaf from a single half-plane
	 * @param h
	 */
	public CSG3D(HalfSpace h)
	{
		hp = new HalfSpace(h);
		op = CSGOp.LEAF;
		c1 = null;
		c2 = null;
		comp = null;
		complexity = 1;
	}
	
	/**
	 * One off constructor for the universal and null sets
	 * @param b
	 */
	private CSG3D(boolean b)
	{
		hp = null;
		if(b)
			op = CSGOp.UNIVERSE;
		else
			op = CSGOp.NULL;
		c1 = null;
		c2 = null;
		comp = null;   // Resist temptation to be clever here
		complexity = 0;
	}
	
	/**
	 * Universal or null set
	 * @return universal or null set
	 */
	public static CSG3D universe()
	{
		return u;
	}
	
	/**
	 * @return nothing/null set
	 */
	public static CSG3D nothing()
	{
		return n;
	}
	
	/**
	 * Deep copy
	 * @param c
	 */
	public CSG3D(CSG3D c)
	{
		if(c == u || c == n)
			Debug.Error("CSG3D deep copy: copying null or universal set.", true);
		
		if(c.hp != null)
			hp = new HalfSpace(c.hp);
		else
			hp = null;
		
		if(c.c1 != null)
			c1 = new CSG3D(c.c1);
		else
			c1 = null;
		
		if(c.c2 != null)
			c2 = new CSG3D(c.c2);
		else
			c2 = null;
		
		comp = null;  // This'll be built if it's needed
		
		op = c.op;
		complexity = c.complexity;
	}
	
	/**
	 * Get children, operator etc
	 * @return children
	 */
	public CSG3D c_1() { return c1; }
	public CSG3D c_2() { return c2; }
	public CSGOp operator() { return op; }
	public HalfSpace hSpace() { return hp; }
	public int complexity() { return complexity; }
	
	/**
	 * Convert to a string
	 * @param result
	 * @param white
	 * @return string representation
	 */
	private String toString_r(String result, String white)
	{
		switch(op)
		{
		case LEAF:
			result = result + white + hp.toString() + "\n";
			break;
			
		case NULL:
			result = result + white + "0\n";
			break;
			
		case UNIVERSE:
			result = result + white + "U\n";
			break;
			
		case UNION:
			result = result + white + "+\n";
			white = white + " ";
			result = c1.toString_r(result, white);
			result = c2.toString_r(result, white);
			break;
			
		case INTERSECTION:
			result = result + white + "&\n";
			white = white + " ";
			result = c1.toString_r(result, white);
			result = c2.toString_r(result, white);
			break;
			
		default:
			Debug.Error("toString_r(): invalid operator.", true);
		}
		return result;
	}
	
	/* (non-Javadoc)
	 * @see java.lang.Object#toString()
	 */
	public String toString()
	{
		String result = "3D CSG3D: complexity = " + 
			Integer.toString(complexity) + "\n";
		result = toString_r(result, " ");
		return result;
	}
	
	/**
	 * Private constructor for common work setting up booleans
	 * @param a
	 * @param b
	 */
	private CSG3D(CSG3D a, CSG3D b)
	{
		hp = null;
		comp = null;
		if(a.complexity <= b.complexity) // So we know the 1st child is the simplest
		{
			c1 = a;
			c2 = b;
		} else
		{
			c1 = b;
			c2 = a;
		}
		complexity = c1.complexity + c2.complexity;
	}
	
	/**
	 * Boolean operations, with de Morgan simplifications
	 * @param a
	 * @param b
	 * @return union of passed CSG objects a and b
	 */
	public static CSG3D union(CSG3D a, CSG3D b)
	{
		if(a == b)
			return a;
		if(a.op == CSGOp.NULL)
			return b;
		if(b.op == CSGOp.NULL)
			return a;
		if((a.op == CSGOp.UNIVERSE) || (b.op == CSGOp.UNIVERSE))
			return universe();
		
		if(a.comp != null && b.comp != null)
			if(a.comp == b)
				return universe();
		
		CSG3D r = new CSG3D(a, b);
		r.op = CSGOp.UNION;
		return r;
	}
	
	/**
	 * Boolean operation to perform an intersection
	 * @param a
	 * @param b
	 * @return intersection of passed CSG objects a and b
	 */
	public static CSG3D intersection(CSG3D a, CSG3D b)
	{
		if(a == b)
			return a;
		if(a.op == CSGOp.UNIVERSE)
			return b;
		if(b.op == CSGOp.UNIVERSE)
			return a;
		if((a.op == CSGOp.NULL) || (b.op == CSGOp.NULL))
			return nothing();
		
		if(a.comp != null && b.comp != null)
			if(a.comp == b)
				return nothing();
		
		CSG3D r = new CSG3D(a, b);
		r.op = CSGOp.INTERSECTION;
		return r;
	}
	
	/**
	 * Lazy evaluation for complement.
	 * @return complement
	 */
	public CSG3D complement()
	{		
		if(comp != null)
			return comp;
		
		CSG3D result;
		
		switch(op)
		{
		case LEAF:
			result = new CSG3D(hp.complement());
			break;
			
		case NULL:
			return universe();
			
		case UNIVERSE:
			return nothing();
			
		case UNION:
			result = intersection(c1.complement(), c2.complement());
			break;
			
		case INTERSECTION:
			result = union(c1.complement(), c2.complement());
			break;
			
		default:
			Debug.Error("complement(): invalid operator.", true);
			return nothing();
		}
		
		// Remember, so we don't have to do it again.
		// (I do hope that the Java garbage collector is up to 
		// spotting this deadly embrace, or we - I mean it - has
		// a memory leak.)
		
		comp = result;
		result.comp = this;
		
		return comp;
	}
	
	/**
	 * Move somewhere else.  Note - this expects the inverse transform
	 * @return
	 */
	private CSG3D xform(Matrix4d iM)
	{				
		CSG3D result;
		
		switch(op)
		{
		case LEAF:
			result = new CSG3D(hp.transform(iM));
			break;
			
		case NULL:
			return nothing();
			
		case UNIVERSE:
			return universe();
			
		case UNION:
			result = union(c1.xform(iM), c2.xform(iM));
			break;
			
		case INTERSECTION:
			result = intersection(c1.xform(iM), c2.xform(iM));
			break;
			
		default:
			Debug.Error("transform(): invalid operator.", true);
			return nothing();
		}
		
		return result;
	}
	
	/**
	 * Go somewhere else
	 * @param m
	 * @return
	 */
	public CSG3D transform(Matrix4d m)
	{
		Matrix4d iM = new Matrix4d(m);
		iM.invert();
		return xform(iM);
	}
	
	/**
	 * Set difference is intersection with complement
	 * @param a
	 * @param b
	 * @return set difference as CSG object based on input CSG objects a and b
	 */		
	public static CSG3D difference(CSG3D a, CSG3D b)
	{
		return intersection(a, b.complement());
	}
	

	

	
	/**
	 * Run through a GSG expression looking at its leaves and return
	 * a list of the distinct leaves.  Note: leaf and leaf.complement() are
	 * not considered distinct.
	 * Recursive internal call.
	 * @param expression
	 * @return
	 */
	private void uniqueList_r(ArrayList<CSG3D> list)
	{
		switch(op)
		{
		case LEAF:
			CSG3D entry;
			for(int i = 0; i < list.size(); i++)
			{
				entry = list.get(i);
				if(this == entry || complement() == entry)
					return;
			}
			list.add(this);
			break;

		case NULL:			
		case UNIVERSE:
			Debug.Error("uniqueList_r: null or universe at a leaf.", true);
			break;
			
		case UNION:
		case INTERSECTION:
			c1.uniqueList_r(list);
			c2.uniqueList_r(list);
			break;
			
		default:
			Debug.Error("uniqueList_r: invalid operator.", true);
		}		
		
		return;
	}
	
	/**
	 * Run through a GSG expression looking at its leaves and return
	 * a list of the distinct leaves.  Note: leaf and leaf.complement() are
	 * not considered distinct.
	 * @param expression
	 * @return
	 */
	private ArrayList<CSG3D> uniqueList()
	{
		ArrayList<CSG3D> result = new ArrayList<CSG3D>();
		uniqueList_r(result);
		return result;
	}	
	
	
	/**
	 * Replace duplicate of leaf with leaf itself
	 * TODO: this should also use known complements
	 * @param leaf
	 * @param tolerance
	 */		
	private void replaceAllSameLeaves(CSG3D leaf, double tolerance)
	{	
		int same;
		switch(op)
		{
		case LEAF:
		case NULL:   
		case UNIVERSE:
			//System.out.println("replace_all_same_leaves(): at a leaf!");
			break;
			
		case UNION:
		case INTERSECTION:    
			HalfSpace hp = leaf.hp;
			if(c1.op == CSGOp.LEAF)
			{
				same = HalfSpace.same(hp, c1.hp, tolerance);
				if(same == 0)
					c1 = leaf;
				if(same == -1)
					c1 = leaf.complement();
			} else
				c1.replaceAllSameLeaves(leaf, tolerance);
			
			if(c2.op == CSGOp.LEAF)
			{
				same = HalfSpace.same(hp, c2.hp, tolerance);
				if(same == 0)
					c2 = leaf;
				if(same == -1)
					c2 = leaf.complement();
			} else
				c2.replaceAllSameLeaves(leaf, tolerance);
			break;
			
		default:
			Debug.Error("replace_all_same(): invalid operator.", true);		
		}
	}
	
	/**
	 * Replace duplicate of all leaves with the last instance of each; also
	 * link up complements.
	 * @param root
	 * @param tolerance
	 * @return simplified CSG object
	 */		
	private void simplify_r(CSG3D root, double tolerance)
	{
		switch(op)
		{
		case LEAF:
			root.replaceAllSameLeaves(this, tolerance);
			break;
			
		case NULL:   
		case UNIVERSE:
			//System.out.println("simplify_r(): at a leaf!");
			break;
			
		case UNION:
		case INTERSECTION:    
			c1.simplify_r(root, tolerance);
			c2.simplify_r(root, tolerance);
			break;
			
		default:
			Debug.Error("simplify_r(): invalid operator.", true);
		
		}
	}
	
	/**
	 * Replace duplicate of all leaves with the last instance of each
	 * @param tolerance
	 * @return simplified CSG object
	 */		
	public CSG3D simplify(double tolerance)
	{
		if(this == u || this == n)
			return this;
		
		CSG3D root = new CSG3D(this);
		simplify_r(root, tolerance);
		return root;
	}
	

	
	/**
	 * Offset by a distance (+ve or -ve)
	 * TODO: this should keep track of complements
	 * @param d
	 * @return offset CSG object by distance d
	 */
	public CSG3D offset(double d)
	{
		CSG3D result;
		
		switch(op)
		{
		case LEAF:
			result = new CSG3D(hp.offset(d));
			break;
			
		case NULL:
		case UNIVERSE:
			result = this;
			break;
			
		case UNION:
			result = union(c1.offset(d), c2.offset(d));
			break;
			
		case INTERSECTION:
			result = intersection(c1.offset(d), c2.offset(d));
			break;
			
		default:
			Debug.Error("offset(): invalid operator.", true);
			result = nothing();
		}
		return result;
	}
	
	
	/**
	 * leaf find the half-plane that generates the value for a point
	 * @param p
	 * @return leaf?
	 */
	public CSG3D leaf(Point3D p)
	{
		CSG3D result, r1, r2;
		
		switch(op)
		{
		case LEAF:
			result = this;
			break;
			
		case NULL:
			result = this;
			break;
			
		case UNIVERSE:
			result = this;
			break;
			
		case UNION:
			r1 = c1.leaf(p);
			r2 = c2.leaf(p);
			if(r1.value(p) < r2.value(p))
				return r1;
			else
				return r2;
			
		case INTERSECTION:
			r1 = c1.leaf(p);
			r2 = c2.leaf(p);
			if(r1.value(p) > r2.value(p))
				return r1;
			else
				return r2;
			
		default:
			Debug.Error("leaf(Rr2Point): invalid operator.", true);
			result = nothing();
		}
		return result;
	}
	
	/**
	 * "Potential" value of a point; i.e. a membership test
	 * -ve means inside; 0 means on the surface; +ve means outside
	 *
	 * @param p
	 * @return potential value of a point
	 */
	public double value(Point3D p)
	{
		double result = 1;
//		CSG3D c = leaf(p);
		switch(op)
		{
		case LEAF:
			result = hp.value(p);
			break;
			
		case NULL:
			result = 1;
			break;
			
		case UNIVERSE:
			result = -1;
			break;
			
		case UNION:
			result = Math.min(c1.value(p), c2.value(p));
			break;
			
		case INTERSECTION:
			result = Math.max(c1.value(p), c2.value(p));
			break;
			
		default:
			Debug.Error("CSG3D.value(): dud operator.", true);
		}
		return result;
	}
	

	
	/**
	 * The interval value of a box (analagous to point)
	 * @param b
	 * @return value of a box
	 */
	public Interval value(Box b)
	{
		Interval result;
		
		switch(op)
		{
		case LEAF:
			result = hp.value(b);
			break;
			
		case NULL:
			result = new Interval(1, 1.01);  // Is this clever?  Or dumb?
			break;
			
		case UNIVERSE:
			result = new Interval(-1.01, -1);  // Ditto.
			break;
			
		case UNION:
			result = Interval.min(c1.value(b), c2.value(b));
			break;
			
		case INTERSECTION:
			result = Interval.max(c1.value(b), c2.value(b));
			break;
			
		default:
			Debug.Error("value(RrBox): invalid operator.", true);
			result = new Interval();
		}
		
		return result;
	}
	
	/**
	 * Prune the set to a box
	 * @param b
	 * @return pruned box as new CSG object
	 */
	public CSG3D prune(Box b)
	{
		CSG3D result = this;
		
		switch(op)
		{
		case LEAF:            
			Interval i = hp.value(b);
			if (i.empty())
				Debug.Error("CSG3D.prune(RrBox): empty interval!", true);
			else if(i.neg())
				result = universe();
			else if (i.pos())
				result = nothing();
			break;
			
		case NULL:
		case UNIVERSE:
			break;
			
		case UNION:
			result =  union(c1.prune(b), c2.prune(b));
			break;
			
		case INTERSECTION:
			result = intersection(c1.prune(b), c2.prune(b));
			break;
			
		default:
			Debug.Error("CSG3D.prune(RrBox): dud op value!", true);
		}
		
		return result;
	}


}
