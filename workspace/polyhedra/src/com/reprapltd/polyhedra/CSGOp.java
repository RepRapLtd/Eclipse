package com.reprapltd.polyhedra;

/**
 * Set operators, operands, and the universal and null sets
 * @author ensab
 *
 */
 
public enum CSGOp 
{
    LEAF("LEAF SET"), 
    NULL("NULL SET"), 
    UNIVERSE("UNIVERSAL SET"), 
    UNION("UNION"), 
    INTERSECTION("INTERSECTION"),
    DIFFERENCE("DIFFERENCE");  // ONLY ever used by CSGReader.  Should NEVER appear in a CSG tree
    
    /**
     * 
     */
    private String name;
    
    /**
     * @param name
     */
    CSGOp(String name)
    {
        this.name = name;
    }
    
    /* (non-Javadoc)
     * @see java.lang.Enum#toString()
     */
    public String toString() { return name; }
}

