package com.reprapltd.polyhedra;

public class Debug 
{
	private static boolean debugging = true;

	public Debug() 
	{
	}
	
	public static void Error(String s, boolean serious)
	{
		if(!serious && !debugging)
			return;
		System.out.println(s);
	}
	
	public static void Warning(String s)
	{
		if(!debugging)
			return;
		System.out.println(s);
	}

}
