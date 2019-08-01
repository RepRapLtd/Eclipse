/*
 * pde.cpp
 *
 * Solves Laplace's equation for a potential field given source and sink voltages at points on the periphery.
 *
 * The region of solution is a disk, and the sources and sinks are specified as a list of potentials.
 *
 * The sources and sinks can be moved between several solutions and the cumulated charge moved through each pixel
 * calculated.
 *
 *  Created on: 26 Jul 2019
 *
 *      Author: Adrian Bowyer
 *              RepRap Ltd
 *              https://reprapltd.com
 *
 *     Licence: GPL
 *
 *
 * To plot the results:
 *
 *  $ gnuplot
 *  gnuplot> set hidden3d
 *  gnuplot> splot 'potential.dat' with lines
 *
 * The output files are:
 *
 *   potential.dat - the electric potentials
 *   field.dat - the magnitude of the electric field
 *   charge.dat - the cumulated charge at each node
 *
 */

#include <stdio.h>
#include <math.h>

// Set true for progress reports etc.

bool debug = false;

// Gauss-Seidel convergence criterion

const double convergence = 0.000001;

// Size of the grid

const int n = 50;
const int m = 50;

// The centre and radius of the disc

const int xc = 25;
const int yc = 25;
const int radius = 22;

// The source and sink

const int sources = 2;

int fixed[sources][2];

// Stop after this many iterations if there's no convergence

const int maxIterations = 3000;

// v[][] is the potential field, V. lastV is a copy of it from the last iteration.
// e[][] is used to store the magnitude of the field vectors, computed at the end of one solution.
// c[][] is the accumulated charge that has flowed through each node for all solutions.

double v[n+2][m+2], lastV[n+2][m+2], e[n+2][m+2], c[n+2][m+2];

// Run the simulation this many times, incrementing the angle of the electrodes each time.

const int angles = 20;

// True in the active region. This could be computed on the fly; but it's faster
// to store it.  What's memory for?

bool inside[n+2][m+2];

//**********************************************************************************************

// Solve the PDE at a single node [i][j]

double PDE(int i, int j)
{
	// Do nothing outside the disc

	if(!inside[i][j])
		return v[i][j];

	// Don't mess with the sources and sinks

	for(int k = 0; k < sources; k++)
	{
		if( (i == fixed[k][0]) && (j == fixed[k][1]) )
			return v[i][j];
	}

	// Make a reflective boundary (i.e. one that does not conduct, so has 0 gradient)

	double vxm = v[i-1][j];
	if(!inside[i-1][j])
		vxm = v[i+1][j];

	double vxp = v[i+1][j];
	if(!inside[i+1][j])
		vxp = v[i-1][j];

	double vym = v[i][j-1];
	if(!inside[i][j-1])
		vym = v[i][j+1];

	double vyp = v[i][j+1];
	if(!inside[i][j+1])
		vyp = v[i][j-1];

	// The actual PDE

	return 0.25*(vxm + vxp + vym + vyp);
}

// Do a single pass of the Gauss-Seidel iteration.

double GaussSeidelOnePass()
{
	double rms = 0.0;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
		{
			v[i][j] = PDE(i, j);
			double r = (v[i][j] - lastV[i][j]);
			rms = rms + r*r;
			lastV[i][j] = v[i][j];
		}
	}
	rms = sqrt(rms/( (double)m*(double)n) );
	return rms;
}


// Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
// of the differences between the last pass and the current one is less than the value
// of the constant convergence.

void GausSeidelIteration()
{
	double rms = 100.0*convergence;
	int k = 0;
	while(k < maxIterations && rms > convergence)
	{
		rms = GaussSeidelOnePass();
		if(debug)
			printf("Iteration: %i, rms: %f\n", k, rms);
		k++;
	}
	if(k >= maxIterations)
		printf("No convergence!, rms: %f\n", rms);
}


// Compute the magnitudes of the gradient vectors in e[][].  This is
// the electric field, E.

void GradientMagnitudes()
{
	double xd, yd;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
		{

			// At the edges use linear gradients; parabolas elsewhere

			if(!inside[i+1][j])
				xd = v[i][j] - v[i-1][j];
			else if(!inside[i-1][j])
				xd = v[i+1][j] - v[i][j];
			else
				xd = 0.5*(v[i+1][j] - v[i-1][j]);

			if(!inside[i][j+1])
				yd = v[i][j] - v[i][j-1];
			else if(!inside[i][j-1])
				yd = v[i][j+1] - v[i][j];
			else
				yd = 0.5*(v[i][j+1] - v[i][j-1]);

			e[i][j] = sqrt(xd*xd + yd*yd);
			c[i][j] += e[i][j];
		}
	}
}


// Initialise the charges at the nodes

void ChargeSetUp()
{
	for(int i=1;i<n;i++)
	{
		for(int j=1;j<m;j++)
		{
			c[i][j] = 0.0;
		}
	}
}


// Set the boundary conditions and initialise one solution

void BoundaryConditions(double angle)
{
	// Set up the active area and initialise the solution to 0.

	for(int i=1;i<n;i++)
	{
		int xd = i - xc;
		for(int j=1;j<m;j++)
		{
			int yd = j - yc;
			inside[i][j] = xd*xd + yd*yd < radius*radius;
			v[i][j] = 0.0;
			lastV[i][j] = 0.0;
		}
	}

	// Source and sink

	fixed[0][0] = xc + round((double)(radius - 1)*cos(angle));
	fixed[0][1] = yc + round((double)(radius - 1)*sin(angle));
	fixed[1][0] = xc + round((double)(radius - 1)*cos(angle + M_PI));
	fixed[1][1] = yc + round((double)(radius - 1)*sin(angle + M_PI));

	v[fixed[0][0]][fixed[0][1]] = 1.0;
	v[fixed[1][0]][fixed[1][1]] = -1.0;
}


// Output for GNUplot.  If activeR is positive, just output that
// radius of the disc for close-ups of the middle.

void Output(char* name, double a[n+2][m+2], int activeR)
{
	// Find the most negative value in the mesh and use that
	// as the values outside the disc.

	double negValue = a[xc][yc];
	for(int i = 0; i <= n; i++)
		for(int j = 0; j <= m ; j++)
		{
			if(activeR <= 0)
			{
				if(a[i][j] < negValue)
					negValue = a[i][j];
			} else
			{
				int xd = i - xc;
				int yd = j - yc;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					if(a[i][j] < negValue)
						negValue = a[i][j];
				}
			}
		}

	// Stick the data in a file so that GNUPlot can
	// plot it.

	FILE *fp;
	fp=fopen(name,"w");
	for(int i = 0; i <= n; i++)
	{
		for(int j = 0; j <= m ; j++)
		{
			double val = negValue;
			if(activeR <= 0)
			{
				if(inside[i][j])
					val = a[i][j];
			} else
			{
				int xd = i - xc;
				int yd = j - yc;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					val = a[i][j];
				}
			}
			fprintf(fp,"%f\n", val);
		}
		fprintf(fp,"\n");
	}
	fclose(fp);
}


// Self-explanatory, I hope.

int main()
{
	ChargeSetUp();

	double angle = 0.0;
	double aInc = 2.0*M_PI/((double)angles);
	for(int a = 0; a < angles; a++)
	{
		printf("Angle: %f\n", angle);
		BoundaryConditions(angle);
		GausSeidelIteration();
		GradientMagnitudes();
		angle += aInc;
	}

	Output("potential.dat", v, -1);
	Output("field.dat", e, -1);
	Output("charge.dat", c, -1);
}





