/*
 * pde.cpp
 *
 * Solves Laplace's equation for a potential field given source and sink voltages at points on the periphery.
 *
 *  Created on: 26 Jul 2019
 *
 *      Author: Adrian Bowyer
 *              RepRap Ltd
 *              https://reprapltd.com
 *
 *     Licence: GPL
 *
 */

#include <stdio.h>
#include <math.h>

// Set true for progress reports etc.

bool debug = true;

// Gauss-Seidel convergence criterion

const float convergence = 0.000001;

// Size of the grid

const int n = 50;
const int m = 50;

// Stop after this many iterations if there's no convergence

const int maxIterations = 2000;

// v[][] is the potential field, V. lastV is a copy of it from the last iteration.
// lastV is also used to store the magnitude of the field vectors, computed at the end.

float v[n+2][m+2], lastV[n+2][m+2];

//**********************************************************************************************

// Do a single pass of the Gauss-Seidel iteration.

float GaussSeidelOnePass()
{
	float rms = 0.0;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
		{
			v[i][j] = (v[i-1][j] + v[i+1][j] + v[i][j-1] + v[i][j+1])/4;
			float r = (v[i][j] - lastV[i][j]);
			rms = rms + r*r;
			lastV[i][j] = v[i][j];
		}
	}
	rms = sqrt(rms/( (float)m*(float)n) );
	return rms;
}


// Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
// of the differences between the last pass and the current one is less than the value
// of the variable convergence.

void GausSeidelIteration()
{
	float rms = 100.0*convergence;
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


// Compute the magnitudes of the gradient vectors in b.  This is
// the electric field, E.

void GradientMagnitudes()
{
	float xd, yd;
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
		{
			xd = 0.5*(v[i+1][j] - v[i-1][j]);
			yd = 0.5*(v[i][j+1] - v[i][j-1]);
			lastV[i][j] = sqrt(xd*xd + yd*yd);
		}
	}
}


// Set the boundary conditions and initialise

void BoundaryConditions()
{
	// Boundary conditions

	for(int i = 0;i <= n; i++)
	{
		v[i][1] = v[i][m] = 0.0;
	}
	for(int j = 0; j <= m; j++)
	{
		v[1][j] = v[n][j] = 0.0;
	}

	// One source and one sink at the edges

	v[0][m/2] = 10.0;
	v[n][m/2] = -10.0;

	// Initialise the rest to 0

	for(int i=1;i<n;i++)
		for(int j=1;j<m;j++)
		{
			v[i][j] = 0.0;
			lastV[i][j] = 0.0;
		}
}


// Outputs for GNUplot

void Output()
{
	FILE *fp;
	fp=fopen("potential.dat","w");
	for(int i = 0;i <= n; i++)
	{
		for(int j = 0; j <= m ; j++)
			fprintf(fp,"%f\n",v[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	// NB the field grid doesn't include the boundaries and so is 2 smaller in each direction.

	fp=fopen("field.dat","w");
	for(int i = 1; i < n; i++)
	{
		for(int j = 1; j < m; j++)
			fprintf(fp,"%f\n",lastV[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

}


// Self-explanatory, I hope.

int main()
{
	BoundaryConditions();
	GausSeidelIteration();
	GradientMagnitudes();
	Output();
}





