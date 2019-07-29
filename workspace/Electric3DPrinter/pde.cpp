/*
 * pde.cpp
 *
 * Solves Laplace's equation
 *
 *  Created on: 26 Jul 2019
 *      Author: Adrian Bowyer
 *     Licence: GPL
 *
 */
#include <stdio.h>
#include <math.h>
#define n 50
#define m 50
#define CONVERGENCE 0.000001

int main()
{
	int i, j, k;
	float a[n+2][m+2], b[n+2][m+2];

	// Boundary conditions

	for(i = 0;i <= n; i++)
	{
		a[i][1] = a[i][m] = 0.0;
	}
	for(j = 0; j <= m; j++)
	{
		a[1][j] = a[n][j] = 0.0;
	}

	// One source and one sink at the edges

	a[0][m/2] = 10.0;
	a[n][m/2] = -10.0;

	// Initialise the rest to 0

	for(i=1;i<n;i++)
		for(j=1;j<m;j++)
		{
			a[i][j] = 0.0;
			b[i][j] = 0.0;
		}

	// Gauss-Seidel iteration

	double r = 0.0;
	double rms = 100.0*CONVERGENCE;
	k = 0;
	while(k < 2000 && rms > CONVERGENCE)
	{
		rms = 0.0;
		for(i = 1; i < n; i++)
		{
			for(j = 1; j < m; j++)
			{
				a[i][j] = (a[i-1][j] + a[i+1][j] + a[i][j-1] + a[i][j+1])/4;
				r = (a[i][j] - b[i][j]);
				rms = rms + r*r;
				b[i][j] = a[i][j];
			}
		}
		rms = sqrt(rms/( (float)m*(float)n) );
		printf("Iteration: %i, rms: %f\n", k, rms);
		k++;
	}

	// Compute gradient magnitudes in b

	double xd, yd;
	for(i = 1; i < n; i++)
	{
		for(j = 1; j < m; j++)
		{
			xd = 0.5*(a[i+1][j] - a[i-1][j]);
			yd = 0.5*(a[i][j+1] - a[i][j-1]);
			b[i][j] = sqrt(xd*xd + yd*yd);
		}
	}

	// Outputs for GNUplot

	FILE *fp;
	fp=fopen("potential.dat","w");
	for(i = 0;i <= n; i++)
	{
		for(j = 0; j <= m ; j++)
			fprintf(fp,"%f\n",a[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

	// NB the field grid doesn't include the boundaries and so is 2 smaller in each direction.

	fp=fopen("field.dat","w");
	for(i = 1; i < n; i++)
	{
		for(j = 1; j < m; j++)
			fprintf(fp,"%f\n",b[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);

}





