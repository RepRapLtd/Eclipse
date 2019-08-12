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

//#include <stdio.h>
#include <iostream>
#include <fstream>
#include <math.h>

using namespace std;

// Set true for progress reports etc.

bool debug = false;

// Gauss-Seidel convergence criterion

const double convergence = 0.00001;

// Size of the grid

const int nodes = 50;

// The centre and radius of the disc

const int xCentre = 25;
const int yCentre = 25;
const int radius = 22;

// The source and sink
// NB sources must be 2*N where N is odd.

const int sources = 2;

int source[sources][3];

// List of nodes on the boundary in cyclic order

const int maxBoundary = 4*nodes;
int boundaryCount = 0;
int boundaryNodes[maxBoundary][2];

// Stop after this many iterations if there's no convergence

const int maxIterations = 3000;

// potential[][] is the potential field, V. lastPotential[][] is a copy of it from the last iteration.
// field[][] is used to store the magnitude of the field vectors, computed at the end of one solution.
// chargeIntegral[][] is the accumulated charge that has flowed through each node for all solutions.
// thresholdedChargeIntegral[][] is a thresholded version of c[][] subject to, for example, a sigmoid function.

double potential[nodes+2][nodes+2][nodes+2], lastPotential[nodes+2][nodes+2][nodes+2], field[nodes+2][nodes+2][nodes+2],
       chargeIntegral[nodes+2][nodes+2][nodes+2], thresholdedChargeIntegral[nodes+2][nodes+2][nodes+2];

// True in the active region. This could be computed on the fly; but it's faster
// to store it.  What's memory for?

bool inside[nodes+2][nodes+2][nodes+2];

//**********************************************************************************************

// Solve the PDE at a single node [i][j]

double PDE(int i, int j, int k)
{
	// Do nothing outside the disc

	if(!inside[i][j][k])
		return potential[i][j][k];

	// Don't mess with the sources and sinks

	for(int l = 0; l < sources; l++)
	{
		if( (i == source[l][0]) && (j == source[l][1]) && (k == source[l][2]))
			return potential[i][j][k];
	}

	// Make a reflective boundary (i.e. one that does not conduct, so has 0 gradient)

	double vxm = potential[i-1][j][k];
	if(!inside[i-1][j][k])
		vxm = potential[i+1][j][k];

	double vxp = potential[i+1][j][k];
	if(!inside[i+1][j][k])
		vxp = potential[i-1][j][k];

	double vym = potential[i][j-1][k];
	if(!inside[i][j-1][k])
		vym = potential[i][j+1][k];

	double vyp = potential[i][j+1][k];
	if(!inside[i][j+1][k])
		vyp = potential[i][j-1][k];

	double vzm = potential[i][j][k-1];
	if(!inside[i][j][k-1])
		vzm = potential[i][j][k+1];

	double vzp = potential[i][j][k+1];
	if(!inside[i][j][k+1])
		vzp = potential[i][j][k-1];

	// The actual PDE

	return (vxm + vxp + vym + vyp + vzm + vzp)/6.0;
}

// Do a single pass of the Gauss-Seidel iteration.

double GaussSeidelOnePass()
{
	double rms = 0.0;
	for(int i = 1; i < nodes; i++)
	{
		for(int j = 1; j < nodes; j++)
		{
			for(int k = 1; k < nodes; k++)
			{
				potential[i][j][k] = PDE(i, j, k);
				double r = (potential[i][j][k] - lastPotential[i][j][k]);
				rms = rms + r*r;
				lastPotential[i][j][k] = potential[i][j][k];
			}
		}
	}
	rms = sqrt( rms/(double)(nodes*nodes) );
	return rms;
}


// Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
// of the differences between the last pass and the current one is less than the value
// of the constant convergence.

void GausSeidelIteration()
{
	double rms = 100.0*convergence;
	int l = 0;
	while(l < maxIterations && rms > convergence)
	{
		rms = GaussSeidelOnePass();
		if(debug)
			cout << "Iteration: " << l << ", rms:" << rms << endl;
		l++;
	}
	if(l >= maxIterations)
		cout << "No convergence!, rms: " << rms << endl;
}


// Compute the magnitudes of the gradient vectors in e[][].  This is
// the electric field, E.

void GradientMagnitudes()
{
	double xd, yd, zd;
	for(int i = 1; i < nodes; i++)
	{
		for(int j = 1; j < nodes; j++)
		{
			for(int k = 1; k < nodes; k++)
			{

				// At the edges use linear gradients; parabolas elsewhere

				if(!inside[i+1][j][k])
					xd = potential[i][j][k] - potential[i-1][j][k];
				else if(!inside[i-1][j][k])
					xd = potential[i+1][j][k] - potential[i][j][k];
				else
					xd = 0.5*(potential[i+1][j][k] - potential[i-1][j][k]);

				if(!inside[i][j+1][k])
					yd = potential[i][j][k] - potential[i][j-1][k];
				else if(!inside[i][j-1][k])
					yd = potential[i][j+1][k] - potential[i][j][k];
				else
					yd = 0.5*(potential[i][j+1][k] - potential[i][j-1][k]);

				if(!inside[i][j][k+1])
					zd = potential[i][j][k] - potential[i][j][k-1];
				else if(!inside[i][j][k-1])
					zd = potential[i][j][k+1] - potential[i][j][k];
				else
					zd = 0.5*(potential[i][j][k+1] - potential[i][j][k-1]);

				field[i][j][k] = sqrt(xd*xd + yd*yd + zd*zd);
				chargeIntegral[i][j][k] += field[i][j][k];
			}
		}
	}
}

void PrintChargeRange()
{
	double low = chargeIntegral[xCentre][yCentre][nodes/2];
	double high = low;
	for(int i=1;i<nodes;i++)
	{
		for(int j=1;j<nodes;j++)
		{
			for(int k=1;k<nodes;k++)
			{
				if(inside[i][j][k])
				{
					if(chargeIntegral[i][j][k] < low)
						low = chargeIntegral[i][j][k];
					if(chargeIntegral[i][j][k] > high)
						high = chargeIntegral[i][j][k];
				}
			}
		}
	}

	cout << "Charge range: " << low << " to " << high << endl;
}

// The sigmoid function that decides if a given charge integral will be solid

double Sigmoid(double a, double sigmoidOffset, double sMultiplier)
{
	return 1.0 - exp(sMultiplier*(a - sigmoidOffset))/(exp(sMultiplier*(a - sigmoidOffset)) + 1.0);
}


// Threshold the charges.  Note - this INVERTS the charges high gives 0; low gives 1.

void SigmoidCharge(double sigmoidOffset, double sMultiplier)
{
	for(int i=1;i<nodes;i++)
	{
		for(int j=1;j<nodes;j++)
		{
			for(int k=1;k<nodes;k++)
			{
				if(inside[i][j][k])
				{
					thresholdedChargeIntegral[i][j][k] = Sigmoid(chargeIntegral[i][j][k], sigmoidOffset, sMultiplier);
				} else
					thresholdedChargeIntegral[i][j][k] = 0.0;
			}
		}
	}
}


// Initialise the charges at the nodes

void ChargeSetUp()
{
	for(int i=0;i<=nodes;i++)
	{
		for(int j=0;j<=nodes;j++)
		{
			for(int k=0;k<=nodes;k++)
			{
				chargeIntegral[i][j][k] = 0.0;
			}
		}
	}
}

bool OnBoundary(int i, int j, int k)
{
	if(!inside[i][j][k])
		return false;

	for(int ii = -1; ii < 2; ii++)
	{
		for(int jj = -1; jj < 2; jj++)
		{
			for(int kk = -1; kk < 2; kk++)
			{
				if(!inside[i+ii][j+jj][k+kk])
					return true;
			}
		}
	}

	return false;
}

void FindBoundary(int k)
{
	// Find the first quadrant

	boundaryCount = 0;
	for(int i = 0; i <= xCentre; i++)
	{
		for(int j = yCentre; j >= 0; j--)
		{
			if(OnBoundary(i,j, k))
			{
				if(boundaryCount >= maxBoundary)
				{
					cout << "Maximum boundary node count (Q0) exceeded!" << endl;
					return;
				}

				boundaryNodes[boundaryCount][0] = i;
				boundaryNodes[boundaryCount][1] = j;
				boundaryCount++;
			}
		}
	}

	// Copy that to the other three

	int bc = boundaryCount;
	for(int i = bc - 2; i >= 0; i--)
	{
		if(boundaryCount >= maxBoundary)
		{
			cout << "Maximum boundary node count (Q1) exceeded!" << endl;
			return;
		}
		boundaryNodes[boundaryCount][0] = 2*xCentre - boundaryNodes[i][0];
		boundaryNodes[boundaryCount][1] = boundaryNodes[i][1];
		boundaryCount++;
	}

	bc = boundaryCount;
	for(int i = bc - 2; i > 0; i--)
	{
		if(boundaryCount >= maxBoundary)
		{
			cout << "Maximum boundary node count (Q34) exceeded!" << endl;
			return;
		}
		boundaryNodes[boundaryCount][0] = boundaryNodes[i][0];
		boundaryNodes[boundaryCount][1] = 2*yCentre - boundaryNodes[i][1];
		boundaryCount++;
	}

	if(boundaryCount%4)
		cout << "Number of boundary nodes is not a multiple of 4! " << boundaryCount << endl;
}


// Set the boundary conditions and initialise one solution

void BoundaryConditions(int b, int k)
{
	// Set up the active area and initialise the solution to 0.

	for(int i = 0; i <= nodes; i++)
	{
		int xd = i - xCentre;
		for(int j = 0; j <= nodes; j++)
		{
			int yd = j - yCentre;
			for(int kk = 1; kk < nodes; kk++)
			{
				inside[i][j][kk] = xd*xd + yd*yd < radius*radius;
				potential[i][j][kk] = 0.0;
				lastPotential[i][j][kk] = 0.0;
			}

			// Bottom and top

			inside[i][j][0] = false;
			potential[i][j][0] = 0.0;
			lastPotential[i][j][0] = 0.0;
			inside[i][j][nodes] = false;
			potential[i][j][nodes] = 0.0;
			lastPotential[i][j][nodes] = 0.0;
		}
	}

	FindBoundary(nodes/2);

	// Sources and sinks

//	source[0][0] = xc + round((double)(radius - 1)*cos(angle));
//	source[0][1] = yc + round((double)(radius - 1)*sin(angle));
//	source[1][0] = xc + round((double)(radius - 1)*cos(angle + M_PI));
//	source[1][1] = yCentre + round((double)(radius - 1)*sin(angle + M_PI));

	source[0][0] = boundaryNodes[b][0];
	source[0][1] = boundaryNodes[b][1];
	double angle = atan2(yCentre - source[0][1], xCentre - source[0][0]);
	int opposite = (b + boundaryCount/2)%boundaryCount;
	source[1][0] = boundaryNodes[opposite][0];
	source[1][1] = boundaryNodes[opposite][1];

//	potential[source[0][0]][source[0][1]] = 2.0 + sin(4.0*angle);
//	potential[source[1][0]][source[1][1]] = 2.0 + sin(4.0*angle + M_PI);

	potential[source[0][0]][source[0][1]][k] = 1.0;
	potential[source[1][0]][source[1][1]][k] = -1.0;
}


// Output one disc at z=k for GNUplot.  If activeR is positive, just output that
// radius of the disc for close-ups of the middle.

void Output(char* name, double a[nodes+2][nodes+2][nodes+2], int activeR, int k)
{
	// Find the most negative value in the mesh and use that
	// as the values outside the disc.

	double negValue = a[xCentre][yCentre][k];
	for(int i = 0; i <= nodes; i++)
		for(int j = 0; j <= nodes ; j++)
		{
			if(activeR <= 0)
			{
				if(inside[i][j][k] && a[i][j][k] < negValue)
					negValue = a[i][j][k];
			} else
			{
				int xd = i - xCentre;
				int yd = j - yCentre;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					if(a[i][j][k] < negValue)
						negValue = a[i][j][k];
				}
			}
		}

	// Stick the data in a file so that GNUPlot can
	// plot it.

	ofstream outputFile;
	outputFile.open(name);
	for(int i = 0; i <= nodes; i++)
	{
		for(int j = 0; j <= nodes ; j++)
		{
			double val = negValue;
			if(activeR <= 0)
			{
				if(inside[i][j][k])
					val = a[i][j][k];
			} else
			{
				int xd = i - xCentre;
				int yd = j - yCentre;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					val = a[i][j][k];
				}
			}
			outputFile << val << '\n';
		}
		outputFile << '\n';
	}
	outputFile.close();
}

void OutputTensor(char* name, double a[nodes+2][nodes+2][nodes+2])
{
	double minValue = a[0][0][0];
	double maxValue = minValue;
	for(int i = 0; i <= nodes; i++)
	{
		for(int j = 0; j <= nodes; j++)
		{
			for(int k = 1; k < nodes; k++)
			{
				if(a[i][j][k] < minValue)
					minValue = a[i][j][k];
				if(a[i][j][k] > maxValue)
					maxValue = a[i][j][k];
			}
		}
	}

	cout << "Tensor minimum and maximum: " << minValue << ", " << maxValue << endl;

	ofstream outputFile;
	outputFile.open(name);
	outputFile << nodes << ' ' << nodes << ' ' << nodes << ' ' << minValue << ' ' << maxValue;

	for(int k = 0; k <= nodes; k++)
	{
		for(int j = 0; j <= nodes; j++)
		{
			for(int i = 1; i < nodes; i++)
			{
				outputFile << ' ' << a[i][j][k];
			}
		}
	}
	outputFile.close();
}


// Function to plot the boundary to test that it's properly set-up.
// Not normally called.

void TestBoundary()
{
	int k = nodes/2;
	BoundaryConditions(0.0, k);
	for(int i = 0; i <= nodes; i++)
	{
		for(int j = 0; j <= nodes ; j++)
			inside[i][j][k] = true;
	}
	potential[1][1][k] = 0.5;
	for(int i = 0; i < boundaryCount; i++)
	{
		potential[boundaryNodes[i][0]][boundaryNodes[i][1]][k] += 1;
		cout << i << ": (" << boundaryNodes[i][0] << ", " << boundaryNodes[i][1] << ")" << endl;
	}
	Output("boundary.dat", potential, -1, k);
}


// Self-explanatory, I hope.

int main()
{
//	TestBoundary();

	ChargeSetUp();
	int k = nodes/2;
	BoundaryConditions(0, k);

	for(int i = 0; i < boundaryCount/4; i++)
	{
		BoundaryConditions(0, k);
		GausSeidelIteration();
		GradientMagnitudes();
	}
	for(int i = boundaryCount/2; i < (3*boundaryCount)/4; i++)
	{
		BoundaryConditions(i, k);
		GausSeidelIteration();
		GradientMagnitudes();
	}

    //OutputTensor("potentialTensor.dat", potential);

	Output("potential.dat", potential, -1, k);
	Output("field.dat", field, -1, k);
	Output("charge.dat", chargeIntegral, -1, k);

	PrintChargeRange();

	double s;
	do
	{
		cout << "Sigmoid value for 0.5 point (-ve to exit): ";
		cin >> s;
		if(s > 0.0)
		{
			SigmoidCharge(s, 50);
		    OutputTensor("thresholdTensor.dat", thresholdedChargeIntegral);
			Output("threshold.dat", thresholdedChargeIntegral, -1, k);
		}
	} while(s > 0.0);

}





