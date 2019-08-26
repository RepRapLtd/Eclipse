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

#include <iostream>
#include <fstream>
#include <math.h>
#include <string>

using namespace std;

// Set true for progress reports etc.

bool debug = false;

// Gauss-Seidel convergence criterion

const double convergence = 0.0001;

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

// List of nodes on the 2D boundary in cyclic order. These are the same for all Z.

const int maxBoundary = 4*nodes;
int boundaryCount = 0;
int boundaryNodes[maxBoundary][2];

// Stop Gauss-Seidel after this many iterations if there's no convergence

const int maxIterations = 3000;

// potential[][][] is the potential field, V. lastPotential[][][] is a copy of it from the last iteration.
// field[][][] is used to store the magnitude of the field vectors, computed at the end of one solution.
// chargeIntegral[][][] is the accumulated charge that has flowed through each node for all solutions.
// thresholdedChargeIntegral[][][] is a thresholded version of chargeIntegral[][][] from a sigmoid function.

double potential[nodes+2][nodes+2][nodes+2], lastPotential[nodes+2][nodes+2][nodes+2], field[nodes+2][nodes+2][nodes+2],
       chargeIntegral[nodes+2][nodes+2][nodes+2], thresholdedChargeIntegral[nodes+2][nodes+2][nodes+2];

// True in the active region. This could be computed on the fly; but it's faster
// to store it.  What's memory for?

bool inside[nodes+2][nodes+2][nodes+2];

// Multiplication is faster than division...

const double oneSixth = 1.0/6.0;

//**********************************************************************************************

// Solve the PDE at a single node [i][j]

double PDE(const int x, const int y, const int z)
{
	// Do nothing outside the active region

	if(!inside[x][y][z])
		return potential[x][y][z];

	// Don't mess with the sources and sinks

	for(int l = 0; l < sources; l++)
	{
		if( (x == source[l][0]) && (y == source[l][1]) && (z == source[l][2]))
			return potential[x][y][z];
	}

	// Make a reflective boundary (i.e. one that does not conduct, so has 0 gradient)

	double vxm = potential[x-1][y][z];
	if(!inside[x-1][y][z])
		vxm = potential[x+1][y][z];

	double vxp = potential[x+1][y][z];
	if(!inside[x+1][y][z])
		vxp = potential[x-1][y][z];

	double vym = potential[x][y-1][z];
	if(!inside[x][y-1][z])
		vym = potential[x][y+1][z];

	double vyp = potential[x][y+1][z];
	if(!inside[x][y+1][z])
		vyp = potential[x][y-1][z];

	double vzm = potential[x][y][z-1];
	if(!inside[x][y][z-1])
		vzm = potential[x][y][z+1];

	double vzp = potential[x][y][z+1];
	if(!inside[x][y][z+1])
		vzp = potential[x][y][z-1];

	// The actual PDE

	return (vxm + vxp + vym + vyp + vzm + vzp)*oneSixth;
}

// Do a single pass of the Gauss-Seidel iteration.

double GaussSeidelOnePass()
{
	double rms = 0.0;
	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{
				potential[x][y][z] = PDE(x, y, z);
				double r = (potential[x][y][z] - lastPotential[x][y][z]);
				rms = rms + r*r;
				lastPotential[x][y][z] = potential[x][y][z];
			}
		}
	}
	rms = sqrt( rms/(double)(nodes*nodes*nodes) );
	return rms;
}


// Iterate the Gauss-Seidel until convergence.  Convergence is when the root-mean-square
// of the differences between the last pass and the current one is less than the value
// of the constant convergence.

void GausSeidelIteration()
{
	double rms = 100.0*convergence;
	int i = 0;
	while(i < maxIterations && rms > convergence)
	{
		rms = GaussSeidelOnePass();
		if(debug)
			cout << "Iteration: " << i << ", rms:" << rms << endl;
		i++;
	}
	if(i >= maxIterations)
		cout << "No convergence!, rms: " << rms << endl;
}


// Compute the magnitudes of the gradient vectors in field[][][].  This is
// the electric field, E.

void GradientMagnitudes()
{
	double xd, yd, zd;
	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{

				// At the edges use linear gradients; parabolas elsewhere

				if(!inside[x+1][y][z])
					xd = potential[x][y][z] - potential[x-1][y][z];
				else if(!inside[x-1][y][z])
					xd = potential[x+1][y][z] - potential[x][y][z];
				else
					xd = 0.5*(potential[x+1][y][z] - potential[x-1][y][z]);

				if(!inside[x][y+1][z])
					yd = potential[x][y][z] - potential[x][y-1][z];
				else if(!inside[x][y-1][z])
					yd = potential[x][y+1][z] - potential[x][y][z];
				else
					yd = 0.5*(potential[x][y+1][z] - potential[x][y-1][z]);

				if(!inside[x][y][z+1])
					zd = potential[x][y][z] - potential[x][y][z-1];
				else if(!inside[x][y][z-1])
					zd = potential[x][y][z+1] - potential[x][y][z];
				else
					zd = 0.5*(potential[x][y][z+1] - potential[x][y][z-1]);

				field[x][y][z] = sqrt(xd*xd + yd*yd + zd*zd);
				chargeIntegral[x][y][z] += field[x][y][z];
			}
		}
	}
}


// Find and print the maximum and minimum values of the charge integral.

void PrintChargeRange()
{
	// Assume the very mid point is not outside the active region...

	double low = chargeIntegral[xCentre][yCentre][nodes/2];
	double high = low;

	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{
				if(inside[x][y][z])
				{
					if(chargeIntegral[x][y][z] < low)
						low = chargeIntegral[x][y][z];
					if(chargeIntegral[x][y][z] > high)
						high = chargeIntegral[x][y][z];
				}
			}
		}
	}

	cout << "Charge range: " << low << " to " << high << endl;
}

// The sigmoid function that decides if a given charge integral will be solid

inline double Sigmoid(const double a, const double sigmoidOffset, const double sigmoidMultiplier)
{
	return 1.0 - exp(sigmoidMultiplier*(a - sigmoidOffset))/(exp(sigmoidMultiplier*(a - sigmoidOffset)) + 1.0);
}


// Threshold the charges.  Note - this INVERTS the charges high gives 0; low gives 1.

void SigmoidCharge(const double sigmoidOffset, const double sigmoidMultiplier)
{
	for(int x = 1; x < nodes; x++)
	{
		for(int y = 1; y < nodes; y++)
		{
			for(int z = 1; z < nodes; z++)
			{
				if(inside[x][y][z])
				{
					thresholdedChargeIntegral[x][y][z] = Sigmoid(chargeIntegral[x][y][z], sigmoidOffset, sigmoidMultiplier);
				} else
					thresholdedChargeIntegral[x][y][z] = 0.0;
			}
		}
	}
}


// Initialise the charges at the nodes.  Set them all 0 inside and outside the active region.

void ChargeSetUp()
{
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes; y++)
		{
			for(int z = 0; z <= nodes; z++)
			{
				chargeIntegral[x][y][z] = 0.0;
			}
		}
	}
}

// Is the node [xx][yy][zz] on the boundary of the active region?
// Tiny inefficiency - it checks [xx][yy][zz] twice. 1/27th unnecessary...

bool OnBoundary(const int x, const int y, const int z)
{
	if(!inside[x][y][z])
		return false;

	for(int xDel = -1; xDel < 2; xDel++)
	{
		for(int yDel = -1; yDel < 2; yDel++)
		{
			for(int zDel = -1; zDel < 2; zDel++)
			{
				if(!inside[x+xDel][y+yDel][z+zDel])
					return true;
			}
		}
	}

	return false;
}

// Find the 2D boundary pattern.  This is a disc in the middle of the Z range,
// and is the same for all Z values.  This uses the 3D OnBoundary() function above, rather
// than bother with an extra 2D version.

void FindBoundary()
{
	// Find the first quadrant

	boundaryCount = 0;
	int z = nodes/2;

	for(int x = 0; x <= xCentre; x++)
	{
		for(int y = yCentre; y >= 0; y--)
		{
			if(OnBoundary(x, y, z))
			{
				if(boundaryCount >= maxBoundary)
				{
					cout << "Maximum boundary node count (Q0) exceeded!" << endl;
					return;
				}

				boundaryNodes[boundaryCount][0] = x;
				boundaryNodes[boundaryCount][1] = y;
				boundaryCount++;
			}
		}
	}

	// Copy that to the other three

	int bc = boundaryCount;
	for(int x = bc - 2; x >= 0; x--)
	{
		if(boundaryCount >= maxBoundary)
		{
			cout << "Maximum boundary node count (Q1) exceeded!" << endl;
			return;
		}
		boundaryNodes[boundaryCount][0] = 2*xCentre - boundaryNodes[x][0];
		boundaryNodes[boundaryCount][1] = boundaryNodes[x][1];
		boundaryCount++;
	}

	bc = boundaryCount;
	for(int x = bc - 2; x > 0; x--)
	{
		if(boundaryCount >= maxBoundary)
		{
			cout << "Maximum boundary node count (Q34) exceeded!" << endl;
			return;
		}
		boundaryNodes[boundaryCount][0] = boundaryNodes[x][0];
		boundaryNodes[boundaryCount][1] = 2*yCentre - boundaryNodes[x][1];
		boundaryCount++;
	}

	if(boundaryCount%4)
		cout << "Number of boundary nodes is not a multiple of 4! " << boundaryCount << endl;
}


// Initialise all the tensors for potential, field etc to 0, and also
// set up the active region.

void Initialise()
{
	// Set up the active area and initialise the solution to 0.

	for(int x = 0; x <= nodes; x++)
	{
		int xd = x - xCentre;
		for(int y = 0; y <= nodes; y++)
		{
			int yd = y - yCentre;
			bool in = xd*xd + yd*yd < radius*radius;
			for(int z = 0; z <= nodes; z++)
			{
				inside[x][y][z] = in;
				potential[x][y][z] = 0.0;
				lastPotential[x][y][z] = 0.0;
			}

			// Bottom and top

			inside[x][y][0] = false;
			inside[x][y][nodes] = false;
		}
	}
}


// Set the boundary conditions and initialise one solution with potentials applied at at z=k.
// b is the index into the boundary array that decides how far round the circle voltages will
// be applied.  v is the potential.

void BoundaryConditions(const int b, const int z, const double v)
{

	Initialise();
	FindBoundary();

	// Sources and sinks

//	source[0][0] = xc + round((double)(radius - 1)*cos(angle));
//	source[0][1] = yc + round((double)(radius - 1)*sin(angle));
//	source[1][0] = xc + round((double)(radius - 1)*cos(angle + M_PI));
//	source[1][1] = yCentre + round((double)(radius - 1)*sin(angle + M_PI));

	source[0][0] = boundaryNodes[b][0];
	source[0][1] = boundaryNodes[b][1];
	//double angle = atan2(yCentre - source[0][1], xCentre - source[0][0]);
	int opposite = (b + boundaryCount/2)%boundaryCount;
	source[1][0] = boundaryNodes[opposite][0];
	source[1][1] = boundaryNodes[opposite][1];

//	potential[source[0][0]][source[0][1]] = 2.0 + sin(4.0*angle);
//	potential[source[1][0]][source[1][1]] = 2.0 + sin(4.0*angle + M_PI);

	potential[source[0][0]][source[0][1]][z] = v;
	potential[source[1][0]][source[1][1]][z] = -v;
}


// Output one disc at z=k for gnuplot.  If activeR is positive, just output that
// radius of the disc for close-ups of the middle.

void Output(const char* name, const double a[nodes+2][nodes+2][nodes+2], const int activeR, const int z)
{
	// Find the most negative value in the mesh and use that
	// as the values outside the disc.

	double negValue = a[xCentre][yCentre][z];
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes ; y++)
		{
			if(activeR <= 0)
			{
				if(inside[x][y][z] && a[x][y][z] < negValue)
					negValue = a[x][y][z];
			} else
			{
				int xd = x - xCentre;
				int yd = y - yCentre;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					if(a[x][y][z] < negValue)
						negValue = a[x][y][z];
				}
			}
		}
	}

	// Stick the data in a file so that GNUPlot can
	// plot it.

	ofstream outputFile;
	outputFile.open(name);
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes ; y++)
		{
			double val = negValue;
			if(activeR <= 0)
			{
				if(inside[x][y][z])
					val = a[x][y][z];
			} else
			{
				int xd = x - xCentre;
				int yd = y - yCentre;
				if(xd*xd + yd*yd < activeR*activeR)
				{
					val = a[x][y][z];
				}
			}
			outputFile << val << '\n';
		}
		outputFile << '\n';
	}
	outputFile.close();
}


// Output the tensor that is the entire mesh so that a 3D iso-surface STL file
// can be generated from it.

void OutputTensor(const char* name, const double a[nodes+2][nodes+2][nodes+2])
{
	// Assume the central point is active...
	double minValue = a[xCentre][yCentre][nodes/2];
	double maxValue = minValue;

	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes; y++)
		{
			for(int z = 0; z <= nodes; z++)
			{
				if(inside[x][y][z])
				{
					if(a[x][y][z] < minValue)
						minValue = a[x][y][z];
					if(a[x][y][z] > maxValue)
						maxValue = a[x][y][z];
				}
			}
		}
	}

	//if(debug)
		cout << "Tensor minimum and maximum: " << minValue << ", " << maxValue << endl;

	// Shift the minimum down a bit for everything that's outside

	//minValue = minValue - 0.1*(maxValue - minValue);
//	if(debug)
//		cout << "New minimum: " << minValue << endl;

	ofstream outputFile;
	outputFile.open(name);

	// Need an extra blank layer so the marching cubes can see the top.

	outputFile << nodes+1 << ' ' << nodes+1 << ' ' << nodes+2 << ' ' << minValue << ' ' << maxValue;

	for(int z = 0; z <= nodes; z++)
	{
		for(int y = 0; y <= nodes; y++)
		{
			for(int x = 0; x <= nodes; x++)
			{
				outputFile << ' ';
				if(!inside[x][y][z])
				{
					outputFile << minValue;
				} else
				{
					outputFile << a[x][y][z];
				}
			}
		}
	}

	// Blank layer

	for(int y = 0; y <= nodes; y++)
	{
		for(int x = 0; x <= nodes; x++)
		{
			outputFile << ' ';
			outputFile << minValue;
		}
	}
	outputFile.close();
}


// Function to plot the boundary to test that it's properly set-up.
// Not normally called.

void TestBoundary()
{
	int z = nodes/2;
	BoundaryConditions(0.0, z, 1.0);
	for(int x = 0; x <= nodes; x++)
	{
		for(int y = 0; y <= nodes ; y++)
			inside[x][y][z] = true;
	}
	potential[1][1][z] = 0.5;
	for(int x = 0; x < boundaryCount; x++)
	{
		potential[boundaryNodes[x][0]][boundaryNodes[x][1]][z] += 1;
		cout << x << ": (" << boundaryNodes[x][0] << ", " << boundaryNodes[x][1] << ")" << endl;
	}
	Output("boundary.dat", potential, -1, z);
}

// Create a cylinder pattern in thresholdedChargeIntegral[][][] and write it out as
// a tensor to test tensor output.  Not normally called.

void TestCylinder(const int r, const int z0, const int z1)
{
	Initialise();
	for(int x = 0; x <= nodes; x++)
	{
		int xd = x - xCentre;
		for(int y = 0; y <= nodes; y++)
		{
			int yd = y - yCentre;
			for(int z = 0; z <= nodes; z++)
			{
				if(z < z0 || z > z1 || (xd*xd + yd*yd) > r*r)
				{
					thresholdedChargeIntegral[x][y][z] = 0.0;
				} else
				{
					thresholdedChargeIntegral[x][y][z] = 1.0;
				}
			}
		}
	}
	OutputTensor("testCylinder.tns", thresholdedChargeIntegral);
}

void Prompt()
{
	cout << endl << "Commands:" << endl;
	cout << " s: - set sigmoid threshold and output a slice" << endl;
	cout << " t: - create the tensor file" << endl;
	cout << " q: quit" << endl<< endl;
}

// Decide how to process the results.

void Control()
{
	cout << "Type h for help." << endl;
	while(1)
	{
		cout << "Command: ";
		char c;
		cin >> c;

		switch(c)
		{
		case 't':
			OutputTensor("thresholdTensor.tns", thresholdedChargeIntegral);
			break;

		case 's':
			double s;
			cout << "Sigmoid value for 0.5 point: ";
			cin >> s;
			SigmoidCharge(s, 50);
			Output("threshold.dat", thresholdedChargeIntegral, -1, nodes/2);
			break;

		case 'q':
			return;

		default:
			cout << endl << "Unrecognised command - " << c << endl;
		case 'h':
			Prompt();
		}
	}
}


// Self-explanatory, I hope.

int main()
{
	//	TestBoundary();
	//  TestCylinder(15, 10, 40);


	BoundaryConditions(0, nodes/2, 1.0);


	//for(int vv = 1; vv < 21; vv++)
	//{
	int vv = 3;
		ChargeSetUp();
		cout << "Z for " << vv << " volts: ";
		for(int z = 1; z < nodes; z++)
		{
			cout << z << ' ';
			cout.flush();

			double v = (double)vv;

			for(int angle = 0; angle < boundaryCount/2; angle++)
			{
				BoundaryConditions(angle, z, v);
				GausSeidelIteration();
				GradientMagnitudes();
			}
		}
		cout << endl;
		SigmoidCharge(0.0, 50);
		string fileName = "thresholdTensor";
		char str[20];
		sprintf(str,"-%d.tns",vv);
		fileName += str;
		OutputTensor(fileName.c_str(), thresholdedChargeIntegral);
	//}

//	Output("potential.dat", potential, -1, nodes/2);
//	Output("field.dat", field, -1, nodes/2);
//	Output("charge.dat", chargeIntegral, -1, nodes/2);
//
//	PrintChargeRange();
//
//	Control();
}






