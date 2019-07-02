/*
 * segment.cxx
 *
 * A C++/Imagemagick program to divide one image into many
 *
 * Adrian Bowyer
 * 25 May 2003
 *
 * A.Bowyer@bath.ac.uk
 *
 */

#include <Magick++.h>
#include <iostream.h>
#include <strings.h>
#include <strstream.h>
#include <math.h>
 
using namespace std;
using namespace Magick;

// String array length (slight hack...)

#define STL 500

// Set to 1 to debug

int debug = 0;


// Function to extract the extension from a file name (fn)

char* get_extn(char* fn)
{
  // Go to the end

  int count = strlen(fn);

  // Look back for the .

  while((count > 0) && (fn[count] != '.')) count--;

  // Forget the .

  count++;

  // How long is the extension?  Get a new char array that long + a bit

  int len = strlen(&fn[count]);
  char* r = new char[len+1];

  // Copy the extension into the new array

  len = 0;
  do
  {
    r[len] = fn[count];
    len++;
    count++;
  } while(fn[count]);

  // Terminate the string

  r[len] = 0;

  // That's it.

  return r;
}

// Function to generate an output file name given the root, the extension,
// and the x and y counts of an image segment.

char* op_name(char* root, char* extn, int x, int y)
{

  // New string to contain the file name; the 50 is a hack...

  int rl = strlen(root);
  int el = strlen(extn);
  char* r = new char[rl + el + 50];

  // Stick the root on the front

  strcpy(r,root);

  // add on _<x>_<y>; the +1s are so the counts dont start at 0

  ostrstream ost(&r[rl], el+49);
  ost << "_" << x+1 << "_" << y+1 << '\0';

  // add on .extn

  strcat(r,".");
  strcat(r,extn);

  // That's it.

  if(debug) cout << r << endl;
  return r;
}
 

// Main prog - no argc & argv stuff makes it more protable.

int main()
{

  // xdiv, ydiv - number of columns and rows of sub-images
  // xl, yl - size of the original image
  // xi, yi - sizes of the sub images
  // xr, yr - remainder pixels if xl and yl are not perfect multiples of xi and yi

  int xdiv, ydiv, xl, yl, xi, yi, xr, yr;
    
  // Imagemagick catch and throw for errors
  
  try 
  {

    // Ask the user the name of the input file and load it into image

    cout << "Input file name (including extension): ";
    char fn[STL];
    cin >> fn;
    Image image(fn);

    // How big is it?

    xl = image.columns();
    yl = image.rows();

    // The extension, and a pointer to the output file name as will be

    char* extn = get_extn(fn);
    char* op_n;

    // What's to be done?

    cout << "Divide horizontally into (number of image columns): ";
    cin >> xdiv;
    cout << "Divide vertically into (number of image rows): ";
    cin >> ydiv;
    cout << "Output file name root: ";
    cin >> fn;

    // Some division and remainder sums to get image sizes

    xi = xl/xdiv;
    xi--;
    yi = yl/ydiv;
    yi--;
    yr = yl%ydiv;

    // Some numbers to increment as we chop the image up

    int x0;      // Left and...
    int y0 = 0;  // ...top
    int x1, y1;  // Right and bottom

    // The sub images, as will be

    Image* image_xy;

    // X-windows/Imagemagick geometry specifier string

    char i_format[STL];

    // Chop image up; Y is the outer loop...

    for(int y = 0; y < ydiv; y++)
    {
      // Y at the bottom of this sub image

      y1 = y0 + yi;

      // If any Y remainder pixels left, absorb one of them

      if(yr > 0)
      {
	y1++;
	yr--;
      }

      // X remainder and left-hand edge

      xr = xl%xdiv;
      x0 = 0;

      // X is the inner loop

      for(int x = 0; x < xdiv; x++)
      {
	// X at the right of this sub image

	x1 = x0 + xi;

        // If any X remainder pixels left, absorb one of them 

	if(xr > 0)
	{
	  x1++;
	  xr--;
	}

	// Put together the X/Imagemagick geometry descriptor string

	ostrstream ost(i_format, STL);
        ost << x1-x0+1 << "x" << y1-y0+1 << "+" << x0 << "+" << y0 << '\0';
	if(debug) cout << i_format << endl;

	// Copy the original...

	image_xy = new Image(image);

	// ...and crop it to this sub-image

	image_xy->crop(i_format);

	// Cook up the output file name and write out the sub-image

	op_n = op_name(fn, extn, x, y);
	image_xy->write(op_n);

	// No memory leaks please...

	delete image_xy;
	delete [] op_n;

	// Left edge of next image

	x0 = x1 + 1;
      }

      // Top edge of next image

      y0 = y1 + 1;
    }
  }
      
  // Something wicked this way comes...
 
  catch( Exception &error_ )
  {
     cerr << "Caught exception: " << error_.what() << endl;
     return 1;
  }

  // Phew.  Made it.

  return 0;
}                
