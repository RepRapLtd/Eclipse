/*
 * draw.h
 *
 * A C++/Imagemagick program to draw from a pre-processed graphics file
 * Adrian Bowyer
 * 8 June 2003
 *
 * A.Bowyer@bath.ac.uk
 *
 */

#ifndef DRAW_H
#define DRAW_H

#include <Magick++.h>
#include <iostream.h>

using namespace std;
using namespace Magick;

#define Colour Color
 
class picture
{
private:

  Image *im;
  PixelPacket *impix;
  int xl, yl;

public:

  picture(char* fn)
  {
    im = new Image(fn);
    //xl = im->columns();
    xl = 100;
    yl = im->rows();
    im->classType( DirectClass );
    im->modifyImage();
    impix = im->getPixels(0, 0, xl, yl);
  }

  int outside(int x, int y)
  {
    if (x < 0) return 1;
    if (y < 0) return 1;
    if (x >= xl) return 1;
    if (y >= yl) return 1;
    return 0;
  }

  int get_pixel(int x, int y)
  {
    if (outside(x, y)) return -1;
    Colour c = *(impix + y*xl + x);
    if( c.redQuantum() > MaxRGB/2 ) return 1;
    return 0;
  }

  void set_pixel(int x, int y, int col)
  {
    if (col)
      *(impix + y*xl + x) = Colour(MaxRGB, MaxRGB, MaxRGB);
    else
      *(impix + y*xl + x) = Colour(0, 0, 0);
    im->syncPixels();
  }

  void erase()
  {
    for(int i = 0; i < xl; i++)
      for(int j = 0; j < yl; j++)
	*(impix + j*xl + i) = Colour(0, 0, 0);
    im->syncPixels();
  }

  int width() { return xl; }
  int height() { return yl; }

  Image* image() { return im; }

};

#endif
