/*
 * draw.cxx
 *
 * A C++/Imagemagick program to draw from a pre-processed graphics file
 * Adrian Bowyer
 * 8 June 2003
 *
 * A.Bowyer@bath.ac.uk
 *
 */

#include <draw.h>
#include <cartesian.h>


// Main prog

int main(int argc, char **argv)
{
  int x_pix, y_pix;
  char resp;
  char file_name[200];
  cartesian* c;

      c = new cartesian(0);
      cout << "Debug robot? ";
      cin >> resp;
      if(resp == 'y' || resp == 'Y')
      {
	cout << "Debug file name: ";
	cin >> file_name;
	c->debug(file_name);
      }
      long count = 0;
      while(1)
      {
	count++;
	c->quick(50, 100, 50);
	c->go(50, 100, 35);
	c->go(50,50, 35);
	c->go(100, 50, 35);
	c->quick(100, 50, 50);
	c->quick(50, 100, 50);
	c->go(50, 100, 35);
	c->go(100, 100, 35);
	c->go(100, 50, 35);
	c->quick(100, 50, 50);
	cout << count << "                    \r";
	cout.flush();
      }

      c->safe();

      delete c;

  return 0;
}                
