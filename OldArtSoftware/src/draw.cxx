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

#define Colour Color

double scl;
double safe; 
double down;
int robot = 0;
int db = 1;

picture *input;
picture *output;
picture *work;

cartesian* c;

// Return 1 if (x_pix, y_pix) is 1 and has at least one 0 neighbour

int b_and_w(int x_pix, int y_pix)
{
  if(input->get_pixel(x_pix, y_pix) < 1) return 0;
  for(int i = -1; i <= 1; i++)
    for(int j = -1; j <= 1; j++)
      if(i || j)
	if(input->get_pixel(x_pix + i, y_pix + j) == 0) return 1;
  return 0;
}

int find_start_r(int &x_pix, int &y_pix)
{
  if(b_and_w(x_pix, y_pix)) return 1;

  work->set_pixel(x_pix, y_pix, 1);

  int x, y, i, j;

  // Breadth-first search - examine all neighbours

  for(i = -1; i <= 1; i++)
  {
    x = x_pix + i;
    for(j = -1; j <= 1; j++)
    {
      if(i || j)
      {
	y = y_pix + j;
	if(work->get_pixel(x, y) == 0)
	{
	  if(b_and_w(x, y)) 
	  {
	    x_pix = x;
	    y_pix = y;
	    return 1;
	  }
	}
      }
    }
  }

  // No neighbour is satisfactory - recurse into them

  for(i = -1; i <= 1; i++)
  {
    x = x_pix + i;
    for(j = -1; j <= 1; j++)
    {
      if(i || j)
      {
	y = y_pix + j;
	if(work->get_pixel(x, y) == 0)
	{
	  if(find_start_r(x, y))
	  {
	    x_pix = x;
	    y_pix = y;
	    return 1;
	  }
	}
      }
    }
  }

  return 0;
}

int find_start(int &x_pix, int &y_pix)
{
  work->erase();
  int i = find_start_r(x_pix, y_pix);
  if(db)
    cout << "find_start(...): x_pix = " << x_pix << ", y_pix = " << y_pix <<
      ", i = " << i << endl;
  return i;
}

void trace_line(int &x_pix, int &y_pix)
{
  if(robot)
  {
    c->quick(x_pix*scl, y_pix*scl, safe);
    c->go(x_pix*scl, y_pix*scl, down);
  } else
    output->set_pixel(x_pix, y_pix, 1);

  int old_x = 0, old_y = -1;
  int i, j, ir, jr, sp, score, x, y;

  do
  {
    score = -3;
    for(i = -1; i <= 1; i++)
    {
      x = x_pix + i;
      for(j = -1; j <= 1; j++)
      {
	if(i || j)
	{
	  y = y_pix + j;
	  if(b_and_w(x, y))
	  {
	    sp = i*old_x + j*old_y;
	    if(sp >=score)
	    {
	      score = sp;
	      ir = i;
	      jr = j;
	    }
	  }
	}
      }
    }

    input->set_pixel(x_pix, y_pix, 0);

    if(score > -3)
    {
      x_pix = x_pix + ir;
      y_pix = y_pix + jr;

      if(robot)
	c->go(x_pix*scl, y_pix*scl, down);
      else
	output->set_pixel(x_pix, y_pix, 1);

      old_x = ir;
      old_y = jr;
      input->set_pixel(x_pix, y_pix, 0);
    }
  } while(score > -3);

  if(robot)
    c->quick(x_pix*scl, y_pix*scl, safe);
}

// Main prog

int main(int argc, char **argv)
{
  int x_pix, y_pix;

  if(argc != 2)
  {
    cerr << "Usage: draw image_file\n";
    return(1);
  }
    
  // Imagemagick catch and throw for errors
  
  try 
  {
    input = new picture(argv[1]);
    work = new picture(argv[1]);

    x_pix = input->width();
    y_pix = input->height();

    cout << "Image is " << x_pix << " pixels by " << y_pix << " pixels.\n" <<
      "What X length do you want in mm? ";
    cin >> scl;
    scl = scl/(double)x_pix;

    cout << "Safe height (mm): ";
    cin >> safe;
    cout << "Down height (mm): ";
    cin >> down;

    if(robot)
    {
      c = new cartesian(0);
      c->quick(0, 0, safe);
      cout << "Bottom left; type any character: ";
      char dummy;
      cin >> dummy;
      c->quick(scl*x_pix, scl*y_pix, safe);
      cout << "Top right; type any character: ";
      cin >> dummy;
    } else
    {
      output = new picture(argv[1]);
      output->erase();
    }

    while(find_start(x_pix, y_pix)) trace_line(x_pix, y_pix);

    if(robot) 
      c->safe();
    else
      output->image()->write("draw_op.gif");

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
