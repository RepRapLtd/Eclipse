package com.reprapltd.paintrobot;

import java.awt.image.BufferedImage;
import java.io.ByteArrayInputStream;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.imageio.ImageIO;
import javax.swing.ImageIcon;
import javax.swing.JFrame;
import javax.swing.JLabel;

import org.opencv.core.CvType;
import org.opencv.core.Mat;
import org.opencv.core.MatOfByte;
import org.opencv.core.Size;
import org.opencv.imgcodecs.Imgcodecs;
import org.opencv.imgproc.Imgproc;
import org.opencv.videoio.VideoCapture;
import org.opencv.core.Core;
import org.opencv.core.MatOfPoint;
import org.opencv.core.Scalar;
import org.opencv.core.Point;

/**
 * Painting Robot
 * 
 * @author <a href="https://reprapltd.com/contact-us/">Adrian Bowyer</a>
 * @since 2016-09-27
 * 
 */
public class PaintRobot
{
	public static int frameCount;
	public static boolean debug = false;
	public static double[] traceColour = new double[3];
	public static double[] black = new double[3];
	public static double[] white = new double[3];
	
	public static void showResult(Mat img, String label) 
	{
	    Imgproc.resize(img, img, new Size(640, 480));
	    MatOfByte matOfByte = new MatOfByte();
	    Imgcodecs.imencode(".jpg", img, matOfByte);
	    byte[] byteArray = matOfByte.toArray();
	    BufferedImage bufImage = null;
	    try 
	    {
	        InputStream in = new ByteArrayInputStream(byteArray);
	        bufImage = ImageIO.read(in);
	        JFrame frame = new JFrame("Frame: " + frameCount + ", Label: " + label);
	        frame.getContentPane().add(new JLabel(new ImageIcon(bufImage)));
	        frame.pack();
	        frame.setVisible(true);
	    } catch (Exception e) 
	    {
	        e.printStackTrace();
	    }
	    frameCount++;
	}
	
	
	public static double Percentile(double[] array, double percentile)
	{
		Arrays.sort(array);
		return array[(int)Math.round(((double)array.length*percentile))];		
	}
	
	public static double Percentile(Mat img, double percentile)
	{
		int count = 0;
		int pixelCount = img.width()*img.height();
		double[] array = new double[pixelCount];
		for(int i = 0; i < img.width(); i++)
			for(int j = 0; j < img.height(); j++)
			{
				double pixel[] = img.get(j, i);
				//double magnitude = Math.sqrt(pixel[0]*pixel[0] + pixel[1]*pixel[1] + pixel[2]*pixel[2]);
				array[count] = pixel[0];
				count++;
			}
		return Percentile(array, percentile);
	}
	
	/***
	 * Get an image from the webcam.
	 * 
	 * @return
	 */
	public static Mat GetPicture()
	{
		String cameraLocation = "/dev/video0";
	    VideoCapture camera = new VideoCapture(cameraLocation); // open the default camera 
	    if(!camera.isOpened())  // check if we succeeded
	    {
	    	System.out.println("Cannot open " + cameraLocation);
	    	System.exit(1);
	    }
	    	    
	    
	    Mat img = new Mat();
	    
	    // Throw away rubbish frames at the start
	    
	    for(int i = 0; i < 5; i++)
	    {
	    	camera.read(img);
	    }

	    camera.release();
	    
	    return img;
	}
	
	/***
	 * Set an entire image to a single colour
	 * 
	 * @param img
	 */
	public static void GoColour(Mat img, double[] c)
	{
		for(int i = 0; i < img.width(); i++)
			for(int j = 0; j < img.height(); j++)
				img.put(j, i, c);
	}
	
	/***
	 * Threshold an image at percentile.
	 * 
	 * @param image
	 * @param percentile
	 */
	public static void Threshold(Mat image, double percentile)
	{
		Imgproc.threshold(image, image, Percentile(image, percentile), 255.0, 0);
	}
	
	/***
	 * Paint a disc.
	 * 
	 * @param image
	 * @param pixel
	 * @param colour
	 * @param d
	 */
	public static void Blob(Mat image, Point pixel, double[] colour, double d)
	{
		double r2 = d*d/4;
		for(int y = 0; y < (int)Math.round(d/2); y++)
		{
			int xm = (int)Math.round(Math.sqrt(r2 - y*y));
			for(int x = 0; x < xm; x++)
			{
				int px = (int)Math.round(pixel.x);
				int py = (int)Math.round(pixel.y);
				image.put(py + y, px + x, colour);
				image.put(py + y, px - x, colour);
				image.put(py - y, px + x, colour);
				image.put(py - y, px - x, colour);
			}
		}
	}
	
	/**
	 * Paint a stroke along the contour pixels[] for count of them (may not go to the end).
	 * 
	 * @param image
	 * @param pixels
	 * @param count
	 * @param colour
	 * @param d
	 */
	public static void Stroke(Mat image, Point[] pixels, int count, double[] colour, double d)
	{
		for(int i = 0; i < count; i++)
			Blob(image, pixels[i], colour, d);
//		if(debug)
//		{
//			for(int i = 0; i < count; i++)
//			{
//				int px = (int)Math.round(pixels[i].x);
//				int py = (int)Math.round(pixels[i].y);
//				image.put(py, px, traceColour);				
//			}
//		}
	}
	
	/***
	 * The squared difference between two pixels
	 * 
	 * @param a
	 * @param b
	 * @return
	 */
	public static double PixelDifference2(double[] a, double[] b)
	{
		double result = 0;
		for(int k = 0; k < 3; k++)
		{
			double diff = a[k] - b[k];
			result += diff*diff;
		}
		return result;		
	}
	
	/***
	 * The difference between two pixels normalised to lie in [0, 255] (the maximum 
	 * squared difference is 3 x 255^2, the square root of which is 441.6...).
	 * @param a
	 * @param b
	 * @return
	 */
	public static double PixelDifference(double[] a, double[] b)
	{
		return Math.sqrt(PixelDifference2(a, b))*255.0/441.67295593;
	}
	
	public static Mat ImageDifference(Mat a, Mat b)
	{
		if(a.width() != b.width() || a.height() != b.height())
		{
			System.out.println("Subtracted images not the same size!");
			return null;
		}
		Mat result = new Mat();
		Imgproc.cvtColor(a, result, Imgproc.COLOR_BGR2GRAY);
		double[] pixel = new double[1];
		for(int i = 0; i < result.width(); i++)
			for(int j = 0; j < result.height(); j++)
			{
				double dataA[] = a.get(j, i);
				double dataB[] = b.get(j, i);
				pixel[0] = PixelDifference(dataA, dataB);
				result.put(j, i, pixel);
			}
		return result;
	}
	
	/**
	 * Work out the average colour of an entire image
	 * 
	 * @param image
	 * @return
	 */
	public static double[] AverageColour(Mat image)
	{
		double[] a = new double[3];
		a[0] = 0.0;
		a[1] = 0.0;
		a[2] = 0.0;
		for(int i = 0; i < image.width(); i++)
			for(int j = 0; j < image.height(); j++)
			{
				double pixel[] = image.get(j, i);
				a[0] += pixel[0];
				a[1] += pixel[1];
				a[2] += pixel[2];
			}
		double inv = 1.0/(image.width()*image.height());
		a[0] = a[0]*inv;
		a[1] = a[1]*inv;
		a[2] = a[2]*inv;		
		return a;
	}
	
	public static double[] TrackAverageColour(Mat original, Point[] pixels)
	{
		double[] result = new double[3];
		double[] cols;
		int k;
		result[0] = result[1] = result[2] = 0;
		int x, y, xOld, yOld;
		xOld = yOld = -1;
		for(int pix = 0; pix < pixels.length; pix++)
		{
			x = (int)Math.round(pixels[pix].x);
			y = (int)Math.round(pixels[pix].y);
//			if(xOld > -1)
//			{
//				xOld = Math.abs(x - xOld);
//				yOld = Math.abs(y - yOld);
//				if(xOld > 1 || yOld > 1)
//					System.out.println("Pixel break: " + pix + " " + xOld + " " + yOld);
//			}
			cols = original.get(y, x);
			for(k = 0; k < 3; k++)
				result[k] += cols[k];
			xOld = x;
			yOld = y;
		}
		for(k = 0; k < 3; k++)
			result[k] = result[k]/(double)pixels.length;
		
		return result;
	}
	
	
	
	public static void ShowTrack(Mat d, Point[] pixels)
	{
		for(int pix = 0; pix < pixels.length; pix++)
		{
			int px = (int)Math.round(pixels[pix].x);
			int py = (int)Math.round(pixels[pix].y);
			d.put(py, px, traceColour);
		}
		showResult(d, "paint track");
	}
	
	/***
	 * Set an entire image to a single colour
	 * 
	 * @param img
	 */
	public static void CopyBWToColour(Mat source, Mat dest)
	{
		for(int x = 0; x < source.width(); x++)
			for(int y = 0; y < source.height(); y++)
			{
				double pixel[] = source.get(y, x);
				if(pixel[0] > 0)
					dest.put(y, x, white);
				else
					dest.put(y, x, black);
			}
	}
	
	/***
	 * Paint the picture.
	 * @param original
	 */
	public static void Paint(Mat original)
	{	
		double blob = original.width()/30;
		double threshold = 0.7;
		double brushReduction = 0.7;
		
		Mat painting = new Mat();
		original.copyTo(painting); // Deep copy
		double[] average = AverageColour(painting);
		GoColour(painting, average);
		
		int strokes = 0;
		boolean keepGoing = true;
		
		Mat result = new Mat();

		//Mat rCopy = new Mat();
		MatOfPoint biggest;
		
		Mat d = new Mat();
		Mat tracks = new Mat();
		Mat kernel = new Mat();
		Point[] pixels;
		double big;
		double b;
		int track;
		int k;
		Mat hierarchy = new Mat();
		List<MatOfPoint> contours = new ArrayList<MatOfPoint>();
		
		if(debug)
			original.copyTo(tracks); // Deep copy
		
		while(strokes < 10 && keepGoing)
		{
			d.release();
			d = ImageDifference(original, painting);
			if(d == null)
				System.out.println("d null!");
			if(debug)
				showResult(d, "original - painting");
			Threshold(d, threshold);
			if(debug)
				showResult(d, "original - painting thresholded");
			kernel.release();
			kernel = Imgproc.getStructuringElement(Imgproc.MORPH_ELLIPSE, new Size(blob, blob));
			if(kernel == null)
				System.out.println("kernel null!");
			result.release();
			Imgproc.erode(d, result, kernel);
			if(result == null)
				System.out.println("result null!");
			if(debug)
				showResult(result, "original - painting thresholded & eroded");

			if(debug)
			{
				tracks.release();
				CopyBWToColour(result, tracks);
			}

			hierarchy.release();
			hierarchy = new Mat();
			contours = new ArrayList<MatOfPoint>();
			Imgproc.findContours(result, contours, hierarchy, Imgproc.RETR_LIST, Imgproc.CHAIN_APPROX_NONE);

			if(contours == null)
				System.out.println("contours null!");
			if(debug)
				System.out.println("contour count: " + contours.size());
			big = 0;

			biggest = null;
			for(MatOfPoint c : contours)
			{
				if((b = Imgproc.contourArea(c)) > big)
				{
					big = b;
					biggest = c;
				}
			}
			if(biggest == null)
				System.out.println("biggest null!");

			//List<MatOfPoint> bigContour = new ArrayList<MatOfPoint>();
			//bigContour.add(biggest);
			//Imgproc.drawContours(painting, bigContour, 0, new Scalar(0, 255, 0), -1);
			//showResult(result);
			//showResult(original);

			pixels = biggest.toArray();
			if(pixels == null)
				System.out.println("pixels null!");
			
			if(debug)
				ShowTrack(tracks, pixels);

			//System.out.println("contour: " + biggest.size());
			//System.out.println("pixels: " + pixels.length);

			track = pixels.length;
			if(track > 20)
			{
				average = TrackAverageColour(original, pixels);
				Stroke(painting, pixels, track, average, blob);
			}else
			{
				System.out.println("Short contour: " + track + ", blob: " + blob);
				blob = blob*brushReduction;
				if(blob < 3)
					keepGoing = false;
			}
			strokes++;
			//if(debug || strokes%10 == 0)
			//	showResult(painting, "painting at" + strokes + " strokes");
			System.gc();
		}
		System.out.println("Strokes: " + strokes);
		showResult(painting, "final painting");
	}

	public static void main(String[] args)
	{
		traceColour[0] = 0.0;
		traceColour[1] = 0.0;
		traceColour[2] = 255.0;
		black[0] = 0.0;
		black[1] = 0.0;
		black[2] = 0.0;
		white[0] = 255.0;
		white[1] = 255.0;
		white[2] = 255.0;
		// load the OpenCV native library
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		frameCount = 0;

		//Mat original = GetPicture();
		//Mat original = Imgcodecs.imread("resources/test.png");
		Mat original = Imgcodecs.imread("resources/sunset.jpg");
		showResult(original, "original");
		Paint(original);

        
		// write the new image on disk
        
		//Imgcodecs.imwrite("resources/result.jpg", result);
		
		// display the result
		
		
		

		
		
		System.out.println("Done!");
		//System.exit(0);
	}
	
}