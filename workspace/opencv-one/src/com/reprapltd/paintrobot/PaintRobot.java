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
	
	public static void showResult(Mat img) 
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
	        JFrame frame = new JFrame("Frame: " + frameCount);
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
	
	public static void GoWhite(Mat img)
	{
		double[] data = new double[3];
		data[0] = 255;
		data[1] = 255;
		data[2] = 255;
		for(int i = 0; i < img.width(); i++)
			for(int j = 0; j < img.height(); j++)
				img.put(j, i, data);
	}
	
	
	public static void Threshold(Mat image, double percentile)
	{
		Imgproc.threshold(image, image, Percentile(image, percentile), 255.0, 0);
	}
	
	
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
	
	public static void Stroke(Mat image, Point[] pixels, int count, double[] colour, double d)
	{
		for(int i = 0; i < count; i++)
			Blob(image, pixels[i], colour, d);
	}
	
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
	
	public static void Paint(Mat original)
	{	
		double blob = original.width()/30;
		double threshold = 0.7;
		double brushReduction = 0.7;
		
		Mat painting = new Mat();
		original.copyTo(painting); // Deep
		GoWhite(painting);
		
		int strokes = 0;
		boolean keepGoing = true;
		
		Mat result = new Mat();

		//Mat rCopy = new Mat();
		MatOfPoint biggest;
		double[] average = new double[3];
		List<MatOfPoint> contours = new ArrayList<MatOfPoint>();
		Mat hierarchy = new Mat();
		
		Mat d;
		Mat kernel;
		Point[] pixels;
		double[] cols;
		double big;
		double b;
		int track;
		int pix;
		int k;
		
		while(strokes < 6 && keepGoing)
		{
			d = ImageDifference(original, painting);
			if(d == null)
				System.out.println("d null!");
			//showResult(d);
			Threshold(d, threshold);
			//if(strokes%10 == 0)
				//showResult(d);
			kernel = Imgproc.getStructuringElement(Imgproc.MORPH_ELLIPSE, new Size(blob, blob));
			if(kernel == null)
				System.out.println("kernel null!");
			Imgproc.erode(d, result, kernel);
			if(result == null)
				System.out.println("result null!");
			//showResult(result);

			//result.copyTo(rCopy);
			Imgproc.findContours(result, contours, hierarchy, Imgproc.RETR_LIST, Imgproc.CHAIN_APPROX_SIMPLE);
			if(contours == null)
				System.out.println("contours null!");
			//System.out.println("contour count: " + contours.size());
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

			//System.out.println("contour: " + biggest.size());
			//System.out.println("pixels: " + pixels.length);

			track = pixels.length;
			if(track > 20)
			{
				track = track/2;

				average[0] = average[1] = average[2] = 0;
				for(pix = 0; pix < track; pix++)
				{
					cols = original.get((int)Math.round(pixels[pix].y), (int)Math.round(pixels[pix].x));
					for(k = 0; k < 3; k++)
						average[k] += cols[k];
				}
				for(k = 0; k < 3; k++)
					average[k] = average[k]/(double)track;
				Stroke(painting, pixels, track, average, blob);
			}else
			{
				System.out.println("Short contour: " + track + ", blob: " + blob);
				blob = blob*brushReduction;
				if(blob < 3)
					keepGoing = false;
			}
			strokes++;
			//if(strokes%20 == 0)
				//showResult(painting);
		}
		showResult(painting);
	}

	public static void main(String[] args)
	{
		// load the OpenCV native library
		System.loadLibrary(Core.NATIVE_LIBRARY_NAME);
		frameCount = 0;

		//Mat original = GetPicture();
		Mat original = Imgcodecs.imread("resources/test.png");
		showResult(original);
		Paint(original);

        
		// write the new image on disk
        
		//Imgcodecs.imwrite("resources/result.jpg", result);
		
		// display the result
		
		
		

		
		
		System.out.println("Done!");
		//System.exit(0);
	}
	
}