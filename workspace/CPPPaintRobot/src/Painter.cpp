/*
 * Painter.cpp
 *
 *  Created on: 24 Jun 2019
 *      Author: ensab
 */


//============================================================================
// Name        : Painter.cpp
// Author      : Adrian Bowyer
// Version     :
// Copyright   : GPL
// Description : Painting robot control program
//============================================================================

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <algorithm>
using namespace std;

#include "Painter.h"

uchar warmGrey[3] = {79, 64, 64};
Scalar red = Scalar(0, 0, 255);
Vec3b black = Vec3b(0, 0, 0);
Vec3b white = Vec3b(255, 255, 255);

//void Painter::ShowResult(Mat img, char* title)
//{
//	//char text[100];
//	//sprintf(text,"%s %d %s %d", title, frameCount, ", strokes: ", strokeCount);
//
//	namedWindow(title, CV_WINDOW_AUTOSIZE);
//	//moveWindow(title, 100+frameXOff, 100+frameYOff);
//
//	//frameXOff += 50;
//	//frameYOff += 50;
//
//	Mat copy = img;
//
//	imshow(title, copy);
//	waitKey(10);
//
//	// release the image
//	//cvReleaseImage(&img );
//	frameCount++;
//}



uchar Painter::Percentile(Mat img, double percentile)
{
	int histogram[MAX_PIXEL];
	for(int i = 0; i < MAX_PIXEL; i++)
		histogram[i] = 0;
	uchar p;

	for(int x = 0; x < imageWidth; x++)
		for(int y = 0; y < imageHeight; y++)
		{
			p = GetPixelMono(img, x, y);
			histogram[p]++;
		}

	//for(int i = 0; i < MAX_PIXEL; i++)
	//	cout << histogram[i] << ' ';
	//cout << endl;

	double frac = percentile*pixelCount;
	//cout << "frac: " << frac << endl;

	double total = 0;

	for(int i = 0; i < MAX_PIXEL; i++)
	{
		total += (double)histogram[i];
		if(total >= frac)
			return i;
	}

	cerr << "Can't find percentile in histogram!" << endl;

	return MAX_PIXEL;
}


void Painter::GoColour(Mat img, uchar p[])
{
	for(int x = 0; x < imageWidth; x++)
		for(int y = 0; y < imageHeight; y++)
			PutPixel(img, x, y, p);
}


void Painter::Threshold(Mat image, double percentile)
{
	double t = (double)Percentile(image, percentile);
	if(debug == STRONG_DEBUG)
		cout << "Threshold: " << t << endl;
	threshold(image, thresholded, t, 255.0, THRESH_BINARY);
}




void Painter::Blob(Mat image, Point pixel, uchar colour[], double d)
{
	double r2 = d*d/4;
	for(int y = 0; y < round(d/2); y++)
	{
		int xm = round(sqrt(r2 - y*y));
		for(int x = 0; x < xm; x++)
		{
			PutPixel(image, pixel.x + x, pixel.y + y, colour);
			PutPixel(image, pixel.x - x, pixel.y + y, colour);
			PutPixel(image, pixel.x + x, pixel.y - y, colour);
			PutPixel(image, pixel.x - x, pixel.y - y, colour);
		}
	}
}


void Painter::Stroke(Mat image, int count, uchar colour[], double d)
{
	for(int i = 0; i < count; i++)
		Blob(image, pixels[i], colour, d+1);
}


uint Painter::PixelDifference2(uchar a[], uchar b[])
{
	uint result = 0;
	for(int k = 0; k < 3; k++)
	{
		int diff = a[k] - b[k];
		result += diff*diff;
	}
	return result;
}


uchar Painter::PixelDifference(uchar a[], uchar b[])
{
	uint d2 = PixelDifference2(a, b);
	d2 = sqrt((double)d2)*255.0/441.67295593;  // Scaling factor from 3D to 1D
	return (uchar)d2;
}


void Painter::ImageDifference(Mat a, Mat b)
{
	uchar dataA[3];
	uchar dataB[3];
	uchar pixel;
	for(int x = 0; x < imageWidth; x++)
		for(int y = 0; y < imageHeight; y++)
		{
			GetPixel(a, x, y, dataA);
			GetPixel(b, x, y, dataB);
			pixel = PixelDifference(dataA, dataB);
			PutPixelMono(difference, x, y, pixel);
		}
}

void Painter::CopyBWToColour(Mat source, Mat dest)
{
	for(int x = 0; x < source.cols; x++)
		for(int y = 0; y < source.rows; y++)
		{
			Vec3b colour = source.at<Vec3b>(Point(x,y));
			if(colour[0] > 0)
				dest.at<Vec3b>(Point(x,y)) = white;
			else
				dest.at<Vec3b>(Point(x,y)) = black;
		}
}


void Painter::RefreshWindows()
{
	CopyBWToColour(eroded, lastTrack);
	drawContours(lastTrack, contours, largestContourIndex, red);

	imshow(originalWindow, original);
	waitKey(10); //???
	imshow(paintingWindow, painting);
	waitKey(10); //???
	imshow(lastTrackWindow, lastTrack);
	waitKey(10); //???
	imshow(differenceWindow, difference);
	waitKey(10); //???
	imshow(thresholdedWindow, thresholded);
	waitKey(10); //???
	imshow(erodedWindow, eroded);
	waitKey(10); //???
}


Painter::Painter(const Mat org, uchar db, double tv, double br, double bl, uchar underCoat[])
{
	debug = db;

	frameXOff = 0;
	frameYOff = 0;

	org.copyTo(original);
	originalWindow = "Original";
	namedWindow(originalWindow, CV_WINDOW_AUTOSIZE);
	imshow(originalWindow, original);
	waitKey(10); //???

	imageWidth = original.cols;
	imageHeight = original.rows;
	pixelCount = imageWidth*imageHeight;

	original.copyTo(painting);
	paintingWindow = "Painting";
	namedWindow(paintingWindow, CV_WINDOW_AUTOSIZE);
	imshow(paintingWindow, painting);
	waitKey(10); //???

	cvtColor(original, difference, CV_BGR2GRAY);
	differenceWindow = "Difference";
	namedWindow(differenceWindow, CV_WINDOW_AUTOSIZE);
	imshow(differenceWindow, difference);
	waitKey(10); //???

	difference.copyTo(thresholded);
	thresholdedWindow = "thresholded";
	namedWindow(thresholdedWindow, CV_WINDOW_AUTOSIZE);
	imshow(thresholdedWindow, thresholded);
	waitKey(10); //???

	difference.copyTo(eroded);
	erodedWindow = "eroded";
	namedWindow(erodedWindow, CV_WINDOW_AUTOSIZE);
	imshow(erodedWindow, eroded);
	waitKey(10); //???

	original.copyTo(lastTrack);
	lastTrackWindow = "Last Track";
	namedWindow(lastTrackWindow, CV_WINDOW_AUTOSIZE);
	imshow(lastTrackWindow, lastTrack);
	waitKey(10); //???

	for(int k = 0; k < 3; k++)
		average[k] = 0;

	if(underCoat)
	{
		GoColour(painting, underCoat);
	} else
	{
		for(int x = 0; x < imageWidth; x++)
			for(int y = 0; y < imageHeight; y++)
			{
				GetPixel(original, x, y, pixel);
				for(int k = 0; k < 3; k++)
					average[k] += (double)pixel[k];
			}

		for(int k = 0; k < 3; k++)
			pixel[k] = (uchar)(average[k]/(double)pixelCount);

		GoColour(painting, pixel);
	}

	frameCount = 0;
	brushDiameter = imageWidth/30;
	thresholdValue = tv;
	brushReduction = br;
	bleed = bl;

	strokeCount = 0;
	totalStrokes = 0;
	keepGoing = true;
}


void Painter::Report()
{
	cout << endl << "Painting report - current state." << endl;
	cout << " Image size: " << imageWidth << " x " << imageHeight << endl;
	cout << " stroke count: " << totalStrokes << endl;
	cout << " frame count: " << frameCount << endl;
	cout << " brush diameter: " << brushDiameter << endl;
	cout << " threshold value: " << thresholdValue << endl;
	cout << " brush reduction: " << brushReduction << endl;
	cout << " bleed: " << bleed << endl;
	cout << " track length: " << track << endl;
	cout << " debug setting: ";
	if(debug == WEAK_DEBUG)
		cout << "weak";
	else if(debug == STRONG_DEBUG)
		cout << "strong";
	else
		cout << "none";
	cout << endl;
	cout << " smallest contour: " << smallestContour << endl;
	cout << " biggest contour: " << biggestContour << endl;
	cout << " mean contour: " << meanContour << endl;
	cout << " contour standard deviation: " << sdContour << endl;
	cout << endl;

}

void Painter::PaintABit(int strokes)
{
	strokeCount = 0;
	keepGoing = true;

	Size *brushSize = new Size(brushDiameter, brushDiameter);
	brush = getStructuringElement(MORPH_ELLIPSE, *brushSize);

	while(strokeCount < strokes && keepGoing)
	{
		ImageDifference(original, painting);
		Threshold(difference, thresholdValue);
		erode(thresholded, eroded, brush);

		findContours(eroded, contours, hierarchy, RETR_LIST, CHAIN_APPROX_NONE);

		largestContourIndex = -1;

		smallestContour = DBL_MAX;
		biggestContour = 0;
		meanContour = 0.0;
		sdContour = 0.0;
		for(int i = 0; i < contours.size(); i++)
		{
			double area = contourArea(contours[i], false);
			if(area > biggestContour)
			{
				biggestContour = area;
				largestContourIndex=i;
			}
			if(area < smallestContour)
			{
				smallestContour = area;
			}
			meanContour += area;
			sdContour += area*area;
		}
		double s = 1.0/(double)contours.size();
		meanContour = meanContour*s;
		sdContour = sqrt(sdContour*s - meanContour*meanContour);

		if(largestContourIndex < 0)
		{
			cout << "biggest null!" << endl;
			return;
		}


		pixels = contours[largestContourIndex];

		track = pixels.size();
		if(track > 20)
		{
			track = track/2;

			average[0] = average[1] = average[2] = 0;
			for(pix = 0; pix < track; pix++)
			{
				GetPixel(original, pixels[pix].x, pixels[pix].y, pixel);
				for(int k = 0; k < 3; k++)
					average[k] += (double)pixel[k];
			}
			for(int k = 0; k < 3; k++)
				pixel[k] = (uchar)(average[k]/(double)track);
			Stroke(painting, track, pixel, brushDiameter*bleed);
		}else
		{
			if(debug >= WEAK_DEBUG)
				cout << "Short contour: " << track << ", brush diameter: " << brushDiameter << endl;

			brushDiameter = brushDiameter*brushReduction;

			if(brushDiameter < MIN_BRUSH)
			{
				keepGoing = false;
				cout << "Stopping, because brush diameter < " << MIN_BRUSH << endl;
			} else
			{
				delete brushSize;
				brushSize = new Size(brushDiameter, brushDiameter);
				brush = getStructuringElement(MORPH_ELLIPSE, *brushSize);
			}
		}
		strokeCount++;
	}
	totalStrokes += strokeCount;
}

void Painter::Prompt()
{
	cout << "Commands: " << endl;
}

void Painter::Control()
{
	cout << "Type h for help." << endl;
	while(1)
	{
		RefreshWindows();
		cout << "Command: ";
		char c;
		cin >> c;

		switch(c)
		{
		case 'p':
			cout << "Number of strokes to paint: ";
			int num;
			cin >> num;
			PaintABit(num);
			break;

		case 'r':
			Report();
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

int main(int argc, char *argv[])
{
	Mat org = cvLoadImage("/home/ensab/rrlOwncloud/RepRapLtd/Engineering/Software/Eclipse/workspace/JavaPaintRobot/resources/sunset.jpg");
	Painter* p = new Painter(org, WEAK_DEBUG, 0.7, 0.7, 1.1, warmGrey);
	p->Control();
	p->Show();

	return 0;
}




