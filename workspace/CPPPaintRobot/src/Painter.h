/*
 * Painter.h
 *
 *  Created on: 3 Oct 2016
 *  This version: about 24 Jun 2019
 *
 *
 *      Author: Adrian Bowyer
 *      		RepRap Ltd
 *      		http://reprapltd.com
 *
 *      Licence: GPL
 */

#ifndef PAINTER_H_
#define PAINTER_H_

#include <cv.hpp>
#include <core/core.hpp>
#include <highgui/highgui.hpp>
#include <imgproc/imgproc.hpp>


using namespace cv;

#define MAX_PIXEL 255
#define STRONG_DEBUG 2
#define WEAK_DEBUG 1
#define MIN_BRUSH 2

enum SnakeState
{
	sdTooBig = 1,
	mayNotContinue = 2,
	mayContinue = 3
};

class Painter
{
public:
	Painter(const Mat original, uchar db, double tv, double br, double bl, int lt, double as, uchar underCoat[]);
	void Control();
	void Show();
	void ImageToSlink(vector< vector<float> > &result, Mat image);

private:

	uchar GetPixelMono(const Mat img, int x, int y);
	void CopyBWToColour(Mat src, Mat dest);
	void PutPixelMono(Mat img, int x, int y, uchar p);
	void GetPixel(const Mat img, int x, int y, uchar p[]);
	void PutPixel(Mat img, int x, int y, uchar p[]);
	uchar Percentile(const Mat img, double percentile);
	void GoColour(Mat img, uchar p[]);
	void Threshold(Mat image, double percentile);
	void Blob(Mat image, Point pixel, uchar colour[], double d);
	void Stroke(Mat image, int count, uchar colour[], double d);
	uint PixelDifference2(uchar a[], uchar b[]);
	uchar PixelDifference(uchar a[], uchar b[]);
	void ImageDifference(Mat a, Mat b);
	void ConvexQuadrilateral(Point corner[]);
	void Report();
	void PaintABit(int strokes);
	//void ShowResult(const Mat img, char* title);
	void ChangeDebug(uchar db);
	void Prompt();
	void RefreshWindows();
	void FindPixelColourAlongTrack();
	void FindBiggestContour();
	bool TrackLength();
	void SmoothSnake();
	SnakeState GrowHeadOnePixel();
	SnakeState EatTailOnePixel();
	void Coarsen();

	uchar coarsenMask;
	int pixelCount;
	int imageWidth;
	int imageHeight;
	double brushDiameter;
	double thresholdValue;
	double brushReduction;
	double bleed;
	double aspect;
	int strokeCount;
	int totalStrokes;
	bool keepGoing;
	int largestContourIndex;
	int track;
	int startTrack;
	int endTrack;
	int pix;
	int trackTooShort;
	uchar debug;
	double smallestContour;
	double biggestContour;
	double meanContour;
	double sdContour;

	Mat original;
	string originalWindow;
	Mat difference;
	string differenceWindow;
	Mat thresholded;
	string thresholdedWindow;
	Mat painting;
	string paintingWindow;
	Mat eroded;
	string erodedWindow;
	Mat lastTrack;
	string lastTrackWindow;
	Mat brush;
	vector<vector<Point> > contours;
	vector<Vec4i> hierarchy;
	double average[3];
	double sd[3];
	uchar pixel[3];
	vector<Point> pixels;
};

inline uchar Painter::GetPixelMono(Mat img, int x, int y)
{
	if(x < 0 || y < 0 || x >= imageWidth || y >= imageHeight)
	{
		if(debug == STRONG_DEBUG)
			cout << "Reading outside the frame." << endl;
		return 0;
	}
	int pos = img.channels()*(img.cols*y + x);
	return img.data[pos];
}

inline void Painter::PutPixelMono(Mat img, int x, int y, uchar p)
{
	if(x < 0 || y < 0 || x >= imageWidth || y >= imageHeight)
	{
		if(debug == STRONG_DEBUG)
			cout << "Writing outside the frame." << endl;
		return;
	}
	int pos = img.channels()*(img.cols*y + x);
	img.data[pos] = p;
}

inline void Painter::GetPixel(Mat img, int x, int y, uchar p[])
{
	if(x < 0 || y < 0 || x >= imageWidth || y >= imageHeight)
	{
		if(debug == STRONG_DEBUG)
			cout << "Reading outside the frame." << endl;
		p[0] = p[1] = p[2] = 0;
		return;
	}
	int pos = img.channels()*(img.cols*y + x);
	p[0] = img.data[pos];
	p[1] = img.data[pos + 1];
	p[2] = img.data[pos + 2];
}

inline void Painter::PutPixel(Mat img, int x, int y, uchar p[])
{
	if(x < 0 || y < 0 || x >= imageWidth || y >= imageHeight)
	{
		if(debug == STRONG_DEBUG)
			cout << "Writing outside the frame." << endl;
		return;
	}
	int pos = img.channels()*(img.cols*y + x);
	img.data[pos] = p[0];
	img.data[pos + 1] = p[1];
	img.data[pos + 2] = p[2];
}

inline void Painter::Show()
{
	imshow(paintingWindow, painting);
	waitKey(10); //???
}

inline void Painter::ChangeDebug(uchar db)
{
	debug = db;
}


#endif /* PAINTER_H_ */




