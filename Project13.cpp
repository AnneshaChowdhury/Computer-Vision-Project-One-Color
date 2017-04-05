#include "opencv2/highgui.hpp"
#include <iostream>
using namespace cv;
using namespace std;

double invGamma(double v) {

	if (v < 0.03928)
	{
		return v / 12.92;
	}
	else
	{
		return pow(((v + 0.055) / 1.055), 2.4);
	}
}

double gamma(double D) {

	if (D < 0.00304)
	{
		return D * 12.92;
	}
	else
	{
		return (1.055 * pow(D, (1 / 2.4)) - 0.055);
	}
}

void runOnWindow(int W1, int H1, int W2, int H2, Mat inputImage, char *outName) {
	int rows = inputImage.rows;
	int cols = inputImage.cols;

	vector<Mat> i_planes;
	split(inputImage, i_planes);
	Mat iB = i_planes[0];
	Mat iG = i_planes[1];
	Mat iR = i_planes[2];

	// dynamically allocate RGB arrays of size rows x cols
	double** L = new double*[rows];
	double** u = new double*[rows];
	double** v = new double*[rows];
	int** R = new int*[rows];
	int** G = new int*[rows];
	int** B = new int*[rows];

	for (int i = 0; i < rows; i++) {
		L[i] = new double[cols];
		u[i] = new double[cols];
		v[i] = new double[cols];
		R[i] = new int[cols];
		G[i] = new int[cols];
		B[i] = new int[cols];
	}

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			double r = iR.at<uchar>(i, j);
			double g = iG.at<uchar>(i, j);
			double b = iB.at<uchar>(i, j);

			double rPrime = r / 255;
			double gPrime = g / 255;
			double bPrime = b / 255;

			double R = invGamma(rPrime);
			double G = invGamma(gPrime);
			double B = invGamma(bPrime);

			double X = 0.412453 * R + 0.35758 * G + 0.180423 * B;
			double Y = 0.212671 * R + 0.71516 * G + 0.072169 * B;
			double Z = 0.019334 * R + 0.119193 * G + 0.950227 * B;

			double Xw = 0.95;
			double Yw = 1.0;
			double Zw = 1.09;

			double uw = (4 * Xw) / (Xw + 15 * Yw + 3 * Zw);
			double vw = (9 * Yw) / (Xw + 15 * Yw + 3 * Zw);
			double t = Y / Yw;

			if (t > 0.008856) {
				L[i][j] = 116 * pow(t, (1 / 3)) - 16;
			}
			else {
				L[i][j] = 903.3 * t;
			}

			double d = X + 15 * Y + 3 * Z;
			double uPrime = 4 * X / d;
			double vPrime = 9 * Y / d;
			u[i][j] = 13 * L[i][j] * (uPrime - uw);
			v[i][j] = 13 * L[i][j] * (vPrime - vw);

		}

	double Lmin = 1000;
	double Lmax = -1;
	int numI = 0;
	int *histogram = new int[256];
	for (int i = 0; i < 256; i++) {
		histogram[i] = 0;
	}

	for (int i = H1; i <= H2; i++)
		for (int j = W1; j <= W2; j++) {
			if (L[i][j] < Lmin) {
				Lmin = L[i][j];
			}
			else if (L[i][j] > Lmax) {
				Lmax = L[i][j];
			}

			histogram[(int)L[i][j]]++;
			numI++;
		}

	int *histogram_map = new int[(int)Lmax + 1];
	for (int i = 0; i < ((int)Lmax + 1); i++) {
		histogram_map[i] = 0;
	}

	// equalization
	int prev_f = 0;
	for (int i = 0; i <= (int)Lmax; i++) {
		int cumulative_f = prev_f + histogram[i];
		histogram_map[i] = (prev_f + cumulative_f) * 101 / (2 * numI);
		prev_f = cumulative_f;
	}

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			double l = L[i][j];
			if (L[i][j] <= Lmin) {
				l = 0;
			}
			else if (L[i][j] >= Lmax) {
				l = 100;
			}
			else {
				l = histogram_map[(int)l];
			}

			double Xw = 0.95;
			double Yw = 1.0;
			double Zw = 1.09;

			double uw = (4 * Xw) / (Xw + 15 * Yw + 3 * Zw);
			double vw = (9 * Yw) / (Xw + 15 * Yw + 3 * Zw);

			double uPrime = (u[i][j] + (13 * uw * L[i][j])) / (13 * L[i][j]);
			double vPrime = (v[i][j] + (13 * vw * L[i][j])) / (13 * L[i][j]);

			double X;
			double Y;
			double Z;

			if (L[i][j] > 7.9996)
			{
				Y = pow(((L[i][j] + 16) / (int)116), 3) * Yw;

			}
			else
			{
				Y = L[i][j] / 903.3 * Yw;
			}

			if (vPrime != 0)
			{
				X = Y * 2.25 * uPrime / vPrime;
				Z = Y * (3 - 0.75 * uPrime - 5 * vPrime) / vPrime;
			}
			else
			{
				X = 0;
				Z = 0;
			}

			/* XYZ to linear sRGB */
			double R2 = 3.240479 * X - 1.53715 * Y - 0.498535 * Z;
			double G2 = -0.969256 * X + 1.875991 * Y + 0.041556 * Z;
			double B2 = 0.055648 * X - 0.204043 * Y + 1.057311 * Z;

			if (R2 > 1) { R2 = 1; }
			else if (R2 < 0) { R2 = 0; }
			if (G2 > 1) { G2 = 1; }
			else if (G2 < 0) { G2 = 0; }
			if (B2 > 1) { B2 = 1; }
			else if (B2 < 0) { B2 = 0; }

			R2 = gamma(R2);
			G2 = gamma(G2);
			B2 = gamma(B2);

			int r, g, b;
			/* XYZ to non - linear sRGB */
			r = (int)(R2 * 255);
			g = (int)(G2 * 255);
			b = (int)(B2 * 255);

			R[i][j] = r;
			G[i][j] = g;
			B[i][j] = b;

		}

	Mat oR(rows, cols, CV_8UC1);
	Mat oG(rows, cols, CV_8UC1);
	Mat oB(rows, cols, CV_8UC1);
	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			oR.at<uchar>(i, j) = R[i][j];;
			oG.at<uchar>(i, j) = G[i][j];;
			oB.at<uchar>(i, j) = B[i][j];;
		}

	Mat o_planes[] = { oB, oG, oR };
	Mat outImage;
	merge(o_planes, 3, outImage);

	namedWindow("output", CV_WINDOW_AUTOSIZE);
	imshow("output", outImage);
	imwrite(outName, outImage);
}

int main(int argc, char** argv) {
	if (argc != 7) {
		cerr << argv[0] << ": "
			<< "got " << argc - 1
			<< " arguments. Expecting six: w1 h1 w2 h2 ImageIn ImageOut."
			<< endl;
		cerr << "Example: proj1b 0.2 0.1 0.8 0.5 fruits.jpg out.bmp" << endl;
		return(-1);
	}
	double w1 = atof(argv[1]);
	double h1 = atof(argv[2]);
	double w2 = atof(argv[3]);
	double h2 = atof(argv[4]);
	char *inputName = argv[5];
	char *outputName = argv[6];

	if (w1<0 || h1<0 || w2 <= w1 || h2 <= h1 || w2>1 || h2>1) {
		cerr << " arguments must satisfy 0 <= w1 < w2 <= 1"
			<< " ,  0 <= h1 < h2 <= 1" << endl;
		return(-1);
	}

	Mat inputImage = imread(inputName, CV_LOAD_IMAGE_UNCHANGED);
	if (inputImage.empty()) {
		cout << "Could not open or find the image " << inputName << endl;
		return(-1);
	}

	string windowInput("input: ");
	windowInput += inputName;

	namedWindow(windowInput, CV_WINDOW_AUTOSIZE);
	imshow(windowInput, inputImage);

	if (inputImage.type() != CV_8UC3) {
		cout << inputName << " is not a standard color image  " << endl;
		return(-1);
	}

	int rows = inputImage.rows;
	int cols = inputImage.cols;
	int W1 = (int)(w1*(cols - 1));
	int H1 = (int)(h1*(rows - 1));
	int W2 = (int)(w2*(cols - 1));
	int H2 = (int)(h2*(rows - 1));

	runOnWindow(W1, H1, W2, H2, inputImage, outputName);

	waitKey(0); // Wait for a keystroke
	return(0);
}
