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
	double** x = new double*[rows];
	double** y = new double*[rows];
	double** capY = new double*[rows];
	int** R = new int*[rows];
	int** G = new int*[rows];
	int** B = new int*[rows];

	for (int i = 0; i < rows; i++) {
		x[i] = new double[cols];
		y[i] = new double[cols];
		capY[i] = new double[cols];
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

			x[i][j] = X / (X + Y + Z);
			y[i][j] = Y / (X + Y + Z);
			capY[i][j] = Y;

		}

	double Ymin = 1000;
	double Ymax = -1;

	for (int i = H1; i <= H2; i++)
		for (int j = W1; j <= W2; j++) {
			if (capY[i][j] < Ymin) {
				Ymin = capY[i][j];
			}
			else if (capY[i][j] > Ymax) {
				Ymax = capY[i][j];
			}
		}

	for (int i = 0; i < rows; i++)
		for (int j = 0; j < cols; j++) {
			double newY = capY[i][j];
			if (capY[i][j] <= Ymin) {
				newY = 0;
			}
			else if (capY[i][j] >= Ymax) {
				newY = 1;
			}
			else {
				newY = (capY[i][j] - Ymin) / (Ymax - Ymin) * 1.0;
			}

			double X = x[i][j] / y[i][j] * newY;
			double Z = (1 - x[i][j] - y[i][j]) / y[i][j] * newY;

			/* XYZ to linear sRGB */
			double R2 = 3.240479 * X - 1.53715 * newY - 0.498535 * Z;
			double G2 = -0.969256 * X + 1.875991 * newY + 0.041556 * Z;
			double B2 = 0.055648 * X - 0.204043 * newY + 1.057311 * Z;

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
