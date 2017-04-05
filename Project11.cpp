#include "opencv2/highgui.hpp"
#include <iostream>
using namespace cv;
using namespace std;

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

int main(int argc, char** argv) {
	if (argc != 3) {
		cout << argv[0] << ": "
			<< "got " << argc - 1 << " arguments. Expecting two: width height."
			<< endl;
		return(-1);
	}

	int width = atoi(argv[1]);
	int height = atoi(argv[2]);
	int** RED1 = new int*[height];
	int** GREEN1 = new int*[height];
	int** BLUE1 = new int*[height];
	int** RED2 = new int*[height];
	int** GREEN2 = new int*[height];
	int** BLUE2 = new int*[height];

	for (int i = 0; i < height; i++) {
		RED1[i] = new int[width];
		GREEN1[i] = new int[width];
		BLUE1[i] = new int[width];
		RED2[i] = new int[width];
		GREEN2[i] = new int[width];
		BLUE2[i] = new int[width];
	}

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++)
		{
			int r1, g1, b1;
			int r2, g2, b2;

			double x = (double)j / (double)width;
			double y = (double)i / (double)height;
			double Y = 1.0;

			double L = 90;
			double u = x * 512 - 255;
			double v = y * 512 - 255;

			/* xyY to XYZ : X = x/y * Y, Y =Y, Z = (1-x-y)/y * Y */
			double X1 = x / y * Y;
			double Z1 = (1 - x - y) / y * Y;

			/* XYZ to linear sRGB */
			double R1 = 3.240479 * X1 - 1.53715 * Y - 0.498535 * Z1;
			double G1 = -0.969256 * X1 + 1.875991 * Y + 0.041556 * Z1;
			double B1 = 0.055648 * X1 - 0.204043 * Y + 1.057311 * Z1;

			if (R1 > 1) { R1 = 1; }
			else if (R1 < 0) { R1 = 0; }
			if (G1 > 1) { G1 = 1; }
			else if (G1 < 0) { G1 = 0; }
			if (B1 > 1) { B1 = 1; }
			else if (B1 < 0) { B1 = 0; }

			/* linear sRGB to non-linear sRGB*/
			r1 = (int)(R1 * 255);
			g1 = (int)(G1 * 255);
			b1 = (int)(B1 * 255);

			/* Luv to XYZ */

			double Xw = 0.95;
			double Yw = 1.0;
			double Zw = 1.09;

			double uw = (4 * Xw) / (Xw + 15 * Yw + 3 * Zw);
			double vw = (9 * Yw) / (Xw + 15 * Yw + 3 * Zw);
			
			double uPrime = (u + 13 * uw * L) / (13 * L);
			double vPrime = (v + 13 * vw * L) / (13 * L);

			double X2;
			double Y2;
			double Z2;

			if (L > 7.9996)
			{
				Y2 = pow(((L + 16) / (int) 116), 3) * Yw;

			}
			else
			{
				Y2 = L / 903.3 * Yw;
			}

			if (vPrime != 0)
			{
				X2 = Y2 * 2.25 * uPrime / vPrime;
				Z2 = Y2 * (3 - 0.75 * uPrime - 5 * vPrime) / vPrime;
			}
			else
			{
				X2 = 0;
				Z2 = 0;
			}

			/* XYZ to linear sRGB */
			double R2 = 3.240479 * X2 - 1.53715 * Y2 - 0.498535 * Z2;
			double G2 = -0.969256 * X2 + 1.875991 * Y2 + 0.041556 * Z2;
			double B2 = 0.055648 * X2 - 0.204043 * Y2 + 1.057311 * Z2;

			if (R2 > 1) { R2 = 1; }
			else if (R2 < 0) { R2 = 0; }
			if (G2 > 1) { G2 = 1; }
			else if (G2 < 0) { G2 = 0; }
			if (B2 > 1) { B2 = 1; }
			else if (B2 < 0) { B2 = 0; }

			R2 = gamma(R2);
			G2 = gamma(G2);
			B2 = gamma(B2);

			/* XYZ to non - linear sRGB */
			r2 = (int)(R2 * 255);
			g2 = (int)(G2 * 255);
			b2 = (int)(B2 * 255);


			// this is the end of your code

			RED1[i][j] = r1;
			GREEN1[i][j] = g1;
			BLUE1[i][j] = b1;
			RED2[i][j] = r2;
			GREEN2[i][j] = g2;
			BLUE2[i][j] = b2;
		}


	Mat R1(height, width, CV_8UC1);
	Mat G1(height, width, CV_8UC1);
	Mat B1(height, width, CV_8UC1);

	Mat R2(height, width, CV_8UC1);
	Mat G2(height, width, CV_8UC1);
	Mat B2(height, width, CV_8UC1);

	for (int i = 0; i < height; i++)
		for (int j = 0; j < width; j++) {

			R1.at<uchar>(i, j) = RED1[i][j];
			G1.at<uchar>(i, j) = GREEN1[i][j];
			B1.at<uchar>(i, j) = BLUE1[i][j];

			R2.at<uchar>(i, j) = RED2[i][j];
			G2.at<uchar>(i, j) = GREEN2[i][j];
			B2.at<uchar>(i, j) = BLUE2[i][j];
		}

	Mat xyY;
	Mat xyY_planes[] = { B1, G1, R1 };
	merge(xyY_planes, 3, xyY);
	namedWindow("xyY", CV_WINDOW_AUTOSIZE);
	imshow("xyY", xyY);

	Mat Luv;
	Mat Luv_planes[] = { B2, G2, R2 };
	merge(Luv_planes, 3, Luv);
	namedWindow("Luv", CV_WINDOW_AUTOSIZE);
	imshow("Luv", Luv);
	waitKey(0); // Wait for a keystroke
	return(0);
}
