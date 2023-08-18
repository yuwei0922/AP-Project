#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <string>
#include "SPoint.h"
#include "Ppoint.h"

using namespace std;
using namespace cv;


class SpaceIntersection
{
public:
	static void read_file(string file1, string file2, vector<double>& inter1, vector<double>& exter1,
		vector<Ppoint>& P1, vector<double>& inter2, vector<double>& exter2, vector<Ppoint>& P2, vector<int>& ids);
	static void calculate_rotation_matrix(Mat_<double>& R, const double phi, const double omega, const double kappa);
	static Mat coordinate_change(Mat_<double>& R, double x, double y, double f);
	static SPoint point_projection(vector<double> exter1, vector<double> exter2, double N1, double N2, Mat_<double>& Aux, Mat_<double>& Aux_);
	static void calculate_A_matrix(Mat_<double>& A, Mat_<double>& R, Mat_<double>& R_, vector<Ppoint>& P1, vector<Ppoint>& P2,
		double Z1, double Z2, int i, double f, double f_);
	static void calculate_L_matrix(Mat_<double>& L, const vector<Ppoint>& P1, const vector<Ppoint>& P2, const vector<Ppoint>& Appro, int i);
	static void update(double& Xt, double& Yt, double& Zt, Mat_<double>& X);
	static bool ifconvergent(Mat_<double>& X, double limit);
	static void collinearity_equation(vector<double> inter1, vector<double> inter2, vector<double> exter1, vector<double> exter2,
		vector<Ppoint> P1, vector<Ppoint> P2, vector<SPoint> S1, vector<SPoint>& S2, double& accuracy);
	static void space_intersection(string file1, string file2, string outfile);
};
