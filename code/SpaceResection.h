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

class SpaceResection
{
public:
	static void read_file(string infile1, vector<SPoint>& G, vector<Ppoint>& P);
	static void init(double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa, const vector<SPoint>& G, double m, double f);
	static void calculate_rotation_matrix(Mat_<double>& R, const double phi, const double omega, const double kappa);
	static void calculate_A_matrix(Mat_<double>& A, Mat_<double>& R, vector<Ppoint>& P, double Z, int i, double phi, double omega, double kappa, double f);
	static void calculate_L_matrix(Mat_<double>& L, const vector<Ppoint>& P, const vector<Ppoint>& A, int i);
	static void update(double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa, Mat_<double>& X);
	static bool ifconvergent(Mat_<double>& X, double limit);
	static void space_resection(string infile1, string outfile1);
};