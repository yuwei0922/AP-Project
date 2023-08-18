#pragma once
#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <string>
#include "SPoint.h"
#include "PPoint.h"

using namespace std;
using namespace cv;


class RelativeOrientation
{
public:
	static void read_file(string file, string img1, string img2, vector<Ppoint>& P1, vector<Ppoint>& P2, vector<double>& inter1, vector<double>& inter2, vector<int>& pname);
	static void calculate_relarotation_matrix(Mat_<double>& R2, Mat_<double>& Para);
	static void calculate_A_matrix(int i, Mat_<double>& A, SPoint& P, double N1, double N2, double Bx);
	static void calculate_L_matrix(int i, Mat_<double>& L, double Q);
	static void update(Mat_<double>& Para, Mat_<double>& X);
	static void relative_orientation(string infile4, string outfile4);
};