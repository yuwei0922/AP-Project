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


class AbsoluteOrientation
{
public:
	static void read_file(string file, vector<string>& pname, vector<SPoint>& Pmodel, vector<SPoint>& Pspace);
	static void calculate_rotation_matrix(Mat_<double>& R, Mat_<double>& Para);
	static void calculate_A_matrix(int i, Mat_<double>& A, vector<SPoint>& Pmodel, Mat_<double>& Para);
	static void calculate_L_matrix(int i, Mat_<double>& L, Mat_<double>& mtp, Mat_<double>& mp, Mat_<double>& Para);
	static void update(Mat_<double>& Para, Mat_<double>& X);
	static void absolute_orientation(string infile5, string outfile5);
};

