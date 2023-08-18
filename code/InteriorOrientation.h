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


class InteriorOrientation
{
public:
	static void read_file(string file, vector<Ppoint>& P1, vector<Ppoint>& P2);
	static void calculate_A_matrix(Mat_<double>& A, vector<Ppoint>& Pi, double pixel);
	static void calculate_L_matrix(Mat_<double>& L, vector<Ppoint>& Pf);
	static void interior_orientation(string infile3, string outfile3);
};