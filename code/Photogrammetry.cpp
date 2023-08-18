#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <opencv2/opencv.hpp>
#include <string>
#include "SPoint.h"
#include "Ppoint.h"
#include "SpaceResection.h"
#include "SpaceIntersection.h"
#include "InteriorOrientation.h"
#include "RelativeOrientation.h"
#include "AbsoluteOrientation.h"

using namespace std;
using namespace cv;

int main()
{
    int task = 0;
    cout << "Here are five tasks: " << endl;
    cout << "  1(Space Resection)" << endl;
    cout << "  2(Space Intersection)" << endl;
    cout << "  3(Interior Orientation)" << endl;
    cout << "  4(Relative Orientation)" << endl;
    cout << "  5(Absolute Orientation)" << endl;
    cout << "Choose a task to run: " << endl;
    cin >> task;
    //Task1 Space Resection
    if (task == 1)
    {
        string infile1 = "..\\SpaceResection.txt";
        string outfile1 = "..\\SpaceResection_Result.txt";
        SpaceResection::space_resection(infile1, outfile1);
    }

    //Task2 Space Intersection
    else if (task == 2)
    {
        string infile2_1= "..\\SpaceIntersection_319.txt";
        string infile2_2 = "..\\SpaceIntersection_320.txt";
        string outfile2 = "..\\SpaceIntersection_Result.txt";
        SpaceIntersection::space_intersection(infile2_1, infile2_2, outfile2);
    }

    //Task3 Interior Orientation
    else if (task == 3)
    {
        string infile3 = "..\\InteriorOrientation.txt";
        string outfile3 = "..\\InteriorOrientation_Result.txt";
        InteriorOrientation::interior_orientation(infile3, outfile3);
    }

    //Task4 Relative Orientation
    else if (task == 4)
    {
        string infile4 = "..\\RelativeOrientation.txt";
        string outfile4 = "..\\RelativeOrientation_Result.txt";
        RelativeOrientation::relative_orientation(infile4, outfile4);
    }

    //Task5 Absolute Orientation
    else if (task == 5)
    {
        string infile5 = "..\\AbsoluteOrientation.txt";
        string outfile5 = "..\\AbsoluteOrientation_Result.txt";
        AbsoluteOrientation::absolute_orientation(infile5, outfile5);
    }

    return 0;

}