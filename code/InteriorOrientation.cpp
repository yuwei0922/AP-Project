#include "InteriorOrientation.h"

using namespace std;
using namespace cv;


void InteriorOrientation::read_file(string file, vector<Ppoint>& Pf, vector<Ppoint>& Pi)
{
    ifstream infile;
    infile.open(file, ios::in);
    if (!infile)
    {
        cout << "File doesn't exist" << endl;
        return;
    }
    string a = "";
    while (!infile.eof())
    {
        getline(infile, a, '\n');
        istringstream str(a);
        string split[4];
        while (str >> split[0] >> split[1] >> split[2] >> split[3])
        {
            double mx = atof(split[0].c_str());
            double my = atof(split[1].c_str());
            double ix = atof(split[2].c_str());
            double iy = atof(split[3].c_str());
            Pf.push_back(Ppoint(mx, my));
            Pi.push_back(Ppoint(ix, iy));
        }
    }
}

void InteriorOrientation::calculate_A_matrix(Mat_<double>& A, vector<Ppoint>& Pi, double pixel)
{
    int point_num = Pi.size();
    for (int i = 0; i < point_num; i++)
    {
        A.at<double>(i * 2, 0) = 1;
        A.at<double>(i * 2, 1) = Pi[i].x * pixel;
        A.at<double>(i * 2, 2) = Pi[i].y * pixel;
        A.at<double>(i * 2, 3) = 0;
        A.at<double>(i * 2, 4) = 0;
        A.at<double>(i * 2, 5) = 0;

        A.at<double>(i * 2 + 1, 0) = 0;
        A.at<double>(i * 2 + 1, 1) = 0;
        A.at<double>(i * 2 + 1, 2) = 0;
        A.at<double>(i * 2 + 1, 3) = 1;
        A.at<double>(i * 2 + 1, 4) = Pi[i].x * pixel;
        A.at<double>(i * 2 + 1, 5) = Pi[i].y * pixel;
    }
}

void InteriorOrientation::calculate_L_matrix(Mat_<double>& L, vector<Ppoint>& Pf)
{
    int point_num = Pf.size();
    for (int i = 0; i < point_num; i++)
    {
        L.at<double>(i * 2, 0) = Pf[i].x;
        L.at<double>(i * 2 + 1, 0) = Pf[i].y;
    }
}

void InteriorOrientation::interior_orientation(string infile3, string outfile3)
{
    //Define and Initiate parameters
    vector<Ppoint> Pf;
    vector<Ppoint> Pi;
    InteriorOrientation::read_file(infile3, Pf, Pi);
    double pixel = 0.021;
    double x0 = 0.011;
    double y0 = 0.002;
    int point_num = Pf.size();
    vector<Ppoint> result;
    Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
    Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);
    Mat_<double> X = Mat::zeros(6, 1, CV_32F);

    //Form normal equation and calculate accuracy
    calculate_A_matrix(A, Pi, pixel);
    calculate_L_matrix(L, Pf);
    X = (A.t() * A).inv() * A.t() * L;
    Mat_<double> V = A * X - L;
    Mat_<double> V2 = V.t() * V;
    double accuracy = sqrt(V2.at<double>(0, 0) / (point_num * 2 - 6));

    double m0 = X.at<double>(0, 0);
    double m1 = X.at<double>(1, 0);
    double m2 = X.at<double>(2, 0);
    double n0 = X.at<double>(3, 0);
    double n1 = X.at<double>(4, 0);
    double n2 = X.at<double>(5, 0);

    //Output result
    cout << "--------------------------------------------" << endl;
    cout << "Interior Orientation Result" << endl;
    cout << "m0 = " << fixed << setprecision(5) << m0 << endl;
    cout << "m1 = " << fixed << setprecision(5) << m1 << endl;
    cout << "m2 = " << fixed << setprecision(5) << m2 << endl;
    cout << "n0 = " << fixed << setprecision(5) << n0 << endl;
    cout << "n1 = " << fixed << setprecision(5) << n1 << endl;
    cout << "n2 = " << fixed << setprecision(5) << n2 << endl;
    cout << "Root Mean Square Error: " << fixed << setprecision(5) << accuracy << endl;
    ofstream outfile;
    outfile.open(outfile3, ios::out);
    outfile << "m0 = " << fixed << setprecision(5) << m0 << endl;
    outfile << "m1 = " << fixed << setprecision(5) << m1 << endl;
    outfile << "m2 = " << fixed << setprecision(5) << m2 << endl;
    outfile << "n0 = " << fixed << setprecision(5) << n0 << endl;
    outfile << "n1 = " << fixed << setprecision(5) << n1 << endl;
    outfile << "n2 = " << fixed << setprecision(5) << n2 << endl;
    outfile << "Root Mean Square Error: " << fixed << setprecision(5) << accuracy << endl;

    cout << "--------------------------------------------" << endl;
    cout << "Input Scanning Coordinates of Point (pixel): " << endl;
    double i, j;
    cin >> i >> j;
    double x = m0 + m1 * i * pixel + m2 * j * pixel - x0;
    double y = n0 + n1 * i * pixel + n2 * j * pixel - y0;
    cout << "Image Coordinates of This Point: " << "(" << x << "," << y << ")" << endl;
}