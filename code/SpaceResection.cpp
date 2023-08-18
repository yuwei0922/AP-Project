#include "SpaceResection.h"

using namespace std;
using namespace cv;

void SpaceResection::read_file(string infile1, vector<SPoint>& S, vector<Ppoint>& P)
{
    ifstream file;
    file.open(infile1, ios::in);
    if (!file)
    {
        cout << "File doesn't exist" << endl;
        return;
    }

    string a = "";
    while (!file.eof())
    {
        getline(file, a, '\n');
        istringstream str(a);
        string split[6];
        while (str >> split[0] >> split[1] >> split[2] >> split[3] >> split[4] >> split[5])
        {
            int id = atoi(split[0].c_str());
            double ix = atof(split[1].c_str()) / 1000.0; //Unit Conversion mm->m
            double iy = atof(split[2].c_str()) / 1000.0;
            double gx = atof(split[3].c_str());
            double gy = atof(split[4].c_str());
            double gz = atof(split[5].c_str());
            P.push_back(Ppoint(ix, iy));
            S.push_back(SPoint(gx, gy, gz));
        }
    }
    cout << "--------------------------------------------" << endl;
    cout <<"Number of Points: " << P.size() << endl;
    cout << endl;
}

void SpaceResection::init(double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa, const vector<SPoint>& S, double m, double f)
{
    double Z = 0;
    int point_num = S.size();
    //Iterate through container
    for (vector<SPoint>::const_iterator it = S.cbegin(); it != S.cend(); it++) 
    {
        Xs += (*it).x;
        Ys += (*it).y;
        Z += (*it).z;
    }
    Xs /= point_num;
    Ys /= point_num;
    Zs /= m * f + Z / point_num;
    phi = 0;
    omega = 0;
    kappa = 0;
}

void SpaceResection::calculate_rotation_matrix(Mat_<double>& R, const double phi, const double omega, const double kappa)
{
    R.at<double>(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
    R.at<double>(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
    R.at<double>(0, 2) = -sin(phi) * cos(omega);
    R.at<double>(1, 0) = cos(omega) * sin(kappa);
    R.at<double>(1, 1) = cos(omega) * cos(kappa);
    R.at<double>(1, 2) = -sin(omega);
    R.at<double>(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
    R.at<double>(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
    R.at<double>(2, 2) = cos(phi) * cos(omega);
}

void SpaceResection::calculate_A_matrix(Mat_<double>& A, Mat_<double>& R, vector<Ppoint>& P, double Z, int i, double phi, double omega, double kappa, double f)
{
    A.at<double>(i * 2, 0) = (R.at<double>(0, 0) * f + R.at<double>(0, 2) * P[i].x) / Z; //a11
    A.at<double>(i * 2, 1) = (R.at<double>(1, 0) * f + R.at<double>(1, 2) * P[i].x) / Z; //a12
    A.at<double>(i * 2, 2) = (R.at<double>(2, 0) * f + R.at<double>(2, 2) * P[i].x) / Z; //a13
    A.at<double>(i * 2, 3) = P[i].y * sin(omega) - (P[i].x * (P[i].x * cos(kappa) - P[i].y * sin(kappa)) / f + f * cos(kappa)) * cos(omega); //a14
    A.at<double>(i * 2, 4) = -f * sin(kappa) - P[i].x * (P[i].x * sin(kappa) + P[i].y * cos(kappa)) / f; //a15
    A.at<double>(i * 2, 5) = P[i].y; //a16
    A.at<double>(i * 2 + 1, 0) = (R.at<double>(0, 1) * f + R.at<double>(0, 2) * P[i].y) / Z; //a21
    A.at<double>(i * 2 + 1, 1) = (R.at<double>(1, 1) * f + R.at<double>(1, 2) * P[i].y) / Z; //a22
    A.at<double>(i * 2 + 1, 2) = (R.at<double>(2, 1) * f + R.at<double>(2, 2) * P[i].y) / Z; //a23
    A.at<double>(i * 2 + 1, 3) = -P[i].x * sin(omega) - (P[i].y * (P[i].x * cos(kappa) - P[i].y * sin(kappa)) / f - f * sin(kappa)) * cos(omega); //a24
    A.at<double>(i * 2 + 1, 4) = -f * cos(kappa) - (P[i].y * (P[i].x * sin(kappa) + P[i].y * cos(kappa))) / f; //a25
    A.at<double>(i * 2 + 1, 5) = -P[i].x; //a26
}

void SpaceResection::calculate_L_matrix(Mat_<double>& L, const vector<Ppoint>& P, const vector<Ppoint>& Appro, int i)
{
    L.at<double>(i * 2, 0) = P[i].x - Appro[i].x;
    L.at<double>(i * 2 + 1, 0) = P[i].y - Appro[i].y;
}

void SpaceResection::update(double& Xs, double& Ys, double& Zs, double& phi, double& omega, double& kappa, Mat_<double>& X)
{
    Xs += X.at<double>(0, 0);
    Ys += X.at<double>(1, 0);
    Zs += X.at<double>(2, 0);
    phi += X.at<double>(3, 0);
    omega += X.at<double>(4, 0);
    kappa += X.at<double>(5, 0);
}

bool SpaceResection::ifconvergent(Mat_<double>& X, double limit)
{
    if (fabs(X.at<double>(3, 0)) < limit && fabs(X.at<double>(4, 0)) < limit && fabs(X.at<double>(5, 0)) < limit)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void SpaceResection::space_resection(string infile1, string outfile1)
{
    
    //Get known data and image coordinates of points
    double f = (153.24 / 1000.0); //Unit Conversion mm->m
    double m = 15000;
    Mat_<double> R = Mat::zeros(3, 3, CV_32F); //Rotation matrix
    Mat_<double> X = Mat::zeros(6, 1, CV_32F); //Solution of error equation
    vector<Ppoint> P;
    vector<SPoint> S;
    SpaceResection::read_file(infile1, S, P);

    //Define variable
    double limit = 0.00003; //Correction of angle element: 0.1'=3.0e-5
    int iteration = 0;
    int point_num = P.size();
    vector<Ppoint> Appro(point_num);
    Mat_<double> A = Mat::zeros(point_num * 2, 6, CV_32F);
    Mat_<double> L = Mat::zeros(point_num * 2, 1, CV_32F);

    //Initialize parameters
    double Xs = 0, Ys = 0, Zs = 0;
    double phi = 0, omega = 0, kappa = 0;
    init(Xs, Ys, Zs, phi, omega, kappa, S, m, f);

    //Update parameters with iterations
    do
    {
        //Calculate rotation matrix
        calculate_rotation_matrix(R, phi, omega, kappa);
        //Calculate approximation of x&y of each point
        //Collinearity Condition Equations
        for (int i = 0; i < point_num; i++)
        {
            //Approximation
            double Z = R.at<double>(0, 2) * (S[i].x - Xs) + R.at<double>(1, 2) * (S[i].y - Ys) + R.at<double>(2, 2) * (S[i].z - Zs);
            Appro[i].x = -f * (R.at<double>(0, 0) * (S[i].x - Xs) + R.at<double>(1, 0) * (S[i].y - Ys) + R.at<double>(2, 0) * (S[i].z - Zs)) / Z;
            Appro[i].y = -f * (R.at<double>(0, 1) * (S[i].x - Xs) + R.at<double>(1, 1) * (S[i].y - Ys) + R.at<double>(2, 1) * (S[i].z - Zs)) / Z;
            //Calculate matrix A and L
            calculate_A_matrix(A, R, P, Z, i, phi, omega, kappa, f);
            calculate_L_matrix(L, P, Appro, i);
        }
        //Calculate corrections of exterior elements
        X = (A.t() * A).inv() * A.t() * L;
        //Update exterior elements
        update(Xs, Ys, Zs, phi, omega, kappa, X);
        iteration += 1;
    } while (!ifconvergent(X, limit));

    //Calculate the mean square error of unit weight from V
    Mat_<double> V = A * X - L;
    Mat_<double> V2 = V.t() * V;
    double accuracy = sqrt(V2.at<double>(0, 0) / (point_num * 2 - 6));

    //Output final result
    cout << "Space Resection Result" << endl;
    cout << "Iterations: " << iteration << endl;
    cout << "Xs = " << fixed << setprecision(2) << Xs << endl;
    cout << "Ys = " << fixed << setprecision(2) << Ys << endl;
    cout << "Zs = " << fixed << setprecision(2) << Zs << endl;
    cout << "phi = " << fixed << setprecision(5) << phi << endl;
    cout << "omega = " << fixed << setprecision(5) << omega << endl;
    cout << "kappa = " << fixed << setprecision(5) << kappa << endl;
    cout << "Rotation Matrix:" << endl;
    cout << R << endl;
    cout << "Root Mean Square Error: " << fixed << setprecision(8) << accuracy << endl;
    cout << endl;

    ofstream outfile;
    outfile.open(outfile1, ios::out);
    outfile << "Xs = " << fixed << setprecision(2) << Xs << endl;
    outfile << "Xs = " << fixed << setprecision(2) << Xs << endl;
    outfile << "Ys = " << fixed << setprecision(2) << Ys << endl;
    outfile << "Zs = " << fixed << setprecision(2) << Zs << endl;
    outfile << "phi = " << fixed << setprecision(5) << phi << endl;
    outfile << "omega = " << fixed << setprecision(5) << omega << endl;
    outfile << "kappa = " << fixed << setprecision(5) << kappa << endl;
    outfile << "Rotation Matrix:" << endl;
    outfile << fixed << setprecision(5) << R << endl;
    outfile << "Root Mean Square Error: " << fixed << setprecision(8) << accuracy << endl;
    outfile << endl;
}