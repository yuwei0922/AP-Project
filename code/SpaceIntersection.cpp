#include "SpaceIntersection.h"

using namespace std;
using namespace cv;


void SpaceIntersection::read_file(string file1, string file2, vector<double>& inter1, vector<double>& exter1,
     vector<Ppoint>& P1, vector<double>& inter2, vector<double>& exter2, vector<Ppoint>& P2, vector<int>& ids)
{
    ifstream infile1;
    infile1.open(file1, ios::in);
    if (!infile1)
    {
        cout << "File1 doesn't exist" << endl;
        return;
    }

    //Read the first image: interior & exterior elements and image coordinates of points
    string a = "";
    int i = 0;
    cout << "--------------------------------------------" << endl;
    cout << "Get Parameters and Image Coordinates" << endl <<endl;
    cout << "The First Image: " << endl;
    while (!infile1.eof())
    {
        i++;
        getline(infile1, a, '\n');
        istringstream str(a);
        string interele[4];
        string exterele[6];
        string split[3];
        if (i == 1)
        {
            cout << " Interior Elements: " << a << endl;
            str >> interele[0] >> interele[1] >> interele[2] >> interele[3];
            double f = atof(interele[0].c_str()) / 1000.0;
            double x0 = atof(interele[1].c_str()) / 1000.0;
            double y0 = atof(interele[2].c_str()) / 1000.0;
            double m = atof(interele[3].c_str());
            inter1.push_back(f);
            inter1.push_back(x0);
            inter1.push_back(y0);
            inter1.push_back(m);
        }
        else if (i == 2)
        {
            cout << " Exterior Elements: " << a << endl;
            str >> exterele[0] >> exterele[1] >> exterele[2] >> exterele[3] >> exterele[4] >> exterele[5];
            double Ys = atof(exterele[0].c_str());
            double Xs = atof(exterele[1].c_str());
            double Zs = atof(exterele[2].c_str());
            double phi = atof(exterele[3].c_str());
            double omega = atof(exterele[4].c_str());
            double kappa = atof(exterele[5].c_str());
            phi = phi * CV_PI / 180;
            omega = omega * CV_PI / 180;
            kappa = kappa * CV_PI / 180;
            exter1.push_back(Xs);
            exter1.push_back(Ys);
            exter1.push_back(Zs);
            exter1.push_back(phi);
            exter1.push_back(omega);
            exter1.push_back(kappa);
            cout << " Image Coordinates: " << endl;
        }
        else
        {
            cout << a << endl;
            while (str >> split[0] >> split[1] >> split[2])
            {
                int id = atoi(split[0].c_str());
                double x = atof(split[1].c_str()) / 1000.0;
                double y = atof(split[2].c_str()) / 1000.0;
                ids.push_back(id);
                P1.push_back(Ppoint(x, y));
            }
        }
    }
    cout << " Number of Points: " << P1.size() << endl;

    ifstream infile2;
    infile2.open(file2, ios::in);
    if (!infile2)
    {
        cout << "File2 doesn't exist" << endl;
        return;
    }

    //Read the second image: interior & exterior elements and image coordinates of points
    string b = "";
    int j = 0;
    cout << endl << "The Second Image: " << endl;
    while (!infile2.eof())
    {
        j ++;
        getline(infile2, b, '\n');
        istringstream str(b);
        string interele[4];
        string exterele[6];
        string split[3];
        if (j == 1)
        {
            cout << " Interior Elements: " << b << endl;
            str >> interele[0] >> interele[1] >> interele[2] >> interele[3];
            double f = atof(interele[0].c_str()) / 1000.0;
            double x0 = atof(interele[1].c_str()) / 1000.0;
            double y0 = atof(interele[2].c_str()) / 1000.0;
            double m = atof(interele[3].c_str());
            inter2.push_back(f);
            inter2.push_back(x0);
            inter2.push_back(y0);
            inter2.push_back(m);
        }
        else if (j == 2)
        {
            cout << " Exterior Elements: " << b << endl;
            str >> exterele[0] >> exterele[1] >> exterele[2] >> exterele[3] >> exterele[4] >> exterele[5];
            double Ys = atof(exterele[0].c_str());
            double Xs = atof(exterele[1].c_str());
            double Zs = atof(exterele[2].c_str());
            double phi = atof(exterele[3].c_str());
            double omega = atof(exterele[4].c_str());
            double kappa = atof(exterele[5].c_str());
            phi = phi * CV_PI / 180;
            omega = omega * CV_PI / 180;
            kappa = kappa * CV_PI / 180;
            exter2.push_back(Xs);
            exter2.push_back(Ys);
            exter2.push_back(Zs);
            exter2.push_back(phi);
            exter2.push_back(omega);
            exter2.push_back(kappa);
            cout << " Image Coordinates: " << endl;
        }
        else
        {
            cout << b << endl;
            while (str >> split[0] >> split[1] >> split[2])
            {
                int id = atoi(split[0].c_str());
                double x = atof(split[1].c_str()) / 1000.0;
                double y = atof(split[2].c_str()) / 1000.0;
                P2.push_back(Ppoint(x, y)); //IDs of points have been recorded
            }
        }
    }
    cout << " Number of Points: " << P2.size() << endl;
}

void SpaceIntersection::calculate_rotation_matrix(Mat_<double>& R, const double phi, const double omega, const double kappa)
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

Mat SpaceIntersection::coordinate_change(Mat_<double>& R, double x, double y, double f)
{
    Mat_<double> P(3, 1);
    P.at<double>(0, 0) = x;
    P.at<double>(1, 0) = y;
    P.at<double>(2, 0) = -f;
    Mat_<double> A = R * P;
    return A;
}

SPoint SpaceIntersection::point_projection(vector<double> exter1, vector<double> exter2, double N1, double N2, Mat_<double>& Aux, Mat_<double>& Aux_)
{
    double A1[3] = {};
    double A2[3] = {};
    double AA[3] = {};
    for (int i = 0; i < 3; i++)
    {
        A1[i] = exter1[i] + N1 * Aux.at<double>(i, 0);
        A2[i] = exter2[i] + N2 * Aux_.at<double>(i, 0);
        AA[i] = (A1[i] + A2[i]) / 2; //Average Y ground coordinates
    }
    return SPoint(AA[0], AA[1], AA[2]);
}

void SpaceIntersection::calculate_A_matrix(Mat_<double>& A, Mat_<double>& R, Mat_<double>& R_, vector<Ppoint>& P1, vector<Ppoint>& P2, 
    double Z1, double Z2,int i, double f, double f_)
{
    A.at<double>(i * 4, 0) = - (R.at<double>(0, 0) * f + R.at<double>(0, 2) * P1[i].x) / Z1; //a11
    A.at<double>(i * 4, 1) = - (R.at<double>(1, 0) * f + R.at<double>(1, 2) * P1[i].x) / Z1; //a12
    A.at<double>(i * 4, 2) = - (R.at<double>(2, 0) * f + R.at<double>(2, 2) * P1[i].x) / Z1; //a13
    A.at<double>(i * 4 + 1, 0) = - (R.at<double>(0, 1) * f + R.at<double>(0, 2) * P1[i].y) / Z1; //a21
    A.at<double>(i * 4 + 1, 1) = - (R.at<double>(1, 1) * f + R.at<double>(1, 2) * P1[i].y) / Z1; //a22
    A.at<double>(i * 4 + 1, 2) = - (R.at<double>(2, 1) * f + R.at<double>(2, 2) * P1[i].y) / Z1; //a23
    A.at<double>(i * 4 + 2, 0) = - (R_.at<double>(0, 0) * f_ + R_.at<double>(0, 2) * P2[i].x) / Z2; //a31
    A.at<double>(i * 4 + 2, 1) = - (R_.at<double>(1, 0) * f_+ R_.at<double>(1, 2) * P2[i].x) / Z2; //a32
    A.at<double>(i * 4 + 2, 2) = - (R_.at<double>(2, 0) * f_ + R_.at<double>(2, 2) * P2[i].x) / Z2; //a33
    A.at<double>(i * 4 + 3, 0) = - (R_.at<double>(0, 1) * f_ + R_.at<double>(0, 2) * P2[i].y) / Z2; //a41
    A.at<double>(i * 4 + 3, 1) = - (R_.at<double>(1, 1) * f_ + R_.at<double>(1, 2) * P2[i].y) / Z2; //a42
    A.at<double>(i * 4 + 3, 2) = - (R_.at<double>(2, 1) * f_ + R_.at<double>(2, 2) * P2[i].y) / Z2; //a43
}

void SpaceIntersection::calculate_L_matrix(Mat_<double>& L, const vector<Ppoint>& P1, const vector<Ppoint>& P2, const vector<Ppoint>& Appro, int i)
{
    L.at<double>(i * 4, 0) = P1[i].x - Appro[i].x;
    L.at<double>(i * 4 + 1, 0) = P1[i].y - Appro[i].y;
    L.at<double>(i * 4 + 2, 0) = P2[i].x - Appro[i * 2].x;
    L.at<double>(i * 4 + 3, 0) = P2[i].y - Appro[i * 2].y;
}

void SpaceIntersection::update(double& Xt, double& Yt, double& Zt, Mat_<double>& X)
{
    Xt += X.at<double>(0, 0);
    Yt += X.at<double>(1, 0);
    Zt += X.at<double>(2, 0);
}

bool SpaceIntersection::ifconvergent(Mat_<double>& X, double limit)
{
    if (fabs(X.at<double>(0, 0)) < limit && fabs(X.at<double>(1, 0)) < limit && fabs(X.at<double>(2, 0)) < limit)
    {
        return true;
    }
    else
    {
        return false;
    }
}

void SpaceIntersection::collinearity_equation(vector<double> inter1, vector<double> inter2, vector<double> exter1, vector<double> exter2, 
    vector<Ppoint> P1, vector<Ppoint> P2, vector<SPoint> S1, vector<SPoint>& S2, double& accuracy)
{
    //Define variable
   double limit = 0.000001;
   int point_num = P1.size();
   vector<Ppoint> Appro(point_num *2);
   Mat_<double> R = Mat::zeros(3, 3, CV_32F); //Rotation matrix
   Mat_<double> R_ = Mat::zeros(3, 3, CV_32F); //Rotation matrix
   Mat_<double> A = Mat::zeros(point_num * 4, 6, CV_32F);
   Mat_<double> L = Mat::zeros(point_num * 4, 1, CV_32F);
   Mat_<double> X = Mat::zeros(3, 1, CV_32F); //Solution of error equation

   double f = inter1[0];
   double x0 = inter1[1];
   double y0 = inter1[2];
   double m = inter1[3];
   double Xs = exter1[0];
   double Ys = exter1[1];
   double Zs = exter1[2];
   double phi = exter1[3];
   double omega = exter1[4];
   double kappa = exter1[5];
   double f_ = inter2[0];
   double x0_ = inter2[1];
   double y0_ = inter2[2];
   double m_ = inter2[3];
   double Xs_ = exter2[0];
   double Ys_ = exter2[1];
   double Zs_ = exter2[2];
   double phi_ = exter2[3];
   double omega_ = exter2[4];
   double kappa_ = exter2[5];

   //Calculate rotation matrix
   calculate_rotation_matrix(R, phi, omega, kappa);
   calculate_rotation_matrix(R_, phi_, omega_, kappa_);

   //Update X,Y,Z with iterations
   for (int i = 0; i < point_num; i++)
   {
       do
       {
           //Calculate approximation of x&y of each point
           double Z1 = R.at<double>(0, 2) * (S1[i].x - Xs) + R.at<double>(1, 2) * (S1[i].y - Ys) + R.at<double>(2, 2) * (S1[i].z - Zs);
           Appro[i * 2].x = -f * (R.at<double>(0, 0) * (S1[i].x - Xs) + R.at<double>(1, 0) * (S1[i].y - Ys) + R.at<double>(2, 0) * (S1[i].z - Zs)) / Z1;
           Appro[i * 2].y = -f * (R.at<double>(0, 1) * (S1[i].x - Xs) + R.at<double>(1, 1) * (S1[i].y - Ys) + R.at<double>(2, 1) * (S1[i].z - Zs)) / Z1;
           double Z2 = R_.at<double>(0, 2) * (S1[i].x - Xs_) + R_.at<double>(1, 2) * (S1[i].y - Ys_) + R_.at<double>(2, 2) * (S1[i].z - Zs_);
           Appro[i * 2 + 1].x = -f_ * (R_.at<double>(0, 0) * (S1[i].x - Xs_) + R_.at<double>(1, 0) * (S1[i].y - Ys_) + R_.at<double>(2, 0) * (S1[i].z - Zs_)) / Z2;
           Appro[i * 2 + 1].y = -f_ * (R_.at<double>(0, 1) * (S1[i].x - Xs_) + R_.at<double>(1, 1) * (S1[i].y - Ys_) + R_.at<double>(2, 1) * (S1[i].z - Zs_)) / Z2;
           //Calculate matrix A and L
           calculate_A_matrix(A, R, R_, P1, P2, Z1, Z2, i, f, f_);
           calculate_L_matrix(L, P1, P2, Appro, i);
           //Calculate corrections of X,Y,Z
           X = (A.t() * A).inv() * A.t() * L;
           //Update X,Y,Z
           update(S1[i].x, S1[i].y, S1[i].z, X);
       } while (!ifconvergent(X, limit));
       S2.push_back(SPoint(S1[i].x, S1[i].y, S1[i].z));
       //Calculate the mean square error of unit weight from V
       Mat_<double> V = A * X - L;
       Mat_<double> V2 = V.t() * V;
       double acy= sqrt(V2.at<double>(0, 0));
       accuracy += acy;
   }
   accuracy /= point_num;

}

void SpaceIntersection::space_intersection(string file1, string file2, string outfile2)
{
    //Get known data and image coordinates of points
    vector<double> inter1, inter2;
    vector<double> exter1, exter2;
    vector<Ppoint> P1, P2;
    vector<SPoint> S1;
    vector<SPoint> S2;
    vector<int> ids;
    SpaceIntersection::read_file(file1, file2, inter1, exter1, P1, inter2, exter2, P2, ids);

    //Define and Initiate parameters
    double f = inter1[0];
    double x0 = inter1[1];
    double y0 = inter1[2];
    double m = inter1[3];
    double Xs = exter1[0];
    double Ys = exter1[1];
    double Zs = exter1[2];
    double phi = exter1[3];
    double omega = exter1[4];
    double kappa = exter1[5];
    double f_ = inter2[0];
    double x0_ = inter2[1];
    double y0_ = inter2[2];
    double m_ = inter2[3];
    double Xs_ = exter2[0];
    double Ys_ = exter2[1];
    double Zs_ = exter2[2];
    double phi_ = exter2[3];
    double omega_ = exter2[4];
    double kappa_ = exter2[5];

    //Rotation Matrixs
    Mat_<double> R = Mat::zeros(3, 3, CV_32F);
    Mat_<double> R_ = Mat::zeros(3, 3, CV_32F);

    //Baseline
    double Bx = Xs_ - Xs;
    double By = Ys_ - Ys;
    double Bz = Zs_ - Zs;

    //Calculate rotation matrixs
    calculate_rotation_matrix(R, phi, omega, kappa);
    calculate_rotation_matrix(R_, phi_, omega_, kappa_);

    //Image auxiliary coordinates
    Mat_<double> Aux = Mat::zeros(3, 1, CV_32F);
    Mat_<double> Aux_ = Mat::zeros(3, 1, CV_32F);

    //Spaceintersection by projection coefficient method
    int point_num = P1.size();
    for (int i = 0; i < point_num; i++)
    {
        //Transfer to image auxiliary coordinates
        Aux = coordinate_change(R, P1[i].x, P1[i].y, f);
        Aux_ = coordinate_change(R_, P2[i].x, P2[i].y, f_);
        //Point projection coefficient
        double N1 = (Bx * Aux_.at<double>(2, 0) - Bz * Aux_.at<double>(0, 0)) / (Aux.at<double>(0, 0) * Aux_.at<double>(2, 0) - Aux.at<double>(2, 0) * Aux_.at<double>(0, 0));
        double N2 = (Bx * Aux.at<double>(2, 0) - Bz * Aux.at<double>(0, 0)) / (Aux.at<double>(0, 0) * Aux_.at<double>(2, 0) - Aux.at<double>(2, 0) * Aux_.at<double>(0, 0));
        //Calculate ground coordinates
        SPoint A = point_projection(exter1, exter2, N1, N2, Aux, Aux_);
        S1.push_back(A);
    }

    //Spaceintersection by collinearity equation
    double accuracy = 0.0;
    collinearity_equation(inter1, inter2, exter1, exter2, P1, P2, S1, S2, accuracy);

    //Output final result
    cout << "--------------------------------------------" << endl;
    cout << "Space Intersection Result of Point Projection Coefficient Method " << endl << endl;
    for (int i = 0; i < S1.size(); i++)
    {
        cout << fixed << setprecision(5) << ids[i] << " " << S1[i].y << " " << S1[i].x << " " << S1[i].z << endl;
    }
    cout << endl;


    cout << "--------------------------------------------" << endl;
    cout << "Space Intersection Result of Collinearity Equation " << endl << endl;
    for (int i = 0; i < S2.size(); i++)
    {
        cout << fixed << setprecision(5) << ids[i] << " " << S2[i].y << " " << S2[i].x << " " << S2[i].z << endl;
    }
    cout << endl;
    cout << "Root Mean Square Error: " << fixed << setprecision(8) << accuracy << endl;


    ofstream outfile;
    outfile.open(outfile2, ios::out);
    outfile << "Space Intersection Result of Point Projection Coefficient Method " << endl << endl;
    for (int i = 0; i < S1.size(); i++)
    {
        outfile << fixed << setprecision(5) << ids[i] << " " << S1[i].y << " " << S1[i].x << " " << S1[i].z << endl;
    }


    outfile << "--------------------------------------------" << endl;
    outfile << "Space Intersection Result of Collinearity Equation " << endl << endl;
    for (int i = 0; i < S2.size(); i++)
    {
        outfile << fixed << setprecision(5) << ids[i] << " " << S2[i].y << " " << S2[i].x << " " << S2[i].z << endl;
    }
    outfile << "Root Mean Square Error: " << fixed << setprecision(8) << accuracy << endl;
}