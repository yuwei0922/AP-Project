#include "AbsoluteOrientation.h"

using namespace std;
using namespace cv;


void AbsoluteOrientation::read_file(string file, vector<string>& pname, vector<SPoint>& Pmodel, vector<SPoint>& Pspace)
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
        string split[7];
        while (str >> split[0] >> split[1] >> split[2] >> split[3] >> split[4] >> split[5] >> split[6])
        {
            string name = split[0];
            double xmodel = atof(split[1].c_str());
            double ymodel = atof(split[2].c_str());
            double zmodel = atof(split[3].c_str());
            double xspace = atof(split[4].c_str());
            double yspace = atof(split[5].c_str());
            double zspace = atof(split[6].c_str());
            pname.push_back(name);
            Pmodel.push_back(SPoint(xmodel, ymodel, zmodel));
            Pspace.push_back(SPoint(xspace, yspace, zspace));
        }
    }
}

void AbsoluteOrientation::calculate_rotation_matrix(Mat_<double>& R, Mat_<double>& Para)
{
    double phi = Para.at<double>(3, 0);
    double omega = Para.at<double>(4, 0);
    double kappa = Para.at<double>(5, 0);
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

void AbsoluteOrientation::calculate_A_matrix(int i, Mat_<double>& A, vector<SPoint>& Pmodel, Mat_<double>& Para)
{
    double phi = Para.at<double>(3, 0);
    double omega = Para.at<double>(4, 0);
    double kappa = Para.at<double>(5, 0);
    double s = Para.at<double>(6, 0);
   // Mat_<double> A = Mat::zeros(3, 7, CV_32F);
    A.at<double>(3 * i, 0) = 1; 
    A.at<double>(3 * i, 1) = 0;
    A.at<double>(3 * i, 2) = 0;
    A.at<double>(3 * i, 3) = Pmodel[i].x;
    A.at<double>(3 * i, 4) = -s * Pmodel[i].z; //phi
    A.at<double>(3 * i, 5) = -s * Pmodel[i].y * sin(phi); //omega
    A.at<double>(3 * i, 6) = -s * Pmodel[i].y * cos(phi) * cos(omega) - s * Pmodel[i].z * sin(omega); //kappa;
    A.at<double>(3 * i + 1, 0) = 0;
    A.at<double>(3 * i + 1, 1) = 1;
    A.at<double>(3 * i + 1, 2) = 0;
    A.at<double>(3 * i + 1, 3) = Pmodel[i].y;
    A.at<double>(3 * i + 1, 4) = 0; //phi
    A.at<double>(3 * i + 1, 5) = s * Pmodel[i].x * sin(phi) - s * Pmodel[i].z * cos(phi); //omega
    A.at<double>(3 * i + 1, 6) = s * Pmodel[i].x * cos(phi) * cos(omega) + s * Pmodel[i].z * sin(phi) * cos(omega); //kappa
    A.at<double>(3 * i + 2, 0) = 0;
    A.at<double>(3 * i + 2, 1) = 0;
    A.at<double>(3 * i + 2, 2) = 1;
    A.at<double>(3 * i + 2, 3) = Pmodel[i].z;
    A.at<double>(3 * i + 2, 4) = s * Pmodel[i].x; //phi
    A.at<double>(3 * i + 2, 5) = s * Pmodel[i].y * cos(phi); //omega
    A.at<double>(3 * i + 2, 6) = s * Pmodel[i].x * sin(omega) - s * Pmodel[i].y * sin(phi) * cos(omega); //kappa
}

void AbsoluteOrientation::calculate_L_matrix(int i, Mat_<double>& L, Mat_<double>& mtp, Mat_<double>& mp, Mat_<double>& Para)
{
    Mat_<double> l = Mat::zeros(3, 1, CV_32F);
    double s = Para.at<double>(6, 0);
    Mat_<double> p0 = Mat::zeros(3, 1, CV_32F);
    p0.at<double>(0, 0) = Para.at<double>(0, 0);
    p0.at<double>(1, 0) = Para.at<double>(1, 0);
    p0.at<double>(2, 0) = Para.at<double>(2, 0);
    l = mtp - s * mp - p0;
    L(3 * i, 0) = l.at<double>(0, 0);
    L(3 * i + 1, 0) = l.at<double>(1, 0);
    L(3 * i + 2, 0) = l.at<double>(2, 0);
}

void AbsoluteOrientation::update(Mat_<double>& Para, Mat_<double>& X)
{
    Para.at<double>(0, 0) += X.at<double>(0, 0); //x
    Para.at<double>(1, 0) += X.at<double>(1, 0); //y
    Para.at<double>(2, 0) += X.at<double>(2, 0); //z
    Para.at<double>(6, 0) *= (1 + X.at<double>(3, 0)); //s
    Para.at<double>(3, 0) += X.at<double>(4, 0); //phi
    Para.at<double>(4, 0) += X.at<double>(5, 0); //omega
    Para.at<double>(5, 0) += X.at<double>(6, 0); //kappa
}

void AbsoluteOrientation::absolute_orientation(string infile5, string outfile5)
{
    //Define and Initiate parameters
    vector<string> Pname;
    vector<SPoint> Pmodel;
    vector<SPoint> Pspace;
    AbsoluteOrientation::read_file(infile5, Pname, Pmodel, Pspace);
    vector<SPoint> Pspace1=Pspace; //For verification
    int point_num = Pmodel.size();
    Mat_<double> Para = Mat::zeros(7, 1, CV_32F); //x, y, z, phi, omega, kappa, s
    Para.at<double>(6, 0) = 1; // initiate scale = 1
    Mat_<double> A = Mat::zeros(3 * point_num, 7, CV_32F);
    Mat_<double> L = Mat::zeros(3 * point_num, 1, CV_32F);
    Mat_<double> X = Mat::zeros(7, 1, CV_32F);
    //Calculate Scale
    double d1 = sqrt(pow((Pmodel[0].x - Pmodel[1].x), 2) + pow((Pmodel[0].y - Pmodel[1].y), 2) + pow((Pmodel[0].z - Pmodel[1].z), 2));
    double d2 = sqrt(pow((Pspace[0].x - Pspace[1].x), 2) + pow((Pspace[0].y - Pspace[1].y), 2) + pow((Pspace[0].z - Pspace[1].z), 2));
    double m = d2 / d1;
    //m = 1;
    for (int i = 0; i < point_num; i++)
    {
        Pspace[i].x /= m;
        Pspace[i].y /= m;
        Pspace[i].z /= m;
    }
    //Calculate barycentric coordinates
    //Calculate center coordinates
    double  xmodel_center = 0;
    double  ymodel_center = 0;
    double  zmodel_center = 0;
    double  xspace_center = 0;
    double  yspace_center = 0;
    double  zspace_center = 0;
    for (int i = 0; i < point_num; i++)
    {
        xmodel_center += Pmodel[i].x;
        ymodel_center += Pmodel[i].y;
        zmodel_center += Pmodel[i].z;
        xspace_center += Pspace[i].x;
        yspace_center += Pspace[i].y;
        zspace_center += Pspace[i].z;
    }
    xmodel_center /= point_num;
    ymodel_center /= point_num;
    zmodel_center /= point_num;
    xspace_center /= point_num;
    yspace_center /= point_num;
    zspace_center /= point_num;

    SPoint model_center(xmodel_center, ymodel_center, zmodel_center);
    SPoint space_center(xspace_center, yspace_center, zspace_center);

    //Barycenterization
    for (int i = 0; i < point_num; i++)
    {
        Pmodel[i].x -= model_center.x;
        Pmodel[i].y -= model_center.y;
        Pmodel[i].z -= model_center.z;
        Pspace[i].x -= space_center.x;
        Pspace[i].y -= space_center.y;
        Pspace[i].z -= space_center.z;
    }

    //Calculate parameters with iteration
    int iteration = 0;
    while (true)
    {
        iteration += 1;
        Mat_<double> R = Mat::zeros(3, 3, CV_32F);
        AbsoluteOrientation::calculate_rotation_matrix(R, Para);
        //Calculate coefficient of each point
        for (int i = 0; i < point_num; i++)
        {
            Mat_<double> mtp = Mat::zeros(3, 1, CV_32F);
            Mat_<double> mp = Mat::zeros(3, 1, CV_32F);
            Mat_<double> t = Mat::zeros(3, 1, CV_32F);
            mtp.at<double>(0, 0) = Pspace[i].x;
            mtp.at<double>(1, 0) = Pspace[i].y;
            mtp.at<double>(2, 0) = Pspace[i].z;
            t.at<double>(0, 0) = Pmodel[i].x;
            t.at<double>(1, 0) = Pmodel[i].y;
            t.at<double>(2, 0) = Pmodel[i].z;
            mp = R * t;
            //Calculate A, L to form norm equation
            AbsoluteOrientation::calculate_A_matrix(i, A, Pmodel, Para);
            AbsoluteOrientation::calculate_L_matrix(i, L, mtp, mp, Para);
        }
        //Calculate corrections of parameters
        X = (A.t() * A).inv() * A.t() * L;
        //Update value of parameters
        AbsoluteOrientation::update(Para, X);

        //Ouput result if convergent
        if (abs(X(0, 0)) < 0.00003 && abs(X(1, 0)) < 0.00003 && abs(X(2, 0)) < 0.00003 && abs(X(3, 0)) < 0.00003 && abs(X(4, 0)) < 0.00003 && abs(X(5, 0)) < 0.00003 && abs(X(6, 0)) < 0.00003)
        {

            cout << "--------------------------------------------" << endl;
            cout << "Absolute Orientation Result: " << endl;
            cout << "Iteration: " << iteration << endl;
            cout << "Residual:" << endl;
            Mat_<double>V = A * X - L;
            cout << V << endl;
            cout << "Seven Parameters of Relative Orientation(x y z phi omega kappa s): " << endl;
            cout << Para.at<double>(0, 0) + space_center.x << " " << Para.at<double>(1, 0) + space_center.y << " " << Para.at<double>(2, 0) + space_center.z << " " << Para.at<double>(3, 0)
                << "  " << Para.at<double>(4, 0) << " " << Para.at<double>(5, 0) << " " << Para.at<double>(6, 0) * m << endl;
            Mat_<double>V_ = V.t() * V;
            double accuracy = sqrt(V_.at<double>(0, 0) / (3 * point_num - 7));
            cout << endl;
            cout << "Root Mean Square Error£º" << accuracy << endl;
            cout << endl;
            ofstream outfile;
            outfile.open(outfile5, ios::out);
            outfile << "--------------------------------------------" << endl;
            outfile << "Absolute Orientation Result: " << endl;
            outfile << "Iteration: " << iteration << endl;
            outfile << "Residual:" << endl;
            outfile << V << endl;
            outfile << "Seven Parameters of Relative Orientation(x y z phi omega kappa s): " << endl;
            outfile << Para.at<double>(0, 0) + space_center.x << " " << Para.at<double>(1, 0) + space_center.y << " " << Para.at<double>(2, 0) + space_center.z << " " << Para.at<double>(3, 0)
                << "  " << Para.at<double>(4, 0) << " " << Para.at<double>(5, 0) << " " << Para.at<double>(6, 0) * m << endl;
            outfile << endl;
            outfile << "Root Mean Square Error£º" << accuracy << endl;
            outfile << endl;

            Mat_<double> Rfinal = Mat::zeros(3, 3, CV_32F);
            Mat_<double> S = Mat::zeros(3, 1, CV_32F);
            Mat_<double> result = Mat::zeros(3, 1, CV_32F);
            AbsoluteOrientation::calculate_rotation_matrix(Rfinal, Para);
            S.at<double>(0, 0) = space_center.x * m;
            S.at<double>(1, 0) = space_center.y * m;
            S.at<double>(2, 0) = space_center.z * m;
            cout << "Ground Coordinate and Residual: " << endl;
            outfile << "Ground Coordinate and Residual: " << endl;
            cout << "name X              Y                Z            dX          dY          dZ " << endl;
            outfile << "name X              Y                Z            dX          dY          dZ " << endl;
            for (int i = 0; i < point_num; i++)
            {
                Mat_<double> m0 = Mat::zeros(3, 1, CV_32F);
                m0.at<double>(0, 0) = Pmodel[i].x;
                m0.at<double>(1, 0) = Pmodel[i].y;
                m0.at<double>(2, 0) = Pmodel[i].z;
                result = Para.at<double>(6, 0) * m * Rfinal * m0 + S;
                cout << Pname[i] << fixed << setprecision(6) << "   " << result.at<double>(0, 0) << "   " << result.at<double>(1, 0) << "   " << result.at<double>(2, 0) << "   " << Pspace1[i].x- result.at<double>(0, 0) << "   " << Pspace1[i].y - result.at<double>(1, 0) << "   " << Pspace1[i].z - result.at<double>(2, 0) << endl;
                outfile << Pname[i] << fixed << setprecision(6) << "   " << result.at<double>(0, 0) << "   " << result.at<double>(1, 0) << "   " << result.at<double>(2, 0) << "   " << Pspace1[i].x - result.at<double>(0, 0) << "   " << Pspace1[i].y - result.at<double>(1, 0) << "   " << Pspace1[i].z - result.at<double>(2, 0) << endl;
            }
            break;
        }
    }
}