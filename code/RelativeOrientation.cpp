#include "RelativeOrientation.h"

using namespace std;
using namespace cv;


void RelativeOrientation::read_file(string file, string img1, string img2, vector<Ppoint>& P1, vector<Ppoint>& P2, vector<double>& inter1, vector<double>& inter2, vector<int>& pname)
{
    ifstream infile;
    infile.open(file, ios::in);
    if (!infile)
    {
        cout << "File doesn't exist" << endl;
        return;
    }
    string a = "";
    int i = 0;
    while (!infile.eof())
    {
        getline(infile, a, '\n');
        i += 1;
        istringstream str(a);
        if (i == 1)
        {
            string name[2];
            while (str >> name[0] >> name[1])
            {
                img1 = name[0];
                img2 = name[1];
            }
        }
        else if (i == 2)
        {
            string p[4];
            while (str >> p[0] >> p[1] >> p[2] >> p[3])
            {
                inter1.push_back(atof(p[0].c_str()) / 1000.0);
                inter1.push_back(atof(p[1].c_str()) / 1000.0);
                inter1.push_back(atof(p[2].c_str()) / 1000.0);
                inter1.push_back(atof(p[3].c_str()));
            }
        }
        else if (i == 3)
        {
            string p[4];
            while (str >> p[0] >> p[1] >> p[2] >> p[3])
            {
                inter2.push_back(atof(p[0].c_str()) / 1000.0);
                inter2.push_back(atof(p[1].c_str()) / 1000.0);
                inter2.push_back(atof(p[2].c_str()) / 1000.0);
                inter2.push_back(atof(p[3].c_str()));
            }
        }
        else
        {
            string split[5];
            while (str >> split[0] >> split[1] >> split[2] >> split[3] >> split[4])
            {
                int id = atoi(split[0].c_str());
                double x1 = atof(split[1].c_str()) / 1000.0;
                double y1 = atof(split[2].c_str()) / 1000.0;
                double x2 = atof(split[3].c_str()) / 1000.0;
                double y2 = atof(split[4].c_str()) / 1000.0;
                pname.push_back(id);
                P1.push_back(Ppoint(x1, y1));
                P2.push_back(Ppoint(x2, y2));
            }
        }
    }
    cout << "Model: " << img1 << "-" << img2 << endl;
    cout << "Number of Points: " << P1.size() << endl;
}

void RelativeOrientation::calculate_relarotation_matrix(Mat_<double>& R2, Mat_<double>& Para)
{
    double phi = Para.at<double>(0, 0);
    double omega = Para.at<double>(1, 0);
    double kappa = Para.at<double>(2, 0);
    R2.at<double>(0, 0) = cos(phi) * cos(kappa) - sin(phi) * sin(omega) * sin(kappa);
    R2.at<double>(0, 1) = -cos(phi) * sin(kappa) - sin(phi) * sin(omega) * cos(kappa);
    R2.at<double>(0, 2) = -sin(phi) * cos(omega);
    R2.at<double>(1, 0) = cos(omega) * sin(kappa);
    R2.at<double>(1, 1) = cos(omega) * cos(kappa);
    R2.at<double>(1, 2) = -sin(omega);
    R2.at<double>(2, 0) = sin(phi) * cos(kappa) + cos(phi) * sin(omega) * sin(kappa);
    R2.at<double>(2, 1) = -sin(phi) * sin(kappa) + cos(phi) * sin(omega) * cos(kappa);
    R2.at<double>(2, 2) = cos(phi) * cos(omega);
}

void RelativeOrientation::calculate_A_matrix(int i, Mat_<double>& A, SPoint& P2, double N1, double N2, double Bx)
{
    double a[5] = {};
    a[0] = -P2.x * P2.y * N2 / P2.z; //¦Õ
    a[1] = -(P2.z + P2.y * P2.y / P2.z) * N2; //¦Ø
    a[2] = P2.x * N2; //¦Ê
    a[3] = Bx; //u
    a[4] = -P2.y * Bx / P2.z; //v
    for (int j = 0; j < 5; j++)
    {
        A.at<double>(i, j) = a[j];
    }
}

void RelativeOrientation::calculate_L_matrix(int i, Mat_<double>& L, double Q)
{
    L.at<double>(i, 0) = Q;
}

void RelativeOrientation::update(Mat_<double>& Para, Mat_<double>& X)
{
    for (int i = 0; i < 5; i++)
    {
        Para.at<double>(i, 0) += X.at<double>(i, 0);
    }
}

void RelativeOrientation::relative_orientation(string infile4, string outfile4)
{
    //Define and initiate parameter
    int iteration = 0;
    string img1, img2;
    vector<Ppoint> P1, P2;
    vector<double> inter1, inter2;
    vector<int> pname;
    RelativeOrientation::read_file(infile4, img1, img2, P1, P2, inter1, inter2, pname);
    double f, f_;
    double x0, x0_, y0, y0_;
    double m1, m2;
    double pixel = 0.021/1000.0;
    f = inter1[0];
    f_ = inter2[0];
    x0 = inter1[1];
    x0_ = inter2[1];
    y0 = inter1[2];
    y0_ = inter2[2];
    m1 = inter1[3];
    m2 = inter2[3];
    //Initiate
    int point_num = P1.size();
    Mat_<double> Para = Mat::zeros(5, 1, CV_32F);
    Mat_<double> X = Mat::zeros(5, 1, CV_32F);
    Mat_<double> A = Mat::zeros(point_num, 5, CV_32F);
    Mat_<double> L = Mat::zeros(point_num, 1, CV_32F);
    vector<SPoint> Aux1(point_num);
    vector<SPoint> Aux2(point_num);
    vector<SPoint> model(point_num);
    double Bx = P1[0].x - P2[0].x;
    double By = 0, Bz = 0;

    //Calculate parameter with iteration
    while (true)
    {
        iteration += 1;
        model.clear();
        Mat_<double> R2 = Mat::zeros(3, 3, CV_32F);
        RelativeOrientation::calculate_relarotation_matrix(R2, Para);
        By = Bx * Para.at<double>(3, 0);
        Bz = Bx * Para.at<double>(4, 0);
        //Calculate model point
        for (int i = 0; i < point_num; i++)
        {
            //Calculate Auxiliary coordinate
            Aux1[i].x = P1[i].x;
            Aux1[i].y = P1[i].y;
            Aux1[i].z = -f;
            Mat_<double> b = Mat::zeros(3, 1, CV_32F);
            b.at<double>(0, 0) = P2[i].x;
            b.at<double>(1, 0) = P2[i].y;
            b.at<double>(2, 0) = -f_;
            Mat_<double> aux2 = R2 * b;
            Aux2[i].x = aux2.at<double>(0, 0);
            Aux2[i].y = aux2.at<double>(1, 0);
            Aux2[i].z = aux2.at<double>(2, 0);
            //Calculate N1, N2
            double N1 = 0, N2 = 0;
            N1 = (Bx * Aux2[i].z - Bz * Aux2[i].x) / (Aux1[i].x * Aux2[i].z - Aux2[i].x * Aux1[i].z);
            N2 = (Bx * Aux1[i].z - Bz * Aux1[i].x) / (Aux1[i].x * Aux2[i].z - Aux2[i].x * Aux1[i].z);
            //Calculate Q of each point
            double Q = N1 * Aux1[i].y - N2 * Aux2[i].y - By;
            //Calculate A, L
            RelativeOrientation::calculate_A_matrix(i, A, Aux2[i], N1, N2, Bx);
            RelativeOrientation::calculate_L_matrix(i, L, Q);
            //Calculate model coordinates
            double mx = m1 * N1 * Aux1[i].x;
            double my = 0.5 * m1 * (N1 * Aux1[i].y + N2 * Aux2[i].y + By);
            double mz = m1 * f + m1 * N1 * Aux1[i].z;
            model.push_back(SPoint(mx, my, mz));
        }
        //Update value of Parameters
        X = (A.t() * A).inv() * A.t() * L;
        RelativeOrientation::update(Para, X);

        //Ouput result if convergent
        if (abs(X.at<double>(0, 0)) < 0.00003 && abs(X.at<double>(1, 0)) < 0.00003 && abs(X.at<double>(2, 0)) < 0.00003 && abs(X.at<double>(3, 0)) < 0.00003 && abs(X.at<double>(4, 0)) < 0.00003)
        {
            cout << "--------------------------------------------" << endl;
            cout << "Relative Orientation Result " << endl;
            cout << "Iteration: " << iteration << endl;
            cout << "Residual(pixel): " << endl;
            Mat_<double>V = A * X - L;
            cout << V/pixel << endl;
            cout << "Five Parameters of Relative Orientation(phi, omega, kappa, u, v): " << endl;
            cout << Para.at<double>(0, 0) << " " << Para.at<double>(1, 0) << " " << Para.at<double>(2, 0) << " " << Para.at<double>(3, 0) << "  " << Para.at<double>(4, 0) << endl;
            Mat_<double>V_ = V.t() * V;
            double accuracy = sqrt(V_.at<double>(0, 0) / (point_num - 5));
            cout << "Root Mean Square Error£º" << accuracy << endl;
            cout << "Model Coordinates(id x y z):" << endl;
            for (int i = 0; i < point_num; i++)
            {
                cout << pname[i] << " " << model[i].x << " " << model[i].y << " " << model[i].z << endl;
            }
            ofstream outfile;
            outfile.open(outfile4, ios::out);
            outfile << "Number of Point: " << point_num << endl;
            outfile << "Iteration: " << iteration << endl;
            outfile << "Residual: " << endl;
            outfile << V << endl;
            outfile << "Five Parameters of Relative Orientation(phi, omega, kappa, u, v): " << endl;
            outfile << Para.at<double>(0, 0) << " " << Para.at<double>(1, 0) << " " << Para.at<double>(2, 0) << " " << Para.at<double>(3, 0) << "  " << Para.at<double>(4, 0) << endl;
            outfile << "Root Mean Square Error£º" << accuracy << endl;
            outfile << "Model coordinates(id x y z):" << endl;
            for (int i = 0; i < point_num; i++)
            {
                outfile << pname[i] << " " << model[i].x << " " << model[i].y << " " << model[i].z << endl;
            }
            break;
        }
    }
}