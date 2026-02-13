#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <vector>
#include <iomanip>

using namespace std;

class ErrorNoFile {
    string filename;
public:
    ErrorNoFile(string fn) : filename(fn) {}
    void Message() {
        cout << "File error: " << filename << " not found. Switching to Algorithm 4." << endl;
    }
};

class SwitchToAlg2 {};
class SwitchToAlg3 {};

void createFiles() {
    ofstream f1("dat_X_1_1.dat");
    f1 << "-1.0 -4.935 1.935\n-0.9 -3.013 0.464\n-0.8 -2.316 1.327\n-0.7 -1.819 1.976\n";
    f1 << "-0.6 -1.425 2.502\n-0.5 -1.097 2.951\n-0.4 -0.816 3.344\n-0.3 -0.571 3.695\n";
    f1 << "-0.2 -0.357 4.013\n-0.1 -0.167 4.303\n0.0 0.000 4.571\n0.1 0.147 4.618\n";
    f1 << "0.2 0.276 4.645\n0.3 0.386 4.652\n0.4 0.477 4.636\n0.5 0.548 4.596\n";
    f1 << "0.6 0.597 4.524\n0.7 0.617 4.412\n0.8 0.597 4.240\n0.9 0.505 3.956\n1.0 0.000 3.000\n";
    f1.close();

    ofstream f2("dat_X_1_00.dat");
    f2 << "0.000 -4.935 1.935\n0.050 -2.663 1.885\n0.100 -1.618 1.834\n0.150 -0.773 1.784\n";
    f2 << "0.200 -0.034 1.732\n0.250 0.635 1.679\n0.300 1.253 1.625\n0.350 1.829 1.570\n";
    f2 << "0.400 2.369 1.512\n0.450 2.877 1.452\n0.500 3.356 1.388\n0.550 3.806 1.322\n";
    f2 << "0.600 4.228 1.251\n0.650 4.622 1.175\n0.700 4.987 1.093\n0.750 5.320 1.003\n";
    f2 << "0.800 5.618 0.905\n0.850 5.876 0.796\n0.900 6.080 0.675\n0.950 6.199 0.536\n1.000 5.890 0.377\n";
    f2.close();

    ofstream f3("dat_X00_1.dat");
    f3 << "0.000 -4.935 1.935\n-0.050 -4.435 1.835\n-0.100 -3.936 1.735\n-0.150 -3.440 1.636\n";
    f3 << "-0.200 -2.948 1.537\n-0.250 -2.461 1.440\n-0.300 -1.980 1.344\n-0.350 -1.506 1.249\n";
    f3 << "-0.400 -1.041 1.156\n-0.450 -0.585 1.065\n-0.500 -0.141 0.976\n-0.550 0.292 0.889\n";
    f3 << "-0.600 0.712 0.806\n-0.650 1.117 0.724\n-0.700 1.507 0.646\n-0.750 1.882 0.572\n";
    f3 << "-0.800 2.239 0.500\n-0.850 2.578 0.432\n-0.900 2.898 0.368\n-0.950 3.199 0.308\n-1.000 3.480 0.252\n";
    f3.close();
}

void GetTU(double x, double& T_val, double& U_val) {
    string filename;
    double search_x = x;

    if (abs(x) <= 1) {
        filename = "dat_X_1_1.dat";
    }
    else if (x < -1) {
        search_x = 1.0 / x;
        filename = "dat_X00_1.dat";
    }
    else {
        search_x = 1.0 / x;
        filename = "dat_X_1_00.dat";
    }

    ifstream is(filename);
    if (!is) throw ErrorNoFile(filename);

    double xi, Ti, Ui, xi1, Ti1, Ui1;
    bool found = false;

    is >> xi >> Ti >> Ui;

    while (is >> xi1 >> Ti1 >> Ui1) {
        if ((search_x >= xi && search_x <= xi1) || (search_x <= xi && search_x >= xi1)) {
            T_val = Ti + (Ti1 - Ti) * (search_x - xi) / (xi1 - xi);
            U_val = Ui + (Ui1 - Ui) * (search_x - xi) / (xi1 - xi);
            found = true;
            break;
        }
        xi = xi1; Ti = Ti1; Ui = Ui1;
    }

    if (!found && search_x == xi) {
        T_val = Ti; U_val = Ui; found = true;
    }

    is.close();
    if (!found) { T_val = 0; U_val = 0; }
}

double Srs(double x, double y, double z);
double Srz(double x, double y, double z);
double Srs1(double x, double y, double z);

double Srz(double x, double y, double z) {
    double Tx, Ux, Ty, Uy, Tz, Uz;
    GetTU(x, Tx, Ux);
    GetTU(y, Ty, Uy);
    GetTU(z, Tz, Uz);

    if (x > y) {
        return Tx + Uz - Ty;
    }
    else {
        return Ty + Uy - Uz;
    }
}

double Srs(double x, double y, double z) {
    double val_Srz = Srz(x, y, z);

    if (z > y) {
        if (z * z + x * y <= 0) throw SwitchToAlg2();
        return val_Srz + y * sqrt(z * z + x * y);
    }
    else {
        if (x * x + z * y <= 0) throw SwitchToAlg3();
        return y + Srz(z, x, y) * sqrt(x * x + z * y);
    }
}

double Qrz(double x, double y) {
    if (abs(x) < 1) return x * Srs(x, y, x);
    else return y * Srs1(y, x, y);
}

double Rrz(double x, double y, double z) {
    if (x > y) return x * z * Qrz(y, z);
    else return y * x * Qrz(x, y);
}

double Srs1(double x, double y, double z) {
    double val_Srz = Srz(x, y, z);
    if (z > y) return val_Srz + 1.44 * y * z;
    else return y + 1.44 * Srz(z, x, y);
}

double Qrz1(double x, double y) {
    if (abs(y) < 1) return x * Srs1(x, y, x);
    else return y * Srs1(y, x, y);
}

double Rrz_Alg2(double x, double y, double z) {
    if (x > y) return x * y * Qrz1(y, z);
    else return x * z * Qrz1(x, y);
}

double Srs2(double x, double y, double z) {
    double val_Srz = Srz(x, y, z);
    if (z > y) return val_Srz + y * x;
    else return y * z + Srz(z, x, y);
}

double Qrz2(double x, double y) {
    if (abs(x) < 1) return x * Srs2(x, y, x);
    else return y * Srs2(y, x, y);
}

double Rrz_Alg3(double x, double y, double z) {
    if (x > y) return x * y * Qrz2(y, z);
    else return y * z * Qrz2(x, y);
}

double Grs(double x, double y, double z) {
    double val_Rrz1, val_Rrz2;

    try {
        val_Rrz1 = Rrz(x, y, z);
    }
    catch (SwitchToAlg2) {
        val_Rrz1 = Rrz_Alg2(x, y, z);
    }
    catch (SwitchToAlg3) {
        val_Rrz1 = Rrz_Alg3(x, y, z);
    }

    try {
        val_Rrz2 = Rrz(x - y, z, y);
    }
    catch (SwitchToAlg2) {
        val_Rrz2 = Rrz_Alg2(x - y, z, y);
    }
    catch (SwitchToAlg3) {
        val_Rrz2 = Rrz_Alg3(x - y, z, y);
    }

    return 0.1389 * val_Rrz1 + 1.8389 * val_Rrz2;
}

double Fun_Alg4(double x, double y, double z) {
    return 1.3498 * x + 2.2362 * y * z - 2.348 * x * y;
}

double fun(double x, double y, double z) {
    return x * Grs(x, y, z) + y * Grs(x, z, y);
}

int main() {
    createFiles();

    double x, y, z, result;

    cout << "Input x, y, z: ";
    if (!(cin >> x >> y >> z)) {
        cout << "Invalid input." << endl;
        return 1;
    }

    try {
        result = fun(x, y, z);
        cout << "Result fun(x,y,z) = " << result << endl;
    }
    catch (ErrorNoFile& e) {
        e.Message();
        result = Fun_Alg4(x, y, z);
        cout << "Result (Alg 4) = " << result << endl;
    }
    catch (...) {
        cout << "Unknown error occurred." << endl;
    }

    return 0;
}