
#include <iostream>
#include <vector>
#include <rapidcsv.h>

using namespace std;

int main() {
    rapidcsv::Document doc("data/dipole_amplitude_with_IP_dependence.csv");

    vector<double> r = doc.GetColumn<double>("r");

    for (int i=0; i<10; i++) {
        cout << r[i] << endl;
    }
}