
#include "dipole_amp_reader.h"
#include <chrono>

int main() {

    static array<array<array<array<double, 4>, 81>, 900>, 30> table;
    string filename = "data/dipole_amplitude_with_IP_dependence.csv";
    
    load_dipole_amplitudes(table, filename);
    
    for (int i=0; i<30; i++) {
        for (int j=0; j<900; j++) {
            for (int k=0; k<81; k++) {
                cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
                cin;
            }
        }
    }
    /*
    auto start = chrono::high_resolution_clock::now();
    double res = get_dipole_amplitude(table, 121.77, 11.84, 1e-05);
    auto stop = chrono::high_resolution_clock::now();

    auto duration = chrono::duration_cast<chrono::nanoseconds>(stop-start);
    cout << "duration: " << duration.count() << endl;

    cout << "amplitude: " << res << endl;
    */
    //0.000155246, 2.43256e-05, 1.12535e-09, 4.95126e-06
}
