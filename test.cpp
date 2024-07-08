
#include "dipole_amp_reader.h"


int main() {

    static array<array<array<array<double, 4>, 81>, 900>, 30> table;
    string filename = "data/dipole_amplitude_with_IP_dependence.csv";

    load_dipole_amplitude(table, filename);

    for (int i=0; i<3; i++) {
        for (int j=0; j<5; j++) {
            for (int k=0; k<5; k++) {
                cout << table[i][j][k][0] << ", " << table[i][j][k][1] << ", " << table[i][j][k][2] << ", " << table[i][j][k][3] << endl;
            }
        }
    }
}
