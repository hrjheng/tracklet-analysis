#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include "pdgidfunc.h"

int main(int argc, char *argv[])
{
    // int pid = 2112;
    int pid = 9010211;
    // int pid = 3312;
    // int pid = 25;
    // int pid = 10543;
    // int pid = 3000211;

    std::cout << "is quark = " << is_quark(pid) << std::endl
              << "is meson = " << is_meson(pid) << std::endl
              << "is baryon = " << is_baryon(pid) << std::endl
              << "is hadron = " << is_hadron(pid) << std::endl
              << "is pentaquark = " << is_pentaquark(pid) << std::endl
              << "is gauge boson or higgs = " << is_gauge_boson_or_higgs(pid) << std::endl
              << "charge = " << charge(pid) << std::endl
              << "is charged hadron = " << is_chargedHadron(pid) << std::endl;
    return 0;
}
