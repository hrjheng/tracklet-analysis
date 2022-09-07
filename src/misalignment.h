
#ifndef MISALIGNMENT_H
#define MISALIGNMENT_H

#include <map>
#include <vector>

#include <TMath.h>

vector<float> misalignment(int c)
{
    float shift_1 = 0.001, shift_2 = 0.005, shift_3 = 0.01;

    std::map<int, std::vector<float>> Map_MisliagnCase = {
        {0, {0, 0, 0}},
        {1, {shift_1, 0, 0}}, // 10 micron
        {2, {-1*shift_1, 0, 0}},
        {3, {shift_2, 0, 0}}, // 50 micron
        {4, {-1*shift_2, 0, 0}},
        {5, {shift_3, 0, 0}}, // 100 micron
        {6, {-1*shift_3, 0, 0}},
        {7, {0, shift_1, 0}}, // angle = 0
        {8, {0, shift_1, TMath::Pi()/6.}}, // angle = 30
        {9, {0, shift_1, TMath::PiOver4()}}, // angle = 45
        {10, {0, shift_1, TMath::Pi()/3.}}, // angle = 60
        {11, {0, shift_1, TMath::Pi()/2.}}, // angle = 90
        {12, {0, shift_1, 3 * TMath::PiOver4()}},
        {13, {0, shift_1, 5 * TMath::PiOver4()}},
        {14, {0, shift_1, 7 * TMath::PiOver4()}},
        {15, {0, shift_2, 0}}, // angle = 0
        {16, {0, shift_2, TMath::Pi()/6.}}, // angle = 30
        {17, {0, shift_2, TMath::PiOver4()}}, // angle = 45
        {18, {0, shift_2, TMath::Pi()/3.}}, // angle = 60
        {19, {0, shift_2, TMath::Pi()/2.}}, // angle = 90
        {20, {0, shift_2, 3 * TMath::PiOver4()}},
        {21, {0, shift_2, 5 * TMath::PiOver4()}},
        {22, {0, shift_2, 7 * TMath::PiOver4()}},
        {23, {0, shift_3, 0}}, // angle = 0
        {24, {0, shift_3, TMath::Pi()/6.}}, // angle = 30
        {25, {0, shift_3, TMath::PiOver4()}}, // angle = 45
        {26, {0, shift_3, TMath::Pi()/3.}}, // angle = 60
        {27, {0, shift_3, TMath::Pi()/2.}}, // angle = 90
        {28, {0, shift_3, 3 * TMath::PiOver4()}},
        {29, {0, shift_3, 5 * TMath::PiOver4()}},
        {30, {0, shift_3, 7 * TMath::PiOver4()}},
        {31, {shift_1, shift_1, TMath::PiOver4()}},
        {32, {shift_1, shift_2, TMath::PiOver4()}},
        {33, {shift_1, shift_3, TMath::PiOver4()}},
        {34, {shift_2, shift_1, TMath::PiOver4()}},
        {35, {shift_2, shift_2, TMath::PiOver4()}},
        {36, {shift_2, shift_3, TMath::PiOver4()}},
        {37, {shift_3, shift_1, TMath::PiOver4()}},
        {38, {shift_3, shift_2, TMath::PiOver4()}},
        {39, {shift_3, shift_3, TMath::PiOver4()}},
    };

    return Map_MisliagnCase[c];
}

#endif