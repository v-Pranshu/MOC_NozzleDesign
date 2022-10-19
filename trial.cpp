#include <stdio.h>
#include <iostream>
#include <stdlib.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <cmath>


float gamma = 1.4;        // Gas constant
double Me = 2;          //Exit Mach number

double ThroatRad = 1;
int nChar = 5;          //Number of characteristic curves


// Prandlt meyer function
double prandtl_meyer(double M){

    double A = sqrt((gamma + 1)/(gamma - 1));
    double B = sqrt((M*M) - 1);
    double nu = A*atan(B/A) - atan(B);           // Arctan returns value in radians

    return nu;

}

//Function to obtain mach number from Prandtl Meyer angle
double inverse_pm(double nu){ 
    //Nu in radians

    double diff = 0.01;
    double Mach = 0.5;
    double nu_calc;
    double smallest  = 100;
    double M = 0;

    while(Mach <= 10){

        nu_calc = prandtl_meyer(Mach);
        double error = nu_calc - nu;

        if(abs(error) < smallest){
            smallest = abs(error);
            M = Mach;
        }

        Mach = Mach + diff;
    }

    return M;
}


int main(){

    double mach = inverse_pm(0.087);
    std::cout << mach << std::endl;

}