#include <stdio.h>
#include <iostream>
//#include <stdlib.h>
//#define _USE_MATH_DEFINES
#include <math.h>


float gamm = 1.4;       // Gas constant
double Me = 3;          //Exit Mach number

double ThroatRad = 1;
int nChar = 200;          //Number of characteristic curves


// Prandlt meyer function
double prandtl_meyer(double M){

    double A = sqrt((gamm + 1)/(gamm - 1));
    double B = sqrt((M*M) - 1);
    double nu = A*atan(B/A) - atan(B);           // Arctan returns value in radians

    return nu;

}

//Function to obtain mach number from Prandtl Meyer angle
double inverse_pm(double nu){
    //nu must be in radians

    double diff = 0.001;
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

double MachAngle(double M){
    return asin(1/M);
}


double nu_m = prandtl_meyer(Me);
double thetaMax = nu_m*0.5;


//Class for A1, A2, A3 ..... 
class OriChar{

    public:
        double theta;
        double nu = theta;
        double M;
        double mu = asin(1.0/M);
        double k_minus = 2.0*theta;
        double k_plus = 0;

};


class CplusPoint{

        public:
            double theta;
            double nu;
            double M;
            double mu = asin(1.0/M);
            double k_minus;
            double k_plus;

            bool isCentrepoint = false;
            bool isWallpoint = false;
            
            struct coord{

                double x;
                double y;

            }coord;

};



int main(){

    std::cout << " -------- STARTING MOC -------- " << std::endl;

    OriChar orichar[nChar];

    int totalPts = (nChar + 3)*(nChar)/2;
    CplusPoint cplusPoint[totalPts];

    for(int i = 0; i<nChar; i++){

        orichar[i].theta = (i+1)*(thetaMax/nChar);
        orichar[i].nu = orichar[i].theta;
        orichar[i].M = inverse_pm(orichar[i].nu);
        orichar[i].k_minus = 2*orichar[i].theta;
        orichar[i].k_plus = 0;


        std::cout << "orichar: " << orichar[i].k_minus*180*(1/3.14) << std::endl;
    }

    int cp[nChar];      //Array to store indices of Centrepoints
    int wp[nChar];      //Array to store indices of Wall Points

    cp[0] = 0;
    wp[0] = nChar;

    cplusPoint[cp[0]].isCentrepoint = true;
    cplusPoint[wp[0]].isWallpoint = true;


    // Obtain indices of centre and wall points
    for(int rr = 1; rr<nChar; rr++){

        cp[rr] = cp[rr-1] + (nChar - rr + 2);
        wp[rr] = (cp[rr] + nChar - rr);

        cplusPoint[cp[rr]].isCentrepoint = true;
        cplusPoint[wp[rr]].isWallpoint = true;

    }


    // Computations for points in first RR

    for(int i = 0; i<nChar+1; i++){

        if(cplusPoint[i].isCentrepoint){

            cplusPoint[i].theta = 0;
            cplusPoint[i].k_minus = orichar[i].k_minus;
            cplusPoint[i].nu = cplusPoint[i].k_minus;
            cplusPoint[i].k_plus = -cplusPoint[i].nu;

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            cplusPoint[i].coord.y = 0;

            //std::cout << "CPlusChar: " << cplusPoint[i].nu << std::endl;

            //Calculate x coordinate

        }

        else if(cplusPoint[i].isWallpoint){

            cplusPoint[i].theta = cplusPoint[i-1].theta;
            //std::cout << "CPlusCharWall: " << cplusPoint[i].theta*180/3.14 << std::endl;

        }

        else{

            cplusPoint[i].k_minus = orichar[i].k_minus;
            cplusPoint[i].k_plus = cplusPoint[i-1].k_plus;
            cplusPoint[i].theta = (0.5)*(cplusPoint[i].k_minus + cplusPoint[i].k_plus);
            cplusPoint[i].nu = (0.5)*(cplusPoint[i].k_minus - cplusPoint[i].k_plus);

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            //std::cout << "theta: " << cplusPoint[i].M << std::endl;

        }

    }

    

    int rr = 0;     //Current RR number
    int rr_next = 0;

    //Computations for other RR

    for(int i = nChar+1; i<totalPts; i++){

        if(cplusPoint[i].isCentrepoint){

            rr+=1;

            cplusPoint[i].theta = 0;
            cplusPoint[i].k_minus = cplusPoint[i - nChar + (rr-1)].k_minus;
            cplusPoint[i].nu = cplusPoint[i].k_minus;
            cplusPoint[i].k_plus = -cplusPoint[i].nu;

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            cplusPoint[i].coord.y = 0;

            //Calculate x coordinate

            //update rr_next
            rr_next++;

        }

        else if(cplusPoint[i].isWallpoint){

            cplusPoint[i].theta = cplusPoint[i-1].theta;

            //Calculate coordinates

        }

        else{

            cplusPoint[i].k_minus = cplusPoint[i - nChar + (rr-1)].k_minus;
            cplusPoint[i].k_plus = cplusPoint[i - 1].k_plus;
            cplusPoint[i].theta = (0.5)*(cplusPoint[i].k_minus + cplusPoint[i].k_plus);
            cplusPoint[i].nu = (0.5)*(cplusPoint[i].k_minus - cplusPoint[i].k_plus);

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            //Calculate coordinates

        }

        //Condition to update rr

    }

    for(int i = 0; i<totalPts; i++){

        if(cplusPoint[i].isWallpoint){
            std::cout << cplusPoint[i].theta*180/3.14159 <<std::endl;

        }

        //std::cout << cplusPoint[i].theta*180/3.14159 <<std::endl;

    }



}


