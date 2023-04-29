#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>


float gamm = 1.414;       // Gas constant
double Me = 2;          //Exit Mach number

double ThroatRad = 1;
int nChar = 150;          //Number of characteristic curves


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

void save_coordinates(char filename[128], double coord);


double nu_m = prandtl_meyer(Me);
double thetaMax = nu_m*0.5;




//Class for A1, A2, A3 ..... 
class OriChar{

    public:
        double theta;
        double nu = theta;
        double M;
        double mu;
        double k_minus = 2.0*theta;
        double k_plus = 0;

};


class CplusPoint{

        public:
            double theta;
            double nu;
            double M;
            double mu;
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

    std::cout << " ############ STARTING NOZZLE DESIGN ########### " << std::endl;

    std::cout << "MAXIMUM WALL ANGLE:  " << thetaMax*180/3.14<< std::endl;
    std::cout << "EXIT MACH NUMBER: " << Me <<std::endl;
    std::cout << "NUMBER OF CHARACTERISTIC CURVES: " << nChar << std::endl;


    OriChar orichar[nChar];

    int totalPts = (nChar + 3)*(nChar)/2;
    CplusPoint cplusPoint[totalPts];

    for(int i = 0; i<nChar; i++){
    
        orichar[i].theta = (i+1)*(thetaMax)/(nChar);
        orichar[i].nu = orichar[i].theta;
        orichar[i].M = inverse_pm(orichar[i].nu);
        orichar[i].k_minus = 2*orichar[i].theta;
        orichar[i].k_plus = 0;
        orichar[i].mu = MachAngle(orichar[i].M);

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

        if(i==0){

            cplusPoint[i].theta = 0;
            cplusPoint[i].k_minus = orichar[i].k_minus;
            cplusPoint[i].nu = cplusPoint[i].k_minus;
            cplusPoint[i].k_plus = -cplusPoint[i].nu;

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            cplusPoint[i].coord.y = 0;

            //Calculate x coordinate
            double t1 = (0.5)*(cplusPoint[i].theta + orichar[i].theta);
            double m1 = (0.5)*(cplusPoint[i].mu + orichar[i].mu);

            cplusPoint[i].coord.x = -(ThroatRad)/tan(t1-m1);

        }

        else if(cplusPoint[i].isWallpoint){

            cplusPoint[i].theta = cplusPoint[i-1].theta;

            double A = (cplusPoint[i-1].theta + cplusPoint[i-1].mu);
            double B = (0.5)*(thetaMax + cplusPoint[i-1].theta);

            double Nr = cplusPoint[i-1].coord.y - ThroatRad - (cplusPoint[i-1].coord.x*tan(A));
            double Dr = tan(B) - tan(A);

            cplusPoint[i].coord.x = Nr/Dr;
            cplusPoint[i].coord.y = ThroatRad + (cplusPoint[i].coord.x*tan(B));

        }

        else{

            cplusPoint[i].k_minus = orichar[i].k_minus;
            cplusPoint[i].k_plus = cplusPoint[i-1].k_plus;
            cplusPoint[i].theta = (0.5)*(cplusPoint[i].k_minus + cplusPoint[i].k_plus);
            cplusPoint[i].nu = (0.5)*(cplusPoint[i].k_minus - cplusPoint[i].k_plus);

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            double A = (0.5)*(cplusPoint[i].theta + cplusPoint[i-1].theta + cplusPoint[i].mu + cplusPoint[i-1].mu);
            double B = (0.5)*(cplusPoint[i].theta + orichar[i].theta - cplusPoint[i].mu - orichar[i].mu);

            cplusPoint[i].coord.x = ((ThroatRad - cplusPoint[i-1].coord.y)+ (cplusPoint[i-1].coord.x*tan(A)))/(tan(A) - tan(B));
            cplusPoint[i].coord.y = ThroatRad + (cplusPoint[i].coord.x*tan(B));

        }

    }

    

    int rr = 0;     //Current RR number
    int rr_next = 0;

    //Computations for other RR

    for(int i = nChar+1; i<totalPts; i++){

        if(cplusPoint[i].isCentrepoint){

            rr+=1;

            /*if(i == totalPts-2){
                //Incase of last centre point
                cplusPoint[i].M = Me;
                cplusPoint[i].nu  = prandtl_meyer(Me);
                cplusPoint[i].theta = 0;
                cplusPoint[i].k_minus = cplusPoint[i].nu;
                cplusPoint[i].k_plus = -cplusPoint[i].nu;
                cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            }
            */

            //else{

                cplusPoint[i].theta = 0;
                cplusPoint[i].k_minus = cplusPoint[i - nChar + (rr-1)].k_minus;
                cplusPoint[i].nu = cplusPoint[i].k_minus;
                cplusPoint[i].k_plus = -cplusPoint[i].nu;

                cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
                cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            //}

            cplusPoint[i].coord.y = 0;

            //Calculate x coordinate
            double t1 = (0.5)*(cplusPoint[i].theta + cplusPoint[i - nChar + (rr-1)].theta);
            double m1 = (0.5)*(cplusPoint[i].mu + cplusPoint[i - nChar + (rr-1)].mu);

            cplusPoint[i].coord.x = cplusPoint[i- nChar + (rr-1)].coord.x - ((cplusPoint[i - nChar + (rr-1)].coord.y)/tan(t1-m1));

            //update rr_next
            rr_next++;

        }

        else if(cplusPoint[i].isWallpoint){

            cplusPoint[i].theta = cplusPoint[i-1].theta;

            double A = (cplusPoint[i-1].theta + cplusPoint[i-1].mu);
            double B = (0.5)*(cplusPoint[i - nChar + (rr-1)].theta + cplusPoint[i].theta);

            //Calculate coordinates

            double y_ = cplusPoint[i - nChar + (rr-1)].coord.y;
            double x_ = cplusPoint[i - nChar + (rr-1)].coord.x;

            cplusPoint[i].coord.x = ((y_ - cplusPoint[i-1].coord.y) + (cplusPoint[i-1].coord.x*tan(A)) - (x_*tan(B)))/(tan(A) - tan(B));
            cplusPoint[i].coord.y = y_ + ((cplusPoint[i].coord.x - x_)*tan(B));

        }

        else{

            cplusPoint[i].k_minus = cplusPoint[i - nChar + (rr-1)].k_minus;
            cplusPoint[i].k_plus = cplusPoint[i - 1].k_plus;
            cplusPoint[i].theta = (0.5)*(cplusPoint[i].k_minus + cplusPoint[i].k_plus);
            cplusPoint[i].nu = (0.5)*(cplusPoint[i].k_minus - cplusPoint[i].k_plus);

            cplusPoint[i].M = inverse_pm(cplusPoint[i].nu);
            cplusPoint[i].mu = MachAngle(cplusPoint[i].M);

            //Calculate coordinates
            double A = (0.5)*(cplusPoint[i].theta + cplusPoint[i-1].theta + cplusPoint[i].mu + cplusPoint[i-1].mu);
            double B = (0.5)*(cplusPoint[i].theta + cplusPoint[i-nChar+(rr-1)].theta - cplusPoint[i].mu - cplusPoint[i-nChar+(rr-1)].mu);

            double yo = cplusPoint[i-nChar+(rr-1)].coord.y;
            double xo = cplusPoint[i-nChar+(rr-1)].coord.x;
            double yprev = cplusPoint[i-1].coord.y;
            double xprev = cplusPoint[i-1].coord.x;

            cplusPoint[i].coord.x = (yo-yprev-(xo*tan(B)) + (xprev*tan(A)))/(tan(A) - tan(B));
            cplusPoint[i].coord.y = yprev + (cplusPoint[i].coord.x - xprev)*tan(A);


        }

    }

    //Save Wall Coordinates in a text file to be used for post processing

    // assume reasonably-sized file names
    char filename_x[128];
    char filename_y[128];
    const char *name_x = "x";
    const char *name_y = "y";


    sprintf(filename_x, name_x);
    sprintf(filename_y, name_y);
    FILE *output_x = fopen(filename_x, "w");
    FILE *output_y = fopen(filename_y, "w");

    fprintf(output_x, "%.10f\t", 0);
    fprintf(output_y, "%.10f\t", ThroatRad);

    for(int i = 0; i<totalPts; i++){

        if(cplusPoint[i].isWallpoint){
            fprintf(output_x, "%.10f\t", cplusPoint[i].coord.x);
            fprintf(output_y, "%.10f\t", cplusPoint[i].coord.y);
        }

    }

    fclose(output_x);
    fclose(output_y);

    std::cout << "NUMBER OF GRID POINTS COMPUTED: " << totalPts << std::endl;
    std::cout << "NOZZLE LENGTH: " << cplusPoint[totalPts-1].coord.x << " (times the throat radius)" << std::endl;
    std::cout << "Please run 'Nozzle_Contour.m' to visualize the nozzle's contour" << std::endl;;
    std::cout << "--------------------------------------------------------------" <<std::endl;

}

