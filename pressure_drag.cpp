#include <stdio.h>
#include <iostream>
#include <math.h>
#include <cmath>
#include <vector>

const int num_points = 80;  //Put even number
double A_throat = 1.0;
double gamma = 1.41;
double Mach_a = 0.4;
double pa = 1;    // in kPa

class Control{
    public:
        double x;
        double y;

};

double ex = (gamma + 1) / (2 * (gamma - 1));

double K = pow(0.5 * (gamma + 1), -ex);

double P_0a = pa * pow((1 + 0.5*(gamma  - 1)*pow(Mach_a, 2)), gamma/(gamma - 1));

class P{

    public:
        double x;
        double y;

        double A ;

        double ratio;

        double alpha;
        double M;
        double pressure;

};

double Bezier_x(Control P0,  Control P1, Control P2, Control P3, double t){

    double t1 = P0.x*( - pow(t, 3) + 3*pow(t, 2) - 3*t + 1);
    double t2 = P1.x * (3*pow(t, 3) - 6 * pow(t, 2) + 3 * t);
    double t3 = P2.x * (-3*pow(t, 3) + 3 * pow(t, 2));
    double t4 = P3.x * (pow(t, 3));

    double P_x = t1+t2+t3+t4;

    return P_x;

}

double Bezier_y(Control P0, Control P1, Control P2, Control P3, double t)
{

    double t1 = P0.y * (-pow(t, 3) + 3 * pow(t, 2) - 3 * t + 1);
    double t2 = P1.y * (3 * pow(t, 3) - 6 * pow(t, 2) + 3 * t);
    double t3 = P2.y * (-3 * pow(t, 3) + 3 * pow(t, 2));
    double t4 = P3.y * (pow(t, 3));

    double P_y = t1 + t2 + t3 + t4;

    return P_y;
}

int main(){

        Control P0, P1, P2, P3;
        P0.x = 0.0;
        P0.y = 0.0;

        P1.x = 1.0;
        P1.y = 2.0;

        P2.x = 1.0;
        P2.y = 4.0;

        P3.x = 0.0;
        P3.y = 6.0;

        P point[num_points];

        //Find Coordinates of all the necessary evaluation points
        for(int i = 0; i<num_points; i++){

            double t = (double)i/ (double)(num_points-1);
            point[i].x = Bezier_x(P0, P1, P2, P3, t);
            point[i].y = Bezier_y(P0, P1, P2, P3, t);

        }

        //Compute pressure at all these points

        for (int i = 0; i < num_points; i++)
        {            
            if(i < num_points/2){
                point[i].A = point[i].y + A_throat;
                point[i].ratio = point[i].A / A_throat;
            }

            else{
                point[i].A = A_throat + (P3.y - point[i].y);
                point[i].ratio = point[i].A / A_throat;
            }

            point[i].alpha = pow(point[i].ratio / K, (1 / ex));
            point[i].M = (point[i].alpha - sqrt(point[i].alpha * point[i].alpha - 2 * (gamma - 1))) / (gamma - 1);
            point[i].pressure = P_0a / (pow((1 + 0.5 * (gamma - 1) * pow(point[i].M, 2)), gamma / (gamma - 1)));

        }

        //Compute Drag by summing pressure*Dy
        double drag = 0.0;

        for(int i = 0; i<(num_points); i++){
            drag += (point[i].pressure) * (point[i].y - point[i - 1].y);
        }

        std::cout << "DRAG FORCE TOTAL : "<< drag << "kN"<<std::endl;
        
}