#pragma once
/************************************************************************************************/
/*** Topic: Dynamic Runge-Kutta 4th Order Method application in the class for each step       ***/
/***         whitout time dependent							                                  ***/
/***        solved numerically using RungeKutta 4th order method with adaptive time step size ***/
/***           Explosive synchronization in interlayer phase-shifted Kuramoto oscillators on  ***/
/***             multiplex networks     --Sarika Jalan--                                      ***/
/*** Version Release 17.12 rev 11256                                                Ali-Seif  ***/
/*** Date: 8/19/2022                                                                          ***/
/*** Code implemented in Microsoft Visual Studio Enterprise 2019 C++ compiler                 ***/
/*** MSI: PX60 6QD/ DDR4                                                                      ***/
/*** Run under a Intel® Core™ i7-6700HQ CPU @ 2.60GHz × 64 based processor with 16 GB RAM     ***/
/*** in the computational server manual parallel                                              ***/
/************************************************************************************************/

/*
#include"RK4.h"                                                            //import library Kuramoto

int main() {                                                                    //Beginning main
    double t_last = 4.0, dt_min = 0.001, y_tol = 0.001, dy_min = 0.008, dy_max = 0.01;
    RK4 RK(t_last, dt_min, y_tol, dy_min, dy_max);
    double t = 0.0, y = 0.0, dt = 0.01;


    //y = RK.DRK4(dt, y);//double (dt,y) //Dynamic Runge-Kutta 4th Order
    //cout << y << endl;


    y = RK.CRK4(dt, y);//double (dt,y)  //Constant Runge-Kutta 4th Order
    cout << y << endl;


    return 0;                                                                   //run program was correct
}                                                                               //

*/

//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
#include<iostream>                                                              //for cout                              //$$$$
#include<fstream>                                                               //infile /ofstream                      //$$$$
using namespace std;                                                                                                    //$$$$
double mfun(double y)                                                                                                   //$$$$
{                                                                                                                       //$$$$
    //double mf =y-y+5; //for this example t_last=40                                                                    //$$$$
    double mf = exp(y) - y*exp(y); //for this example t_last=4                                                          //$$$$
    return mf;                                                                                                          //$$$$
}                                                                                                                       //$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$



//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@                                                                                @@@@
//@@@                                class Kuramoto                                  @@@@
//@@@                                                                                @@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
class RK4 {                                                                //create and define class 
private:                                                                        //private values
    double t_last;                                                                  //step lenth
    double dt_min;                                                              //time final
    double y_tol;                                                            //stringh cupling
    double dy_min;
    double dy_max;
public:                                                                         //public values
    RK4                                                                    //
    (double t_last_, double dt_min_, double y_tol_, double dy_min_, double dy_max_) ://input data in class
        t_last(t_last_), dt_min(dt_min_), y_tol(y_tol_), dy_min(dy_min_), dy_max(dy_max_) {//change name input to privet


    }
    double  DRK4(double,double);      //Dynamic Runge-Kutta 4th Order
    double  CRK4(double, double);     //Constant Runge-Kutta 4th Order

};



//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$                                                                                                                  $$$$$
//$$$$$                                                    DRK4                                                          $$$$$
//$$$$$                                                                                                                  $$$$$
double RK4::DRK4(double dt,double y) {                                                                                  //$$$$
    //ofstream Print_DRK4("DRK4.txt");                                                                                  //$$$$
        double k1, k2, k3, k4, sy, hsy, dsy, ny;                                                                        //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@                                    Nurmal step                                 @@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        k1 = mfun(y);                                                                                                   //$$$$
        k2 = mfun(y + dt * k1 / 2);                                                                                     //$$$$
        k3 = mfun(y + dt * k2 / 2);                                                                                     //$$$$
        k4 = mfun(y + dt * k3);                                                                                         //$$$$
        sy = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);                                                                  //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@                                    Half step                                   @@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        k2 = mfun(y + dt * k1 / 4);                                                                                     //$$$$
        k3 = mfun(y + dt * k2 / 4);                                                                                     //$$$$
        k4 = mfun(y + dt * k3 / 2);                                                                                     //$$$$
        hsy = y + dt / 12 * (k1 + 2 * k2 + 2 * k3 + k4);                                                                //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@                                   Double step                                  @@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        k2 = mfun(y + dt * k1);                                                                                         //$$$$
        k3 = mfun(y + dt * k2);                                                                                         //$$$$
        k4 = mfun(y + dt * k3 * 2);                                                                                     //$$$$
        dsy = y + dt / 3 * (k1 + 2 * k2 + 2 * k3 + k4);                                                                 //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@                                  compareand use                                @@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        if (abs(sy) < y_tol) {                                                              //@@@                       //$$$$
            if (dt = !dt_min) {                                                             //@@@                       //$$$$
                cout << "New step size=\t" << dt_min<<endl;                                 //@@@                       //$$$$
                dt = dt_min;                                                                //@@@                       //$$$$
            }                                                                               //@@@                       //$$$$
            ny = sy;                                                                        //@@@                       //$$$$
        }else {                                                                             //@@@                       //$$$$
            if (abs(sy) > y_tol && (abs(sy - hsy) / abs(sy)) > dy_max) {                    //@@@                       //$$$$
                dt = dt / 2;                                                                //@@@                       //$$$$
                cout << "New step size=\t" << dt << endl;                                   //@@@                       //$$$$
                ny = hsy;                                                                   //@@@                       //$$$$
            }                                                                               //@@@                       //$$$$
            else if (abs(sy) > y_tol && (abs(sy - dsy) / abs(sy)) < dy_min) {               //@@@                       //$$$$
                dt = dt * 2;                                                                //@@@                       //$$$$
                cout << "New step size=\t" << dt << endl;                                   //@@@                       //$$$$
                ny = dsy;                                                                   //@@@                       //$$$$
            }else {                                                                         //@@@                       //$$$$
                ny = sy;                                                                    //@@@                       //$$$$
            }                                                                               //@@@                       //$$$$
        }                                                                                   //@@@                       //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        return ny;                                                                                                      //$$$$
        //Print_DRK4 << t << '\t' << y << endl;                                                                         //$$$$
}                                                                                                                       //$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$


//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
//$$$$$                                                                                                                  $$$$$
//$$$$$                                                    CRK4                                                          $$$$$
//$$$$$                                                                                                                  $$$$$
double  RK4::CRK4(double dt, double y) {                                                                                //$$$$
    //ofstream Print_CRK4("CRK4.txt");                                                                                  //$$$$
        double k1, k2, k3, k4, sy;                                                                                      //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@                                    Nurmal step                                 @@@@                       //$$$$
        //@@@                                                                                @@@@                       //$$$$
        //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@                       //$$$$
        k1 = mfun(y);                                                                                                   //$$$$
        k2 = mfun(y + dt * k1 / 2);                                                                                     //$$$$
        k3 = mfun(y + dt * k2 / 2);                                                                                     //$$$$
        k4 = mfun(y + dt * k3);                                                                                         //$$$$
        sy = y + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);                                                                  //$$$$
        return sy;                                                                                                      //$$$$
       // Print_CRK4 << t << '\t' << y << endl;                                                                         //$$$$
}                                                                                                                       //$$$$
//$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$