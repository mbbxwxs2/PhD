//
//  main.cpp
//  Fokker_Planck
//
//  Created by 宋晓东 on 30/06/2020.
//  Copyright © 2020 Xiaodong Song. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string.h>
#include<iomanip>
using namespace std;
//Lagrange Interpolation
double lagrangeInterpolation(const vector<double>& y, const vector<double>& x, double x0, int n)
{
    if (x.size()<n) return lagrangeInterpolation(y, x, x0, int(x.size()));
    if (n == 0) throw;
    int nHalf = n / 2;
    int jStar;
    double dx = x[1] - x[0];
    if (n % 2 == 0)
        jStar = int((x0 - x[0]) / dx) - (nHalf - 1);
    else
        jStar = int((x0 - x[0]) / dx + 0.5) - (nHalf);
    jStar = std::max(0, jStar);
    jStar = std::min(int(x.size() - n), jStar);
    if (n == 1)return y[jStar];
    double temp = 0.;
    for (unsigned int i = jStar; i<jStar + n; i++){
        double  int_temp;
        int_temp = y[i];
        for (unsigned int j = jStar; j<jStar + n; j++){
            if (j == i){ continue; }
            int_temp *= (x0 - x[j]) / (x[i] - x[j]);
        }
        temp += int_temp;
    }
    // end of interpolate
    return temp;
}
//Check the conservation law
double monitor(int n,vector<double> uNew)
{
    double result=0;
    for(int j=0;j<=n;j++) result+=uNew[j];
    return result;
}

//PDF of normal distribution
double normal_PDF(double x, double mu, double sigma)
{
    return(exp(-0.5*pow((x-mu)/sigma,2))/sigma/sqrt(2*M_PI));
}

//Find the lower and upper matrix

vector<double> thomasSolve(const vector<double> &a,const vector<double> &b,const vector<double> &c,vector<double> &d,int jMax)
{
    vector<double> rho(jMax+1),tau(jMax+1),temp(jMax+1);
    // initial first value of b
    rho[0]=d[0]/b[0];
    tau[0]=c[0]/b[0];
    for(int j=1;j<=jMax;j++)
    {
        tau[j]=c[j]/(b[j]-a[j]*tau[j-1]);
        rho[j]=(d[j]-a[j]*rho[j-1])/(b[j]-a[j]*tau[j-1]);
    }
    // calculate solution
    temp[jMax]=rho[jMax];
    for(int j=jMax-1;j>=0;j--)
        temp[j]=rho[j]-tau[j]*temp[j+1];
    return temp;
}

int main()
{
    auto start1 = std::chrono::steady_clock::now();
    // declare and initialise parameters
    //double theta=1.5,mu=0.5,sigma=1,gamma=0.8, T=10, x0=0.1,t_0=0.01;//con1
    //double theta=1,mu=8,sigma=1.2,gamma=0.7, T=10, x0=5, t_0=0.01;//con2
    //double theta=0.5,mu=5,sigma=0.5,gamma=0.6, T=10, x0=4,t_0=0.01;//con3
    //double theta=1.5,mu=0.5,sigma=1,gamma=0.3, T=10, x0=2,t_0=0.01;//con4
    //double theta=1,mu=9,sigma=1.2,gamma=0.2, T=10, x0=7, t_0=0.01;//con5
    //double theta=0.5,mu=5,sigma=0.5,gamma=0.1, T=10, x0=4,t_0=0.01;//con6
    double theta=1.5,mu=0.5,sigma=1,gamma=0.2, T=10, x0=0.2,t_0=0.01;//con7
    // declare and initialise grid paramaters
    int m=1000,n=500;
    double s=10;
    double dx=s/(n+1);
    double dt=(T-t_0)/m;
    double Single_Test=0.09;
    cout<<"dt= "<<dt<<endl;
    cout<<"dx= "<<dx<<endl;
    vector<double> conservation(m+1);
    vector<double> X(n+1);
    // setup and initialise x
    for(int j=0;j<=n;j++)
    {
        X[j] = (0.5+j)*dx;
    }
    //Declare nd initialize A, where A*U=D
    vector<double> integral(n+1);
    //Decomposite A=LU
    vector<vector<double>> lower(n+1,vector<double> (n+1,0)), upper(n+1,vector<double> (n+1,0));
    
    //Declare nd initialize u(x,t)
    vector<double> uOld(n+1),uNew(n+1);
    for (int j=0;j<=n;j++){
        uOld[j]=exp(-pow(X[j]-x0, 2)/(4*t_0))/(2*sqrt(M_PI*t_0));
        uNew[j]=exp(-pow(X[j]-x0, 2)/(4*t_0))/(2*sqrt(M_PI*t_0));
    }
    conservation[0]=1-monitor(n,uNew)*dx;
    cout<<"Energy Loss= "<<conservation[0]<<endl;
    //reset energy loss to zero
    for (int j=0;j<=n;j++){
        uOld[j]=uOld[j]/(1-conservation[0]);
        uNew[j]=uNew[j]/(1-conservation[0]);
    }
    conservation[0]=1-monitor(n,uNew)*dx;
    cout<<"Energy Loss= "<<conservation[0]<<endl;
    //Single point tests
    vector<double> single_point(m+1);
    single_point[0]=lagrangeInterpolation(uNew,X,Single_Test,4);
    // start looping through time levels
    vector<double> a_tile(n+1),b_tile(n+1),c_tile(n+1),d_tile(n+1);
    double para_a, para_b, para_c, para_d;
    for(int i=1;i<=m;i++)
    {
        // declare vectors for matrix equations
        para_a=-theta*(mu-X[0]-0.5*dx)+pow(sigma,2)*gamma*pow(X[0]+0.5*dx,2*gamma-1);
        para_b=pow(sigma,2)/dx*pow(X[0]+0.5*dx,2*gamma);
        a_tile[0]=0.;
        b_tile[0]=1-dt/4/dx*(para_a-para_b);
        c_tile[0]=-dt/4/dx*(para_a+para_b);
        d_tile[0]=dt/4/dx*(para_a+para_b)*uOld[1]+(1+dt/4/dx*(para_a-para_b))*uOld[0];
        for(int j=1;j<n;j++)
        {
            para_a=-theta*(mu-X[j]-0.5*dx)+pow(sigma,2)*gamma*pow(X[j]+0.5*dx,2*gamma-1);
            para_b=pow(sigma,2)/dx*pow(X[j]+0.5*dx,2*gamma);
            para_c=-theta*(mu-X[j]+0.5*dx)+pow(sigma,2)*gamma*pow(X[j]-0.5*dx,2*gamma-1);
            para_d=pow(sigma,2)/dx*pow(X[j]-0.5*dx,2*gamma);
            a_tile[j]=dt/4/dx*(para_c-para_d);
            b_tile[j]=1-dt/4/dx*(para_a-para_b-para_c-para_d);
            c_tile[j]=-dt/4/dx*(para_a+para_b);
            d_tile[j]=dt/4/dx*(para_a+para_b)*uOld[j+1]+(1+dt/4/dx*(para_a-para_b-para_c-para_d))*uOld[j]-dt/4/dx*(para_c-para_d)*uOld[j-1];
            
        }
        para_c=-theta*(mu-X[n]+0.5*dx)+pow(sigma,2)*gamma*pow(X[n]-0.5*dx,2*gamma-1);
        para_d=pow(sigma,2)/dx*pow(X[n]-0.5*dx,2*gamma);
        a_tile[n]=dt/4/dx*(para_c-para_d);
        b_tile[n]=1+dt/4/dx*(para_c+para_d);
        c_tile[n]=0;
        d_tile[n]=(1-dt/4/dx*(para_c+para_d))*uOld[n]-dt/4/dx*(para_c-para_d)*uOld[n-1];
        //Solve X from LUu=d
        uNew=thomasSolve(a_tile,b_tile,c_tile,d_tile,n);
        //Check the conservation law
        conservation[i]=1-monitor(n,uNew)*dx;
        //Set old=new
        uOld = uNew;
        single_point[i]=lagrangeInterpolation(uNew,X,Single_Test,4);
    }
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start1);
    cout<<"u(x*,T)= "<<setprecision(10)<<single_point[m]<<"\n";
    cout<<"Total time= "<<elapsed.count()<<endl;
    // finish looping through time levels
    ofstream output;
    //output.open("/Users/sxd/Downloads/FP_SOR.csv");
    output.open("/Users/sxd/Downloads/FP_mid.csv");
    if (!output.is_open())
    {
        cout<<"File not opened \n";
        throw;
    }
    // output the estimated u(x,T)
    for (int j=0;j<=n;j++)
    {
        output <<X[j] <<"," << uNew[j]<<endl;
        //    cout << "x: "<<X[j] <<", u(x,T)= " << uNew[j]<< "\n";
    }
    
    cout<<"The final conservation is "<<conservation[m]<<endl;
    output.close();
    output.open("/Users/sxd/Downloads/mid_conservation.csv");
    for (int i=0;i<=m;i++) output<<i*dt<<","<<conservation[i]<<"\n";
    output.close();
    
    output.open("/Users/sxd/Downloads/mid_single.csv");
    for (int i=0;i<=m;i++) output<<i*dt<<","<<single_point[i]<<"\n";
    output.close();
}
