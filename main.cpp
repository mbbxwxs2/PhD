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
#include <iomanip>
using namespace std;


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


//Check the conservation law (trapezium rule)
double monitor(int n,vector<double> vNew,vector<double> &A)
{
    double result=vNew[0]*A[0]*0.5;
    for(int j=1;j<n;j++) result+=vNew[j]*A[j];
    result+=vNew[n]*A[n]*0.5;
    return result;
}

double calculation(int n, int m, double kappa, double theta, double sigma, double gamma, double T, double x0, double t_0, double s, double r, double Single_Test, vector<double> &vOld, vector<double> &vNew, vector<double> &conservation, vector<double> &single_point,vector<double> &X)
{
    double dx=(s-r)/m;
    double dt=(T-t_0)/n;
    vector<double> A(m+1),W(m),C(m),omega(m);
    // setup and initialise X
    for (int j=0;j<m;j++)
    {
        X[j] =r+j*dx;
        A[j]=exp((1-2*gamma)*X[j]);
        double B=-kappa*theta*exp(-2*gamma*(X[j]+0.5*dx))+kappa*exp((1-2*gamma)*(X[j]+0.5*dx));
        C[j]=pow(sigma,2)/2*exp(-X[j]-0.5*dx); //C[0 - m] stand for C_0, 1, ..., m
        omega[j]=dx*B/C[j];
        W[j]=omega[j]/(exp(omega[j])-1); //W[0 - m] stand for W_0, ..., m
        //_______________________________________________________________________________//
        //double delta=1/omega[j]-1/(exp(omega[j])-1);
        //if (delta>0.5 || delta<0) cout<<"j= "<<j<<", delta= "<<delta<<endl;
    }
    X[m]=s;
    A[m]=exp((1-2*gamma)*X[m]);
    //Declare nd initialize v(X,t)
    for (int j=0;j<=m;j++){
        vOld[j]=exp(2*gamma*X[j])*exp(-pow(exp(X[j])-x0, 2)/(4*t_0))/(2*sqrt(M_PI*t_0));
        vNew[j]=exp(2*gamma*X[j])*exp(-pow(exp(X[j])-x0, 2)/(4*t_0))/(2*sqrt(M_PI*t_0));
    }
    conservation[0]=1-monitor(m,vNew,A)*dx;
    //cout<<"Energy Loss= "<<conservation[0]<<endl;
    //reset energy loss to zero
    for (int j=0;j<=m;j++){
        vOld[j]=vOld[j]/(1-conservation[0]);
        vNew[j]=vNew[j]/(1-conservation[0]);
    }
    conservation[0]=1-monitor(m,vNew,A)*dx;
    //cout<<"Energy Loss= "<<conservation[0]<<endl;
    // start looping through time levels
    vector<double> a_tile(m+1),b_tile(m+1),c_tile(m+1),d_tile(m+1);
    
    //Single point tests
    double value=lagrangeInterpolation(vNew,X,Single_Test,4);
    single_point[0]=exp(-2*gamma*Single_Test)*value;
    for(int i=1;i<=n;i++)
    {
        // declare vectors for matrix equations
        a_tile[0]=0.;
        b_tile[0]=dt/A[0]/pow(dx,2)*C[0]*W[0]+1;
        c_tile[0]=-dt/A[0]/pow(dx,2)*C[0]*W[0]*exp(omega[0]);
        d_tile[0]=vOld[0];
        for (int j=1;j<m;j++)
        {
            a_tile[j]=-dt/A[j]/pow(dx,2)*C[j-1]*W[j-1];
            b_tile[j]=dt/A[j]/pow(dx,2)*(C[j]*W[j]+C[j-1]*W[j-1]*exp(omega[j-1]))+1;
            c_tile[j]=-dt/A[j]/pow(dx,2)*C[j]*W[j]*exp(omega[j]);
            d_tile[j]=vOld[j];
        }
        a_tile[m]=-dt/A[m]/pow(dx,2)*C[m-1]*W[m-1];
        b_tile[m]=dt/A[m]/pow(dx,2)*C[m-1]*W[m-1]*exp(omega[m-1])+1;
        c_tile[m]=0;
        d_tile[m]=vOld[m];
        vNew=thomasSolve(a_tile,b_tile,c_tile,d_tile,m);
        //Check the conservation law
        conservation[i]=1-monitor(m,vNew,A)*dx;
        //Set old=new
        vOld = vNew;
        
        // this is a cubic interpolation of V(S_0) using the 4 closest points (n=4)
        double value=lagrangeInterpolation(vNew,X,Single_Test,4);
        single_point[i]=exp(-2*gamma*Single_Test)*value;
    }
    return dt;
}


int main()
{
    auto start = std::chrono::steady_clock::now();
    //double kappa=1.5,theta=0.5,sigma=1,gamma=0.8, T=10, x0=0.1,t_0=0.01;//con1
    //double kappa=1,theta=8,sigma=1.2,gamma=0.7, T=10, x0=5,t_0=0.01;//con2
    //double kappa=0.5,theta=5,sigma=0.5,gamma=0.6, T=10, x0=4,t_0=0.01;//con3
    //double kappa=1.5,theta=0.5,sigma=1,gamma=0.3, T=10, x0=2,t_0=0.01;//con4
    //double kappa=1,theta=9,sigma=1.2,gamma=0.2, T=10, x0=7,t_0=0.01;//con5
    double kappa=0.5,theta=5,sigma=0.5,gamma=0.1, T=10, x0=4,t_0=0.01;//con6
    //double kappa=1.5,theta=0.5,sigma=1,gamma=0.2, T=10, x0=0.2,t_0=0.01;//con7
    // declare and initialise grid paramaters
    int n=6400,m=3200;
    double s=log(10);
    double r=-10;
    double Single_Test=5;
    Single_Test=log(Single_Test);
    vector<double> X_1(m+1), X_2(m/2+1),X_3(m/4+1);
    vector<double> vOld_1(m+1,0),vNew_1(m+1,0),vOld_2(m/2+1,0),vNew_2(m/2+1,0),vOld_3(m/4+1,0),vNew_3(m/4+1,0);
    vector<double> conservation_1(n+1), conservation_2(n+1),conservation_3(n+1);
    vector<double> single_point_1(n+1),single_point_2(n+1),single_point_3(n+1),single_point(n+1);
    double dt=calculation(n, m, kappa, theta, sigma, gamma, T, x0,t_0, s, r, Single_Test, vOld_1, vNew_1, conservation_1, single_point_1,X_1);
    auto finish_1 = std::chrono::steady_clock::now();
    auto elapsed_1 = std::chrono::duration_cast<std::chrono::duration<double> >(finish_1 - start);
    dt=calculation(n, m/2, kappa, theta, sigma, gamma, T, x0,t_0, s, r, Single_Test, vOld_2, vNew_2, conservation_2, single_point_2,X_2);
    for (int i=0;i<=n;i++){
        single_point[i]=(single_point_1[i]*4-single_point_2[i])/3;
    }
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
    dt=calculation(n, m/4, kappa, theta, sigma, gamma, T, x0,t_0, s, r, Single_Test, vOld_3, vNew_3, conservation_3, single_point_3,X_3);
    cout<<"n= "<<n<<" , "<<"dt= "<<(T-t_0)/n<<endl;
    cout<<"m= "<<m<<" , "<<"dx= "<<(s-r)/m<<endl;
    cout<<"u(x*,T)= "<<setprecision(10)<<single_point_1[n]<<"\n";
    cout<<"Single time= "<<elapsed_1.count()<<endl;
    cout<<"u_e(x*,T)= "<<setprecision(10)<<single_point[n]<<"\n";
    cout<<"Total time= "<<elapsed.count()<<endl;
    
    cout<<"u(x*,T)= "<<setprecision(10)<<single_point_1[n]<<"\n";
    cout<<"u(x*,T)= "<<setprecision(10)<<single_point_2[n]<<"\n";
    cout<<"u(x*,T)= "<<setprecision(10)<<single_point_3[n]<<"\n";
    cout<<"R= "<<setprecision(6)<<(single_point_2[n]-single_point_3[n])/(single_point_1[n]-single_point_2[n])<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[3200]-single_point_3[3200])/(single_point_1[3200]-single_point_2[3200])<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[400]-single_point_3[400])/(single_point_1[400]-single_point_2[400])<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[2]-single_point_3[2])/(single_point_1[2]-single_point_2[2])<<endl;
    // finish looping through time levels
    ofstream output;
    //output.open("/Users/sxd/Downloads/FP_SOR.csv");
    output.open("/Users/sxd/Downloads/Chang_cooper.csv");
    if (!output.is_open())
    {
        cout<<"File not opened \n";
        throw;
    }
    // output the estimated u(x,T)
    vector<double> uNew(m+1,0);
    for (int j=0;j<=m;j++)
    {
        uNew[j]=exp((-2*gamma)*X_1[j])*vNew_1[j];
        output <<exp(X_1[j]) <<"," << uNew[j]<<endl;
        //output <<X[j] <<"," << vNew[j]<<endl;
        //cout << "x: "<<X[j] <<", u(x,T)= " << uNew[j]<< "\n";
    }
    cout<<"The final conservation is "<<conservation_1[n]<<endl;
    output.close();
    output.open("/Users/sxd/Downloads/Chang_cooper_conservation.csv");
    for (int i=0;i<=n;i++) output<<i*dt<<","<<conservation_1[i]<<"\n";
    output.close();
    
    output.open("/Users/sxd/Downloads/Chang_cooper_single.csv");
    for (int i=0;i<=n;i++) output<<i*dt<<","<<setprecision(10)<<single_point[i]<<"\n";
    output.close();
}

