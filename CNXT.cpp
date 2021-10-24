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
//-------------------------------------------------------------------------------------------------------------------//
//Check the conservation law (trapezium rule)
double monitor(int n,vector<double> vNew,const vector<double> &X,double gamma)
{
    double result=vNew[0]*exp((1-2*gamma)*X[0])*0.5;
    for(int j=1;j<n;j++) result+=vNew[j]*exp((1-2*gamma)*X[j]);
    result+=vNew[n]*exp((1-2*gamma)*X[n])*0.5;
    return result;
}
//-------------------------------------------------------------------------------------------------------------------//
//Find the lower and upper matrix
void LU_find(vector<vector<double>> &A,vector<vector<double>> &lower, vector<vector<double>> &upper, int n)
{
    // Decomposing matrix into Upper and Lower
    // triangular matrix
    for (int i = 0; i <=n; i++)
    {
        // Upper Triangular
        for (int k = i; k <=n; k++)
        {
            // Summation of L(i, j) * U(j, k)
            double sum = 0;
            for (int j = 0; j < i; j++)
                sum += (lower[i][j] * upper[j][k]);
            // Evaluating U(i, k)
            upper[i][k] = A[i][k] - sum;
        }
        // Lower Triangular
        for (int k = i; k <=n; k++)
        {
            if (i == k)
                lower[i][i] = 1; // Diagonal as 1
            else {
                
                // Summation of L(k, j) * U(j, i)
                double sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (lower[k][j] * upper[j][i]);
                // Evaluating L(k, i)
                lower[k][i] = (A[k][i] - sum) / upper[i][i];
            }
        }
    }
}
//-------------------------------------------------------------------------------------------------------------------//
//LU matrix solver
vector<double> LU_solver(vector<vector<double>> &lower,vector<vector<double>> &upper,const vector<double> &d,int n)
{
    vector<double> y(n+1),output(n+1);
    y[0]=d[0];
    for (int i=1;i<=n;i++)
    {
        y[i]=d[i];
        for (int j=0;j<i;j++) y[i]-=y[j]*lower[i][j];
    }
    for (int i=n;i>=0;i--)
    {
        output[i]=y[i];
        for (int j=n;j>i;j--) output[i]-=upper[i][j]*output[j];
        output[i]=output[i]/upper[i][i];
    }
    return(output);
}
//-------------------------------------------------------------------------------------------------------------------//
double Calculation(int m, int n, double s, double r, double T, double t_0, double kappa, double theta, double gamma,double sigma, double x0, double Single_Test, vector<double> &vOld, vector<double> &vNew, vector<double> &X, vector<double> &single_point, vector<double> &conservation)
{
    double dx=(s-r)/m;
    double dt=(T-t_0)/n;
    vector<double> A_tile(m+1),B(m),C(m);
    //setup and initialise X
    for (int j=0;j<m;j++)
    {
        X[j]=r+j*dx;
        A_tile[j]=exp((1-2*gamma)*X[j]);
        B[j]=-kappa*theta*exp(-2*gamma*(X[j]+0.5*dx))+kappa*exp((1-2*gamma)*(X[j]+0.5*dx));
        C[j]=0.5*pow(sigma,2)*exp(-X[j]-0.5*dx);
    }
    X[m]=r+m*dx;
    A_tile[m]=exp((1-2*gamma)*X[m]);
    //Declare 2-dimensional matrix A, where Av=D
    vector<vector<double>> A(m+1,vector<double> (m+1,0));
    A[0][0]=0.5*exp((1-2*gamma)*X[0])*dx;
    for (int j1=1;j1<m;j1++)
    {
        A[0][j1]=exp((1-2*gamma)*X[j1])*dx;
    }
    A[0][m]=0.5*exp((1-2*gamma)*X[m])*dx;
    for (int j=1;j<m;j++)
    {
        A[j][j-1]=-1/(2*A_tile[j]*dx)*(-B[j-1]/2+C[j-1]/dx);
        A[j][j]=1/dt-1/(2*A_tile[j]*dx)*(B[j]/2-C[j]/dx-B[j-1]/2-C[j-1]/dx);
        A[j][j+1]=-1/(2*A_tile[j]*dx)*(B[j]/2+C[j]/dx);
    }
    //Boundary
    A[m][m-1]=1/(2*A_tile[m]*dx)*(B[m-1]/2-C[m-1]/dx);
    A[m][m]=1/dt+1/(2*A_tile[m]*dx)*(B[m-1]/2+C[m-1]/dx);
    
    //Decomposite A=LU
    vector<vector<double>> lower(m+1,vector<double> (m+1,0)), upper(m+1,vector<double> (m+1,0));
    LU_find(A, lower, upper, m);
    //Declare nd initialize v(X,t)
    
    for (int j=0;j<=m;j++){
        vOld[j]=exp(2*gamma*X[j])*exp(-pow(exp(X[j])-x0, 2)/(4*t_0))/(2*sqrt(M_PI*t_0));
        vNew[j]=exp(2*gamma*X[j])*exp(-pow(exp(X[j])-x0, 2)/(4*t_0))/(2*sqrt(M_PI*t_0));
    }
    conservation[0]=1-monitor(m,vNew,X,gamma)*dx;
    //cout<<"Energy Loss= "<<conservation[0]<<endl;
    //reset energy loss to zero
    for (int j=0;j<=m;j++){
        vOld[j]=vOld[j]/(1-conservation[0]);
        vNew[j]=vNew[j]/(1-conservation[0]);
    }
    conservation[0]=1-monitor(m,vNew,X,gamma)*dx;
    // finish looping through time levels
    vector<double> d(m+1);
    vector<double> a_tile_2(m+1),b_tile(m+1),c_tile(m+1),d_tile(m+1);
    //Single point tests
    double value=lagrangeInterpolation(vNew,X,Single_Test,4);
    single_point[0]=exp(-2*gamma*Single_Test)*value;
    for (int i=1;i<=n;i++)
    {
        d[0]=1;
        for (int j=1;j<m;j++)
        {
            d[j]=1/(2*A_tile[j]*dx)*(-B[j-1]/2+C[j-1]/dx)*vOld[j-1]+(1/dt+1/(2*A_tile[j]*dx)*(B[j]/2-C[j]/dx-B[j-1]/2-C[j-1]/dx))*vOld[j]+1/(2*A_tile[j]*dx)*(B[j]/2+C[j]/dx)*vOld[j+1];
        }
        d[m]=-1/(2*A_tile[m]*dx)*(B[m-1]/2-C[m-1]/dx)*vOld[m-1]+(1/dt-1/(2*A_tile[m]*dx)*(B[m-1]/2+C[m-1]/dx))*vOld[m];
        vNew=LU_solver(lower,upper, d,m);
        //Check the conservation law
        conservation[i]=1-monitor(m,vNew,X,gamma)*dx;
        //Set old=new
        vOld = vNew;
        // this is a cubic interpolation of V(S_0) using the 4 closest points (n=4)
        double value=lagrangeInterpolation(vNew,X,Single_Test,4);
        single_point[i]=exp(-2*gamma*Single_Test)*value;
    }
    return dt;
}
//-------------------------------------------------------------------------------------------------------------------//
int main()
{
    auto start = std::chrono::steady_clock::now();
    double kappa=0.5,theta=5,sigma=0.5,gamma=0.1, T=10, x0=4,t_0=0.01;//con6
    // declare and initialise grid paramaters
    int n=6400,m=400;
    double s=log(10);
    double r=-10;
    double Single_Test=5;
    Single_Test=log(Single_Test);
    vector<double> vOld(m+1,0),vNew(m+1,0),vOld_2(m/2+1,0),vNew_2(m/2+1,0),vOld_3(m/4+1,0),vNew_3(m/4+1,0);
    vector<double> single_point_1(n+1), single_point_2(n+1),single_point_3(n+1),single_point(n+1), conservation_1(n+1),conservation_2(n+1),conservation_3(n+1);
    vector<double> X_1(m+1),X_2(m/2+1),X_3(m/4+1);
    double dt=Calculation(m,n,s,r,T,t_0, kappa, theta, gamma, sigma, x0, Single_Test, vOld, vNew, X_1,single_point_1,conservation_1);
    auto finish_1 = std::chrono::steady_clock::now();
    auto elapsed_1 = std::chrono::duration_cast<std::chrono::duration<double> >(finish_1 - start);
    dt=Calculation(m/2,n,s,r,T,t_0, kappa, theta, gamma, sigma, x0, Single_Test, vOld_2, vNew_2, X_2,single_point_2,conservation_2);
    for (int i=0;i<=n;i++){
        single_point[i]=(single_point_1[i]*4-single_point_2[i])/3;
    }
    auto finish = std::chrono::steady_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::duration<double> >(finish - start);
    dt=Calculation(m/4,n,s,r,T,t_0, kappa, theta, gamma, sigma, x0, Single_Test, vOld_3, vNew_3, X_3,single_point_3,conservation_3);
    cout<<"n= "<<n<<" , "<<"dt= "<<(T-t_0)/n<<endl;
    cout<<"m= "<<m<<" , "<<"dx= "<<(s-r)/m<<endl;
    cout<<"u(x*,T)= "<<setprecision(10)<<single_point_1[n]<<"\n";
    cout<<"Single time= "<<elapsed_1.count()<<endl;
    cout<<"u_e(x*,T)= "<<setprecision(10)<<single_point[n]<<"\n";
    cout<<"Total time= "<<elapsed.count()<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[n]-single_point_3[n])/(single_point_1[n]-single_point_2[n])<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[3200]-single_point_3[3200])/(single_point_1[3200]-single_point_2[3200])<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[400]-single_point_3[400])/(single_point_1[400]-single_point_2[400])<<endl;
    cout<<"R= "<<setprecision(6)<<(single_point_2[2]-single_point_3[2])/(single_point_1[2]-single_point_2[2])<<endl;

    ofstream output;
    output.open("/Users/sxd/Downloads/FDM_2.csv");
    if (!output.is_open())
    {
        cout<<"File not opened \n";
        throw;
    }
    // output the estimated u(x,T)
    vector<double> uNew(m+1,0);
    for (int j=0;j<=m;j++)
    {
        //Transform v(X,T) to u(x,T)
        uNew[j]=exp(-2*gamma*X_1[j])*vNew[j];
        output <<exp(X_1[j]) <<"," << uNew[j]<<endl;
    }
    cout<<"The final conservation is "<<conservation_1[n]<<endl;
    output.close();
    output.open("/Users/sxd/Downloads/FDM_conservation_2.csv");
    for (int i=0;i<=n;i++) output<<i*dt<<","<<conservation_1[i]<<"\n";
    output.close();
    output.open("/Users/sxd/Downloads/FDM_single_2.csv");
    for (int i=0;i<=n;i++) output<<i*dt<<","<<setprecision(10)<<single_point[i]<<"\n";
    output.close();
}




