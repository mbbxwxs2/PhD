//
//  main.cpp
//  Fokker_Planck
//
//  Created by 宋晓东 on 30/06/2020.
//  Copyright © 2020 Xiaodong Song. All rights reserved.
//
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <algorithm>
#include <string.h>
#include <cstdio>
using namespace std;

//Check the conservation law
double monitor(int n,vector<double> uNew)
{
    double result=0;
    for(int j=0;j<=n;j++) result+=uNew[j];
    return result;
}
//-------------------------------------------------------------------------------------------------------------------//
//PDF of normal distribution
double lognormal_PDF(double f, double mu, double sigma)
{
    if (f>0)
        return(exp(-0.5*pow((log(f)-mu)/sigma,2))/f/sigma/sqrt(2*M_PI));
    else
        return(0.);
}
//Probability of positive jump (logistic)
double logistic(double f,double beta_0,double beta_1)
{
    return(exp(beta_0+beta_1*f)/(1+exp(beta_0+beta_1*f)));
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
double Finite_difference(double theta, double mu, double sigma, double T, double lambda, double mu_pos, double sigma_pos, double mu_neg,double sigma_neg, double beta_0, double beta_1, int m, int n, double s, const vector<double> &X, vector<double> &u,int regime, int part)
{
    vector<double> conservation(m+1);
    vector<double> uOld(n+1),uNew(n+1);
    for(int j=0;j<=n;j++)
    {
        uOld[j] = u[j];
        uNew[j] = u[j];
    }
    double dx=s/(n+1);
    double dt=T/m;
    conservation[0]=1-monitor(n,uNew)*dx;
    //Declare nd initialize A, where A*U=D
    vector<vector<double>> A(n+1,vector<double> (n+1));
    vector<double> integral_pos(n+1),integral_neg(n+1);
    //The jump size f_j(y) (truncated lognormal distribution)
    //Denominator part (integral)
    vector<double> P(n+1,0);
    if (lambda!=0){
        for (int j1=0;j1<=n;j1++)
        {
            integral_pos[j1]=0;
            integral_neg[j1]=0;
            for (int j2=0;j2<=j1;j2++) integral_neg[j1]+=lognormal_PDF(X[j1]-X[j2],mu_neg,sigma_neg);
            for (int j2=j1;j2<=n;j2++) integral_pos[j1]+=lognormal_PDF(X[j2]-X[j1],mu_pos,sigma_pos);
            integral_pos[j1]=integral_pos[j1]*dx;
            integral_neg[j1]=integral_neg[j1]*dx;
        }
        //Probability of positive jump for X_j
        if (mu_pos!=-100 && sigma_pos!=-100 && mu_neg==-100 && sigma_neg==-100){
            for (int j=0;j<=n;j++)
            {
                P[j]=1;
            }
        } else if (mu_pos==-100 && sigma_pos==-100 && mu_neg!=-100 && sigma_neg!=-100) {
            for (int j=0;j<=n;j++)
            {
                P[j]=0;
            }
        } else {
            P[0]=1;
            for (int j=1;j<n;j++)
            {
                P[j]=logistic(X[j], beta_0, beta_1);
            }
            P[n]=0;
        }
    }
    //Construct the matrix A
    for (int j1=0;j1<=n;j1++)
    {
        A[j1][j1]=0;
        //Jump part
        if (lambda!=0){
            for (int j2=0;j2<=j1;j2++)
                if (P[j2]!=0) A[j1][j2]=-0.5*lambda*dx*dt*P[j2]*lognormal_PDF(X[j1]-X[j2],mu_pos,sigma_pos)/integral_pos[j2];
            for (int j2=j1;j2<=n;j2++)
                if (P[j2]!=1) A[j1][j2]=-0.5*lambda*dx*dt*(1-P[j2])*lognormal_PDF(X[j2]-X[j1],mu_neg,sigma_neg)/integral_neg[j2];
        }
    }
    //Diffusion part
    if (regime<part)
    {
        A[0][0]+=1+0.5*lambda*dt+dt/2/dx*(0.5*theta*(mu-X[0]-0.5*dx)+pow(sigma,2)/2/dx);
        A[0][1]+=dt/2/dx*(0.5*theta*(mu-X[0]-0.5*dx)-pow(sigma,2)/2/dx);
        for (int j1=1;j1<n;j1++)
        {
            A[j1][j1-1]+=-dt/2/dx*(pow(sigma,2)/2/dx+0.5*theta*(mu-X[j1]+0.5*dx));
            A[j1][j1]+=1+0.5*lambda*dt+dt/2/dx*(pow(sigma,2)/dx+0.5*theta*(mu-X[j1]-0.5*dx)-0.5*theta*(mu-X[j1]+0.5*dx));
            A[j1][j1+1]+=dt/2/dx*(0.5*theta*(mu-X[j1]-0.5*dx)-pow(sigma,2)/2/dx);
        }
        A[n][n-1]+=-dt/2/dx*(0.5*theta*(mu-X[n]+0.5*dx)+pow(sigma,2)/2/dx);
        A[n][n]+=1+0.5*lambda*dt-dt/2/dx*(0.5*theta*(mu-X[n]+0.5*dx)-pow(sigma,2)/2/dx);
    } else {
        A[0][0]+=1+0.5*lambda*dt+dt/2/dx*(0.5*mu+pow(sigma,2)/2/dx);
        A[0][1]+=dt/2/dx*(0.5*mu-pow(sigma,2)/2/dx);
        for (int j1=1;j1<n;j1++)
        {
            A[j1][j1-1]+=-dt/2/dx*(pow(sigma,2)/2/dx+0.5*mu);
            A[j1][j1]+=1+0.5*lambda*dt+dt/2/dx*pow(sigma,2)/dx;
            A[j1][j1+1]+=dt/2/dx*(0.5*mu-pow(sigma,2)/2/dx);
        }
        A[n][n-1]+=-dt/2/dx*(0.5*mu+pow(sigma,2)/2/dx);
        A[n][n]+=1+0.5*lambda*dt-dt/2/dx*(0.5*mu-pow(sigma,2)/2/dx);
    }
    vector<vector<double>> lower(n+1,vector<double> (n+1)), upper(n+1,vector<double> (n+1));
    for (int j1=0;j1<=n;j1++)
        for (int j2=0;j2<=n;j2++)
        {
            lower[j1][j2]=0;
            upper[j1][j2]=0;
        }
    LU_find(A, lower, upper, n);
    // start looping through time levels
    for(int i=1;i<=m;i++)
    {
        // declare vectors for matrix equations
        vector<double> d(n+1);
        //Jump part
        for (int j=0;j<=n;j++)
        {
            d[j]=0;
            if (lambda!=0){
                for (int j2=0;j2<=j;j2++)
                    if (P[j2]!=0) d[j]+=0.5*lambda*dx*dt*P[j2]*lognormal_PDF(X[j]-X[j2], mu_pos,sigma_pos)/integral_pos[j2]*uOld[j2];
                for (int j2=j;j2<=n;j2++)
                    if (P[j2]!=1) d[j]+=0.5*lambda*dx*dt*(1-P[j2])*lognormal_PDF(X[j2]-X[j],mu_neg,sigma_neg)/integral_neg[j2]*uOld[j2];
            }
        }
        //Diffusion part
        if (regime<part)
        {
            d[0]+=dt/2/dx*(-0.5*theta*(mu-X[0]-0.5*dx)+pow(sigma,2)/2/dx)*uOld[1]+(1-0.5*lambda*dt-dt/2/dx*(0.5*theta*(mu-X[0]-0.5*dx)+pow(sigma,2)/2/dx))*uOld[0];
            for(int j=1;j<n;j++)
            {
                double para_a=dt/2/dx*(-0.5*theta*(mu-X[j]-0.5*dx)+pow(sigma,2)/2/dx);
                double para_b=1-0.5*lambda*dt-dt/2/dx*(pow(sigma,2)/dx+0.5*theta*(mu-X[j]-0.5*dx)-0.5*theta*(mu-X[j]+0.5*dx));
                double para_c=dt/2/dx*(pow(sigma,2)/2/dx+0.5*theta*(mu-X[j]+0.5*dx));
                d[j]+=para_a*uOld[j+1]+para_b*uOld[j]+para_c*uOld[j-1];
            }
            d[n]+=dt/2/dx*(0.5*theta*(mu-X[n]+0.5*dx)+pow(sigma,2)/2/dx)*uOld[n-1]+(1-0.5*lambda*dt+dt/2/dx*(0.5*theta*(mu-X[n]+0.5*dx)-pow(sigma,2)/2/dx))*uOld[n];
        } else {
            d[0]+=dt/2/dx*(-0.5*mu+pow(sigma,2)/2/dx)*uOld[1]+(1-0.5*lambda*dt-dt/2/dx*(0.5*mu+pow(sigma,2)/2/dx))*uOld[0];
            for(int j=1;j<n;j++)
            {
                double para_a=dt/2/dx*(-0.5*mu+pow(sigma,2)/2/dx);
                double para_b=1-0.5*lambda*dt-dt/2/dx*pow(sigma,2)/dx;
                double para_c=dt/2/dx*(pow(sigma,2)/2/dx+0.5*mu);
                d[j]+=para_a*uOld[j+1]+para_b*uOld[j]+para_c*uOld[j-1];
            }
            d[n]+=dt/2/dx*(0.5*mu+pow(sigma,2)/2/dx)*uOld[n-1]+(1-0.5*lambda*dt+dt/2/dx*(0.5*mu-pow(sigma,2)/2/dx))*uOld[n];
        }
        //Solve X from LUu=d
        uNew = LU_solver(lower,upper,d,n);
        
        //Check the conservation law
        conservation[i]=1-monitor(n,uNew)*dx;
        //cout<<conservation[i]<<endl;
        //Set old=new
        uOld = uNew;
    }
    for (int i=0;i<=n;i++)
    {
        u[i]=uNew[i];
    }
    return (conservation[m]);
}
//-------------------------------------------------------------------------------------------------------------------//

int main()
{
    // declare and initialise parameters
    vector<double> lambda,mu,theta,sigma,mu_pos,sigma_pos,mu_neg,sigma_neg;
    ifstream file_para("/Users/sxd/Downloads/parameters.csv");
    if (file_para.is_open())
    {
        double l1,l2,l3,l4,l5,l6,l7,l8;
        char delimiter;
        while (file_para>>l1>>delimiter>>l2>>delimiter>>l3>>delimiter>>l4>>delimiter>>l5>>delimiter>>l6>>delimiter>>l7>>delimiter>>l8)
        {
            lambda.push_back(l1);
            mu.push_back(l2);
            theta.push_back(l3);
            sigma.push_back(l4);
            mu_pos.push_back(l5);
            sigma_pos.push_back(l6);
            mu_neg.push_back(l7);
            sigma_neg.push_back(l8);
        }
        file_para.close();
    }
    vector<int> chain_number;
    ifstream file_chain("/Users/sxd/Downloads/chain.csv");
    if (file_chain.is_open())
    {
        int l1;
        while (file_chain>>l1)
        {
            chain_number.push_back(l1-1);
        }
        file_chain.close();
    }
//-------------------------------------------------------------------------------------------------------------------//
// declare and initialise grid paramaters
// Update:
//    part:         Total number of regimes 2*part
//    length_day:   Length of the CI
//    x0:           The first CI value
//    T:            The total time of a regime
//    beta_0:       Para_1 in Logistic model of jump sign
//    beta_1:       Para_2 in Logistic model of jump sign
//    n:            The grid number of CI (dx)
//    Period:       The total period in the day
//    s:            Maximum value of CI
//    day_number:   the number of day in the whole data set
//-------------------------------------------------------------------------------------------------------------------//
    int part=4;
    double T=1, x0=0.0671642;
    int length_day=619;
    int n=1599,period=16,TH=20;
    int day_number=1;
    int slice=round(length_day/(double) period);
    int slice_end=length_day-slice*(period-1);
    cout<<"slice="<<slice<<endl;
    cout<<"end="<<slice_end<<endl;
    double beta_0=4.406,beta_1=-7.533;
    double s=1.5;
    double dx=s/(n+1);
    cout<<"dt= "<<double(T/slice/TH)<<endl;
    cout<<"dx= "<<dx<<endl;
    vector<int> day_chain(period);
    cout<<"Regime number= ";
    for (int i=0;i<period;i++){
        day_chain[i]=chain_number[(day_number-1)*period+i];
        cout<<day_chain[i]<<" ";
    }
    cout<<'\n';
    int m;
    vector<double> conservation_total(period+1);
    vector<double> X(n+1);
    // setup and initialise x
    for(int j=0;j<=n;j++)
    {
        X[j] = (0.5+j)*dx;
    }
    //Declare nd initialize u(dx,dt)
    vector<double> u(n+1);
    for(int j=0;j<=n;j++) u[j]=0;
    int start=round(x0/dx);
    cout<<"X0= "<<start*dx<<endl;
    u[start]=1/dx;
    conservation_total[0]=0;
    for (int time=0;time<period;time++)
    {
        if (time<period-1) m=slice*TH; else m=slice_end*TH;
        double energy=Finite_difference(theta[day_chain[time]],mu[day_chain[time]],sigma[day_chain[time]],T,lambda[day_chain[time]],mu_pos[day_chain[time]],sigma_pos[day_chain[time]],mu_neg[day_chain[time]],sigma_neg[day_chain[time]],beta_0,beta_1,m,n,s,X,u,day_chain[time],part);
        conservation_total[time+1]=conservation_total[time]+energy;
        char data_name[200];
        sprintf(data_name, "%s%i%s", "/Users/sxd/Downloads/FP_log_period_", time,".csv");
        ofstream output(data_name);
        if (!output.is_open())
        {
            cout << "error" << endl;
            return 0;
        }
        // output the estimated u(x,T)
        for (int j=0;j<=n;j++)
        {
            output<<X[j] <<"," << u[j]<<"\n";
        }
        output.close();
        cout<<"Period "<<time+1<<" energy loss= "<<energy<<endl;
        //cout<<"Period "<<time+1<<" energy loss= "<<conservation_total[time+1]<<endl;
    }
    //cout<<"The final conservation is "<<conservation[m]<<endl;
    //output.close();
    ofstream output_energy("/Users/sxd/Downloads/FP_log_period_conservation.csv");
    if (!output_energy.is_open())
    {
        cout << "error" << endl;
        return 0;
    }
    for (int i=0;i<=period;i++) output_energy<<i<<","<<conservation_total[i]<<endl;
    output_energy.close();
    cout<<"Energy losses= "<<conservation_total[period]<<endl;
}

