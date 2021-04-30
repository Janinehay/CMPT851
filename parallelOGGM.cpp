#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <omp.h>
#include <numeric>
#include <fstream>
#include <chrono>

int main(){
    //Define constants and variables
    int Tmelt = -1;
    double Tliquid = 2.0;
    int i;
    int y;
    int e;
    int Months=12;
    int Years = 1;
    //double mu_star=0;
    //double micalc;
    //double epsilon;
    
    auto start = std::chrono::high_resolution_clock::now();
    
    //import data
    std::ifstream massbalance("massbalance.dat");
    std::vector<double>MB;
    double MBdata;
    while (massbalance >> MBdata)
    {
    	MB.push_back(MBdata);
    }
    massbalance.close();
    
    std::ifstream Temp("Temp.dat");
    double Ti[Years][Months];
    for(y=0; y<Years; y++)
    {
    	for(i=0; i<Months; i++)
    	{
    		Temp >> Ti[y][i];
    	}
    }
    
    std::ifstream Precip("Precip.dat");
    double Pi[Years][Months];
    for (y=0; y<Years; y++)
    {
    	for(i=0; i<Months; i++)
    	{
    		Precip >> Pi[y][i];
    	}
    }
    
    //Allocate Memory
    
    double* AnnualMB = (double*) malloc(Years*sizeof(double));
    double* mu_star= (double*) malloc(Years*sizeof(double));
    double* epsilon= (double*) malloc(Years*sizeof(double));
    double Pisolid [Years][Months];
    
    
    //Run over multiple years
    #pragma omp parallel for shared (Pisolid, Ti) private(y, i)
        for(y=0; y<Years; y++){
    		
    		//Calculate Solid Precipitation
    		for(i=0; i<Months; i++){
    			if(Ti[y][i]<0){
    			Pisolid[y][i]=Pi[y][i];
    		}
    		else if (Ti[y][i]>=0 && Ti[y][i] <= Tliquid){
    			double Percentsolid=(0-1)/(Tliquid-0)*Ti[y][i]+1;
    			Pisolid[y][i]=Pi[y][i]*Percentsolid;
    		}
    		else{
    			Pisolid[y][i]=0;
    		}
    		}
		}
    //Calculate mustar
	#pragma omp parallel for shared (Pisolid, mu_star, Ti) private(y, i)
    	for(y=0; y<Years; y++){
    		std::vector<double> mu;		
		for(i=0; i<Months; i++){
    			double T = Ti[y][i]-Tmelt;
    			double maxT = std::max(T,0.0);
    			double mustar;
    		
			if(maxT != 0.0){
    				mustar = Pisolid[y][i]/maxT;
    				mu.push_back(mustar);
    			}
    		
			else{
    				mustar = 1;
    				mu.push_back(mustar);
    			}
   			}
   			mu_star[y] = accumulate(mu.begin(), mu.end(), 0);	
		}
    	//Calculate epsilon
	#pragma omp parallel for shared (Pisolid, mu_star, Ti, epsilon) private(y, i)
    	for(y=0; y<Years; y++){	
			double micalc=0;
    		for(i=0; i<Months; i++){
    			double T = Ti[y][i]-Tmelt;
    			double maxT = std::max(T, 0.0);
    			micalc += Pisolid[y][i]-mu_star[y]*maxT;	
    		}
    
    		epsilon[y] = MB[y]-micalc;
		}
    		//Calculate mass Balance
    #pragma omp parallel for shared (AnnualMB, Pisolid, mu_star, Ti, epsilon) private(y, i)
    	for(y=0; y<Years; y++){	
		double sum=0;
		double mi[Months];	
    		for(i=0; i<Months; i++){
    			double T = Ti[y][i] -Tmelt;
    			double maxT = std::max(T, 0.0);
    			double mihere = Pisolid[y][i] - mu_star[y]*maxT+epsilon[y];
    			mi[i]=mihere;
    			sum+=mihere;
    		}
    		AnnualMB[y]=sum;
    	}
    	
    	
    
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
    
    std::cout << duration.count() <<std::endl;
    
    std::ofstream AnnualMBout;
    AnnualMBout.open("AnnualMBout.dat");
    for(y=0; y<Years; y++){
	AnnualMBout << AnnualMB[y]<< std::endl;
	}
    	
    return 0;
    
}
