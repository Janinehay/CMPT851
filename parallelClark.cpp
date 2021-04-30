#include <stdio.h>
#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <omp.h>
#include <chrono>
#include <random>

int main(){
	
	auto start = std::chrono::high_resolution_clock::now();
	
	//constants
	double fm = 0.83881;
	double asnow = 0.01022;
	int I = 1362;
	int ti = 30;
	double sigmai=2;
	int month = 12;
	int year = 1;
	int y;
	int i;
	int j;
	int sigma = 2;
	double c = 2.01;
	double L = 3.35;
	int d = 2;
	double Tmp = 0;
	double pw = 997;
	double p = 917;
	
	
	//declare
	//std::vector<double> D;
	double D[year][month];
	
	//import data
	std::ifstream Temp("Temp.dat");
	std::ifstream Precip("Precip.dat");
	double Ti[year][month];
	double Pi[year][month];
	
	for(y=0; y<year; y++)
	{
		for(i=0; i<month; i++)
		{
			Temp >> Ti[y][i];
			Precip >> Pi[y][i];
		}
	}
	
	
	//Calculate Positive degree days
	for(y=0; y<year; y++)
	{

		for(i=0; i<month; i++)
		{
			
			double Dmonth=0;
			std::default_random_engine generator;
			std::normal_distribution<double> distribution(Ti[y][i], sigma);
			for(j=0; j<ti; j++)
			{
				double number = distribution(generator);
				Dmonth +=number;
			}
			D[y][i] = Dmonth;	
		}
		
	}
	
	//Calculate degreeday factor
	double fdd = fm+asnow*I;
	
	
	//Calculate Ablation
	//Allocate Memory
	double M[year][month];
	
	#pragma omp parallel for shared(M, D, fdd) private(i)
		for(y=0; y<year; y++)
		{
			for(i=0; i<month; i++)
			{
				M[y][i] = fdd*D[y][i];
			}		
		}

	
	//Calculate Refreeze
	
	//Calculate Ta
	std::vector<double> Ta;
	for(y=0; y<year; y++)
	{
		double annualT = 0;
		
		for(i=0; i<month; i++)
		{
			annualT+=Ti[y][i];
		}
		
		Ta.push_back(annualT/12);	
	}
	
	//calculate h
	double R[year];
	#pragma omp parallel for shared(R) private(i)
		for(y=0; y<year; y++)
		{
			double h=0;
			for(i=0; i<month; i++)
			{
				if(Ti[y][i]<0){
					h+=Pi[y][i];
				}
				else h+=0;
			}
			h = h/12;
			double maxt = std::max(Tmp-Ta[y], 0.0);
			double r = c/L *(d/h)*maxt;
			double eta = std::min(r, 1.0);	
		
			R[y]=eta*h;
		}
	
	//Calculate Accumulation
	double fsnow[year][month];
	
	//calculate fsnow
	for (y=0; y<year; y++)
	{
		for(i=0; i<month; i++)
		{
			if(Ti[y][i]<0)
			{
				fsnow[y][i]=1;
			}
			else if(Ti[y][i]>=0 && Ti[y][i] <= 2)
			{
				double percentsolid = (0-1)/(2-0)*Ti[y][i]+1;
				fsnow[y][i] = percentsolid;
			}
			else
			{
				fsnow[y][i] = 0;
			}
		}
	}
	
	double acc[year][month];
	#pragma omp parallel for shared(acc, fsnow) private(y, i)
		for(y=0; y<year; y++)
		{
			for(i=0; i<month; i++)
			{
				acc[y][i] = pw/p *fsnow[y][i]* Pi[y][i];
			}
		}
	
	//Calculate annual Mass balance
	double MB[year];
	for(y=0; y<year; y++)
	{
		double My = 0;
		double accy = 0;
		for(i=0; i<month; i++)
		{
			My+=M[y][i];
			accy+=acc[y][i];
		}
		
		MB[y] = My+accy+R[y];
	}
	
	auto stop = std::chrono::high_resolution_clock::now();
	
	auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop-start);
	
	std::cout<< "Time take by function: "<<duration.count() << " microseconds"<< std::endl;
	
	std::ofstream MBout;
	MBout.open("MBout.dat");
	for(y=0; y<year; y++)
	{
		MBout << MB[y] <<std::endl;
	}	
	
	
	return 0;
}
