#include<iostream>
using namespace std;
#include"eignwcpp.h"

int main()
{
	//eigenfilter e=eigenfilter();
	double **a,**egvctor,**egvalu,**c;
	int m,n;
	cout<<"Enter order of the matrix\n"
	cout<<"row= ";
	cin>>m;
	cout<<"col= ";
	cin>>n;
	a=new double*[m];
	egvctor=new double*[m];
	egvalu=new double*[m];
	cout<<"Enter the user data whose eigen value and eigen vector is to be created";
	for(int i=0;i<m;i++)
	{
		a[i]=new double[n];
		egvctor[i]=new double[n];
		egvalu[i]=new double[n];
		for(int j=0;j<n;j++)
		{
			cin>>a[i][j];
			
		}
	}
	cout<<"*****Matrix A*****\n";
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			cout<<" "<<a[i][j];
			
		}
		cout<<"\n";
	}
	
	
	eignwcpp e=eignwcpp(a,egvctor,egvalu,m,n); 
	e.eigen_mem_alloc(a,egvctor,egvalu,m,n);
	cout<<"\n\n*****Eigen Vector of Matrix A*****\n";
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			cout<<" "<<egvctor[i][j];
			
		}
		cout<<"\n";
	}
	
	cout<<"\n\n*****Eigen Vector of Matrix A*****\n";
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			cout<<" "<<egvalu[i][j];
			
		}
		cout<<"\n";
	}
	
	
	
	return 0;
}