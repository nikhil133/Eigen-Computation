#pragma once
#include<iostream>
using namespace std;
#include "geteigen.c"
#include "geteigen_emxAPI.c"
#include "geteigen_emxutil.c"
#include "geteigen_initialize.c"
#include "geteigen_rtwutil.c"
#include "geteigen_terminate.c"
#include "geteigen_types.h"
#include "rt_nonfinite.c"
#include "rtGetInf.c"
#include "rtGetNaN.c"
#include "rtwtypes.h"
#include "sqrt.c"
#include "eignwcpp.h"


eignwcpp::eignwcpp(double **A,double **egvctr,double **eigval,int m,int n)
{
}
void eignwcpp::eigen_mem_alloc(double **A,double **egvctr,double **eigval,int m,int n)
{
	
	emxArray_real32_T *x;
	emxArray_real32_T *y;
	emxArray_real32_T *z;

	emxInit_real32_T(&x,2);
	emxInit_real32_T(&y,2);
	emxInit_real32_T(&z,2);

	x->size[0]=m;
	x->size[1]=n;
	y->size[0]=m;
	y->size[1]=n;
	z->size[0]=m;
	z->size[1]=n;

	emxEnsureCapacity((emxArray__common *)x,0,(int32_T)sizeof(real_T));
	emxEnsureCapacity((emxArray__common *)y,0,(int32_T)sizeof(real_T));
	emxEnsureCapacity((emxArray__common *)z,0,(int32_T)sizeof(real_T));
	
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			x->data[i*n+j]=A[i][j];
		}
	}
	geteigen(x,y,z);
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			egvctr[i][j]=y->data[i*n+j];
		}
	}
	for(int i=0;i<m;i++)
	{
		for(int j=0;j<n;j++)
		{
			eigval[i][j]=z->data[i*n+j];
		}
	}
}
