#ifndef MATPROPLE_H_
#define MATPROPLE_H_

#include "definitions.h"

class MatPropLE
{
public:
	void preCompute()
	{
		_Cmat.resize(6,6) ; 

		_Cmat.setZero();

		double mu	= _Emat / (2*(1 + _nu));
		double lambda	= _Emat * _nu / ((1+_nu)*(1-2*_nu)) ;
		
		_Cmat(0,0) = lambda + 2*mu ; 
		_Cmat(1,1) = lambda + 2*mu ; 
		_Cmat(2,2) = lambda + 2*mu ; 

		_Cmat(3,3) = mu ;   
		_Cmat(4,4) = mu ;   
		_Cmat(5,5) = mu ;   

		_Cmat(0,1) = lambda ;
		_Cmat(0,2) = lambda ;
		_Cmat(1,0) = lambda ;
		_Cmat(2,0) = lambda ;
		_Cmat(2,1) = lambda ;
		_Cmat(1,2) = lambda ;
		
		_Cs	=	sqrt(mu/_rho);
		_Cp	=	sqrt((lambda+2*mu)/_rho);
	}
	
	void readMaterialProperties(ifstream & fid)
	{
		size_t dumy ; 
		fid >> dumy >> _Emat >> _nu >> _rho ; 
	}

	const MatrixXd & getCmat()
	{
		return _Cmat ;
	}

	const double & getDensity()
	{
		return _rho ;
	}

	void getBndryData_solid(double &density, MatrixXd & Cmat, double & Cp, 
				double &Cs)
	{
		density = _rho ; 
		Cmat	= _Cmat ; 
		Cs	= _Cs ; 
		Cp	= _Cp ;
	}

private:
	double		_Emat ; // Youngs modulus
	double		_nu ;   // Poisson ratio 
	double		_rho ;  // density 

	MatrixXd	_Cmat ;
	double		_Cs ; // shear wave speed
	double		_Cp ; // pressure wave speed
};

#endif
