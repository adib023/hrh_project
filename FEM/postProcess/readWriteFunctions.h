#ifndef READWRITEFUNCTIONS_H_
#define READWRITEFUNCTIONS_H_

#include "definitions.h"



extern std::string readString(ifstream & inputFile)
{
	string line ;

	getline(inputFile,line) ; 
	while ((line[0] == '%')||(line.length() == 0))
		getline(inputFile,line) ; 

	return line ; 
}


extern double readDouble(ifstream & inputFile)
{
	string line ;

	getline(inputFile,line) ; 
	while ((line[0] == '%')||(line.length() == 0))
		getline(inputFile,line) ; 

	return atof(line.c_str()) ; 
}

extern int readInteger(ifstream & inputFile)
{
	string line ;

	getline(inputFile,line) ; 
	while ((line[0] == '%')||(line.length() == 0))
		getline(inputFile,line) ; 

	return atoi(line.c_str()) ; 
}

extern Vector3d readVector(ifstream & inputFile)
{
	string line ;

	getline(inputFile,line) ; 
	while ((line[0] == '%')||(line.length() == 0))
		getline(inputFile,line) ; 

	Vector3d vec ; 
	istringstream iss ; 
	string partial;
	
	iss.str(line) ;

	for (int p = 0; p<3 ; ++p)
	{
		getline(iss,partial,' ');
		vec(p) = atof(partial.c_str());
	}

	iss.clear();

	return vec ;
}

extern Vector3i readIntVector(ifstream & inputFile)
{
	string line;

	getline(inputFile,line) ; 
	while ((line[0] == '%')||(line.length() == 0))
		getline(inputFile,line) ; 

	Vector3i vec ; 
	istringstream iss ; 
	string partial;
	
	iss.str(line) ;

	for (int p = 0; p<3 ; ++p)
	{
		getline(iss,partial,' ');
		vec(p) = atoi(partial.c_str());
	}

	iss.clear();

	return vec ;
}

#endif
