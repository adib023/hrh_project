#ifndef DEFINITIONS_H_
#define DEFINITIONS_H_

// C libraries
#include <float.h>
#include <stdio.h>
#include <stdlib.h>

// C++ libraries
#include <cmath>

#include <cstdlib> 
using std::size_t ;

#include <string>
using std::string;

#include <iomanip>
using std::setprecision ;

#include <iostream>
using std::cout;
using std::cerr;
using std::endl;

#include <fstream>
using std::stringstream ;
using std::ifstream ; 
using std::getline ; 
using std::istringstream ; 
using std::ofstream ; 

#include <vector>
using std::vector ; 

#include <Eigen/Dense>
using Eigen::MatrixXd ; 
using Eigen::Matrix3d ; 
using Eigen::Vector2d ; 
using Eigen::Vector3d ; 
using Eigen::Vector4d ; 
using Eigen::VectorXd ;
using Eigen::VectorXi ;
using Eigen::Vector3i ; 
using Eigen::Vector4i ; 
using Eigen::Matrix3Xd ; 

#define space "\t"
#define comma ","
#define PI 3.141592653589793

#define BASIS_TYPE 1 // Bloch mode basis
//#define BASIS_TYPE 2 // polynomial basis

#define TRANSIENT 1

#endif
