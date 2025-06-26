#ifndef ARCHSOLID3D_H_
#define ARCHSOLID3D_H_

#include "definitions.h"
#include "MatPropLE.h"
#include "readWriteFunctions.h"


class ArchSolid3D
{
public:
	ArchSolid3D()
	{
		_nFacet = 0 ;
	}

	void preProcess(string comsolMeshFile)
	{
		readUnitCellDescription(comsolMeshFile);
	}

	size_t & getNumElem()
	{
		return _nElem ;
	}


	vector<VectorXi> & getElem()
	{
		return _elem ; 
	}
// number of nodes
	size_t & getNumNodes()
	{
		return _nNode ;
	}

// nodal coords
	vector<Vector3d> & getNodalCoords()
	{
		return _xCoord ; 
	}	

	void getDisp(MatrixXd & disp,ifstream &datFile)
	{ 	
		disp.resize(_nNode,3) ; 

		disp.setZero() ;

		for (size_t p = 0 ; p < _nNode ; ++p)
			datFile >> 
				disp(p,0) >>  
				disp(p,1) >>  
				disp(p,2) ; 

	}



	void getDefXcoord(MatrixXd & def_xCoord, ifstream & defXcoordFile)
	{
	
	def_xCoord.resize(_nNode,3);
	def_xCoord.setZero();
	
	for (size_t p = 0 ; p < _nNode ; ++p)
			defXcoordFile >> 
				def_xCoord(p,0) >>  
				def_xCoord(p,1) >>  
				def_xCoord(p,2) ;


	}

	void readUnitCellDescription(string unitCellFile)
	{
		ifstream fid(unitCellFile.c_str());   
		if ((bool)fid == false)
		{
			cerr << "unit cell file:" << unitCellFile << " does not exist" << endl ; 
			exit(0);
		}
		string elemType ; 
		fid >> elemType ; 


		size_t nMat, nProp ; 
		fid >> nMat >> nProp; 
		
		_matProp.resize(nMat);

		for (size_t p = 0 ; p < nMat ; ++p)
			_matProp[p].readMaterialProperties(fid);

		fid >> _nNode ;
		
		Vector3d xc; 
		for (size_t p = 0 ; p < _nNode ; ++p)
		{
			fid >> xc(0) >> xc(1) >> xc(2) ; 
			_xCoord.push_back(xc);
		}

		size_t elemOrder ;

		fid >> _nElem >> elemOrder ;

		if (elemType == "TET1")
			_nNodeElem = 4 ; 
		else if (elemType == "TET2")
			_nNodeElem = 10 ; 

		else if (elemType == "HEX2")
			_nNodeElem = 27 ; 
		else
		{
			cerr << "invalid elem " << elemOrder << endl ;
			exit(0);
		} 

		_matIndxElem.resize(_nElem);
		_elem.resize(_nElem);

		size_t matindx1;

		for (size_t p = 0 ; p < _nElem ; ++p)
		{
			_elem[p].resize(_nNodeElem) ; 

			for (size_t q = 0 ; q < _nNodeElem ; ++q)
				fid >> _elem[p](q) ; 

			fid >> matindx1 ;
			_matIndxElem[p] = matindx1 - 1;
		}
	
		fid.close();

		_elemType = elemType ;
	}

	string 			_elemType ; 
		

private:
	// total number of elements, facets, elem and nodes in 
	// the entire lattice
	size_t			_nElem , _nFacet, _nNode ;

	size_t 			_nNodeElem ;

	vector<VectorXi>	_elem ;
	vector<Vector3d>	_xCoord ; 
	vector<Vector3d>	_def_xCoord ;

	vector<MatPropLE> 		_matProp ;
	vector<size_t>			_matIndxElem ;


	double			_dt, _finalTime ; 
	size_t			_nOutput;
	ifstream		_rawSolnFile ;

};


#endif
