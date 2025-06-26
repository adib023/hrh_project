#ifndef VTK_WRITE_H_
#define VTK_WRITE_H_

#include "definitions.h"
#include "ArchSolid3D.h"

extern void writeVTK_disp(string & fname,MatrixXd & dval, 
		ArchSolid3D& archSolid)
{
		ofstream fid(fname.c_str());

		fid << "# vtk DataFile Version 2.0\n" << "mode shapes\n" ;
		fid << "ASCII\n" << "DATASET UNSTRUCTURED_GRID\n";

		fid << "POINTS " << archSolid.getNumNodes() << " float\n";

		size_t nNode = archSolid.getNumNodes() ; 

		vector<Vector3d> xc  ; 
		xc = archSolid.getNodalCoords() ; 

		double dfactr = 1.0 ;

		vector<Vector3d> xd ; xd.resize(nNode); 

		/*for (size_t p = 0 ; p < nNode ; ++p)
			for (size_t q = 0 ; q < 3 ; ++q)
				xd[p](q) = xc[p](q) + dfactr * dval(p,q) ;*/ 

		cout << "nNode " << nNode<< endl; 
		for (size_t p = 0 ; p < nNode ; ++p) 
			fid << xc[p](0) + dval(p,0) << " " << xc[p](1) + dval(p,1) << " " << xc[p](2) + dval(p,2) << "\n" ; 

		size_t ntets = archSolid.getNumElem() ; 

		vector<VectorXi> tetID = archSolid.getElem() ; 

		fid << "CELLS " << ntets << " " << 28 * ntets << "\n" ;

		for (size_t p = 0; p < ntets; ++p) {
			/*fid << "27 " < tetID[p](0) << " " << tetID[p](1) <<
				" " << tetID[p](2) << " " << tetID[p](3) <<
				"\n" ;*/

			fid << "27";
			for (size_t q = 0; q < 27; q++) {
				fid << " " << tetID[p](q);
			}
			fid << "\n";
		}
		fid << "\n" << "CELL_TYPES " << ntets << "\n" ;
	
		for (size_t p = 0 ; p < ntets ; ++p)
			fid << "29\n" ;

		fid << "\n\n" << "POINT_DATA " << nNode << "\n" ;
		fid << "SCALARS scalars float" << "\n" ; 
		fid << "LOOKUP_TABLE default" << "\n" ; 


		for (size_t p = 0 ; p < nNode ; ++p)
			fid << 1.0e2*sqrt( dval(p,0)*dval(p,0) + 
				 dval(p,1)*dval(p,1) + 
				 dval(p,2)*dval(p,2))  << "\n" ; 

		fid.close();
}

extern void writeVTK_disp2(string & fname,MatrixXd & dval, 
		ArchSolid3D& archSolid, double dfactr)
{
		ofstream fid(fname.c_str());

		fid << "# vtk DataFile Version 2.0\n" << "mode shapes\n" ;
		fid << "ASCII\n" << "DATASET UNSTRUCTURED_GRID\n";

		fid << "POINTS " << archSolid.getNumNodes() << " float\n";

		size_t nNode = archSolid.getNumNodes() ; 

		vector<Vector3d> xc  ; 
		xc = archSolid.getNodalCoords() ; 


		vector<Vector3d> xd ; xd.resize(nNode); 

		for (size_t p = 0 ; p < nNode ; ++p)
			for (size_t q = 0 ; q < 3 ; ++q)
				xd[p](q) = xc[p](q) + dfactr * dval(p,q) ; 


		for (size_t p = 0 ; p < nNode ; ++p) 
			fid << xd[p](0) << " " << xd[p](1) << " " << xd[p](2) << "\n" ; 

		size_t ntets = archSolid.getNumElem() ; 

		vector<VectorXi> tetID = archSolid.getElem() ; 

		fid << "CELLS " << ntets << " " << 11 * ntets << "\n" ;

		for (size_t p = 0 ; p < ntets ; ++p)
			fid << "10 " << tetID[p](0) << " " << tetID[p](1) << 
				" " << tetID[p](2) << " " << tetID[p](3) << 
				" " << tetID[p](4) << " " << tetID[p](7) << 
				" " << tetID[p](5) << " " << tetID[p](6) << 
				" " << tetID[p](9) << " " << tetID[p](8) << 
				"\n" ;

		fid << "\n" << "CELL_TYPES " << ntets << "\n" ;
	
		for (size_t p = 0 ; p < ntets ; ++p)
			fid << "24\n" ;

		fid << "\n\n" << "POINT_DATA " << nNode << "\n" ;
		fid << "SCALARS scalars float" << "\n" ; 
		fid << "LOOKUP_TABLE default" << "\n" ; 


		for (size_t p = 0 ; p < nNode ; ++p) //// this is place where I need to change is if I want to plot differernt components
			fid << dfactr * sqrt( dval(p,0)*dval(p,0) + 
				 dval(p,1)*dval(p,1) + 
				 dval(p,2)*dval(p,2))  << "\n" ; 

		fid.close();
}



extern void writeVTK_disp3(string & fname, MatrixXd & dval, 
		ArchSolid3D& archSolid, double dfactr)
{
		ofstream fid(fname.c_str());

		fid << "# vtk DataFile Version 2.0\n" << "mode shapes\n" ;
		fid << "ASCII\n" << "DATASET UNSTRUCTURED_GRID\n";

		fid << "POINTS " << archSolid.getNumNodes() << " float\n";

		size_t nNode = archSolid.getNumNodes() ; 

		vector<Vector3d> xc  ; 
		xc = archSolid.getNodalCoords() ; 
 		
		//MatrixXd xc;
		//archSolid.getDefXcoord(xc, defXcoordFile) ; 

		 vector<Vector3d> xd ; xd.resize(nNode); 
		
		string finalUURfileName =  "../Results/final_UUR.txt";

		ifstream finalUURfile(finalUURfileName);

		MatrixXd finalUUR; 

		archSolid.getDefXcoord(finalUUR, finalUURfile);

		//cout << finalUUR << endl;

		
		for (size_t p = 0 ; p < nNode ; ++p)
			for (size_t q = 0 ; q < 3 ; ++q)
				xd[p](q) = xc[p](q) + finalUUR(p,q) ;
	
		for (size_t p = 0 ; p < nNode ; ++p)
			for (size_t q = 0 ; q < 3 ; ++q)
				xd[p](q) = xd[p](q) + dfactr * dval(p,q) ; 


		for (size_t p = 0 ; p < nNode ; ++p) 
			fid << xd[p](0) << " " << xd[p](1) << " " << xd[p](2) << "\n" ; 

		size_t ntets = archSolid.getNumElem() ; 

		vector<VectorXi> tetID = archSolid.getElem() ; 

		fid << "CELLS " << ntets << " " << 11 * ntets << "\n" ;

		for (size_t p = 0 ; p < ntets ; ++p)
			fid << "10 " << tetID[p](0) << " " << tetID[p](1) << 
				" " << tetID[p](2) << " " << tetID[p](3) << 
				" " << tetID[p](4) << " " << tetID[p](7) << 
				" " << tetID[p](5) << " " << tetID[p](6) << 
				" " << tetID[p](9) << " " << tetID[p](8) << 
				"\n" ;

		fid << "\n" << "CELL_TYPES " << ntets << "\n" ;
	
		for (size_t p = 0 ; p < ntets ; ++p)
			fid << "24\n" ;

		fid << "\n\n" << "POINT_DATA " << nNode << "\n" ;
		fid << "SCALARS scalars float" << "\n" ; 
		fid << "LOOKUP_TABLE default" << "\n" ; 


		for (size_t p = 0 ; p < nNode ; ++p) //// this is place where I need to change is if I want to plot differernt components
			fid << dfactr * sqrt( dval(p,0)*dval(p,0) + 
				 dval(p,1)*dval(p,1) + 
				 dval(p,2)*dval(p,2))  << "\n" ; 

		fid.close();
}



#endif
