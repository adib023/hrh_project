#include "definitions.h"
#include "ArchSolid3D.h"
#include "vtkWrite.h"

int main(int argc, char * argv[])
{	
	// give config file name here: 
	string configFileName = "../hrh_new_long_fem.txt" ;
	int caseNum =  atoi(argv[1]); 

 	string datFileName, outFile;

	if(caseNum  == 1 ) {
		cout << "Here I am" << endl;
		datFileName = "../Results/output.txt";
		outFile = "../vtkDeformed/defFile_" ; 
	}
	else if(caseNum  == 2 ){
		datFileName = "../Results/vtkDefMode.txt";
		outFile = "../vtkModeDeformed/defModeFile_" ;
		
	}
	else{
		datFileName = "../Results/vtkMode.txt";
		outFile = "../vtkModeUndeformed/modeFile_" ;
	}



	// give number of steps here: 
	size_t nSteps = atoi(argv[2]) ;


	double dfactr;
	ArchSolid3D archSolid ;

	archSolid.preProcess(configFileName) ; 
	
	ifstream datFile(datFileName) ;
	//datFile>>nSteps;
	MatrixXd dval ;
	
	//string outFile = "../vtkFiles/file_" ; 

	for (size_t p = 1 ; p <= nSteps ; ++p)
	{	
		std::stringstream tstr ; tstr << p ;  
		string outT = outFile  + tstr.str() + ".vtk" ; 

		archSolid.getDisp(dval, datFile) ; 

		if (archSolid._elemType == "TET2")

			if(caseNum == 1){ 
			dfactr = 1;
			writeVTK_disp2(outT,dval,archSolid, dfactr);
			}
			
			else if	(caseNum == 2){
			dfactr = 0.001;
			writeVTK_disp3(outT,dval,archSolid,dfactr);
			}
			
			else if	(caseNum == 3){
			dfactr = 0.0005;
			writeVTK_disp2(outT,dval,archSolid,dfactr);
			}

		else
			writeVTK_disp(outT,dval,archSolid);
			
	}

	return 0;
}
