#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>

#include <stdlib.h> 
#include <iostream>
#include <fstream> 
#include <algorithm> 
#include <list>
#include <map> 
#include <sstream>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Periodic_3_triangulation_traits_3<K> GT;
typedef CGAL::Periodic_3_Delaunay_triangulation_3<GT>       PDT;

typedef PDT::Cell_handle                    Cell_handle;
typedef PDT::Point                          Point;
typedef PDT::Iso_cuboid                     Iso_cuboid;
typedef PDT::Periodic_tetrahedron_iterator  Periodic_tetr_iterator;
typedef PDT::Iterator_type                  Iterator_type;   

typedef GT::Point_3          Point;
typedef GT::Tetrahedron_3    Tetrahedron;
               
class cellData {
	public:
		int      vertex[4]  ;
		int      neighbor[4];
		bool     OnBoundary ;
		int      Img[4][3]  ;
	
		// get_ID uses Cantor Pairing function three consecutive times 
		// to generate unique 1-1 ID's for a given tetrahedron, based 
		// on vertex indices
		long get_ID() {
			long ID1;
			long ID2;
		
			ID1 =  CantorPair(vertex[1],vertex[2]);
			ID2 =  CantorPair(vertex[3],vertex[4]);
			return CantorPair(      ID1,      ID2);
		}
	private:
		long CantorPair(long ii, long jj){
			return ( (ii+jj)*(ii+jj+1) ) / 2  +  jj ;
		}
}; 

void readInput(int iFileType, std::list<Point>    &L     , double      H[3]    , 
			   std::map<Point,int> &vIndex, const char*  filename);

int findNeighbor(int &ICell, int &i , int &IPoint , std::list<int>);

void  Get_Image_FromPos(Point &p, double CellDim[3], int &i, int &d);

// Global Variables
cellData*  Cell      ;
int        NCells    ;
int        NPoints   ;
int        IFrame    ;

// Variables for reading input
bool FirstEntry  = true;
std::ifstream      ifs ;
int iFr_previous = -1  ;

//*******************************************************************************
// ***This program is designed to perform 3D Periodic Delaunay Triangulations
// for input files containing many frames. The readme file contains  all the 
// necessary information about usage and design choices.
//
// argument 1 = file name and path
// argument 2 = file type (4 -> type 4 file, 11 -> lammps dump file)
//*******************************************************************************
int main(int argc, char *argv[]) {
	bool           LFound        ;
	int            iFileType     ;
	int            ICell , JCell ;
	int            IPoint, JPoint;
	int            ii, jj, kk    ;
	double         vertex[3]     ;
	double         H[3]          ;
	Point          p             ;
	Tetrahedron    tetr          ;
	
	int NSurfCells = 0;
	
	const char       *inputFile = argv[1];
	std::stringstream streamArg2(argv[2]);
	streamArg2 >> iFileType;  // 1 - type 4, 2 - lammps dump
	
	std::map<Point,int>               vIndex   ;
	std::map<int,std::list<int> >     nnTetr   ;
	std::map<int,std::list<int> >     nnVert   ;
	std::list<Point>                  DummyList;	
	
	// The cube for the periodic domain
	Iso_cuboid domain;

	// construction from a list of points :
	std::list<Point> L;
	
	//output files
	FILE* oFileT    = fopen ("output.tvt","w"); // tetr-vert-tetr
	FILE* oFileVert = fopen ("output.vv" ,"w"); // vert-vert
	FILE* oFileTetr = fopen ("output.vt" ,"w"); // vert-tetr 

	// file type to use
	// 1 - type 4, 2 - lammps dump
	iFileType     =  1;

	std::cout << inputFile << " " << iFileType;

	if ((iFileType != 4) && (iFileType != 11)) { 
		std::cout << "ERROR - Filetype not supported"       << std::endl;
		std::cout << "        > Second argument can be:"    << std::endl;
		std::cout << "                  4 - type 4 file"    << std::endl;
		std::cout << "                 11 - lammps dump"    << std::endl;
		return                                                          ; 
	}

	for(;;) {
		// Read data from file - break if reached end fo file
		if (ifs.eof()) break;
		readInput(iFileType, L, H, vIndex, inputFile);

		// fix an end-of-file bug
		if (iFr_previous == IFrame) break;
		iFr_previous = IFrame;
		
		std::cout << "Analyzing Frame #  " << IFrame << std::endl;

		// construct periodic domain
		domain = Iso_cuboid(0.0, 0.0, 0.0, H[0], H[1], H[2]);
	
		// Put the domain with the constructor
		PDT T(DummyList.begin(), DummyList.end(), domain); 
		
		// Insert Points in the triangulatiom 
		// (true -> Iterator range insertion using spatial sorting and dummy point heuristic)
		T.insert(L.begin(), L.end(), true); 
				
		NCells  =  T.number_of_cells(); 	
		Cell    =  (cellData*)malloc(NCells*sizeof(cellData));	

		// Loops over all periodic tetrahedra and keep ONLY those with a
		// circumcenter inside the original domain. Also prepare lists of
		// tetrahedra that contain a specific Point.
		ICell = 0;
		for (Periodic_tetr_iterator tit  = T.periodic_tetrahedra_begin(PDT::UNIQUE) ;
									tit != T.periodic_tetrahedra_end  (PDT::UNIQUE) ; ++tit) {
			for (int i = 0; i < 4; i++ ) {
				IPoint = vIndex [ tit -> at(i).first ];
				Cell[ICell].vertex[i] =  IPoint       ;
				nnTetr[IPoint].push_back(ICell)       ;    
			}
			
			Cell[ICell].OnBoundary  =  false;
			
			// Determine if the Cell is located on the surface of the
			// original domain
			tetr = T.construct_tetrahedron( tit->at(0).first , tit->at(1).first , 
			                                tit->at(2).first , tit->at(3).first ,
			                                tit->at(0).second, tit->at(1).second,
			                                tit->at(2).second, tit->at(3).second );
			
			// Locate Cells intersecting the boundary of the original domain
			for (int i = 0; i < 4; i++) {
				p = tetr[i];
				for (int j = 0; j < 3; j++ ) {
					Cell[ICell].Img[i][j] = 0;
					if ( (p[j] >= H[j]) || (p[j] <= 0.0) ) {
						Cell[ICell].OnBoundary = true  ;
                        Get_Image_FromPos( p, H, ii, j);
                        Cell[ICell].Img[i][j]  = ii    ;
					}
				}
			}
			ICell  +=  1;
		}
		
		
		// Locate neighbor Cells & count surface Cells
		for (ICell = 0; ICell < NCells; ICell++) {
			for (int i = 0; i < 4; i++) {
				//std::cout << "searching for neighbor:  " << i+1 
				//          << "for Cell:  " << ICell+1 << std::endl;
				if (i==3) {IPoint = Cell[ICell].vertex[0];   }
				else      {IPoint = Cell[ICell].vertex[i+1]; }
			 
				Cell[ICell].neighbor[i] = findNeighbor(ICell, i, IPoint, 
													    nnTetr[IPoint] );
			}
			if (Cell[ICell].OnBoundary) NSurfCells += 1;
		}
		
		// Map connectivity between vertices
		for (ICell = 0; ICell < NCells ; ICell++){
			for (int i = 0; i < 4; i++ ) {
				IPoint  =  Cell[ICell].vertex[i];
				for (int j = 0; j < 4; j++) {
					JPoint =  Cell[ICell].vertex[j];
					if (IPoint==JPoint) continue   ;
					LFound =  (std::find(nnVert[IPoint].begin(), nnVert[IPoint].end(), JPoint) 
							          != nnVert[IPoint].end());
					if (!LFound) {
						nnVert[IPoint].push_back(JPoint);
					}
				}
			}
		}
		
		// Output IFrame, number of total Points and number of cells in all files
		fprintf (oFileT   , "%7d  %7d  %7d  %7d\n"  , IFrame, NPoints   , 
		                                              NCells, NSurfCells);
		fprintf (oFileVert, "%7d  %7d  %7d  %7d\n"  , IFrame, NPoints   , 
		                                              NCells, NSurfCells);
		fprintf (oFileTetr, "%7d  %7d  %7d  %7d\n"  , IFrame, NPoints   , 
		                                              NCells, NSurfCells);                   

		// output domain size in all files
		for (int i = 0; i < 3; i++) {fprintf ( oFileT   , "%lf  " , H[i] );} 
		for (int i = 0; i < 3; i++) {fprintf ( oFileVert, "%lf  " , H[i] );} 
		for (int i = 0; i < 3; i++) {fprintf ( oFileTetr, "%lf  " , H[i] );} 
		fprintf ( oFileT   , "\n" );
		fprintf ( oFileVert, "\n" );
		fprintf ( oFileTetr, "\n" );
		
		//output Cell triangulation data
		for (ICell = 0; ICell < NCells; ICell++ ) {
			fprintf ( oFileT, "%7d     " , ICell+1 );
		
			for (int i = 0; i < 4; i++) {
				fprintf ( oFileT, " %7d " , Cell[ICell].vertex[i] );
			}
			fprintf ( oFileT, "    " );
			for (int i = 0; i < 4; i++) {
				fprintf ( oFileT, " %7d " , Cell[ICell].neighbor[i] );
			}
			
			if (Cell[ICell].OnBoundary) {fprintf ( oFileT, "%9d" , 1 );}
			else                        {fprintf ( oFileT, "%9d" , 0 );}
			
			fprintf ( oFileT, "\n" );
			if (Cell[ICell].OnBoundary) {
				fprintf ( oFileT, "           " );
				for (int i=0;i<4;i++) {
					fprintf ( oFileT, "%2d %2d %2d" , Cell[ICell].Img[i][0], 
					                                  Cell[ICell].Img[i][1], 
					                                  Cell[ICell].Img[i][2] );
					fprintf ( oFileT, "  " );
				} 
				fprintf ( oFileT, "\n" );
			}
		}
		
		// output vertex connectivity
		for (int i = 0 ; i < NPoints ; i++) {
			fprintf (oFileVert, "%7d  %7lu\n", i+1 , nnVert[i+1].size() );
			std::list<int>::const_iterator IVert;
			for (IVert = nnVert[i+1].begin(); IVert != nnVert[i+1].end(); ++IVert) {
				fprintf (oFileVert, "%7d  ", *IVert);
			}
			fprintf (oFileVert, "\n");
		}
	   
		// output reverse vertex-tetrahedra connectivity
		for (int i = 0 ; i < NPoints ; i++) {
			fprintf (oFileTetr, "%7d  %7lu\n", i+1 , nnTetr[i+1].size() );
			std::list<int>::const_iterator ITetr;
			for (ITetr = nnTetr[i+1].begin(); ITetr != nnTetr[i+1].end(); ++ITetr) {
				fprintf (oFileTetr, "%7d  ", *ITetr+1);
			}
			fprintf (oFileTetr, "\n");
		}

		// Clear Variables before next frame
		L.clear()      ; free(Cell)     ; T.clear()     ; 
		vIndex.clear() ; nnTetr.clear() ; nnVert.clear();
		NSurfCells = 0 ;
	}
	fclose(oFileT);
}

//*******************************************************************************
// ***Function readInput is used to access the input file and read frame data***
// INPUT:
// Type of input File: 1 - type 4, 2 - lammps dump
// Empty List to read Point coords         - Length of the box in each direction
// Empty Map to be used for Point indexing - Input Filename
//*******************************************************************************
void readInput(int iFileType, std::list<Point>  &L, double H[3]    , 
			   std::map<Point,int> &vIndex, const char*  filename) {
	double   Lx, Ly, Lz      ;
	int      IPoint          ;
	int      NAtoms          ;
	double   R[3]         	 ;
	double   box[3][2]       ;
	double   *xr, *yr, *zr   ;
	Point    p               ;
	std::string  line        ;


	int    NChains, IDummy;  
	float  Dummy; 
	
	if (FirstEntry) {
		ifs.open      (filename, std::ifstream:: in)  ;
		std::cout  << "Opened Input File" << std::endl; 
		FirstEntry  = false                           ;
		IFrame      = 0                               ;
	}

	// Type 4 file
	if (iFileType == 1) {
	
		ifs >> IFrame >> NChains >> IDummy >> NPoints >> IDummy  
			>> IDummy >> IDummy  >> IDummy >> Dummy; 
			
		for (int i=0; i<NChains; i++) { ifs >> IDummy; }
			
		ifs >> Lx >> Ly >> Lz;
		H[0] = Lx ; H[1] = Ly ; H[2] = Lz;

		// Read points data from file.  
		for (int IPoint=1; IPoint<NPoints+1; IPoint++) { 
			ifs >> R[0] >> R[1] >> R[2];
		
			// if p is on boundary, displace in the periodic domain
			for (int i=0 ; i<3 ; i++) {
				if ( abs(R[i] - 0.   ) < 1E-5 ) R[i] += 1E-3;
				if ( abs(R[i] - H[i] ) < 1E-5 ) R[i] -= 1E-3;
			}
			p = Point(R[0], R[1], R[2]);
			L.push_back(p)      ;
			vIndex[p]  =  IPoint;
		}	
	} else {
		IFrame += 1;

		for (int i=0; i<3; i++) { getline(ifs,line); } //skip lines
		
		getline(ifs,line);
		std::stringstream ss(line);
		ss >> NAtoms;

		xr = (double*) malloc(NAtoms*sizeof(double));
	    yr = (double*) malloc(NAtoms*sizeof(double));
        zr = (double*) malloc(NAtoms*sizeof(double));

		getline(ifs,line); // skip line

		// input box length. Also fix box to (0,0,0)
		for (int i=0; i<3; i++) {
			getline(ifs,line)                  ; 
			std::stringstream ss(line)         ;
			ss >> box[i][0] >> box[i][1]       ;
			H[i]    = box[i][1]  - box[i][0]   ;
		}

		getline(ifs,line); //skip line

		// read coordinates and displace them properly
		for (int i=0; i<NAtoms; i++) {

			getline(ifs,line)             ;
			std::stringstream ss(line)    ;
			ss >> IPoint       >> IDummy  ;
			ss << std::setprecision(6)    ; // stream with 6 decimal accuracy
            ss >> xr[IPoint-1] >> yr[IPoint-1]
			   >> zr[IPoint-1]; 

			// convert fractional coordinates to real coordinates
			xr[IPoint-1] = H[0]*xr[IPoint-1];
			yr[IPoint-1] = H[1]*yr[IPoint-1];
			zr[IPoint-1] = H[2]*zr[IPoint-1];
		}

		// put coordinates in the triangulation sites list
		for (int ip=0; ip<NAtoms; ip++) { 

			R[0] = xr[ip+1];
			R[1] = yr[ip+1];
			R[2] = zr[ip+1];

			// if p is on boundary, displace in the periodic domain
			for (int i=0 ; i<3 ; i++) {
				if ( abs(R[i] - 0.   ) < 1E-6 ) R[i] += 1E-6;
				if ( abs(R[i] - H[i] ) < 1E-6 ) R[i] -= 1E-6;
			}

			p = Point(R[0], R[1], R[2]);
			L.push_back(p)      ;
			vIndex[p]    =  ip+1;
		}
		free(xr) ; free(yr) ; free(zr) ;
	}
}


//*******************************************************************************
// ***Function findNeighbor returns the index of Neighbor Cell i of ICell***
// INPUT:
// 1) Index of ICell -
// 2) Neighbor i is the Cell that is located opposite of vertex i of ICell
// 3) One of each of the 3 common vertices between ICell and each neighbor
// 4) List of Tetrahedra that contain IPoint as a vertex
// NOTE: If the neighboring Cell is not located the function returns ZERO.
//*******************************************************************************
int findNeighbor ( int &ICell, int &i , int &IPoint , std::list<int> nnTetr ) {
	int  it, j  ;
	int  JCell  ;
	int  face[2];
	bool flag[2];
	
	j = 0;
	for (it = 0; it < 4; it++) {
		if ( (it == i) || (Cell[ICell].vertex[it] == IPoint) )  continue;
		face[j]   =   Cell[ICell].vertex[it];
		j        +=   1;
	}
	
	std::list<int>::const_iterator itetr;
	for ( itetr = nnTetr.begin(); itetr != nnTetr.end(); ++itetr) {
		if (*itetr == ICell) continue;
		
		for (int i = 0; i < 3; i++) {flag[i] = false;}
		for (j = 0; j < 2; j++) {
			for (it = 0; it < 4; it++){
				if ( face[j] == Cell[*itetr].vertex[it] ) {
					 flag[j] = true;
					 break;
				}
			}
		}
		
		if ( flag[0] && flag[1] ) {
			JCell = *itetr + 1;
			return  JCell    ;
		}
	}
	return 0;
}

//*******************************************************************************
// ***Function Get_Image_FromPos returns the Image of coordinates R, at dimension d.
// We assume that the simulation Cell extends between (0.,0.,0.) and 'CellDim(:)'
//*******************************************************************************
void Get_Image_FromPos(Point &p, double CellDim[3], int &i, int &d) {	
	double r;
	i  =  0 ;
	r = p[d];
	
	if (p[d] > CellDim[d])  i = (int)(p[d]/CellDim[d]);
	if (p[d] < 0.)          i = (int)(p[d]/CellDim[d]) - 1;
}
