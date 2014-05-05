// OptimizedPitRemoval.cpp : main project file.

//============== OPTIMIZED PIT REMOVAL V1.5.1 ================================
//Stephen Jackson, Center for Research in Water Resources, University of Texas at Austin; 
//srj9@utexas.edu
//Acknowledgements: Soille, P. (2004), Optimal removal of spurious pits in grid digital elevation models, Water Resour. Res., 40, W12509; 
//Dr. David Maidment, University of Texas at Austin; 
//Dr. David Tarboton, Utah State University; 
//TauDEM (Terrain Analysis Using Digital Elevation Models), Hydrology Research Group, Utah State University

#include "stdafx.h"
#include <vector> //Used for vector arrays
#include <fstream> //Used for reading ASCII
#include <iostream> 
#include <sstream> //Used for converting numbers and strings
#include <stdlib.h> //Used for converting numbers and strings
#include <string>
#include <queue>
using namespace std;
using namespace System;
using std::vector;
#include <map> //Used to create cost functions

//Using vectors rather than traditional arrays as they seem to be more intuitive without a large performance loss
//When declaring vectors, use syntax: std::vector<int> v; // declares a vector of integers
//This saves from having to bring in the whole std namespace (using namespace std;) so keeps the final product from getting bloated.
//http://www.codeguru.com/cpp/cpp/cpp_mfc/stl/article.php/c4027/C-Tutorial-A-Beginners-Guide-to-stdvector-Part-1.htm

// ====== FILES FROM TAUDEM ================================
#include "commonLib.h" //I'm using a simplified version of commonLib.h that does not reference as much. On a computer with working TauDEM, the standard commonLib.h should be used instead.

// ====== CLASSES AND STRUCTURES ================================
struct point
{
	//A point contains a raster cell location in 1-D array format (0-based, starting from upper left corner) and the associated elevation
	int id;	//0-based grid cell number in a 1-D array
	double elev; //Elevation associated with the point
};
class ComparePoint
{
	public:
    bool operator()(point& p1, point& p2) // Returns true if p1 is greater than p2. This makes the priority queue return low elevations first
    {
       if (p1.elev > p2.elev) return true;
       return false;
    }
};
class CompareElevation
{
	public:
    bool operator()(double e1, double e2) // Returns true if p1 is greater than p2. This makes the priority queue return low elevations first
    {
       if (e1 > e2) return true;
       return false;
    }
};

// ====== FUNCTION PROTOTYPES ================================

void PitRemove(char *demfile, char *newfile, bool use_hybrid, double cut_coef, bool balance, double step_size);
void InitializeMainQueue();
void IterateMainQueue();
void PitRemoveHybrid(int PitID);
void CutToElevation(int PitID);

bool ParseInputs(int argc, char **argv);
void ImportTerrain(char *demfile);

void GetAsciiHeader(char *demfile);
void GetAsciiDoubleVector(char *inputpath,vector<double>& Vect);
void SaveAscii(char *outputpath, vector<double> Vect);

bool IsBorder(int ID);
bool NeighborNoValue(int ID);
void GetNeighbors(int ID);
void AddToMainQueue(int ID, bool ConfirmDescend);
void GetDryNeighbors(int ID);
void SetFlowDirection(int FromID, int ToID);

double GetIdealFillLevel(double CrestElev);

int TraceFlow(int FromID, int FlowDirection);
void GetFillSteps(int PitID, double CrestElev, priority_queue<point, vector<point>, ComparePoint>& FillStepQueue);
void FillToElevation(int PitID, double FillElev);
void GetDepressionExtent(int PitID, double CrestElev);
double GetCrestElevation(int PitID);
void CreateFillFunction(int PitID, double CrestElev);
void CreateCutFunction(int PitID, double CrestElev);
bool CheckCell(int ID, int Direction, int& CurNeighborID, double CrestElev);
bool IsLocalMinimum(int CurID);

// ====== GLOBAL VARIABLES  ================================
//Global Variables tend to be strongly frowned upon. There may be a better option.
//These are very large arrays (10^6 or larger). Keeping them private and passing as through the function call was very time consuming at runtime.
vector<double> Terrain(0); //This begins as the input DEM and is modified by the algorithm.
vector<int> Direction(0); //8-direction indicator of which cell caused the current cell to become flooded
vector<int> Flooded(0); //Value of 0=unflooded; Value of 1=flooded; Value of 2=flooded and has confirmed descending path to an outlet
vector<bool> Checked(0); //Used to determine the extent of a depression. Reset after each depression is identified
vector<bool> BlankBool(0);  //Used for clearing the contents of a vector-bool of Terrain size
vector<int> Depression(0); //Stores the extent of 
vector<int> BlankInt(0); //Used for clearing the contents of a vector-int  of Zero size
vector<int> Neighbors(8); //Stores the ID of the current cell's eight neighbors
vector<double> IsPit(0); //Optional binary record of which cells were identified as local minima.

map<double,double> CutFunction; //Stores paired elevation/cost data
map<double,double> FillFunction; //Stores paired elevation/cost data
map<double,double> BlankMap; //Used to reset Cut and Fill Functions
priority_queue<point, vector<point>, ComparePoint> MainQueue;
priority_queue<point, vector<point>, ComparePoint> NeighborQueue;
priority_queue<point, vector<point>, ComparePoint> DepressionExtent;
priority_queue<point, vector<point>, ComparePoint> EmptyQueue;
priority_queue<point, vector<point>, ComparePoint> FillStepQueue;

//Common variables which would otherwise be passed back and forth frequently
double PitElev;
int numCols, numRows; // numCols/numRows are 1-based grid size from ASCII definition
int Size;
double Xll, Yll, cell_size, no_data;
enum ModeType {CUT, BAL, MIN_COST, MIN_CELL};
int Mode;
double step_size;
char demfile[MAXLN], newfile[MAXLN], pitfile[MAXLN]; //MAXLN defined in commonLib.h as 4096
bool SavePits;

// ===== MAINLINE CONTROL ================================

int main(int argc,char **argv)
{
    //Inputs: -z: Input DEM file path (as ASCII). -fel: Output filepath (as ASCII). -pit: Optional Output pit location filepath (as ASCII)
	//Inputs: -mode Optional mode value ('cut' for cut only, 'bal' for minimize sum(delta z), 'mincost' for minimize sum(abs(delta z)), 'mincell' for minimize number of cells modified)
	//Inputs: -step: Optional step_size value (double, default 0.1)

	//MINCELL has not yet been implemented, so the ArcMap tool does not have that option

	//Initialize variables
	Mode = MIN_COST;
	step_size = 0.1; //Vertical step size for optimization calculations
	SavePits = false;
		
	//Parse Inputs
	if (ParseInputs(argc, argv));
	else goto errexit;

	//Import Terrain
	printf("\nImporting Input DEM: ");
	printf(demfile);
	ImportTerrain(demfile);

	//Initialize Main Queue
	printf("\nInitializing Main Queue...");
	InitializeMainQueue();

	//Iterate Main Queue, removing all pits
	printf("\nIterating Main Queue...");
	IterateMainQueue();
	
	//Save the output files
	printf("\nSaving output dem: ");
	printf(newfile);
	SaveAscii(newfile, Terrain);
	if (SavePits == true)
	{
		printf("\nSaving pit locations: ");
		printf(pitfile);
		SaveAscii(pitfile, IsPit);

	}
	printf("\nProcess Completed");

    return 0;

	errexit:
	printf("ERROR!");
}

bool ParseInputs(int argc, char **argv)
{
	bool Parse = true;
	char inputstring[MAXLN];
	int i = 1;
	if(argc < 4)
    {  
       printf("Error: Inputs required");
	   Parse = false;
	   i = argc;
    }
	while(argc > i)
	{
		if(strcmp(argv[i], "-z") == 0)
		{
			i++;
			if(argc > i)
			{
				strcpy_s (demfile, MAXLN, argv[i]);
				i++;
			}
			else Parse = false;
		}
		else if(strcmp(argv[i],"-fel")==0)
		{
			i++;
			if(argc > i)
			{
				strcpy_s (newfile, MAXLN, argv[i]);
				i++;
			}
			else Parse = false;
		}
		else if(strcmp(argv[i],"-pit")==0)
		{
			i++;
			if(argc > i)
			{
				strcpy_s (pitfile, MAXLN, argv[i]);
				SavePits = true;
				i++;
			}
			else Parse = false;
		}
		else if(strcmp(argv[i],"-mode")==0)
		{
			i++;
			if(argc > i)
			{
				if (strcmp(argv[i],"cut")==0)
				{
					Mode = CUT;
				}
				else if (strcmp(argv[i],"bal")==0)
				{
					Mode = BAL;
				}
				else if (strcmp(argv[i],"mincost")==0)
				{
					Mode = MIN_COST;
				}
				else if (strcmp(argv[i],"mincell")==0)
				{
					Mode = MIN_CELL;
				}
				i++;
			}
			else Parse = false;
		}
		else if(strcmp(argv[i],"-step")==0)
		{
			i++;
			if(argc > i)
			{
				strcpy_s (inputstring, argv[i]);
				step_size = atof(inputstring);
				i++;
			}
			else Parse = false;
		}
		else 
		{
			Parse = false;
		}
	}

	return Parse;
}
void ImportTerrain(char *demfile)
{
	//Get Header Information
	GetAsciiHeader(demfile);

	//Resize vectors
	Terrain.resize(numCols*numRows); //Values input from file, modified throughout program
	Direction.resize(numCols*numRows);  //Starts empty
	Flooded.resize(numCols*numRows); //Starts empty
	Checked.resize(numCols*numRows);
	BlankBool.resize(numCols*numRows);
	if (SavePits == true)
		IsPit.resize(numCols*numRows);
	
	Size = Terrain.size();
		
	//Import Terrain values
	GetAsciiDoubleVector(demfile, Terrain);
}
void InitializeMainQueue()
{
	//The entire DEM is scanned and all outlets are added to the Main Queue
	//An outlet is defined as a cell that is either on the border of the grid or has a neighbor with no_data
	//This allows internal points of no data to be used as outlets
	size_t i=0;
	bool AddToQueue = 0;

	for (i=0; i < Size; i++)
	{
		//Test if cell is on border
		if(IsBorder(i)) 
			AddToQueue = 1; 
		//Test if cell has a neighbor with no data
		else if (NeighborNoValue(i)) 
			AddToQueue = 1;

		//Add cells to queue
		if (AddToQueue == 1)
		{
			AddToMainQueue(i, true); //All outlet cells by definition have a path to the outlet, thus ConfirmDescend = true.
			Direction.at(i)=0;
		}
		AddToQueue = 0;
	}
}
void IterateMainQueue()
{
	//MainQueue will gradually have all cells in grid added to it, starting from the borders and those next to no_data_value cells and 
	//progressing upward with a rising flood. 
	//Each cell is removed from the queue in order of lowest elevation and checked to see if it is a pit minimum.
	//If it is a minimum, the pit removal procedures are run.
	//After each cell is removed from the queue, each of its unflooded neighbors is added.
	//The vector Direction keeps track of which neighboring cell caused that cell to become flooded. This can be used to trace a natural path to an outlet.
	//This continues until all cells have been flooded and the Main Queue is empty, and thus all pits removed.
	//Note: Preprocessing to identify all pits does not seem a good strategy, as the number of pits which need to be handled is dependant on the step size. 
	//(A coarser step size will cause more pits to be filled before they are removed from the queue.)

	point CurCell;
	point CurNeighbor;
	int CurPitNum=0;
	bool ConfirmDescend=false;
			
	while(!MainQueue.empty())
	{
		CurCell = MainQueue.top();
		MainQueue.pop();

		//Check Cell to determine if it is a Pit Minimum
		if(IsLocalMinimum(CurCell.id))
		{
			//printf("\nRemoving Pit: %i", CurPitNum);
			//User-selected method for removing pit
			if(Mode == CUT)
			{
				CutToElevation(CurCell.id);
			}
			else
			{
				PitRemoveHybrid(CurCell.id);
			}
			//printf(".");
			CurPitNum++;
		}
		else //Some cells within a depression may still be classified as Flooded=1 after pit has been removed. Need to correct this
		{
			if(Flooded.at(CurCell.id)==1)
			{
				GetNeighbors(CurCell.id);
				int i;
				for (i=0;i<8;i++)
				{
					if (Neighbors.at(i)>-1)
					{
						if ((Flooded.at(Neighbors.at(i))==2)&&(Terrain.at(Neighbors.at(i))<=Terrain.at(CurCell.id)))
						{
							Flooded.at(CurCell.id)=2;
							break;
						}
					}
				}
			}
		}

		//Add unflooded neighbors to Main Queue and identify the direction the flooding came from
		GetDryNeighbors(CurCell.id);
		while(!NeighborQueue.empty())
		{
			CurNeighbor = NeighborQueue.top();
			NeighborQueue.pop();

			ConfirmDescend=false;
			if(Terrain.at(CurNeighbor.id)>=Terrain.at(CurCell.id)&&(Flooded.at(CurCell.id)==2)) 
				ConfirmDescend=true;
			AddToMainQueue(CurNeighbor.id, ConfirmDescend);
			SetFlowDirection(CurNeighbor.id, CurCell.id);
		}
	}
}

void PitRemoveHybrid(int PitID)
{
	//This handles the cases where Mode = BAL or MIN_COST
	//Cost is defined as the difference between the original terrain and the modified terrain for each point, summed across all modified points.
	//BAL minimizes the sum of cost across all cells modified for each depression (i.e., tries to get Cut Volume = Fill Volume)
	//MIN_COST minimizes the sum of the absolute value of cost across all cells modified for each depression (i.e. tries to have the least disturbance to the original terrain)
	//BAL will *always* have a mix of cut and fill (for suitably small step size)
	//MIN_COST will *sometimes* have a mix of cut and fill, but some pits may be best removed using only cut or only fill
	//Filling modifies a 2-D region (depression). Cutting modifies a 1-D path to an outlet. Thus MIN_COST will often result in more cut than fill.

	double CrestElev, IdealFillLevel;
	PitElev = Terrain.at(PitID);
	
	//Get Crest Elevation
	CrestElev = GetCrestElevation(PitID);
	
	//Get Depression Extents and set PitElev
	GetDepressionExtent(PitID, CrestElev);

	//Create Cut Function
	CreateCutFunction(PitID, CrestElev);

	//Create Fill Function
	CreateFillFunction(PitID, CrestElev);

	if(PitElev < CrestElev)
	{
		//Find desired Fill Elevation
		IdealFillLevel = GetIdealFillLevel(CrestElev);

		//Modify Terrain to remove Pit
		FillToElevation(PitID, IdealFillLevel);
		CutToElevation(PitID);
	}

	Flooded.at(PitID) = 2; //Confirm that cell has descending path to outlet (i.e. Pit has been removed)
	
	//Reset global variables
	CutFunction = BlankMap;
	FillFunction = BlankMap;
	Depression = BlankInt;
}


double GetCrestElevation(int PitID)
{
	//Find the highest point along the path to an outlet

	bool ReachedOutlet = 0;
	double Crest;
	int CurID, NextID;

	//Initialize Crest as the Pit Minimum elevation
	CurID = PitID;
	Crest = PitElev; 

	
	while(!ReachedOutlet)
	{
		NextID = TraceFlow(CurID, Direction.at(CurID));
		if(NextID<0) //CurID is a border cell
			ReachedOutlet = 1;
		else if (Terrain.at(NextID) == no_data) //CurID is next to an internal outlet
			ReachedOutlet = 1;
		else if((Terrain.at(NextID) < Terrain.at(PitID))&&(Flooded.at(NextID)==2)) //NextID is lower than Pit and NextID on confirmed descending path to outlet
			ReachedOutlet = 1;
		else 
		{
			if((Terrain.at(NextID) > Crest) && (Terrain.at(NextID) != no_data)) Crest = Terrain.at(NextID);
		}
		CurID = NextID;
	}
	return Crest;
}
void GetDepressionExtent(int PitID, double CrestElev)
{
	//Makes an elevation ordered list of every cell in the depression
	//A compound depression (neighboring pit with separating ridge lower than crest elevation) is treated as a separate depression. That pit will be removed later.
	int CurNeighborID;
	int CurID;
	point CurPoint;
	point NeighborPoint;
	priority_queue<point, vector<point>, ComparePoint> DepressionQueue;
		
	CurPoint.id = PitID;
	CurPoint.elev = Terrain.at(PitID);
	DepressionQueue.push(CurPoint);
	Depression.push_back(CurPoint.id);

	while(!DepressionQueue.empty())
	{
		CurPoint = DepressionQueue.top();
		DepressionQueue.pop();
		CurID = CurPoint.id;

		//Get each neighbor cell that is lower than the Crest elevation and with elevation greater than or equal to the present cell
		int i;
		for (i=0;i<8;i++)
		{
			if (CheckCell(CurID, i, CurNeighborID, CrestElev))
			{
				NeighborPoint.id=CurNeighborID;
				NeighborPoint.elev=Terrain.at(CurNeighborID);

				DepressionQueue.push(NeighborPoint);
				Depression.push_back(CurNeighborID);
			}
		}
	}
	Checked = BlankBool;
}
void CreateCutFunction(int PitID, double CrestElev)
{
	//Create a list of elevations from the Pit to the Crest, incremented by step size
	//For each elevation in list, sum all positive difference between Terrain elev and List Elev for all cells along path from pit to outlet
	//CutFunction represents the (Cut Volume / Cell Area)
	bool ReachedOutlet = 0;
	int CurID, NextID;
	double CellElev = Terrain.at(PitID);
	double CurStep;
	double OldCut;
	CurID = PitID;

	//Initialize all values
	CurStep = PitElev;
	while (CurStep < CrestElev)
	{
		CutFunction[CurStep] = 0;
		CurStep = CurStep + step_size;
	}
	CutFunction[CrestElev] = 0; //Make sure an entry for the Crest is added, but that cut will always be zero

	//Calculate Cost
	while(!ReachedOutlet)
	{
		NextID = TraceFlow(CurID, Direction.at(CurID));
		if(NextID<0) //CurID is a border cell
			ReachedOutlet = 1;
		else if (Terrain.at(NextID) == no_data) //CurID is next to an internal outlet
			ReachedOutlet = 1;
		else if((Terrain.at(NextID) < PitElev)&&(Flooded.at(NextID)==2)) //NextID is lower than Pit and NextID on confirmed descending path to outlet
			ReachedOutlet = 1;
		else 
		{
			CellElev = Terrain.at(NextID);
			CurStep = PitElev;
			while (CurStep < CellElev)
			{
				if (CutFunction.count(CurStep)) //check if entry exists
					OldCut = CutFunction.find(CurStep)->second;
				else
					OldCut = 0;
				CutFunction[CurStep] = OldCut + CellElev - CurStep;
				CurStep = CurStep + step_size;
			}
		}
		CurID = NextID;
	}	
}
void CreateFillFunction(int PitID, double CrestElev)
{
	//Create a list of elevations from the Pit to the Crest, incremented by step size
	//For each elevation in list, sum all positive difference between List elev and Terrain Elev for all cells within the depression
	//FillFunction represents the (Fill Volume / Cell Area)
	int CurID;
	double CurGroundElev;
	double OldFill;
	double CurStep;

	//Initialize
	CurStep = PitElev;
	while (CurStep < CrestElev)
	{
		FillFunction[CurStep] = 0;
		CurStep = CurStep + step_size;
	}
	FillFunction[CrestElev] = 0; //Make sure an entry for the Crest is added

	//Calculate Cost
	int i;
	for (i=0;i<Depression.size();i++)
	{
		CurID = Depression.at(i);
		CurGroundElev = Terrain.at(CurID);

		CurStep = PitElev;
		while (CurStep < CrestElev)
		{
			if (CurStep > CurGroundElev)
			{
				OldFill = FillFunction.find(CurStep)->second;
				FillFunction[CurStep] = OldFill + CurStep - CurGroundElev;
			}
			
			CurStep = CurStep + step_size;

			//Exit Conditions - make sure the last step gets accounted for
			if (CurStep > CrestElev)
			{
				OldFill = FillFunction.find(CrestElev)->second;
				FillFunction[CrestElev] = OldFill + CrestElev - CurGroundElev;
			}
		}
	}
}
double GetIdealFillLevel(double CrestElev)
{
	//If balancing, find point of minimum difference between cut and fill
	//If not balancing, find point of minimum total cost
	//NOTE: CutFunction and FillFunction are currently both defined with only positive values (rather than cut being negative, so the calculations below account for this.

	double FillLevel;
	double CurTotalCost, CurCutCost, CurFillCost, MinCost, CurDifference, MinDifference;
	double CurStep;

	switch (Mode)
	{	
	case BAL: 
		{
			//Initialize
			CurStep = PitElev;
			CurCutCost = CutFunction.find(CurStep)->second;
			CurFillCost = FillFunction.find(CurStep)->second;
			MinDifference = fabs(CurFillCost - CurCutCost);
			FillLevel = CurStep;

			//Optimize
			while (CurStep < CrestElev)
			{
				CurCutCost = CutFunction.find(CurStep)->second;
				CurFillCost = FillFunction.find(CurStep)->second;
				CurDifference = fabs(CurFillCost - CurCutCost);

				if (CurDifference < MinDifference)
				{
					MinDifference = CurDifference;
					FillLevel = CurStep;
				}
				CurStep = CurStep + step_size;
			}
		}
		break;
	case MIN_COST:
		{
			//Initialize
			CurStep = PitElev;
			MinCost = CutFunction.find(CurStep)->second;
			MinCost = MinCost;
			FillLevel = CurStep;

			if (MinCost > 0)
			{
				//Optimize
				while (CurStep < CrestElev)
				{
					CurCutCost = CutFunction.find(CurStep)->second;
					CurFillCost = FillFunction.find(CurStep)->second;
					CurTotalCost = CurCutCost + CurFillCost;
					if (CurTotalCost < MinCost)
					{
						MinCost = CurTotalCost;
						FillLevel = CurStep;
					}
					CurStep = CurStep + step_size;
				}
			}
		}
		break;
	}
	return FillLevel;
}

bool IsBorder(int ID)
{
	//Tests if cell is on the outer edge of the grid, and thus is an outlet
	bool border = 0;

	if (ID < numCols) 
		border=1;
	else if (Size-ID <numCols + 1)
		border=1;
	else if (!(ID % numCols))
		border=1;
	else if (!((ID+1) % numCols))
		border = 1;
	else
		border=0;

	return border;
}
bool NeighborNoValue(int ID)
{
	//Tests if cell is next to a cell with no_data, and thus is an outlet
	bool novalue = 0;

	GetNeighbors(ID);

	int i;
	for (i=0;i<8;i++)
	{
		if (Neighbors.at(i)==-1) 
		{
			novalue=1;
			break;
		}
		else if (Terrain.at(Neighbors.at(i))==no_data)
		{
			novalue=1;
			SetFlowDirection(ID, Neighbors.at(i));
			break;
		}
	}
	return novalue;
}
void AddToMainQueue(int ID, bool ConfirmDescend)
{
	//Adds cell to main queue. Input ConfirmDescend = true if there is a known continuously descending path to an outlet.
	//Cells are added to Main Queue at the moment they become flooded. 
	//The vector Flooded is a permanent record of whether each cell has been flooded and whether there is a known path to an outlet
	
	point CurCell;
	if (Flooded.at(ID)==0)
	{
		CurCell.id = ID;
		CurCell.elev = Terrain.at(ID);
		MainQueue.push(CurCell);

		if(ConfirmDescend)
		{
			Flooded.at(ID)=2;
		}
		else
		{
			Flooded.at(ID) = 1;
		}
	}
}

bool IsLocalMinimum(int CurID)
{
	//A local Minimum has been found if the cell does not have a confirmed path to an outlet, and all neighbors are either higher or the same elevation but flooded.
	//This means for a wide pit with a flat bottom, the last flooded cell (roughly farthest from the outlet) is considered the minimum.

	bool IsMinimum = true;
	if (Flooded.at(CurID)==2)  //Cell is on confirmed path to outlet
	{
		IsMinimum = false;
	}
	else //Flooded = 1
	{
		GetNeighbors(CurID);
		int i = 0;
		while ((i<8)&&(IsMinimum==true))
		{
			if (Neighbors.at(i)>-1)
			{
				if( Terrain.at(Neighbors.at(i)) < Terrain.at(CurID))
					IsMinimum = false;
				if((Terrain.at(Neighbors.at(i)) == Terrain.at(CurID)) && (Flooded.at(Neighbors.at(i)) == 0))
					IsMinimum = false;
			}
			i++;
		}		
	}
	if (SavePits == true)
		if (IsMinimum==true) IsPit.at(CurID) = 1;
	return IsMinimum;
}
void GetNeighbors(int ID)
{
	//Neighbors is a 0-7 vector, defined clockwise from Northwest
	//Returns the ID value for the eight neighbors, with -1 for cells off the grid

	//Northwest
	if(ID<numCols)
		Neighbors.at(0)=-1;
	else if (!(ID % numCols))
		Neighbors.at(0)=-1;
	else 
		Neighbors.at(0)=ID-1-numCols;

	//North
	if(ID<numCols)
		Neighbors.at(1)=-1;
	else 
		Neighbors.at(1)=ID-numCols;

	//Northeast
	if(ID<numCols)
		Neighbors.at(2)=-1;
	else if (!((ID+1) % numCols))
		Neighbors.at(2)=-1;
	else 
		Neighbors.at(2)=ID+1-numCols;

	//East
	if (!((ID+1) % numCols))
		Neighbors.at(3)=-1;
	else 
		Neighbors.at(3)=ID+1;

	//Southeast
	if (Size-ID <numCols + 1)
		Neighbors.at(4)=-1;
	else if (!((ID+1) % numCols))
		Neighbors.at(4)=-1;
	else 
		Neighbors.at(4)=ID+1+numCols;

	//South
	if (Size-ID <numCols + 1)
		Neighbors.at(5)=-1;
	else
		Neighbors.at(5)=ID+numCols;

	//Southwest
	if (Size-ID <numCols + 1)
		Neighbors.at(6)=-1;
	else if (!(ID % numCols))
		Neighbors.at(6)=-1;
	else 
		Neighbors.at(6)=ID-1+numCols;

	//West
	if (!(ID % numCols))
		Neighbors.at(7)=-1;
	else 
		Neighbors.at(7)=ID-1;
}
void GetDryNeighbors(int ID)
{
	//Adds all neighbors which have not yet been flooded to NeighborQueue.

	point DryNeighbor;
	GetNeighbors(ID);

	int i;
	for (i=0;i<8;i++)
	{
		if (Neighbors.at(i)>-1)
		{
			if (!Flooded.at(Neighbors.at(i)))
			{
				DryNeighbor.id = Neighbors.at(i);
				DryNeighbor.elev = Terrain.at(Neighbors.at(i));
				NeighborQueue.push(DryNeighbor);
			}
		}
	}
}

void SetFlowDirection(int FromID, int ToID)
{
	 //Flow direction is FROM current cell TO cell which caused it to flood, clockwise from East (1,2,4,8,16,32,64,128)
     //If two cells are not neighbors or if a neighbor is off the grid/ has no_data, direction set to 0.
	
	if(ToID == FromID + 1)
		//Flow is to East
		Direction.at(FromID)=1;
	else if(ToID == FromID + 1 + numCols)
		//Flow is to Southeast
		Direction.at(FromID)=2;
	else if(ToID == FromID + numCols)
		//Flow is to South
		Direction.at(FromID) = 4;
	else if(ToID == FromID - 1 + numCols)
		//Flow is to Southwest
		Direction.at(FromID) = 8;
	else if(ToID == FromID - 1)
		//Flow is to West
		Direction.at(FromID) = 16;
	else if(ToID == FromID - 1 - numCols)
		//Flow is to Northwest
		Direction.at(FromID) = 32;
	else if(ToID == FromID - numCols)
		//Flow is to North
		Direction.at(FromID) = 64;
	else if(ToID == FromID + 1 - numCols)
		//Flow is to Northeast
		Direction.at(FromID) = 128;
	else
		//Cells are not neighbors
		Direction.at(FromID) = 0;
}
int TraceFlow(int FromID, int FlowDirection)
{
	//Returns the cell pointed to by the direction grid at the given location.
	//If flow direction equals 0, This is a border cell. Return -1.
	int ToID;

	if(FlowDirection == 1)
		ToID = FromID + 1;
	else if(FlowDirection == 2)
		ToID = FromID + 1 + numCols;
	else if(FlowDirection == 4)
		ToID = FromID + numCols;
	else if(FlowDirection == 8)
		ToID = FromID - 1 + numCols;
	else if(FlowDirection == 16)
		ToID = FromID - 1;
	else if(FlowDirection == 32)
		ToID = FromID - 1 - numCols;
	else if(FlowDirection == 64)
		ToID = FromID - numCols;
	else if(FlowDirection == 128)
		ToID = FromID + 1 - numCols;
	else 
		ToID = -1;

	return ToID;
}

void FillToElevation(int PitID, double FillElev)
{
	//Fills all cells within a depression to the specified elevation
	int CurID;
		
	if ((Terrain.at(PitID) < FillElev) && (Terrain.at(PitID) != no_data)) Terrain.at(PitID) = FillElev;

	int i;
	for (i=0;i<Depression.size();i++)
	{
		CurID = Depression.at(i);
		if ((Terrain.at(CurID) < FillElev) && (Terrain.at(CurID) != no_data)) 
		{
			Terrain.at(CurID) = FillElev;
		}
	}
}
void CutToElevation(int PitID)
{
	//This is used both when Mode = CUT_ONLY and by the Hybrid process after Filling has occured
	//Starting at the Pit, Backtrack along flow direction and set each cell equal to pit elevation
	//Stop when either an equal or lower elevation is found, or when an outlet is reached
	PitElev = Terrain.at(PitID);
	bool ReachedOutlet = 0;
	int CurID, NextID;

	CurID = PitID;

	while(!ReachedOutlet)
	{
		NextID = TraceFlow(CurID, Direction.at(CurID));
		if(NextID<0) //CurID is a border cell
			ReachedOutlet = 1;
		else if (Terrain.at(NextID) == no_data) //CurID is next to an internal outlet
			ReachedOutlet = 1;
		else if((Terrain.at(NextID) < PitElev)&&(Flooded.at(NextID)==2)) //NextID is lower than Pit and NextID on confirmed descending path to outlet
			ReachedOutlet = 1;
		else 
		{
			if ((Terrain.at(NextID) > PitElev) && (Terrain.at(NextID) != no_data))
				Terrain.at(NextID) = PitElev;
			Flooded.at(NextID) = 2; //Confirm that cell has descending path to outlet
		}
		CurID = NextID;
	}
}


bool CheckCell(int ID, int Direction, int& CurNeighborID, double CrestElev)
{
	bool Check = false;
	int Size = Terrain.size();
	switch (Direction)
	{	
	case 0: //Northwest
		{
			if(!(ID<numCols) && (ID % numCols))
			{
				CurNeighborID = ID-1-numCols;
				Check = true;
			}
		}
		break;
	case 1: //North
		{
			if(!(ID<numCols))
			{
				CurNeighborID = ID-numCols;
				Check = true;
			}
		}
		break;
	case 2: //Northeast
		{
			if(!(ID<numCols) && ((ID+1) % numCols))
			{
				CurNeighborID = ID+1-numCols;
				Check = true;
			}
		}
		break;
	case 3: //East
			
		{
			if (((ID+1) % numCols))
			{
				CurNeighborID = ID+1;
				Check = true;
			}
		}
		break;
	case 4: //Southeast
		{
			if (!(Size-ID <numCols + 1) && ((ID+1) % numCols))
			{
				CurNeighborID = ID+1+numCols;
				Check = true;
			}
		}
		break;
	case 5: //South
		{
			if (!(Size-ID <numCols + 1))
			{
				CurNeighborID = ID+numCols;
				Check = true;
			}
		}
		break;
	case 6: //Southwest
		{
			if (!(Size-ID <numCols + 1) && (ID % numCols))
			{
				CurNeighborID = ID-1+numCols;
				Check = true;
			}
		}
		break;
	case 7: //West
		{
			if ((ID % numCols))
			{
				CurNeighborID = ID-1;
				Check = true;
			}
		}
		break;
	}
	
	if (Check==true)
	{
		if (!( (Checked.at(CurNeighborID)==0) && (Terrain.at(CurNeighborID) < CrestElev) && (Terrain.at(CurNeighborID)>=Terrain.at(ID))))
		{
			Check = false;
		}
		Checked.at(CurNeighborID) = 1;
	}

	return Check;
}


// ===== FILE I/O ================================
void GetAsciiHeader(char *inputpath)
{
	ifstream inputfile;
	inputfile.open(inputpath);
	if(inputfile)  //if statement used to handle errors
	{
		std::string inputstring;
		inputfile >> inputstring >> inputstring;
		stringstream(inputstring) >> numCols;
		inputfile >> inputstring >> inputstring;
		stringstream(inputstring) >> numRows;

		inputfile >> inputstring >> inputstring;
		Xll = atof(inputstring.c_str());
		inputfile >> inputstring >> inputstring;
		Yll = atof(inputstring.c_str());
		inputfile >> inputstring >> inputstring;
		cell_size = atof(inputstring.c_str());
		inputfile >> inputstring >> inputstring;
		no_data = atof(inputstring.c_str());
	}
	inputfile.close();
}
void GetAsciiDoubleVector(char *inputpath, vector<double>& Vect)
{
	ifstream inputfile;
	inputfile.open(inputpath);
	if(inputfile)  //if statement used to handle errors
	{
		std::string inputstring;
		size_t i=0;

		//Skip the header
		size_t numheader = 12;
		while (i < numheader)
		{
			inputfile >> inputstring;
			i++;
		}

		//Input the data
		i=0;
		for (i=0; i < Vect.size(); i++)
		{
			inputfile >> inputstring;
			Vect.at(i) = atof(inputstring.c_str());
		}
	}
	inputfile.close();
}
void SaveAscii(char *outputpath, vector<double> Vect)
{
	ofstream outputfile;
	outputfile.open(outputpath);
	//Set ouput precision. This is the maximum number of digits that will be printed (before & after decimal). Currently set to a number that should be conservative; there may be a better explicit approach.
	int prec = 10;
	outputfile.unsetf(ios::floatfield);  
	outputfile.precision(prec);

	//Write the header
	outputfile << "ncols     " << numCols << "\n";
	outputfile << "nrows     " << numRows << "\n";
	outputfile << "xllcorner     " << Xll << "\n";
	outputfile << "yllcorner     " << Yll << "\n";
	outputfile << "cellsize      " << cell_size << "\n";
	outputfile << "NODATA_value      " << no_data << "\n";

	//Write data
	size_t i=0;
	int curcol=1;
	while(i<Size)
	{
		for(curcol=1;curcol<numCols+1;curcol++)
		{

			outputfile << Vect.at(i) << " ";
			i++;
		}
		outputfile << "\n";
	}
	outputfile.flush();
	outputfile.close();
}

