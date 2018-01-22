/******************************************
 * headfile containing most global variables
 ******************************************/
#include <stdio.h>

#define MAX_NumberOfParticles 50000
//#define MAX_NumberOfNeighbors 500 // max number of neighboring particles in the Verlet list
#define MAX_drBins 10000 // for g(r)

#define MAX_NumberOfCells 10000
#define MAX_NumberOfNeighborCells 30
//#define MAX_Dimension 3

#define SQR(x) ((x)*(x))
#define CUBIC(x) ((x)*(x)*(x))
#define MAX(x,y) ((x)>(y) ? (x) : (y))
#define MIN(x,y) ((x)<(y) ? (x) : (y))

/************************/
extern int CellSwitch; // use or not use cell list
struct LinkedList{
 int WhichCell;
 int Prev;
 int Next;
};
extern struct LinkedList CellTrack[MAX_NumberOfParticles];
extern int HeadOfChain[MAX_NumberOfCells];
extern double rcell; // size of cell rcell > sigma
extern int LCell,NumberOfCells;
extern int NumberOfNeighborCells;
extern double rc;// rcutoff
extern int NeighborCellList[MAX_NumberOfCells][MAX_NumberOfNeighborCells];
extern int CellCoordinate[MAX_NumberOfCells][8];
/************************/
 
extern double randomseed;

extern int CompressSwitch;
extern double collisionfrequency; 
extern double maxfrequency;
extern double timewindow;

extern double compresstime;// compress every such time interval 
extern double compressfactor;// 0~1
extern int icompress,jcompress;// particle i,j which gives dmin
extern double compressratio; // particle size increase ratio
extern double tolerance,timetol;


extern double TimeBig; // collision time upper bound
extern double CollisionTime[MAX_NumberOfParticles];
extern int CollisionPartner[MAX_NumberOfParticles];
extern double EscapeTime[MAX_NumberOfParticles]; // time needed to escape a cell
extern int EscapeID;
//extern int Escapexyz[MAX_NumberOfParticles]; // 1 x 2 y 3 z
//extern int OnBoundx[MAX_NumberOfParticles],OnBoundy[MAX_NumberOfParticles],OnBoundz[MAX_NumberOfParticles]; // xyz: 0 off; 1 on boundary
//extern int Cellold[MAX_NumberOfParticles],Cellnew[MAX_NumberOfParticles];

extern int JobIndex;
extern double  rate; // cooling rate

extern int rInitialType,vInitialType;
extern int PackType; // SC, BCC or FCC
extern int PotentialType; //1: shifted-force L-J
extern int EnsembleType; // 0: NVE 1: isokinetic 2: NVT 3:NpT
extern int VType; // volume input type

extern int NumberOfParticles; //N
extern int Nf; // degrees of freedom 3N-3 if velocity center of mass fixed
extern int NumberOfLatticeSites;
extern int NA,NB,NC;
extern double fA,fB,fC;// fraction

//extern int NumberOfCollisions; // number of collisions
//extern int NumberOfSteps; // number of collisions or transer
extern int NumberOfInitialSteps;
extern int SampleMultiplier,MovieMultiplier;
//extern int TimeMultiplier,TimeFrequency;


/**************************************************/
extern double kB; // Boltzmann constant, set to 1
extern double T,T0; // temperature
extern double Tinstant; // temperature
extern double Kinstant,Upotential,Virial;
extern double dT; // temperature increment
extern double beta; // 1/(kB*T)
extern double rho,rhoA,rhoB,rhoC; // number density rho = N/V
extern double packingfraction; // phi = pi/6 rho in 3D
extern double phimax; // final packing fraction
extern double V,Vo;  // Volume
extern double L,Lo; //simulation box length V = L^3
extern double latticeconst; // lattice constant
extern double P; // pressure
/**************************************************/

extern double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
extern double mass[MAX_NumberOfParticles]; // particle mass  m_i
extern double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
extern double sigmaA,sigmaB,sigmaC,sigmaAB,sigmaBC,sigmaCA;
extern double massA,massB,massC;
extern double g[MAX_drBins],gAA[MAX_drBins],gBB[MAX_drBins],gCC[MAX_drBins];
extern double gAB[MAX_drBins],gBC[MAX_drBins],gCA[MAX_drBins];
extern double dradial;
extern int drBins;

typedef struct
{
	double x;
	double y;
	double z;
} VECTOR;

extern VECTOR position[MAX_NumberOfParticles]; // t
extern VECTOR position_old[MAX_NumberOfParticles]; // t-dt
extern VECTOR position_new[MAX_NumberOfParticles]; // t+dt
extern VECTOR velocity[MAX_NumberOfParticles];
extern VECTOR velocity_old[MAX_NumberOfParticles];
extern VECTOR velocity_new[MAX_NumberOfParticles];

void ReadInput(void);
void Initialization(int job);
double BoxMuller(double mm, double ss);

double Distance(int i,int j); // return rij^2
double SigmaIJ(int iID,int jID);
void PBC(double *xx);
void MinimumImage(double *xx);

void Kinetic(void);
void RadialDis(void);
void CollisionInfo(int i);
void Collision(int i,int j);
void Writemovie(FILE *FilePtr);
void MinRatio(void); // update i,j,dmin
void Compress(void); // compress volume

// cell list
int CellDetermine(int i);
void NeighborCell(void); // generate neighbor cell list
void CoordinateTrans(int i);
void MakeCell(void);
void AddToCell(int i,int iCell);
void RemoveFromCell(int i,int iCell);
void UpdateCell(int i);
void EscapeInfo(int i); 
void OverlapCheck(void); 
void CollisionUpdate(int i, int iCell); 
