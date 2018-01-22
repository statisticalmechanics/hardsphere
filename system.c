#include "system.h"
int CellSwitch; // use or not use cell list
struct LinkedList CellTrack[MAX_NumberOfParticles];
int HeadOfChain[MAX_NumberOfCells];
double rcell; // size of cell rcell > sigma
int LCell,NumberOfCells;
int NumberOfNeighborCells;
double rc;// rcutoff
int NeighborCellList[MAX_NumberOfCells][MAX_NumberOfNeighborCells];
int CellCoordinate[MAX_NumberOfCells][8];

double randomseed;

int CompressSwitch;
double collisionfrequency; 
double maxfrequency;
double timewindow;

double compresstime;// compress every such time interval 
double compressfactor;// 0~1
int icompress,jcompress;// particle i,j which gives dmin
double compressratio; // particle size increase ratio
double tolerance,timetol;

double TimeBig; // collision time upper bound
double CollisionTime[MAX_NumberOfParticles];
int CollisionPartner[MAX_NumberOfParticles];
double EscapeTime[MAX_NumberOfParticles]; // time needed to escape a cell
int EscapeID;
//int Escapexyz[MAX_NumberOfParticles]; // 1 x 2 y 3 z
//int OnBoundx[MAX_NumberOfParticles],OnBoundy[MAX_NumberOfParticles],OnBoundz[MAX_NumberOfParticles]; // xyz: 0 off; 1 on boundary
//int Cellold[MAX_NumberOfParticles],Cellnew[MAX_NumberOfParticles];

int JobIndex;
double  rate; // cooling rate

int rInitialType,vInitialType;
int PackType; // SC, BCC or FCC
int PotentialType; //1: shifted-force L-J
int EnsembleType; // 0: NVE 1: isokinetic 2: NVT 3:NpT
int VType; // volume input type

int NumberOfParticles; //N
int Nf; // degrees of freedom 3N-3 if velocity center of mass fixed
int NumberOfLatticeSites;
int NA,NB,NC;
double fA,fB,fC;// fraction

//int NumberOfCollisions; // number of collisions
//int NumberOfSteps; // number of collisions or transer
int NumberOfInitialSteps;
int SampleMultiplier,MovieMultiplier;
//int TimeMultiplier,TimeFrequency;


double kB; // Boltzmann constant, set to 1
double T,T0; // temperature
double Tinstant; // temperature
double Kinstant,Upotential,Virial;
double dT; // temperature increment
double beta; // 1/(kB*T)
double rho,rhoA,rhoB,rhoC; // number density rho = N/V
double packingfraction; // phi = pi/6 rho in 3D
double phimax; // final packing fraction
double V,Vo;  // Volume
double L,Lo; //simulation box length V = L^3
double latticeconst; // lattice constant
double P; // pressure

double sigma[MAX_NumberOfParticles]; // hardcore diameter  sigma_i
double mass[MAX_NumberOfParticles]; // particle mass  m_i
double identity[MAX_NumberOfParticles]; // particle identity 1:A 2:B
double sigmaA,sigmaB,sigmaC,sigmaAB,sigmaBC,sigmaCA;
double massA,massB,massC;
double g[MAX_drBins],gAA[MAX_drBins],gBB[MAX_drBins],gCC[MAX_drBins];
double gAB[MAX_drBins],gBC[MAX_drBins],gCA[MAX_drBins];
double dradial;
int drBins;

VECTOR position[MAX_NumberOfParticles]; // t
VECTOR position_old[MAX_NumberOfParticles]; // t-dt
VECTOR position_new[MAX_NumberOfParticles]; // t+dt
VECTOR velocity[MAX_NumberOfParticles];
VECTOR velocity_old[MAX_NumberOfParticles];
VECTOR velocity_new[MAX_NumberOfParticles];

