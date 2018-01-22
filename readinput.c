/*********************************************************************
 * input simulation parameters from file "input"
 *********************************************************************/
#include <stdio.h>
#include <math.h>
#include "system.h"
#include "ran_uniform.h"

void ReadInput(void)
{
 int i;
 FILE *fp;
 
 fp=fopen("input","r"); 
 fscanf(fp,"%d",&EnsembleType);
 fscanf(fp,"%d %d %d %d",&rInitialType,&PackType,&NumberOfLatticeSites,&vInitialType);
 fscanf(fp,"%d",&PotentialType);
 fscanf(fp,"%d %d %d %d",&NumberOfParticles,&NA,&NB,&NC);
 fscanf(fp,"%lf %lf %lf",&massA,&massB,&massC);
 fscanf(fp,"%lf %lf %lf",&sigmaA,&sigmaB,&sigmaC);
 fscanf(fp,"%lf %lf %lf",&sigmaAB,&sigmaBC,&sigmaCA);
 fscanf(fp,"%lf",&T);
 fscanf(fp,"%lf",&V);
 fscanf(fp,"%d",&NumberOfInitialSteps);
 fscanf(fp,"%d %d",&SampleMultiplier,&MovieMultiplier);
 fscanf(fp,"%d",&drBins);
 fscanf(fp,"%lf",&randomseed);
 fscanf(fp,"%d",&CellSwitch);
 fscanf(fp,"%d %lf",&CompressSwitch,&compressfactor);
 fscanf(fp,"%lf %lf",&timewindow,&maxfrequency);
 fscanf(fp,"%lf",&phimax);
 fscanf(fp,"%d",&VType);
 fclose(fp);

 if(phimax < 0.01) phimax = 1.0;


 rc = MAX(sigmaA,sigmaB);
 rc = MAX(sigmaC,rc);
 rc = MAX(sigmaAB,rc);
 rc = MAX(sigmaBC,rc);
 rc = MAX(sigmaCA,rc);

 TimeBig = 1.0E10;
 tolerance = 1.0E-10;
 timetol = 1.0E-9;

 fA = 1.*NA/NumberOfParticles;
 fB = 1.*NB/NumberOfParticles;
 fC = 1.*NC/NumberOfParticles;

 Nf = 3*NumberOfParticles-3;
// Nf = 3*NumberOfParticles;

 //InitializeRandomNumberGenerator(time(0l)); //time(0L)
 //InitializeRandomNumberGenerator(randomseed); //time(0L)

 //if(PotentialType == 0) printf("Lennard-Jones potential\n");

 if(EnsembleType == 0) printf("NVE ensemble\n");
 if(EnsembleType == 1) printf("isokinetic ensemble by velocity rescaling\n");
 if(EnsembleType == 2) printf("NVT ensemble via Nose-Hoover thermostat\n");
 if(EnsembleType == 3) printf("NVT ensemble via Brown-Clarke\n");

 T0 = T; // initial temperature
 kB = 1.;
 beta = 1.0 / T / kB;

 //T = T + dT*JobIndex;
 //rho = rho + dT*JobIndex;

 rho = NumberOfParticles/V*CUBIC(sigmaA);
 rhoA = NA/V*CUBIC(sigmaA);
 rhoB = NB/V*CUBIC(sigmaA);
 rhoC = NC/V*CUBIC(sigmaA);
 L = pow(V,1.0/3.0); 
 dradial = L/2./drBins; // PBC is used, L/2

 //packingfraction = M_PI/6.*rho*CUBIC(sigmaA);
 packingfraction = M_PI/6.*rho*(fA*CUBIC(sigmaA)+fB*CUBIC(sigmaB)+fC*CUBIC(sigmaC));

 //V = NumberOfParticles / rho;

 printf("\n");
 printf("Number of equilibrium MD steps: %d\n",NumberOfInitialSteps);
 printf("Sample multiplier (sampling frequency): %d\n",SampleMultiplier);
 printf("Movie multiplier (draw snapshot frequency): %d\n",MovieMultiplier);
 printf("\n");

 printf("\n");
 printf("T = %lf\tbeta = %lf\n",T,beta);
 if(VType == 0)
 {
 printf("rho = %lf\n",rho);
 printf("packing fraction = %lf\n",packingfraction);
 printf("V = %lf\n",V);
 printf("L = %lf\n",L);
 }

 printf("\n");
 printf("N = %d\tNA:NB:NC = %d:%d:%d\tfA:fB:fC = %lf:%lf:%lf\n",NumberOfParticles,NA,NB,NC,fA,fB,fC);
 printf("massA = %lf\tmassB = %lf\tmassC = %lf\n",massA,massB,massC);
 printf("sigmaA = %lf\tsigmaB = %lf\tsigmaC = %lf\n",sigmaA,sigmaB,sigmaC);
 printf("sigmaAB = %lf\tsigmaBC = %lf\tsigmaCA = %lf\n",sigmaAB,sigmaBC,sigmaCA);
 printf("\n");

 if(CompressSwitch == 0)
 printf("no compression\n");
 else
 printf("compression is applied\n");

 printf("compression factor = %lf\n",compressfactor);
 printf("compresstion time window = %lf\n",timewindow);
 printf("maximum collision frequency = %lf\n",maxfrequency);
 printf("\n");
 
 return;
}
