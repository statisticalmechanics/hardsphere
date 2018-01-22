/***************************************************
* 3D Molecular Dynamics (MD) simulation of
* Ternary mixtures: particle A B C
* with hard sphere interactions 
* Kai Zhang, Yale University, 2013
****************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "system.h"
#include "ran_uniform.h"

//int main(void)
int main(int argc, char *argv[])
{
 //int step;
 int coln;
 int compressnumber;
 int event; // compress collide transfer
 double time; // real time
 double sampletime; // sampling time 
 double tij; // collision time
 double tescape; // time needs to escape a cell
 int I,J;// collision particle index
 int cello,celln;
 int i;
 char filename[20];
 FILE *fp;
 FILE *fpmovie;
 FILE *fpvelocity;
 FILE *fpsample;
 FILE *fpcompress;
 FILE *fpr; // position
 FILE *fpv; // velocity
 int onboundx,onboundy,onboundz;

 double Sum_K,Sum_U,Sum_E,Sum_W,Sum_T;
 double Count;

 double maxL;

 sscanf(argv[1],"%d",&JobIndex);
 
 printf("**************** 3D Hard Sphere Molecular Dynamics simulation ********************");
 printf("\n");

 ReadInput();

 if(randomseed == 0.)  randomseed = (double)(JobIndex);
 printf("randomseed = %lf\n",randomseed);
 InitializeRandomNumberGenerator(randomseed); //time(0L)

 Initialization(JobIndex);

 if(VType == 1)
 {
  sprintf(filename,"V0_%d",JobIndex);
  fp = fopen(filename,"r");
  fscanf(fp,"compress = %*d\tt = %*lf\ttsample = %*lf\trho = %*lf\tp = %*lf\tcolns = %*d\tfrequency = %*lf\tphi = %*lf\tV = %lf\n",&V);
  V = ((int)(V*1000.)+1)/1000.;
  L = pow(V,1.0/3.0); 
  fclose(fp);
  /*********
  maxL = 0.0;
  for(i=0;i<NumberOfParticles;i++)
  {
   if(position[i].x > maxL) maxL = position[i].x;
   if(position[i].y > maxL) maxL = position[i].y;
   if(position[i].z > maxL) maxL = position[i].z;
  }
  L = maxL;
  V = CUBIC(L);
  ***********/
  rho = NumberOfParticles/V*CUBIC(sigmaA);
  rhoA = NA/V*CUBIC(sigmaA);
  rhoB = NB/V*CUBIC(sigmaA);
  rhoC = NC/V*CUBIC(sigmaA);
  dradial = L/2./drBins; // PBC is used, L/2
  packingfraction = M_PI/6.*rho*(fA*CUBIC(sigmaA)+fB*CUBIC(sigmaB)+fC*CUBIC(sigmaC));
 
  printf("volume input from current configuration\n");
  printf("rho = %lf\n",rho);
  printf("packing fraction = %lf\n",packingfraction);
  printf("V = %lf\n",V);
  printf("L = %lf\n",L);
 }

 OverlapCheck(); 

 MakeCell();

 fflush(stdout);

for(i=0;i<NumberOfParticles;i++)
{
 CollisionTime[i] = TimeBig;
 CollisionPartner[i] = NumberOfParticles-1;
 EscapeTime[i] = TimeBig;
// Escapexyz[i] = 0;
}
tescape = TimeBig;

//initialize collision and cell escaping time
for(i=0;i<NumberOfParticles;i++)
{
 CollisionInfo(i); //return collision time and partner

 if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)//use Cell list
 EscapeInfo(i); // return excape time
}
   
fflush(stdout);

 /***********************************************************/
 sprintf(filename,"position0_%d.dat",JobIndex);
 fpr=fopen(filename,"w");
 sprintf(filename,"velocity0_%d.dat",JobIndex);
 fpv=fopen(filename,"w");
 for(i=0;i<NumberOfParticles;i++)
 {
  fprintf(fpr,"%lf\t%lf\t%lf\n",position[i].x,position[i].y,position[i].z);
  fprintf(fpv,"%lf\t%lf\t%lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
 }
 fclose(fpr);
 fclose(fpv);
 /***********************************************************/

 //sprintf(filename,"movie_%d.xyz",JobIndex);
 //fpmovie=fopen(filename,"w");
// sprintf(filename,"sample_%d.dat",JobIndex);
// fpsample=fopen(filename,"w");
// sprintf(filename,"compress_%d.dat",JobIndex);
// fpcompress=fopen(filename,"w");
// sprintf(filename,"velocity_%d.xyz",JobIndex);
 //fpvelocity=fopen(filename,"w");

 printf("start MD loop ......\n");
 printf("\n");

compressnumber = 0;
collisionfrequency = 0.0;  

 if(MovieMultiplier < 10000) 
 {
  sprintf(filename,"movie_%d.xyz",JobIndex);
  fpmovie=fopen(filename,"a+");
  Writemovie(fpmovie);
  fclose(fpmovie);
 }

 Sum_K = 0.;
 Sum_U = 0.;
 Sum_E = 0.;
 Sum_W = 0.;
 Sum_T = 0.;
 Count = 0.;
 for(i=0;i<drBins;i++)
 {
  g[i] = 0.;
  gAA[i] = 0.;
  gBB[i] = 0.;
  gCC[i] = 0.;
  gAB[i] = 0.;
  gBC[i] = 0.;
  gCA[i] = 0.;
 }
 
time = 0.;
tij = TimeBig;
sampletime = 0.;
coln = 0; // collision numbers in timewindow

 sprintf(filename,"sample_%d.dat",JobIndex);
 fpsample=fopen(filename,"a+");
   fprintf(fpsample,"collision = %d\tdt = %lf\tt = %lf\tK = %lf\tVirial = %lf\tTinstant = %lf\tV = %lf\tphi = %lf\tncell = %d\n", \
   coln,tij,time,Kinstant/NumberOfParticles,Virial/V/3.,Tinstant,V,packingfraction,NumberOfCells);
 fclose(fpsample);


while(collisionfrequency < maxfrequency && packingfraction < phimax) // collision frequency is below threshod
{
 Sum_K = 0.;
 Sum_U = 0.;
 Sum_E = 0.;
 Sum_W = 0.;
 Sum_T = 0.;
 Count = 0.;
 for(i=0;i<drBins;i++)
 {
  g[i] = 0.;
  gAA[i] = 0.;
  gBB[i] = 0.;
  gCC[i] = 0.;
  gAB[i] = 0.;
  gBC[i] = 0.;
  gCA[i] = 0.;
 }
 
time = 0.;
sampletime = 0.;
coln = 0; // collision numbers in timewindow
event = 100;

while(1) // relaxation loop
{
  // next collision
  /*********************************/
   tij = TimeBig;
   for(i=0;i<NumberOfParticles;i++)
   {
    if(CollisionTime[i] < tij)
    {
     tij = CollisionTime[i];
     I = i; // collision particle 1
    }
   }
   J = CollisionPartner[I]; // collision particle 2
  /*********************************/
   
  // next transfer
  /*********************************/
  if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells) // if use cell list 
  {
   tescape = TimeBig;
   for(i=0;i<NumberOfParticles;i++)
   {
    if(EscapeTime[i] < tescape)
    {
     tescape = EscapeTime[i];
     EscapeID = i;
    }
  
   }
  }
/*********************************/
  
if(time+MIN(tij,tescape) > timewindow) // evolve dt/2, prepare for compression
{
  event = 2;
  if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells) tij =  MIN(tij,tescape);
  time += tij/2.; // evolve dt/2 so that not particle on contact
  collisionfrequency = coln/time; // collision frequency in timewindow

  for(i=0;i<NumberOfParticles;i++)
  {
    position[i].x = position[i].x + velocity[i].x * tij/2.0;
    position[i].y = position[i].y + velocity[i].y * tij/2.0;
    position[i].z = position[i].z + velocity[i].z * tij/2.0;
    PBC(&(position[i].x)); // periodic boundary condition
    PBC(&(position[i].y)); // periodic boundary condition
    PBC(&(position[i].z)); // periodic boundary condition
    CollisionTime[i] = CollisionTime[i] - tij/2.0;
    
    if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)
     EscapeTime[i] = EscapeTime[i] - tij/2.0;
  }
  
  break; // getting out of relaxation loop
}
/************  transfer or collision ********************/
else // do not compress. transfer or collision
{
  if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells && tescape < tij) // if transfer 
  {
   event = 1;
   tij = tescape;
   if(coln>=NumberOfInitialSteps) sampletime += tij;
   time += tij; // increase time
   for(i=0;i<NumberOfParticles;i++)
   {
    cello = CellTrack[i].WhichCell;

    position[i].x = position[i].x + velocity[i].x * tij;
    position[i].y = position[i].y + velocity[i].y * tij;
    position[i].z = position[i].z + velocity[i].z * tij;
    PBC(&(position[i].x)); // periodic boundary condition
    PBC(&(position[i].y)); // periodic boundary condition
    PBC(&(position[i].z)); // periodic boundary condition
 
    onboundx = 0;
    onboundy = 0;
    onboundz = 0;
 
    if((fabs(position[i].x/rcell - (double)((int)(position[i].x/rcell))) < tolerance) && (velocity[i].x < 0.) \
     ||(fabs(position[i].x/rcell - (double)((int)(position[i].x/rcell))+1) < tolerance) && (velocity[i].x > 0.))
    onboundx = 1;
    
    if((fabs(position[i].y/rcell - (double)((int)(position[i].y/rcell))) < tolerance) && (velocity[i].y < 0.) \
     ||(fabs(position[i].y/rcell - (double)((int)(position[i].y/rcell))+1) < tolerance) && (velocity[i].y > 0.))
    onboundy = 1;
    
    if((fabs(position[i].z/rcell - (double)((int)(position[i].z/rcell))) < tolerance) && (velocity[i].z < 0.) \
     ||(fabs(position[i].z/rcell - (double)((int)(position[i].z/rcell))+1) < tolerance) &&(velocity[i].z > 0.))
    onboundz = 1;
   
    if(i == EscapeID || onboundx==1 || onboundy==1 || onboundz==1  ) 
    {
     if(onboundx == 1 ) position[i].x = round(position[i].x/rcell)*rcell;
     if(onboundy == 1 ) position[i].y = round(position[i].y/rcell)*rcell;
     if(onboundz == 1 ) position[i].z = round(position[i].z/rcell)*rcell;

     celln = CellDetermine(i);

     RemoveFromCell(i,cello);
     AddToCell(i,celln);
     CellTrack[i].WhichCell = celln;
     
     EscapeInfo(i);

    }// end if particle on boundary
    else
     EscapeTime[i] = EscapeTime[i] - tij; //only for particles not on boundary

    CollisionTime[i] = CollisionTime[i] - tij;
   }//end loop i

   CollisionUpdate(EscapeID,celln);
  }// end if transfer
  else // collision
  {
   event = 0;
   time += tij; // increase time
   coln++;
 
   for(i=0;i<NumberOfParticles;i++)
   {
    CollisionTime[i] = CollisionTime[i] - tij;
    position[i].x = position[i].x + velocity[i].x * tij;
    position[i].y = position[i].y + velocity[i].y * tij;
    position[i].z = position[i].z + velocity[i].z * tij;
    PBC(&(position[i].x)); // periodic boundary condition
    PBC(&(position[i].y)); // periodic boundary condition
    PBC(&(position[i].z)); // periodic boundary condition
    
    if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)
     EscapeTime[i] = EscapeTime[i] - tij;
   }

   Collision(I,J); // collision dynamics

   if(coln>=NumberOfInitialSteps)
   {
    Sum_W += Virial;
    sampletime += tij;
   }

/***************************update*****************************/
  if(CellSwitch == 1 && NumberOfCells > NumberOfNeighborCells)//use Cell list
  {
    EscapeInfo(I);
    EscapeInfo(J);
   for(i=0;i<NumberOfParticles;i++)
   {
    if(i == I || CollisionPartner[i] == I || i == J || CollisionPartner[i] == J)
     CollisionInfo(i); // O(27*Nc), Nc~1
   }
  }
  else // not use cell list
  {
  for(i=0;i<NumberOfParticles;i++)
  {
    if(i == I || CollisionPartner[i] == I || i == J || CollisionPartner[i] == J)
    CollisionInfo(i); // O(N) computation
  }
  }
/************************end update****************************/
  } // end if collision
  
 fflush(stdout);
} // end if no compress
  
 /************** sample *******************/
  if(event == 0 && coln%SampleMultiplier==0)
  {
   Kinetic();
   
   if(CompressSwitch == 0)
   {
   sprintf(filename,"sample_%d.dat",JobIndex);
   fpsample=fopen(filename,"a+");
   fprintf(fpsample,"collision = %d\tdt = %lf\tt = %lf\tK = %lf\tVirial = %lf\tTinstant = %lf\tV = %lf\tphi = %lf\tncell = %d\n", \
   coln,tij,time,Kinstant/NumberOfParticles,Virial/V/3.,Tinstant,V,packingfraction,NumberOfCells);
   fclose(fpsample);
   }
   fflush(stdout);

  if(coln>=NumberOfInitialSteps)
  {
   Sum_K += Kinstant;
   Sum_T += Tinstant;
   Count += 1.;
   RadialDis(); // g(r)
  }//endif after equilibration
  }//endif every ? steps
  /************** end of sample *******************/

 if(event == 0 && coln%MovieMultiplier==0 && CompressSwitch == 0) 
 {
  sprintf(filename,"movie_%d.xyz",JobIndex);
  fpmovie=fopen(filename,"a+");
  Writemovie(fpmovie);
  fclose(fpmovie);
 }

} // end relaxation while loop
 

 Sum_K /= Count;
 Sum_U /= Count;
 Sum_E /= Count;
 Sum_T /= Count;
 Sum_W /= sampletime;
 Sum_W /= 3.;
 
 sprintf(filename,"compress_%d.dat",JobIndex);
 fpcompress=fopen(filename,"a+");
 fprintf(fpcompress,"compress = %d\tt = %lf\ttsample = %lf\trho = %lf\tp = %lf\tcolns = %d\tfrequency = %lf\tphi = %lf\tV = %lf\n", \
 compressnumber,compressnumber*timewindow+time,sampletime,rho,rho*kB*Sum_T+Sum_W/V,coln,collisionfrequency,packingfraction,V);
 fclose(fpcompress);
 
 /***********************************************************/
 sprintf(filename,"position1_%d.dat",JobIndex);
 fpr=fopen(filename,"w");
 sprintf(filename,"velocity1_%d.dat",JobIndex);
 fpv=fopen(filename,"w");
 for(i=0;i<NumberOfParticles;i++)
 {
  fprintf(fpr,"%lf\t%lf\t%lf\n",position[i].x,position[i].y,position[i].z);
  fprintf(fpv,"%lf\t%lf\t%lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
 }
 fclose(fpr);
 fclose(fpv);
 /***********************************************************/
   
 if(CompressSwitch == 0)
 break;
  
/*************** compression ******************/
  MinRatio(); // O(N2) or O(N*27*Nc), Nc~1
  Compress(); // O(N)

  compressnumber++;

  for(i=0;i<NumberOfParticles;i++)
  {
    CollisionInfo(i);
   if(CellSwitch == 1 && NumberOfCells>NumberOfNeighborCells) // if use cell list 
    EscapeInfo(i);
  }
   
 if(CompressSwitch == 1 && compressnumber%MovieMultiplier==0) 
 {
  sprintf(filename,"movie_%d.xyz",JobIndex);
  fpmovie=fopen(filename,"a+");
  Writemovie(fpmovie);
  fclose(fpmovie);
 }
 fflush(stdout);
/************end of compression********************/


}// end while collision frequency below threshod

  
 sprintf(filename,"movie_%d.xyz",JobIndex);
 fpmovie=fopen(filename,"a+");
 Writemovie(fpmovie);
 fclose(fpmovie);

 //fclose(fpmovie);
 //fclose(fpsample);
 //fclose(fpcompress);
 //fclose(fpvelocity);


 /***********************************************************/
 sprintf(filename,"position1_%d.dat",JobIndex);
 fpr=fopen(filename,"w");
 sprintf(filename,"velocity1_%d.dat",JobIndex);
 fpv=fopen(filename,"w");
 for(i=0;i<NumberOfParticles;i++)
 {
  fprintf(fpr,"%lf\t%lf\t%lf\n",position[i].x,position[i].y,position[i].z);
  fprintf(fpv,"%lf\t%lf\t%lf\n",velocity[i].x,velocity[i].y,velocity[i].z);
 }
 fclose(fpr);
 fclose(fpv);
 /***********************************************************/
 sprintf(filename,"gr_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  g[i] /= Count;
  g[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*rho;
  fprintf(fp,"r = %lf\tg(r) = %lf\n",(i+0.5)*dradial,g[i]/NumberOfParticles);
 }
 fclose(fp);
 
 sprintf(filename,"grAA_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gAA[i] /= Count;
  gAA[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*rhoA;
  fprintf(fp,"r = %lf\tgAA(r) = %lf\n",(i+0.5)*dradial,gAA[i]/NA);
 }
 fclose(fp);

 sprintf(filename,"grBB_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gBB[i] /= Count;
  gBB[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*rhoB;
  fprintf(fp,"r = %lf\tgBB(r) = %lf\n",(i+0.5)*dradial,gBB[i]/NB);
 }
 fclose(fp);
 
 sprintf(filename,"grCC_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gCC[i] /= Count;
  gCC[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial)*rhoC;
  fprintf(fp,"r = %lf\tgCC(r) = %lf\n",(i+0.5)*dradial,gCC[i]/NC);
 }
 fclose(fp);

 sprintf(filename,"grAB_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gAB[i] /= Count;
  gAB[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial);
  fprintf(fp,"r = %lf\tgAB(r) = %lf\n",(i+0.5)*dradial,gAB[i]*V/(NA*NB)/2.);
 }
 fclose(fp);
 
 sprintf(filename,"grBC_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gBC[i] /= Count;
  gBC[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial);
  fprintf(fp,"r = %lf\tgBC(r) = %lf\n",(i+0.5)*dradial,gBC[i]*V/(NB*NC)/2.);
 }
 fclose(fp);
 
 sprintf(filename,"grCA_%d.dat",JobIndex);
 fp=fopen(filename,"w");
 for(i=0;i<drBins;i++)
 {
  gCA[i] /= Count;
  gCA[i] /= 4./3.*M_PI*(CUBIC(i+1)-CUBIC(i))*CUBIC(dradial);
  fprintf(fp,"r = %lf\tgCA(r) = %lf\n",(i+0.5)*dradial,gCA[i]*V/(NC*NA)/2.);
 }
 fclose(fp);
 
 printf("\n");
 printf("MD simulation is finished\n");
 printf("final packing fraction is %lf\n",packingfraction);
 
 printf("\n");

 printf("****************************** the end *******************************");


 return 0;
}
