/*******************************************************
* divide simulation box into small cells
* of size rcell (rcell > sigma)
*******************************************************/

#include <stdio.h>
#include <math.h>
#include "system.h"

void MakeCell(void)
{
 int i,j;

 LCell = (int)(pow(NumberOfParticles,1.0/3.0));
 rcell = L/LCell;

 if(rcell < rc) // rc is the max of sigmaAA, sigmaBB,sigmaCC, sigmaAB, sigmaBC,sigmaCA
 {
 LCell = (int)(L/rc);
 rcell = L/LCell;
 }
 
 NumberOfNeighborCells = 27;
 NumberOfCells = CUBIC(LCell);
 
 printf("\n");
 printf("rcell = %lf\tLCell = %d\tNCell = %d\n",rcell,LCell,NumberOfCells);
 printf("\n");

 if(CellSwitch == 0)
 printf("cell list is NOT used\n");
 printf("\n");


 NeighborCell(); // assign each cell i of all its neighboring cells

 //group particles into cells

 for(i=0;i<NumberOfCells;i++)
   HeadOfChain[i] = -1;  

 for(i=0;i<NumberOfParticles;i++){

   CellTrack[i].WhichCell = CellDetermine(i);

   CellTrack[i].Next = HeadOfChain[CellTrack[i].WhichCell];

   if(CellTrack[i].Next != -1) // avoid array[-1] which is segmentation  fault
    CellTrack[CellTrack[i].Next].Prev = i;

   HeadOfChain[CellTrack[i].WhichCell] = i;
 }// end loop particle i

return;
}

