/*************************************
* Add Particle i to Cell iCell by 
* making i Head of Chain
*************************************/
#include <stdio.h>
#include "system.h"

void AddToCell(int i,int iCell)
{
     CellTrack[i].Next = HeadOfChain[iCell];
     if(CellTrack[i].Next!=-1) // avoid array[-1] which is segmentation 
        CellTrack[CellTrack[i].Next].Prev = i;
     HeadOfChain[iCell] = i;

 return;
}
/*************************************
* Remove Particle i from Cell iCell
*************************************/

void RemoveFromCell(int i,int iCell)
{

     if(i == HeadOfChain[iCell]) // if i is the first in the cell
      HeadOfChain[iCell] = CellTrack[i].Next;
     else // i is NOT the first in the cell eg, ? -> i -> ? -> -1
     {
      CellTrack[CellTrack[i].Prev].Next = CellTrack[i].Next;
      if(CellTrack[i].Next != -1) // if i is not last in the old cell
      CellTrack[CellTrack[i].Next].Prev = CellTrack[i].Prev;
     }

 return;

}

/*************************************
* check and update the cell of particle i
*************************************/

void UpdateCell(int i)
{
 int OldCell,NewCell;

 OldCell = CellTrack[i].WhichCell;
 NewCell = CellDetermine(i);

 if(OldCell != NewCell)
 {
  CellTrack[i].WhichCell = NewCell;
  RemoveFromCell(i,OldCell);
  AddToCell(i,NewCell);
 }
 return;
}
