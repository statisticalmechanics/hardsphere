/***************************************
 * periodic boundary condition
 ***************************************/

#include <stdio.h>
#include <stdlib.h>
#include "system.h"
#include <math.h>

/************************/
void PBC(double *xx)
{
  if(*xx < 0.0)
   *xx += L;
  else if(*xx >= L)
   *xx -= L; 

 return;
}
/************************/

void MinimumImage(double *xx)
{
 *xx = *xx - L*round(*xx/L);

 return;
}

/***************
void PBC(void)
{
    if(position[i].x < 0.0)
     position[i].x += L;
    else if(position[i].x >= L)
     position[i].x -= L; 
    
    if(position[i].y < 0.0)
     position[i].y += L;
    else if(position[i].y >= L)
     position[i].y -= L; 
    
    if(position[i].z < 0.0)
     position[i].z += L;
    else if(position[i].z >= L)
     position[i].z -= L; 

return;
}
***************/
