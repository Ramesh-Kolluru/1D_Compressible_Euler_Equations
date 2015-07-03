#include "1D_Euler_FVM.h"



void Slope_Limiter(vector<double> & r)
{
  double epsilon = 0.0;
    r[0] = (r[0]>1) ? 1 : ( (r[0]<0) ? epsilon :r[0]); 
    r[1] = (r[1]>1) ? 1 : ( (r[1]<0) ? epsilon :r[1]);
    r[2] = (r[2]>1) ? 1 : ( (r[2]<0) ? epsilon :r[2]);
    
//      cout<<r[0]<<"\t"<<r[1]<<"\t"<<r[2]<<endl;

}


void Min_Mod(int & Cell_No,vector<double> & Phi,int & Face_No)
{
  double inv_epsilon = 1e6;
  if(Cell_No == Total_Cells-2)
  {Face_No = 0;}
  switch(Face_No)
  {
    case 0:
        (U[Cell_No + 1][0]-U[Cell_No][0]) == 0 ? Phi[0]= (U[Cell_No][0]-U[Cell_No-1][0])*inv_epsilon: (Phi[0]= (U[Cell_No][0]-U[Cell_No-1][0])/(U[Cell_No + 1][0]-U[Cell_No][0]));
	(U[Cell_No + 1][1]-U[Cell_No][1]) == 0 ? Phi[1]= (U[Cell_No][1]-U[Cell_No-1][1])*inv_epsilon: (Phi[1]= (U[Cell_No][1]-U[Cell_No-1][1])/(U[Cell_No + 1][1]-U[Cell_No][1]));
	(U[Cell_No + 1][2]-U[Cell_No][2]) == 0 ? Phi[2]= (U[Cell_No][2]-U[Cell_No-1][2])*inv_epsilon: (Phi[2]= (U[Cell_No][2]-U[Cell_No-1][2])/(U[Cell_No + 1][2]-U[Cell_No][2]));
      break;
    case 1:
        (U[Cell_No + 2][0]-U[Cell_No+1][0]) == 0 ? Phi[0]= (U[Cell_No + 1][0]-U[Cell_No][0])*inv_epsilon: (Phi[0]= (U[Cell_No+1][0]-U[Cell_No][0])/(U[Cell_No + 2][0]-U[Cell_No+1][0]));
	(U[Cell_No + 2][1]-U[Cell_No + 1][1]) == 0 ? Phi[1]= (U[Cell_No + 1][1]-U[Cell_No][1])*inv_epsilon: (Phi[1]= (U[Cell_No + 1][1]-U[Cell_No][1])/(U[Cell_No + 2][1]-U[Cell_No + 1][1]));
	(U[Cell_No + 2][2]-U[Cell_No + 1][2]) == 0 ? Phi[2]= (U[Cell_No + 1][2]-U[Cell_No][2])*inv_epsilon: (Phi[2]= (U[Cell_No + 1][2]-U[Cell_No][2])/(U[Cell_No + 2][2]-U[Cell_No + 1][2]));
      break;
  }
  if(isnan(Phi[0])|isnan(Phi[1])|isnan(Phi[2]))
  {
    Plot_Solution();
    exit(0);
  }
  Slope_Limiter(Phi);
}


