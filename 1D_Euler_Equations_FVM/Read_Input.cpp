#include "1D_Euler_FVM.h"

void Read_Exact_Solution_Input(string & File_Input)
{
  int Cell_Index=0,number=0;
  double delx=0.0,density =0.0,P=0.0,v=0.0,e=0.0;
 
  ifstream ipfile(File_Input.c_str(),ios::in);
 
  if(ipfile.is_open())
  {
   do
   {
     ipfile>>delx>>density>>v>>P>>e;
     
     if(number%10==0)
     {
       Exact_Density[Cell_Index] = density;
       Exact_Pressure[Cell_Index] = P;
       Exact_Velocity[Cell_Index] = v;
       Exact_e[Cell_Index] = e;
       Cell_Index++;
       number++;
//        cout<<Cell_Index<<"\t"<<delx<<"\t"<<density<<"\t"<<P<<"\t"<<v<<"\t"<<e<<endl;
     }
     else
     {
      number++; 
      }
   }while(!ipfile.eof());
   }
  else
  {
   cout<<"Could not open file to read......... Please check the file path\n";
  }
//   cout<<Exact_Density.size()<<endl;
}


