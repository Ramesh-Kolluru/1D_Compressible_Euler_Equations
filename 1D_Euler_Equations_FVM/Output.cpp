#include "1D_Euler_FVM.h"


void Write_Solution(string & File_Name)
{
  ofstream outputfile(File_Name.c_str(),ios::out);
  if(outputfile.is_open())
  {
     outputfile<<"Domain Length"<<"\t"<<"Pressure"<<"\t"<<"Density"<<"\t"<<"Velocity"<<"\t"<<"Temperature"<<"\t"<<"Sound Speed"<<endl;
    for(int Cell_Index=0;Cell_Index<Total_Cells;Cell_Index++)
    {
      outputfile<<Cell_del_X[Cell_Index]<<"\t"<<Pressure[Cell_Index]<<"\t"<<Density[Cell_Index]<<"\t"<<Velocity[Cell_Index]<<"\t"<<Temperature[Cell_Index]<<"\t"<<C[Cell_Index]<<"\t"<<Internal_Energy[Cell_Index]<<endl;
    }
  }
  else
  {
    cout<<"Could not open file to write output.........Please check the path"<<endl;
  }
}





void Directory_Name()
{
  switch(Dissipation_Type)
  {
    case 0:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/LLF/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/LLF/";
      break;
   case 1:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/PVU/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/PVU/";
      break;
   case 2:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/PVU_P/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/PVU_P/";
      break;
  case 3:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/MOVERS/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/MOVERS/";
      break;
  case 4:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/ROE/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/ROE/";
      break;
  case 5:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/KFVS/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/KFVS/";
      break;
  case 6:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Steger_Warming/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Steger_Warming/";
  case 7:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Vanleer/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Vanleer/";
      break;
  case 8:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/AUSM/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/AUSM/";
      break;
   case 9:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/ROE/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/ROE/";
      break;
   case 10:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/HLL/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/HLL/";
      break;
   case 11:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/HLLC/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/HLLC/";
      break;
   case 12:
        Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Rusanov/";
	Initial_Solution_Op_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Rusanov/";
      break;


  }
}


void File_Name()
{
  switch(Dissipation_Type)
  {
    case 0:
        Solution_Op_File +="LLF";
	Initial_Solution_Op_File +="Initial_LLF";
	Plot_Density = "set output '../1D_Euler_Solutions/LLF/";
	Plot_Pressure = "set output '../1D_Euler_Solutions/LLF/";
	Plot_Velocity = "set output '../1D_Euler_Solutions/LLF/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/LLF/";
	Plot_e = 	"set output '../1D_Euler_Solutions/LLF/";
      break;
   case 1:
        Solution_Op_File +="PVU";
	Initial_Solution_Op_File +="Initial_PVU";
	Plot_Density = "set output '../1D_Euler_Solutions/PVU/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/PVU/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/PVU/";
	Plot_e = "set output '../1D_Euler_Solutions/PVU/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/PVU/";
      break;
   case 2:
        Solution_Op_File +="PVU_P";
	Initial_Solution_Op_File +="Initial_PVU_P";
	Plot_Pressure ="set output '../1D_Euler_Solutions/PVU_P/";
	Plot_Density = "set output '../1D_Euler_Solutions/PVU_P/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/PVU_P/";
	Plot_e = "set output '../1D_Euler_Solutions/PVU_P/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/PVU_P/";

      break;
  case 3:
        Solution_Op_File +="MOVERS";
	Initial_Solution_Op_File +="Initial_MOVERS";
	Plot_Density = "set output '../1D_Euler_Solutions/MOVERS/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/MOVERS/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/MOVERS/";
	Plot_e = "set output '../1D_Euler_Solutions/MOVERS/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/MOVERS/";

      break;
  case 4:
        Solution_Op_File +="ROE";
	Initial_Solution_Op_File +="Initial_ROE";
	Plot_Density = "set output '../1D_Euler_Solutions/ROE/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/ROE/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/ROE/";
	Plot_e = "set output '../1D_Euler_Solutions/ROE/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/ROE/";
      break;
  case 5:
        Solution_Op_File +="KFVS";
	Initial_Solution_Op_File +="Initial_KFVS";
	Plot_Density = "set output '../1D_Euler_Solutions/KFVS/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/KFVS/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/KFVS/";
	Plot_e = "set output '../1D_Euler_Solutions/KFVS/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/KFVS/";
      break;
      case 6:
        Solution_Op_File +="Steger_Warming";
	Initial_Solution_Op_File +="Initial_Steger_Warming";
	Plot_Density = "set output '../1D_Euler_Solutions/Steger_Warming/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/Steger_Warming/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/Steger_Warming/";
	Plot_e = "set output '../1D_Euler_Solutions/Steger_Warming/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/Steger_Warming/";
      break;
      case 7:
        Solution_Op_File +="Vanleer";
	Initial_Solution_Op_File +="Initial_Vanleer";
	Plot_Density = "set output '../1D_Euler_Solutions/Vanleer/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/Vanleer/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/Vanleer/";
	Plot_e = "set output '../1D_Euler_Solutions/Vanleer/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/Vanleer/";
      break;
       case 8:
        Solution_Op_File +="AUSM";
	Initial_Solution_Op_File +="Initial_AUSM";
	Plot_Density = "set output '../1D_Euler_Solutions/AUSM/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/AUSM/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/AUSM/";
	Plot_e = "set output '../1D_Euler_Solutions/AUSM/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/AUSM/";
      break;
       case 9:
        Solution_Op_File +="ROE";
	Initial_Solution_Op_File +="Initial_ROE";
	Plot_Density = "set output '../1D_Euler_Solutions/ROE/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/ROE/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/ROE/";
	Plot_e = "set output '../1D_Euler_Solutions/ROE/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/ROE/";
      break;
  case 10:
        Solution_Op_File +="HLL";
	Initial_Solution_Op_File +="Initial_HLL";
	Plot_Density = "set output '../1D_Euler_Solutions/HLL/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/HLL/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/HLL/";
	Plot_e = "set output '../1D_Euler_Solutions/HLL/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/HLL/";
      break;
  case 11:
        Solution_Op_File +="HLLC";
	Initial_Solution_Op_File +="Initial_HLLC";
	Plot_Density = "set output '../1D_Euler_Solutions/HLLC/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/HLLC/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/HLLC/";
	Plot_e = "set output '../1D_Euler_Solutions/HLLC/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/HLLC/";
      break;
 case 12:
        Solution_Op_File +="Rusanov";
	Initial_Solution_Op_File +="Initial_Rusanov";
	Plot_Density = "set output '../1D_Euler_Solutions/Rusanov/";
	Plot_Pressure ="set output '../1D_Euler_Solutions/Rusanov/";
	Plot_Velocity ="set output '../1D_Euler_Solutions/Rusanov/";
	Plot_e = "set output '../1D_Euler_Solutions/Rusanov/";
	Plot_Temperature ="set output '../1D_Euler_Solutions/Rusanov/";
      break;

      
  }
  
  
  switch(Test_Case)
  {
    case 1:
        Solution_Op_File +="_Toro_Test_Case_1.txt";
	Initial_Solution_Op_File +="Initial_Toro_Test_Case_1.txt";
	Plot_Density +="Density_Toro_Test_Case_1.png'";
	Plot_Pressure +="Pressure_Toro_Test_Case_1.png'";
	Plot_Velocity +="Velocity_Toro_Test_Case_1.png'";
	Plot_Temperature +="Temperature_Toro_Test_Case_1.png'";
	Plot_e +="Internal_Energy_Toro_Test_Case_1.png'";
	Exact_Solution_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Exact_Solutions/prob1.dat";
      break;
   case 2:
        Solution_Op_File +="_Toro_Test_Case_2.txt";
	Initial_Solution_Op_File +="Initial_Toro_Test_Case_2.txt";
	Plot_Density +="Density_Toro_Test_Case_2.png'";
	Plot_Pressure +="Pressure_Toro_Test_Case_2.png'";
	Plot_Velocity +="Velocity_Toro_Test_Case_2.png'";
	Plot_Temperature +="Temperature_Toro_Test_Case_2.png'";
	Plot_e +="Internal_Energy_Toro_Test_Case_2.png'";
	Exact_Solution_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Exact_Solutions/prob2.dat";

      break;
   case 3:
        Solution_Op_File +="_Toro_Test_Case_3.txt";
	Initial_Solution_Op_File +="Initial_Toro_Test_Case_3.txt";
	Plot_Density +="Density_Toro_Test_Case_3.png'";
	Plot_Pressure +="Pressure_Toro_Test_Case_3.png'";
	Plot_Velocity +="Velocity_Toro_Test_Case_3.png'";
	Plot_Temperature +="Temperature_Toro_Test_Case_3.png'";
	Plot_e +="Internal_Energy_Toro_Test_Case_3.png'";
	Exact_Solution_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Exact_Solutions/prob3.dat";

	break;
  case 4:
        Solution_Op_File +="_Toro_Test_Case_4.txt";
	Initial_Solution_Op_File +="Initial_Toro_Test_Case_4.txt";
	Plot_Density +="Density_Toro_Test_Case_4.png'";
	Plot_Pressure +="Pressure_Toro_Test_Case_4.png'";
	Plot_Velocity +="Velocity_Toro_Test_Case_4.png'";
	Plot_Temperature +="Temperature_Toro_Test_Case_4.png'";
	Plot_e +="Internal_Energy_Toro_Test_Case_4.png'";
	Exact_Solution_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Exact_Solutions/prob4.dat";

	break;
  case 5:
        Solution_Op_File +="_Toro_Test_Case_5.txt";
	Initial_Solution_Op_File +="Initial_Toro_Test_Case_5.txt";
	Plot_Density +="Density_Toro_Test_Case_5.png'";
	Plot_Pressure +="Pressure_Toro_Test_Case_5.png'";
	Plot_Velocity +="Velocity_Toro_Test_Case_5.png'";
	Plot_Temperature +="Temperature_Toro_Test_Case_5.png'";
	Plot_e +="Internal_Energy_Toro_Test_Case_5.png'";
	Exact_Solution_File ="/home/ramesh/Desktop/My_Codes/1D_Codes/1D_Euler_Solutions/Exact_Solutions/prob5.dat";

	break;
  case 6:
        Solution_Op_File +="_Steady_Shock.txt";
	Initial_Solution_Op_File +="Initial_Steady_Shock.txt";
	Plot_Density +="Density_Steady_Shock.png'";
	Plot_Pressure +="Pressure_Steady_Shock.png'";
	Plot_Velocity +="Velocity_Steady_Shock.png'";
	Plot_Temperature +="Temperature_Steady_Shock.png'";
	Plot_e +="Internal_Energy_Steady_Shock.png'";

	break;
   case 7:
        Solution_Op_File +="_Steady_Contact.txt";
	Initial_Solution_Op_File +="Initial_Stead_Contact.txt";
	Plot_Density +="Density_Steady_Contact.png'";
	Plot_Pressure +="Pressure_Steady_Contact.png'";
	Plot_Velocity +="Velocity_Steady_Contact.png'";
	Plot_Temperature +="Temperature_Steady_Contact.png'";
	Plot_e +="Internal_Energy_Steady_Contact.png'";

	break;
   case 8:
        Solution_Op_File +="_SOD_Test_Case_1.txt";
	Initial_Solution_Op_File +="Initial_SOD_Test_Case_1.txt";
	Plot_Density +="Density_SoD_Test_Case_1.png'";
	Plot_Pressure +="Pressure_SoD_Test_Case_1.png'";
	Plot_Velocity +="Velocity_SoD_Test_Case_1.png'";
	Plot_Temperature +="Temperature_SoD_Test_Case_1.png'";
	Plot_e +="Internal_Energy_SoD_Test_Case_1.png'";

	break;
  case 9:
        Solution_Op_File +="_SOD_Test_Case_2.txt";
	Initial_Solution_Op_File +="Initial_SOD_Test_Case_2.txt";
	Plot_Density +="Density_SoD_Test_Case_2.png'";
	Plot_Pressure +="Pressure_SoD_Test_Case_2.png'";
	Plot_Velocity +="Velocity_SoD_Test_Case_2.png'";
	Plot_Temperature +="Temperature_SoD_Test_Case_2.png'";
	Plot_e +="Internal_Energy_SoD_Test_Case_1.png'";
	break;
  }
}