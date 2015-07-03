#include "1D_Euler_FVM.h"

void Plot_Solution()
{
	  Gnuplot g1=Gnuplot("linespoints");
	  Gnuplot g2=Gnuplot("linespoints");
          Gnuplot g3=Gnuplot("linespoints");
  	  Gnuplot g4=Gnuplot("linespoints");
	  
	  g1.cmd("set terminal png size 1280,960 enhanced");
	  g1.cmd(Plot_Density.c_str());
 	  g1.cmd("set multiplot");
	  g1.cmd("unset key");
	  g1.cmd("set origin 0.0,0.0");
	  g1.plot_xy(Cell_del_X,Density,"Density");
	  g1.cmd("clear");
	  g1.plot_xy(Cell_del_X,Exact_Density,"Exact_Solution");

	  g2.cmd("set terminal png size 1280,960 enhanced");
	  g2.cmd(Plot_Pressure.c_str());
   	  g2.cmd("set multiplot");
	  g2.cmd("unset key");
	  g2.cmd("set origin 0.0,0.0");
	  g2.plot_xy(Cell_del_X,Pressure,"Pressure");
	  g2.cmd("clear");
	  g2.plot_xy(Cell_del_X,Exact_Pressure,"Exact Solution for Pressure");
	  
	  g3.cmd("set terminal png size 1280,960 enhanced");
	  g3.cmd(Plot_Velocity.c_str());
   	  g3.cmd("set multiplot");
	  g3.cmd("unset key");
	  g3.cmd("set origin 0.0,0.0");
          g3.plot_xy(Cell_del_X,Velocity,"Velocity");
	  g3.cmd("clear");
	  g3.plot_xy(Cell_del_X,Exact_Velocity,"Exact Solution for Velocity");
	  
	  g4.cmd("set terminal png size 1280,960 enhanced");
	  g4.cmd(Plot_e.c_str());
   	  g4.cmd("set multiplot");
	  g4.cmd("unset key");
	  g4.cmd("set origin 0.0,0.0");
          g4.plot_xy(Cell_del_X,Internal_Energy,"Internal Energy");
	  g4.cmd("clear");
	  g4.plot_xy(Cell_del_X,Exact_e,"Exact Solution for Internal_Energy");
   	  
	  
	  sleep(sleeplength);
	  
	  g1.reset_plot();
	  g2.reset_plot();
	  g3.reset_plot();
	  g4.reset_plot();
}

void Plot_Solution(vector< vector<double>  > & u)
{
}
