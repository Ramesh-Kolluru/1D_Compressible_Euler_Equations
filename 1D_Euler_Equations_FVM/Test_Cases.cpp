#include "1D_Euler_FVM.h"


void Sod_Test_Case_1()
{
    Rho_L = 1.0;
    u_L   = 0.0;
    P_L   = 100000.0;
    
    Rho_R = 0.125;
    u_R   = 0.0;
    P_R   = 10000.0;
    
  Terminating_Time = 0.035;

    
    for(int Cell_Index =0;Cell_Index<0.5*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.5*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    Directory_Name();
    File_Name();

   Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

void Sod_Test_Case_2()
{
    Rho_L = 1.0;
    u_L   = 0.0;
    P_L   = 100000.0;
    
    Rho_R = 0.010;
    u_R   = 0.0;
    P_R   = 1000.0;

    Terminating_Time = 0.035;
    
    for(int Cell_Index =0;Cell_Index<0.5*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.5*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }

    Directory_Name();
    File_Name();

   Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

// For each test we select a convenient position x0 of the initial discontinuity and an output time. These are stated in the legend of each figure displaying computational results.

void Toro_Test_Case_1()
{

  // Mild test case Solution consists of Left Rarefraction, a contact and Right shock time to stopped = t = 0.25 unitsTest 1 is a modified version of the popular Sod’s test [453]; the solution consists of a right shock wave, a right travelling contact wave and a left sonic rarefaction wave; this test is very useful in assessing the entropy satisfaction property of numerical methods. 

  // x_0 = 0.3, t = 0.2 units
  Rho_L = 1.0;
    u_L   = 0.75;
    P_L   = 1.0;
    
    Rho_R = 0.125;
    u_R   = 0.0;
    P_R   = 0.1;
    
    Terminating_Time = 0.2;
    
    for(int Cell_Index =0;Cell_Index<0.3*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.3*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }

    Directory_Name();
    File_Name();

    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

void Toro_Test_Case_2()
{
  
  //Test 2 has solution consisting of two symmetric rarefaction waves and a trivial contact wave of zero speed; the Star Region between the non–linear waves is close to vacuum, which makes this problem a suitable test for assessing the performance of numerical methods for low–density flows; this is the so called 123 problem introduced in chapter

    Rho_L = 1.0;
    u_L   = -2.0;
    P_L   = 0.4;
    
    Rho_R = 1.0;
    u_R   = 2.0;
    P_R   = 0.4;
    
    // x_0 = 0.5, t = 0.15 units
    
    Terminating_Time = 0.15;
    
    for(int Cell_Index =0;Cell_Index<0.5*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.5*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    Directory_Name();
    File_Name();

    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

void Toro_Test_Case_3()
{
  
  // Chap. 4. Test 3 is designed to assess the robustness a,Exact_Solution_Filend accuracy of numerical methods; its solution consists of a strong shock wave, a contact surface and a left rarefaction wave. 

    Rho_L = 1.0;
    u_L   = 0.0;
    P_L   = 1000.0;
    
    Rho_R = 1.0;
    u_R   = 0.0;
    P_R   = 0.01;
    
    // x_0 = 0.5, t = 0.012 units
    
    Terminating_Time = 0.012;
    
    for(int Cell_Index =0;Cell_Index<0.5*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.5*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    Directory_Name();
    File_Name();

    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

void Toro_Test_Case_4()
{
  
  // Test 4 is also designed to test robustness of numerical methods; the solution consists of three strong discontinuities travelling to the right. See Sect. 4.3.3 of Chap. 4 for more details on the exact solution of these test problems.

    Rho_L = 5.99924;
    u_L   = 19.5975;
    P_L   = 460.894;
    
    Rho_R = 5.99242;
    u_R   = -6.19633;
    P_R   = 46.0950;
    
    // x_0 = 0.4, t = 0.035 units
    Terminating_Time = 0.035;
    
    for(int Cell_Index =0;Cell_Index<0.4*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.4*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
  
    Directory_Name();
    File_Name();


    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

void Toro_Test_Case_5()
{
  
  // Test 5 is also designed to test the robustness of numerical methods but the main reason for devising this test is to assess the ability of the numerical methods to resolve slowly– moving contact discontinuities. The exact solution of Test 5 consists of a left rarefaction wave, a right–travelling shock wave and a stationary contact discontinuity.

    Rho_L = 1.0;
    u_L   = -19.59745;
    P_L   = 1000.0;
    
    Rho_R = 1.0;
    u_R   = -19.59745;
    P_R   = 0.01;
    // x_0 = 0.8, t = 0.012 units
    
    Terminating_Time = 0.012;
    for(int Cell_Index =0;Cell_Index<0.8*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }
    
    for(int Cell_Index =0.8*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
    }

    Directory_Name();
    File_Name();


    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}


void Steady_Shock()
{
  double M_inf  = 2.0;
    Rho_L = 1.0;
    u_L   = 1.0;
    P_L   = 1.0/(gamma*M_inf*M_inf);
    
    P_R   = P_L*(2.0*gamma*M_inf*M_inf-(gamma-1.0))/(gamma+1.0);
    Rho_R = ((gamma+1.0)/(gamma-1.0)*(P_R/P_L)+1.0)/((gamma+1.0)/(gamma-1.0)+(P_R/P_L));
    u_R   =sqrt(gamma*(((2.0+(gamma-1.0)*M_inf*M_inf)*P_R)/((2.0*gamma*M_inf*M_inf+(1.0-gamma))*Rho_R)));
    
    
    // x_0 = 0.8, t = 0.012 units
    
    Terminating_Time = 10.0;
    for(int Cell_Index =0;Cell_Index<0.5*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
	Exact_Density[Cell_Index] = Rho_L;
	Exact_Pressure[Cell_Index] = P_L;
	Exact_Velocity[Cell_Index] = u_L;
	Exact_e[Cell_Index]= P_L/(Rho_L*(gamma-1.0));
    }
    
    for(int Cell_Index =0.5*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
	Exact_Density[Cell_Index] = Rho_R;
	Exact_Pressure[Cell_Index] = P_R;
	Exact_Velocity[Cell_Index] = u_R;
	Exact_e[Cell_Index]= P_R/(Rho_R*(gamma-1.0));

    }

    Directory_Name();
    File_Name();

    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
  
}


void Steady_Contact()
{
    Rho_L = 1.4;
    u_L   = 0.0;
    P_L   = 0.4;
    
    Rho_R = 1.0;
    u_R   = 0.0;
    P_R   = 0.4;
    
    // x_0 = 0.8, t = 0.012 units
    
    Terminating_Time = 5.0;
    g1.cmd("set xrange[0:1]");
    for(int Cell_Index =0;Cell_Index<0.5*(Total_Cells-1);Cell_Index++)
    {
        Calculate_Conserved_Variables(P_L,u_L,Rho_L);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
	Exact_Density[Cell_Index] = Rho_L;
	Exact_Pressure[Cell_Index] = P_L;
	Exact_Velocity[Cell_Index] = u_L;
	Exact_e[Cell_Index]= P_L/(Rho_L*(gamma-1.0));

    }
    
    for(int Cell_Index =0.5*(Total_Cells);Cell_Index<Total_Cells;Cell_Index++)
    {
        Calculate_Conserved_Variables(P_R,u_R, Rho_R);
        Calculate_Primitive_Variables(U_Cell,Cell_Index);
	Exact_Density[Cell_Index] = Rho_R;
	Exact_Pressure[Cell_Index] = P_R;
	Exact_Velocity[Cell_Index] = u_R;
	Exact_e[Cell_Index]= P_R/(Rho_R*(gamma-1.0));

    }

    Directory_Name();
    File_Name();


    Plot_Solution();
    Write_Solution(Initial_Solution_Op_File);
}

