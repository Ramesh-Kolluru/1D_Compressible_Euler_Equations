#include "1D_Euler_FVM.h"

vector<vector<double > >U,U_New,Eigen_Values;
vector<double> U_Cell,U_Cell_New,Left_Face_Average_Flux,Right_Face_Average_Flux,Pressure,Temperature,Density,Velocity,C,Internal_Energy;
vector<double> Left_Face_Dissipative_Flux,Right_Face_Dissipative_Flux,Cell_del_X,Cell_del_t,Net_Cell_Flux;
vector<double> G_P_L_Plus,G_P_L_Minus,G_P_R_Plus,G_P_R_Minus,G_C_L_Plus,G_C_L_Minus,G_C_R_Plus,G_C_R_Minus,Phi;
int Physical_Cells,Ghost_Cells,Total_Cells,Boundary_Cell_No,Interior_Cell_No,Test_Case;
double CFL,del_X,del_T,Eigen_Value_1,Eigen_Value_2,Eigen_Value_3,u_L,P_L,Rho_L,c_L,T_L,V_Mag_L,E_L,u_R,P_R,Rho_R,c_R,T_R,V_Mag_R;
double rho_error,rhou_error,rhoE_error,E_R,u,P,Rho,c,T,V_Mag,E,Length,Lambda_Max,Lambda_Min;
double R_I_1,R_I_2,R_I_3,Min_delta_t,Min_delta_X,Sigma,Terminating_Time;
int Boundary_Type,Implement_Type,Dissipation_Type,Flux_Scheme_Type;
Gnuplot g1,g2,g3,g4,g5;
vector<double> Exact_Density,Exact_Velocity,Exact_Pressure,Exact_e;

void Allocate_Memory()
{
//  Arrays to store values in each cell    
    U_Cell.resize(3,0.0);                                // Conserved variable vector
    U_Cell_New.resize(3,0.0);                           // Updated Conserved variable vector
    Left_Face_Average_Flux.resize(3,0.0);                   // Average Flux Vector on Left Face
    Right_Face_Average_Flux.resize(3,0.0);                  // Average Flux Vector on Right Face
    Left_Face_Dissipative_Flux.resize(3,0.0);                  // Dissipation Vector on Left Face
    Right_Face_Dissipative_Flux.resize(3,0.0);                 // Dissipation Vector on Right Face
    Net_Cell_Flux.resize(3,0.0);            // Net Flux from a cell 
    G_P_L_Plus.resize(3,0.0);
    G_P_L_Minus.resize(3,0.0);
    G_P_R_Plus.resize(3,0.0);
    G_P_R_Minus.resize(3,0.0);
    G_C_L_Plus.resize(3,0.0);
    G_C_L_Minus.resize(3,0.0);
    G_C_R_Plus.resize(3,0.0);
    G_C_R_Minus.resize(3,0.0);
    Phi.resize(3,0.0);
    
// Arrays for global values 

    for( int Cell_No = 0; Cell_No<Total_Cells;Cell_No++)
    {
        U.push_back(U_Cell);                        // Conserved Variable Vector for all Cells
        U_New.push_back(U_Cell);                     // Updated Consered Variable Vector for all Cells 
        Eigen_Values.push_back(U_Cell);             // Eigen values of each cell
        Pressure.push_back(0.0);                    // Pressure data for all the cells
        Temperature.push_back(0.0);                 // Temperature for all the cells
        Density.push_back(0.0);                     // Density for all the cells 
        Velocity.push_back(0.0);                    // Velocity for all the cells 
        C.push_back(0.0);                           // Acoustic velocity for all the cells
	Cell_del_X.push_back(Cell_No*del_X);                  // Width of the cell
	Exact_Density.push_back(0.0);
	Exact_e.push_back(0.0);
	Exact_Pressure.push_back(0.0);
	Exact_Velocity.push_back(0.0);
	Internal_Energy.push_back(0.0);
    }
    
    for( int Cell_No = 1; Cell_No<Total_Cells-1;Cell_No++)
    {
        Cell_del_t.push_back(0.0);                  // Timestep of the cell
    }
        
}


void Reset_Variables()
{
    u_L=0;P_L=0.0;Rho_L=0.0;c_L=0.0;T_L=0.0;V_Mag_L=0.0,E_L=0.0;
    u_R=0;P_R=0.0;Rho_R=0.0;c_R=0.0;T_R=0.0;V_Mag_R=0.0,E_R=0.0;
    u=0;P=0.0;Rho=0.0;c=0.0;T=0.0;V_Mag=0.0,E=0.0;

    Left_Face_Average_Flux[0] = 0.0;
    Left_Face_Average_Flux[1] = 0.0;
    Left_Face_Average_Flux[2] = 0.0;

// Average Average flux on Right Face     
    Right_Face_Average_Flux[0] = 0.0;
    Right_Face_Average_Flux[1] = 0.0;
    Right_Face_Average_Flux[2] = 0.0;

}



void Initalize()
{
    Total_Cells = Physical_Cells + Ghost_Cells;
    del_X = Length/Total_Cells;
    
    Allocate_Memory();
    
// Test case refers to the cases to be solved on Euler equations     
    cout<<"Test case being used \t"<<Test_Case<<endl;
    switch(Test_Case)
    {
        case 8: Sod_Test_Case_1();
            break;
        case 9: Sod_Test_Case_2();
            break;
        case 1:Toro_Test_Case_1();
            break;
	case 2:Toro_Test_Case_2();
            break;
        case 3:Toro_Test_Case_3();
            break;
        case 4:Toro_Test_Case_4();
            break;
        case 5:Toro_Test_Case_5();
            break;
	case 6: Steady_Shock();
	    break;
	case 7: Steady_Contact();
	    break;
        default: Toro_Test_Case_1();
            break;
    }
}




// void Initalize(string & IpFile)
// {
//  Initalize();
// }