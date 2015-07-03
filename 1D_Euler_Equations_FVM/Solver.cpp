#include "1D_Euler_FVM.h"


void Evaluate_Average_Flux(int & Cell_No)
{
// Reset_Variables();
// Left state Variables------------------ i-1
    u_L=Velocity[Cell_No - 1];
    P_L=Pressure[Cell_No - 1];
    Rho_L=Density[Cell_No - 1];
    V_Mag_L = u_L*u_L;
    E_L = U[Cell_No-1][2];

// Right State Variables---------------------- i+1
    u_R=Velocity[Cell_No + 1];
    P_R=Pressure[Cell_No + 1];
    Rho_R=Density[Cell_No + 1];
    V_Mag_R = u_R*u_R;
    E_R = U[Cell_No + 1][2];

// Current Cell Variables-------------------------- i
    u=Velocity[Cell_No];
    P=Pressure[Cell_No];
    Rho=Density[Cell_No];
    V_Mag = u*u;
    E = U[Cell_No][2];
    
// Average Convective flux on left Face  between i and i -1 
    Left_Face_Average_Flux[0] = 0.5*(Rho*u + Rho_L*u_L);
    Left_Face_Average_Flux[1] = 0.5*((Rho*V_Mag + P)  + (Rho_L*V_Mag_L + P_L));
    Left_Face_Average_Flux[2] = 0.5*((E + P)*u  + (E_L + P_L)*u_L);

// Average Convective flux on Right Face between i and i+1
    Right_Face_Average_Flux[0] = 0.5*(Rho*u + Rho_R*u_R);
    Right_Face_Average_Flux[1] = 0.5*((Rho*V_Mag + P)  + (Rho_R*V_Mag_R + P_R));
    Right_Face_Average_Flux[2] = 0.5*((E + P)*u + (E_R + P_R)*u_R);
}



void Evaluate_Dissipation_Flux(int & Cell_No, int & Dissipation_Type)
{
    switch(Dissipation_Type)
    {
        case 0: LLF_Dissipation(Cell_No);
		break;
        case 1: PVU_Dissipation(Cell_No);
		break;
        case 2: PVU_P_Dissipation(Cell_No);
		break;
        case 3: MOVERS_Dissipation(Cell_No);
		break;
	case 9:	ROE(Cell_No);// Roe scheme for evaluating the dissipation	  
		break;	  
	case 12:Rusanov_Dissipation(Cell_No);
		break;
        default: LLF_Dissipation(Cell_No);
            break;
    }
}

//Dissipation_Type is a global Variable
void Evaluate_Net_Cell_Flux(int & Cell_No)
{
  switch(Flux_Scheme_Type)
  {
    case 0:
          // Calculates the Average Flux from Right and left faces 
	  Evaluate_Average_Flux(Cell_No); 
        // Calculate the Dissipation flux based on the type
	Evaluate_Dissipation_Flux(Cell_No,Dissipation_Type);
      break;
    case 1:      KFVS(Cell_No);
      break;
    case 2:      Steger_Warming(Cell_No);
      break;
    case 3:      Vanleer(Cell_No);
      break;
    case 4:      AUSM(Cell_No);
      break;
    case 5:      HLL(Cell_No);
      break;
    case 6:      HLLC(Cell_No);
      break;
    default: Evaluate_Average_Flux(Cell_No);           // Calculates the Average Flux from Right and left faces 
	     Evaluate_Dissipation_Flux(Cell_No,Dissipation_Type);// Calculate the Dissipation flux based on the type
   break;
  }
    Net_Cell_Flux[0] = 0.0;
    Net_Cell_Flux[1] = 0.0;
    Net_Cell_Flux[2] = 0.0;
    
    // Net Mass Flux from a cell
    Net_Cell_Flux[0] = 	(Right_Face_Average_Flux[0] - Left_Face_Average_Flux[0])
			-(Right_Face_Dissipative_Flux[0] - Left_Face_Dissipative_Flux[0]);
    
    //Net Momentum Flux from a cell
    Net_Cell_Flux[1] = (Right_Face_Average_Flux[1] - Left_Face_Average_Flux[1])
		      -(Right_Face_Dissipative_Flux[1] - Left_Face_Dissipative_Flux[1]);
    
    // Net Energy Flux from a cell
    Net_Cell_Flux[2] = (Right_Face_Average_Flux[2] - Left_Face_Average_Flux[2])
			-(Right_Face_Dissipative_Flux[2] - Left_Face_Dissipative_Flux[2]);
}

void Evaluate_Sigma()
{
//     cout<<"CFL\t"<<CFL<<endl;
    for(int Cell_No =1; Cell_No<Total_Cells;Cell_No++)
    {
//         cout<<Cell_del_X[Cell_No]<<"\t"<<Eigen_Values[Cell_No][0]<<endl;
      
      switch(Is_Kinetic_Scheme)  
      {
	case true:
	  Cell_del_t[Cell_No-1] = CFL*del_X/(fabs(Velocity[Cell_No]) + 3.0*sqrt(R*Temperature[Cell_No]));
	  break;
	case false:
	  Cell_del_t[Cell_No-1] = CFL*del_X/Eigen_Values[Cell_No][0];
	  break;
      }
//          cout<<Cell_del_t[Cell_No-1]<<endl;
    }
// Finds the minimum number in the list 
    Min_delta_t = *min_element(Cell_del_t.begin(),Cell_del_t.end());
//       cout<<Cell_del_t.size()<<"\tMinimum Delta t\t"<<Min_delta_t<<endl;
//      cout<<Min_delta_X<<"\t"<<del_X<<endl;
    Sigma = Min_delta_t/del_X;
//     cout<<Min_delta_t<<"\t"<<Min_delta_X<<"\t"<<Sigma<<endl;
}


void Evaluate_Time_Evolution()
{
// This function returns sigma     
    Evaluate_Sigma();
//     cout<< "value of sigma\t"<<Sigma<<endl;
    for(int Cell_No =1; Cell_No<Total_Cells-1;Cell_No++)
    {
        Evaluate_Net_Cell_Flux(Cell_No);
        // Update of Mass 
        U_New[Cell_No][0] = U[Cell_No][0] - Sigma*Net_Cell_Flux[0];
        //Update of Momentum
        U_New[Cell_No][1] = U[Cell_No][1] - Sigma*Net_Cell_Flux[1];
        //Update of Energy
        U_New[Cell_No][2] = U[Cell_No][2] - Sigma*Net_Cell_Flux[2];
    }
}


// void Evaluate_Time_Evolution()
// {
// // This function returns sigma     
//     Evaluate_Sigma();
// //     cout<< "value of sigma\t"<<Sigma<<endl;
//     for(int Cell_No =1; Cell_No<Total_Cells-1;Cell_No++)
//     {
//         Evaluate_Net_Cell_Flux(Cell_No);
// //  	Sigma = Cell_del_t[Cell_No-1]/del_X;
//         // Update of Mass 
//         U_New[Cell_No][0] = 0.5*(U[Cell_No-1][0] + U[Cell_No+1][0]) - Sigma*Net_Cell_Flux[0];
//         //Update of Momentum
//         U_New[Cell_No][1] = 0.5*(U[Cell_No-1][1] + U[Cell_No+1][1]) - Sigma*Net_Cell_Flux[1];
//         //Update of Energy
//         U_New[Cell_No][2] = 0.5*(U[Cell_No-1][2] + U[Cell_No+1][2]) - Sigma*Net_Cell_Flux[2];
//     }
// }

void Find_Error_and_Update()
{
    rho_error = 0.0,rhou_error=0.0,rhoE_error=0.0;
    double temp1 =0.0,temp2=0.0,temp3=0.0;
    for(int Cell_No=1;Cell_No<Total_Cells-1;Cell_No++)
    {
 
// Error in Density        Interior_Cell_No
        temp1 = (U_New[Cell_No][0]-U[Cell_No][0]);///U[Cell_No][0];
        rho_error += temp1*temp1;
// Error in Momentum
        temp2 = (U_New[Cell_No][1]-U[Cell_No][1]);///U[Cell_No][1];
        rhou_error += temp2*temp2;
// Error in Energy
        temp3 = (U_New[Cell_No][2]-U[Cell_No][2]);///U[Cell_No][2];
        rhoE_error += temp3*temp3; 
        
// Updating the New values with old values for density, Momentum and Energy    
        U[Cell_No][0] = U_New[Cell_No][0];
        U[Cell_No][1] = U_New[Cell_No][1];
        U[Cell_No][2] = U_New[Cell_No][2];
	  if(isnan(U[Cell_No][0])|isnan(U[Cell_No][1])|isnan(U[Cell_No][2]))
	  {
	    cout<<"NaN occured exiting\n";
	    Write_Solution(Solution_Op_File);
	    
	    Plot_Solution();
	    exit(0);
	  }
    }

        rho_error = sqrt(rho_error);
        rhou_error = sqrt(rhou_error);
        rhoE_error = sqrt(rhoE_error);
}


void Solver()
{
    int iterations =0;
    double Total_Time =0.0;
    do
    {
        Boundary_Conditions(Boundary_Type,Implement_Type);
        Evaluate_Time_Evolution();
        Find_Error_and_Update();
        Calculate_Primitive_Variables();
        if(iterations%1 ==0)
        {
//             Write_Solution(Solution_Op_File);
            cout<<Terminating_Time<<"\t"<<Min_delta_t << "\t" <<Total_Time<<"\t"<<iterations<<"\t"<<rho_error<<"\t"<<rhou_error<<"\t"<<rhoE_error<<endl;
//             Plot_Solution();
        }
        iterations++;
        Total_Time+=Min_delta_t;
//         Plot_Solution();
//         cout<<iterations<<"\t"<<rho_error<<"\t"<<rhou_error<<"\t"<<rhoE_error<<"\t"<<Total_Time<<endl;
    }while(Total_Time<Terminating_Time);
    
}