#include "1D_Euler_FVM.h"

void Error_Function(double & x, double & f6)
{
  double f1;//ARG
	double f3;//E
	double f4;//VB
	double f5;//T
// 	double f6;//ERR
	
	f1 = x*x;
	f3 = 0.0;

	if (f1 < 20.0)
	{
		f3 = exp(-f1);
	}

	f4 = fabs(x);
	f5 = 1.0/(1.0+0.3275911*f4);
	f6 = 1.061405429*f5;
	f6 = (f6-1.453152027)*f5;
	f6 = (f6+1.421413741)*f5;
	f6 = (f6-0.284496736)*f5;
	f6 = (f6+0.254829592)*f5;
	f6 = 1 - f6*f3;

	if (x < 0.0)
	{
		f6 = -f6;
	}

	
}

void KFVS(int & Cell_No)
{
  double Beta =0.0, Term1 = 0.0,A_Plus,A_Minus,B,s=0.0,alpha=0.0;
 vector<double> G_1_Plus(3,0.0),G_1_Minus(3,0.0),G_2_Plus(3,0.0),G_2_Minus(3,0.0),G_3_Plus(3,0.0),G_3_Minus(3,0.0);
  
  u_L = 0.0; P_L = 0.0;Rho_L = 0.0; E_L=0.0;T_L =0.0;
 u_R = 0.0; P_R = 0.0;Rho_R = 0.0; E_R=0.0;T_R =0.0;
 
//  cout<<"In PVU_Dissipation Function\n";  
  
// Upwind flux at i-1/2 interface between i and i-1   
    
// values at i -1  Left state values
    u_L = Velocity[Cell_No - 1];
    P_L = Pressure[Cell_No - 1];
    Rho_L = Density[Cell_No - 1];
    E_L =  U[Cell_No - 1][2];
    c_L = C[Cell_No - 1];

    Beta = 0.5*gamma/(c_L*c_L);     	// beta = 1.0/(2*R*T) = 0.5*gamma / a*a; a*a = gamma*R*T
    s = u_L*sqrt(Beta);			//s = u*sqrt(Beta)
//     Error_Function(s,alpha);
    A_Plus = 0.5*(1.0 + erf(s));
    A_Minus = 0.5*(1.0 - erf(s));
    B = exp(-s*s);
    Term1 = 0.5*(B/sqrt(M_PI*Beta)); 	// Term1 = B/(2*sqrt(pi*beta))
    
    
    G_1_Plus[0] = (Rho_L*u_L)*A_Plus  		+	Rho_L*Term1;
    G_1_Plus[1] = (P_L + Rho_L*u_L*u_L)*A_Plus	+ 	Rho_L*u_L*Term1;
    G_1_Plus[2] = ((P_L + E_L)*u_L)*A_Plus 	+ 	(0.5*P_L + E_L)*Term1;
    
    G_1_Minus[0] = (Rho_L*u_L)*A_Minus  		-	Rho_L*Term1;
    G_1_Minus[1] = (P_L + Rho_L*u_L*u_L)*A_Minus	- 	Rho_L*u_L*Term1;
    G_1_Minus[2] = ((P_L + E_L)*u_L)*A_Minus 		- 	(0.5*P_L + E_L)*Term1;
    
    
// values at i Left state values
    
    u = Velocity[Cell_No];
    P = Pressure[Cell_No];
    Rho = Density[Cell_No];
    E =  U[Cell_No][2];
    c = C[Cell_No];
    
    Beta = 0.5*gamma/(c*c);
    s = u*sqrt(Beta);
//     Error_Function(s,erf(s));
    A_Plus = 0.5*(1.0 + erf(s));
    A_Minus = 0.5*(1.0 - erf(s));
    B = exp(-s*s);
    Term1 = 0.5*(B/sqrt(M_PI*Beta)); 
    
    
    G_2_Plus[0] = (Rho*u)*A_Plus  		+	Rho*Term1;
    G_2_Plus[1] = (P + Rho*u*u)*A_Plus		+ 	Rho*u*Term1;
    G_2_Plus[2] = ((P + E)*u)*A_Plus 		+ 	(0.5*P + E)*Term1;
    
    G_2_Minus[0] = (Rho*u)*A_Minus  		-	Rho*Term1;
    G_2_Minus[1] = (P + Rho*u*u)*A_Minus	- 	Rho*u*Term1;
    G_2_Minus[2] = ((P + E)*u)*A_Minus 		- 	(0.5*P + E)*Term1;


// values at i+1 Left state values
    u_R = Velocity[Cell_No + 1];
    P_R = Pressure[Cell_No + 1];
    Rho_R = Density[Cell_No + 1];
    E_R =  U[Cell_No + 1][2];
    c_R = C[Cell_No + 1];
    
    Beta = 0.5*gamma/(c_R*c_R);
    s = u_R*sqrt(Beta);
//     Error_Function(s,erf(s));
    A_Plus = 0.5*(1.0 + erf(s));
    A_Minus = 0.5*(1.0 - erf(s));
    B = exp(-s*s);
    Term1 = 0.5*(B/sqrt(M_PI*Beta)); 
    
    
    G_3_Plus[0] = (Rho_R*u_R)*A_Plus	  	+	Rho_R*Term1;
    G_3_Plus[1] = (P_R + Rho_R*u_R*u_R)*A_Plus	+ 	Rho_R*u_R*Term1;
    G_3_Plus[2] = ((P_R + E_R)*u_R)*A_Plus 	+ 	(0.5*P_R + E_R)*Term1;
    
    G_3_Minus[0] = (Rho_R*u_R)*A_Minus  	-	Rho_R*Term1;
    G_3_Minus[1] = (P_R + Rho*u_R*u_R)*A_Minus	- 	Rho_R*u_R*Term1;
    G_3_Minus[2] = ((P_R + E_R)*u_R)*A_Minus 	- 	(0.5*P_R + E_R)*Term1;
    
//    Flux at i+1/2 interface = G+|i + G-|i+1
    Right_Face_Average_Flux[0] = (G_2_Plus[0] + G_3_Minus[0]);
    Right_Face_Average_Flux[1] = (G_2_Plus[1] + G_3_Minus[1]);
    Right_Face_Average_Flux[2] = (G_2_Plus[2] + G_3_Minus[2]);
    
    Right_Face_Dissipative_Flux[0] = 0.0;
    Right_Face_Dissipative_Flux[1] = 0.0;
    Right_Face_Dissipative_Flux[2] = 0.0;

//    Flux at i -1/2 interface = G+|i-1 + G-|i
    Left_Face_Average_Flux[0] = (G_1_Plus[0] + G_2_Minus[0]);
    Left_Face_Average_Flux[1] = (G_1_Plus[1] + G_2_Minus[1]);
    Left_Face_Average_Flux[2] = (G_1_Plus[2] + G_2_Minus[2]);
    
    Left_Face_Dissipative_Flux[0] = 0.0;
    Left_Face_Dissipative_Flux[1] = 0.0;
    Left_Face_Dissipative_Flux[2] = 0.0;

}