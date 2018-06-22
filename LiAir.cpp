/*
Modeling and Sim of Semiconductors
Grant Brown
*/

#include "toolkit.hpp"

//--------------------Variables-------------------------------------start
double precision = 1e-6;    // epsilon

const int N = 30;           // discretized intervals
const int n = 5*(N+1);      // number of unknowns
double tau = 3600;             // time step                    (s)
double t = 0;

double L = 0.075;          // width of cathode region      (cm)
double Vt = 0.025;          // thermal voltage              (V)
double h = L/(N-1);         // width step                   (cm)

double F = 96485;           // Faraday constant             (C/mol)
double I = 0.1e-3;            // discharge current            (A)
double A = 1;               // cross-sectional area         (cm2)
double e0 = 0.75;           // initial porosity             (%)
double B = 0.5;             // factor for Butler-Volmer equation
double r0 = 20e-7;          // initial pore radius          (cm)

double rho_Li2O2 = 2.14;    // mass density of Li2O2        (g/cm3)
double rho_C = 2.26;        // mass density of carbon       (g/cm3)
double M_Li2O2 = 45.88;     // molecular mass of Li2O2      (g/mol)
double c0 = 3.26e-6;        // external O2 concentration    (mol/cm3)

double sigma = 1.2;         // electron conductivity        (S/cm)
double K = 0.0124;          // lithium conductivity         (S/cm)
double D_Li = 1e-5;         // lithium diffusion coefficient (cm2/s)
double t_plus = 0.25;        // transference number
double K_D = Vt*K*(1-t_plus);   // lithium diffusional conductivity
double D_O2 = 7e-6;         // oxygen diffusion coefficient

double k = 1.76e-10;        // reaction rate at cathode     (mol cm/s)
double I0 = 1.3e-8;         // reaction rate at anode       (mol/s/cm2)

//--------------------Variables - States & Jacobian----------------start
vector<double> xo(n,0), x(n), y(n,0), line(n+1,0);
vector<vector<double>> AugMatrix(n,line); // J | rhs

double *phi = &x[0];
double *phi_Li = &x[(N+1)];
double *c_Li = &x[2*(N+1)];
double *c_O2 = &x[3*(N+1)];
double *porosity = &x[4*(N+1)];

double *ophi = &xo[0];
double *ophi_Li = &xo[(N+1)];
double *oc_Li = &xo[2*(N+1)];
double *oc_O2 = &xo[3*(N+1)];
double *oporosity = &xo[4*(N+1)];

//--------------------Variables - Parameters------------------------end
//--------------------Variables-------------------------------------end

//--------------------Functions Declarations----------------------start
double specificCap(double t);
double minporosity(vector<double> x);
void ComputeFunct(vector < vector <double>> &rhs);
void ComputeJacob(vector < vector <double>> &J);

// Reaction Rate of Cathode        (mol/s/cm3)
double rrC(double cO2, double porosity, double phi_Li, double phi);
// Derivative of Reaction Rate at Cathode
double drrC_dphi(double cO2, double porosity, double phi_Li, double phi);
double drrC_dphiLi(double cO2, double porosity, double phi_Li, double phi);
double drrC_dcO2(double cO2, double porosity, double phi_Li, double phi);
double drrC_dporosity(double cO2, double porosity, double phi_Li, double phi);
// Reaction Rate of Anode          (mol/s/cm2)
double rrA(double phi_Li);
double drrA_dphiLi(double phi_Li);

//--------------------Function Declarations-------------------------end

int main()
{
    FILE *fp = fopen("voltage.csv", "w");//save the data in this file
    FILE *fp1 = fopen("VoltageSpecificCap.csv", "w");
    FILE *fp2 = fopen("VoltageHalfSpecificCap.csv", "w");
	cout << "SIMULATION STARTED...\n";
	int epoch = 0, NewtonCnt;
	double normY, specCap;
	bool nanFlag = false, specCapFlag = true;

        
    for(int i = 0; i <= N; i++)
    {
        phi[i] = 2.95;
        phi_Li[i] = 2.95;
        c_Li[i] = 1e-3;
        c_O2[i] = c0;
        porosity[i] = 0.75;
    }
    while(minporosity(x) > 0 && !nanFlag)
    {
        xo = x;
        NewtonCnt=0;
        cout << "Epoch: "<< ++epoch << endl;
        do
        {
            ComputeFunct(AugMatrix);
            ComputeJacob(AugMatrix);
            
            y = solveGauss(AugMatrix);
            
            for(int i = 0; i < n; i++)
                x[i] -= y [i];
            
            cout << "NewtonCnt: " << ++NewtonCnt << endl;
            normY = norm(y);
            cout << "Norm of Y: " << normY << endl;
            cout << "---------------------------------------------------------------------------" << endl;
            
            if (std::isnan(normY))
                nanFlag = true;
        } while(!nanFlag && normY > precision);

        epoch++;                                            // epoch counter
        t += tau;                                           // update time by tau (time step)
        printf("T = %f [s], %f [hr], %f [days]\n", t, t/3600, t/(24*3600));
        fprintf(fp, "%f, %f\n", t/3600, phi[N]);            // prints voltage with time in file "voltage.csv"
        specCap = specificCap(t/3600);
        fprintf(fp1, "%f, %f\n", specCap, phi[N]);          // prints voltage with specCap in file "voltageSpecificCap.csv"
        if(specCap <= 1.11975 && specCapFlag)
        {// 2.2395 / 2 = 1.11975
            fprintf(fp2, "%f, %f\n", specCap, phi[N]);
            exportdata(x, "x_halfSpecCap.csv");
        }
    }
    // saves the values of o old into x because if not they are nan values.
    x=xo;
    exportdata(x,"x_file.csv");
    
    // Print state variable when baterry dies
    cout << "---------------------------------------------------------------------------\n"
        << "Cell discharge completed after:\n"
        << "---------------------------------------------------------------------------\n";
    printf("T = %f [s], %f [hr], %f [days]\n", t, t/3600, t/(24*3600));
    int width = 15;
    cout << setprecision(7);
    cout << left << setw(10) << "Position"
        << setw(width+2) << "| Phi"
        << setw(width+2) << "| Phi Li"
        << setw(width+2) << "| C Li"
        << setw(width+2) << "| C O2"
        << setw(width+2) << "| Porosity" << endl;
    for(int i = 0; i <= N; i++)
    {
        cout << left << setw(10) << i << "| "
            << setw(width) << phi[i] << "| "
            << setw(width) << phi_Li[i] << "| "
            << setw(width) << c_Li[i] << "| "
            << setw(width) << c_O2[i] << "| "
            << setw(width) << porosity[i] << endl;
    }
    
    cout << "...SIMULATION FINSIHED!\n";
    return 0;
}
double specificCap(double t)
{
    return t*I/(rho_C*(1-e0)*A*L);
}
double minporosity(vector<double> x)
{
    double r = 1;
    // minporisity less than 1
    for (int i = 0; i <= N; i++)
        if (porosity[i] < r)
            r = porosity[i];
            
    return r;
    // return min_element(porosity, porosity[N]);
}
void ComputeFunct(vector < vector <double>> &rhs)
{
    // ------------------Equation 1-----------------------------------------start
    // electrostatic potential of electrons ---- phi                (V)
    rhs[0][n] = (phi[1] - phi[0])/h;
    
    for(int i = 1; i < N; i++)
        rhs[i][n] = sigma * (phi[i+1]+phi[i-1]-2*phi[i])/pow(h,2)
            + F*rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        
    rhs[N][n] = -sigma*A*(phi[N]-phi[N-1])/h - I;
    // ------------------Equation 1-------------------------------------------end               check
    
    // ------------------Equation 2-----------------------------------------start
    // electrostatic potential of Li ---- phi li                      (V)
    rhs[(N+1)][n] = -K * (phi_Li[1]-phi_Li[0])/h + K_D * (log(c_Li[1])-log(c_Li[0]))/h
        - F * rrA(phi_Li[0]);
        
    for(int i = 1; i < N; i++)
        rhs[(N+1)+i][n] = K * (phi_Li[i+1]+phi_Li[i-1]-2*phi_Li[i])/pow(h,2)
            - K_D * (log(c_Li[i+1])+log(c_Li[i-1])-2*log(c_Li[i]))/pow(h,2)
            - F * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    
    rhs[(N+1)+N][n] = -K * (phi_Li[N] - phi_Li[N-1])/h
        + K_D * (log(c_Li[N])-log(c_Li[N-1]))/h;
    // ------------------Equation 2-------------------------------------------end               check
    
    // ------------------Equation 3-----------------------------------------start
    // concentration of lithium ions ---- c li                        (mol/cm3)
    rhs[2*(N+1)][n] = c_Li[0] - 1e-3;
    
    for(int i = 1; i < N; i++)
        rhs[2*(N+1)+i][n] = (porosity[i]*c_Li[i]-oporosity[i]*oc_Li[i])/tau
            - D_Li * (c_Li[i+1]+c_Li[i-1]-2*c_Li[i])/pow(h,2)
            + (1-t_plus) * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
            
    rhs[2*(N+1)+N][n] = c_Li[N]-c_Li[N-1];
    // ------------------Equation 3-------------------------------------------end               check
    
    // ------------------Equation 4-----------------------------------------start
    // concentration of oxygen ions ---- c o2                        (mol/cm3)
    rhs[3*(N+1)][n] = c_O2[1] - c_O2[0];
    
    for(int i = 1; i < N; i++)
        rhs[3*(N+1)+i][n] = (porosity[i]*c_O2[i]-oporosity[i]*oc_O2[i])/tau
            - D_O2 *(c_O2[i+1]+c_O2[i-1]-2*c_O2[i])/pow(h,2)
            + 0.5 * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
            
    rhs[3*(N+1)+N][n] = c_O2[N] - c0;
    // ------------------Equation 4-------------------------------------------end               check
    
    // ------------------Equation 5-----------------------------------------start
    // porosity ------------------------------------------------------- (%)
    for(int i = 0; i <= N; i++)
        rhs[4*(N+1)+i][n] = (porosity[i]-oporosity[i])/tau
            + M_Li2O2/(2*rho_Li2O2) * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    // ------------------Equation 5-------------------------------------------end
    
}
void ComputeJacob(vector < vector <double>> &J)
{
    // ------------------Equation 1-----------------------------------------start
    /*  electrostatic potential of electrons ---- phi                (V)
    rhs[0][n] = (phi[1] - phi[0])/h;
    
    for(int i = 1; i < N; i++)
        rhs[i][n] = sigma * (phi[i+1]+phi[i-1]-2*phi[i])/pow(h,2)
            + F*rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        
    rhs[N][n] = -sigma*A*(phi[N]-phi[N-1])/h - I;
    */
    //double * phi_J = &J;
    J[0][1] = 1/h;
    J[0][0] = -1/h;
    
    for(int i = 1; i < N; i++)
    {
        J[i][i+1] = sigma/pow(h,2);
        J[i][i-1] = sigma/pow(h,2);
        J[i][i] = -2*sigma/pow(h,2) + F * drrC_dphi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[i][(N+1)+i] = F * drrC_dphiLi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[i][3*(N+1)+i] = F * drrC_dcO2(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[i][4*(N+1)+i] = F * drrC_dporosity(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    }
    
    J[N][N-1] = sigma*A/h;
    J[N][N] = -sigma*A/h;
    // ------------------Equation 1-------------------------------------------end               check
    
    
    // ------------------Equation 2-----------------------------------------start
    /*  electrostatic potential of Li ---- phi li                      (V)
    rhs[(N+1)][n] = -K * (phi_Li[1]-phi_Li[0])/h + K_D * (log(c_Li[1])-log(c_Li[0]))/h
        - F * rrA(phi_Li[0]);
        
    for(int i = 1; i < N; i++)
        rhs[(N+1)+i][n] = K * (phi_Li[i+1]+phi_Li[i-1]-2*phi_Li[i])/pow(h,2)
            - K_D * (log(c_Li[i+1])+log(c_Li[i-1])-2*log(c_Li[i]))/pow(h,2)
            - F * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    
    rhs[(N+1)+N][n] = -K * (phi_Li[N] - phi_Li[N-1])/h
        + K_D * (log(c_Li[N])-log(c_Li[N-1]))/h;
    */
    J[(N+1)][(N+1)] = K/h - F * drrA_dphiLi(phi_Li[0]);
    J[(N+1)][(N+1)+1] = -K/h;
    J[(N+1)][2*(N+1)] = -K_D/(h*c_Li[0]);
    J[(N+1)][2*(N+1)+1] = K_D/h*(1/c_Li[1]);
    
    for(int i = 1; i < N; i++)
    {
        J[(N+1)+i][i] = -F * drrC_dphi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[(N+1)+i][(N+1)+i+1] = K/pow(h,2);
        J[(N+1)+i][(N+1)+i-1] = K/pow(h,2);
        J[(N+1)+i][(N+1)+i] = -2*K/pow(h,2) - F * drrC_dphiLi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[(N+1)+i][2*(N+1)+i-1] = -K_D/pow(h,2) * (1/c_Li[i-1]);
        J[(N+1)+i][2*(N+1)+i+1] = -K_D/pow(h,2) * (1/c_Li[i+1]);
        J[(N+1)+i][2*(N+1)+i] = 2*K_D/pow(h,2) * (1/c_Li[i]);
        J[(N+1)+i][3*(N+1)+i] = -F * drrC_dcO2(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[(N+1)+i][4*(N+1)+i] = -F * drrC_dporosity(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    }
    
    J[(N+1)+N][(N+1)+N] = -K/h;
    J[(N+1)+N][(N+1)+N-1] = K/h;
    J[(N+1)+N][2*(N+1)+N] = K_D/h * (1/c_Li[N]);
    J[(N+1)+N][2*(N+1)+N-1] = -K_D/h * (1/c_Li[N-1]);
    // ------------------Equation 2-------------------------------------------end               check
    
    
    // ------------------Equation 3-----------------------------------------start
    /*  concentration of lithium ions ---- c li                        (mol/cm3)
    rhs[2*(N+1)][n] = c_Li[0] - 1e-3;
    
    for(int i = 1; i < N; i++)
        rhs[2*(N+1)+i][n] = (porosity[i]*c_Li[i]-oporosity[i]*oc_Li[i])/tau
            - D_Li * (c_Li[i+1]+c_Li[i-1]-2*c_Li[i])/pow(h,2)
            + (1-t_plus) * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
            
    rhs[2*(N+1)+N][n] = c_Li[N]-c_Li[N-1];
    */
    J[2*(N+1)][2*(N+1)] = 1;
    
    for(int i = 1; i < N; i++)
    {
        J[2*(N+1)+i][i] = (1-t_plus) * drrC_dphi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[2*(N+1)+i][(N+1)+i] = (1-t_plus) * drrC_dphiLi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[2*(N+1)+i][2*(N+1)+i-1] = -D_Li/pow(h,2);
        J[2*(N+1)+i][2*(N+1)+i+1] = -D_Li/pow(h,2);
        J[2*(N+1)+i][2*(N+1)+i] = porosity[i]/tau+2*D_Li/pow(h,2);
        J[2*(N+1)+i][3*(N+1)+i] = (1-t_plus) * drrC_dcO2(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[2*(N+1)+i][4*(N+1)+i] = c_Li[i]/tau-(1-t_plus)*drrC_dporosity(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    }
    
    J[2*(N+1)+N][2*(N+1)+N] = 1/h;
    J[2*(N+1)+N][2*(N+1)+N-1] = -1/h;
    // ------------------Equation 3-------------------------------------------end               check
    
    
    // ------------------Equation 4-----------------------------------------start
    /*  concentration of oxygen ions ---- c o2                        (mol/cm3)
    rhs[3*(N+1)][n] = c_O2[1] - c_O2[0];
    
    for(int i = 1; i < N; i++)
        rhs[3*(N+1)+i][n] = (porosity[i]*c_O2[i]-oporosity[i]*oc_O2[i])/tau
            - D_O2 *(c_O2[i+1]+c_O2[i-1]-2*c_O2[i])/pow(h,2)
            + 1/2 * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
            
    rhs[3*(N+1)+N][n] = c_O2[N] - c0;
    */
    J[3*(N+1)][3*(N+1)] = -1/h;
    J[3*(N+1)][3*(N+1)+1] = 1/h;
    
    for(int i = 1; i < N; i++)
    {
        J[3*(N+1)+i][i] = 0.5 * drrC_dphi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[3*(N+1)+i][(N+1)+i] = 0.5 * drrC_dphiLi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[3*(N+1)+i][3*(N+1)+i-1] = -D_O2/pow(h,2);
        J[3*(N+1)+i][3*(N+1)+i] = porosity[i]/tau + 2*D_O2/pow(h,2)+ 0.5*drrC_dcO2(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[3*(N+1)+i][3*(N+1)+i+1] = -D_O2/pow(h,2);
        J[3*(N+1)+i][4*(N+1)+i] = c_O2[i]/tau + 0.5 * drrC_dporosity(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    }
    
    J[3*(N+1)+N][3*(N+1)+N] = 1;
    // ------------------Equation 4-------------------------------------------end               check
    
    
    // ------------------Equation 5-----------------------------------------start
    /*  Porosity ---- E                                             (%)
    for(int i = 0; i < N; i++)
        rhs[4*(N+1)+i][n] = (porosity[i]-oporosity[i])/tau
            + M_Li2O2/(2*rho_Li2O2) * rrC(c_O2[i], porosity[i], phi_Li[i], phi[i]);
    */
    for(int i = 0; i <= N; i++)
    {
        J[4*(N+1)+i][i] = M_Li2O2/(2*rho_Li2O2) * drrC_dphi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[4*(N+1)+i][(N+1)+i] = M_Li2O2/(2*rho_Li2O2) * drrC_dphiLi(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[4*(N+1)+i][3*(N+1)+i] = M_Li2O2/(2*rho_Li2O2) * drrC_dcO2(c_O2[i], porosity[i], phi_Li[i], phi[i]);
        J[4*(N+1)+i][4*(N+1)+i] = 1/tau + M_Li2O2/(2*rho_Li2O2) * drrC_dporosity(c_O2[i], porosity[i], phi_Li[i], phi[i]);

    }
    // ------------------Equation 5-------------------------------------------end               check
    
}
double rrC(double cO2, double porosity, double phi_Li, double phi)//------------------------    check?
{
    return (k*2*cO2*sqrt(porosity*e0)/r0)*(exp(B*(phi_Li-phi)/Vt)-exp(-B*(phi_Li-phi)/Vt));
} // check
double drrC_dphi(double cO2, double porosity, double phi_Li, double phi)
{
    return (k*2*cO2*sqrt(porosity*e0)/r0)*(-B/Vt)*(exp(B*(phi_Li-phi)/Vt)+exp(-B*(phi_Li-phi)/Vt));
} // check
double drrC_dphiLi(double cO2, double porosity, double phi_Li, double phi)
{
    return (k*2*cO2*sqrt(porosity*e0)/r0)*(B/Vt)*(exp(B*(phi_Li-phi)/Vt)+exp(-B*(phi_Li-phi)/Vt));
} // check
double drrC_dcO2(double cO2, double porosity, double phi_Li, double phi)
{
    return k*(2*sqrt(porosity*e0)/r0)*(exp(B*(phi_Li-phi)/Vt)-exp(-B*(phi_Li-phi)/Vt));
} // check
double drrC_dporosity(double cO2, double porosity, double phi_Li, double phi)
{
    return (k*cO2*sqrt(e0)/(sqrt(porosity)*r0))*(exp(B*(phi_Li-phi)/Vt)-exp(-B*(phi_Li-phi)/Vt));
} // check
double rrA(double phi_Li)
{
    return I0*(exp(B*(-phi_Li+2.95)/Vt)-exp(-B*(-phi_Li+2.95)/Vt));
} // check
double drrA_dphiLi(double phi_Li)
{
    return I0*(-B/Vt)*(exp(B*(-phi_Li+2.95)/Vt) + exp(-B*(-phi_Li+2.95)/Vt));
} // check
