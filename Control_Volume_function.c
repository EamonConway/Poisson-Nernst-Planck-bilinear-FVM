#include "mex.h"
#include "math.h"
#include "matrix.h"
void mexFunction (int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{ 
    
double *element = mxGetPr(prhs[0]);    
double *nodes = mxGetPr(prhs[1]);   
double *Control_Volume;
int num_nodes = mxGetNumberOfElements(prhs[1])/3;
int num_elements = mxGetNumberOfElements(prhs[0])/5;
int i,j,ii;
int node_numbers[4];
double Sub_Volume;
double shift_xi[4] = {-1, 0, 0, -1};
double shift_eta[4] = {-1, -1, 0, 0};
double x_coord[4], y_coord[4];
const double gauss_points[2] = {-sqrt(3)/3,sqrt(3)/3};

        
double xi_1, eta_1, N1_1, N2_1, N3_1, N4_1;
double xi_2, eta_2, N1_2, N2_2, N3_2, N4_2;
double xi_3, eta_3, N1_3, N2_3, N3_3, N4_3;
double xi_4, eta_4, N1_4, N2_4, N3_4, N4_4;
double r_c_1, r_c_2;
double r_c_3, r_c_4;


double diff_N1_1[2], diff_N2_1[2], diff_N3_1[2], diff_N4_1[2];
double diff_N1_2[2], diff_N2_2[2], diff_N3_2[2], diff_N4_2[2];
double diff_N1_3[2], diff_N2_3[2], diff_N3_3[2], diff_N4_3[2];
double diff_N1_4[2], diff_N2_4[2], diff_N3_4[2], diff_N4_4[2];

double diff_x_1[2], diff_y_1[2];
double diff_x_2[2], diff_y_2[2];
double diff_x_3[2], diff_y_3[2];
double diff_x_4[2], diff_y_4[2];
        
double Volume_map1, Volume_map2, Volume_map3, Volume_map4;
double F_1, F_2, F_3, F_4;

plhs[0] = mxCreateDoubleMatrix(num_nodes,1,mxREAL);
Control_Volume = mxGetPr(plhs[0]);


    for(i=0;i<num_elements;i++){
     
        for(ii=0;ii<4;ii++){            
            node_numbers[ii]    = element[num_elements*(ii+1) +i]-1; 
            y_coord[ii]         = nodes[node_numbers[ii] + 2*num_nodes];
            x_coord[ii]         = nodes[node_numbers[ii] + num_nodes]; 
        }


                for(j = 0;j<4;j++){ 
                    xi_1 = 0.5*(gauss_points[0] + 1) + shift_xi[j];
                    eta_1= 0.5*(gauss_points[0] + 1) + shift_eta[j];

                    xi_2 = 0.5*(gauss_points[1] + 1) + shift_xi[j];
                    eta_2= 0.5*(gauss_points[0] + 1) + shift_eta[j];

                    xi_3 = 0.5*(gauss_points[1] + 1) + shift_xi[j];
                    eta_3= 0.5*(gauss_points[1] + 1) + shift_eta[j];

                    xi_4 = 0.5*(gauss_points[0] + 1) + shift_xi[j];
                    eta_4= 0.5*(gauss_points[1] + 1) + shift_eta[j];


                    N1_1 = 0.25*(1-xi_1)*(1-eta_1);
                    N2_1 = 0.25*(1+xi_1)*(1-eta_1);
                    N3_1 = 0.25*(1+xi_1)*(1+eta_1);
                    N4_1 = 0.25*(1-xi_1)*(1+eta_1);

                    N1_2 = 0.25*(1-xi_2)*(1-eta_2);
                    N2_2 = 0.25*(1+xi_2)*(1-eta_2);
                    N3_2 = 0.25*(1+xi_2)*(1+eta_2);
                    N4_2 = 0.25*(1-xi_2)*(1+eta_2);
                    
                    N1_3 = 0.25*(1-xi_3)*(1-eta_3);
                    N2_3 = 0.25*(1+xi_3)*(1-eta_3);
                    N3_3 = 0.25*(1+xi_3)*(1+eta_3);
                    N4_3 = 0.25*(1-xi_3)*(1+eta_3);
                    
                    N1_4 = 0.25*(1-xi_4)*(1-eta_4);
                    N2_4 = 0.25*(1+xi_4)*(1-eta_4);
                    N3_4 = 0.25*(1+xi_4)*(1+eta_4);
                    N4_4 = 0.25*(1-xi_4)*(1+eta_4);
                    
                    diff_N1_1[0] = -0.25*(1-eta_1);
                    diff_N1_1[1] =  -0.25*(1-xi_1);
                    diff_N2_1[0] = 0.25*(1-eta_1);
                    diff_N2_1[1] =  -0.25*(1+xi_1);
                    diff_N3_1[0] = 0.25*(1+eta_1);
                    diff_N3_1[1] =  0.25*(1+xi_1);
                    diff_N4_1[0] = -0.25*(1+eta_1);
                    diff_N4_1[1] =  0.25*(1-xi_1);

                    diff_N1_2[0] = -0.25*(1-eta_2);
                    diff_N1_2[1] =  -0.25*(1-xi_2);
                    diff_N2_2[0] = 0.25*(1-eta_2);
                    diff_N2_2[1] =  -0.25*(1+xi_2);
                    diff_N3_2[0] = 0.25*(1+eta_2);
                    diff_N3_2[1] =  0.25*(1+xi_2);
                    diff_N4_2[0] = -0.25*(1+eta_2);
                    diff_N4_2[1] =  0.25*(1-xi_2);
                    
                    diff_N1_3[0] = -0.25*(1-eta_3);
                    diff_N1_3[1] =  -0.25*(1-xi_3);
                    diff_N2_3[0] = 0.25*(1-eta_3);
                    diff_N2_3[1] =  -0.25*(1+xi_3);
                    diff_N3_3[0] = 0.25*(1+eta_3);
                    diff_N3_3[1] =  0.25*(1+xi_3);
                    diff_N4_3[0] = -0.25*(1+eta_3);
                    diff_N4_3[1] =  0.25*(1-xi_3);

                    diff_N1_4[0] = -0.25*(1-eta_4);
                    diff_N1_4[1] =  -0.25*(1-xi_4);
                    diff_N2_4[0] = 0.25*(1-eta_4);
                    diff_N2_4[1] =  -0.25*(1+xi_4);
                    diff_N3_4[0] = 0.25*(1+eta_4);
                    diff_N3_4[1] =  0.25*(1+xi_4);
                    diff_N4_4[0] = -0.25*(1+eta_4);
                    diff_N4_4[1] =  0.25*(1-xi_4);
                    
                    diff_x_1[0] = x_coord[0]*diff_N1_1[0] + x_coord[1]*diff_N2_1[0] + x_coord[2]*diff_N3_1[0] + x_coord[3]*diff_N4_1[0]; 
                    diff_x_1[1] = x_coord[0]*diff_N1_1[1] + x_coord[1]*diff_N2_1[1] + x_coord[2]*diff_N3_1[1] + x_coord[3]*diff_N4_1[1];          
    
                    diff_y_1[0] = y_coord[0]*diff_N1_1[0] + y_coord[1]*diff_N2_1[0] + y_coord[2]*diff_N3_1[0] + y_coord[3]*diff_N4_1[0]; 
                    diff_y_1[1] = y_coord[0]*diff_N1_1[1] + y_coord[1]*diff_N2_1[1] + y_coord[2]*diff_N3_1[1] + y_coord[3]*diff_N4_1[1]; 
                                        
                    diff_x_2[0] = x_coord[0]*diff_N1_2[0] + x_coord[1]*diff_N2_2[0] + x_coord[2]*diff_N3_2[0] + x_coord[3]*diff_N4_2[0]; 
                    diff_x_2[1] = x_coord[0]*diff_N1_2[1] + x_coord[1]*diff_N2_2[1] + x_coord[2]*diff_N3_2[1] + x_coord[3]*diff_N4_2[1];          
                
                    diff_y_2[0] = y_coord[0]*diff_N1_2[0] + y_coord[1]*diff_N2_2[0] + y_coord[2]*diff_N3_2[0] + y_coord[3]*diff_N4_2[0]; 
                    diff_y_2[1] = y_coord[0]*diff_N1_2[1] + y_coord[1]*diff_N2_2[1] + y_coord[2]*diff_N3_2[1] + y_coord[3]*diff_N4_2[1]; 

                    diff_x_3[0] = x_coord[0]*diff_N1_3[0] + x_coord[1]*diff_N2_3[0] + x_coord[2]*diff_N3_3[0] + x_coord[3]*diff_N4_3[0]; 
                    diff_x_3[1] = x_coord[0]*diff_N1_3[1] + x_coord[1]*diff_N2_3[1] + x_coord[2]*diff_N3_3[1] + x_coord[3]*diff_N4_3[1];          
    
                    diff_y_3[0] = y_coord[0]*diff_N1_3[0] + y_coord[1]*diff_N2_3[0] + y_coord[2]*diff_N3_3[0] + y_coord[3]*diff_N4_3[0]; 
                    diff_y_3[1] = y_coord[0]*diff_N1_3[1] + y_coord[1]*diff_N2_3[1] + y_coord[2]*diff_N3_3[1] + y_coord[3]*diff_N4_3[1]; 
                                        
                    diff_x_4[0] = x_coord[0]*diff_N1_4[0] + x_coord[1]*diff_N2_4[0] + x_coord[2]*diff_N3_4[0] + x_coord[3]*diff_N4_4[0]; 
                    diff_x_4[1] = x_coord[0]*diff_N1_4[1] + x_coord[1]*diff_N2_4[1] + x_coord[2]*diff_N3_4[1] + x_coord[3]*diff_N4_4[1];          
                
                    diff_y_4[0] = y_coord[0]*diff_N1_4[0] + y_coord[1]*diff_N2_4[0] + y_coord[2]*diff_N3_4[0] + y_coord[3]*diff_N4_4[0]; 
                    diff_y_4[1] = y_coord[0]*diff_N1_4[1] + y_coord[1]*diff_N2_4[1] + y_coord[2]*diff_N3_4[1] + y_coord[3]*diff_N4_4[1]; 

                    
                    r_c_1					= y_coord[0]*N1_1 + y_coord[1]*N2_1 + y_coord[2]*N3_1 + y_coord[3]*N4_1; 
                    r_c_2					= y_coord[0]*N1_2 + y_coord[1]*N2_2 + y_coord[2]*N3_2 + y_coord[3]*N4_2; 
                    r_c_3					= y_coord[0]*N1_3 + y_coord[1]*N2_3 + y_coord[2]*N3_3 + y_coord[3]*N4_3; 
                    r_c_4					= y_coord[0]*N1_4 + y_coord[1]*N2_4 + y_coord[2]*N3_4 + y_coord[3]*N4_4; 

                    F_1 = r_c_1;        
                    F_2 = r_c_2;  
                    F_3 = r_c_3;  
                    F_4 = r_c_4;
                    
                    Volume_map1 = diff_x_1[0]*diff_y_1[1] - diff_x_1[1]*diff_y_1[0];
                    Volume_map2 = diff_x_2[0]*diff_y_2[1] - diff_x_2[1]*diff_y_2[0];
                    Volume_map3 = diff_x_3[0]*diff_y_3[1] - diff_x_3[1]*diff_y_3[0];
                    Volume_map4 = diff_x_4[0]*diff_y_4[1] - diff_x_4[1]*diff_y_4[0];
                    
                    Sub_Volume = 0.25*(Volume_map1*F_1 + Volume_map2*F_2 + Volume_map3*F_3 + Volume_map4*F_4);
                    Control_Volume[node_numbers[j]] = Control_Volume[node_numbers[j]] + Sub_Volume;
           
          
                }
    }


}
