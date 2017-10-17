#include "mex.h"
#include "math.h"


void mexFunction (int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[])
{ 
    
// Code is structured such that
//	- t (time)
//	- y (solution) C1 = 1:3:end, C2 = 2:3:end, Phi = 3:3:end
//	- Parameters 
//	- element
//	- nodes
//  - DL
double *y = mxGetPr(prhs[1]);
double *element = mxGetPr(prhs[3]);    
double *nodes = mxGetPr(prhs[4]);    
double *Dl_mat = mxGetPr(prhs[5]);

// Define Pointers for Parameters

double *D = mxGetPr(mxGetField(prhs[2],0,"D")); // Diffusion coefficients [D1,D2]
double *R = mxGetPr(mxGetField(prhs[2],0,"R")); // Gas Constant
double *T = mxGetPr(mxGetField(prhs[2],0,"T")); // Temperature **Must be scalar** Functionality could be incorporated
double *vphi = mxGetPr(mxGetField(prhs[2],0,"varphi")); // Dimensionless coefficent phi_0*F/(R*T)
double *mu = mxGetPr(mxGetField(prhs[2],0,"mu")); // Dimensionless coefficient T_0/L_0^2
double *permit = mxGetPr(mxGetField(prhs[2],0,"permit")); // Permitivity
double *z = mxGetPr(mxGetField(prhs[2],0,"z")); // Valence of ions [z1, z2]
double *Flux;
int num_nodes = mxGetNumberOfElements(prhs[4])/3; // nodes must be inputed as nx3 matrix. [node number , x_location, y_location]
int num_elements = mxGetNumberOfElements(prhs[3])/5; // elements inputed as mx5 matrix and counter clockwise [element_number, first_node, second_node, third_node, fourth_node]
int i,j,ii;
int node_numbers[4];


const double face_xi[4] = {0, 0.5, 0, -0.5};
const double face_eta[4] = {-0.5, 0, 0.5, 0};
const double unit_normal_x[4] = {1, 0, -1, 0};
const double unit_normal_y[4] = {0, 1, 0, -1};
double temp_unit_normal_x, temp_unit_normal_y;
double temp_unit_normal;
        
double dl, xi, eta, N1, N2, N3, N4;
double diff_N1[2], diff_N2[2], diff_N3[2], diff_N4[2];
double x_coord[4], y_coord[4];
double C_1[4], C_2[4], Phi[4];
double C_1face, C_2face, Phiface;
double Grad_C1[2], Grad_C2[2], Grad_Phi[2];
double diff_x[2], diff_y[2], r_c;
double det_J, jacobian_inverse[4];
double fluxC_1, fluxC_2, fluxPhi;
double scale;

// Output Information 
plhs[0] = mxCreateDoubleMatrix(3*num_nodes,1,mxREAL);
Flux = mxGetPr(plhs[0]);

    for(i=0;i<num_elements;i++){

        // Nodes in the element.
        for(ii=0;ii<4;ii++){            
            node_numbers[ii]    = element[num_elements*(ii+1) +i]-1;   // -1 for  referencing . 
            x_coord[ii]         = nodes[node_numbers[ii] + num_nodes]; 
            y_coord[ii]         = nodes[node_numbers[ii] + 2*num_nodes];
            C_1[ii]             = y[3*node_numbers[ii]];
            C_2[ii]             = y[3*node_numbers[ii]+1];
            Phi[ii]             = y[3*node_numbers[ii]+2];
     
        }
  
       
        
                // Calculate flux by looping over elements
                for(j = 0;j<4;j++){ 
                    // Length of face
                    dl = Dl_mat[(j+1)*num_elements + i];                    

					temp_unit_normal_x = unit_normal_x[j];
					temp_unit_normal_y = unit_normal_y[j];
                    // Natural Coordinates
                    xi = face_xi[j];
                    eta= face_eta[j];

                    // Shape Functions                    
                    N1 = 0.25*(1-xi)*(1-eta);
                    N2 = 0.25*(1+xi)*(1-eta);
                    N3 = 0.25*(1+xi)*(1+eta);
                    N4 = 0.25*(1-xi)*(1+eta);
                    
                    // Derivatives of Shape Functions     
                    diff_N1[0] = -0.25*(1-eta);
                    diff_N1[1] =  -0.25*(1-xi);
                    diff_N2[0] = 0.25*(1-eta);
                    diff_N2[1] =  -0.25*(1+xi);
                    diff_N3[0] = 0.25*(1+eta);
                    diff_N3[1] =  0.25*(1+xi);
                    diff_N4[0] = -0.25*(1+eta);
                    diff_N4[1] =  0.25*(1-xi);
                    
                    
                    // Value on faces and gradients;
                    C_1face = C_1[0]*N1 + C_1[1]*N2 + C_1[2]*N3 + C_1[3]*N4;
                    C_2face = C_2[0]*N1 + C_2[1]*N2 + C_2[2]*N3 + C_2[3]*N4;



                    diff_x[0] = x_coord[0]*diff_N1[0] + x_coord[1]*diff_N2[0] + x_coord[2]*diff_N3[0] + x_coord[3]*diff_N4[0]; 
                    diff_x[1] = x_coord[0]*diff_N1[1] + x_coord[1]*diff_N2[1] + x_coord[2]*diff_N3[1] + x_coord[3]*diff_N4[1];          
                
                    diff_y[0] = y_coord[0]*diff_N1[0] + y_coord[1]*diff_N2[0] + y_coord[2]*diff_N3[0] + y_coord[3]*diff_N4[0]; 
                    diff_y[1] = y_coord[0]*diff_N1[1] + y_coord[1]*diff_N2[1] + y_coord[2]*diff_N3[1] + y_coord[3]*diff_N4[1]; 

                    // Evaluate Unit Normals 
                    temp_unit_normal = temp_unit_normal_x;     
                    temp_unit_normal_x = diff_x[0]*temp_unit_normal_x + diff_x[1]*temp_unit_normal_y;
                    temp_unit_normal_y = diff_y[0]*temp_unit_normal + diff_y[1]*temp_unit_normal_y;

                    
                    scale = sqrt(pow(temp_unit_normal_x,2) + pow(temp_unit_normal_y,2));

                    temp_unit_normal_x = (1/scale)*temp_unit_normal_x;
                    temp_unit_normal_y = (1/scale)*temp_unit_normal_y;

                    Grad_C1[0] = C_1[0]*diff_N1[0] + C_1[1]*diff_N2[0] + C_1[2]*diff_N3[0] + C_1[3]*diff_N4[0]; 
                    Grad_C1[1] = C_1[0]*diff_N1[1] + C_1[1]*diff_N2[1] + C_1[2]*diff_N3[1] + C_1[3]*diff_N4[1]; 
                    
                    Grad_C2[0] = C_2[0]*diff_N1[0] + C_2[1]*diff_N2[0] + C_2[2]*diff_N3[0] + C_2[3]*diff_N4[0]; 
                    Grad_C2[1] = C_2[0]*diff_N1[1] + C_2[1]*diff_N2[1] + C_2[2]*diff_N3[1] + C_2[3]*diff_N4[1]; 
                    
                    Grad_Phi[0] = Phi[0]*diff_N1[0] + Phi[1]*diff_N2[0] + Phi[2]*diff_N3[0] + Phi[3]*diff_N4[0]; 
                    Grad_Phi[1] = Phi[0]*diff_N1[1] + Phi[1]*diff_N2[1] + Phi[2]*diff_N3[1] + Phi[3]*diff_N4[1]; 
                    
                    r_c         = y_coord[0]*N1 + y_coord[1]*N2 + y_coord[2]*N3 + y_coord[3]*N4; 

                    // Evaluate jacobian
                    det_J               = diff_x[0]*diff_y[1] - diff_x[1]*diff_y[0];
                    jacobian_inverse[0] = diff_y[1]/det_J;
                    jacobian_inverse[1] = -diff_x[1]/det_J;
                    jacobian_inverse[2] = -diff_y[0]/det_J;
                    jacobian_inverse[3] = diff_x[0]/det_J;

                    // Evaluate Flux
                    fluxC_1     = dl*r_c*mu[0]*temp_unit_normal_x*(jacobian_inverse[0]*(-D[0]*Grad_C1[0] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[0]) + 
                                        jacobian_inverse[2]*(-D[0]*Grad_C1[1] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[1])) +
                                    dl*r_c*mu[0]*temp_unit_normal_y*(jacobian_inverse[1]*(-D[0]*Grad_C1[0] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[0]) + 
                                            jacobian_inverse[3]*(-D[0]*Grad_C1[1] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[1]));

                    fluxC_2     = dl*r_c*mu[0]*temp_unit_normal_x*(jacobian_inverse[0]*(-D[1]*Grad_C2[0] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[0]) + 
                                        jacobian_inverse[2]*(-D[1]*Grad_C2[1] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[1])) +
                                    dl*r_c*mu[0]*temp_unit_normal_y*(jacobian_inverse[1]*(-D[1]*Grad_C2[0] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[0]) + 
                                            jacobian_inverse[3]*(-D[1]*Grad_C2[1] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[1]));
//                                  //           mexPrintf(" fluxC_2 = %.15e\n",  fluxC_2);
//                       mexPrintf(" fluxC_1 = %.15e\n",  fluxC_2);
                    fluxPhi     = dl*permit[0]*r_c*temp_unit_normal_x*(jacobian_inverse[0]*(Grad_Phi[0]) + jacobian_inverse[2]*(Grad_Phi[1])) +
                                    dl*permit[0]*r_c*temp_unit_normal_y*(jacobian_inverse[1]*(Grad_Phi[0]) + jacobian_inverse[3]*(Grad_Phi[1]));
                   
                    
//                       fluxC_1     = r_c*mu[0]*unit_normal_x[j]*(jacobian_inverse[0]*(-D[0]*Grad_C1[0] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[0]) + 
//                                         jacobian_inverse[2]*(-D[0]*Grad_C1[1] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[1])) +
//                                     r_c*mu[0]*unit_normal_y[j]*(jacobian_inverse[1]*(-D[0]*Grad_C1[0] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[0]) + 
//                                             jacobian_inverse[3]*(-D[0]*Grad_C1[1] - z[0]*vphi[0]*D[0]*C_1face*Grad_Phi[1]));
// //                        mexPrintf(" fluxC_1 = %e\n",  fluxC_1);
//                     fluxC_2     = r_c*mu[0]*unit_normal_x[j]*(jacobian_inverse[0]*(-D[1]*Grad_C2[0] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[0]) + 
//                                         jacobian_inverse[2]*(-D[1]*Grad_C2[1] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[1])) +
//                                   r_c*mu[0]*unit_normal_y[j]*(jacobian_inverse[1]*(-D[1]*Grad_C2[0] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[0]) + 
//                                             jacobian_inverse[3]*(-D[1]*Grad_C2[1] - z[1]*vphi[0]*D[1]*C_2face*Grad_Phi[1]));
// //                                            mexPrintf(" fluxC_2 = %.15f\n",  fluxC_2);
// //                        mexPrintf(" fluxC_2 = %e\n",  fluxC_2);
//                     fluxPhi     = permit[0]*r_c*unit_normal_x[j]*(jacobian_inverse[0]*Grad_Phi[0] + jacobian_inverse[2]*Grad_Phi[1]) +
//                                   permit[0]*r_c*unit_normal_y[j]*(jacobian_inverse[1]*Grad_Phi[0] + jacobian_inverse[3]*Grad_Phi[1]);
//                    mexPrintf(" fluxPhi = %e\n",  fluxPhi);
//                              mexPrintf(" fluxPhi1 = %.15f\n",  jacobian_inverse[0]*Grad_Phi[0] + jacobian_inverse[2]*Grad_Phi[1]);
//                     mexPrintf(" fluxPhi2 = %.15f\n",  jacobian_inverse[1]*Grad_Phi[0] + jacobian_inverse[3]*Grad_Phi[1]);
                    
//                     mexPrintf("f %.15f\n",unit_normal_x[j]);
//                     mexPrintf("times %.15f\n",permit[0]*(r_c*unit_normal_x[j]*(jacobian_inverse[0]*Grad_Phi[0] + jacobian_inverse[2]*Grad_Phi[1]) +
//                             r_c*unit_normal_y[j]*(jacobian_inverse[1]*Grad_Phi[0] + jacobian_inverse[3]*Grad_Phi[1])));
                    // Flux Concentration 1
                    Flux[3*node_numbers[j]] =  Flux[3*node_numbers[j]] + fluxC_1;
                    Flux[3*node_numbers[(j+1)%4]] =  Flux[3*node_numbers[(j+1)%4]] - fluxC_1;
               
                    // Flux Concentration 2 
                    Flux[3*node_numbers[j]+1] =  Flux[3*node_numbers[j]+1] + fluxC_2;
                    Flux[3*node_numbers[(j+1)%4]+1] =  Flux[3*node_numbers[(j+1)%4]+1] - fluxC_2;
                    
                    // Flux Potential 
                    Flux[3*node_numbers[j]+2] =  Flux[3*node_numbers[j]+2] + fluxPhi;
                    Flux[3*node_numbers[(j+1)%4]+2] =  Flux[3*node_numbers[(j+1)%4]+2] - fluxPhi;
//                     
//                      
               
                }
//                     
    }
// End Mex file


}