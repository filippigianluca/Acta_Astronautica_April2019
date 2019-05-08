//  Action_calculator_mex.c
//  
//
//  Created by Lorenzo Teofili on 14/10/15.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// Action_calculator_forces2_torques2_massive_point computes
// the DIMENTIONAL gravitational actions acting on a body due to the
// presence of an other body. The irregular shape of both is taken into
// consideration. Consider as a reference the model outlined by Eugene
// G. Fahnestock & Daniel J. Scheeres in the paper:
//
// "Simulation of the full two rigid body problem using polyhedral
// mutual potential and potential derivatives approach".
//
// The approximation is carried out up to the SECOND ORDER both for
// the torques and for the forces.
//
// This particular function considers one of the bodies as a massive
// point, simplifying the problem.
//
// G is expressed in Km^3, thus :  6.67384*10^(-20);
//
//
//
// [F, E_m] = action_calculator( R, r1a, r2a, r3a, Ja, rho_a, M_b, nfaces_a)
//
// INPUT:
// R            - A vector that indicates the relative position of the two asteroids.
//                It have to be expressed in one reference frame fixed with one of the
//                bodies, centered in its baricenter and rotating with the selected asteroid,
//                one can call it reference A or reference B.
//                HERE R MUST BE EXPRESSED IN THE REFERENCE FRAME OF THE MASSIVE POINT B, attitude
//                in that case is not important.
//
// r1a,r2a,r3a  - Matrices obtained by the function "load_model_2bp", they take into account for
//                the irregular shape of the bodies. They have to be expressed coerently with
//                R, for instance if R is given in the A reference r1a r2a and r3a must be provided in
//                the A reference (it means they remain the same). In the other hand, if R is
//                given in the B reference r1a r2a and r3a must be provided in
//                the B reference (that is they must be rotate!).
//
// Ja           - Vector obtained by the function "Asteroid_Inertia_matrix_2".
//                It represents the vector of the determinant of the polihedra
//                of the asteroids.
//
// rho_a        - Scalar, asteroid A density, the asteroid is suppose to have
//                constant density.
//
// M_b          - The mass of the massive point which represents the B body.
//
// nfaces_a     - Scalar obtained by the function "load_model_2bp", it needs for the for cycle.
//
// OUTPUT:
// F            - The force vector expressed in the selected body fixed reference frame.
//
// E_m          - Matrix to compute the torques acting on the the selected body fixed reference frame,
//                see the reference above for further details.
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


#include <math.h>
#include <matrix.h>
#include <mex.h>

//using namespace std;


void action_calculator(double *R,double *r1a,double *r2a,double *r3a,double *Ja,double rho_a,double M_b,int nfaces_a,double *Fa,double *E_m){

    //mexPrintf("after function calling\n");
    //mexPrintf("rho_a=%f\n",rho_a);
    //mexPrintf("M_b=%f\n",M_b);
    //mexPrintf("R=%f\n",R[0]);
    //mexPrintf("nfaces_a=%d\n",nfaces_a);
    //mexPrintf("nfaces_bf=%f\n",nfaces_b);
    //mexPrintf("r1a=%20.3f\n",r1a[0]);
    //mexPrintf("r1b1=%20.3f\n",r1b[0]);
    //mexPrintf("r2b1=%20.3f\n",r2b[0]);
    //mexPrintf("r3b1=%20.3f\n",r3b[0]);
    
    
    
    //do something
    
    const double Grav_con = 6.67384*pow(10,-20);
    
    const double Q0 = 1.0/36.0;
    const double Q1 = 1.0/144.0;
    
    const double Q2[36] = {
        8.0/2880.0, 4.0/2880.0, 4.0/2880.0, 5.0/2880.0, 5.0/2880.0, 5.0/2880.0,
        4.0/2880.0, 8.0/2880.0, 4.0/2880.0, 5.0/2880.0, 5.0/2880.0, 5.0/2880.0,
        4.0/2880.0, 4.0/2880.0, 8.0/2880.0, 5.0/2880.0, 5.0/2880.0, 5.0/2880.0,
        5.0/2880.0, 5.0/2880.0, 5.0/2880.0, 8.0/2880.0, 4.0/2880.0, 4.0/2880.0,
        5.0/2880.0, 5.0/2880.0, 5.0/2880.0, 4.0/2880.0, 8.0/2880.0, 4.0/2880.0,
        5.0/2880.0, 5.0/2880.0, 5.0/2880.0, 4.0/2880.0, 4.0/2880.0, 8.0/2880.0
    };
    
    //const double Q3[216]=
    //{
    //    12.0/8640.0,4.0/8640.0,4.0/8640.0,6.0/8640.0,6.0/8640.0,6.0/8640.0,
    //    4.0/8640.0,4.0/8640.0,2.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,
    //    4.0/8640.0,2.0/8640.0,4.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,
    //    6.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,
    //    6.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,
    //    6.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,
    //    
    //    4.0/8640.0,4.0/8640.0,2.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,
   //     4.0/8640.0,12.0/8640.0,4.0/8640.0,6.0/8640.0,6.0/8640.0,6.0/8640.0,
    //    2.0/8640.0,4.0/8640.0,4.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,
     //   3.0/8640.0,6.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,
      //  3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,
       // 3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,
        
  //      4.0/8640.0,2.0/8640.0,4.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,
   //     2.0/8640.0,4.0/8640.0,4.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,
    //    4.0/8640.0,4.0/8640.0,12.0/8640.0,6.0/8640.0,6.0/8640.0,6.0/8640.0,
     //   3.0/8640.0,3.0/8640.0,6.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,
      //  3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,
       // 3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,
        
 //       6.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,
//        3.0/8640.0,6.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,
 //       3.0/8640.0,3.0/8640.0,6.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,
  //      6.0/8640.0,6.0/8640.0,6.0/8640.0,12.0/8640.0,4.0/8640.0,4.0/8640.0,
   //     3.0/8640.0,3.0/8640.0,3.0/8640.0,4.0/8640.0,4.0/8640.0,2.0/8640.0,
    //    3.0/8640.0,3.0/8640.0,3.0/8640.0,4.0/8640.0,2.0/8640.0,4.0/8640.0,
        
 //       6.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,
  //      3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,
   //     3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,
    //    3.0/8640.0,3.0/8640.0,3.0/8640.0,4.0/8640.0,4.0/8640.0,2.0/8640.0,
     //   6.0/8640.0,6.0/8640.0,6.0/8640.0,4.0/8640.0,12.0/8640.0,4.0/8640.0,
       // 3.0/8640.0,3.0/8640.0,3.0/8640.0,2.0/8640.0,4.0/8640.0,4.0/8640.0,
        
//        6.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,
  //      3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,
  //      3.0/8640.0,3.0/8640.0,6.0/8640.0,3.0/8640.0,3.0/8640.0,6.0/8640.0,
  //      3.0/8640.0,3.0/8640.0,3.0/8640.0,4.0/8640.0,2.0/8640.0,4.0/8640.0,
  //      3.0/8640.0,3.0/8640.0,3.0/8640.0,2.0/8640.0,4.0/8640.0,4.0/8640.0,
  //      6.0/8640.0,6.0/8640.0,6.0/8640.0,4.0/8640.0,4.0/8640.0,12.0/8640.0
  //      
  //  };
    
    const double arg_sqrt = pow(R[0],2)+pow(R[1],2)+pow(R[2],2);
    const double n_R = sqrt(arg_sqrt);
    const double R2 = arg_sqrt;
    const double R3 = n_R*n_R*n_R;
    const double R5 = R3*n_R*n_R;
    const double R7 = R5*n_R*n_R;
    const double R9 = R7*n_R*n_R;
    double J_prod;
    
    // forces constants
    
    const double a_R1 = -3*Q1/R5*R[0];
    const double a_R2 = -3*Q1/R5*R[1];
    const double a_R3 = -3*Q1/R5*R[2];
    
    const double b = Q1/R3;
    
    
    // Torques constant
    
    
    double a_t_R1 = - Q1*R[0];
    double a_t_R2 = - Q1*R[1];
    double a_t_R3 = - Q1*R[2];

    
    
    double dU_dAt_x;
    double dU_dAt_y;
    double dU_dAt_z;
    
    double dU0_dAt_x = Q0/R3*R[0];
    double dU0_dAt_y = Q0/R3*R[1];
    double dU0_dAt_z = Q0/R3*R[2];
    
    double dU1_dAt_x;
    double dU1_dAt_y;
    double dU1_dAt_z;
    
    double dU2_dAt_x;
    double dU2_dAt_y;
    double dU2_dAt_z;
    
    double dU3_dAt_x;
    double dU3_dAt_y;
    double dU3_dAt_z;
    
    double dU1_dT_11;
    double dU1_dT_12;
    double dU1_dT_13;
    double dU1_dT_21;
    double dU1_dT_22;
    double dU1_dT_23;
    double dU1_dT_31;
    double dU1_dT_32;
    double dU1_dT_33;
    
    
    double dU2_dT_11;
    double dU2_dT_12;
    double dU2_dT_13;
    double dU2_dT_21;
    double dU2_dT_22;
    double dU2_dT_23;
    double dU2_dT_31;
    double dU2_dT_32;
    double dU2_dT_33;
    
    
    double dU_dT_11;
    double dU_dT_12;
    double dU_dT_13;
    double dU_dT_21;
    double dU_dT_22;
    double dU_dT_23;
    double dU_dT_31;
    double dU_dT_32;
    double dU_dT_33;
    
    
    double sum_v_r1;
    double sum_v_r2;
    double sum_v_r3;
    double sum_w;
    double p_sum_v_r1;
    double p_sum_v_r2;
    double p_sum_v_r3;
    
    
    
    
    double dU2_dAt_x_non_simmetric_term,dU2_dAt_y_non_simmetric_term,dU2_dAt_z_non_simmetric_term;
    double dU2_dAt_simmetric_term;
    
    double w[6],r_[9];
    
    
    
    double v[9];
    
    int l_m =3;
    int l_v = 3; // length matrix v
    int l_vt2 = 2*l_v;
    
    
    int i,ii,iii,j,k;
    
    double cdR2 = 5/R2;
    double mtmR5 = -3/(2*R5);
    double u_R3 = 1/R3;
    double R1t3dR2 = 3*R[0]/R2;
    double R2t3dR2 = 3*R[1]/R2;
    double R3t3dR2 = 3*R[2]/R2;
    
    double v_012tQ_012;
    double v_012tQ_678;
    double v_012tQ_121314;
    double v_012tQ_181920;
    double v_012tQ_242526;
    double v_012tQ_303132;
    
    double v_678tQ_012;
    double v_678tQ_678;
    double v_678tQ_121314;
    double v_678tQ_181920;
    double v_678tQ_242526;
    double v_678tQ_303132;
    
    double v_121314tQ_012;
    double v_121314tQ_678;
    double v_121314tQ_121314;
    double v_121314tQ_181920;
    double v_121314tQ_242526;
    double v_121314tQ_303132;
    
    double a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c2,c3,c4,c5,c6;
    
    
    
        
        
        for( iii = 0; iii < nfaces_a; iii = iii + 1 ) // nfaces_a
            
        {
            
            
            dU1_dAt_x = 0;
            dU1_dAt_y = 0;
            dU1_dAt_z = 0;
            
            dU2_dAt_x = 0;
            dU2_dAt_y = 0;
            dU2_dAt_z = 0;
            
            
            dU2_dT_11=0;
            dU2_dT_12=0;
            dU2_dT_13=0;
            dU2_dT_21=0;
            dU2_dT_22=0;
            dU2_dT_23=0;
            dU2_dT_31=0;
            dU2_dT_32=0;
            dU2_dT_33=0;
            

            v[0] = -r1a[iii*l_m];
            v[1] = -r2a[iii*l_m];
            v[2] = -r3a[iii*l_m];
            v[3] = -r1a[iii*l_m+1];
            v[4] = -r2a[iii*l_m+1];
            v[5] = -r3a[iii*l_m+1];
            v[6] = -r1a[iii*l_m+2];
            v[7] = -r2a[iii*l_m+2];
            v[8] = -r3a[iii*l_m+2];
            
            
            //v[3] = r1b[ii*l_m];
            //v[3] = r2b[ii*l_m];
            //v[4] = r3b[ii*l_m];
            //v[9] = r1b[ii*l_m+1];
            //v[10] = r2b[ii*l_m+1];
            //v[11] = r3b[ii*l_m+1];
            //v[15] = r1b[ii*l_m+2];
            //v[16] = r2b[ii*l_m+2];
            //v[17] = r3b[ii*l_m+2];
            
            //mexPrintf("V=%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n",v[0],v[1],v[2],v[3],v[3],v[4],v[3],v[4],v[5],v[9],v[10],v[11],v[6],v[7],v[8],v[15],v[16],v[17]);
            
            sum_v_r1 = v[0] + v[1] + v[2];
            sum_v_r2 = v[3] + v[4] + v[5];
            sum_v_r3 = v[6] + v[7] + v[8];
            p_sum_v_r1 = v[0] + v[1] + v[2];
            p_sum_v_r2 = v[3] + v[4] + v[5];
            p_sum_v_r3 = v[6] + v[7] + v[8];
           
            
            //for( j = 0; j < 6; j = j + 1 )
            //{
                //w[j]= R[0]*v[j]+R[1]*v[l_v+j]+R[2]*v[l_vt2+j];
            
                w[0]= R[0]*v[0]+R[1]*v[l_v]+R[2]*v[l_vt2];
                w[1]= R[0]*v[1]+R[1]*v[l_v+1]+R[2]*v[l_vt2+1];
                w[2]= R[0]*v[2]+R[1]*v[l_v+2]+R[2]*v[l_vt2+2];
            
 
            //};
            
            sum_w = w[0]+ w[1]+ w[2];
            
            
             //mexPrintf("R=%f\n%f\n%f",R[0],R[1],R[2]);
            //mexPrintf("w=%f\n%f\n%f\n%f\n%f\n%f",w[0],w[1],w[2],w[3],w[4],w[5]);
            
            
            // Torques
            
            dU1_dT_11 =  a_t_R1*p_sum_v_r1;
            dU1_dT_12 =  a_t_R1*p_sum_v_r2;
            dU1_dT_13 =  a_t_R1*p_sum_v_r3;
            dU1_dT_21 =  a_t_R2*p_sum_v_r1;
            dU1_dT_22 =  a_t_R2*p_sum_v_r2;
            dU1_dT_23 =  a_t_R2*p_sum_v_r3;
            dU1_dT_31 =  a_t_R3*p_sum_v_r1;
            dU1_dT_32 =  a_t_R3*p_sum_v_r2;
            dU1_dT_33 =  a_t_R3*p_sum_v_r3;
            
            //mexPrintf("dU1_dt_11=%30.29f\n",dU1_dT_11);
            //mexPrintf("dU1_dt_12=%30.29f\n",dU1_dT_12);
            //mexPrintf("dU1_dt_13=%30.29f\n",dU1_dT_13);
            //mexPrintf("dU1_dt_21=%30.29f\n",dU1_dT_21);
            //mexPrintf("dU1_dt_22=%30.29f\n",dU1_dT_22);
            //mexPrintf("dU1_dt_23=%30.29f\n",dU1_dT_23);
            //mexPrintf("dU1_dt_31=%30.29f\n",dU1_dT_31);
            //mexPrintf("dU1_dt_32=%30.29f\n",dU1_dT_32);
            //mexPrintf("dU1_dt_33=%30.29f\n",dU1_dT_33);

            
            
            r_[0]=   v[0]*v[0]+v[l_v]*v[l_v]    +v[l_vt2]*v[l_vt2];
            r_[1]=   v[0]*v[1]+v[l_v]*v[l_v+1]  +v[l_vt2]*v[l_vt2+1];
            r_[2]=   v[0]*v[2]+v[l_v]*v[l_v+2]  +v[l_vt2]*v[l_vt2+2];
           
            
            r_[3]=   v[1]*v[0]+v[l_v+1]*v[l_v]  +v[l_vt2+1]*v[l_vt2];
            r_[4]=   v[1]*v[1]+v[l_v+1]*v[l_v+1]+v[l_vt2+1]*v[l_vt2+1];
            r_[5]=   v[1]*v[2]+v[l_v+1]*v[l_v+2]+v[l_vt2+1]*v[l_vt2+2];
            
            
            r_[6]=  v[2]*v[0]+v[l_v+2]*v[l_v]  +v[l_vt2+2]*v[l_vt2];
            r_[7]=  v[2]*v[1]+v[l_v+2]*v[l_v+1]+v[l_vt2+2]*v[l_vt2+1];
            r_[8]=  v[2]*v[2]+v[l_v+2]*v[l_v+2]+v[l_vt2+2]*v[l_vt2+2];
            
            
            //mexPrintf("V=%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n%f\n",r_[0],r_[1],r_[2],r_[3],r_[4],r_[5],r_[3],r_[4],r_[5],r_[9],r_[10],r_[11],r_[6],r_[7],r_[8],r_[15],r_[16],r_[17],r_[18],r_[19],r_[20],r_[21],r_[22],r_[23],r_[24],r_[25],r_[26],r_[27],r_[28],r_[29],r_[30],r_[31],r_[32],r_[33],r_[34],r_[35]);
            
            
            
            dU2_dAt_x_non_simmetric_term = 2*(w[0]*(Q2[0]*v[0] +Q2[1]*v[1] +Q2[2]*v[2])+
                                              w[1]*(Q2[6]*v[0] +Q2[7]*v[1] +Q2[8]*v[2])+
                                              w[2]*(Q2[12]*v[0]+Q2[13]*v[1]+Q2[14]*v[2]));
                                              
            dU2_dAt_y_non_simmetric_term = 2*(w[0]*(Q2[0]*v[3] +Q2[1]*v[4] +Q2[2]*v[5] )+
                                              w[1]*(Q2[6]*v[3] +Q2[7]*v[4] +Q2[8]*v[5] )+
                                              w[2]*(Q2[12]*v[3]+Q2[13]*v[4]+Q2[14]*v[5]));
                                              
            dU2_dAt_z_non_simmetric_term = 2*(w[0]*(Q2[0]*v[6] +Q2[1]*v[7] +Q2[2]*v[8] )+
                                              w[1]*(Q2[6]*v[6] +Q2[7]*v[7] +Q2[8]*v[8] )+
                                              w[2]*(Q2[12]*v[6]+Q2[13]*v[7]+Q2[14]*v[8]));
                                           
            
            
            dU2_dAt_simmetric_term = Q2[0]*(r_[0]-cdR2*w[0]*w[0])+ //diagonal terms
                                     Q2[7]*(r_[4]-cdR2*w[1]*w[1])+
                                     Q2[14]*(r_[8]-cdR2*w[2]*w[2])+
                                    
            
                                   2*(Q2[1]*(r_[1]-cdR2*w[0]*w[1])+
                                      Q2[2]*(r_[2]-cdR2*w[0]*w[2])+
                                     
               
                                      Q2[8]*(r_[5]-cdR2*w[1]*w[2]));
            
            
            dU1_dAt_x = a_R1*sum_w+b*sum_v_r1;
            dU1_dAt_y = a_R2*sum_w+b*sum_v_r2;
            dU1_dAt_z = a_R3*sum_w+b*sum_v_r3;
            
            dU2_dAt_x = mtmR5*(dU2_dAt_x_non_simmetric_term + R[0]* dU2_dAt_simmetric_term);
            dU2_dAt_y = mtmR5*(dU2_dAt_y_non_simmetric_term + R[1]* dU2_dAt_simmetric_term);
            dU2_dAt_z = mtmR5*(dU2_dAt_z_non_simmetric_term + R[2]* dU2_dAt_simmetric_term);
            
            // second order torques
            
            v_012tQ_012 =    (v[0]*Q2[0] +v[1]*Q2[1] +v[2]*Q2[2]);
            v_012tQ_678 =    (v[0]*Q2[6] +v[1]*Q2[7] +v[2]*Q2[8]);
            v_012tQ_121314 = (v[0]*Q2[12]+v[1]*Q2[13]+v[2]*Q2[14]);
            v_012tQ_181920 = (v[0]*Q2[18]+v[1]*Q2[19]+v[2]*Q2[20]);
            v_012tQ_242526 = (v[0]*Q2[24]+v[1]*Q2[25]+v[2]*Q2[26]);
            v_012tQ_303132 = (v[0]*Q2[30]+v[1]*Q2[31]+v[2]*Q2[32]);
            
            v_678tQ_012 =    (v[3]*Q2[0] +v[4]*Q2[1] +v[5]*Q2[2]);
            v_678tQ_678 =    (v[3]*Q2[6] +v[4]*Q2[7] +v[5]*Q2[8]);
            v_678tQ_121314 = (v[3]*Q2[12]+v[4]*Q2[13]+v[5]*Q2[14]);
            v_678tQ_181920 = (v[3]*Q2[18]+v[4]*Q2[19]+v[5]*Q2[20]);
            v_678tQ_242526 = (v[3]*Q2[24]+v[4]*Q2[25]+v[5]*Q2[26]);
            v_678tQ_303132 = (v[3]*Q2[30]+v[4]*Q2[31]+v[5]*Q2[32]);
            
            v_121314tQ_012 =    (v[6]*Q2[0] +v[7]*Q2[1] +v[8]*Q2[2]);
            v_121314tQ_678 =    (v[6]*Q2[6] +v[7]*Q2[7] +v[8]*Q2[8]);
            v_121314tQ_121314 = (v[6]*Q2[12]+v[7]*Q2[13]+v[8]*Q2[14]);
            v_121314tQ_181920 = (v[6]*Q2[18]+v[7]*Q2[19]+v[8]*Q2[20]);
            v_121314tQ_242526 = (v[6]*Q2[24]+v[7]*Q2[25]+v[8]*Q2[26]);
            v_121314tQ_303132 = (v[6]*Q2[30]+v[7]*Q2[31]+v[8]*Q2[32]);
            
            a1 = R1t3dR2*w[0]-v[0];
            a2 = R1t3dR2*w[1]-v[1];
            a3 = R1t3dR2*w[2]-v[2];
                       
            b1 = R2t3dR2*w[0]-v[3];
            b2 = R2t3dR2*w[1]-v[4];
            b3 = R2t3dR2*w[2]-v[5];
           
            
            c1 = R3t3dR2*w[0]-v[6];
            c2 = R3t3dR2*w[1]-v[7];
            c3 = R3t3dR2*w[2]-v[8];
           


            
            dU2_dT_11 =  v_012tQ_012*(a1)+v_012tQ_678*(a2)+v_012tQ_121314*(a3);
            
            dU2_dT_12 =  v_678tQ_012*(a1)+v_678tQ_678*(a2)+v_678tQ_121314*(a3);
            
            dU2_dT_13 = v_121314tQ_012*(a1)+v_121314tQ_678*(a2)+v_121314tQ_121314*(a3);
                                            
            dU2_dT_21 =  v_012tQ_012*(b1)+v_012tQ_678*(b2)+v_012tQ_121314*(b3);
                          
            dU2_dT_22 =  v_678tQ_012*(b1)+v_678tQ_678*(b2)+v_678tQ_121314*(b3);
                          
            dU2_dT_23 =  v_121314tQ_012*(b1)+v_121314tQ_678*(b2)+v_121314tQ_121314*(b3);
                          
            dU2_dT_31 =  v_012tQ_012*(c1)+v_012tQ_678*(c2)+v_012tQ_121314*(c3);
                          
            dU2_dT_32 =  v_678tQ_012*(c1)+v_678tQ_678*(c2)+v_678tQ_121314*(c3);
                          
            dU2_dT_33 =  v_121314tQ_012*(c1)+v_121314tQ_678*(c2)+v_121314tQ_121314*(c3);
                          
            

            
            //mexPrintf("dU2_dt_11=%30.29f\n",dU2_dT_11);
            //mexPrintf("dU2_dt_12=%30.29f\n",dU2_dT_12);
            //mexPrintf("dU2_dt_13=%30.29f\n",dU2_dT_13);
            //mexPrintf("dU2_dt_21=%30.29f\n",dU2_dT_21);
            //mexPrintf("dU2_dt_22=%30.29f\n",dU2_dT_22);
            //mexPrintf("dU2_dt_23=%30.29f\n",dU2_dT_23);
            //mexPrintf("dU2_dt_31=%30.29f\n",dU2_dT_31);
            //mexPrintf("dU2_dt_32=%30.29f\n",dU2_dT_32);
            //mexPrintf("dU2_dt_33=%30.29f\n",dU2_dT_33);
            
                              
                              
                              
            

            dU_dAt_x = dU0_dAt_x + dU1_dAt_x + dU2_dAt_x; //+ dU3_dAt_x;
            dU_dAt_y = dU0_dAt_y + dU1_dAt_y + dU2_dAt_y; //+ dU3_dAt_y;
            dU_dAt_z = dU0_dAt_z + dU1_dAt_z + dU2_dAt_z; //+ dU3_dAt_z;
            
            dU_dT_11= dU1_dT_11 + dU2_dT_11;
            dU_dT_12= dU1_dT_12 + dU2_dT_12;
            dU_dT_13= dU1_dT_13 + dU2_dT_13;
            dU_dT_21= dU1_dT_21 + dU2_dT_21;
            dU_dT_22= dU1_dT_22 + dU2_dT_22;
            dU_dT_23= dU1_dT_23 + dU2_dT_23;
            dU_dT_31= dU1_dT_31 + dU2_dT_31;
            dU_dT_32= dU1_dT_32 + dU2_dT_32;
            dU_dT_33= dU1_dT_33 + dU2_dT_33;
            
            //mexPrintf("dU_dt_11=%30.29f\n",dU_dT_11);
            //mexPrintf("dU_dt_12=%30.29f\n",dU_dT_12);
            //mexPrintf("dU_dt_13=%30.29f\n",dU_dT_13);
            //mexPrintf("dU_dt_21=%30.29f\n",dU_dT_21);
            //mexPrintf("dU_dt_22=%30.29f\n",dU_dT_22);
            //mexPrintf("dU_dt_23=%30.29f\n",dU_dT_23);
            //mexPrintf("dU_dt_31=%30.29f\n",dU_dT_31);
            //mexPrintf("dU_dt_32=%30.29f\n",dU_dT_32);
            //mexPrintf("dU_dt_33=%30.29f\n",dU_dT_33);
            
            J_prod = Ja[iii];


            
            
            Fa[0] = Fa[0] +  J_prod*dU_dAt_x;
            Fa[1] = Fa[1] +  J_prod*dU_dAt_y;
            Fa[2] = Fa[2] +  J_prod*dU_dAt_z;
            
            E_m[0] = E_m[0] +  J_prod*dU_dT_11;  // it works, it is the traspost
            E_m[3] = E_m[3] +  J_prod*dU_dT_12;
            E_m[6] = E_m[6] +  J_prod*dU_dT_13;
            E_m[1] = E_m[1] +  J_prod*dU_dT_21;
            E_m[4] = E_m[4] +  J_prod*dU_dT_22;
            E_m[7] = E_m[7] +  J_prod*dU_dT_23;
            E_m[2] = E_m[2] +  J_prod*dU_dT_31;
            E_m[5] = E_m[5] +  J_prod*dU_dT_32;
            E_m[8] = E_m[8] +  J_prod*dU_dT_33;
            
            
            
        };
    
    
    Fa[0] = Grav_con*6*rho_a*M_b*Fa[0];
    Fa[1] = Grav_con*6*rho_a*M_b*Fa[1];
    Fa[2] = Grav_con*6*rho_a*M_b*Fa[2];
    
    E_m[0] = Grav_con*6*rho_a*M_b/R3*E_m[0];  // it works, it is the traspost
    E_m[3] = Grav_con*6*rho_a*M_b/R3*E_m[3];
    E_m[6] = Grav_con*6*rho_a*M_b/R3*E_m[6];
    E_m[1] = Grav_con*6*rho_a*M_b/R3*E_m[1];
    E_m[4] = Grav_con*6*rho_a*M_b/R3*E_m[4];
    E_m[7] = Grav_con*6*rho_a*M_b/R3*E_m[7];
    E_m[2] = Grav_con*6*rho_a*M_b/R3*E_m[2];
    E_m[5] = Grav_con*6*rho_a*M_b/R3*E_m[5];
    E_m[8] = Grav_con*6*rho_a*M_b/R3*E_m[8];
    
    
};


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray  *prhs[] ) {
    
    
    
    //declare variables
    mxArray *R_in_m,*r1a_in_m,*r2a_in_m,*r3a_in_m,*Ja_in_m,*rho_a_in_m,*M_b_in_m,*nfaces_a_in_m,*Fa_out_m,*E_out_m;
    //const mwSize *dims;
    double *R,*r1a,*r2a,*r3a,*Ja;
    double rho_a,M_b;
    double *Fa,*E_m;
    int nfaces_a,dimx, dimy, numdims;;
    
    
    //figure out dimensions
    //dims = mxGetDimensions(prhs[0]);
    //numdims = mxGetNumberOfDimensions(prhs[0]);
    //dimy = (int)dims[0]; dimx = (int)dims[1];
    
    //associate outputs
    plhs[0] = mxCreateDoubleMatrix(3,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3,3,mxREAL);
    Fa=mxGetPr(plhs[0]);
    E_m = mxGetPr(plhs[1]);
    
    //associate pointers at the outputs and inputs
    
    // Input
    R=mxGetPr(prhs[0]);
    r1a=mxGetPr(prhs[1]);
    r2a=mxGetPr(prhs[2]);
    r3a=mxGetPr(prhs[3]);
    Ja=mxGetPr(prhs[4]);
    rho_a=mxGetScalar(prhs[5]);
    M_b=mxGetScalar(prhs[6]);
    nfaces_a=mxGetScalar(prhs[7]);
   
    
    // check print
    
    //mexPrintf("rho_a=%f\n",rho_a);
    //mexPrintf("M_b=%f\n",M_b);
    //mexPrintf("R=%f\n",R[0]);
    //mexPrintf("nfaces_a=%d\n",nfaces_a);
    //mexPrintf("nfaces_bf=%f\n",nfaces_b);
    //mexPrintf("r1a=%20.3f\n",r1a[0]);
    //mexPrintf("r1b1=%20.3f\n",r1b[0]);
    //mexPrintf("r2b1=%20.3f\n",r2b[0]);
    //mexPrintf("r3b1=%20.3f\n",r3b[0]);
    
    action_calculator(R,r1a,r2a,r3a,Ja,rho_a,M_b,nfaces_a,Fa,E_m);
    
    return;
    }
    
