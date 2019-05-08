//  ast_gravity_fun_grad.c
//
//
//  Created by Lorenzo Teofili on 14/10/15.
//
//
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// ast_gravity_fun_grad computes the ADIMENTIONAL gravitational
// force acting on a point in space due to the
// presence of a body of whatever shape. Further, the gradient
// of such a force is computed.
// Consider as a reference the model outlined by
// ROBERT A. WERNER and DANIEL J. SCHEERES in the paper:
//
// "Exterior gravitation of a polyhedron derived and compared
//  with harmonic and mascon gravitation representations of asteroid
//  4769 Castalia".
//
// The model is exact up to the precision of the mesh rapresentation of the asteorid
// and the machine precision.
//
//
// [F, E_m] = ast_gravity_fun(rho_adim, vert_adim, faces, rf_adim, nfaces, R_adim)
//
// INPUT:
// rho          - Scalar, asteroid density, the asteroid is suppose to have
//                constant density.
//
// vert         - A vector obtained by the function "load_model_2bp", it takes into account for
//                the shape of the bodies.
//
// faces        - A vector obtained by the function "load_model_2bp", it takes into account for
//                the shape of the bodies.
//
// rf_          - A matrix (nfaces x 3) obtained by the function "load_model_2bp", it takes into account for
//                the shape of the bodies.
//
//
// Ja           - Vector obtained by the function "Asteroid_Inertia_matrix_2".
//                It represents the vector of the determinant of the polihedra
//                of the asteroids.
//
// R            - Vector which indicates the position of the point where the forces and its gradient
//                have to be computed. It must be expressed in a reference frame fixed with the asteroid,
//                centered in its baricenter.
//
// OUTPUT:
// out          - This is a 12 components vector, the first 9 components are the gravity force gradient, the last
//                3 components are the force acting on the selected point.
//
// N.B.           TO MAKE THE OUTPUT DIMENTIONAL THE FIRST 9 COMPONENT HAVE TO BE MULTIPLIED FOR THE SCALE DENSITY
//                (rho_ref) THE LAST 3 ONES FOR THE PRODUCT BETWEEN THE SCALE DENSITY TIMES THE SCALE DISTANCE (rho_ref*r_ref).
//
//
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




#include "mex.h"
#include "matrix.h"
#include "math.h"
		
double doto(double *b, double *c){
    double a;
    a = b[0]*c[0]+b[1]*c[1]+b[2]*c[2];
    return a;
}
void array_sum(double *b, double *c,double *a){
    a[0] = b[0]+c[0];
    a[1] = b[1]+c[1];
    a[2] = b[2]+c[2];
    return;
}

void array_diff(double *b, double *c,double *a){
    a[0] = b[0]-c[0];
    a[1] = b[1]-c[1];
    a[2] = b[2]-c[2];
    return;
}



void assign_vect(double *y, mwSize m, mwSize n, double *a){
    {
        mwSize i;
        for (i=0; i<=(n-m); i++) {
            a[i] = y[i+m];

        }
        
        return;
    }
}

double norma(double *b){
    double a;
    a = sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);
    return a;
}

void matrix_trans(double *r1, double *a){
    {
        a[0] = r1[0];
        a[1] = r1[3];
        a[2] = r1[6];
        a[3] = r1[1];
        a[4] = r1[4];
        a[5] = r1[7];
        a[6] = r1[2];
        a[7] = r1[5];
        a[8] = r1[8];
        return;
    }
}


void matrix_vect(double *r1, double *r2, double *a){
    {
        mwSize j, i;
        for (j=0; j<=2; j++) {
            for (i=0; i<=2; i++) {
                a[j] = a[j] + r1[i+3*(j)] * r2[i];
                
            }
        }
        return;
    }
}

void arrayProduct(double x, double *y, double *z, mwSize n){
    mwSize i;

    for (i=0; i<n; i++) {
        z[i] = x * y[i];
    }
}

void crossi(double *h, double *r2, double *a){
    {   
        double rb_m[9];
        rb_m[0] = 0.;
        rb_m[1] =  -h[2];
        rb_m[2] =  h[1];
        rb_m[3] = h[2];
        rb_m[4]	= 0.;
        rb_m[5] = -h[0];
        rb_m[6] = -h[1];
        rb_m[7]	= h[0];
        rb_m[8] = 0;
        matrix_vect(rb_m, r2, a);
        return;
    }
}


void ast_gravity_fun(double rho, double *vert, double *faces, double *rf_, int nfaces,double r[3],double out[3]){

    double rr[3]={0.};
    int *gg;               
                          
	int     i, j, m, n;  
 //   const double  G0= 6.67259*pow(10.,-11.);
                        
    double  uno[3]={0.}; 
	double  r1_[3]={0.};
	double  r2_[3]={0.};
	double  r3_[3]={0.}; 
    double  r1[3]={0.};
	double  r2[3]={0.};
	double  r3[3]={0.}; 
	double  aa[3]={0.};
	double  bb[3]={0.};
	double  cc[3]={0.};
	double adot;
	double bdot;
    double	cdot;
	double	ddot;
	double Le, lij, r1_norm, r2_norm, r3_norm,ri_norm,rj_norm;
    double nf[3]={0.}, rf[3]={0.};
    double  g[3]={0.},ri[3]={0.},rj[3]={0.},ri_[3]={0.},rj_[3]={0.},nijf[3]={0.},dlij[3]={0.}; 
    double ss[3], vv[3],re[3];	
	double  ri_v [9]={0.}; 
	double  riv [9], rf_i[3]={0.}; 
	double Ff[9]={0.},Eef[9]={0.};
	double volume = 0.;
    double  wf=0., aux1[3]={0.};

    mwSize nr, nc, ii1, ii2, ii3;      
   double gigi1=0.,gigi2=0.,gigi3=0.;    
   
    rr[0] = r[0];
	rr[1] = r[1];
	rr[2] = r[2];

 for (i=1; i<=nfaces; i++){

     ii1 = 3*(i-1);
     ii2 = ii1+1;
     ii3 = ii1+2;
     
     ii1 = faces[ii1]-1;
    nc = 3*(ii1);

    r1[0] = vert[nc];
    nc = nc+1;

    r1[1] =	vert[nc];
    nc = nc+1;

	r1[2] =vert[nc];

    ii2 = faces[ii2]-1;
    nc = 3*(ii2);
    r2[0] = vert[nc];
    nc = nc+1;
	r2[1] = vert[nc];
	nc = nc+1;
    r2[2] = vert[nc];
    ii3 = faces[ii3]-1;
    nc = 3*(ii3);
    r3[0] = vert[nc];
	nc = nc +1;
    r3[1] = vert[nc];
	nc = nc+1;
    r3[2] =vert[nc];

    r1_[0] = r1[0] - rr[0];
	r1_[1] = r1[1] - rr[1];
	r1_[2] = r1[2] - rr[2];
    
    r2_[0] = r2[0] - rr[0];
	r2_[1] = r2[1] - rr[1];
	r2_[2] = r2[2] - rr[2];
    
    r3_[0] = r3[0] - rr[0];
	r3_[1] = r3[1] - rr[1];
	r3_[2] = r3[2] - rr[2];
    
    rf[0] = rf_[3*(i-1)]-rr[0];
	rf[1] = rf_[3*(i-1)+1]-rr[1];
	rf[2] = rf_[3*(i-1)+2]-rr[2] ; 
    
	aa[0] = r1[1]*r2[2]-r1[2]*r2[1];
    aa[1] = r1[2]*r2[0]-r1[0]*r2[2];
    aa[2] = r1[0]*r2[1]-r1[1]*r2[0];
    
	bb[0] = r2[1]*r3[2]-r2[2]*r3[1];
    bb[1] = r2[2]*r3[0]-r2[0]*r3[2];
    bb[2] = r2[0]*r3[1]-r2[1]*r3[0];
    
	cc[0] = r3[1]*r1[2]-r3[2]*r1[1];
    cc[1] = r3[2]*r1[0]-r3[0]*r1[2];
    cc[2] = r3[0]*r1[1]-r3[1]*r1[0];

    nf[0] = aa[0] + bb[0] + cc[0]; 
	nf[1] = aa[1] + bb[1] + cc[1];
	nf[2] = aa[2] + bb[2] + cc[2];

	aa[0] = rf[0]+rr[0];
	aa[1] = rf[1]+rr[1];
	aa[2] = rf[2]+rr[2];
	adot =doto(nf,aa);

	adot =norma(nf);

	nf[0] = nf[0]/adot;
    nf[1] = nf[1]/adot;
	nf[2] = nf[2]/adot;

    if (doto(nf,r1)<0){
        nf[0] = -nf[0];
        nf[1] = -nf[1];
        nf[2] = -nf[2];
    }

    ri_v[0] = r1_[0];
	ri_v[1] = r1_[1];
	ri_v[2] = r1_[2];
	ri_v[3] = r2_[0];
	ri_v[4] = r2_[1];
	ri_v[5] = r2_[2];
    ri_v[6] = r3_[0];
    ri_v[7] = r3_[1];
	ri_v[8] = r3_[2];
    riv[0] = r1[0];
	riv[1] = r1[1];
	riv[2] = r1[2];
	riv[3] = r2[0];
	riv[4] = r2[1];
	riv[5] = r2[2];
    riv[6] = r3[0];
    riv[7] = r3[1];
	riv[8] = r3[2];

    r1_norm= norma(r1_);
    r2_norm=norma(r2_);
    r3_norm=norma(r3_); 

	ss[0] = r2_[0];
	ss[1] = r2_[1];
	ss[2] = r2_[2];
	vv[0] = r3_[0];
	vv[1] = r3_[1];
	vv[2] = r3_[2];
	aa[0] = ss[1]*vv[2]-ss[2]*vv[1];
    aa[1] = ss[2]*vv[0]-ss[0]*vv[2];
    aa[2] = ss[0]*vv[1]-ss[1]*vv[0];
	adot = doto(r1_,aa);
    bdot = doto(r2_,r3_);
	cdot = doto(r3_,r1_);
	ddot = doto(r1_,r2_);

    wf = 2.*atan(adot/((r1_norm*r2_norm*r3_norm) + (r1_norm*bdot) + (r2_norm*cdot) + (r3_norm*ddot)));

    Ff[0] = nf[0]*nf[0];
    Ff[1] = nf[0]*nf[1];
    Ff[2] = nf[0]*nf[2];
	Ff[3] = nf[1]*nf[0]; 
	Ff[4] = nf[1]*nf[1]; 
	Ff[5] = nf[1]*nf[2]; 
	Ff[6] = nf[2]*nf[0]; 
	Ff[7] = nf[2]*nf[1]; 
	Ff[8] = nf[2]*nf[2];

	aa[0] = Ff[0]*rf[0] + Ff[1]*rf[1] + Ff[2]*rf[2];
    aa[1] = Ff[3]*rf[0] + Ff[4]*rf[1] + Ff[5]*rf[2];
    aa[2] = Ff[6]*rf[0] + Ff[7]*rf[1] + Ff[8]*rf[2];
    
	gigi1=gigi1+ wf*aa[0];
	gigi2= gigi2+wf*aa[1];
	gigi3 = gigi3 +wf*aa[2];
     
     // computing the gravitational gradient acting on a point (first part)
     
     out[0]= out[0]+ wf*Ff[0];
     out[1]= out[1]+ wf*Ff[1];
     out[2]= out[2]+ wf*Ff[2];
     out[3]= out[3]+ wf*Ff[3];
     out[4]= out[4]+ wf*Ff[4];
     out[5]= out[5]+ wf*Ff[5];
     out[6]= out[6]+ wf*Ff[6];
     out[7]= out[7]+ wf*Ff[7];
     out[8]= out[8]+ wf*Ff[8];

	for (j=1; j<=3; j++){
        
            ri[0] = riv[3*(j-1)]; 
            ri[1] = riv[3*(j-1)+1]; 
            ri[2] = riv[3*(j-1)+2]; 
			ri_[0] = ri_v[3*(j-1)];
            ri_[1] = ri_v[3*(j-1)+1];
            ri_[2] = ri_v[3*(j-1)+2];

            if (j ==3){
            rj[0] = riv[0];
            rj[1] = riv[1];
            rj[2] = riv[2];
			rj_[0] = ri_v[0];
            rj_[1] = ri_v[1];
            rj_[2] = ri_v[2];}
            else{
             rj[0] = riv[3*(j)]; 
             rj[1] = riv[3*(j)+1];
             rj[2] = riv[3*(j)+2];
			 rj_[0] = ri_v[3*(j)];
             rj_[1] = ri_v[3*(j)+1];
             rj_[2] = ri_v[3*(j)+2];
            }

        ri_norm = norma(ri_);
        rj_norm = norma(rj_);

        dlij[0] = rj[0] - ri[0];
        dlij[1] = rj[1] - ri[1];
        dlij[2] = rj[2] - ri[2];
        lij = norma(dlij);

        	ss[0] =dlij[0];
	ss[1] = dlij[1];
	ss[2] = dlij[2];
	vv[0] = nf[0];
	vv[1] = nf[1];
	vv[2] = nf[2];

	aa[0] = ss[1]*vv[2]-ss[2]*vv[1];
    aa[1] = ss[2]*vv[0]-ss[0]*vv[2];
    aa[2] = ss[0]*vv[1]-ss[1]*vv[0];
    adot = norma(aa);
    nijf[0] = aa[0]/adot;
    nijf[1] = aa[1]/adot;
    nijf[2] = aa[2]/adot;

    Eef[0] = nijf[0]*nf[0];
    Eef[1] = nijf[0]*nf[1];
    Eef[2] = nijf[0]*nf[2];
	Eef[3] = nijf[1]*nf[0]; 
	Eef[4] = nijf[1]*nf[1]; 
	Eef[5] = nijf[1]*nf[2]; 
	Eef[6] = nijf[2]*nf[0]; 
	Eef[7] = nijf[2]*nf[1]; 
	Eef[8] = nijf[2]*nf[2];     

        Le = log((ri_norm + rj_norm + lij)/(ri_norm + rj_norm - lij));

        re[0] = (ri_[0]  + rj_[0] )/2.;
        re[1]  = (ri_[1]  + rj_[1] )/2.;
        re[2]  = (ri_[2]  + rj_[2] )/2.;


        gigi1=gigi1 - Le*(Eef[0]*re[0] + Eef[1]*re[1] + Eef[2]*re[2]);
        gigi2=gigi2 - Le*(Eef[3]*re[0] + Eef[4]*re[1] + Eef[5]*re[2]);
        gigi3= gigi3 - Le*(Eef[6]*re[0] + Eef[7]*re[1] + Eef[8]*re[2]);
        
        // Second Part of the Gravitational Gradient
        
        out[0]= out[0] - Le*Eef[0];
        out[1]= out[1] - Le*Eef[1];
        out[2]= out[2] - Le*Eef[2];
        out[3]= out[3] - Le*Eef[3];
        out[4]= out[4] - Le*Eef[4];
        out[5]= out[5] - Le*Eef[5];
        out[6]= out[6] - Le*Eef[6];
        out[7]= out[7] - Le*Eef[7];
        out[8]= out[8] - Le*Eef[8];

    }
 }

    out[0]= out[0]*rho;
    out[1]= out[1]*rho;
    out[2]= out[2]*rho;
    out[3]= out[3]*rho;
    out[4]= out[4]*rho;
    out[5]= out[5]*rho;
    out[6]= out[6]*rho;
    out[7]= out[7]*rho;
    out[8]= out[8]*rho;
    out[9]= -rho*gigi1;
    out[10]= -rho*gigi2;
    out[11] = -rho*gigi3;

return;
		}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{

  
    
    double rho;       
	double *faces;
    double *vert;          
	double *rf_;
	int nfaces;     
    double *r;
    int *gg;         
    double *out;        
    if(nrhs != 6)
        mexErrMsgTxt("gravity_force requires 6 and only 6 inputs.");
    if(nlhs > 1)
        mexErrMsgTxt("Too many outputs. gravity_force returns only 1 output.");
	
	if( !mxIsDouble(prhs[0]) || 
         mxIsComplex(prhs[0]) ||
         mxGetNumberOfElements(prhs[0])!=1 ) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notScalar","Input multiplier must be a scalar.");
    }
    rho = mxGetScalar(prhs[0]);     
	vert = mxGetPr(prhs[1]);
	faces = mxGetPr(prhs[2]);

	rf_ = mxGetPr(prhs[3]);
	nfaces = mxGetScalar(prhs[4]);
	r = mxGetPr(prhs[5]);

   plhs[0] = mxCreateDoubleMatrix(1,12,mxREAL);
	
	 
    out = mxGetPr(plhs[0]);
    ast_gravity_fun(rho, vert, faces, rf_, nfaces,r,out);
return;
		}

