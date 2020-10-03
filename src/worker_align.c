#include <stdio.h>
#include <string.h>
#include <math.h>
#include "defines.h"
#include "worker_align.h"
#include "module_linux_simple.h"
#include "system_cl.h"
#include "module_cl.h"
#include "device_cl.h"
#include "driver_cl.h"

#define CORR 0
#if CORR
		data_size rmean =0.0,tmean = 0.0;
		data_size std = 0.0;
#endif

static int worker_align_start(struct worker*);
static void worker_align_process(struct worker*, struct module*, struct frame*); 
static int worker_align_stop(struct worker*);
static int worker_align_module_cl_bld(struct worker*, struct module*, struct system*);
static int worker_align_module_cl_dty(struct worker*, struct module*); 
static int worker_align_module_linux_simple_bld(struct worker*, struct module*, 
		struct system*);
static int worker_align_module_linux_simple_dty(struct worker*, struct module*); 
static void worker_align_prepare_buffers(struct worker*, struct module*, struct frame**);
static void worker_align_process_buffers(struct worker*);
struct worker_interface worker_align_interface = 
{
	.start = worker_align_start,
	.stop = worker_align_stop,
	.prepare_buffers = worker_align_prepare_buffers,
	.process_buffers = worker_align_process_buffers,
	.dty = worker_align_dty,
	.module_linux_simple_bld = worker_align_module_linux_simple_bld,
	.module_linux_simple_dty = worker_align_module_linux_simple_dty,
	.module_cl_bld = worker_align_module_cl_bld,
	.module_cl_dty = worker_align_module_cl_dty
};

int worker_align_bld(struct worker* self, struct device* host)
{
	int ret = worker_bld(self, &worker_align_interface, host);
	self->properties = malloc(sizeof(struct worker_align_properties));
	((struct worker_align_properties*)self->properties)->frmsd = fopen("qalign.out","w");
	return ret;
}

void worker_align_dty(struct worker* self)
{
	fclose( ((struct worker_align_properties*)self->properties)->frmsd );
	((struct worker_align_properties*)self->properties)->frmsd = NULL;
	int in_used[] = {2,3,-1};
	int out_used[] = {1,-1};
	worker_clean_buffers(self,self->input,in_used);
	worker_clean_buffers(self,self->output,out_used);

}
static int worker_align_start(struct worker* self) { return WORKER_UNIMPLEMENTED; }
static void worker_align_process(struct worker* self, struct module* module, struct frame* frame) 
{
	return; 
}
static int worker_align_stop(struct worker* self) { return WORKER_UNIMPLEMENTED; }

int
FastCalcRMSDAndRotation(double *rot, double *A, double *rmsd, double E0, int len, double minScore)
{
    double Sxx, Sxy, Sxz, Syx, Syy, Syz, Szx, Szy, Szz;
    double Szz2, Syy2, Sxx2, Sxy2, Syz2, Sxz2, Syx2, Szy2, Szx2,
           SyzSzymSyySzz2, Sxx2Syy2Szz2Syz2Szy2, Sxy2Sxz2Syx2Szx2,
           SxzpSzx, SyzpSzy, SxypSyx, SyzmSzy,
           SxzmSzx, SxymSyx, SxxpSyy, SxxmSyy;
    double C[4];
    int i;
    double mxEigenV; 
    double oldg = 0.0;
    double b, a, delta, rms, qsqr;
    double q1, q2, q3, q4, normq;
    double a11, a12, a13, a14, a21, a22, a23, a24;
    double a31, a32, a33, a34, a41, a42, a43, a44;
    double a2, x2, y2, z2; 
    double xy, az, zx, ay, yz, ax; 
    double a3344_4334, a3244_4234, a3243_4233, a3143_4133,a3144_4134, a3142_4132; 
    double evecprec = 1e-6;
    double evalprec = 1e-11;

    Sxx = A[0]; Sxy = A[1]; Sxz = A[2];
    Syx = A[3]; Syy = A[4]; Syz = A[5];
    Szx = A[6]; Szy = A[7]; Szz = A[8];

    Sxx2 = Sxx * Sxx;
    Syy2 = Syy * Syy;
    Szz2 = Szz * Szz;

    Sxy2 = Sxy * Sxy;
    Syz2 = Syz * Syz;
    Sxz2 = Sxz * Sxz;

    Syx2 = Syx * Syx;
    Szy2 = Szy * Szy;
    Szx2 = Szx * Szx;

    SyzSzymSyySzz2 = 2.0*(Syz*Szy - Syy*Szz);
    Sxx2Syy2Szz2Syz2Szy2 = Syy2 + Szz2 - Sxx2 + Syz2 + Szy2;

    C[2] = -2.0 * (Sxx2 + Syy2 + Szz2 + Sxy2 + Syx2 + Sxz2 + Szx2 + Syz2 + Szy2);
    C[1] = 8.0 * (Sxx*Syz*Szy + Syy*Szx*Sxz + Szz*Sxy*Syx - Sxx*Syy*Szz - Syz*Szx*Sxy - Szy*Syx*Sxz);

    SxzpSzx = Sxz + Szx;
    SyzpSzy = Syz + Szy;
    SxypSyx = Sxy + Syx;
    SyzmSzy = Syz - Szy;
    SxzmSzx = Sxz - Szx;
    SxymSyx = Sxy - Syx;
    SxxpSyy = Sxx + Syy;
    SxxmSyy = Sxx - Syy;
    Sxy2Sxz2Syx2Szx2 = Sxy2 + Sxz2 - Syx2 - Szx2;

    C[0] = Sxy2Sxz2Syx2Szx2 * Sxy2Sxz2Syx2Szx2
         + (Sxx2Syy2Szz2Syz2Szy2 + SyzSzymSyySzz2) * (Sxx2Syy2Szz2Syz2Szy2 - SyzSzymSyySzz2)
         + (-(SxzpSzx)*(SyzmSzy)+(SxymSyx)*(SxxmSyy-Szz)) * (-(SxzmSzx)*(SyzpSzy)+(SxymSyx)*(SxxmSyy+Szz))
         + (-(SxzpSzx)*(SyzpSzy)-(SxypSyx)*(SxxpSyy-Szz)) * (-(SxzmSzx)*(SyzmSzy)-(SxypSyx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzpSzy)+(SxzpSzx)*(SxxmSyy+Szz)) * (-(SxymSyx)*(SyzmSzy)+(SxzpSzx)*(SxxpSyy+Szz))
         + (+(SxypSyx)*(SyzmSzy)+(SxzmSzx)*(SxxmSyy-Szz)) * (-(SxymSyx)*(SyzpSzy)+(SxzmSzx)*(SxxpSyy-Szz));

    /* Newton-Raphson */
    mxEigenV = E0;
    for (i = 0; i < 128; ++i)
    {
        oldg = mxEigenV;
        x2 = mxEigenV*mxEigenV;
        b = (x2 + C[2])*mxEigenV;
        a = b + C[1];
        delta = ((a*mxEigenV + C[0])/(2.0*x2*mxEigenV + b + a));
        mxEigenV -= delta;
         printf("\n diff[%3d]: %16g %16g %16g", i, mxEigenV - oldg, evalprec*mxEigenV, mxEigenV); 
        if (fabs(mxEigenV - oldg) < fabs(evalprec*mxEigenV))
            break;
    }

    if (i == 128 )
       fprintf(stderr,"\nMore than %d iterations needed!\n", i);

    /* the fabs() is to guard against extremely small, but *negative* numbers due to floating point error */
    rms = sqrt(fabs(2.0 * (E0 - mxEigenV)));
    (*rmsd) = rms;
     printf("\n\n %.4f %16g %16g \n", rms, E0, 2.0 * (E0 - mxEigenV)/len); 

    if (minScore > 0) 
        if (rms < minScore)
            return (-1); // Don't bother with rotation. 

    a11 = SxxpSyy + Szz-mxEigenV; a12 = SyzmSzy; a13 = - SxzmSzx; a14 = SxymSyx;
    a21 = SyzmSzy; a22 = SxxmSyy - Szz-mxEigenV; a23 = SxypSyx; a24= SxzpSzx;
    a31 = a13; a32 = a23; a33 = Syy-Sxx-Szz - mxEigenV; a34 = SyzpSzy;
    a41 = a14; a42 = a24; a43 = a34; a44 = Szz - SxxpSyy - mxEigenV;
    a3344_4334 = a33 * a44 - a43 * a34; a3244_4234 = a32 * a44-a42*a34;
    a3243_4233 = a32 * a43 - a42 * a33; a3143_4133 = a31 * a43-a41*a33;
    a3144_4134 = a31 * a44 - a41 * a34; a3142_4132 = a31 * a42-a41*a32;
    q1 =  a22*a3344_4334-a23*a3244_4234+a24*a3243_4233;
    q2 = -a21*a3344_4334+a23*a3144_4134-a24*a3143_4133;
    q3 =  a21*a3244_4234-a22*a3144_4134+a24*a3142_4132;
    q4 = -a21*a3243_4233+a22*a3143_4133-a23*a3142_4132;

    qsqr = q1 * q1 + q2 * q2 + q3 * q3 + q4 * q4;

/* The following code tries to calculate another column in the adjoint matrix when the norm of the 
   current column is too small.
   Usually this block will never be activated.  To be absolutely safe this should be
   uncommented, but it is most likely unnecessary.
*/
    if (qsqr < evecprec)
    {
        q1 =  a12*a3344_4334 - a13*a3244_4234 + a14*a3243_4233;
        q2 = -a11*a3344_4334 + a13*a3144_4134 - a14*a3143_4133;
        q3 =  a11*a3244_4234 - a12*a3144_4134 + a14*a3142_4132;
        q4 = -a11*a3243_4233 + a12*a3143_4133 - a13*a3142_4132;
        qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

        if (qsqr < evecprec)
        {
            double a1324_1423 = a13 * a24 - a14 * a23, a1224_1422 = a12 * a24 - a14 * a22;
            double a1223_1322 = a12 * a23 - a13 * a22, a1124_1421 = a11 * a24 - a14 * a21;
            double a1123_1321 = a11 * a23 - a13 * a21, a1122_1221 = a11 * a22 - a12 * a21;

            q1 =  a42 * a1324_1423 - a43 * a1224_1422 + a44 * a1223_1322;
            q2 = -a41 * a1324_1423 + a43 * a1124_1421 - a44 * a1123_1321;
            q3 =  a41 * a1224_1422 - a42 * a1124_1421 + a44 * a1122_1221;
            q4 = -a41 * a1223_1322 + a42 * a1123_1321 - a43 * a1122_1221;
            qsqr = q1*q1 + q2 *q2 + q3*q3+q4*q4;

            if (qsqr < evecprec)
            {
                q1 =  a32 * a1324_1423 - a33 * a1224_1422 + a34 * a1223_1322;
                q2 = -a31 * a1324_1423 + a33 * a1124_1421 - a34 * a1123_1321;
                q3 =  a31 * a1224_1422 - a32 * a1124_1421 + a34 * a1122_1221;
                q4 = -a31 * a1223_1322 + a32 * a1123_1321 - a33 * a1122_1221;
                qsqr = q1*q1 + q2 *q2 + q3*q3 + q4*q4;
                
                if (qsqr < evecprec)
                {
                    /* if qsqr is still too small, return the identity matrix. */
                    rot[0] = rot[4] = rot[8] = 1.0;
                    rot[1] = rot[2] = rot[3] = rot[5] = rot[6] = rot[7] = 0.0;

                    return(0);
                }
            }
        }
    }

    normq = sqrt(qsqr);
    q1 /= normq;
    q2 /= normq;
    q3 /= normq;
    q4 /= normq;

    a2 = q1 * q1;
    x2 = q2 * q2;
    y2 = q3 * q3;
    z2 = q4 * q4;

    xy = q2 * q3;
    az = q1 * q4;
    zx = q4 * q2;
    ay = q1 * q3;
    yz = q3 * q4;
    ax = q1 * q2;

    rot[0] = a2 + x2 - y2 - z2;
    rot[1] = 2 * (xy + az);
    rot[2] = 2 * (zx - ay);
    rot[3] = 2 * (xy - az);
    rot[4] = a2 - x2 + y2 - z2;
    rot[5] = 2 * (yz + ax);
    rot[6] = 2 * (zx + ay);
    rot[7] = 2 * (yz - ax);
    rot[8] = a2 - x2 - y2 + z2;

    return (1);
}



const int 
xx = 0, xy = 1, xz = 2,
yx = 3, yy = 4, yz = 5,
zx = 6, zy = 7, zz = 8;

/*!
 * Calculate the determinant of a 3x3 matrix
 * \param[in] m Pointer to the matrix.
 * \param[out] det The calculated determinant.
 */
void determinant3by3(data_size *m, data_size *det)
{
	*det = 0.0;
	*det += ( m[xx] * m[yy] * m[zz] );
	*det += ( m[xy] * m[yz] * m[zx] );
	*det += ( m[xz] * m[yx] * m[zy] );
	*det -= ( m[zx] * m[yy] * m[xz] );
	*det -= ( m[zy] * m[yz] * m[xx] );
	*det -= ( m[zz] * m[yx] * m[xy] );
}

/*!
 * Calculates the 3x3 rotation matrix between to structures using quaternions.
 * \param[in] target The structure to align.
 * \param[in] model The reference structure to align from.
 * \param[in] atoms The number of atoms in the structures.
 * \param[in] S The 3x3 covariance matrix between the target and model.
 * \param[out] rot The 3x3 rotation matrix.
 */
void quaternion(data_size *target, data_size *model,const int atoms,
		double *S,data_size *rot)
{
	double Syz_minus_zy = S[yz] - S[zy];
	double Szx_minus_xz = S[zx] - S[xz];
	double Sxy_minus_yx = S[xy] - S[yx];
	double Syz_plus_zy  = S[yz] + S[zy];
	double Sxz_plus_zx  = S[xz] + S[zx];
	double Sxy_plus_yx  = S[xy] + S[yx];

	double S2[9] = 
	{
		S[xx] * S[xx], S[xy] * S[xy], S[xz] * S[xz],
		S[yx] * S[yx], S[yy] * S[yy], S[yz] * S[yz],
		S[zx] * S[zx], S[zy] * S[zy], S[zz] * S[zz] 
	};

	double sigma = 0.0;
	size_t i;
	for(i = 0; i < atoms; i++)
	{
		sigma += target[i*3    ] * target[i*3    ] / atoms; 
		sigma += target[i*3 + 1] * target[i*3 + 1] / atoms; 
		sigma += target[i*3 + 2] * target[i*3 + 2] / atoms; 
		//sigma += (target[i*3] + target[i*3+1] + target[i*3+2]) *(target[i*3] + target[i*3+1] + target[i*3+2]); 
	}
	double inners = sigma;
	sigma = 0.0;
	for(i = 0; i < atoms; i++)
	{
		sigma += model[i*3    ] * model[i*3    ] / atoms; 
		sigma += model[i*3 + 1] * model[i*3 + 1] / atoms; 
		sigma += model[i*3 + 2] * model[i*3 + 2] / atoms; 
		//sigma += (model[i*3] + model[i*3+1] + model[i*3+2]) * (model[i*3] + model[i*3+1] + model[i*3+2]); 
	}
	inners += sigma;;
	//xxx+=sqrt(sigma);
	//sigma = xx/atoms;
	sigma = inners * 0.5;//S[0] + S[4] + S[8];
	//sigma /= 2.0;
	//sigma = 3.0;
    
    //inners = 2.0*((S[xx] + S[yy] + S[zz]) * std);
	data_size rmsd = 0.0;
	//printf("****RMSD IS %f ??? \n",sqrt( 
		//FastCalcRMSDAndRotation(rot, S, &rmsd , sigma, atoms, 1e-7);
	//rot[0] = rmsd;
    //return;
	double C[3] = {0.0,0.0,0.0};
	C[2] = 0.0;
	for(i = 0; i < 9; i++)
		C[2] += S2[i];

	C[2] = (-2.0) * C[2];

	C[1] = 8.0 * (   ( S[xx] * S[yz] * S[zy] )   + 
			         ( S[yy] * S[zx] * S[xz] )   + 
			         ( S[zz] * S[xy] * S[yx] ) ) - 
		    8.0 * (  ( S[xx] * S[yy] * S[zz] )   + 
		    		 ( S[yz] * S[zx] * S[xy] )   + 
		    		 ( S[zy] * S[yx] * S[xz] ) );

	double scratch = ( S2[xy] + S2[xz] - S2[yx] - S2[zx] );
	double D = scratch * scratch;
	scratch = S2[yy] - S2[xx] + S2[zz] + S2[yz] + S2[zy];

	double scratch2 = S[yy] * S[zz] - S[yz] * S[zy];
	double E = ( scratch - 2.0 * scratch2) * ( scratch + 2.0 * scratch2 );

	scratch =  ( ( S[xy] - S[yx] ) * ( S[xx] - S[yy] - S[zz] ) ) - 
		(( S[xz] + S[zx] ) * (S[yz] - S[zy] )) ;
	scratch2 =  ( ( S[xy] - S[yx] ) * ( S[xx] - S[yy] + S[zz] ) ) - 
		(( S[xz] + S[zx] ) * (S[yz] - S[zy] )) ;
	
	double F = scratch * scratch2;

	scratch =  ( -1.0 * (( S[xz] + S[zx] ) * ( S[yz] + S[zy] ) )) -                   
						 ( ( S[xy] + S[yx] ) * ( S[xx] + S[yy] - S[zz] ) );
	
	scratch2 = ( -1.0 * (( S[xz] - S[zx] ) * ( S[yz] - S[zy] )) ) -                   
						 ( ( S[xy] + S[yx] ) * ( S[xx] + S[yy] + S[zz] ) );
	double G = scratch * scratch2;

	scratch =  (      ( S[xy] + S[yx] ) * ( S[yz] + S[zy] ) ) +                   
						 ( ( S[xz] + S[zx] ) * ( S[xx] - S[yy] + S[zz] ) );
	
	scratch2 = ( -1.0 * (( S[xy] - S[yx] ) * ( S[yz] - S[zy] )) ) +                   
						 ( ( S[xz] + S[zx] ) * ( S[xx] + S[yy] + S[zz] ) );
	double H = scratch * scratch2;

	scratch =  ( 	   ( S[xy] + S[yx] ) * ( S[yz] - S[zy] ) ) +                   
						 ( ( S[xz] - S[zx] ) * ( S[xx] - S[yy] - S[zz] ) );
	
	scratch2 = ( -1.0 * (( S[xy] - S[yx] ) * ( S[yz] + S[zy] ) )) +                   
						 ( ( S[xz] - S[zx] ) * ( S[xx] + S[yy] - S[zz] ) );
	double I = scratch * scratch2;

	C[0] = D + E + F + G + H + I;
    
	//printf("WORKERRMS: C0 is %f\n",C[0]);
	//printf("WORKERRMS: C1 is %f\n",C[1]);
	//printf("WORKERRMS: C2 is %f\n",C[2]);
	//printf("WORKERRMS: sigma is %f\n",sigma);
	//printf("C is %f %f %f\n",C[0],C[1],C[2]);
	double P,dP,sigma2;
	const double eps = 1e-8;
	double delta = 0.0;
	double delta_old = 0.0;
	double sigma_old;
	int max_iter = 128;
	int conv = 0;
	int iter = 0;
	while ( iter++ < max_iter )
	{
		sigma_old = sigma;
		sigma2 = sigma * sigma;
		P = ( sigma2 * sigma2 ) +
			 ( C[2]   * sigma2 ) +
			 ( C[1]   * sigma  ) +
			 ( C[0] 				 );
		dP = ( 4.0 * ( sigma2 * sigma ) ) + ( 2.0 * C[2] * sigma ) + C[1];

		sigma = sigma - ( P / dP );
		delta = sigma - sigma_old;
		//delta = sigma - inners;
		//printf("i = %d SIGMA IS %8.6f D = %8.6f P = %8.6f dP = %8.6f P/dP = %8.6f\n",iter,2.0*sigma,delta,P,dP,P/dP);
		if ( fabs( delta ) < eps)
		{
			break;
		}
	}
	if ( iter == max_iter && conv == 0) 
	{
		printf("WORKERRMS: Maximum (%d) steps exceeded and d = %e. Not rotating.\n", max_iter,delta);
		return;
	}
    //sigma *= std;
	rot[0] = sqrt(fabs(inners - 2.0*sigma));
	//printf("RMSD SHOULD BE %f s = %f inner = %f diff = %f\n",
	//		sqrt(fabs(inners - (2.0*sigma))),
	//			2.0*sigma,inners,inners - (2.0*sigma));
	return;
#if 0
	double K[16] = 
	{
		S[xx] + S[yy] + S[zz], Syz_minus_zy, Szx_minus_xz, Sxy_minus_yx,
		Syz_minus_zy, S[xx] - S[yy] - S[zz], Sxy_plus_yx, Sxz_plus_zx,
		Szx_minus_xz, Sxy_plus_yx, S[yy] - S[zz] - S[xx], Syz_plus_zy,
		Sxy_minus_yx, Sxz_plus_zx, Syz_plus_zy, S[zz] - S[xx] - S[yy] 
	};
	data_size scratch3[9];
	data_size det[4] = { 0.0, 0.0, 0.0, 0.0 };

	scratch3[xx] = K[ 5] - sigma; scratch3[xy] = K[6 ]; 		 scratch3[xz] = K[7 ];
	scratch3[yx] = K[ 9]; 		  scratch3[yy] = K[10] - sigma;  scratch3[yz] = K[11];
	scratch3[zx] = K[13]; 		  scratch3[zy] = K[14]; 		 scratch3[zz] = K[15] - sigma;
	determinant3by3(scratch3, &det[0]);

	scratch3[xx] = K[ 4]; 
	scratch3[yx] = K[ 8]; 
	scratch3[zx] = K[12]; 
	determinant3by3(scratch3, &det[1]);

	scratch3[xy] = K[ 5] - sigma; 
	scratch3[yy] = K[ 9]; 
	scratch3[zy] = K[13]; 
	determinant3by3(scratch3, &det[2]);

	scratch3[xz] = K[ 6]; 
	scratch3[yz] = K[10] - sigma; 
	scratch3[zz] = K[14]; 
	determinant3by3(scratch3, &det[3]);

	det[1] = -1.0 * det[1];
	det[3] = -1.0 * det[3];

	data_size detsq =    det[0] * det[0] 
		            + det[1] * det[1] 
		            + det[2] * det[2] 
		            + det[3] * det[3];
	data_size detnrm = sqrt(detsq);
	det[0] /= detnrm;
	det[1] /= detnrm;
	det[2] /= detnrm;
	det[3] /= detnrm;
	//printf("l = %f q = %f %f %f %f\n",sigma, det[0],det[1],det[2],det[3]);
	//rot[xx] = det[0] * det[0] + det[1] * det[1] - det[2] * det[2] - det[3] * det[3];
	//rot[yy] = det[0] * det[0] - det[1] * det[1] + det[2] * det[2] - det[3] * det[3];
	//rot[zz] = det[0] * det[0] - det[1] * det[1] - det[2] * det[2] + det[3] * det[3];
	//
	//rot[xy] = 2 * ( det[1] * det[2] + det[0] * det[3]);
	//rot[yx] = 2 * ( det[1] * det[2] - det[0] * det[3]);
	//
	//rot[xz] = 2 * ( det[1] * det[3] - det[0] * det[2]);
	//rot[zx] = 2 * ( det[1] * det[3] + det[0] * det[2]);
	//
	//rot[yz] = 2 * ( det[2] * det[3] + det[0] * det[1]);
	//rot[zy] = 2 * ( det[2] * det[3] - det[0] * det[1]);
	
	rot[0] = det[0] * det[0] + det[1] * det[1] - det[2] * det[2] - det[3] * det[3];
	rot[1] = 2.0 * ( det[1] * det[2] - det[0] * det[3]);
	rot[2] = 2.0 * ( det[0] * det[2] + det[1] * det[3]);
	rot[3] = 2.0 * ( det[0] * det[3] + det[1] * det[2]);
	rot[4] = det[0] * det[0] - det[1] * det[1] + det[2] * det[2] - det[3] * det[3];
	rot[5] = 2.0 * ( det[2] * det[3] - det[0] * det[1]);
	rot[6] = 2.0 * ( det[1] * det[3] - det[0] * det[2]);
	rot[7] = 2.0 * ( det[0] * det[1] + det[2] * det[3]);
	rot[8] = det[0] * det[0] - det[1] * det[1] - det[2] * det[2] + det[3] * det[3];
	//data_size v = 0.0;
	//for(i = 0; i < 3; i++)
	//{
	//	v += S[i*3    ] * rot[0*3 + i]; 
	//	v += S[i*3 + 1] * rot[1*3 + i]; 
	//	v += S[i*3 + 2] * rot[2*3 + i]; 
	//	//sigma += (target[i*3] + target[i*3+1] + target[i*3+2]) *(target[i*3] + target[i*3+1] + target[i*3+2]); 
	//}
	printf("WORKERRMS: det quaternion is %f %f %f %f\n",
			det[0],det[1],det[2],det[3]);
#endif
}

static void worker_align_module_linux_simple_process(struct worker* self, struct module* module, 
		struct buffer** in, struct buffer** out)
{
	int i = 0;
	unsigned int index = 0;
	unsigned int atoms3 = 0;
	unsigned int stride3 = 0;
	unsigned int chunkSize = 0;
	size_t stride = in[0]->dim[0];
	size_t atoms = in[0]->dim[1];
	data_size* data = in[0]->data;
	data_size* mask = in[0]->data;
	size_t maskSize = atoms;
	atoms3 = atoms *3;
	chunkSize = stride *atoms3;
	data_size* offset = NULL;
	data_size* maskRef = in[3]->data;
	size_t counter;

#if STATS
	struct timespec start,finish;
	clock_gettime(CLOCK_REALTIME,&start);
#endif
	data_size* d;

	for(counter = 0; counter < stride; counter++)
	{

		d = ((data_size*)(out[0]->data)) + 9*counter;
		// calculate the 3x3 covariance matrix (correlation matrix)
		data_size *rotation,rot[9] = { 1.0,0.0,0.0,0.0,1.0,0.0,0.0,0.0,1.0};
		double cov[9] = {0.0,0.0,0.0,  0.0,0.0,0.0,  0.0,0.0,0.0};
		
		double buf=0.0;
		size_t k =0,j= 0;
        data_size tstd,rstd;
#if CORR
        std = 0.0;
#endif
		for(i = 0; i < 3; i++)
		{
#if CORR
			for(k=0,rmean=0.0; k < maskSize; k++)
				rmean += maskRef[k*3+i] / maskSize;
			for(k=0,rstd=buf=0.0; k < maskSize; k++)
			{
				buf = (maskRef[k*3+i] - rmean);
				rstd += (buf*buf) / maskSize;
			}
			rstd = sqrt(rstd);
#endif
			for(j = 0; j < 3; j++)
			{
				//printf("SUM START\n");
#if CORR
				for(k=0,tmean=0.0; k < maskSize; k++)
				{
					tmean += mask   [k*3+j] / maskSize;
					//printf("SUM IS %f / %zu = %f \n",mask[k*3+j],maskSize,tmean);
				}
				for(k=0,tstd=buf=0.0; k < maskSize; k++)
				{
					buf = (mask[k*3+j] - tmean);
					tstd += (buf*buf) / maskSize;
				}
				tstd = sqrt(tstd);
#endif
				for(k=0,buf=0.0; k < maskSize; k++)
                {
#if CORR
                    buf +=  ((maskRef[k*3+i] - rmean ) * (mask[k*3+j] - tmean) ) / maskSize;
#else
					buf +=  ((maskRef[k*3+i]) * (mask[k*3+j]) ) / maskSize;
#endif
                }
#if CORR
				cov[i*3+j] = buf / (rstd * tstd);
                if(j == i) std += (rstd*tstd) / 3.0;
#else
				//cov[i*3+j] = cov[j*3+i] = buf ;
				cov[i*3+j] = buf ;
#endif
				//cov[i*3+j] =  (tmean) ;//* tstd);
			}
		}
		
		quaternion(mask,maskRef,maskSize,cov,rot);
		rotation = &rot[0];
		//memcpy(out[0]->data + (counter*9*sizeof(data_size)),rot,out[0]->bytes_per_elem*out[0]->dim[1]);
		//memcpy(out[0]->data + (counter*sizeof(data_size)),rot,out[0]->bytes_per_elem);
		memcpy(out[0]->data + (out[0]->dim[1]*counter*sizeof(data_size)),rot,out[0]->bytes_per_elem*out[0]->dim[1]);
		
		//printf("RMS: ROT is:\t");
		//	for(i=0;i<9;i++) printf(" %4.7f",cov[i]);
		//printf("\n");
		//data_size result[3];
		//data_size *index = data;
		//rotation
		//for(i = buf = 0; i < atoms; i+=3,buf=0.0)
		//{
		//	buf += data[i + 0] * rotation[0];
		//	buf += data[i + 1] * rotation[3];
		//	buf += data[i + 2] * rotation[6];
		//	data[i + 0] = buf;
		//	buf = 0.0;
		//	buf += data[i + 0] * rotation[1];
		//	buf += data[i + 1] * rotation[4];
		//	buf += data[i + 2] * rotation[7];
		//	data[i + 1] = buf;
		//	buf = 0.0;
		//	buf += data[i + 0] * rotation[2];
		//	buf += data[i + 1] * rotation[5];
		//	buf += data[i + 2] * rotation[8];
		//	data[i + 2] = buf;
		//}
	//	int sync = false;
	//	if ( sync )
	//	{
	//		for(i = buf = 0; i < maskSize*3; i+=3,buf=0.0)
	//		{
	//			buf += mask[i] * rotation[0];
	//			buf += mask[i] * rotation[3];
	//			buf += mask[i] * rotation[6];
	//			mask[i] = buf;
	//			buf = 0.0;
	//			buf += mask[i + 1] * rotation[1];
	//			buf += mask[i + 1] * rotation[4];
	//			buf += mask[i + 1] * rotation[7];
	//			mask[i + 1] = buf;
	//			buf = 0.0;
	//			buf += mask[i + 2] * rotation[2];
	//			buf += mask[i + 2] * rotation[5];
	//			buf += mask[i + 2] * rotation[8];
	//			mask[i + 2] = buf;
	//		}
	//	}

		data += atoms 	  * 3;
		mask += maskSize * 3;
	}
#if STATS
	clock_gettime(CLOCK_REALTIME, &finish);
			self->exe_tot += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
	//double mib = (in[0]->size)	/ 1e6;
	//double sec = (finish.tv_sec - start.tv_sec) + ((finish.tv_nsec - start.tv_nsec) / 1e9);
#ifdef ALIGN_BANDWIDTH
	//printf("ALIGN  Kernel Bandwidth is %f MiB at %f = %f MiB/s\n", mib, sec , mib/sec);
#endif

}
static void worker_align_module_cl_process(struct worker* self, struct module* module, 
		struct buffer** in, struct buffer** out)
{
	struct worker_align_properties* props = 
		(struct worker_align_properties*)self->properties;
	struct device_cl_properties* device_props = 
		(struct device_cl_properties*)module->driver[0].device->properties;
	struct driver_cl_properties* driver_props = 
		(struct driver_cl_properties*)module->driver[0].properties;
	
	size_t l = ALIGN_BLK;
	size_t dims[] =  { l , in[0]->dim[0] };
	size_t ldims[] = { l , 1};
	size_t dim[] = { in[0]->dim[0], in[0]->dim[1], in[0]->dim[1] };
	cl_uint num_dims = 2;
	cl_bool blk = CL_ASYNC ? CL_FALSE : CL_TRUE;
	static int first = 0;
#if STATS
	struct timespec start,finish;
	unsigned long long time;
	clock_gettime(CLOCK_REALTIME,&start);
#endif
	
	cl_int ret = 0;
	if(first == 0)
	{
		ret = clEnqueueWriteBuffer(device_props->in,
			in[2]->data,
			blk,
			0,
			sizeof(size_t) * 3,
			dim,
			//0, NULL, NULL);
			0, NULL, (cl_event*)&in[2]->lock);
		//clFlush(device_props->in);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("ALGN: CPY FAIL %d\n",ret);
	}
	ret = clSetKernelArg(driver_props->kernel,
			0,
			sizeof(cl_mem),
			&in[2]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			3,
			sizeof(cl_mem),
			&in[3]->data);
	}
	ret |= clSetKernelArg(driver_props->kernel,
			4,
			sizeof(cl_mem),
			&out[0]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			1,
			sizeof(cl_mem),
			&in[0]->data);
	ret |= clSetKernelArg(driver_props->kernel,
			2,
			sizeof(cl_mem),
			&in[0]->data);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("ARG FAIL %d\n",ret);
	}
	if(!first)
	clWaitForEvents(1, (cl_event*)&in[2]->lock);
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
			self->memin_tot += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
	clock_gettime(CLOCK_REALTIME,&start);
	//device_cl_print_kernel_info(driver_props->kernel,device_props->id);
	register_profile_event(self,"proc_ocl_kern_wait_start",self->frame_cnt,0);
#endif
	ret = clEnqueueNDRangeKernel(device_props->exe,
			driver_props->kernel,
			num_dims,
			NULL,
			dims,
			ldims,
			0,
			NULL,
			//NULL);
			(cl_event*)&in[0]->lock);
		//clFlush(device_props->exe);
	if(unlikely(ret != CL_SUCCESS))
	{
		printf("EXEC ERROR: %d\n",ret);
	}
	clWaitForEvents(1,(cl_event*)&in[0]->lock);
#if STATS
		clock_gettime(CLOCK_REALTIME, &finish);
			self->exe_tot += (finish.tv_sec - start.tv_sec) + 
				((finish.tv_nsec - start.tv_nsec) / 1e9);
	if(!first)
	{
	self->memin_exe += device_cl_get_exe_time_seconds(in[2]->lock);
	self->memin_wait += device_cl_get_queued_time_seconds(in[2]->lock);
	self->memin_delay += device_cl_get_delay_time_seconds(in[2]->lock);
	first = 1;
	}
	self->exe_exe += device_cl_get_exe_time_seconds(in[0]->lock);
	self->exe_wait += device_cl_get_queued_time_seconds(in[0]->lock);
	self->exe_delay += device_cl_get_delay_time_seconds(in[0]->lock);
		time = (start.tv_sec * 1e9) + (start.tv_nsec);
		time +=  (unsigned long long)
			(device_cl_get_queued_time_seconds(in[0]->lock)*1e9);
		register_profile_event(self,"proc_ocl_kern_delay_start",self->frame_cnt,&time);
		time +=  (unsigned long long)
			(device_cl_get_delay_time_seconds(in[0]->lock)*1e9);
		register_profile_event(self,"proc_ocl_kern_exe_start",self->frame_cnt,&time);
		time +=  (unsigned long long)
			(device_cl_get_exe_time_seconds(in[0]->lock)*1e9);
		//time = (finish.tv_sec * 1e9) + (finish.tv_nsec);
		register_profile_event(self,"proc_fin",self->frame_cnt,&time);
	clock_gettime(CLOCK_REALTIME, &start);
	clock_gettime(CLOCK_REALTIME, &finish);
		self->proc_fin += (finish.tv_sec - start.tv_sec) + 
			((finish.tv_nsec - start.tv_nsec) / 1e9);
#endif
#if CLEAN_UP
	if(in[2]->lock) clReleaseEvent(in[2]->lock) ;in[2]->lock = 0;
	if(in[0]->lock) clReleaseEvent(in[0]->lock) ;in[0]->lock = 0;
#endif
}

static int worker_align_module_linux_simple_bld(struct worker* self, 
		struct module* module, struct system* system) 
{ 
	if(! module_linux_simple_bld(module,system) )
	{
		worker_add_module(self,module,&worker_align_module_linux_simple_process);
		return  0;
	}
	else
		return -1;
}

static int worker_align_module_linux_simple_dty(struct worker* self,struct module* module) 
{ 
	return 0; 
}

static int worker_align_module_cl_bld(struct worker* self, struct module* module,
		struct system* system) 
{ 
	char* filename = "align_shm.cl";
	char* kname    = "align";
	if(! module_cl_bld(module,system,filename,kname) )
	{
		worker_add_module(self,module,&worker_align_module_cl_process);
		return 0;
	}
	else
		return -1;
}

static int worker_align_module_cl_dty(struct worker* self,struct module* module) 
{ 
	return 0; 
}

static void worker_align_prepare_buffers(struct worker* self, struct module* module, 
		struct frame** f)
{
	struct frame* frame = f[0];
	size_t dims[] = { 3 };
	size_t num_dims = 1;
	struct device* device = module->driver[0].device;
	if(unlikely (self->input == NULL) ) 
	{
		self->input = malloc(sizeof(struct buffer*) * 5);
		if(likely   (self->input != NULL) )
		{
			self->input[0] = NULL;
			self->input[1] = NULL;
			self->input[2] = NULL;
			self->input[3] = NULL;
			self->input[4] = NULL;
		}
	}

	if(likely   (self->input != NULL) )
	{
		self->input[0] = &frame->buffer[0];
		self->input[1] = &frame->buffer[0];
		if(unlikely(self->input[2] == NULL))
		{
			self->input[2] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[2], 
					'r' + 'p',sizeof(size_t), dims,num_dims);
		}
		if(unlikely(self->input[3] == NULL))
		{
			size_t dims[] = { self->input[0]->dim[1], self->input[0]->dim[2] };
			size_t num_dims = 2;
			self->input[3] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->input[3], 
					'r',sizeof(data_size), dims,num_dims);
			data_size* tdata = device->fn->map(device,self->input[0]);
			data_size* sdata = device->fn->map(device,self->input[3]);
			memcpy(sdata,tdata,self->input[3]->size);
			device->fn->unmap(device,self->input[0]);
			device->fn->unmap(device,self->input[3]);
		}
		self->input[4] = NULL;
	}
	if(unlikely (self->output == NULL) ) 
	{
		self->output = malloc(sizeof(struct buffer*) * 2);
		self->output[0] = NULL;
		self->output[1] = NULL;
	}

	if(likely   (self->output != NULL) )
	{
	//	const int num_dims = 2;
	//	size_t d[] = { 3, frame->buffer[1].dim[0] };
	//	self->output[0] = (struct buffer*)self->host->fn->mem_alloc(self->host,sizeof(struct buffer));
	//	struct device* device = module->driver[module->driver_next_idx].device;
	//	device->fn->buffer_bld(device,self->host,self->output[0], 'w',sizeof(double),
	//			d,num_dims);
		if(unlikely(self->output[0] == NULL))
		{
			size_t dims[] = { self->input[0]->dim[0], 1 };
			size_t num_dims = 2;
			self->output[0] = (struct buffer*)self->host->fn->mem_alloc(self->host,
					sizeof(struct buffer) );
			device->fn->buffer_bld(device,self->host,self->output[0], 
					'w' + 'p',sizeof(data_size), dims,num_dims);
		}
		self->output[1] = NULL;
	}
}
static void worker_align_process_buffers(struct worker* self)
{
	FILE* fout;
	size_t i,j;
	struct device_cl_properties* device_props = self->input[0]->device->properties;
#if 0
	data_size* d;
	fout = fopen("align.dat","a");
	if(fout)
	{
		//data_size*d = malloc(self->output[0]->size / self->output[0]->dim[0]);
		//clEnqueueReadBuffer(device_props->out,
		//		self->output[0]->data,
		//		CL_TRUE,
		//		0,
		//		self->output[0]->size / self->output[0]->dim[0],
		//		d,
		//		0, NULL, NULL);
		//free(d);
	//((struct buffer_cl_properties*)self->input[0]->properties)->map_flags = CL_MAP_READ;
	data_size* d = self->input[0]->device->fn->map(self->input[0]->device,
			self->input[0]);
	for(i = 0; i < self->input[0]->dim[0] * self->input[0]->dim[1] * 3; i+=3)
		//fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%zu: %f %f %f\n",i,d[i],d[i+1],d[i+2]);
	fprintf(fout,"**\n");
	self->input[0]->device->fn->unmap(self->input[0]->device,self->input[0]);
	fclose(fout);
	}
#endif
#if 1
	fout = ((struct worker_align_properties*)self->properties)->frmsd;
	//fout = fopen("qrmsd.out","a");
	if(likely(fout != NULL))
	{
	//	data_size*d = malloc(self->output[0]->size);
	//	clEnqueueReadBuffer(device_props->out,
//				self->output[0]->data,
//				CL_TRUE,
//				0,
//				self->output[0]->size,
//				d,
//				0, NULL, NULL);
#if USE_CL
	//((struct buffer_cl_properties*)self->output[0]->properties)->map_flags = CL_MAP_READ;
	//((struct buffer_cl_properties*)self->output[0]->properties)->async = CL_TRUE;
#endif
	data_size* d = self->output[0]->device->fn->map(self->output[0]->device,
			self->output[0]);
	for(i = 0; i < self->output[0]->dim[0] ; i++)
    {
       //fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
	    for(j = 0; j < self->output[0]->dim[1] ; j++)
		    fprintf(fout,"%8.4f ",d++[0]);
		fprintf(fout,"\n");
    }
	self->output[0]->device->fn->unmap(self->output[0]->device,self->output[0]);
	//	free(d);
	//fclose(fout);
	}
#endif
#if 0
	fout = fopen("rot.dat","a");
	if(fout)
	{
	//	data_size*d = malloc(self->output[0]->size);
	//	clEnqueueReadBuffer(device_props->out,
//				self->output[0]->data,
//				CL_TRUE,
//				0,
//				self->output[0]->size,
//				d,
//				0, NULL, NULL);
	//((struct buffer_cl_properties*)self->output[0]->properties)->map_flags = CL_MAP_READ;
	data_size* d = self->output[0]->device->fn->map(self->output[0]->device,
			self->output[0]);
	for(i = 0; i < self->output[0]->dim[0] * 9 ; i+=9)
		//fprintf(fout,"%p %p %p\n",&d[i],&d[i+1],&d[i+2]);
		//fprintf(fout,"%x %x %x\n",(uint)d[i],(uint)d[i+1],(uint)d[i+2]);
		fprintf(fout,"%zu: %f %f %f %f %f %f %f %f %f\n",i/9,
				d[i+0],d[i+1],d[i+2],
				d[i+3],d[i+4],d[i+5],
				d[i+6],d[i+7],d[i+8]);
	fprintf(fout,"**\n");
	self->output[0]->device->fn->unmap(self->output[0]->device,self->output[0]);
	//	free(d);
	fclose(fout);
	}
#endif

}
