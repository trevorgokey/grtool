#pragma OPENCL EXTENSION cl_khr_fp64 : enable

#define xx  0
#define xy  1
#define xz  2
#define yx  3
#define yy  4
#define yz  5
#define zx  6
#define zy  7
#define zz  8

#ifdef USE_DOUBLE
typedef double data_size;
#else
typedef float data_size;
#endif

#define EPS 1e-11
#define LOOSE 1e-7


#define SHR_SIZE (ALIGN_SHR_SIZE*3)

void reduce_three(__private double* sum, __local double* shm, __private size_t *a)
{
	size_t x = get_local_id(0);
	size_t x_len = get_local_size(0);
	size_t blk = get_global_id(1);
	const uint r = (x_len % (SHR_SIZE/3)) > 0;
	size_t i,j;
	const size_t shr_idx = x % (SHR_SIZE/3);
	
	char write = 0;
	
	barrier(CLK_LOCAL_MEM_FENCE);
	
        //zero out shr mem
	for(i = x; i < SHR_SIZE; i+=x_len)
		shm[i] = 0.0;

	
	//condense into shr mem
	for(j = 1; j <= (x_len/(SHR_SIZE/3)) + r ;  j++ )
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(!write && x < j*(SHR_SIZE/3) && (write = 1))
		{
			shm[3*(shr_idx)    ] += sum[0];
			shm[3*(shr_idx) + 1] += sum[1];
			shm[3*(shr_idx) + 2] += sum[2];
		}
	}
	//reduce the shr mem

	for(i = SHR_SIZE/3 >> 1 ; i > 0; i >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(x < i)
		{
			shm[3*x    ] += shm[3*(x + i)    ];
			shm[3*x + 1] += shm[3*(x + i) + 1];
			shm[3*x + 2] += shm[3*(x + i) + 2];
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	sum[0] = shm[0];
	sum[1] = shm[1];
	sum[2] = shm[2];
}

void reduce(__private double* sum, __local double* shm)
{
	size_t x = get_local_id(0);
	size_t x_len = get_local_size(0);
	size_t blk = get_global_id(1);
	const size_t ss = SHR_SIZE /3 ;
	const uint r = (x_len % (ss)) > 0;
	size_t i,j;
	const size_t shr_idx = x % (ss);
	
	char write = 0;
	
	//printf("MY ID IS %zu MY LEN IS %zu\n",x,x_len);
	barrier(CLK_LOCAL_MEM_FENCE);
	//zero out shr mem
	for(i = x; i < ss; i+=x_len)
		shm[i] = 0.0;

	
	//condense into shr mem
	for(j = 1; j <= (x_len/(ss)) + r ;  j++ )
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(!write && x < j*(ss) && (write = 1))
		{
			shm[(shr_idx)    ] += sum[0];
		}
	}
	//reduce the shr mem
	for(i = ss >> 1 ; i > 0; i >>= 1)
	{
		barrier(CLK_LOCAL_MEM_FENCE);
		if(x < i)
		{
			shm[x    ] += shm[(x + i)    ];
		}
	}
	barrier(CLK_LOCAL_MEM_FENCE);
	sum[0] = shm[0];
}


inline data_size determinant3by3(__private data_size *m)
{
	data_size det = 0.0;
	det += ( m[xx] * m[yy] * m[zz] );
	det += ( m[xy] * m[yz] * m[zx] );
	det += ( m[xz] * m[yx] * m[zy] );
	det -= ( m[zx] * m[yy] * m[xz] );
	det -= ( m[zy] * m[yz] * m[xx] );
	det -= ( m[zz] * m[yx] * m[xy] );
	return det;
}


void calc_cp_coeff(__private double* C, __private double* S)
{
	C[2] = 0.0;
	for(int i = 0; i < 9; i++)
		C[2] += S[i] * S[i];
	C[2] = -2.0 * C[2];
	C[1] = 8.0 * ( ( S[xx] * S[yz] * S[zy] )   + 
                   ( S[yy] * S[zx] * S[xz] )   + 
				   ( S[zz] * S[xy] * S[yx] ) ) - 
		   8.0 * ( ( S[xx] * S[yy] * S[zz] )   + 
				   ( S[yz] * S[zx] * S[xy] )   + 
				   ( S[zy] * S[yx] * S[xz] ) );

	double scratch = ( S[xy]*S[xy] + S[xz]*S[xz] - S[yx]*S[yx] - S[zx]*S[zx] );
	double D = scratch * scratch;
	scratch = S[yy]*S[yy] - S[xx]*S[xx] + S[zz]*S[zz] + S[yz]*S[yz] + S[zy]*S[zy];
	
    double scratch2 = S[yy] * S[zz] - S[yz] * S[zy];
	double E = ( scratch - 2.0 * scratch2) * ( scratch + 2.0 * scratch2 );

	scratch =  ( -1.0 * ( S[xz] + S[zx] ) * ( S[yz] - S[zy] ) ) +
				       ( ( S[xy] - S[yx] ) * ( S[xx] - S[yy] - S[zz] ) );
	scratch2 = ( -1.0 * ( S[xz] - S[zx] ) * ( S[yz] + S[zy] ) ) +
				       ( ( S[xy] - S[yx] ) * ( S[xx] - S[yy] + S[zz] ) );
	double F = scratch * scratch2;

	scratch =  ( -1.0 * ( S[xz] + S[zx] ) * ( S[yz] + S[zy] ) ) -                   
				       ( ( S[xy] + S[yx] ) * ( S[xx] + S[yy] - S[zz] ) );
	scratch2 = ( -1.0 * ( S[xz] - S[zx] ) * ( S[yz] - S[zy] ) ) -                   
				       ( ( S[xy] + S[yx] ) * ( S[xx] + S[yy] + S[zz] ) );
	double G = scratch * scratch2;

	scratch =  (      ( S[xy] + S[yx] ) * ( S[yz] + S[zy] ) ) +                   
				 	    ( ( S[xz] + S[zx] ) * ( S[xx] - S[yy] + S[zz] ) );
	scratch2 = ( -1.0 * ( S[xy] - S[yx] ) * ( S[yz] - S[zy] ) ) +                   
				       ( ( S[xz] + S[zx] ) * ( S[xx] + S[yy] + S[zz] ) );
	double H = scratch * scratch2;
	
	scratch = ( 	   ( S[xy] + S[yx] ) * ( S[yz] - S[zy] ) ) +                   
				 		 ( ( S[xz] - S[zx] ) * ( S[xx] - S[yy] - S[zz] ) );
	scratch2 = ( -1.0 * ( S[xy] - S[yx] ) * ( S[yz] + S[zy] ) ) +                   
				  	    ( ( S[xz] - S[zx] ) * ( S[xx] + S[yy] - S[zz] ) );
	double I = scratch * scratch2;

	C[0] = D + E + F + G + H + I;
	//if(!get_local_id(0))	
	//printf("BLK %zu C is %f %f %f\n",get_group_id(0),C[0],C[1],C[2]);
}

int nr_minimize(__private double* C, __private double* s, __private int max_iter)
{
	double P,dP,sigma2,sigma = *s;
	double delta = 0.0;
	double sigma_old;
	int iter = 0;
	while ( ++iter < (max_iter) )
	{
		//if(!get_local_id(0))
		sigma_old = sigma;
		sigma2 = sigma * sigma;
		P = ( sigma2 * sigma2 ) +
		    ( C[2]   * sigma2 ) +
		    ( C[1]   * sigma  ) +
		    ( C[0] 				 );
		//printf("RMSD P %f\n",C[1]);
		dP = ( 4.0 * ( sigma2 * sigma ) ) + ( 2.0 * C[2] * sigma ) + C[1];
		
		sigma = sigma - ( P / dP );
		delta = sigma - sigma_old;
		//if(!get_global_id(0))
	//	printf("WORKERRMS: sigma: %.16f delta: %18.16e (%d iter), P = %.16f dP = %.16f P/dP = %18.16f\n", 2.0*sigma, delta, iter,P,dP,P/dP);
		
		if ( fabs( delta ) < EPS)
		{
		//	printf("WORKERRMS: converged to %lf (%d iter)\n",sigma,iter);
			*s = sigma;
			return 1;

		}
	}
	return 0;

}
void calc_quaternion(__private data_size* det, __private double* S, __private double sigma)
{
#define    Syz_minus_zy  (S[yz] - S[zy] )
#define    Szx_minus_xz  (S[zx] - S[xz] )
#define    Sxy_minus_yx  (S[xy] - S[yx] )
#define    Syz_plus_zy   (S[yz] +  S[zy])
#define    Sxz_plus_zx   (S[xz] + S[zx] )
#define    Sxy_plus_yx   (S[xy] + S[yx] )
	//data_size K[16] = {
	//S[xx] + S[yy] + S[zz],
	//Syz_minus_zy,
	//Szx_minus_xz, 
	//Sxy_minus_yx,
	//Syz_minus_zy,
	//S[xx] - S[yy] - S[zz], 
	//Sxy_plus_yx,
	//Sxz_plus_zx,
	//Szx_minus_xz, 
	//Sxy_plus_yx,
	//S[yy] - S[zz] - S[xx], 
	//Syz_plus_zy,
	//Sxy_minus_yx,
	//Sxz_plus_zx,
	//Syz_plus_zy,
	//S[zz] - S[xx] - S[yy] };
	
#define K_0   	(S[xx] + S[yy] + S[zz])
#define K_1   	(Syz_minus_zy)
#define K_2   	(Szx_minus_xz) 
#define K_3   	(Sxy_minus_yx)
#define K_4   	(Syz_minus_zy)
#define K_5   	(S[xx] - S[yy] - S[zz]) 
#define K_6   	(Sxy_plus_yx)
#define K_7   	(Sxz_plus_zx)
#define K_8   	(Szx_minus_xz) 
#define K_9   	(Sxy_plus_yx)
#define K_10  	(S[yy] - S[zz] - S[xx]) 
#define K_11  	(Syz_plus_zy)
#define K_12  	(Sxy_minus_yx)
#define K_13  	(Sxz_plus_zx)
#define K_14  	(Syz_plus_zy)
#define K_15  	(S[zz] - S[xx] - S[yy])
	

	double scratch3[9] = {0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0};
	
	//scratch3[xx] = K[ 5] - sigma; scratch3[xy] = K[6 ]; scratch3[xz] = K[7 ];
	//scratch3[yx] = K[ 9]; scratch3[yy] = K[10] - sigma; scratch3[yz] = K[11];
	//scratch3[zx] = K[13]; scratch3[zy] = K[14]; scratch3[zz] = K[15] - sigma;
	//det[0] = determinant3by3(scratch3);

	//scratch3[xx] = K[ 4]; 
	//scratch3[yx] = K[ 8]; 
	//scratch3[zx] = K[12]; 
	//det[1] = determinant3by3(scratch3);

	//scratch3[xy] = K[ 5] - sigma; 
	//scratch3[yy] = K[ 9]; 
	//scratch3[zy] = K[13]; 
	//det[2] = determinant3by3(scratch3);

	//scratch3[xz] = K[ 6]; 
	//scratch3[yz] = K[10] - sigma; 
	//scratch3[zz] = K[14]; 
	//det[3] = determinant3by3(scratch3);
	
	
	scratch3[xx] = K_5 - sigma; scratch3[xy] = K_6; scratch3[xz] = K_7;
	scratch3[yx] = K_9; scratch3[yy] = K_10 - sigma; scratch3[yz] = K_11;
	scratch3[zx] = K_13; scratch3[zy] = K_14; scratch3[zz] = K_15 - sigma;
	det[0] = determinant3by3(scratch3);

	scratch3[xx] = K_4; 
	scratch3[yx] = K_8; 
	scratch3[zx] = K_12; 
	det[1] = -1.0 * determinant3by3(scratch3);

	scratch3[xy] = K_5 - sigma; 
	scratch3[yy] = K_9; 
	scratch3[zy] = K_13; 
	det[2] = determinant3by3(scratch3);

	scratch3[xz] = K_6; 
	scratch3[yz] = K_10 - sigma; 
	scratch3[zz] = K_14; 
	det[3] = -1.0 * determinant3by3(scratch3);
	
	data_size detsq = det[0]*det[0] + det[1]*det[1] + det[2]*det[2] + det[3]*det[3];
	data_size detnrm = sqrt(detsq);
	det[0] /= detnrm;
	det[1] /= detnrm;
	det[2] /= detnrm;
	det[3] /= detnrm;

}

void quaternion_reduce_reg(__global data_size *target, 
				    __global data_size *model,
				    __private size_t *atoms,
				    __private double *S,
				    __private data_size *rot,
					__local double* shm)
{
		 

//	data_size S2[9] = {
//	S[xx] * S[xx], S[xy] * S[xy], S[xz] * S[xz],
//	S[yx] * S[yx], S[yy] * S[yy], S[yz] * S[yz],
//	S[zx] * S[zx], S[zy] * S[zy], S[zz] * S[zz] };

	double sigma = 0.0;
	double C[5];
	size_t len = *atoms * 3;
	double s = 0.0;
	for(int i = get_local_id(0); i< len; i += get_local_size(0))
		s += (double)(target[i] * target[i]) / (*atoms);
	reduce(&s,shm);
	sigma = s;
	s = 0.0;
	for(int i = get_local_id(0); i < len; i += get_local_size(0))
		s += (double)(model[i] * model[i]) / (*atoms);
	reduce(&s,shm);
	s = (s + sigma) ;
	sigma = 0.5 * s;
	
	calc_cp_coeff(C,S);
	if( nr_minimize(C,&sigma,256) == 0)
		return;
	
	rot[0] = sqrt(fabs(s - 2.0*sigma));
    return;
	//if(!get_local_id(0))
	//printf("RMSD is s = %.8f i = %.8f d = %.8f r = %.8f\n",2.0*sigma, s,s - 2.0*sigma,rot[0]); 
	data_size det[4] = { 0.0, 0.0, 0.0, 0.0 };
	calc_quaternion(det,S,sigma);
	
	rot[0] = det[0] * det[0] + det[1] * det[1] - det[2] * det[2] - det[3] * det[3];
	rot[1] = 2 * ( det[1] * det[2] - det[0] * det[3]);
	rot[2] = 2 * ( det[0] * det[2] + det[1] * det[3]);
	rot[3] = 2 * ( det[0] * det[3] + det[1] * det[2]);
	rot[4] = det[0] * det[0] - det[1] * det[1] + det[2] * det[2] - det[3] * det[3];
	rot[5] = 2 * ( det[2] * det[3] - det[0] * det[1]);
	rot[6] = 2 * ( det[1] * det[3] - det[0] * det[2]);
	rot[7] = 2 * ( det[0] * det[1] + det[2] * det[3]);
	rot[8] = det[0] * det[0] - det[1] * det[1] - det[2] * det[2] + det[3] * det[3];
}

void quaternion(__global data_size *target, 
				    __global data_size *model,
				    __private uint *atoms,
				    __private data_size *S,
				    __private data_size *rot,
					__local data_size* shm)
{
		 
#define    Syz_minus_zy  (S[yz] - S[zy] )
#define    Szx_minus_xz  (S[zx] - S[xz] )
#define    Sxy_minus_yx  (S[xy] - S[yx] )
#define    Syz_plus_zy   (S[yz] +  S[zy])
#define    Sxz_plus_zx   (S[xz] + S[zx] )
#define    Sxy_plus_yx   (S[xy] + S[yx] )
	data_size K[16] = {
	S[xx] + S[yy] + S[zz], Syz_minus_zy, Szx_minus_xz, Sxy_minus_yx,
	Syz_minus_zy, S[xx] - S[yy] - S[zz], Sxy_plus_yx, Sxz_plus_zx,
	Szx_minus_xz, Sxy_plus_yx, S[yy] - S[zz] - S[xx], Syz_plus_zy,
	Sxy_minus_yx, Sxz_plus_zx, Syz_plus_zy, S[zz] - S[xx] - S[yy] };

	data_size S2[9] = {
	S[xx] * S[xx], S[xy] * S[xy], S[xz] * S[xz],
	S[yx] * S[yx], S[yy] * S[yy], S[yz] * S[yz],
	S[zx] * S[zx], S[zy] * S[zy], S[zz] * S[zz] };

	double sigma = 0.0;
	for(int i = 0; i < *atoms; i++)
	{
		sigma += target[i*3    ] * target[i*3    ]; 
		sigma += target[i*3 + 1] * target[i*3 + 1]; 
		sigma += target[i*3 + 2] * target[i*3 + 2]; 
	}
	for(int i = 0; i < *atoms; i++)
	{
		sigma += model[i*3    ] * model[i*3    ]; 
		sigma += model[i*3 + 1] * model[i*3 + 1]; 
		sigma += model[i*3 + 2] * model[i*3 + 2]; 
	}
	double inners = sigma;
	sigma *= 0.5;
	
	data_size C[3];
	C[2] = 0.0;
	for(int i = 0; i < 9; i++)
		C[2] += S2[i];
	C[2] = -2.0 * C[2];
	C[1] = 8.0 * ( ( S[xx] * S[yz] * S[zy] ) + 
					 	( S[yy] * S[zx] * S[xz] ) + 
					 	( S[zz] * S[xy] * S[yx] ) ) - 
			 8.0 * ( ( S[xx] * S[yy] * S[zz] ) + 
					 	( S[yz] * S[zx] * S[xy] ) + 
					 	( S[zy] * S[yx] * S[xz] ) );

	data_size scratch = ( S2[xy] + S2[xz] - S2[yx] - S2[zx] );
	data_size D = scratch * scratch;
	
	scratch = S2[yy] - S2[xx] + S2[zz] + S2[yz] + S2[zy];
	data_size scratch2 = S[yy] * S[zz] - S[yz] * S[zy];
	data_size E = ( scratch - 2.0 * scratch2) * ( scratch + 2.0 * scratch2 );

	scratch =  ( -1.0 * ( S[xz] + S[zx] ) * ( S[yz] - S[zy] ) ) +
				       ( ( S[xy] - S[yx] ) * ( S[xx] - S[yy] - S[zz] ) );
	scratch2 = ( -1.0 * ( S[xz] - S[zx] ) * ( S[yz] + S[zy] ) ) +
				       ( ( S[xy] - S[yx] ) * ( S[xx] - S[yy] + S[zz] ) );
	data_size F = scratch * scratch2;

	scratch =  ( -1.0 * ( S[xz] + S[zx] ) * ( S[yz] + S[zy] ) ) -                   
				       ( ( S[xy] + S[yx] ) * ( S[xx] + S[yy] - S[zz] ) );
	scratch2 = ( -1.0 * ( S[xz] - S[zx] ) * ( S[yz] - S[zy] ) ) -                   
				       ( ( S[xy] + S[yx] ) * ( S[xx] + S[yy] + S[zz] ) );
	data_size G = scratch * scratch2;

	scratch =  (      ( S[xy] + S[yx] ) * ( S[yz] + S[zy] ) ) +                   
				 	    ( ( S[xz] + S[zx] ) * ( S[xx] - S[yy] + S[zz] ) );
	scratch2 = ( -1.0 * ( S[xy] - S[yx] ) * ( S[yz] - S[zy] ) ) +                   
				       ( ( S[xz] + S[zx] ) * ( S[xx] + S[yy] + S[zz] ) );
	data_size H = scratch * scratch2;
	
	scratch = ( 	   ( S[xy] + S[yx] ) * ( S[yz] - S[zy] ) ) +                   
				 		 ( ( S[xz] - S[zx] ) * ( S[xx] - S[yy] - S[zz] ) );
	scratch2 = ( -1 * ( S[xy] - S[yx] ) * ( S[yz] + S[zy] ) ) +                   
				  	    ( ( S[xz] - S[zx] ) * ( S[xx] + S[yy] - S[zz] ) );
	data_size I = scratch * scratch2;

	C[0] = D + E + F + G + H + I;
	
	data_size P,dP,sigma2;
	int iter = 0;
	data_size delta = 0.0;
	data_size sigma_old;
	while ( iter++ < 128 )
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
		#ifdef _DEBUG
		//printf("WORKERRMS: sigma: %lf delta: %lf (%d iter)\n",
		//	sigma, delta, iter);
		#endif
		
		if ( fabs( delta ) < EPS)
		{
			#ifdef _DEBUG
		//	printf("WORKERRMS: converged to %lf (%d iter)\n",sigma,iter);
			#endif
			break;
		}
	}
	if ( iter == 128 ) 
	{
		//printf("WORKERRMS: Maximum (%d) steps exceeded\n", 128);
		return;
	}
	rot[0] = inners; //sqrt(fabs(inners - 2.0*sigma));
    return;
	data_size scratch3[9];
	data_size det[4] = { 0.0, 0.0, 0.0, 0.0 };
	
	scratch3[xx] = K[ 5] - sigma; scratch3[xy] = K[6 ]; scratch3[xz] = K[7 ];
	scratch3[yx] = K[ 9]; scratch3[yy] = K[10] - sigma; scratch3[yz] = K[11];
	scratch3[zx] = K[13]; scratch3[zy] = K[14]; scratch3[zz] = K[15] -sigma;
	det[0] = determinant3by3(scratch3);

	scratch3[xx] = K[ 4]; 
	scratch3[yx] = K[ 8]; 
	scratch3[zx] = K[12]; 
	det[1] = determinant3by3(scratch3);

	scratch3[xy] = K[ 5] - sigma; 
	scratch3[yy] = K[ 9]; 
	scratch3[zy] = K[13]; 
	det[2] = determinant3by3(scratch3);

	scratch3[xz] = K[ 6]; 
	scratch3[yz] = K[10] - sigma; 
	scratch3[zz] = K[14]; 
	det[3] = determinant3by3(scratch3);
	
	det[1] = -1.0 * det[1];
	det[3] = -1.0 * det[3];
	
	data_size detsq =   det[0] * det[0] 
						+ det[1] * det[1] 
						+ det[2] * det[2] 
						+ det[3] * det[3];
	data_size detnrm = sqrt(detsq);
	det[0] /= detnrm;
	det[1] /= detnrm;
	det[2] /= detnrm;
	det[3] /= detnrm;
	
	//rot[0] = K[0];
	//rot[1] = K[1];
	//rot[2] = K[2];
	//rot[3] = K[3];
	//rot[4] = K[4];
	//rot[5] = K[5];
	//rot[6] = K[6];
	//rot[7] = K[7];
	//rot[8] = K[8];
	//rot[5] = det[0];
	//rot[6] = det[1];
	//rot[7] = det[2];
	//rot[8] = det[3];
	
	
	#ifdef _DEBUG
	//printf("WORKERRMS: det quaternion is %f %f %f %f\n",
	//	det[0],det[1],det[2],det[3]);
	#endif
	rot[0] = det[0] * det[0] + det[1] * det[1] - det[2] * det[2] - det[3] * det[3];
	rot[1] = 2 * ( det[1] * det[2] - det[0] * det[3]);
	rot[2] = 2 * ( det[0] * det[2] + det[1] * det[3]);
	rot[3] = 2 * ( det[0] * det[3] + det[1] * det[2]);
	rot[4] = det[0] * det[0] - det[1] * det[1] + det[2] * det[2] - det[3] * det[3];
	rot[5] = 2 * ( det[2] * det[3] - det[0] * det[1]);
	rot[6] = 2 * ( det[1] * det[3] - det[0] * det[2]);
	rot[7] = 2 * ( det[0] * det[1] + det[2] * det[3]);
	rot[8] = det[0] * det[0] - det[1] * det[1] - det[2] * det[2] + det[3] * det[3];
}

void rotate_atoms(__global data_size *tar,
				__private size_t *size,
				__private data_size *rot,
				__local data_size* shm)
{
	size_t i = 0, N = *size;
	size_t tid = get_local_id(0);
	data_size x = 0.0, y = 0.0, z = 0.0;
	for(i = tid*3; i < N*3 ; i+=get_local_size(0)*3)
	{
		x = tar[i];
		y = tar[i+1];
		z = tar[i+2];
		
		tar[i] = 
			( x * rot[ 0 ] )
		+	( y * rot[ 3 ] )
		+	( z * rot[ 6 ] );
		tar[i+1] = 
			( x * rot[ 1 ] )
		+	( y * rot[ 4 ] )
		+	( z * rot[ 7 ] );
		tar[i+2] = 
			( x * rot[ 2 ] )
		+	( y * rot[ 5 ] )
		+	( z * rot[ 8 ] );
		//data_size rmsd = sqrt(((tar[i] - x) * (tar[i] - x)) + ((tar[i+1] - y) * (tar[i+1] - y)) +((tar[i+2] - z) * (tar[i+2] - z)));
		//printf("BLK %zu TID %zu OLD IS %f %f %f NEW IS %f %f %f RMSD IS %f\n",blk,tid,x,y,z,tar[i],tar[i+1],tar[i+2],rmsd);
	}

}



void get_cov(    __global  data_size* ref,
				 __global  data_size* tar,
				 __private size_t  atoms,
				 __private double* cov,
				 __local   double* shm)
{
	size_t i = 0;	
	double buffer = 0.0;
	size_t x = get_local_id(0);
	size_t blk = get_global_id(1);
	size_t tid = x; //get_local_size(0)*y + x;
	size_t x_len =  get_local_size(0);
	//data_size d = (data_size)atoms;
	//data_size sum[] =  {0.0, 0.0, 0.0,0.0,0.0,0.0,0.0,0.0,0.0};
#if 0
	data_size rmean[] = {0.0,0.0,0.0};
	data_size rstd[]  = {0.0,0.0,0.0};
	data_size tmean[]   = {0.0,0.0,0.0};
	data_size tstd[]    = {0.0,0.0,0.0};
	for( i = 3*tid ; i < atoms*3  ; i+= x_len*3 )
	{
		data_size3 r = (data_size3)(ref[i],ref[i+1],ref[i+2])  ;
		data_size3 t = (data_size3)(tar[i],tar[i+1],tar[i+2])  ;
		rmean[0] += r.x / atoms;
		rmean[1] += r.y / atoms;
		rmean[2] += r.z / atoms;
		tmean[0] += t.x / atoms;
		tmean[1] += t.y / atoms;
		tmean[2] += t.z / atoms;
	}
	reduce_three(rmean    ,shm,&atoms);
	reduce_three(tmean    ,shm,&atoms);
	
	for( i = 3*tid ; i < atoms*3  ; i+= x_len*3 )
	{
		data_size3 t = (data_size3)(tar[i],tar[i+1],tar[i+2])  ;
		data_size3 r = (data_size3)(ref[i],ref[i+1],ref[i+2])  ;
		rstd[0] += ((r.x - rmean[0]) * (r.x - rmean[0])) / atoms;
		rstd[1] += ((r.y - rmean[1]) * (r.y - rmean[1])) / atoms;
		rstd[2] += ((r.z - rmean[2]) * (r.z - rmean[2])) / atoms;
		tstd[0] += ((t.x - tmean[0]) * (t.x - tmean[0])) / atoms;
		tstd[1] += ((t.y - tmean[1]) * (t.y - tmean[1])) / atoms;
		tstd[2] += ((t.z - tmean[2]) * (t.z - tmean[2])) / atoms;
	}
	reduce_three(rstd    ,shm,&atoms);
	reduce_three(tstd    ,shm,&atoms);
	rstd[0] = sqrt(rstd[0]);
	rstd[1] = sqrt(rstd[1]);
	rstd[2] = sqrt(rstd[2]);
	tstd[0] = sqrt(tstd[0]);
	tstd[1] = sqrt(tstd[1]);
	tstd[2] = sqrt(tstd[2]);
	for( i = 3*tid ; i < atoms*3  ; i+= x_len*3 )
	{
		data_size3 t = (data_size3)(tar[i],tar[i+1],tar[i+2]) ;
		data_size3 r = (data_size3)(ref[i],ref[i+1],ref[i+2]) ;
		t.x = (t.x - tmean[0]) ;
		t.y = (t.y - tmean[1]) ;
		t.z = (t.z - tmean[2]) ;
		r.x = (r.x - rmean[0]) ;
		r.y = (r.y - rmean[1]) ;
		r.z = (r.z - rmean[2]) ;
		
		cov[0] += (r.x * t.x) / atoms  ;
		cov[1] += (r.x * t.y) / atoms  ;
		cov[2] += (r.x * t.z) / atoms  ;
		cov[3] += (r.y * t.x) / atoms  ;
		cov[4] += (r.y * t.y) / atoms  ;
		cov[5] += (r.y * t.z) / atoms  ;
		cov[6] += (r.z * t.x) / atoms  ;
		cov[7] += (r.z * t.y) / atoms  ;
		cov[8] += (r.z * t.z) / atoms  ;
		
		//if(tid == 0)
			//printf("%d i is %zu going until %zu sum is %f %f %f %f %f\n",tid,i,atoms,ref[i],tar[i],sum[0],	sum[1],sum[2]);
			//printf("BLK %zu TID %zu i is %zu going until %zu sum is %f %f %f\n",blk,tid,i/3,atoms,cov[0], cov[1],cov[2]);
	}
	reduce_three(cov    ,shm,&atoms);
	reduce_three(cov + 3,shm,&atoms);
	reduce_three(cov + 6,shm,&atoms);
	cov[0] /= (rstd[0] * tstd[0]) ;
	cov[1] /= (rstd[0] * tstd[1]) ;
	cov[2] /= (rstd[0] * tstd[2]) ;
	cov[3] /= (rstd[1] * tstd[0]) ;
	cov[4] /= (rstd[1] * tstd[1]) ;
	cov[5] /= (rstd[1] * tstd[2]) ;
	cov[6] /= (rstd[2] * tstd[0]) ;
	cov[7] /= (rstd[2] * tstd[1]) ;
	cov[8] /= (rstd[2] * tstd[2]) ;
#else
	for( i = 3*tid ; i < atoms*3  ; i+= x_len*3 )
	{
		double3 t = (double3)(tar[i],tar[i+1],tar[i+2]) ;
		double3 r = (double3)(ref[i],ref[i+1],ref[i+2]) ;
		
		cov[0] += (r.x * t.x) / atoms  ;
		cov[1] += (r.x * t.y) / atoms  ;
		cov[2] += (r.x * t.z) / atoms  ;
		cov[3] += (r.y * t.x) / atoms  ;
		cov[4] += (r.y * t.y) / atoms  ;
		cov[5] += (r.y * t.z) / atoms  ;
		cov[6] += (r.z * t.x) / atoms  ;
		cov[7] += (r.z * t.y) / atoms  ;
		cov[8] += (r.z * t.z) / atoms  ;
		
		//if(tid == 0)
			//printf("%d i is %zu going until %zu sum is %f %f %f %f %f\n",tid,i,atoms,ref[i],tar[i],sum[0],	sum[1],sum[2]);
			//printf("BLK %zu TID %zu i is %zu going until %zu sum is %f %f %f\n",blk,tid,i/3,atoms,cov[0], cov[1],cov[2]);
	}
	reduce_three(cov    ,shm,&atoms);
	reduce_three(cov + 3,shm,&atoms);
	reduce_three(cov + 6,shm,&atoms);
#endif
	
	
//	cov[0] = tmean[0];
//	cov[1] = tmean[1];
//	cov[2] = tmean[2];
//	cov[3] = rmean[0];
//	cov[4] = rmean[1];
//	cov[5] = rmean[2];
	//cov[1] = reduce3(&sum[1],shm,&atoms);
	//cov[2] = reduce3(&sum[2],shm,&atoms);
	//cov[0] = sum[0];
	//cov[1] = sum[1];
	//cov[2] = sum[2];

	//for( i = 0, buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i] * tar[3 * i] + buffer;
	//cov[0] = buffer;
	//
	//for( i = 0,buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i] * tar[3 * i + 1 ] + buffer;
	//cov[1] = buffer;
	//
	//for( i = 0,buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i] * tar[3 * i + 2 ] + buffer;
	//cov[2] = buffer;
	//	
	//for( i = 0, buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i + 1] * tar[3 * i] + buffer;
	//cov[3] = buffer;
	//
	//for( i = 0,buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i + 1] * tar[3 * i + 1 ] + buffer;
	//cov[4] = buffer;
	//
	//for( i = 0,buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i + 1] * tar[3 * i + 2 ] + buffer;
	//cov[5] = buffer;
	//
	//for( i = 0, buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i + 2] * tar[3 * i] + buffer;
	//cov[6] = buffer;
	//
	//for( i = 0,buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i + 2] * tar[3 * i + 1 ] + buffer;
	//cov[7] = buffer;
	//
	//for( i = 0,buffer = 0.0; i < *sz ; i++)
	//	buffer = ref[3 * i + 2] * tar[3 * i + 2 ] + buffer;
	//cov[8] = buffer;
}


__kernel void align(__global size_t *size,
					__global data_size *data,
					__global data_size *mask,
					__global data_size *ref,
					__global data_size *out) 
{ 
	__local double shm[SHR_SIZE];
	size_t tid = get_local_id(0);
	size_t blk = get_global_id(1);
	//size_t stride = get_global_size(1);//size[0];
	data_size rot[9] = { 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0 };
	double cov[9] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };
	//if(tid < 1 && blk < stride)
	//{
		//int i = 0,j = 0;
		size_t anum = size[1];
		size_t mnum = size[2];
		//uint sz = anum*3*stride;
		//uint st = stride;
		//uint ms = mnum;

		//blk is the stride, need the offset multiplied by atoms
		size_t mask_offset = blk * mnum * 3;
		size_t data_offset = blk * anum * 3;

		get_cov(ref,mask+mask_offset,mnum,cov,shm);
		//if(tid == 0)
		//	printf("BLK %zu TID %zu COV IS %f %f %f %f %f %f %f %f %f %p\n",get_global_id(1),get_global_id(0),
		//			cov[0],cov[1],cov[2],
		//			cov[3],cov[4],cov[5],
		//			cov[6],cov[7],cov[8],cov);
		quaternion_reduce_reg(mask+mask_offset,ref,&mnum,cov,rot,shm);
		//quaternion(mask+mask_offset,ref,&mnum,cov,rot,shm);
		if(tid == 0)
			out[blk] = rot[0];
   //         for(size_t p = 0; p < 9; p++) 
   //         {
   //         	out[9*blk+p] = data[p];
   //         }
        return;

//		if(1 && tid == 0)
//			printf("BLK %zu ROT IS %f %f %f %f %f %f %f %f %f\n",blk,
//					rot[0],rot[1],rot[2],
//					rot[3],rot[4],rot[5],
//					rot[6],rot[7],rot[8]);
		//barrier(CLK_GLOBAL_MEM_FENCE);
		//rotate_atoms(data + data_offset,&anum,rot,shm);
		//barrier(CLK_GLOBAL_MEM_FENCE);
	//	rotate_atoms(mask+mask_offset,&mnum,rot,shm);
		//return;
		
		//barrier(CLK_GLOBAL_MEM_FENCE);
		//if(tid == 0)
		//{
        //    int p = 0;
        //    uint idx = 0;
        //    idx = (blk * 9);
        //    //printf("RMSD FOR FRAME IS mask %f ref %f out %f\n",mask[0],ref[0],out[0]);
		//	//printf("SIZES ARE %zu %d %d\n",stride,anum,mnum);
        //    
        //    for(p = 0; p < 9; p++) 
        //    {
        //    	out[idx+p] = rot[p];
        //    }
		//}
	//}
}

// vim:ft=c
