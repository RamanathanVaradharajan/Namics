#include "tools.h"
#include "namics.h"
#include "stdio.h"
#ifdef PAR_MESODYN
	#include <thrust/inner_product.h>
#endif
#define MAX_BLOCK_SZ 512
#define HALF_MAX_BLOCK_SZ 256

Real* SUM_RESULT;

#ifdef CUDA
//cublasHandle_t handle;
//cublasStatus_t stat=cublasCreate(&handle);
const int block_size = 512;
cudaStream_t CUDA_STREAMS[CUDA_NUM_STREAMS];

#if !defined(__CUDA_ARCH__) || __CUDA_ARCH__ >= 600
#else
__device__ void atomicAdd(Real* address, Real val)
{
    unsigned long long int* address_as_ull =
                             (unsigned long long int*)address;
    unsigned long long int old = *address_as_ull, assumed;
    do {
        assumed = old;
        old = atomicCAS(address_as_ull, assumed,
                        __double_as_longlong(val +
                               __longlong_as_double(assumed)));
    } while (assumed != old);
}
#endif

void Propagate_gs_locality(Real* gs, Real* gs_1, Real* G1, int JX, int JY, int JZ, int M) {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	propagate_gs_locality<<<n_blocks,block_size>>>(gs, gs_1, G1, JX, JY, JZ, M);
}

__global__ void propagate_gs_locality(Real* gs, Real* gs_1, Real* G1, int JX, int JY, int JZ, int M) {
	int index = blockIdx.x*blockDim.x+threadIdx.x;

	if (index < M-JX) {
		Real gs_register = gs[index];

		gs_register += gs_1[index-JZ];
		gs_register += gs_1[index+JZ];
		gs_register += gs_1[index-JY];
		gs_register += gs_1[index+JY];
		gs_register += gs_1[index-JX];
		gs_register += gs_1[index+JX];
		gs_register *= 1.0/6.0;
		gs_register *= G1[index];
		gs[index] = gs_register;
	}
}

#define BDIMX    	16			// tile (and threadblock) size in x
#define BDIMY    	16			// tile (and threadblock) size in y
#define WAVEFRONTS 	1			// use multiple wave-fronts
#define radius   	1 			// half of the order in space (k/2)
#define BDIMZ    	2			// tile (and threadblock) size in z

 void Second_order_fd_stencil(Real *g_output, Real *g_input, Real coeff, const int dimx, const int dimy, const int dimz) {

	dim3 dimBlock(BDIMX,BDIMY,BDIMZ);
	dim3 dimGridz(ceil(static_cast<Real>(dimz)/static_cast<Real>(BDIMX)),ceil(static_cast<Real>(dimy)/static_cast<Real>(BDIMY)),1);

	second_order_fd_stencil<<<dimGridz,dimBlock>>>(g_output, g_input, coeff, dimz, dimy, dimx);
}

__global__ void second_order_fd_stencil(Real *g_output, Real *g_input, Real coeff, const int dimx, const int dimy, const int dimz)
{ 

// **********	WARNING: OUR ARRAYS ARE DEPTH MAJOR, SO .x SHOULD BE FILLED WITH ELEMENTS IN Z AND VICE VERSA *********
	// Contains the current 2D slice
	__shared__ Real s_data[BDIMY+2*radius][BDIMX+2*radius];

	int ix  = blockIdx.x*blockDim.x + threadIdx.x;
	int iy  = blockIdx.y*blockDim.y + threadIdx.y;
	
	if (ix < dimx and iy < dimy) {
		int wavefront = blockIdx.z;
		int end = dimz/WAVEFRONTS;
		
		if (wavefront == gridDim.z-1)
			end = dimz-((end-radius)*wavefront)-radius;

		int stride  = dimx*dimy; 											// distance between 2D slices (in elements) 
		int in_idx  = ix + iy*dimx + stride*(dimz/WAVEFRONTS-radius)*wavefront;	// index for reading input
		int out_idx = 0;              										// index for writing output  

		Real infront1;								 	// variables for input “in front of” the current slice: Add more for higher order stencils.
		Real behind1; 									// variables for input “behind” the current slice: Add more for higher order stencils.
		Real current;                           		// input value in the current slice

		int tx = threadIdx.x + radius;  				// thread’s x-index into corresponding shared memory tile (adjusted for halos)
		int ty = threadIdx.y + radius; 					// thread’s y-index into corresponding shared memory tile (adjusted for halos)

		// fill the "in-front" and "behind" data
		// This is read into the GPU's register, tiles are in shared mem
		// Add in_idx += stride; for every extra stencil order register variable
		current  = g_input[in_idx];		out_idx = in_idx; in_idx += stride;
		infront1 = g_input[in_idx];		in_idx += stride;

		for(int i=radius; i < end; i++) { 

			////////////////////////////////////////// // advance the slice (move the thread-front)     
			behind1  = current;     
			current  = infront1;     
			infront1 = g_input[in_idx];

			in_idx  += stride;     
			out_idx += stride;    

			__syncthreads();

			///////////////////////////////////////// // update the data slice in smem 
			if(threadIdx.y<radius) // halo above/below     
			{     
				s_data[threadIdx.y][tx]              = g_input[out_idx-radius*dimx];     
				s_data[threadIdx.y+BDIMY+radius][tx] = g_input[out_idx+BDIMY*dimx];    
			} 
			if(threadIdx.x<radius) // halo left/right     
			{          
				s_data[ty][threadIdx.x]              = g_input[out_idx-radius];          
				s_data[ty][threadIdx.x+BDIMX+radius] = g_input[out_idx+BDIMX];     
			} // update the slice in smem     

			s_data[ty][tx] = current;

			__syncthreads(); 

			///////////////////////////////////////// // compute the output value     

			// Optionally multiply any of the following with a coefficient
			Real div  = 0; //current;  
			div += coeff*(infront1 + behind1     + s_data[ty-1][tx] + s_data[ty+1][tx] + s_data[ty][tx-1] + s_data[ty][tx+1]);
			g_output[out_idx] += div;
		} 
	}
}

void Xr_times_ci(int posi, int k_diis, int k, int m, int nvar, Real* x, Real* xR, Real* Ci)   {
	int n_blocks=(nvar)/block_size + ((nvar)%block_size == 0 ? 0:1);
	xr_times_ci<<<n_blocks,block_size, (k_diis)*sizeof(Real)>>>(posi, k_diis, k, m, nvar, x, xR, Ci);
}

__global__ void xr_times_ci(int posi, int k_diis, int k, int m, int nvar, Real* x, Real* xR, Real* Ci) {
	int idx  = blockIdx.x*blockDim.x + threadIdx.x;

 	extern __shared__ Real s_Ci[];

 	if (threadIdx.x < k_diis)
		s_Ci[threadIdx.x] = Ci[threadIdx.x];

	__syncthreads();
	 
	if (idx < nvar) {
		int thread_posi = posi;

		Real x_register = x[idx];
		x_register += s_Ci[0] * (xR+thread_posi*nvar)[idx];

		for (int i=1; i<k_diis; i++) {
			thread_posi = k-k_diis+1+i;
    		if (thread_posi<0) {
    	  		thread_posi +=m;
			}
			x_register += s_Ci[i] * xR[thread_posi*nvar+idx];
		}

		x[idx] = x_register;
	}
}

__global__ void distributeg1(Real *G1, Real *g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	if (k < Mz && j < My && i < Mx){
		int pM=jx+jy+1+k+jx*i+jy*j;
		int pos_r=JX+JY+1;
		int ii;
		int MXm1 = MX-1;
		int MYm1 = MY-1;
		int MZm1 = MZ-1;
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1) ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX;
			if (By[p]+j > MYm1) ii+=(By[p]+j-MY)*JY; else ii+=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1) ii+=(Bz[p]+k-MZ); else ii+=(Bz[p]+k);
			g1[pM]=G1[pos_r+ii];
			pM+=M;
		}
	}
}


__global__ void collectphi(Real* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int i = blockIdx.x*blockDim.x+threadIdx.x;
	int j = blockIdx.y*blockDim.y+threadIdx.y;
	int k = blockIdx.z*blockDim.z+threadIdx.z;
	int pM=jx+jy+1+jx*i+jy*j+k;
	int ii,jj,kk;
	int MXm1 = MX-1;
	int MYm1 = MY-1;
	int MZm1 = MZ-1;
	if (k < Mz && j < My && i < Mx){
		for (int p = 0;p<n_box;++p){
			if (Bx[p]+i > MXm1)  ii=(Bx[p]+i-MX)*JX; else ii=(Bx[p]+i)*JX;
			if (By[p]+j > MYm1)  jj=(By[p]+j-MY)*JY; else jj=(By[p]+j)*JY;
			if (Bz[p]+k > MZm1)  kk=(Bz[p]+k-MZ); else kk=(Bz[p]+k);
			//__syncthreads(); //will not work when two boxes are idential....
			//phi[pos_r+ii+jj+kk]+=rho[pM+jx*i+jy*j+k]*Inv_GNp;
			atomicAdd(&phi[ii+jj+kk], rho[pM]/GN[p]);
			pM+=M;
		}
	}
}

__global__ void dot(Real *v1, Real *v2, Real *out, int N)
{
	__shared__ Real cache[ MAX_BLOCK_SZ ];
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    cache[ threadIdx.x ] = 0.0;
    while( i < N ) {
        cache[ threadIdx.x ] += v1[ i ] * v2[ i ];
        i += gridDim.x * blockDim.x;
    }
    __syncthreads();
    i = HALF_MAX_BLOCK_SZ;
    while( i > 0 ) {
        if( threadIdx.x < i ) cache[ threadIdx.x ] += cache[ threadIdx.x + i ];
        __syncthreads();
        i /= 2; 
    }
    if( threadIdx.x == 0 ) atomicAdd( out, cache[ 0 ] );
}

__global__ void sum(Real *g_idata, Real *g_odata, int M) {
    __shared__ Real sdata[MAX_BLOCK_SZ]; //thread shared memory
	int tid = threadIdx.x;
    int i=threadIdx.x + blockIdx.x * blockDim.x;
    Real temp = 0;

	sdata[tid] = -DBL_MAX;
    while (i < M) {
        temp += g_idata[i];
        i += blockDim.x * gridDim.x;
    }
    sdata[tid] = temp;
    __syncthreads();

    for (int s=blockDim.x/2; s>0; s>>=1) {
        if (tid < s) {
            sdata[tid] += sdata[tid + s];
        }
        __syncthreads();
    }
    __syncthreads();
    if (tid==0) {
        atomicAdd(g_odata,sdata[0]);
    }
}

__global__ void composition(Real *phi, Real *Gf, Real *Gb, Real* G1, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (G1[idx]>0) phi[idx] += C*Gf[idx]*Gb[idx]/G1[idx];
}
__global__ void times(Real *P, Real *A, Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void times(Real *P, Real *A, int *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=A[idx]*B[idx];
}
__global__ void addtimes(Real *P, Real *A, Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx]*B[idx];
}
__global__ void norm(Real *P, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] *= C;
}
__global__ void zero(Real *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0.0;
}
__global__ void unity(Real *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 1.0;
}

__global__ void zero(int *P, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 0;
}
__global__ void cp (Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = A[idx];
}
__global__ void cp (Real *P, int *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx] = 1.0*A[idx];
}
__global__ void yisaminb(Real *Y, Real *A,Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]-B[idx];
}
__global__ void yisaplusc(Real *Y, Real *A, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+C;
}
__global__ void yisaplusb(Real *Y, Real *A,Real *B, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = A[idx]+B[idx];
}

__global__ void yplusisctimesx(Real *Y, Real *X, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*X[idx];
}
__global__ void yplusisctimesatimesb(Real *Y, Real *A, Real*B, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*A[idx]*B[idx];
}

__global__ void updatealpha(Real *Y, Real *X, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] += C*(X[idx]-1.0);
}
__global__ void picard(Real *Y, Real *X, Real C, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) Y[idx] = C*Y[idx]+(1-C)*X[idx];
}
__global__ void add(Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void add(int *P, int *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]+=A[idx];
}
__global__ void subtract(Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]-=A[idx];
}
__global__ void dubble(Real *P, Real *A, Real norm, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]*=norm/A[idx];
}
__global__ void minlog(Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx]=-log(A[idx]); else P[idx]=0;
}
__global__ void boltzmann(Real *P, Real *A, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) P[idx]=exp(-A[idx]);
}
__global__ void overwritec(Real* P, int* Mask, Real X,int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (Mask[idx]==1) P[idx] = X ; else P[idx]=0;
}
__global__ void overwritea(Real* P, int* Mask, Real* A,int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (Mask[idx]==1) P[idx] = A[idx] ; else P[idx]=0;
}

__global__ void upq(Real* g, Real* q, Real* psi, Real* eps, int jx, int jy, Real C, int* Mask, int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	Real* Px=psi+jx;
	Real* P_x=psi-jx;
	Real* Py=psi+jy;
	Real* P_y=psi-jy;
	Real* Pz=psi+1;
	Real* P_z=psi-1;
	Real* ex=eps+jx;
	Real* e_x=eps-jx;
	Real* ey=eps+jy;
	Real* e_y=eps-jy;
	Real* ez=eps+1;
	Real* e_z=eps-1;
	Real* e=eps;
	if (idx<M && Mask[idx]==1)  {
		q[idx]=(((e_x[idx]+e[idx])*P_x[idx] + (ex[idx]+e[idx])*Px[idx] +
			(e_y[idx]+e[idx])*P_y[idx] + (ey[idx]+e[idx])*Py[idx] +
			(e_z[idx]+e[idx])*P_z[idx] + (ez[idx]+e[idx])*Pz[idx]) -
			(e_x[idx]+ex[idx]+e_y[idx]+ey[idx]+e_z[idx]+ez[idx]+6*e[idx])*psi[idx])/C;
		g[idx] -=q[idx];
	}
}
__global__ void uppsi(Real* q, Real* psi, Real* X, Real* eps, int jx, int jy, Real C, int* Mask, int M) {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	Real* Px=X+jx;
	Real* P_x=X-jx;
	Real* Py=X+jy;
	Real* P_y=X-jy;
	Real* Pz=X+1;
	Real* P_z=X-1;
	Real* ex=eps+jx;
	Real* e_x=eps-jx;
	Real* ey=eps+jy;
	Real* e_y=eps-jy;
	Real* ez=eps+1;
	Real* e_z=eps-1;
	Real* e=eps;
	if (idx<M && Mask[idx]==0)  {
		psi[idx]=((e_x[idx]+e[idx])*P_x[idx] + (ex[idx]+e[idx])*Px[idx] +
			(e_y[idx]+e[idx])*P_y[idx] + (ey[idx]+e[idx])*Py[idx] +
			(e_z[idx]+e[idx])*P_z[idx] + (ez[idx]+e[idx])*Pz[idx] +
			C*q[idx])/(e_x[idx]+ex[idx]+e_y[idx]+ey[idx]+e_z[idx]+ez[idx]+6*e[idx]);
		q[idx] -=psi[idx];
	}
}

__global__ void invert(int *SKAM, int *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) SKAM[idx]=(MASK[idx]-1)*(MASK[idx]-1);
}
__global__ void invert(Real *SKAM, Real *MASK, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) SKAM[idx]=(MASK[idx]-1)*(MASK[idx]-1);
}
__global__ void addgradsquare(Real *EE, Real* X,  Real* Y, Real* Z, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) EE[idx]+=pow(X[idx]-Y[idx],2)+pow(Y[idx]-Z[idx],2);
}
__global__ void putalpha(Real *g,Real *phitot,Real *phi_side,Real chi,Real phibulk,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (phitot[idx]>0) g[idx] = g[idx] - chi*(phi_side[idx]/phitot[idx]-phibulk);
}

__global__ void putalpha(Real *g,Real *phi_side,Real chi,Real phibulk,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx] = g[idx] - chi*(phi_side[idx]-phibulk);
}

__global__ void div(Real *P,Real *A,int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) if (A[idx]>0) P[idx] /=A[idx];
}
__global__ void oneminusphitot(Real *g, Real *phitot, int M)   {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) g[idx]= 1/phitot[idx]-1;
}
__global__ void addg(Real *g, Real *phitot, Real *alpha, int M)    {
	int idx = blockIdx.x*blockDim.x+threadIdx.x;
	if (idx<M) {
		if (phitot[idx]>0) g[idx]= g[idx] -alpha[idx] +1/phitot[idx]-1; else g[idx]=0;
	}
}

//__global__ void computegn(Real *GN, Real* G, int M, int n_box)   {
//	int p= blockIdx.x*blockDim.x+threadIdx.x;
//	if (p<n_box){
//		cublasDasum(handle,M,G+p*M,1,GN+p); //in case of float we need the cublasSasum etcetera.
//	}
//}

//Define types allowed by TransferDataToHost
template void TransferDataToHost<Real>(Real*, Real*, int);
template void TransferDataToHost<int>(int*, int*, int);

template<typename T>
void TransferDataToHost(T *H, T *D, int M)    {
	cudaMemcpy(H, D, sizeof(T)*M,cudaMemcpyDeviceToHost);
}

//Define types allowed by TransferDataToDevice
template void TransferDataToDevice<Real>(Real*, Real*, int);
template void TransferDataToDevice<int>(int*, int*, int);

template<typename T>
void TransferDataToDevice(T *H, T *D, int M)    {
	cudaMemcpy(D, H, sizeof(T)*M,cudaMemcpyHostToDevice);
}

__global__ void bx(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}
__global__ void bx(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx, jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=P[bx1_jx+idx];
		P[jx_mmx+idx]=P[jx_bxm+idx];
	}
}
__global__ void b_x(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int idx, jx_mmx=jx*mmx;
	int yi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (yi<My && zi<Mz) {
		idx=jy*yi+zi;
		P[idx]=0;
		P[jx_mmx+idx]=0;
	}
}
__global__ void by(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=P[jy_by1+idx];
		P[jy_mmy+idx]=P[jy_bym+idx];
	}
}
__global__ void b_y(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int idx, jy_mmy=jy*mmy;
	int xi =blockIdx.x*blockDim.x+threadIdx.x, zi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && zi<Mz) {
		idx=jx*xi+zi;
		P[idx]=0;
		P[jy_mmy+idx]=0;
	}
}
__global__ void bz(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=P[idx+bz1];
		P[idx+mmz]=P[idx+bzm];
	}
}
__global__ void b_z(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int idx, xi =blockIdx.x*blockDim.x+threadIdx.x, yi =blockIdx.y*blockDim.y+threadIdx.y;
	if (xi<Mx && yi<My) {
		idx=jx*xi+jy*yi;
		P[idx]=0;
		P[idx+mmz]=0;
	}
}
#else
void bx(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
	}
}
void b_x(Real *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(Real *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(Real *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}
void bx(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i;
	int jx_mmx=jx*mmx;
	int jx_bxm=jx*bxm;
	int bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=P[bx1_jx+i];
		P[jx_mmx+i]=P[jx_bxm+i];
	}
}
void b_x(int *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}
void by(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}
void b_y(int *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}
void bz(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}
void b_z(int *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}
#endif

#ifdef CUDA
bool GPU_present(int deviceIndex)    {
	cudaDeviceReset();
	int deviceCount =0;
	cuDeviceGetCount(&deviceCount);
   	if (deviceCount ==0) printf("There is no device supporting Cuda.\n");
	else {
		cudaDeviceReset();
		if (deviceIndex > deviceCount-1) {
			cerr << "Error: tried to set device index to " << deviceIndex << " but the highest index number available is " << deviceCount-1
			<< ". Defaulting to 0. Press enter to continue, or CTRL+c to abort." << endl;
			cin.get();
			deviceIndex = 0;
		}
		cudaSetDevice(deviceIndex);
		cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
	}
	//if (deviceCount>0) {
	//	stat = cublasCreate(&handle);
	//	if (stat !=CUBLAS_STATUS_SUCCESS) {printf("CUBLAS handle creation failed \n");}
	//	cublasSetPointerMode(handle, CUBLAS_POINTER_MODE_DEVICE);
	//}

	const int num_streams = 4;

	for (int i = 0; i < num_streams; i++)
    	cudaStreamCreate(&CUDA_STREAMS[i]);

	return deviceCount > 0;
}
#else

bool GPU_present()    {
	return false;
}
#endif


#ifdef CUDA

int* AllIntOnDev(int N) {
  int* X;
	if (cudaSuccess != cudaMalloc((void **) &X, sizeof(int)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}

Real* AllOnDev(int N) {
  Real* X;
	cudaMalloc((void **) &X, sizeof(Real)*N);
	cudaError_t error = cudaPeekAtLastError();
	if (error != cudaSuccess)
		printf("CUDA error: %s\n", cudaGetErrorString(error));
	return X;
}

Real* AllUnpagableOnHost(int N) {
  Real* X;
	cudaHostAlloc((void **) &X, sizeof(Real)*N, cudaHostAllocDefault);
	cudaError_t error = cudaPeekAtLastError();
	if (error != cudaSuccess)
		printf("CUDA error: %s\n", cudaGetErrorString(error));
	return X;
}

int* AllManagedIntOnDev(int N) {
  int* X;
	if (cudaSuccess != cudaMallocManaged((void **) &X, sizeof(int)*N))
	printf("Memory allocation on GPU failed.\n Please reduce size of lattice and/or chain length(s) or turn on CUDA\n");
	return X;
}

Real* AllManagedOnDev(int N) {
  Real* X;
	cudaMallocManaged((void **) &X, sizeof(Real)*N);
	cudaError_t error = cudaPeekAtLastError();
	if (error != cudaSuccess)
		printf("CUDA error: %s\n", cudaGetErrorString(error));
	return X;
}

void Dot(Real &result, Real *x,Real *y, int M)   {
	cudaMemset((void**)SUM_RESULT, 0, sizeof(Real));
	//Use a pre-allocated (member) variable! Allocating memory for every call is way too costly
	//Memcopies can be masked by asynchronous transfer, allocations and frees are blocking.
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dot<<<n_blocks,block_size>>>(x,y,SUM_RESULT,M);
	TransferDataToHost(&result, SUM_RESULT, 1);
}

void Sum(Real& result, Real *x, int M)   {
	cudaMemset((void**)SUM_RESULT, 0, sizeof(Real));
	//Use a pre-allocated (member) variable! Allocating memory for every call is way too costly
	//Memcopies can be masked by asynchronous transfer, allocations and frees are blocking.
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	sum<<<n_blocks,block_size>>>(x,SUM_RESULT,M);
	TransferDataToHost(&result, SUM_RESULT, 1);
}

void Sum(int &result, int *x, int M)   {
	int *H_XXX=(int*) malloc(M*sizeof(int));
	TransferDataToHost(H_XXX, x, M);
	result=H_Sum(H_XXX,M);
	for (int i=0; i<M; i++) if (std::isnan(H_XXX[i])) cout <<" At "  << i << " NaN" << endl;
 	free(H_XXX);
}

void AddTimes(Real *P, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addtimes<<<n_blocks,block_size>>>(P,A,B,M);
}

void Times(Real *P, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}

void Times(Real *P, Real *A, int *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	times<<<n_blocks,block_size>>>(P,A,B,M);
}

void Norm(Real *P, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	norm<<<n_blocks,block_size>>>(P,C,M);
}

void Composition(Real *phi, Real *Gf, Real *Gb, Real* G1, Real C, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	composition<<<n_blocks,block_size>>>(phi,Gf,Gb,G1,C,M);
}

void Unity(Real* P, int M)   {
	cudaMemset((void**)P, 1.0, M*sizeof(Real)); //much faster than a kernel
//int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
//	unity<<<n_blocks,block_size>>>(P,M);
}

void Zero(Real* P, int M)   {
	cudaMemset((void**)P, 0.0, M*sizeof(Real)); //much faster than a kernel
/*  int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size, 0>>>(P,M); */
}

void Zero(int* P, int M)   {
	cudaMemset((void**)P, 0, M*sizeof(int)); //much faster than a kernel
/* 	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	zero<<<n_blocks,block_size, 0>>>(P,M);  */
}

void Cp(Real *P,Real *A, int M)   {
	//cudaMemcpy(P, A, sizeof(Real)*M,cudaMemcpyDeviceToDevice); //much slower than a kernel
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}

void Cp(Real *P,int *A, int M)   {
int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	cp<<<n_blocks,block_size>>>(P,A,M);
}

void YisAminB(Real *Y, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaminb<<<n_blocks,block_size>>>(Y,A,B,M);
}

void YisAplusC(Real *Y, Real *A, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusc<<<n_blocks,block_size>>>(Y,A,C,M);
}

void YisAplusB(Real *Y, Real *A, Real *B, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yisaplusb<<<n_blocks,block_size>>>(Y,A,B,M);
}

void YplusisCtimesX(Real *Y, Real *X, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesx<<<n_blocks,block_size>>>(Y,X,C,M);
}

void YplusisCtimesAtimesB(Real *Y, Real *A, Real *B, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	yplusisctimesatimesb<<<n_blocks,block_size>>>(Y,A,B,C,M);
}


void UpdateAlpha(Real *Y, Real *X, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	updatealpha<<<n_blocks,block_size>>>(Y,X,C,M);
}

  void Picard(Real *Y, Real *X, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	picard<<<n_blocks,block_size>>>(Y,X,C,M);
}

void Add(Real *P, Real *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
}

void Add(int *P, int *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	add<<<n_blocks,block_size>>>(P,A,M);
}

void Subtract(Real *P, Real *A, int M)    {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	subtract<<<n_blocks,block_size>>>(P,A,M);
}

void Dubble(Real *P, Real *A, Real norm,int M)   {
       int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	dubble<<<n_blocks,block_size>>>(P,A,norm,M);
}

void MinLog(Real *P, Real *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	minlog<<<n_blocks,block_size>>>(P,A,M);
}

void Boltzmann(Real *P, Real *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	boltzmann<<<n_blocks,block_size>>>(P,A,M);
}

void OverwriteC(Real *P, int *Mask, Real C, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	overwritec<<<n_blocks,block_size>>>(P,Mask,C,M);
}

void OverwriteA(Real *P, int *Mask, Real* A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	overwritea<<<n_blocks,block_size>>>(P,Mask,A,M);
}

void UpPsi(Real* g, Real* psi, Real* X, Real* eps, int JX, int JY, Real C, int* Mask, int M)  {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	uppsi<<<n_blocks,block_size>>>(g,psi,X,eps,JX,JY,C,Mask,M);
}

void UpQ(Real* g, Real* q, Real* psi, Real* eps, int JX, int JY, Real C, int* Mask, int M)  {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	upq<<<n_blocks,block_size>>>(g,q,psi,eps,JX,JY,C,Mask,M);
}

void Invert(Real *KSAM, Real *MASK, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	invert<<<n_blocks,block_size>>>(KSAM,MASK,M);
}
void Invert(int *KSAM, int *MASK, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	invert<<<n_blocks,block_size>>>(KSAM,MASK,M);
}
void AddGradSquare(Real* EE, Real* X, Real* Y, Real* Z, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addgradsquare<<<n_blocks,block_size>>>(EE,X,Y,Z,M);
}

void PutAlpha(Real *g, Real *phitot, Real *phi_side, Real chi, Real phibulk, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phitot,phi_side,chi,phibulk,M);
}

void PutAlpha(Real *g, Real *phi_side, Real chi, Real phibulk, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	putalpha<<<n_blocks,block_size>>>(g,phi_side,chi,phibulk,M);
}

void Div(Real *P, Real *A, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	div<<<n_blocks,block_size>>>(P,A,M);
}

void AddG(Real *g, Real *phitot, Real *alpha, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	addg<<<n_blocks,block_size>>>(g,phitot,alpha,M);
}

void OneMinusPhitot(Real *g, Real *phitot, int M)   {
	int n_blocks=(M)/block_size + ((M)%block_size == 0 ? 0:1);
	oneminusphitot<<<n_blocks,block_size>>>(g,phitot,M);
}

#ifdef PAR_MESODYN
Real ComputeResidual(Real* array, int size) {
	Real residual{0};

    auto temp_residual = thrust::minmax_element( thrust::device_pointer_cast(array), thrust::device_pointer_cast(array+size) );
    if(abs(*temp_residual.first) > abs(*temp_residual.second) ) {
      residual = abs(*temp_residual.first);
    } else {
      residual = abs(*temp_residual.second);
    }
	return residual;
}

#endif

namespace tools {

void DistributeG1(Real* G1, Real* g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//distributeg1<<<n_blocks,block_size>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
        dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        distributeg1<<<blocks,blockdim>>>(G1,g1,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);

}


void CollectPhi(Real* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	//int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
	//collectphi<<<n_blocks,block_size>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
 	 dim3 blocks(ceil(Mx/8.0),ceil(My/8.0),ceil(Mz/8.0));
        dim3 blockdim(8,8,8);
        collectphi<<<blocks,blockdim>>>(phi,GN,rho,Bx,By,Bz,MM,M,n_box,Mx,My,Mz,MX,MY,MZ,jx,jy,JX,JY);
}

}

//#ifdef CUDA
//void ComputeGN(Real *GN, Real *G, int M, int n_box)   {
//	int n_blocks=(n_box)/block_size + ((n_box)%block_size == 0 ? 0:1);
//	computegn<<<n_blocks,block_size>>>(GN,G,M,n_box);
//}
//#else
//void ComputeGN(Real* GN, Real* G, int M, int n_box)    {
//	for (int p=0; p<n_box; p++) GN[p]=Sum(G+p*M,M);
//}
//#endif

template void SetBoundaries<int>(int*, int, int, int, int, int, int, int, int, int, int, int);
template void SetBoundaries<Real>(Real*, int, int, int, int, int, int, int, int, int, int, int);

template <typename T>
void SetBoundaries(T *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(16,16);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	bx<<<dimGridx,dimBlock, 0, CUDA_STREAMS[0]>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	by<<<dimGridy,dimBlock, 0, CUDA_STREAMS[1]>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	bz<<<dimGridz,dimBlock, 0, CUDA_STREAMS[2]>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}

template void RemoveBoundaries<Real>(Real*, int, int, int, int, int, int, int, int, int, int, int);
template void RemoveBoundaries<int>(int*, int, int, int, int, int, int, int, int, int, int, int);

template <typename T>
void RemoveBoundaries(T *P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz)   {
	dim3 dimBlock(8,8);
	dim3 dimGridz((Mx+dimBlock.x+1)/dimBlock.x,(My+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridy((Mx+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	dim3 dimGridx((My+dimBlock.x+1)/dimBlock.x,(Mz+dimBlock.y+1)/dimBlock.y);
	b_x<<<dimGridx,dimBlock, 0, CUDA_STREAMS[0]>>>(P,Mx+1,My+2,Mz+2,bx1,bxm,jx,jy);
	b_y<<<dimGridy,dimBlock, 0, CUDA_STREAMS[1]>>>(P,Mx+2,My+1,Mz+2,by1,bym,jx,jy);
	b_z<<<dimGridz,dimBlock, 0, CUDA_STREAMS[2]>>>(P,Mx+2,My+2,Mz+1,bz1,bzm,jx,jy);
}

#endif

#define MAX_ITER 100
#define SIGN(a, b) ((b) >= 0.0 ? fabs(a) : -fabs(a))
#define MAX(x,y) ((x)>(y)?(x):(y))

Real PYTHAG(Real a, Real b)
{
    Real at = fabs(a), bt = fabs(b), ct, result;

    if (at > bt)       { ct = bt / at; result = at * sqrt(1.0 + ct * ct); }
    else if (bt > 0.0) { ct = at / bt; result = bt * sqrt(1.0 + ct * ct); }
    else result = 0.0;
    return(result);
}

int svdcmp(Real** a, int m, int n, Real *w, Real** v)
{
    int flag, i, its, j, jj, k, l, nm;
    Real c, f, h, s, x, y, z;
    Real anorm = 0.0, g = 0.0, scale = 0.0;
    Real *rv1;

    if (m < n)
    {
        fprintf(stderr, "#rows must be > #cols \n");
        return(0);
    }

    rv1 = (Real *)malloc((unsigned int) n*sizeof(Real));

/* Householder reduction to bidiagonal form */
    for (i = 0; i < n; i++)
    {
        /* left-hand reduction */
        l = i + 1;
        rv1[i] = scale * g;
        g = s = scale = 0.0;
        if (i < m)
        {
            for (k = i; k < m; k++)
                scale += fabs((Real)a[i][k]);
            if (scale)
            {
                for (k = i; k < m; k++)
                {
                    a[i][k] = (Real)((Real)a[i][k]/scale);
                    s += ((Real)a[i][k] * (Real)a[i][k]);
                }
                f = (Real)a[i][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[i][i] = (Real)(f - g);
                if (i != n - 1)
                {
                    for (j = l; j < n; j++)
                    {
                        for (s = 0.0, k = i; k < m; k++)
                            s += ((Real)a[i][k] * (Real)a[j][k]);
                        f = s / h;
                        for (k = i; k < m; k++)
                            a[j][k] += (Real)(f * (Real)a[i][k]);
                    }
                }
                for (k = i; k < m; k++)
                    a[i][k] = (Real)((Real)a[i][k]*scale);
            }
        }
        w[i] = (Real)(scale * g);

        /* right-hand reduction */
        g = s = scale = 0.0;
        if (i < m && i != n - 1)
        {
            for (k = l; k < n; k++)
                scale += fabs((Real)a[k][i]);
            if (scale)
            {
                for (k = l; k < n; k++)
                {
                    a[k][i] = (Real)((Real)a[k][i]/scale);
                    s += ((Real)a[k][i] * (Real)a[k][i]);
                }
                f = (Real)a[l][i];
                g = -SIGN(sqrt(s), f);
                h = f * g - s;
                a[l][i] = (Real)(f - g);
                for (k = l; k < n; k++)
                    rv1[k] = (Real)a[k][i] / h;
                if (i != m - 1)
                {
                    for (j = l; j < m; j++)
                    {
                        for (s = 0.0, k = l; k < n; k++)
                            s += ((Real)a[k][j] * (Real)a[k][i]);
                        for (k = l; k < n; k++)
                            a[k][j] += (Real)(s * rv1[k]);
                    }
                }
                for (k = l; k < n; k++)
                    a[k][i] = (Real)((Real)a[k][i]*scale);
            }
        }
        anorm = MAX(anorm, (fabs((Real)w[i]) + fabs(rv1[i])));
    }

    /* accumulate the right-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        if (i < n - 1)
        {
            if (g)
            {
                for (j = l; j < n; j++)
                    v[j][i] = (Real)(((Real)a[j][i] / (Real)a[l][i]) / g);
                    /* Real division to avoid underflow */
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < n; k++)
                        s += ((Real)a[k][i] * (Real)v[k][j]);
                    for (k = l; k < n; k++)
                        v[k][j] += (Real)(s * (Real)v[k][i]);
                }
            }
            for (j = l; j < n; j++)
                v[i][j] = v[j][i] = 0.0;
        }
        v[i][i] = 1.0;
        g = rv1[i];
        l = i;
    }

    /* accumulate the left-hand transformation */
    for (i = n - 1; i >= 0; i--)
    {
        l = i + 1;
        g = (Real)w[i];
        if (i < n - 1)
            for (j = l; j < n; j++)
                a[j][i] = 0.0;
        if (g)
        {
            g = 1.0 / g;
            if (i != n - 1)
            {
                for (j = l; j < n; j++)
                {
                    for (s = 0.0, k = l; k < m; k++)
                        s += ((Real)a[i][k] * (Real)a[j][k]);
                    f = (s / (Real)a[i][i]) * g;
                    for (k = i; k < m; k++)
                        a[j][k] += (Real)(f * (Real)a[i][k]);
                }
            }
            for (j = i; j < m; j++)
                a[i][j] = (Real)((Real)a[i][j]*g);
        }
        else
        {
            for (j = i; j < m; j++)
                a[i][j] = 0.0;
        }
        ++a[i][i];
    }

    for (int i=0; i<n; i++)
      for (int j=0; j<n; j++)
        if (a[i][j] != a[i][j])
            throw -3;

    /* diagonalize the bidiagonal form */
    for (k = n - 1; k >= 0; k--)
    {                             /* loop over singular values */
        for (its = 0; its < 30; its++)
        {                         /* loop over allowed iterations */
            flag = 1;
            for (l = k; l >= 0; l--)
            {                     /* test for splitting */
                nm = l - 1;
                if (fabs(rv1[l]) + anorm == anorm)
                {
                    flag = 0;
                    break;
                }
                if (fabs((Real)w[nm]) + anorm == anorm)
                    break;
            }
            if (flag)
            {
                c = 0.0;
                s = 1.0;
                for (i = l; i <= k; i++)
                {
                    f = s * rv1[i];
                    if (fabs(f) + anorm != anorm)
                    {
                        g = (Real)w[i];
                        h = PYTHAG(f, g);
                        w[i] = (Real)h;
                        h = 1.0 / h;
                        c = g * h;
                        s = (- f * h);
                        for (j = 0; j < m; j++)
                        {
                            y = (Real)a[nm][j];
                            z = (Real)a[i][j];
                            a[nm][j] = (Real)(y * c + z * s);
                            a[i][j] = (Real)(z * c - y * s);
                        }
                    }
                }
            }
            z = (Real)w[k];
            if (l == k)
            {                  /* convergence */
                if (z < 0.0)
                {              /* make singular value nonnegative */
                    w[k] = (Real)(-z);
                    for (j = 0; j < n; j++)
                        v[j][k] = (-v[j][k]);
                }
                break;
            }
            if (its >= 30) {
                free((void*) rv1);
                fprintf(stderr, "No convergence after 30! iterations \n");
                return(0);
            }

            /* shift from bottom 2 x 2 minor */
            x = (Real)w[l];
            nm = k - 1;
            y = (Real)w[nm];
            g = rv1[nm];
            h = rv1[k];
            f = ((y - z) * (y + z) + (g - h) * (g + h)) / (2.0 * h * y);
            g = PYTHAG(f, 1.0);
            f = ((x - z) * (x + z) + h * ((y / (f + SIGN(g, f))) - h)) / x;

            /* next QR transformation */
            c = s = 1.0;
            for (j = l; j <= nm; j++)
            {
                i = j + 1;
                g = rv1[i];
                y = (Real)w[i];
                h = s * g;
                g = c * g;
                z = PYTHAG(f, h);
                rv1[j] = z;
                c = f / z;
                s = h / z;
                f = x * c + g * s;
                g = g * c - x * s;
                h = y * s;
                y = y * c;
                for (jj = 0; jj < n; jj++)
                {
                    x = (Real)v[jj][j];
                    z = (Real)v[jj][i];
                    v[jj][j] = (Real)(x * c + z * s);
                    v[jj][i] = (Real)(z * c - x * s);
                }
                z = PYTHAG(f, h);
                w[j] = (Real)z;
                if (z)
                {
                    z = 1.0 / z;
                    c = f * z;
                    s = h * z;
                }
                f = (c * g) + (s * y);
                x = (c * y) - (s * g);
                for (jj = 0; jj < m; jj++)
                {
                    y = (Real)a[j][jj];
                    z = (Real)a[i][jj];
                    a[j][jj] = (Real)(y * c + z * s);
                    a[i][jj] = (Real)(z * c - y * s);
                }
            }
            rv1[l] = 0.0;
            rv1[k] = f;
            w[k] = (Real)x;
        }
    }
    //free((void*) rv1);
    free (rv1);
    return(1);
}
