#ifndef HOST_ONLY_TOOLSxH
#define HOST_ONLY_TOOLSxH

#include <numeric>
#include <algorithm>
#include <functional>
#if __SSE__
	#include <smmintrin.h>
#endif
#include <cmath>

struct saxpy_functor
{
    const double a;

    saxpy_functor(double _a) : a(_a) {}

      double operator()(const double& x, const double& y) const { 
            return a * x + y;
        }
};

struct const_multiply_functor
{
    const double a;

    const_multiply_functor(double _a) : a(_a) {}

      double operator()(const double& x, const double& y) const { 
            return a * x * y;
        }
};


struct order_param_functor
{

    order_param_functor() {}

        double operator()(const double& x, const double& y) const { 
            return pow(x-y,2);
        }
};

struct is_negative_functor
{

  const double tolerance{0};

  is_negative_functor(double _tolerance=0) : tolerance(_tolerance) {}

  bool operator()(const double &x) const
  {
    return x < 0-tolerance || x > 1+tolerance;
  }
};

struct is_not_unity_functor
{
  const double tolerance{0};

  is_not_unity_functor(double _tolerance = 1e-4 ) : tolerance(_tolerance) {}

  bool operator()(const double &x) const
  {
    bool result{0};

    if (x > (1+tolerance) || x < (1-tolerance))
      result = 1;

    return result;
  }
};

typedef double Real;

template <typename T>
inline void Times(Real* P, Real* A, T* B, int M) {
  std::transform(A, A + M, B, P, std::multiplies<Real>());
}

template <typename T>
inline void Add(T* P, T* A, int M) {
  std::transform(P, P + M, A, P, std::plus<T>());
}
template <typename T>
inline void Subtract(T* P, T* A, int M) {
  #pragma GCC ivdep
  for (int i = 0 ; i < M ; ++i)
	P[i] -= A[i];
  //std::transform(P, P + M, A, P, std::plus<T>());
}

template <typename T>
inline void Sum(T &result, T *x, int M)   {
  result = std::accumulate(x, x+M, 0.0);
//  result = 0;
//  for (int i=0; i<M; i++) result +=x[i];
}

template <typename T>
inline void Invert(T* KSAM, T* MASK, int M) {
  std::transform(MASK, MASK + M, KSAM, KSAM, [](Real A, Real B) { if (A==0.0) return 1.0; else return 0.0; });
  //std::transform(MASK, MASK + M, KSAM, KSAM, (1.0 - std::placeholders::_1) *(1.0 - std::placeholders::_1);
}

template <typename T>
inline void Zero(T* P, int M) {
  std::fill(P, P + M, 0);
}

template <typename T>
inline void Cp(Real* P, T* A, int M) {
  std::copy(A, A + M, P);
}

template <typename T>
void bx(T* P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
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

template<typename T>
void b_x(T *P, int mmx, int My, int Mz, int bx1, int bxm, int jx, int jy)   {
	int i, jx_mmx=jx*mmx;// jx_bxm=jx*bxm, bx1_jx=bx1*jx;
	for (int y=0; y<My; y++)
	for (int z=0; z<Mz; z++){
		i=jy*y+z;
		P[i]=0;
		P[jx_mmx+i]=0;
	}
}

template<typename T>
void by(T *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy, jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=P[jy_by1+i];
		P[jy_mmy+i]=P[jy_bym+i];
	}
}

template<typename T>
void b_y(T *P, int Mx, int mmy, int Mz, int by1, int bym, int jx, int jy)   {
	int i, jy_mmy=jy*mmy;// jy_bym=jy*bym, jy_by1=jy*by1;
	for (int x=0; x<Mx; x++)
	for (int z=0; z<Mz; z++) {
		i=jx*x+z;
		P[i]=0;
		P[jy_mmy+i]=0;
	}
}

template<typename T>
void bz(T *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=P[i+bz1];
		P[i+mmz]=P[i+bzm];
	}
}

template<typename T>
void b_z(T *P, int Mx, int My, int mmz, int bz1, int bzm, int jx, int jy)   {
	int i;
	for (int x=0; x<Mx; x++)
	for (int y=0; y<My; y++) {
		i=jx*x+jy*y;
		P[i]=0;
		P[i+mmz]=0;
	}
}

template <typename T>
inline void SetBoundaries(T* P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  bx(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy);
  by(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy);
  bz(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy);
}

template <typename T>
inline void RemoveBoundaries(T* P, int jx, int jy, int bx1, int bxm, int by1, int bym, int bz1, int bzm, int Mx, int My, int Mz) {
  b_x(P, Mx + 1, My + 2, Mz + 2, bx1, bxm, jx, jy);
  b_y(P, Mx + 2, My + 1, Mz + 2, by1, bym, jx, jy);
  b_z(P, Mx + 2, My + 2, Mz + 1, bz1, bzm, jx, jy);
}

#if __SSE__
template<typename T>
void Dot(T &result, T *x,T *y, int M)   {
	T z = 0.0;
	result = 0.0;
	T ftmp[2] = { 0.0, 0.0 };
	__m128d mres;

	if ((M / 2) != 0) {
		mres = _mm_load_sd(&z);
		for (int i = 0; i < M / 2; i++)
			mres = _mm_add_pd(mres, _mm_mul_pd(_mm_loadu_pd(&x[2*i]),
			_mm_loadu_pd(&y[2*i])));

		_mm_store_pd(ftmp, mres);

		result = ftmp[0] + ftmp[1];
}

	if ((M % 2) != 0) {
		for (int i = M - M % 2; i < M; i++)
			result += x[i] * y[i];
	}
}

#else

template< typename T>
void Dot(T &result, T *x,T *y, int M)   {
	result = 0.0;
	for (int i = 0 ; i < M ; i++)
		result += x[i] * y[i];
}

#endif

template<typename T>
void AddTimes(T *P, T *A, T *B, int M)   {
	for (int i=0; i<M; i++) P[i]+=A[i]*B[i];
}

template<typename T>
void Composition(T* phi, T* Gf, T* Gb, T* G1, T C, int M)   {
	for (int i=0; i<M; i++) if (G1[i]>0) phi[i]+=C*Gf[i]*Gb[i]/G1[i];
}

template<typename T, typename D>
void Norm(T *P, D C, int M)   {
     for_each(P, P+M, [C](T& a) { a *= C;} );
	 //for (int i=0; i<M; i++) P[i] *= C;
	}

template<typename T>
void Unity(T* P, int M)   {
  std::fill(P, P+M, 1);
}

template<typename T>
void Assign(T* P, T C, int M)   {
  std::fill(P, P+M, C);
}

template<typename T>
void Assign(T* P, T* C, int M)   {
  std::copy(P, P+M, C);
}

template<typename T>
void YisAminB(T *Y, T *A, T *B, int M)   {
  std::transform(A, A + M, B, Y, std::minus<Real>());
}

template<typename T>
void YisAplusC(T *Y, T *A, T C, int M)   {
	for (int i=0; i<M; i++) Y[i] = A[i]+C;
}

template<typename T>
void YisAplusB(T *Y, T *A, T *B, int M)   {
  std::transform(A, A + M, B, Y, std::plus<Real>());
}

template<typename T>
void YplusisCtimesX(T *Y, T *X, T C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*X[i];
}

template<typename T>
void Xr_times_ci(int posi, int k_diis, int k, int m, int nvar, T* x, T* xR, T* Ci) {
	YplusisCtimesX(x,xR+posi*nvar,Ci[0],nvar); //pv = Ci[0]*xR[0];

	for (int i=1; i<k_diis; i++) {
		posi = k-k_diis+1+i;
    	if (posi<0) {
      		posi +=m;
		}
		YplusisCtimesX(x,xR+posi*nvar,Ci[i],nvar);
	}
}

template<typename T>
void UpdateAlpha(T *Y, T *X, T C, int M)    {
	for (int i=0; i<M; i++) Y[i] += C*(X[i]-1.0);
}

template<typename T>
void Picard(T *Y, T *X, T C, int M)    {
	for (int i=0; i<M; i++) Y[i] = C*Y[i]+(1.0-C)*X[i];
}

template<typename T>
void Dubble(Real *P, T *A, T norm,int M)   {
	for (int i=0; i<M; i++) P[i]*=norm/A[i];
}

template<typename T>
void MinLog(Real *P, T *A, int M)   {
	for (int i=0; i<M; i++) if (A[i]>0) P[i]=-log(A[i]); else P[i]=0;
}

template<typename T>
void Boltzmann(Real *P, T *A, int M)   {
	for (int i=0; i<M; i++) P[i]=exp(-A[i]);
}

template<typename T>
void AddGradSquare(T* EE, T* X, T* Y, T* Z, int M)    {
	for (int i=0; i<M; i++) EE[i] += pow(X[i]-Y[i],2)+pow(Y[i]-Z[i],2);
}

template<typename T>
void PutAlpha(T *g, T *phitot, T *phi_side, T chi, T phibulk, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0) g[i] = g[i] - chi*(phi_side[i]/phitot[i]-phibulk);
}

template<typename T>
void PutAlpha(T *g, T *phi_side, T chi, T phibulk, int M)   {
	std::transform(phi_side, phi_side+M, g, g, std::placeholders::_2 - (chi*std::placeholders::_1-phibulk)) ;
	//for (int i=0; i<M; i++) g[i] = g[i] - chi*(phi_side[i]-phibulk);
}

template<typename T>
void Div(T *P, T *A, int M)   {
	for (int i=0; i<M; i++) if (A[i]!=0) P[i]/=A[i]; else P[i]=0;
}

template<typename T>
void AddG(T *g, T *phitot, T *alpha, int M)   {
	for (int i=0; i<M; i++) if (phitot[i]>0)  g[i]= g[i] -alpha[i] +1/phitot[i]-1.0; else g[i]=0;
}

template<typename T>
void OneMinusPhitot(T *g, T *phitot, int M)   {
	for (int i=0; i<M; i++) g[i]= 1/phitot[i]-1;
}

namespace tools {

template<typename T>
void DistributeG1(T* G1, Real* g1, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	int ii=0,jj=0,kk=0;

	for (int p=0; p<n_box; p++) { pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p];
		for (int i=1; i<Mx+1; i++) { ii+=jx; jj=0; if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					g1[pos_l+ii+jj+kk]=G1[pos_x+pos_y+pos_z];
				}
			}
		}
	}
}

template<typename T>
void CollectPhi(T* phi, Real* GN, Real* rho, int* Bx, int* By, int* Bz, int MM, int M, int n_box, int Mx, int My, int Mz, int MX, int MY, int MZ, int jx, int jy, int JX, int JY) {
	int pos_l=-M;
	int pos_x,pos_y,pos_z;
	int Bxp,Byp,Bzp;
	Real Inv_H_GNp;
	int ii=0,jj=0,kk=0;
	for (int p=0; p<n_box; p++) {pos_l +=M; ii=0; Bxp=Bx[p]; Byp=By[p]; Bzp=Bz[p]; Inv_H_GNp=1.0/GN[p];
		for (int i=1; i<Mx+1; i++) {ii+=jx; jj=0;  if (Bxp+i>MX) pos_x=(Bxp+i-MX)*JX; else pos_x = (Bxp+i)*JX;
			for (int j=1; j<My+1; j++) {jj+=jy;  kk=0; if (Byp+j>MY) pos_y=(Byp+j-MY)*JY; else pos_y = (Byp+j)*JY;
				for (int k=1; k<Mz+1; k++) { kk++; if (Bzp+k>MZ) pos_z=(Bzp+k-MZ); else pos_z = (Bzp+k);
					phi[pos_x+pos_y+pos_z]+=rho[pos_l+ii+jj+kk]*Inv_H_GNp;
				}
			}
		}
	}
}

}

template<typename T, typename D>
void OverwriteC(T *P, D *Mask, T C, int M) {
	for (int i=0; i<M; i++) if (Mask[i]==1) P[i]=C; else P[i]=0;
}

template<typename T, typename D>
void OverwriteA(T *P, D *Mask,T* A,int M) {
	for (int i=0; i<M; i++) if (Mask[i]==1) P[i]=A[i]; else P[i]=0;
}

Real pythag(Real, Real);
int svdcmp(Real**, int, int, Real*, Real**);

#endif
