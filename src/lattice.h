#ifndef LATTICExH
#define LATTICExH

#include "namics.h"
#include "input.h"
#include "tools.h"

class Lattice {
public:
	Lattice(vector<Input*>,string);

virtual ~Lattice();

	string name;
	vector<Input*> In;
	int MX,MY,MZ;
	vector<int> mx;
	vector<int> my;
	vector<int> mz;
	vector<int> m;
	vector<int> jx;
	vector<int> jy;
	vector<int> n_box;
	vector<string> BC;
	int BX1,BY1,BZ1,BXM,BYM,BZM;
	int *B_X1;
	int *B_Y1;
	int *B_Z1;
	int *B_XM;
	int *B_YM;
	int *B_ZM;
	int JX,JY,JZ,M;
	bool all_lattice;
	bool ignore_sites;
	bool fcc_sites;
	bool stencil_full;
	int sub_box_on;
	Real volume;
	Real Accesible_volume;

	int gradients;
	string lattice_type;
	string geometry;
	Real offset_first_layer;
	Real bond_length;
	Real *L;
	int Z;
	Real lambda;
	Real *lambda0;
	Real *fcc_lambda0;
	Real *lambda_1;
	Real *fcc_lambda_1;
	Real *lambda1;
	Real *fcc_lambda1;
	Real *LAMBDA;
	Real l0,l1;
	int fjc, FJC;
	Real *X;
	int VarInitValue;
	string Var_type;
	int Var_target;
	int Var_step;
	int Var_end_value;

	std::vector<string> KEYS;
	std::vector<string> PARAMETERS;
	std::vector<string> VALUES;

	vector<string> ints;
	vector<string> Reals;
	vector<string> bools;
	vector<string> strings;
	vector<Real> Reals_value;
	vector<int> ints_value;
	vector<bool> bools_value;
	vector<string> strings_value;

	void DeAllocateMemory(void);
	void AllocateMemory(void);
	void push(string,Real);
	void push(string,int);
	void push(string,bool);
	void push(string,string);
	void PushOutput();
	Real* GetPointer(string,int&); 
	int* GetPointerInt(string,int&);
	int P(int,int,int); 
	int P(int,int);
	int P(int);

	int GetValue(string,int&,Real&,string&);
	//bool PutMask(int* ,vector<int>,vector<int>,vector<int>,int);
	Real GetValue(Real*,string);
	bool ReadGuess(string, Real* ,string&, vector<string>&, vector<string>&, bool&, int&, int&, int&, int&, int);
	bool StoreGuess(string,Real*,string, vector<string>,vector<string>, bool,int); 
	bool PutVarInfo(string,string,Real);
	bool UpdateVarInfo(int);
	bool ResetInitValue();
	bool CheckInput(int);	
	bool PutSub_box(int,int,int,int);
	void PutParameter(string);
	string GetValue(string);
	int PutVarScan(int, int);
	bool PrepareForCalculations(void);
	void DistributeG1(Real*, Real*, int*, int*, int*, int);
	void CollectPhi(Real*, Real*, Real*, int*, int*, int*, int);
	void ComputeGN(Real*, Real*, int*, int*, int*, int*, int*, int*, int, int);
	bool GuessVar(Real*, Real, string, Real, Real);
	bool GenerateGuess(Real*, string, string, Real, Real);

	virtual void ComputeLambdas(void)=0;
	virtual Real WeightedSum(Real*)=0;
	virtual Real Moment(Real*,Real,int) =0;
	virtual void TimesL(Real*) =0; 
	virtual void DivL(Real*) =0; 
	virtual void vtk(string, Real*,string,bool) =0;
	virtual void PutProfiles(FILE*,vector<Real*>,bool)=0;
	virtual bool PutM(void)=0;
	virtual void propagate(Real*,Real*, int, int,int)=0;
	virtual void Side(Real *, Real *, int) =0;
	virtual bool ReadRange(int*, int*, int&, bool&, string, int, string, string)=0;
	virtual bool ReadRangeFile(string,int* H_p,int&, string, string) =0; 
	virtual bool FillMask(int*, vector<int>, vector<int>, vector<int>, string)=0;
	virtual bool CreateMASK(int*, int*, int*, int, bool) =0; 
	virtual Real ComputeTheta(Real*) =0;
	virtual void UpdateEE(Real*, Real*,Real*) =0;
	virtual void UpdatePsi(Real*, Real*, Real* , Real*, int*,bool,bool)=0;
	virtual void UpdateQ(Real*,Real*,Real*,Real*,int*,bool)=0;
       virtual void outputtest(); //later weghalen	
	virtual void remove_bounds(Real*)=0;
	virtual void set_bounds(Real*)=0;
	virtual void remove_bounds(int*)=0;
	virtual void set_bounds(int*)=0;


};
#endif
