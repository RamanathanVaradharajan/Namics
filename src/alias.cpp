#include "alias.h"

Alias::Alias(vector<Input*> In_,vector<Lattice*> Lat_, string name_) {
	In=In_; name=name_;   Lat=Lat_;
	KEYS.push_back("value");
}
Alias::~Alias() {
	DeAllocateMemory();
}
void Alias::DeAllocateMemory(){
if (debug) cout <<"Destructor for alias " + name << endl;
	free(H_phi);
	if (clamp) free(rho);
#ifdef CUDA
	cudaFree(phi);
	if (clamp) {
		cudaFree(rho);
		cudaFree(phi);
	} else cudaFree(phi);
#else
#endif
}

void Alias::AllocateMemory(int Clamp_nr, int n_box) {
if (debug) cout <<"AllocateMemory in Alias " + name << endl;
	int M=Lat[0]->M;
	clamp = Clamp_nr>0;
	int m=0;
	if (clamp) m=Lat[0]->m[Clamp_nr];
	H_phi = (Real*) malloc(M*sizeof(Real)); H_Zero(H_phi,M);
#ifdef CUDA
	if (clamp) {
		rho=(Real*)AllOnDev(m*n_box);
		phi=(Real*)AllOnDev(M); Zero(phi,M);
	} else {
		phi=(Real*)AllOnDev(M); Zero(phi,M);
		rho=phi;
	}
#else
	if (clamp) {
		rho=(Real*) malloc(m*n_box*sizeof(Real));
		phi=H_phi;
	} else {
		phi=H_phi;
		rho=phi;
	}
#endif

}

void Alias::PrepareForCalculations() {
if (debug) cout <<"PrepareForCalculations in Alias " + name << endl;

	int M= Lat[0]->M;
	Zero(phi,M);
}


void Alias::PutParameter(string new_param) {
if (debug) cout <<"PutParameter in Alias " + name << endl;
	KEYS.push_back(new_param);
}

bool Alias::CheckInput(int start) {
if (debug) cout <<"CheckInput in Alias " + name << endl;
	bool success=true;
	success= In[0]->CheckParameters("alias",name,start,KEYS,PARAMETERS,VALUES);
	if (success) {
		value =0;
		if (GetValue("value").size()>0) {
			int defaultvalue=12345678;
			value=In[0]->Get_int(GetValue("value"),defaultvalue);
			if (value==12345678) {composition=GetValue("value"); value=-1; }
			else {
				if (value < 0|| value >1e7 ) {cout <<"In alias " + name + " the numerical value of 'value' out of range 0 ... 1e7" << endl; success=false; }
			}
		} else {
			success=false;
			cout << "In alias " << name << " the property 'value' is not found. Enter a valid molecular 'fragment' or a positive integer." << endl;
		}
	}
	return success;
}

string Alias::GetValue(string parameter){
if (debug) cout <<"GetValue in Alias " + name << endl;
	int i=0;
	int length = PARAMETERS.size();
	while (i<length) {
		if (parameter==PARAMETERS[i]) {
			return VALUES[i];
		}
		i++;
	}
	return "" ;
}

void Alias::push(string s, Real X) {
if (debug) cout <<"push (Real) in Alias " + name << endl;
	Reals.push_back(s);
	Reals_value.push_back(X);
}
void Alias::push(string s, int X) {
if (debug) cout <<"push (int) in Alias " + name << endl;
	ints.push_back(s);
	ints_value.push_back(X);
}
void Alias::push(string s, bool X) {
if (debug) cout <<"push (boool) in Alias " + name << endl;
	bools.push_back(s);
	bools_value.push_back(X);
}
void Alias::push(string s, string X) {
if (debug) cout <<"push (string) in Alias " + name << endl;
	strings.push_back(s);
	strings_value.push_back(X);
}
void Alias::PushOutput() {
if (debug) cout <<"PushOutput in Alias " + name << endl;
	strings.clear();
	strings_value.clear();
	bools.clear();
	bools_value.clear();
	Reals.clear();
	Reals_value.clear();
	ints.clear();
	ints_value.clear();

	push("value",value);
#ifdef CUDA
	TransferDataToHost(H_phi,phi,Lat[0]->M);
#endif
}

Real* Alias::GetPointer(string s) {
if (debug) cout <<"GetPointer in Alias " + name << endl;
	return NULL;
}


int Alias::GetValue(string prop,int &int_result,Real &Real_result,string &string_result){
if (debug) cout <<"GetValue (long)  in Alias " + name << endl;
	int i=0;
	int length = ints.size();
	while (i<length) {
		if (prop==ints[i]) {
			int_result=ints_value[i];
			return 1;
		}
		i++;
	}
	i=0;
	length = Reals.size();
	while (i<length) {
		if (prop==Reals[i]) {
			Real_result=Reals_value[i];
			return 2;
		}
		i++;
	}
	i=0;
	length = bools.size();
	while (i<length) {
		if (prop==bools[i]) {
			if (bools_value[i]) string_result="true"; else string_result="false";
			return 3;
		}
		i++;
	}
	i=0;
	length = strings.size();
	while (i<length) {
		if (prop==strings[i]) {
			string_result=strings_value[i];
			return 3;
		}
		i++;
	}
	return 0;
}
