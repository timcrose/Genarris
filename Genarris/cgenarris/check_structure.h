#ifndef _CHECK_STRUCTURE_H
#define _CHECK_STRUCTURE_H

int check_structure(crystal random_crystal, float sr);
void vector_cpy(float A[], float B[][3], int index);
float pdist(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );
			
float pdist_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  );

void pdist_2(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  );
			
void pdist_2_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  );
			
int check_pair_tier2( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *atom_vdw,
				int i,
				int j,
				int N,
				float sr);

int check_self_tier2(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *atom_vdw, int i,int  N,float  sr, float *bond_length);

int check_pair_tier3( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *atom_vdw,
				int i,
				int j,
				int N,
				float sr);

int check_self_tier3(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *atom_vdw, int i,int  N,float  sr, float *bond_length);

void convert_atom2atom_vdw(char *atom,float *atom_vdw, int num_atoms);

void calc_bond_length(float *bond_length,float *X, float *Y, float *Z, int N);


int fast_screener(crystal xtal, float sr, float *atom_vdw);


int check_self_vdw(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length);
	
int check_pair_vdw( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms);

int fast_screener_vdw(crystal xtal, float *vdw_matrix);
int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2);

int check_pair_vdw_tier2( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms);

int check_self_vdw_tier2(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length);

int check_pair_vdw_tier3( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms);

int check_self_vdw_tier3(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length);

#endif
