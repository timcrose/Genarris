#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "read_input.h"
#include "spg_generation.h"
#include "crystal_utils.h"
#include "check_structure.h"
#include "algebra.h"

float CONSTANT_TOL = 0.001;
float MOL_SIM_TOL = 0.05;

static inline int int_floor(float x)
{
  int i = (int)x; /* truncate */
  return i - ( i > x ); /* convert trunc to floor */
}



/* The most function of this file. returns the distance between
 * two points under periodic boundary. T and T_inv are lattice
 * vector matrix and its inverse. T is in row major form. x1,x2,x3 are
 * the first point and y's are second in cartesian
 * returns the shortest distance as a float.
 */
float pdist(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  )
{
	float p_dist = pdist_appx(T,
			T_inv,
			x1, 
			x2,
			x3,
			y1,
			y2,
			y3  );
			
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((long)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

	//find cartesian distance vector in the bounding box
	//reduced cartesian distance vector
	float red_cart_distance[3];
	for(int i = 0; i < 3; i++)
		red_cart_distance[i] = T[0][i] * fractional_distance[0] + 
							   T[1][i] * fractional_distance[1] +
							   T[2][i] * fractional_distance[2];
	
	//use pdist as radius of sphere are search
	int search_safety = 2;
	int limx = abs( int_floor ( p_dist/ T[0][0] ) ) + 1 + search_safety;
	int limy = abs( int_floor ( p_dist/ T[1][1] ) ) + 1 + search_safety;
	int limz = abs( int_floor ( p_dist/ T[2][2] ) ) + 1 + search_safety;
	
	float test_dist[3] = {0,0,0};
	
	for (int i = -limx; i <= limx; i++)
	for (int j = -limy; j <= limy; j++)
	for (int k = -limz; k <= limz; k++)
	{
		test_dist[0] = i*T[0][0] + j*T[1][0] + k*T[2][0] - red_cart_distance[0];
		test_dist[1] = j*T[1][1] + k*T[2][1] - red_cart_distance[1];
		test_dist[2] = k*T[2][2] - red_cart_distance[2];
		
		if (vector3_norm(test_dist) < p_dist )
			p_dist = vector3_norm(test_dist);
	}

	return p_dist;
}

/* Almost the same as p_dist, but returns the shortest and the second
 * shortest distance. return value is in argument. used for self image 
 * check
 */
void pdist_2(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  )
{
	*p_dist = 0.0;
	*p_dist_2 = 0.0;
	
	pdist_2_appx(T,
				T_inv,
				x1,
				x2,
				x3,
				y1,
				y2,
				y3,
				p_dist,
				p_dist_2);
	
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((int)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

	float red_cart_distance[3];
	for(int i = 0; i < 3; i++)
		red_cart_distance[i] = T[0][i] * fractional_distance[0] + 
							   T[1][i] * fractional_distance[1] +
							   T[2][i] * fractional_distance[2];
	
	//use pdist as radius of sphere are search
	int search_safety = 2;
	int limx = abs( int_floor ( *p_dist/ T[0][0] ) ) + 1 + search_safety;
	int limy = abs( int_floor ( *p_dist/ T[1][1] ) ) + 1 + search_safety;
	int limz = abs( int_floor ( *p_dist/ T[2][2] ) ) + 1 + search_safety;
	
	float test_dist[3] = {0,0,0};
	
	for (int i = -limx; i <= limx; i++)
	for (int j = -limy; j <= limy; j++)
	for (int k = -limz; k <= limz; k++)
	{
		test_dist[0] = i*T[0][0] + j*T[1][0] + k*T[2][0] - red_cart_distance[0];
		test_dist[1] = i*T[0][1] + j*T[1][1] + k*T[2][1] - red_cart_distance[1];
		test_dist[2] = i*T[0][2] + j*T[1][2] + k*T[2][2] - red_cart_distance[2];
		
		float norm = vector3_norm(test_dist);
		if ( norm + CONSTANT_TOL < *p_dist )
		{/*
            printf("ijl = %d %d %d , norm = %f, red_cart_dist=%f\n", i, j, k, norm, vector3_norm(red_cart_distance));
            printf("red_cart_dist = %f %f %f\n", red_cart_distance[0], red_cart_distance[1], red_cart_distance[2]);
            printf("x =  %f, %f, %f \n", x1, x2, x3);
            printf("y =  %f, %f, %f \n", y1, y2, y3);
            printf("fractional dist =  %f %f %f \n", fractional_distance[0], fractional_distance[1], fractional_distance[2] );
           */
           // print_mat3b3(T);
           // print_mat3b3(T_inv);
			*p_dist_2 = *p_dist;
			*p_dist = norm;
		}
		else if ( norm + CONSTANT_TOL < *p_dist_2 && *p_dist > norm + CONSTANT_TOL)
			*p_dist_2 = norm;
	}
	
}



/*aproximate version of pdist()
 */ 
float pdist_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3  )
{
    //intialising variables
	float p_dist = 0;

	static float Q[8][3]={ {0,0,0},{1,0,0},{0,1,0},{0,0,1},
						   {1,1,0},{0,1,1},{1,0,1},{1,1,1}  };
		
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
	//vector3_subtract(x,y,cartesian_distance);

    /*computes tthe distance in fractinal space and reduces it 
     * inside first unit cube
     * #TODO: make this outside pdist().
     *  Needs to be done only once per crystal
     */
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((long)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

	for (int z = 0; z < 8; z++)
	{
		float A[3] = {Q[z][0], Q[z][1], Q[z][2]};
		float dist_corner[3];

        dist_corner[0] = fractional_distance[0] - A[0];
		dist_corner[1] = fractional_distance[1] - A[1];
        dist_corner[2] = fractional_distance[2] - A[2];

		//distance vector to 8 corners in cartesian.
		float dist_z = 0;
		for (int i = 0; i < 3; i++)
		{
			float sum = 0;
			for (int j = 0; j < 3 ;j++)
				sum = sum + dist_corner[j] * T[j][i];	
			dist_z += sum*sum;
		}
        // length odistnce vector to the zth corner 
        // for finding the minimum distance (min dist_z)
		if (z == 0)
			p_dist = dist_z;
		else if (p_dist > dist_z)
			p_dist = dist_z;
	}
	return sqrt(p_dist); 
}

/* Almost the same as p_dist, but returns the shortest and the second
 * shortest distance. return value is in argument. used for self image 
 * check
 */
void pdist_2_appx(float T[3][3],
			float T_inv[3][3],
			float x1, 
			float x2,
			float x3,
			float y1,
			float y2,
			float y3,
			float *p_dist,
			float *p_dist_2  )
{
    //intialising variables
	*p_dist = 0;
	*p_dist_2 = 0;
	static float Q[8][3]={ {0,0,0},{1,0,0},{0,1,0},{0,0,1},
						   {1,1,0},{0,1,1},{1,0,1},{1,1,1}  };
		
	float cartesian_distance[3] = {x1 - y1, x2 - y2, x3 - y3};
	//vector3_subtract(x,y,cartesian_distance);

    /*computes tthe distance in fractinal space and reduces it 
     * inside first unit cube
     * #TODO: make this outside pdist().
     *  Needs to be done only once per crystal
     */
    float fractional_distance[3];
	for (int i = 0; i < 3; i++)
	{	
		float sum = 0;
		for (int j = 0; j < 3; j++)
			sum = sum + T_inv[j][i] * cartesian_distance[j];	

		float frac_part = sum - ((long)sum); //if -ve add one
		if (frac_part < 0)
			fractional_distance[i] = frac_part + 1;	
		else
			fractional_distance[i] = frac_part;
	}

    //computes the distances to all the corners of the unit cube
	for (int z = 0; z < 8; z++)
	{
		float A[3] = {Q[z][0], Q[z][1], Q[z][2]};
		float dist_corner[3];

        dist_corner[0] = fractional_distance[0] - A[0];
		dist_corner[1] = fractional_distance[1] - A[1];
        dist_corner[2] = fractional_distance[2] - A[2];

		//distance vector to 8 corners in cartesian.
		float dist_z = 0;
		for (int i = 0; i < 3; i++)
		{
			float sum = 0;
			for (int j = 0; j < 3 ;j++)
				sum = sum + dist_corner[j] * T[j][i];	
			dist_z += sum*sum;
		}
        // length odistnce vector to the zth corner 
        // for finding the minimum distance (min dist_z)
		if (z == 0)
		{
			*p_dist = dist_z;
		}
		else if (z == 1)
		{
			if (*p_dist > dist_z)
				{*p_dist_2 = *p_dist; *p_dist = dist_z;}
			else
				*p_dist_2 = dist_z;
		}
		else if (*p_dist > dist_z)
		{
			*p_dist_2 = *p_dist;
			*p_dist = dist_z;
		}
		else if (*p_dist_2 > dist_z)
		{
			*p_dist_2 = dist_z;
		}
		
	}
	*p_dist 	= sqrt(*p_dist);
	*p_dist_2	= sqrt(*p_dist_2);

}



/* convert the atoms char array into vdw radii information for structure checking
 * uses Bondii radii 
 */

void convert_atom2atom_vdw(char *atom,float *atom_vdw, int num_atoms)
{

	for (int i = 0; i < num_atoms; i++)
	{
		if      (atom[2*i] == 'C' && atom[2*i+1] == ' ')
		atom_vdw[i]=1.7;
		else if (atom[2*i] == 'H' && atom[2*i+1] == ' ')
		atom_vdw[i]=1.1;
		else if (atom[2*i] == 'N' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.55 ;
		else if (atom[2*i] == 'O' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.52;
		else if (atom[2*i] == 'F' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.47;
		else if (atom[2*i] == 'P' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.8;
		else if (atom[2*i] == 'S' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.8;
		else if (atom[2*i] == 'C' && atom[2*i+1] == 'l')
		atom_vdw[i] = 1.75;
		else if (atom[2*i] == 'B' && atom[2*i+1] == 'r')
		atom_vdw[i] = 1.85;
		else if (atom[2*i] == 'I' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.98;
		else if (atom[2*i] == 'B' && atom[2*i+1] == ' ')
		atom_vdw[i] = 1.92;
		else if (atom[2*i] == 'H' && atom[2*i+1] == 'e')
		atom_vdw[i] = 1.40;
		else if (atom[2*i] == 'N' && atom[2*i+1] == 'e')
		atom_vdw[i] = 1.54;
		else if (atom[2*i] == 'K' && atom[2*i+1] == 'r')
		atom_vdw[i] = 2.02;
		else if (atom[2*i] == 'S' && atom[2*i+1] == 'i')
		atom_vdw[i] = 2.10;		
		else
			{printf("***ERROR: atom2atom_vdw: atom not found -> %c%c\n", atom[2*i], atom[2*i+1]);exit(0);}
	//	printf("%d --> %f \n", i, atom_vdw[i]);
	}
}

/* used for checking self-image overlap. this has the information if the 
 * atoms are bonded.
 */
void calc_bond_length(float *bond_length,float *X, float *Y, float *Z, int N)
{
	
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			*(bond_length+i*N+j) = sqrt( (X[i] - X[j])*(X[i] - X[j]) + 
			(Y[i] - Y[j])*(Y[i] - Y[j]) + (Z[i] - Z[j])*(Z[i] - Z[j]) );
		//printf("%f \t" , *(bond_length+i*N+j));	
		}
		//printf("\n");
	}
}

/*checks the distance between two pair of molecules
 */
//approximate
int check_pair_tier2( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *atom_vdw,
				int i,
				int j,
				int N,
				float sr)
{
	
	for (int k = i; k < i + N; k++)
	{
		for (int l = j; l < j + N; l++)
		{
		//	printf("k=%d \t l= %d  \n", k, l);
			if (pdist_appx(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]) <
				sr * ( atom_vdw[k]+atom_vdw[l] )					   )
			{
				//printf("atom1 = %f %f %f \n", X[k], Y[k], Z[k]);
				//printf("atom2 = %f %f %f \n", X[l], Y[l], Z[l]);
				//printf("%f  %f  %f", sr , atom_vdw[k], atom_vdw[l]);
				//printf("failed at %d atom and %d atom with distance %f,
				// expected %f \n", k+1%N, l+1%N,
				// pdist(T,T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]),
				//sr*(atom_vdw[k]+atom_vdw[l]) );
				return 0;
			}
		}
	}
//	printf("yes");
	return 1;
}

//rigorous
int check_pair_tier3( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *atom_vdw,
				int i,
				int j,
				int N,
				float sr)
{
	
	for (int k = i; k < i + N; k++)
	{
		for (int l = j; l < j + N; l++)
		{
		//	printf("k=%d \t l= %d  \n", k, l);
			if (pdist(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]) <
				sr * ( atom_vdw[k]+atom_vdw[l] )					   )
			{
				/*
				printf("atom1 = %f %f %f \n", X[k], Y[k], Z[k]);
				printf("atom2 = %f %f %f \n", X[l], Y[l], Z[l]);
				printf("sr = %f  %f  %f\n", sr , atom_vdw[k], atom_vdw[l]);
				printf("failed at %d atom and %d atom with distance %f expected %f \n", k+1%N, l+1%N,\
				 pdist(T,T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]),\
				sr*(atom_vdw[k]+atom_vdw[l]) );

				float p_dist_app = pdist_appx(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]);
				printf("pdist_appx = %f\n", p_dist_app);
				*/
				return 0;
			}
		}
	}
//	printf("yes");
	return 1;
}

//approximate
int check_self_tier2(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *atom_vdw, int i,int  N,float  sr, float *bond_length)
{

	int j,k, modj,modk;
	float pdist_kj, pdist_kj_2;
	for (j = i; j < i+N; j++)
	{
		modj = j % N;
		for (k = j ; k < i + N; k++)
		{
			modk = k % N;
			pdist_2_appx(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &pdist_kj, &pdist_kj_2);
			
			if (k == j && pdist_kj_2 < sr*(atom_vdw[k]+atom_vdw[j]) )
				return 0;
				
			if (pdist_kj - *(bond_length+modj*N+modk) < -MOL_SIM_TOL)
			{
				if (pdist_kj < sr*(atom_vdw[k]+atom_vdw[j]) )
				{
					///printf("failed at %d atom and %d atom with distance %f, expected %f , bondlength of %f\n"
					//, k+1, j+1, pdist_kj,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					//return 0 ;
				}
			}
			
			if(pdist_kj_2 < sr*(atom_vdw[k]+atom_vdw[j]) )
				return 0;
			
		}
	}
	return 1;
}


//rigourous
int check_self_tier3(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *atom_vdw, int i,int  N,float  sr, float *bond_length)
{

	int j,k, modj,modk;
	float pdist_kj, pdist_kj_2;
	for (j = i; j < i+N; j++)
	{
		modj = j % N;
		for (k = j ; k < i + N; k++)
		{
			modk = k % N;
			pdist_2(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &pdist_kj, &pdist_kj_2);
			
			if (k == j && pdist_kj_2 < sr*(atom_vdw[k]+atom_vdw[j]) )
			{	
					//printf("failed at %d atom and %d atom with distance %f, %f, expected %f ,bondlength of %f\n",
					// k+1, j+1, pdist_kj, pdist_kj_2,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					//float p_dist1, p_dist2;
					//pdist_2_appx(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &p_dist1, &p_dist2);
					//printf("appx = %f, %f \n", p_dist1, p_dist2 );
				
				return 0;
			}	
			if (pdist_kj - *(bond_length+modj*N+modk) < -MOL_SIM_TOL)
			{
				if (pdist_kj < sr*(atom_vdw[k]+atom_vdw[j]) )
				{
					/*
					printf("failed at %d atom and %d atom with distance %f, %f, expected %f ,bondlength of %f\n",
					 k+1, j+1, pdist_kj, pdist_kj_2,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					float p_dist1, p_dist2;
					pdist_2_appx(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &p_dist1, &p_dist2);
					printf("appx = %f, %f \n", p_dist1, p_dist2 );
					*/
					return 0 ;
				}
			}
			
			if(pdist_kj_2 < sr*(atom_vdw[k]+atom_vdw[j]) )
			{
				/*
				   printf("failed at %d atom and %d atom with distance %f, %f, expected %f ,bondlength of %f\n",
					 k+1, j+1, pdist_kj, pdist_kj_2,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					float p_dist1, p_dist2;
					pdist_2_appx(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &p_dist1, &p_dist2);
					printf("appx = %f, %f \n", p_dist1, p_dist2 );
				*/
				return 0;
			}
		}
	}
	return 1;
}

/* Precheck if intermolecular distances are too close.
 * Doesn't take into account of periodic boundary condition, but is 
 * much faster than rigourous checking.
 * Doesn't check molecule with its own periodic image
 */
int fast_screener(crystal xtal, float sr, float *atom_vdw)
{
	//number of atoms in a molecule
	int N = xtal.num_atoms_in_molecule;	
	int m = xtal.Z;	//number of molecules in a unit cell;
	// numberof atom in a unit cell = N*m
	int total_atoms = N*m;	

	float mol_len = 11.6;
	float small_number = 0.1;
	
	for(int i = 0; i < total_atoms; i += N)
	{
		float com1[3];
		compute_molecule_COM( xtal, com1, i);
		for(int j = i + N; j < total_atoms; j += N)
		{
			//check if the molecule COM are far. if they are, dont 
			//bother checking distances
			if (j == i + N)
			{	float com2[3];
				compute_molecule_COM( xtal, com2, j);
				if( sqrt ((com1[0] - com2[0])*(com1[0] - com2[0])+
						  (com1[1] - com2[1])*(com1[1] - com2[1])+
						  (com1[2] - com2[2])*(com1[2] - com2[2]))
					<     (mol_len + small_number)                  )
					continue;
			}
			
			//molecule COM are close than molecule length
			else
			{
				for(int k = i; k < i + N; k++)
				{	
					for (int z = j; z < j + N; z++ )
					{
						if(		(xtal.Xcord[k] - xtal.Xcord[z])*
							   (xtal.Xcord[k] - xtal.Xcord[z])+
							   (xtal.Ycord[k] - xtal.Ycord[z])*
							   (xtal.Ycord[k] - xtal.Ycord[z])+
							   (xtal.Zcord[k] - xtal.Zcord[z])*
							   (xtal.Zcord[k] - xtal.Zcord[z])
							   <sr * (atom_vdw[k]+atom_vdw[z])*
								sr * (atom_vdw[k]+atom_vdw[z])		)
						    return 0;
					}
				}
			}
		}
	}

	return 1;
}


//void main()
int check_structure(crystal random_crystal, float sr)
{
	//float T[3][3]; //lattice_vectors
	int N = random_crystal.num_atoms_in_molecule;	//number of atoms in a molecule
	int m = random_crystal.Z;	//number of molecules in a unit cell; numberof atom in a unit cell = N*m	
	int total_atoms = m*N;
	float *atom_vdw= malloc (N*m*sizeof(float));
	//convert atoms array to atom_vdw distance. this array has vdw informtion of each atom
	
	convert_atom2atom_vdw(random_crystal.atoms, atom_vdw, N*m);
	
	//tier 1 check 
	if( fast_screener(random_crystal, sr, atom_vdw) == 0)
		{ free (atom_vdw);  return 0;}
	//end tier 1 check
	
	
	// tier 2 check
	float *bond_length = malloc (N*N*sizeof(float));
	int final_verdict =1;
	int check_val=1;
	float T_inv[3][3];
	int i,j;
	inverse_mat3b3(T_inv, random_crystal.lattice_vectors);

	//start checking each pair of molecule
	for (i = 0; i < total_atoms; i = i + N)
	{
		for(j= i + N; j < total_atoms; j = j + N)
		{
			check_val = check_pair_tier2(random_crystal.lattice_vectors,T_inv,
			random_crystal.Xcord,random_crystal.Ycord, random_crystal.Zcord,
			atom_vdw, i, j, N, sr);
			
			if (check_val == 0) // 0 if check failed , 1 if it is ok
			{
				i = N*m + 1; //to break out of both loops if the check failed
				j = N*m + 1;
				final_verdict = 0;
			}
		}
	}
	/*
    if(final_verdict == 0)
      printf("failed at tier2 pair check");
	*/
    //check self-image
	 calc_bond_length(bond_length, random_crystal.Xcord, 
		random_crystal.Ycord,  random_crystal.Zcord, N);
    if (final_verdict == 1)
    {
		for (i = 0; i < total_atoms; i = i + N)
		{
			check_val = check_self_tier2(random_crystal.lattice_vectors, T_inv,
				random_crystal.Xcord, random_crystal.Ycord, random_crystal.Zcord,
				atom_vdw, i, N, sr, bond_length);
				
			if (check_val == 0)
			{
				i = total_atoms+ 1;
				//break;
			}
		}
	}

	if (check_val == 0)
	{	free(atom_vdw);
		free(bond_length);
        //printf("failed at tier2 self check\n");
		return final_verdict = 0;
	}


	/*
	if (final_verdict == 1)
	printf("Verdict : Pass \n");
	else
	printf("Verdict : fail \n");
	*/	
	//end of tier-2 check

														//start tier-3 check
	
	if (final_verdict == 1)
	{
		for (i = 0; i < total_atoms; i = i + N)
		{
			for(j= i + N; j < total_atoms; j = j + N)
			{
				check_val = check_pair_tier3(random_crystal.lattice_vectors,T_inv,
				random_crystal.Xcord,random_crystal.Ycord, random_crystal.Zcord,
				atom_vdw, i, j, N, sr);
				
				if (check_val == 0) // 0 if check failed , 1 if it is ok
				{
					i = N*m + 1; //to break out of both loops if the check failed
					j = N*m + 1;
					final_verdict = 0;
				}
			}
		}
	}
    /*
    if (final_verdict == 0)
       printf("failed at tier3 pair check\n"); 
	*/
    if (final_verdict == 1)
    {
		for (i = 0; i < total_atoms; i = i + N)
		{
			check_val = check_self_tier3(random_crystal.lattice_vectors, T_inv,
				random_crystal.Xcord, random_crystal.Ycord, random_crystal.Zcord,
				atom_vdw, i, N, sr, bond_length);
				
			if (check_val == 0)
			{
				i = total_atoms+ 1;
				final_verdict = 0;
			}
		}
	}
	/*													//end tier-3
	if (final_verdict == 0)
       printf("failed at tier3 self check\n") ;
    */
	free(atom_vdw);
	free(bond_length);
	return final_verdict;
} 

/*            			END OF STRUCTURE CHECKS 				*/
/*						BEGIN STRUCTURE CHECKS WTIH VDW DISTANCE CUTOFF  */



int check_structure_with_vdw_matrix(crystal random_crystal,
	float *vdw_matrix,
	int dim1,
	int dim2)
{
	int N = random_crystal.num_atoms_in_molecule;	//number of atoms in a molecule
	int m = random_crystal.Z;						//number of molecules in a unit cell;
													// numberof atom in a unit cell = N*m	
	int total_atoms = m*N;
	
	if (dim1 != dim2)
		{printf("ERROR:matrix not square\n"); exit(0);}

	if (dim1 != total_atoms)
		{printf("matrix size doesnt match total atoms\n"); exit(0);}

	if( fast_screener_vdw(random_crystal, vdw_matrix) == 0)
		{return 0;}
	
	int final_verdict =1;
	int check_val=1;
	float T_inv[3][3];
	int i,j;
	inverse_mat3b3(T_inv, random_crystal.lattice_vectors);

	//tier2 checks
	//start checking each pair of molecule
	for (i = 0; i < total_atoms; i = i + N)
	{
		for(j= i + N; j < total_atoms; j = j + N)
		{
			check_val = check_pair_vdw_tier2(random_crystal.lattice_vectors,
									   T_inv,
									   random_crystal.Xcord,
									   random_crystal.Ycord,
									   random_crystal.Zcord,
									   vdw_matrix,
									   i, 
									   j,
									   N,
									   total_atoms);
			
			//if(i == 0)
			//printf("check val = %d \n", check_val );
			if (check_val == 0) // 0 if check failed , 1 if it is ok
			{
				i = N*m + 1; //to break out of both loops if the check failed
				j = N*m + 1;
				final_verdict = 0;
			}
		}
	}
	
	//check self-image
	float *bond_length = malloc (N*N*sizeof(float));
	 calc_bond_length(bond_length,
					  random_crystal.Xcord,
					  random_crystal.Ycord,
					  random_crystal.Zcord,
					  N);
		
    if (final_verdict == 1)
    {
		for (i = 0; i < total_atoms; i = i + N)
		{
			check_val = check_self_vdw_tier2(random_crystal.lattice_vectors,
									T_inv,
									random_crystal.Xcord,
								   random_crystal.Ycord,
								   random_crystal.Zcord,
								   vdw_matrix,
								   i,
								   N,
								   total_atoms,
								   bond_length);
				
			if (check_val == 0)
			{
				i = total_atoms+ 1;
				//break;
			}
		}
	}
	if (check_val == 0)
		final_verdict = 0;

														//start tier-3 check
	
	if (final_verdict == 1)
	{
		for (i = 0; i < total_atoms; i = i + N)
		{
			for(j= i + N; j < total_atoms; j = j + N)
			{
				check_val = check_pair_vdw_tier3(random_crystal.lattice_vectors,
									   T_inv,
									   random_crystal.Xcord,
									   random_crystal.Ycord,
									   random_crystal.Zcord,
									   vdw_matrix,
									   i, 
									   j,
									   N,
									   total_atoms);
				
				if (check_val == 0) // 0 if check failed , 1 if it is ok
				{
					i = N*m + 1; //to break out of both loops if the check failed
					j = N*m + 1;
					final_verdict = 0;
				}
			}
		}
	} 
	
    if (final_verdict == 1)
    {
		for (i = 0; i < total_atoms; i = i + N)
		{
			check_val = check_self_vdw_tier3(random_crystal.lattice_vectors,
									T_inv,
									random_crystal.Xcord,
								   random_crystal.Ycord,
								   random_crystal.Zcord,
								   vdw_matrix,
								   i,
								   N,
								   total_atoms,
								   bond_length);
				
			if (check_val == 0)
			{
				i = total_atoms+ 1;
				final_verdict = 0;
			}
		}
	}
														//end tier-3

	free(bond_length);

	return final_verdict;
} 

int fast_screener_vdw(crystal xtal, float *vdw_matrix)
{
									//number of atoms in a molecule
	int N = xtal.num_atoms_in_molecule;	
	int m = xtal.Z;					//number of molecules in a unit cell;
								// numberof atom in a unit cell = N*m
	int total_atoms = N*m;	
	float mol_len = 11.6;
	float small_number = 0.1;
	
	for(int i = 0; i < total_atoms; i += N)
	{
		float com1[3];
		compute_molecule_COM( xtal, com1, i);
		for(int j = i + N; j < total_atoms; j += N)
		{
			//check if the molecule COM are far. if they are, dont 
			//bother checking distances
			if (j == i + N)
			{	float com2[3];
				compute_molecule_COM( xtal, com2, j);
				if( sqrt ((com1[0] - com2[0])*(com1[0] - com2[0])+
						  (com1[1] - com2[1])*(com1[1] - com2[1])+
						  (com1[2] - com2[2])*(com1[2] - com2[2]))
					<     (mol_len + small_number)                  )
					continue;
			}
			
			//molecule COM are close than molecule length
			else
			{
				for(int k = i; k < i + N; k++)
				{	
					for (int z = j; z < j + N; z++ )
					{
						if(		(xtal.Xcord[k] - xtal.Xcord[z])*
							   (xtal.Xcord[k] - xtal.Xcord[z])+
							   (xtal.Ycord[k] - xtal.Ycord[z])*
							   (xtal.Ycord[k] - xtal.Ycord[z])+
							   (xtal.Zcord[k] - xtal.Zcord[z])*
							   (xtal.Zcord[k] - xtal.Zcord[z])
							   < *(vdw_matrix + total_atoms*k + z) *
								 *(vdw_matrix + total_atoms*k + z)	)
						    return 0;
					}
				}
			}
		}
	}

	return 1;
}


/*checks the distance between two pair of molecules
 */
int check_pair_vdw_tier2( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms)
{
	
	for (int k = i; k < i + N; k++)
	{
		for (int l = j; l < j + N; l++)
		{
		//	printf("k=%d \t l= %d  \n", k, l);
			if (pdist_appx(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]) <
				 *(vdw_matrix + k*total_atoms + l) 					)
			{
				//printf("atom1 = %f %f %f \n", X[k], Y[k], Z[k]);
				//printf("atom2 = %f %f %f \n", X[l], Y[l], Z[l]);
				//printf("%f  %f  %f", sr , atom_vdw[k], atom_vdw[l]);
				//printf("failed at %d atom and %d atom with distance %f,
				// expected %f \n", k+1%N, l+1%N,
				// pdist(T,T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]),
				//sr*(atom_vdw[k]+atom_vdw[l]) );
				return 0;
			}
		}
	}
//	printf("yes");
	return 1;
}

/*check with periodic image of one molecule: approximate but faster
*/
int check_self_vdw_tier2(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length)
{

	int j,k, modj,modk;
	float pdist_kj, pdist_kj_2;
	for (j = i; j < i+N; j++)
	{
		modj = j % N;
		for (k = j ; k < i + N; k++)
		{
			modk = k % N;
			pdist_2_appx(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &pdist_kj, &pdist_kj_2);
			
			if (k == j && pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
				
			if (pdist_kj - *(bond_length+modj*N+modk) < -0.02)
			{
				if (pdist_kj < *(vdw_matrix + k*total_atoms + j) )
				{
					//printf("failed at %d atom and %d atom with distance %f, expected %f , bondlength of %f\n",
					// k+1, j+1, pdist_kj,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					return 0 ;
				}
			}
			
			if(pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
			
		}
	}
	return 1;
}


int check_pair_vdw_tier3( float T[3][3],
				float T_inv[3][3],
				float *X,
				float *Y,
				float *Z,
				float *vdw_matrix,
				int i,
				int j,
				int N,
				int total_atoms)
{
	
	for (int k = i; k < i + N; k++)
	{
		for (int l = j; l < j + N; l++)
		{
		//	printf("k=%d \t l= %d  \n", k, l);
			if (pdist(T, T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]) <
				 *(vdw_matrix + k*total_atoms + l) 					)
			{
				//printf("atom1 = %f %f %f \n", X[k], Y[k], Z[k]);
				//printf("atom2 = %f %f %f \n", X[l], Y[l], Z[l]);
				//printf("%f  %f  %f", sr , atom_vdw[k], atom_vdw[l]);
				//printf("failed at %d atom and %d atom with distance %f,
				// expected %f \n", k+1%N, l+1%N,
				// pdist(T,T_inv, X[k], Y[k], Z[k], X[l], Y[l], Z[l]),
				//sr*(atom_vdw[k]+atom_vdw[l]) );
				return 0;
			}
		}
	}
//	printf("yes");
	return 1;
}

/*check with periodic image of one molecule: rigourous and slow
*/
int check_self_vdw_tier3(float T[3][3], float T_inv[3][3],float *X,float *Y,
	float *Z, float *vdw_matrix, int i,int  N,int total_atoms, float *bond_length)
{

	int j,k, modj,modk;
	float pdist_kj, pdist_kj_2;
	for (j = i; j < i+N; j++)
	{
		modj = j % N;
		for (k = j ; k < i + N; k++)
		{
			modk = k % N;
			pdist_2(T,T_inv, X[k], Y[k], Z[k], X[j], Y[j], Z[j], &pdist_kj, &pdist_kj_2);
			
			if (k == j && pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
				
			if (pdist_kj - *(bond_length+modj*N+modk) < -0.02)
			{
				if (pdist_kj < *(vdw_matrix + k*total_atoms + j) )
				{
					//printf("failed at %d atom and %d atom with distance %f, expected %f , bondlength of %f\n", k+1, j+1, pdist_kj,sr*(atom_vdw[k]+atom_vdw[j]), *(bond_length+modj*N+modk) );
					return 0 ;
				}
			}
			
			if(pdist_kj_2 < *(vdw_matrix + k*total_atoms + j) )
				return 0;
			
		}
	}
	return 1;
}
