#include <stdio.h>
#include <math.h>
#include <string.h>
#include "read_input.h"
#include "spg_generation.h"
#include "crystal_utils.h"
#include "algebra.h"
#include "check_structure.h"
#include "spglib.h"
#include "combinatorics.h"

extern float TOL;

float get_crystal_volume(crystal *xtal)
{
	return xtal->lattice_vectors[0][0]*
		   xtal->lattice_vectors[1][1]*
		   xtal->lattice_vectors[2][2];
}

void print_crystal(crystal* xtal)
{
	
	
	 int N = xtal->num_atoms_in_molecule;
	 int m = xtal->Z;
	 
	 printf("#Z = %d \n", xtal->Z);
	 printf("#napm = %d \n", xtal->num_atoms_in_molecule);
	
	for(int i = 0; i < 3; i++)
	{
		printf("lattice_vector %12f %12f %12f \n",
			xtal->lattice_vectors[i][0], xtal->lattice_vectors[i][1],
			xtal->lattice_vectors[i][2]);
	}
	
	for(int i = 0; i < N*m; i++)
	{
		printf("atom %12f %12f %12f %4c \n", xtal->Xcord[i],
			xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[2*i]);
	}
}

void print_crystal2file(crystal* xtal, FILE* out_file)
{
	
	static int counter = 1;
	int N = xtal->num_atoms_in_molecule;
	int m = xtal->Z;
	 
	fprintf(out_file, "####### BEGIN STRUCTURE #######\n");
	fprintf(out_file, "#structure_number = %d\n", counter);
	fprintf(out_file, "#Z = %d\n", xtal->Z);
	fprintf(out_file, "#number_of_atoms_in_molecule = %d \n",
	 xtal->num_atoms_in_molecule);
    
    fprintf(out_file, "#unit_cell_volume = %f cubic Angstrom\n", get_crystal_volume(xtal));
	fprintf(out_file, "#attempted_spacegroup = %d\n", xtal->spg);
	
	char letter = spg_positions[xtal->spg - 1].wyckoff_letter[xtal->wyckoff_position];
	fprintf(out_file, "#attempted_wyckoff_position = %d%c\n", xtal->Z, letter);

	char site_symm[6];
	strcpy(site_symm, spg_positions[xtal->spg - 1].site_symmetry[xtal->wyckoff_position]);
	fprintf(out_file, "#site_symmetry_group = %s\n",
	 site_symm);
	
	int spglib_spg = detect_spg_using_spglib(xtal);
	fprintf(out_file, "#SPGLIB_detected_spacegroup = %d\n", spglib_spg);
	
	fprintf(out_file, "#\"All distances in Angstroms and using Cartesian coordinate system\"\n");
		
	for(int i = 0; i < 3; i++)
	{
		fprintf(out_file,"lattice_vector %12f %12f %12f \n",
			xtal->lattice_vectors[i][0], xtal->lattice_vectors[i][1],
			xtal->lattice_vectors[i][2]);
	}
	
	for(int i = 0; i < N*m; i++)
	{
		fprintf(out_file,"atom %12f %12f %12f  %c%c \n", xtal->Xcord[i],
			xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[2*i], xtal->atoms[2*i+1]);
	}
	fprintf(out_file, "#######  END  STRUCTURE #######\n\n");
	
	fflush(out_file);
	counter++;
}

void print_crystal_fractional(crystal* xtal)
{
	
	
	 int N = xtal->num_atoms_in_molecule;
	 int m = xtal->Z;
	 
	 printf("#Z = %d \n", xtal->Z);
	 printf("#napm = %d \n", xtal->num_atoms_in_molecule);
	
	for(int i = 0; i < 3; i++)
	{
		printf("lattice_vector %f %f %f \n",
			xtal->lattice_vectors[i][0], xtal->lattice_vectors[i][1],
			xtal->lattice_vectors[i][2]);
	}
	
	for(int i = 0; i < N*m; i++)
	{
		printf("atom_frac %f %f %f %c%c \n", xtal->Xcord[i],
			xtal->Ycord[i],  xtal->Zcord[i],  xtal->atoms[i], xtal->atoms[i+1]);
	}
}


void compute_molecule_COM(crystal xtal, float com[3], int i)
{	/* computes the COM of a molecule (whose first atom is i) in
	a crystal. 
	*/
	com[0] = 0;
	com[1] = 0;
	com[2] = 0;
	
	int N = xtal.num_atoms_in_molecule;
	for(int j = i; j < i + N ; j++)
	{
		com[0] += xtal.Xcord[j]; 
		com[1] += xtal.Ycord[j];
		com[2] += xtal.Zcord[j];
	}
	com[0] /= N;
	com[1] /= N;
	com[2] /= N;
	
	//print_vec3(com);
	
}

void bring_all_molecules_to_first_cell(crystal *xtal)
{
	int N = xtal->num_atoms_in_molecule;
	int Z = xtal->Z;
	
	convert_xtal_to_fractional(xtal);
	//convert to fractional
	for(int i = 0; i < N*Z; i+=N)
	{
		float com[3] = {0,0,0};
		compute_molecule_COM(*xtal,com,i);
		int integer_part1 = (int)com[0];
		int integer_part2 = (int)com[1];
		int integer_part3 = (int)com[2];
		if(com[0] < 0)
		{integer_part1 -= 1;}
		if(com[1] < 0)
		{integer_part2 -= 1;}
		if(com[2] < 0)
		{integer_part3 -= 1;}

		for(int j = i; j < i + N; j++)
		{
			xtal->Xcord[j] -= integer_part1;
			xtal->Ycord[j] -= integer_part2;
			xtal->Zcord[j] -= integer_part3;
		}
		
	}
	
	convert_xtal_to_cartesian(xtal);
}

void convert_xtal_to_cartesian(crystal *xtal)
{
	float lattice_vectors_transpose[3][3];
	copy_mat3b3_mat3b3(lattice_vectors_transpose,
					   xtal->lattice_vectors);
	mat3b3_transpose(lattice_vectors_transpose,
					 lattice_vectors_transpose);
	int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;
	
	for(int i = 0; i < total_atoms; i++)
	{
		float temp[3] = {xtal->Xcord[i],
						 xtal->Ycord[i],
						 xtal->Zcord[i] };
		vector3_mat3b3_multiply(lattice_vectors_transpose, temp, temp );
		xtal->Xcord[i] = temp[0];	
		xtal->Ycord[i] = temp[1];
		xtal->Zcord[i] = temp[2];
	}
}

void remove_close_molecules(crystal* xtal)
{
	/* function to remove close molecules in a cell. Originally written
	 * to remove duplicate overlapping molecules in a cell created by 
	 * performing symmetry operations when the molecule is at an 
	 * inversion center.
	 * remove duplicates if COM is close
	 */
	

	int N = xtal->num_atoms_in_molecule;
	int m = xtal->Z;
	int total_atoms = m * N;
	//to count the number of molecules after duplicate removal
	int molecule_counter = 0;
	 
	 //for computing pdist
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);

	
	 //temp variables for storing unique molecules
	float tempx[total_atoms];
	float tempy[total_atoms];
	float tempz[total_atoms];
	char tempatom[2*total_atoms];
	

	//printf("total atoms = %d \n", total_atoms);
	 
	for(int i = 0; i < total_atoms; i = i + N)
	{
		int same = 0;
		float com1[3] = {0,0,0};
		compute_molecule_COM(*xtal, com1, i);
		//print_vec3(com1);
		for(int j = i + N; j < total_atoms; j = j + N)
		{
			float com2[3] = {0,0,0};
			compute_molecule_COM(*xtal, com2, j);
			/*if(fabs(com1[0] - com2[0]) < tol && fabs(com1[1] - com2[1]) < tol 
				&& fabs(com1[2] - com2[2]) < tol )
			*/
			if ( pdist_appx(lattice_vectors, inv_lattice_vectors, com1[0],
					com1[1], com1[2], com2[0], com2[1], com2[2]) < TOL )	
			{
				same = 1;
				//printf("i = %d j = %d \n ",i/N+1 ,j/N+1);
				break;
			}
		}
		if (same == 0) //unique molecule
		{
			for(int z = 0; z < N; z++)
			{
				tempx[molecule_counter*N + z] = xtal->Xcord[i+z];
				tempy[molecule_counter*N + z] = xtal->Ycord[i+z];
				tempz[molecule_counter*N + z] = xtal->Zcord[i+z]; 
				tempatom[2*(molecule_counter*N +z)] = xtal->atoms[2*(i+z)];
                tempatom[2*(molecule_counter*N +z)+1] = xtal->atoms[2*(i+z)+1];	
			}
			
			molecule_counter++;
		}	
	}
	
	//printf("number of uniqu mol = %d \n", molecule_counter);
	for(int i = 0; i < molecule_counter*N; i++)
	{
		xtal->Xcord[i] = tempx[i];
		xtal->Ycord[i] = tempy[i];
		xtal->Zcord[i] = tempz[i];
		xtal->atoms[2*i] = tempatom[2*i];
        xtal->atoms[2*i+1] = tempatom[2*i+1];
	}
	xtal->Z = molecule_counter;
}

void convert_xtal_to_fractional(crystal *xtal)
{
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);
	mat3b3_transpose(inv_lattice_vectors, inv_lattice_vectors);
	
	int total_atoms = xtal->Z *xtal->num_atoms_in_molecule;
	
	for(int i = 0; i < total_atoms; i++)
	{
		float temp[3] = {xtal->Xcord[i], xtal->Ycord[i], xtal->Zcord[i] };
		vector3_mat3b3_multiply(inv_lattice_vectors, temp, temp );
		xtal->Xcord[i] = temp[0];	
		xtal->Ycord[i] = temp[1];
		xtal->Zcord[i] = temp[2];
	}
}

void combine_close_molecules(crystal* xtal)
{
	/* function to combine close molecules in a cell. 
	 * remove duplicates if COM is close.
	 * average the cordinates of the atom with the nearest atom of the
	 * neighbouring molecule.
	 */

	int N = xtal->num_atoms_in_molecule;
	int Z = xtal->Z;
	int total_atoms = Z * N;
	//to count the number of molecules after duplicate removal
	int molecule_counter = 0;
	int same = 0;
	//for partitioning molecules into different clusters
	int partition_list[Z];
	int averaging_list[Z];
	for(int i =0; i < Z; i++) //initialise
	{	
		partition_list[i] = -1;
		averaging_list[i] = -1;
	}
	 
	 //for computing pdist
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);

	
	 //temp variables for storing unique molecules
	float tempx[total_atoms];
	float tempy[total_atoms];
	float tempz[total_atoms];
	char tempatom[2*total_atoms];
	
	//printf("total atoms = %d \n", total_atoms);
	 
	 //outer loop over molecules in xtal
	for(int i = 0; i < total_atoms; i = i + N)
	{
		if(partition_list[i/N] != -1)
			continue;
			
		float com1[3] = {0,0,0};
		compute_molecule_COM(*xtal, com1, i);
		//print_vec3(com1);
		//inner loop over molecules in xtal
		for(int j = i + N; j < total_atoms; j = j + N)
		{
			if(partition_list[j/N] != -1)
				continue;
			
			float com2[3] = {0,0,0};
			compute_molecule_COM(*xtal, com2, j);
			//if two molecules are overlapping
			if ( pdist_appx(lattice_vectors,
					   inv_lattice_vectors,
					   com1[0],
					   com1[1],
					   com1[2],
					   com2[0],
					   com2[1],
					   com2[2]) < TOL )	
			{
				int index1 = j / N;
				int index2 = i / N;
				partition_list[index1] = j;
				averaging_list[index1] = j;
				partition_list[index2] = i;
				averaging_list[index2] = i;
				same = 1;
				
				/*if a molecule is at corner and the other is at other
				 * corner, move the molecule such that they overlap
				 */
				 if (cart_dist(com1, com2) > TOL) 
				 {
					float displace[3];
					vector3_subtract(com1, com2, displace);
					for(int k = j; k < j + N; k++)
					{
						xtal->Xcord[k] += displace[0];
						xtal->Ycord[k] += displace[1];
						xtal->Zcord[k] += displace[2];  
					}
				}
				
			}	
		}
		
		//average overlapping molecules
		if (same)
		{
			average_positions_using_list(tempx,
									 tempy,
									 tempz,
									 tempatom,
									 molecule_counter,
									 N,
									 Z,
									 xtal,
									 averaging_list);
			//reset averaging list
			for(int i = 0; i < Z; i++ )
				averaging_list[i] = -1;
			
		}
		else
		{
			for(int z = 0; z < N; z++)
			{
				tempx[molecule_counter*N + z] = xtal->Xcord[i+z];
				tempy[molecule_counter*N + z] = xtal->Ycord[i+z];
				tempz[molecule_counter*N + z] = xtal->Zcord[i+z]; 
				tempatom[molecule_counter*N +z] = xtal->atoms[i+z];
                tempatom[2*(molecule_counter*N +z)+1] = xtal->atoms[2*(i+z)+1];	
			}
		}	
		
		molecule_counter++;
	}
	//printf("number of uniqu mol = %d \n", molecule_counter);
	//place averaged molecules back into the xtal structure and update Z
	for(int i = 0; i < molecule_counter*N; i++)
	{
		xtal->Xcord[i] = tempx[i];
		xtal->Ycord[i] = tempy[i];
		xtal->Zcord[i] = tempz[i];
		xtal->atoms[2*i] = tempatom[2*i];
        xtal->atoms[2*i+1] = tempatom[2*i+1];
	}
	xtal->Z = molecule_counter;
}

void average_positions_using_list(float 	*tempx,
								  float 	*tempy,
								  float 	*tempz,
								  char	 	*tempatom,
								  int   	 molecule_counter,
								  int		 N,
								  int	 	 Z,
								  crystal 	*xtal,
								  int 		*averaging_list)
{
	//modify averaging list by moving all elements to front
	//len is now the new length of the list
	int len = 0;
	for(int i = 0; i < Z ; i++ )
	{
		if (averaging_list[i] != -1)
		{
			averaging_list[len] = averaging_list[i];
			len++;
		} 
	}	
	//traverse through average list
	for (int i = 1; i < len; i++)
	{
		int first_mol = averaging_list[0];
		int second_mol = averaging_list[i];
		//loop through atoms in first mol
		for (int j = 0; j < N; j++)
		{

			float min = 0;
			int mink = 0;
			float p1[3] = 
			{ xtal->Xcord[first_mol+j],
			  xtal->Ycord[first_mol+j],
			  xtal->Zcord[first_mol+j] };
			
			//add first mol to tempx and tempatoms
            if (i == 1)
            {
			    tempx[molecule_counter*N+j] = p1[0];
			    tempy[molecule_counter*N+j] = p1[1];
			    tempz[molecule_counter*N+j] = p1[2];
			    tempatom[2*(molecule_counter*N+j)] = xtal->atoms[2*(first_mol+j)];
                tempatom[2*(molecule_counter*N+j)+1] = xtal->atoms[2*(first_mol+j)+1];
            }            

			for(int k = 0; k < N; k++)
			{
				
				float p2[3] = 
				{ xtal->Xcord[second_mol + k],
				  xtal->Ycord[second_mol + k],
				  xtal->Zcord[second_mol + k] };
				
				float temp_min =  (p1[0]-p2[0])*(p1[0]-p2[0]) +
								  (p1[1]-p2[1])*(p1[1]-p2[1]) +
								  (p1[2]-p2[2])*(p1[2]-p2[2]);
				
				if(temp_min < min || k == 0)
				{
					min = temp_min;
					mink = k;
				}
			}
			
			tempx[molecule_counter*N+j] += xtal->Xcord[second_mol+mink];
			tempy[molecule_counter*N+j] += xtal->Ycord[second_mol+mink];
			tempz[molecule_counter*N+j] += xtal->Zcord[second_mol+mink];
			
		}
	}
	
	//take average
	for (int i = 0; i < N; i++)
	{
		tempx[molecule_counter*N+i] /= len;
		tempy[molecule_counter*N+i] /= len;
		tempz[molecule_counter*N+i] /= len;
	}
}

/*
void remove_close_molecules_mirror_plane(crystal *xtal)
{
	function to remove close molecules in mirror plane.
	 * COM closeness is checked first
	 * normal to mirror is calculated from molecular orientation
	 * one molecule is deleted,other is rotated so that axis coincides 
	 * with the normal
	 
	 
	int N = xtal->num_atoms_in_molecule;
	int m = xtal->Z;
	int total_atoms = m * N;
	//to count the number of molecules after duplicate removal
	int molecule_counter = 0;
	 
	 //for computing pdist
	float lattice_vectors[3][3];
	float inv_lattice_vectors[3][3];
	copy_mat3b3_mat3b3(lattice_vectors, xtal->lattice_vectors);
	inverse_mat3b3(inv_lattice_vectors, lattice_vectors);

	
	 //temp variables for storing unique molecules
	float tempx[total_atoms];
	float tempy[total_atoms];
	float tempz[total_atoms];
	char tempatom[total_atoms];
	
	for(int i = 0; i < total_atoms; i = i + N)
	{
		int same = 0;
		int j;
		float com1[3] = {0,0,0};
		compute_molecule_COM(*xtal, com1, i);
		//print_vec3(com1);
		for(j = i + N; j < total_atoms; j = j + N)
		{
			float com2[3] = {0,0,0};
			compute_molecule_COM(*xtal, com2, j);
			if(fabs(com1[0] - com2[0]) < tol && fabs(com1[1] - com2[1]) < tol 
				&& fabs(com1[2] - com2[2]) < tol )
			
			if ( pdist(lattice_vectors, inv_lattice_vectors, com1[0],
					com1[1], com1[2], com2[0], com2[1], com2[2]) < tol )	
			{
				same = 1;
				printf("identical molecules i = %d j = %d detected \n",i/N+1 ,j/N+1);
				break;
			}
		}
		
		if(same == 1)
		{
			
			//at this point i/N+1 and j/N+1 molecules are equivalent
			int mol_i = i/N;
			int mol_j = j/N;
			
			
		
		
		}

	}
	
	printf("number of uniq mol = %d \n", molecule_counter);
	for(int i = 0; i < molecule_counter*N; i++)
	{
		xtal->Xcord[i] = tempx[i];
		xtal->Ycord[i] = tempy[i];
		xtal->Zcord[i] = tempz[i];
		xtal->atoms[i] = tempatom[i];
	}
	
	if(molecule_counter != m/2)
	{
			printf("**WARNING: all equivalent molecules not detected");
			exit(0);
	}
	xtal->Z = molecule_counter; 
}
*/

							/*spglib check*/
int detect_spg_using_spglib(crystal* xtal)
{
	float tol = 1e-4;  
	//print_crystal(xtal); 
	convert_xtal_to_fractional(xtal);
	//variable declarations				
	int num_atoms_in_cell = xtal->Z * xtal-> num_atoms_in_molecule ;
	int types[num_atoms_in_cell];
	double positions[num_atoms_in_cell][3];
	char atom[num_atoms_in_cell*2];
	char symbol[21];
	int num_spg;
	double lattice_vector[3][3];
	//printf("total atoms = %d \n", num_atoms_in_cell);
	for(int i = 0; i < num_atoms_in_cell; i++)
	{
		positions[i][0] = xtal->Xcord[i]; 
		positions[i][1] = xtal->Ycord[i];
		positions[i][2] = xtal->Zcord[i];
		atom[2*i] = xtal->atoms[2*i];
		atom[2*i+1] = xtal->atoms[2*i+1];
		
		// if (xtal->atoms[2*i] == 'C' && xtal->atoms[2*i+1] == ' ' )
		// types[i] = 12;
		// else if (xtal->atoms[2*i]== 'H' && xtal->atoms[2*i+1] == ' ')
		// types[i] = 1;
		// else if (xtal->atoms[2*i] == 'N' && xtal->atoms[2*i+1] == ' ')
		// types[i] = 7;
		// else if (xtal->atoms[2*i] == 'O' && xtal->atoms[2*i+1] == ' ')
		// types[i] = 8;
		// else
		// {printf("error atom not found : %c\n", xtal->atoms[2*i]);exit(0);}

		if      (atom[2*i] == 'C' && atom[2*i+1] == ' ')
		types[i] = 6; 
		else if (atom[2*i] == 'H' && atom[2*i+1] == ' ')
		types[i] = 1;
		else if (atom[2*i] == 'N' && atom[2*i+1] == ' ')
		types[i] = 7;
		else if (atom[2*i] == 'O' && atom[2*i+1] == ' ')
		types[i] = 8;
		else if (atom[2*i] == 'F' && atom[2*i+1] == ' ')
		types[i] = 9;
		else if (atom[2*i] == 'P' && atom[2*i+1] == ' ')
		types[i] = 15;
		else if (atom[2*i] == 'S' && atom[2*i+1] == ' ')
		types[i] = 16;
		else if (atom[2*i] == 'C' && atom[2*i+1] == 'l')
		types[i] = 17;
		else if (atom[2*i] == 'B' && atom[2*i+1] == 'r')
		types[i] = 35;
		else if (atom[2*i] == 'I' && atom[2*i+1] == ' ')
		types[i] = 53;
		else if (atom[2*i] == 'B' && atom[2*i+1] == ' ')
		types[i] = 5;
		else if (atom[2*i] == 'H' && atom[2*i+1] == 'e')
		types[i] = 2;
		else if (atom[2*i] == 'N' && atom[2*i+1] == 'e')
		types[i] = 10;
		else if (atom[2*i] == 'K' && atom[2*i+1] == 'r')
		types[i] = 36;
		else if (atom[2*i] == 'S' && atom[2*i+1] == 'i')
		types[i] = 14;		
		else
		{printf("***ERROR: spglib detector: atom not found -> %c%c\n", atom[2*i], atom[2*i+1]);exit(0);}
		
		
	}
	for(int i = 0; i < 3; i++)
		for(int j = 0; j < 3; j++ )
		{
			//note transpose
			lattice_vector[i][j] = xtal->lattice_vectors[j][i];
		}

	num_spg = spg_get_international(symbol,
									lattice_vector, 
									positions,
									types,
									num_atoms_in_cell,
									tol);
	//SpglibError error;
	//error = spg_get_error_code();
	//printf("#SPGLIB says %s\n", spg_get_error_message(error));
	//printf("#SPGLIB detected space group = %d\n",
	//											num_spg);
	convert_xtal_to_cartesian(xtal);
	return num_spg;

}

void allocate_xtal(crystal* xtal, int Z, int N)
{
	//xtal = (crystal*)malloc(sizeof(crystal));
	xtal->Xcord = (float *)malloc(Z*N*sizeof(float));
	xtal->Ycord = (float *)malloc(Z*N*sizeof(float));
	xtal->Zcord = (float *)malloc(Z*N*sizeof(float));
	xtal->atoms = (char *)malloc(Z*N*2*sizeof(char));
}

void free_xtal(crystal* xtal)
{
	free(xtal->atoms);
	free(xtal->Xcord);
	free(xtal->Ycord);
	free(xtal->Zcord);
	free(xtal);
}

void copy_xtal(crystal* xtal1, crystal* xtal2)
{
	xtal1->lattice_vectors[0][0] = xtal2->lattice_vectors[0][0];
	xtal1->lattice_vectors[0][1] = xtal2->lattice_vectors[0][1];
	xtal1->lattice_vectors[0][2] = xtal2->lattice_vectors[0][2];
	xtal1->lattice_vectors[1][0] = xtal2->lattice_vectors[1][0];
	xtal1->lattice_vectors[1][1] = xtal2->lattice_vectors[1][1];
	xtal1->lattice_vectors[1][2] = xtal2->lattice_vectors[1][2];
	xtal1->lattice_vectors[2][0] = xtal2->lattice_vectors[2][0];
	xtal1->lattice_vectors[2][1] = xtal2->lattice_vectors[2][1];
	xtal1->lattice_vectors[2][2] = xtal2->lattice_vectors[2][2];
	
	int N = xtal2->num_atoms_in_molecule;
	xtal1->num_atoms_in_molecule = N;
	
	int Z = xtal2->Z;
	xtal1->Z = Z;
	
	xtal1->spg = xtal2->spg;
	
	for(int i = 0; i < N*Z; i++)
	{
		xtal1->Xcord[i] = xtal2->Xcord[i];
		xtal1->Ycord[i] = xtal2->Ycord[i];
		xtal1->Zcord[i] = xtal2->Zcord[i];
		
		xtal1->atoms[2*i] = xtal2->atoms[2*i];
		xtal1->atoms[2*i+1] = xtal2->atoms[2*i+1];
	}
		
}


