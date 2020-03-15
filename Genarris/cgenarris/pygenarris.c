#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "read_input.h"
#include "spg_generation.h"
#include "pygenarris.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"
#include "spglib.h"

#define ZMAX 192
#define VOL_ATTEMPT  100000

unsigned int *seed;
unsigned int *seed2;
#pragma omp threadprivate(seed)
#pragma omp threadprivate(seed2)
extern float TOL;

/*

void generate_molecular_crystals(char *filename,int num_structures, int Z,
	double volume_mean1, double volume_std1, double sr1, double tol1, 
	int max_attempts)
{

    float Zp_max=1;
    float volume_mean = volume_mean1;
    float volume_std = volume_std1;
    float sr = sr1;
    TOL = tol1;

    printf("printing inputs...\n");
    printf("num_of_structures_per_spg = %d, Z=%d volume_mean = %f, volume_std = %f, Zp_max= % f, sr = %f, tol = %f \n",
		num_structures,Z, volume_mean,volume_std, Zp_max, sr, TOL);
    // exit(0);
	//random number seeding
	//random number seeding
	srand((unsigned int)time(NULL));
	int seed_shift = rand()% 1061 + 7 ;
	
	int counter = 0 ;
	int spg_rand = 0;
	int stop_flag = 0;

	FILE *out_file;
	out_file = fopen(filename,"w");
	if(!out_file)
	{
		printf("***ERROR: cannot create geometry.out \n");
		exit(0);
	}

	
	#pragma omp parallel shared(counter, spg_rand, seed_shift, out_file, stop_flag) //default(none)
																	{
	seed = (int*)malloc(sizeof(int));								
	seed2 = (unsigned int*)malloc(sizeof(unsigned int));
	int thread_num = omp_get_thread_num() + 1;
	int total_threads = omp_get_num_threads();
	*seed += thread_num*7 + seed_shift*13;
	*seed2 = thread_num*17 + seed_shift*11;
	init_genrand(*seed);
	
	COMPATIBLE_SPG compatible_spg[230]; 
	int num_compatible_spg = 0;
	int num_axes;
	float *mol_axes;
		
	molecule *mol = (molecule*)malloc(sizeof(molecule));
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));
	float volume;
	int spg;
	//read input from file, read molecular structure from geometry,Z, Zp
    #pragma omp critical
	{
        read_geometry(mol);

        //printf("%d %d %f %f %f\n", num_structures, Z, sr, volume_mean, volume_std);
    }

    //recenter molecule to origin
	recenter_molecule(mol);
	
	#pragma omp barrier
	{}
   
    //inititalise
    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
	int N = mol->num_of_atoms;
	allocate_xtal(random_crystal, ZMAX, N);
	random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
	
	if(thread_num == 1)
	{
		printf("Detecting compatible spacegroups"
			   " from molecule symmetries...\n");
	}
	
	find_compatible_spg_positions(mol,
								  Z,
								  compatible_spg,
								  &num_compatible_spg,
								  &mol_axes,
								  &num_axes,
								  thread_num);
			
	if(thread_num == 1)
	{
		printf("\nTotal compatible space groups = %d\n", num_compatible_spg);
		printf("Number of mol axes = %d \n", num_axes);
		printf("Starting generator ...\n\n");
		sleep(1);
	}
	
	#pragma omp barrier
	{}
	
	//find allowed space groups for general position (deprecated)
	//int allowed_spg[230];
	//int num_allowed_spg = 0;
	//find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
	//printf("allowed %d \n", num_compatible_spg);
	
	while( spg_rand < num_compatible_spg)
	{
		//conter counts the number of structures generated
		spg = compatible_spg[spg_rand].spg;
		while( counter < num_structures)
		{
			int verdict = 0;
			int i = 0; 
			// attempts for an spg
			for(; i < max_attempts/total_threads; i++)
			{
				if(counter >= num_structures)
					break;
				
				if(stop_flag)
					break;
					
				//if successful
				if (verdict == 1 )
				{
					#pragma omp critical
					if(counter < num_structures)
					{
						random_crystal->spg = spg;
						//print_crystal(random_crystal);
						print_crystal2file(random_crystal, out_file);
						counter++;
						//printf("#All Structure checks passed \n");
						printf("#thread %d:Generation Successful after %d attempts\n",
																		 thread_num,
																		 i*total_threads);
									
						printf("#Structure number = %d  spg = %d vol = %f\n", 
																	counter, 
																	spg, volume);
						
						int spglib_spg = detect_spg_using_spglib(random_crystal);
						printf("#SPGLIB detected space group = %d\n",
												                     spglib_spg);
						fflush(stdout);
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
					}
									
					break;
				}
				// else generate again
				int result = generate_crystal(random_crystal,
								 mol,
								 volume,
								 Z,
								 Zp_max,
								 spg,
								 compatible_spg,
								 num_compatible_spg,
								 spg_rand,
								 mol_axes,
								 num_axes);
				
				//alignment failure
				if(!result)
					continue;
				
				//print_crystal(random_crystal); //exit(0);		    
				verdict = check_structure(*random_crystal, sr);    
				
				if(i % VOL_ATTEMPT == 0 && i != 0)
				{
					do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
					if(thread_num == 1)
						printf("#thread 1: completed %d attempts\n", i*total_threads);
					*seed = *seed2 + thread_num*rand_r(seed2);
					fflush(stdout);
				}
				
			}//end of attempt loop	
			
			//if max limit is reached or if some thread hit the limit
			if (i >= max_attempts/total_threads || stop_flag == 1)
			{	
				//stop other threads
				#pragma omp critical
					{stop_flag = 1;}
				#pragma omp barrier
				{}
				do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				if(thread_num == 1)
				{	
					printf("**WARNING: generation failed for spg = %d \
	after %d attempts \n\n", spg, max_attempts);
					printf("#counter reset...\n\n");
					counter = num_structures + 1;
					print_crystal(random_crystal);
					fflush(stdout);
				}
				#pragma omp barrier
				{}
			}	
			
		}//end of numof structures llop
	
	#pragma omp barrier
	{}
	if(thread_num == 1 )
	{	
         counter = 0;
         spg_rand++;
         stop_flag = 0;
         printf("#counter reset...\n\n");
         fflush(stdout);
         //sleep(1);
         //usleep(90000);
    } 
	#pragma omp barrier
	{}
	
	}//while loop //end of spg loop
	
		//free memory
																	   }
	fclose(out_file);
	
	printf("Generation completed!\n");
//	free(mol_axes);
}


//Generator with van der waal cutoff distance matrix
void generate_molecular_crystals_with_vdw_cutoff_matrix(char *filename,
	int seedstate,
	float *vdw_matrix,
	int dim1,
	int dim2,
	int num_structures,
	int Z,
	double volume_mean1,
	double volume_std1,
	double tol1, 
	int max_attempts)
{

    float Zp_max=1;
    float volume_mean = volume_mean1;
    float volume_std = volume_std1;
    TOL = tol1;

    if (dim1 != dim2)
	{printf("***ERROR:vdw cutoff matrix is not square***\n"); exit(0);}

    printf("printing inputs...\n");
    printf("num_of_structures_per_spg = %d, Z = %d, volume_mean = %f, volume_std = %f zp_max = %f ,tol= %f\n",
		num_structures,Z, volume_mean, volume_std, Zp_max, TOL);
    printf("The vdw distance cutoff matrix is :\n");
    for (int i = 0; i < dim1; i++)
    {
    	for (int j= 0; j < dim2; j++)
    	{
    		printf("%f ", *(vdw_matrix + i*dim1 + j) );
    	}
    	printf("\n");
    }

    // exit(0);
	//random number seeding
	//random number seeding
	srand((unsigned int)time(NULL));
	int seed_shift = rand() % 1061 + 7*seedstate;
	
	int counter = 0 ;
	int spg_rand = 0;
	int stop_flag = 0;

	FILE *out_file;
	out_file = fopen(filename,"w");
	if(!out_file)
	{
		printf("***ERROR: cannot create geometry.out \n");
		exit(0);
	}

	
	#pragma omp parallel shared(counter, spg_rand, seed_shift, out_file, stop_flag) //default(none)
																	{
	seed = (int*)malloc(sizeof(int));								
	seed2 = (unsigned int*)malloc(sizeof(unsigned int));
	int thread_num = omp_get_thread_num() + 1;
	int total_threads = omp_get_num_threads();
	*seed += thread_num*7 + seed_shift*13;
	*seed2 = thread_num*17 + seed_shift*11;
	init_genrand(*seed);
	
	COMPATIBLE_SPG compatible_spg[230]; 
	int num_compatible_spg = 0;
	int num_axes;
	float *mol_axes;
		
	molecule *mol = (molecule*)malloc(sizeof(molecule));
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));
	float volume;

	int spg;
	//read input from file, read molecular structure from geometry,Z, Zp
    #pragma omp critical
	{
        read_geometry(mol);

        //printf("%d %d %f %f %f\n", num_structures, Z, sr, volume_mean, volume_std);
    }
	
    //recenter molecule to origin
	recenter_molecule(mol);

	#pragma omp barrier
	{}
   
    //inititalise
    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
	int N = mol->num_of_atoms;
	allocate_xtal(random_crystal, ZMAX, N);
	random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
	
	if(thread_num == 1)
	{
		printf("Detecting compatible spacegroups"
			   " from molecule symmetries...\n");
	}
	
	find_compatible_spg_positions(mol,
								  Z,
								  compatible_spg,
								  &num_compatible_spg,
								  &mol_axes,
								  &num_axes,
								  thread_num);
			
	if(thread_num == 1)
	{
		printf("\nTotal compatible space groups = %d\n", num_compatible_spg);
		printf("Number of mol axes = %d \n", num_axes);
		printf("Starting generator ...\n\n");
		sleep(1);
	}
	
	#pragma omp barrier
	{}
	
	//find allowed space groups for general position (deprecated)
	//int allowed_spg[230];
	//int num_allowed_spg = 0;
	//find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
	//printf("allowed %d \n", num_compatible_spg);
	
	while( spg_rand < num_compatible_spg)
	{
		//conter counts the number of structures generated
		spg = compatible_spg[spg_rand].spg;
		while( counter < num_structures)
		{
			int verdict = 0;
			int i = 0; 
			// attempts for an spg
			for(; i < max_attempts/total_threads; i++)
			{
				if(counter >= num_structures)
					break;
				
				if(stop_flag)
					break;
					
				//if successful
				if (verdict == 1 )
				{
					#pragma omp critical
					if(counter < num_structures)
					{
						random_crystal->spg = spg;
						//print_crystal(random_crystal);
						print_crystal2file(random_crystal, out_file);
						counter++;
						//printf("#All Structure checks passed \n");
						printf("#thread %d:Generation Successful after %d attempts\n",
																		 thread_num,
																		 i*total_threads);
									
						printf("#Structure number = %d,  attempted spg = %d, volume = %f\n", 
																	counter, 
																	spg, volume);
						
						int spglib_spg = detect_spg_using_spglib(random_crystal);
						printf("#SPGLIB detected space group = %d\n",
												                     spglib_spg);
						fflush(stdout);
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
					}
									
					break;
				}
				// else generate again
				int result = generate_crystal(random_crystal,
								 mol,
								 volume,
								 Z,
								 Zp_max,
								 spg,
								 compatible_spg,
								 num_compatible_spg,
								 spg_rand,
								 mol_axes,
								 num_axes);
				
				//alignment failure
				if(!result)
					continue;
				
				//print_crystal(random_crystal); //exit(0);		    
				verdict = check_structure_with_vdw_matrix(*random_crystal, vdw_matrix, dim1, dim2);    
				
				if(i % VOL_ATTEMPT == 0 && i != 0)
				{
					do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
					if(thread_num == 1)
						printf("#thread 1: completed %d attempts\n", i*total_threads);
					*seed = *seed2 + thread_num*rand_r(seed2);
					fflush(stdout);
				}
				
			}//end of attempt loop	
			
			//if max limit is reached or if some thread hit the limit
			if (i >= max_attempts/total_threads || stop_flag == 1)
			{	
				//stop other threads
				#pragma omp critical
					{stop_flag = 1;}
				#pragma omp barrier
				{}
				do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
				if(thread_num == 1)
				{	
					printf("**WARNING: generation failed for spg = %d \
	after %d attempts \n\n", spg, max_attempts);
					printf("#counter reset...\n\n");
					counter = num_structures + 1;
					//print_crystal(random_crystal);
					fflush(stdout);
				}
				#pragma omp barrier
				{}
			}	
			
		}//end of numof structures llop
	
	#pragma omp barrier
	{}
	if(thread_num == 1 )
	{	
         counter = 0;
         spg_rand++;
         stop_flag = 0;
         printf("#counter reset...\n\n");
         fflush(stdout);
         //sleep(1);
         //usleep(90000);
    } 
	#pragma omp barrier
	{}
	
	}//while loop //end of spg loop
	
		//free memory
																	   }
	fclose(out_file);
	
	printf("Generation completed!\n");
	//free(mol_axes);
}



*/









/* atoms supported
 */
void create_crystal_from_array(crystal *xtal, double lattice_vector[3][3], double *Xc,int total_atoms1,
		double *Yc,int total_atoms2, double *Zc, int total_atoms3,char* atoms, int total_atoms, int Z, int spg)
{
	
	xtal->spg = spg;
	xtal->Z = Z;
	
	int num_atoms_in_molecule = total_atoms/Z;
	allocate_xtal(xtal, Z, num_atoms_in_molecule);
	xtal->num_atoms_in_molecule = num_atoms_in_molecule;
	
	xtal->lattice_vectors[0][0] = lattice_vector[0][0];
	xtal->lattice_vectors[0][1] = lattice_vector[0][1];
	xtal->lattice_vectors[0][2] = lattice_vector[0][2];
	
	xtal->lattice_vectors[1][0] = lattice_vector[1][0];
	xtal->lattice_vectors[1][1] = lattice_vector[1][1];
	xtal->lattice_vectors[1][2] = lattice_vector[1][2];
	
	xtal->lattice_vectors[2][0] = lattice_vector[2][0];
	xtal->lattice_vectors[2][1] = lattice_vector[2][1];
	xtal->lattice_vectors[2][2] = lattice_vector[2][2];
	
	for(int i = 0; i < total_atoms; i++)
	{
		xtal->Xcord[i] = Xc[i];
		xtal->Ycord[i] = Yc[i];
		xtal->Zcord[i] = Zc[i]; 
		xtal->atoms[2*i] = atoms[2*i];
		xtal->atoms[2*i+1] = atoms[2*i+1];
		
	}
	
	print_crystal(xtal);
}


int c_check_structure(crystal xtal, double sr)
{
	float f_sr = sr;
	return check_structure(xtal, f_sr);
}

/*
crystal generate_one_molecular_crystal(int Z, int spg, double volume1, double sr,
	double tol, int max_attempts, int seed_in)
{
	srand((unsigned int)time(NULL));
	int seed_shift = rand()% 1061 + 7 ;
	*seed = seedin*13+ seed_shift*7;
	int verdict;
	molecule *mol = (molecule*)malloc(sizeof(molecule));
	read_geometry(mol);
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));
	allocate_xtal(random_crystal, ZMAX, N);
	
	for(int i = 0; i < max_attempts; i++)
	{
		
	}
	 
}
*/

int num_compatible_spacegroups(int Z, double tolerance)
{
	//set global variable tolerance
	TOL = tolerance;
    seed2 = (unsigned int*)malloc(sizeof(unsigned int));

	COMPATIBLE_SPG compatible_spg[230]; 
	int num_compatible_spg = 0;
	int num_axes;
	float *mol_axes;
	int thread_num = 1; 
	molecule *mol = (molecule*)malloc(sizeof(molecule));

	//read geometry from geometry.in
	read_geometry(mol);

	find_compatible_spg_positions(mol,
								  Z,
								  compatible_spg,
								  &num_compatible_spg,
								  thread_num);

	free(mol_axes);
	free(mol);
    free(seed2);
	return num_compatible_spg;

}

