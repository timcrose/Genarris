#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <omp.h>
#include "read_input.h"
#include "spg_generation.h"
#include "combinatorics.h"
#include "check_structure.h"
#include "crystal_utils.h"
#include "molecule_utils.h"
#include "lattice_generator.h"
#include "randomgen.h"
#include "algebra.h"

//maximum mulipicity possible
#define ZMAX 192
#define VOL_ATTEMPT  100000

//seeds are private to threads and initialised diffrently
unsigned int *seed;
unsigned int *seed2;
#pragma omp threadprivate(seed)
#pragma omp threadprivate(seed2)

int main(int argc, char **argv)
{
	//Testing
	/*
	float T[3][3] = {27.67178 ,   0.00000 ,   0.00000,  2.56691  ,  7.12630  ,  0.00000 ,4.92891 ,   4.40403   , 3.24262};
	float T_inv[3][3] = { 0.03614   ,0.00000  , 0.00000, -0.01302 ,  0.14033,   0.00000, -0.03725 , -0.19059  , 0.30839};
	float x[3] = { 13.80239 ,   8.44218  ,  0.44715};
	float y[3] = {21.3550   , 4.0151 ,   3.7916};
	float p_dist, p_dist_2;

	//pdist_2(T, T_inv, x[0], x[1], x[2], y[0], y[1], y[2], &p_dist, &p_dist_2);
	p_dist = pdist(T, T_inv, x[0], x[1], x[2], y[0], y[1], y[2]);
	printf("%f \n", p_dist);

	exit(0);
	*/
	//end testing

	//random number seeding
	srand((unsigned int)time(NULL));
	int seed_shift = rand()% 1061 + 7 ;
	//variable declarations
	int stop_flag = 0;	// to stop threads if limit is reached
	int counter = 0;	//counts number of structures
	int spg_rand = 0;	//space group to be generated
	FILE *out_file;		//file to output geometries
	out_file = fopen("geometry.out","w");
	if(!out_file)		//check permissions
	{
		printf("***ERROR: cannot create geometry.out \n");
		exit(0);
	}
	
	#pragma omp parallel shared(counter, spg_rand, seed_shift,\
		out_file, stop_flag, stdout) default(none)						//omp block starts here
																		{
	//random number seeding, different seeds for different threads
	seed = (unsigned int*)malloc(sizeof(unsigned int)); //seed for uniform gen
	seed2 = (unsigned int*)malloc(sizeof(unsigned int)); //seed for random
	int thread_num = omp_get_thread_num() + 1; //get threads
	int total_threads = omp_get_num_threads();
	*seed += thread_num*7 + seed_shift*13; //some random seed private for each threads
	*seed2 = thread_num*17 + seed_shift*11;
	init_genrand(*seed);

	//storing information for compatible space groups
	COMPATIBLE_SPG compatible_spg[230]; 
	int num_compatible_spg = 0;
	int num_axes;	//storing number molecule axes
	float *mol_axes; //storing possible molecular axes
	
	//variable declarartion	
	molecule *mol = (molecule*)malloc(sizeof(molecule));//store molecule
	crystal *random_crystal = (crystal*)malloc(sizeof(crystal));//dummy crystal
	float volume_std;	//standard dev for volumes
	float volume_mean;	//mean volume
	float sr;			//specific radius proportion for structure checks
						//see paper for definition
	float Zp_max;		//Z'' . not implemented
	float volume;		//random volume used of generation
	int Z;				//multiplicity of general position
	int num_structures;	//num of structures
	int spg;			//space group attempted
	long max_attempts;	//max attempts per space group

	//read input from file, read molecular structure from geometry,Z, Zp

	if(thread_num == 1)
	{
		int len = 100;
		char name[len];
		gethostname(name, len);
		printf("PARALLELIZATION INFO:\n");
		printf("---------------------------\n");
		printf("Using OpenMP for parallelization.\n");
		printf("cgenarris created %d threads on %s. \n", total_threads, name);
		printf("NOTE: Number of threads can be changed using the environment variable OMP_NUM_THREADS.\n");
		printf("---------------------------\n");
		printf("\n");
	}



    #pragma omp critical
	{
        read_geometry(mol);				//read molecule from geometry.in
        read_control(&num_structures,
					 &Z,
					 &Zp_max,
					 &volume_mean, 
					 &volume_std,
					 &sr,
					 &max_attempts);	//get settings
    }
	
	//recenter molecule to origin
	recenter_molecule(mol);

	if (thread_num == 1)
	{
		print_input_geometry(mol);
		print_input_settings(&num_structures,
							 &Z,
							 &Zp_max,
							 &volume_mean, 
							 &volume_std,
							 &sr,
							 &max_attempts);
	}


	#pragma omp barrier //wait to get the settings
	{}
   
    //inititalise volume
    do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
	int N = mol->num_of_atoms;
	allocate_xtal(random_crystal, ZMAX, N); //allcate memory
	random_crystal->num_atoms_in_molecule = mol->num_of_atoms;
	
	if(thread_num == 1)
	{
		printf("COMPATIBLE SPACE GROUP INFO:\n");
		printf("-----------------------------------------------------\n");
		printf("Detecting compatible space groups"
			   " from molecule symmetries...\n");
	}

	//compatible space groups are stored in compatible_spg. see the 
	//object definition in spg_generation.h
	//every thread has its own copy
	find_compatible_spg_positions(mol,
								  Z,
								  compatible_spg,
								  &num_compatible_spg,
								  &mol_axes,
								  &num_axes,
								  thread_num);
			
	if(thread_num == 1)
	{
		printf("\n");
		printf("Total compatible space groups = %d\n", num_compatible_spg);
		printf("Number of molecular axes = %d \n", num_axes);
		printf("-----------------------------------------------------\n\n");
		printf("Starting generation...\n\n");
		sleep(1);
	}
	
	#pragma omp barrier // wait for other friends to join
	{}
	
	/*deprecated
	//find allowed space groups for general position (deprecated)
	//int allowed_spg[230];
	//int num_allowed_spg = 0;
	//find_allowed_spg(allowed_spg, &num_allowed_spg, Z);
	//printf("allowed %d \n", num_compatible_spg);
	*/
	
	
	while( spg_rand < num_compatible_spg )
	{
		//counter counts the number of structures generated
		spg = compatible_spg[spg_rand].spg; //pick a spg 

		//print attempted space group 
		if (thread_num == 1)
			{printf("Attempting to generate spacegroup number %d....\n", spg);}

		while( counter < num_structures )
		{
			int verdict = 0; //for structure check
			int i = 0; 		 //counts attempts for spg
			//attempts for an spg. Assume all threads run equally fast
			for(; i < max_attempts/total_threads; i++) 
			{
				if(counter >= num_structures)
					break;
				
				if(stop_flag)
					break;
					
				//if generation is successful
				if (verdict == 1 )
				{
					#pragma omp critical //one at a time please
					if(counter < num_structures)
					{
						random_crystal->spg = spg;
						//print_crystal(random_crystal);
						print_crystal2file(random_crystal, out_file);
						counter++; //one less structure to generate
						printf("#thread %d:Generation Successful after %d attempts\n",
								thread_num,
								i*total_threads);
									
						printf("#Structure number = %d ,\n"
							   "#attempted space group = %d,\n"
							   "#unit cell volume (cubic Angstrom)= %f\n", 
							   counter, 
							   spg,
							   volume);
						fflush(stdout);
						int spglib_spg = detect_spg_using_spglib(random_crystal);
						printf("#SPGLIB detected space group = %d\n\n",
												                     spglib_spg);
						do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
						i = 0;
						verdict = 0;
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
				
				//check if molecules are too close with sr	    
				verdict = check_structure(*random_crystal, sr);    
				
				//reset volume after volume attempts
				if(i % VOL_ATTEMPT == 0 && i != 0)
				{
					do {volume = normal_dist_ab(volume_mean, volume_std);} while(volume < 0.1);
					if(thread_num == 1)
						printf("#thread1:completed %d attempts\n",
							i*total_threads);
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
					printf("**WARNING: generation failed for space group = %d "
							"after %ld attempts. \n",
							spg,
							max_attempts);
					counter = num_structures + 1;
					fflush(stdout);
					//print_crystal(random_crystal);
				}
				#pragma omp barrier
				{}
			}	
			
		}//end of numof structures whileloop
	
		//move to next spacegroup
		#pragma omp barrier
		{}
		if(thread_num == 1 )
		{	
			 counter = 0;
			 spg_rand++;
			 stop_flag = 0;
			 printf("#space group counter reset. Moving to next space group...\n\n");
			 fflush(stdout);
		} 
		#pragma omp barrier
		{}
		
	}//end of spg while loop

																	   } //omp block ends here.
	fclose(out_file);
	printf("Generation completed!\n");

}

