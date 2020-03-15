#include "read_input.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

//returns the atomic number for a given symbol
/*
int convert_atom_to_atomic_numbers(char *atom)
{
	if (strcmp(atom, "C") == 0)
		return 12;
	else if (strcmp(atom, "H") == 0)
		return 1;
}
*/
float TOL;


void read_control(int* num_structures, int* Z, float* Zp_max, 
	float* volume_mean, float* volume_std, float *sr, long *max_attempts)
{
	FILE *fileptr;
	size_t len = 0;
	char *line = NULL;
	char *sub_line = NULL;
	int read;
	fileptr = fopen("control.in","r");
	
	if(!fileptr)
	{
		printf("***ERROR: no control.in file \n");
		exit(0);
	}
	
	//read from control
	while ((read = getline(&line, &len, fileptr)) != -1)
	{
		//printf("Retrieved line of length %zu :\n", read);
		//if comment
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
		//printf("%s \n" , sub_line);
        if(strcmp(sub_line, "Z") == 0)
		{
			sub_line = strtok(NULL," ");
			*Z = atoi(sub_line);
			continue;
		}
		
		if(strcmp(sub_line, "sr") == 0)
		{
			sub_line = strtok(NULL," ");
			*sr = atof(sub_line);
			continue;
		}
		
		if(strcmp(sub_line, "volume_mean") == 0)
		{
			sub_line = strtok(NULL," ");
			*volume_mean = atof(sub_line);
			continue;
		}
		
		if(strcmp(sub_line, "volume_std") == 0)
		{
			sub_line = strtok(NULL," ");
			*volume_std = atof(sub_line);
			continue;
		}
		
		if(strcmp(sub_line, "number_of_structures") == 0)
		{
			sub_line = strtok(NULL," ");
			*num_structures = atoi(sub_line);
			continue;
		}
		
		if(strcmp(sub_line, "tolerance") == 0)
		{
			sub_line = strtok(NULL," ");
			TOL = atof(sub_line);
			continue;
		}	
		
		if(strcmp(sub_line, "max_attempts") == 0)
		{
			sub_line = strtok(NULL," ");
			*max_attempts = atoi(sub_line);
			continue;
		}	
			//printf("char = %c \n",mol->atoms[2*i+1]); exit(0);
            //printf("string = %c l %c \n", *sub_line, *(sub_line+1));	

   	}
	
	/*
    *num_structures = 1;
    *volume_mean = 1200;
    *volume_std  = 100;
	*Z = 6;
	*sr = 0.75;
	*/     
	fclose(fileptr);
    *Zp_max = 192;
}


void read_geometry(molecule* mol)
{
    FILE *fileptr;
	size_t len = 0;
    int read;	
	char *line = NULL;
	char *sub_line = NULL;
	int i = 0;
    int atom_count = 0;

    //find_number_of atoms
	fileptr = fopen("geometry.in","r");
	//check if file exits
	if(!fileptr)
	{
		printf("***ERROR: no geometry.in file \n");
		exit(0);
	}

     while ((read = getline(&line, &len, fileptr)) != -1)
	{
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
    }
    fclose(fileptr);

    //printf("Total number of atoms in molecule = %d\n", atom_count);
    int N = atom_count;
	//memory allocation
	(*mol).atoms = (char *)malloc(2*N*sizeof(char));
	(*mol).X = (float *)malloc(N*sizeof(float));
	(*mol).Y = (float *)malloc(N*sizeof(float));
	(*mol).Z = (float *)malloc(N*sizeof(float));
	
	//creates a file pointer and opens geometry.in
	    
    fileptr = fopen("geometry.in","r");
    while ((read = getline(&line, &len, fileptr)) != -1)
	{
		//printf("Retrieved line of length %zu :\n", read);
		//printf("%s", line);
        if (strstr(line, "#") != NULL)
            continue;

		sub_line=strtok(line," ");
        //printf("%s \n" , sub_line);
        if(strcmp(sub_line, "atom") == 0)
            atom_count++;
        else
            continue;
		sub_line=strtok(NULL," ");
		(*mol).X[i]=atof(sub_line);
		//	printf("%f \t",X[i-3]);	
		sub_line=strtok(NULL," ");
		(*mol).Y[i]=atof(sub_line);
		//	printf("%f \t",Y[i-3]);
		sub_line=strtok(NULL," ");
		(*mol).Z[i]=atof(sub_line);
		//	printf("%f \t",Z[i-3]);	
		sub_line=strtok(NULL," ");
		(*mol).atoms[2*i]=*sub_line;
        if(*(sub_line+1) == '\n' || *(sub_line+1) == ' ' || *(sub_line+1) == '\0' )
            (*mol).atoms[2*i+1]=' ';            
        else
            (*mol).atoms[2*i+1]=*(sub_line+1);
			//printf("char = %c \n",mol->atoms[2*i+1]); exit(0);
            //printf("string = %c l %c \n", *sub_line, *(sub_line+1));	
		i++;

   	}
 

	fclose(fileptr);
   // printf("atoms = %d \n", atom_count);	
	mol->num_of_atoms = N;
	

}


void print_input_settings(int* num_structures, int* Z, float* Zp_max, 
	float* volume_mean, float* volume_std, float *sr, long *max_attempts)
{
    *Zp_max = 192; //useless argument
	printf("INPUT SETTINGS:\n");
	printf("-----------------------------\n");
	printf("Number of structures per space group:         %d \n", *num_structures);
	printf("Number of molecules in the cell:              %d\n", *Z);
	printf("Mean volume of unit cell:                     %f\n", *volume_mean);
	printf("Standard deviation of unit cell volume:       %f\n", *volume_std);
	printf("Specific radius proportion:                   %f\n", *sr );
	printf("Maximum attempts per space group:             %ld\n", *max_attempts);
	printf("Tolerance:                                    %f\n", TOL);
	printf("-----------------------------\n\n");

}


void print_input_geometry(molecule* mol)
{
	printf("MOLECULAR GEOMETRY:\n");
	printf("-----------------------------\n");
	print_molecule(mol);
	printf("-----------------------------\n");
	printf("\n");
}


