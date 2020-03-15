#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "algebra.h"
#include "randomgen.h"

#define PI 3.141592653
extern unsigned int *seed2;


void vector3_add(float a[3], float b[3], float sum[3])
{
	sum[0] = a[0] + b[0];
	sum[1] = a[1] + b[1];
	sum[2] = a[2] + b[2];
	return;
}

/* Both b and c can be the same vector; hence a temporary variable
 * is needed
 */
void vector3_mat3b3_multiply(float a[3][3], float b[3], float c[3])
{
	float temp[3];
	temp[0] = a[0][0] * b[0] + a[0][1] * b[1] + a[0][2] * b[2];
	temp[1] = a[1][0] * b[0] + a[1][1] * b[1] + a[1][2] * b[2];
	temp[2] = a[2][0] * b[0] + a[2][1] * b[1] + a[2][2] * b[2];
	c[0]  = temp[0];
	c[1]  = temp[1];
	c[2]  = temp[2];
	return;
}

void vector3_intmat3b3_multiply(int a[3][3], float b[3], float c[3])
{
	float temp[3];
	for (int i = 0; i < 3; i++)
		temp[i] = a[i][0] * b[0] + a[i][1] * b[1] + a[i][2] * b[2];
	copy_vector3_vector3(c,temp);
	return;
}


//tested
void mat3b3_mat3b3_multiply(float a[3][3], float b[3][3], float c[3][3])
{
	float temp[3][3];
	for (int i = 0; i < 3; i++) 
    { 
        for (int j = 0; j < 3; j++) 
        { 
            temp[i][j] = 0; 
            for (int k = 0; k < 3; k++) 
                temp[i][j] += a[i][k]*b[k][j]; 
        } 
    } 
	
	copy_mat3b3_mat3b3(c,temp);
	return;
}

void rotation_mat_around_axis(float rot[3][3], float axis[3], float psi)
{
	float cosp = cos(psi);
	float vercosp = 1 - cosp; //Thanks S.L. Loney
	float sinp = sin(psi);

	//using Rodrigues rotation formula
	rot[0][0] = cosp + axis[0]*axis[0]*vercosp;
	rot[0][1] = axis[0]*axis[1]*vercosp - axis[2]*sinp;
	rot[0][2] = axis[0]*axis[2]*vercosp + axis[1]*sinp;

	rot[1][0] = axis[0]*axis[1]*vercosp + axis[2]*sinp;
	rot[1][1] = cosp + axis[1]*axis[1]*vercosp;
	rot[1][2] = axis[1]*axis[2]*vercosp - axis[0]*sinp;
	
	rot[2][0] = axis[0]*axis[2]*vercosp - axis[1]*sinp;
	rot[2][1] = axis[1]*axis[2]*vercosp + axis[0]*sinp;
	rot[2][2] = cosp + axis[2]*axis[2]*vercosp;
}	

void print_mat3b3(float mat[3][3])
{
	for(int i = 0; i < 3; i++)
	{
		for(int j = 0; j < 3; j++)
			printf("%f   ", mat[i][j]);
		printf("\n");
	}
}

void print_vec3(float vec[3])
{
	for(int i = 0; i < 3; i++)
		printf("%f   ", vec[i]);
	printf("\n");
}

float det_mat3b3(float a[3][3])
{
	 return a[0][0] * ((a[1][1]*a[2][2]) - (a[2][1]*a[1][2])) -a[0][1] 
	 * (a[1][0] * a[2][2] - a[2][0] * a[1][2]) + a[0][2] * (a[1][0] 
	 * a[2][1] - a[2][0] * a[1][1]);
}

void copy_vector3_vector3(float a[3], float b[3])
{
	a[0] = b[0];
	a[1] = b[1];
	a[2] = b[2];
	return;
}

void copy_mat3b3_mat3b3(float a[3][3], float b[3][3])
{
	a[0][0] = b[0][0];
	a[0][1] = b[0][1];
	a[0][2] = b[0][2];
	
	a[1][0] = b[1][0];
	a[1][1] = b[1][1];
	a[1][2] = b[1][2];
	
	a[2][0] = b[2][0];
	a[2][1] = b[2][1];
	a[2][2] = b[2][2];
	
	return;
}

void copy_mat3b3_mat3b3bN(float a[3][3], float b[][3][3], int index)
{
	a[0][0] = b[index][0][0];
	a[0][1] = b[index][0][1];
	a[0][2] = b[index][0][2];
	
	a[1][0] = b[index][1][0];
	a[1][1] = b[index][1][1];
	a[1][2] = b[index][1][2];
	
	a[2][0] = b[index][2][0];
	a[2][1] = b[index][2][1];
	a[2][2] = b[index][2][2];
	
	return;	
}

void copy_mat3b3_intmat3b3bN(float a[3][3], const int b[][3][3], int index)
{
	a[0][0] = b[index][0][0];
	a[0][1] = b[index][0][1];
	a[0][2] = b[index][0][2];
	
	a[1][0] = b[index][1][0];
	a[1][1] = b[index][1][1];
	a[1][2] = b[index][1][2];
	
	a[2][0] = b[index][2][0];
	a[2][1] = b[index][2][1];
	a[2][2] = b[index][2][2];
	
	return;	
}

void copy_mat3b3bN_mat3b3(float b[][3][3], float a[3][3], int index)
{
	b[index][0][0] = a[0][0];
	b[index][0][1] = a[0][1];
	b[index][0][2] = a[0][2];
	
	b[index][1][0] = a[1][0];
	b[index][1][1] = a[1][1];
	b[index][1][2] = a[1][2];
	
	 b[index][2][0] = a[2][0];
	 b[index][2][1] = a[2][1];
	 b[index][2][2] = a[2][2];
	
	return;	
}
//for coying integer matrices
void copy_intmat3b3_intmat3b3bN(int a[3][3], int b[][3][3], int index)
{
	a[0][0] = b[index][0][0];
	a[0][1] = b[index][0][1];
	a[0][2] = b[index][0][2];
	
	a[1][0] = b[index][1][0];
	a[1][1] = b[index][1][1];
	a[1][2] = b[index][1][2];
	
	a[2][0] = b[index][2][0];
	a[2][1] = b[index][2][1];
	a[2][2] = b[index][2][2];
	
	return;	
}

void copy_vector3_vector3bN(float a[3], const float b[][3], int index)
{
	a[0] = b[index][0];
	a[1] = b[index][1];
	a[2] = b[index][2];
}

//tested
void inverse_mat3b3(float Tinv[3][3], float T[3][3])
{
	float det = det_mat3b3(T);

    if(det == 0)
       {printf("error: det = 0, matrix has no inverse \n"); print_mat3b3(T);}
 
	for(int i = 0;i < 3; i++)
	{
    	for(int j = 0; j < 3; j++)
		{
			Tinv[j][i]=((T[(i+1)%3][(j+1)%3] * T[(i+2)%3][(j+2)%3]) - 
				(T[(i+1)%3][(j+2)%3]*T[(i+2)%3][(j+1)%3]))/ det;
		}
    }
}

int hall_number_from_spg(int spg)
{
		static int spg_to_hall_number[230] =
	{
    1,   2,   3,   6,   9,  18,  21,  30,  39,  57,
   60,  63,  72,  81,  90, 108, 109, 112, 115, 116,
  119, 122, 123, 124, 125, 128, 134, 137, 143, 149,
  155, 161, 164, 170, 173, 176, 182, 185, 191, 197,
  203, 209, 212, 215, 218, 221, 227, 228, 230, 233,
  239, 245, 251, 257, 263, 266, 269, 275, 278, 284,
  290, 292, 298, 304, 310, 313, 316, 322, 334, 335,
  337, 338, 341, 343, 349, 350, 351, 352, 353, 354,
  355, 356, 357, 358, 359, 361, 363, 364, 366, 367,
  368, 369, 370, 371, 372, 373, 374, 375, 376, 377,
  378, 379, 380, 381, 382, 383, 384, 385, 386, 387,
  388, 389, 390, 391, 392, 393, 394, 395, 396, 397,
  398, 399, 400, 401, 402, 404, 406, 407, 408, 410,
  412, 413, 414, 416, 418, 419, 420, 422, 424, 425,
  426, 428, 430, 431, 432, 433, 435, 436, 438, 439,
  440, 441, 442, 443, 444, 446, 447, 448, 449, 450,
  452, 454, 455, 456, 457, 458, 460, 462, 463, 464,
  465, 466, 467, 468, 469, 470, 471, 472, 473, 474,
  475, 476, 477, 478, 479, 480, 481, 482, 483, 484,
  485, 486, 487, 488, 489, 490, 491, 492, 493, 494,
  495, 497, 498, 500, 501, 502, 503, 504, 505, 506,
  507, 508, 509, 510, 511, 512, 513, 514, 515, 516,
  517, 518, 520, 521, 523, 524, 525, 527, 529, 530,
	};

	return spg_to_hall_number[spg-1];
	
}

void copy_floatmat3b3_intmat3b3(float a[3][3], int b[3][3])
{
	a[0][0] = b[0][0];
	a[0][1] = b[0][1];
	a[0][2] = b[0][2];
	
	a[1][0] = b[1][0];
	a[1][1] = b[1][1];
	a[1][2] = b[1][2];
	
	a[2][0] = b[2][0];
	a[2][1] = b[2][1];
	a[2][2] = b[2][2];
	
}
 void mat3b3_transpose(float b_trans[3][3], float b[3][3])
{
		float temp[3][3];
		
		temp[0][0] = b[0][0];
		temp[1][1] = b[1][1];	
		temp[2][2] = b[2][2];
		
		temp[1][0] = b[0][1];
		temp[2][0] = b[0][2];
		temp[2][1] = b[1][2];
		
		temp[0][1] = b[1][0];
		temp[0][2] = b[2][0];
		temp[1][2] = b[2][1];
		
		copy_mat3b3_mat3b3(b_trans, temp);
}

void print_mat3b3bN(float a[][3][3], int N)
{
	for(int i= 0;i < N; i++)
	{
		float temp[3][3];
		copy_mat3b3_mat3b3bN(temp,a,i);
		print_mat3b3(temp);
		printf("\n");
	}
}

void copy_vector3_mat3b3(float a[3], float b[3][3],int index)
{
	a[0] = b[index][0];
	a[1] = b[index][1];
	a[2] = b[index][2];
	
}

void vector3_inverse(float a[3])
{
	a[0] = -a[0];
	a[1] = -a[1];
	a[2] = -a[2];
}

float vector3_norm(float a[3])
{
	return sqrt(a[0]*a[0] + a[1]*a[1] + a[2]*a[2]);
}

void vector3_subtract(float a[3], float b[3], float diff[3])
{
	diff[0] = a[0] - b[0];
	diff[1] = a[1] - b[1];
	diff[2] = a[2] - b[2];
	return;
}

void cross_vector3_vector3(float cross[3], float a[3], float b[3])
{
	cross[0] = a[1]*b[2] - a[2]*b[1];
	cross[1] = a[2]*b[0] - a[0]*b[2];
	cross[2] = a[0]*b[1] - a[1]*b[0];

	return;
}

void normalise_vector3(float a[3])
{
	float norm;
	norm = vector3_norm(a);
	a[0] /= norm;
	a[1] /= norm;
	a[2] /= norm;
}

float dot_vector3_vector3(float a[3], float b[3])
{
	return a[0]*b[0] + a[1]*b[1] + a[2]*b[2];
}

int check_vec3_isNull(float a[3], float tol)
{
	if (fabs(a[0]) < tol &&
		fabs(a[1]) < tol &&
		fabs(a[2]) < tol  )
		return 1;
	else
		return 0;
}

float cart_dist(float p1[3], float p2[3])
{
		float c[3];
		vector3_subtract(p1,p2,c);
		return vector3_norm(c);
}

void generate_random_rotation_matrix( float rotation_matrix[3][3] )
{
		float phi = 2*PI*uniform_dist_01();
		//theta, angle with x-axis. this should acos(u), -1<u<1
		float u = 2*uniform_dist_01() - 1;
		float theta = acos(u);
		//compute random axis with theta and phi
		float axis[] = {cos(phi)*sin(theta),
						sin(phi)*sin(theta),
						cos(theta)};
		float psi = 2*PI*uniform_dist_01();
		rotation_mat_around_axis(rotation_matrix, axis, psi);
}

//fisher-yates shuffle algorithn
void array_shuffler_2(float a[][3], int len )
{
	for(int i = len - 1; i > 0; i--)
	{
		int j = rand_r(seed2) %(i + 1);
		
		//swap i and j element
		float temp[3] = {a[i][0], a[i][1], a[i][2]};
		a[i][0] = a[j][0];
		a[i][1] = a[j][1];
		a[i][2] = a[j][2];
		
		a[j][0] = temp[0];
		a[j][1] = temp[1];
		a[j][2] = temp[2];
	}
}

//fisher-yates shuffle algorithn
void array_shuffler_1(float *a, int len )
{
	for(int i = len - 1; i > 0; i--)
	{
		int j = rand_r(seed2) %(i + 1);
		
		//swap i and j element
		float temp[3] = {*(a+3*i+0), *(a+3*i+1), *(a+3*i+2)};
		*(a+3*i + 0) = *(a+3*j+0);
		*(a+3*i + 1) = *(a+3*j+1);
		*(a+3*i + 2) = *(a+3*j+2);
		
		*(a+3*j+0) = temp[0];
		*(a+3*j+1) = temp[1];
		*(a+3*j+2) = temp[2];
	}
}
