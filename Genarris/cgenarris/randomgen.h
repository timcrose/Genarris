#ifndef RANDOMGEN_H
#define RANDOMGEN_H

void init_genrand(unsigned int s);
unsigned long genrand_int32(void);
float uniform_dist_01 (void);
float normal_dist_01 (void);
float normal_dist_ab ( float a, float b);

#endif
