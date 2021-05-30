//
//  utils.h
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#ifndef utils_h
#define utils_h

#include <float.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
void **malloc2D(int w, int h, int s);
void free2D(void **p, int w);
void swap2D(float ***a, float ***b);
void print2D(float **matrix, int w, int h);
void print2DInt(int **matrix, int w, int h);

#endif /* utils_h */
