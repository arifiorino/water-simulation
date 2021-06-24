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
#include <time.h>
#include <limits.h>
#include <vector>
#include <algorithm>
#include <thread>
#include <functional>
void **malloc2D(int w, int h, int s);
void free2D(void **p, int w);
double rand2(void);
double timestamp(void);
void cross_prod(float ax, float ay, float az, float bx, float by, float bz,
                float *sx, float *sy, float *sz);
void parallel_for(int n, std::function<void (int start, int end)> functor);

#endif /* utils_h */
