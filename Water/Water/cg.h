//
//  cg.h
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#ifndef cg_h
#define cg_h

void mallocCG(int n, float max_row);
void initCG(void);
void write_A(int i, int j, float x);
void write_b(int i, float x);
float *cg(void);
void testCG(void);

#endif /* cg_h */
