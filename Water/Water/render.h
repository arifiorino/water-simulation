//
//  render.h
//  Water
//
//  Created by Ari Fiorino on 6/5/21.
//

#ifndef render_h
#define render_h

#include <stdio.h>

int n_indices;
uint *indices;

int n_vertices;
float *vertices;
float *normals;

void init_render(void);
void render(void);

#endif /* render_h */
