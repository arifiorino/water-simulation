//
//  animate.h
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#ifndef animate_h
#define animate_h

#define N 16
extern int n_particles;

typedef struct particle{
  float x,y,z;
  struct particle *next;
} particle_t;
extern particle_t *particles;
extern particle_t *particles_hash[N][N][N];

void init_animation(void);
void animate(void);

#endif /* animate_h */
