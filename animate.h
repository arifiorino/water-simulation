//
//  animate.h
//  Water
//
//  Created by Ari Fiorino on 5/28/21.
//

#ifndef animate_h
#define animate_h

extern int N;
extern int n_particles;

typedef struct particle{
  float x,y;
  struct particle *next;
} particle_t;
extern particle_t *particles;
extern particle_t ***particles_hash;

void initAnimation(void);
void animate(void);

#endif /* animate_h */
