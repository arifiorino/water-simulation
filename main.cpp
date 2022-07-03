#include <GL/glut.h>
#include "render.h"
#include "animate.h"
 
int refreshMills = 15;

void display() {
  glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  glClear(GL_COLOR_BUFFER_BIT);

  animate();
  render();
  float *triangles2=(float*)malloc(sizeof(float)*n_triangles*3*3);
  for (int i=0;i<n_triangles*3;i++){
    triangles2[i*3]=triangles[i*2];
    triangles2[i*3+1]=triangles[i*2+1];
    triangles2[i*3+2]=0.0f;
  }
  
  glEnableClientState(GL_VERTEX_ARRAY);
  glColor3f(0.0f, 0.0f, 1.0f);
  glVertexPointer(3, GL_FLOAT, 0, triangles2);
  glDrawArrays(GL_TRIANGLES, 0, n_triangles*3);
  glDisableClientState(GL_VERTEX_ARRAY);
  
  glFlush();
  free(triangles2);
}

void timer(int value) {
  glutPostRedisplay();
  glutTimerFunc(refreshMills, timer, 0);
}
 
int main(int argc, char** argv) {
  initAnimation();
  glutInit(&argc, argv);
  glutInitWindowSize(1000, 1000);
  glutInitWindowPosition(50, 50);
  glutCreateWindow("Water");
  glutDisplayFunc(display);
  glutTimerFunc(0, timer, 0);
  glutMainLoop();
  return 0;
}
