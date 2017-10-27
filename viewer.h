#include <string>
#include "particle.h"
#include "parameter.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef VIEWER

class Viewer {
private:
	FILE *file_p;
	FILE *file_c;
  char file_name_p[200];
  char file_name_c[200];
private:
	Physical_particle *pp;
	Computational_particle *cp;

public:
	Viewer(Physical_particle *physical_particle1, Computational_particle *computational_particle1);
	int Set_File_Name(char * p_name, char * c_name){strcpy(file_name_p,p_name);strcpy(file_name_c,c_name);return 0;}
	void View_Physical_Particles(double t);
	void View_Physical_Particles_Limit(double t, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax);
	void View_Computational_Particles(double t);
};

#define VIEWER
#endif
