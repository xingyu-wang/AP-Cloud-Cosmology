#include "viewer.h"

Viewer::Viewer(Physical_particle *physical_particle1, Computational_particle *computational_particle1)
: pp(physical_particle1), cp(computational_particle1)
{
	Set_File_Name(OutputFile,OutputFile);
	strcat (file_name_p,"/physical_");
	strcat (file_name_c,"/computational_");
}


void Viewer::View_Physical_Particles(double t)
{
	int i;
  char file_name_temp[200],temp[200];
	strcpy(file_name_temp,file_name_p);
	sprintf(temp,"%f",t);
	strcat (temp,".vtk");
	strcat (file_name_temp,temp);
	file_p = fopen(file_name_temp,"w");

	if (file_p == NULL) perror ("Error opening file");
	fprintf(file_p,"# vtk DataFile Version 3.0\n");
	fprintf(file_p,"Data on physical particles for cosmology problem\n"); 
	fprintf(file_p,"ASCII\n");
	fprintf(file_p,"DATASET POLYDATA\n");
	fprintf(file_p,"POINTS %d double\n\n",pp->size);
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%f %f %f\n",pp->x[i],pp->y[i],pp->z[i]);
	fprintf(file_p,"\nPOINT_DATA %d",pp->size);
	fprintf(file_p,"\nVECTORS velocity double\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e %e %e\n",pp->v_x[i],pp->v_y[i],pp->v_z[i]);
	fprintf(file_p,"\nSCALARS v_x double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e\n",pp->v_x[i]);
	fprintf(file_p,"\nSCALARS v_y double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e\n",pp->v_y[i]);
	fprintf(file_p,"\nSCALARS v_z double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e\n",pp->v_z[i]);
	fprintf(file_p,"\nVECTORS acceleration double\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e %e %e\n",pp->a_x[i],pp->a_y[i],pp->a_z[i]);
	fprintf(file_p,"\nSCALARS a_x double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e\n",pp->a_x[i]);
	fprintf(file_p,"\nSCALARS a_y double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e\n",pp->a_y[i]);
	fprintf(file_p,"\nSCALARS a_z double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	for (i=0;i<pp->size;i++)
		fprintf(file_p,"%e\n",pp->a_z[i]);
	fclose (file_p);
}

void Viewer::View_Computational_Particles(double t)
{
	int i;
  char file_name_temp[200],temp[200];
	strcpy(file_name_temp,file_name_c);
	sprintf(temp,"%f",t);
	strcat (temp,".vtk");
	strcat (file_name_temp,temp);

	file_c = fopen(file_name_temp,"w");
	if (file_c == NULL) perror ("Error opening file");
	fprintf(file_c,"# vtk DataFile Version 3.0\n");
	fprintf(file_c,"Data on computational particles for cosmology problem\n"); 
	fprintf(file_c,"ASCII\n");
	fprintf(file_c,"DATASET POLYDATA\n");
	fprintf(file_c,"POINTS %d double\n\n",cp->size);
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%f %f %f\n",cp->x[i],cp->y[i],cp->z[i]);
	fprintf(file_c,"\nPOINT_DATA %d\n",cp->size);
	fprintf(file_c,"SCALARS n_p int\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%d\n",cp->NumberOfPhysicalParticles[i]);
	fprintf(file_c,"\nSCALARS type int\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%d\n",cp->Type[i]);
	fprintf(file_c,"\nSCALARS first_id int\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%d\n",cp->IndexOfFirstPhysicalParticle[i]);
	fprintf(file_c,"\nSCALARS density double\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%e\n",cp->Density[i]);
	fprintf(file_c,"\nSCALARS potential double\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%e\n",cp->Potential[i]);
	fprintf(file_c,"\nSCALARS p_x double\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%e\n",cp->p_x[i]);
	fprintf(file_c,"\nSCALARS p_y double\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%e\n",cp->p_y[i]);
	fprintf(file_c,"\nSCALARS p_z double\n");
	fprintf(file_c,"LOOKUP_TABLE default\n");
	for (i=0;i<cp->size;i++)
		fprintf(file_c,"%e\n",cp->p_z[i]);
	fclose (file_c);
}

	void Viewer::View_Physical_Particles_Limit(double t, double xmin, double ymin, double zmin, double xmax, double ymax, double zmax)
{
	int i;
  char file_name_temp[200],temp[200];
	strcpy(file_name_temp,file_name_p);
	sprintf(temp,"%f",t);
	strcat (temp,".vtk");
	strcat (file_name_temp,temp);
	file_p = fopen(file_name_temp,"w");

	if (file_p == NULL) perror ("Error opening file");
	fprintf(file_p,"# vtk DataFile Version 3.0\n");
	fprintf(file_p,"Data on physical particles for cosmology problem\n"); 
	fprintf(file_p,"ASCII\n");
	fprintf(file_p,"DATASET POLYDATA\n");
	fprintf(file_p,"POINTS %d double\n\n",2);
	fprintf(file_p,"%f %f %f\n",xmin,ymin,zmin);
	fprintf(file_p,"%f %f %f\n",xmax,ymax,zmax);
	fprintf(file_p,"\nPOINT_DATA %d",2);
	fprintf(file_p,"\nVECTORS velocity double\n");
	for (i=0;i<2;i++)
		fprintf(file_p,"%f %f %f\n",0.0,0.0,0.0);
	fprintf(file_p,"\nSCALARS v_x double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"\nSCALARS v_y double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"\nSCALARS v_z double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"\nVECTORS acceleration double\n");
	for (i=0;i<2;i++)
		fprintf(file_p,"%f %f %f\n",0.0,0.0,0.0);
	fprintf(file_p,"\nSCALARS a_x double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"\nSCALARS a_y double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"\nSCALARS a_z double\n");
	fprintf(file_p,"LOOKUP_TABLE default\n");
	fprintf(file_p,"%f\n",0.0);
	fprintf(file_p,"%f\n",0.0);
	fclose (file_p);
}
