#include "apcloud_macro.h"
#include "initialization.h"

int initialization(Physical_particle *physical_particle)
{
	int errorflag;
	switch (ICFileType)
	{
	case 1://ascii .vtk file
		errorflag=initialization_vtk(physical_particle);
		break;
	default:
		errorflag=1;
		printf("ICFileType invalid!\n");
	}
	if(errorflag)
		exit(1);
	return errorflag;
}

int initialization_vtk(Physical_particle *physical_particle)
{
	int errorflag=0;
	int i;
	FILE *file;

	file = fopen(ICFile,"r");
	if(file==NULL)
	{
		printf("Error in opening input file %s\n",ICFile);
		errorflag=1;
		exit(1);
	}
	char buf[200];
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fscanf(file,"%s %d",buf, &(physical_particle->size));

  physical_particle->x = new double[physical_particle->size];
  physical_particle->y = new double[physical_particle->size];
  physical_particle->z = new double[physical_particle->size];
  physical_particle->v_x = new double[physical_particle->size];
  physical_particle->v_y = new double[physical_particle->size];
  physical_particle->v_z = new double[physical_particle->size];
  physical_particle->a_x = new double[physical_particle->size];
  physical_particle->a_y = new double[physical_particle->size];
  physical_particle->a_z = new double[physical_particle->size];
  physical_particle->ID = new int[physical_particle->size];

  fgets(buf, 200, file);
	for(i=0;i<physical_particle->size;i++)
	{
		fscanf(file,"%lf %lf %lf",&(physical_particle->x[i]),&(physical_particle->y[i]),&(physical_particle->z[i]));
	}
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);

	for(i=0;i<physical_particle->size;i++)
		fscanf(file,"%lf",&(physical_particle->v_x[i]));
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	for(i=0;i<physical_particle->size;i++)
		fscanf(file,"%lf",&(physical_particle->v_y[i]));
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	for(i=0;i<physical_particle->size;i++)
		fscanf(file,"%lf",&(physical_particle->v_z[i]));
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	fgets(buf, 200, file);
	for(i=0;i<physical_particle->size;i++)
		fscanf(file,"%d",&(physical_particle->ID[i]));

	fclose(file);

	for(i=0;i<physical_particle->size;i++)
	{
		physical_particle->v_x[i]=physical_particle->v_x[i]*1.0/UnitVelocity_in_cm_per_s;
		physical_particle->v_y[i]=physical_particle->v_y[i]*1.0/UnitVelocity_in_cm_per_s;
		physical_particle->v_z[i]=physical_particle->v_z[i]*1.0/UnitVelocity_in_cm_per_s;
	}
#ifdef COSMOLOGY_UNIT
	for(i=0;i<physical_particle->size;i++)
	{
		physical_particle->v_x[i]=physical_particle->v_x[i]*1.0/(1.0+Redshift)/sqrt(1.0+Redshift);
		physical_particle->v_y[i]=physical_particle->v_y[i]*1.0/(1.0+Redshift)/sqrt(1.0+Redshift);
		physical_particle->v_z[i]=physical_particle->v_z[i]*1.0/(1.0+Redshift)/sqrt(1.0+Redshift);
	}
#endif
	return errorflag;
}
