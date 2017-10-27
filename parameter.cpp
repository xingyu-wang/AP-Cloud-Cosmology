#include "parameter.h"

double Omega;
double OmegaLambda;
double Redshift;
double HubbleParam;
double GravitionalConstant;
double UnitLength_in_cm;
double UnitMass_in_g;
double UnitVelocity_in_cm_per_s;
char ICFile[200], FileBase[200], OutputFile[200];
int ICFileType;
double DeltaA;
double EndRedShift;
int MaxOctreeLevel;
double CellSizeParameter;
int NumberOfNeighbors;
int DensityEstimator;
double BoundingBox_min_x;
double BoundingBox_max_x;
double BoundingBox_min_y;
double BoundingBox_max_y;
double BoundingBox_min_z;
double BoundingBox_max_z;
int Periodic;
double MassPerParticle;

void read_parameterfile(char *fname)
{
#define FLOAT 1
#define STRING 2
#define INT 3
#define MAXTAGS 300

  FILE *fd;
  char buf[200], buf1[200], buf2[200], buf3[200];
  int i, j, nt;
  int id[MAXTAGS];
  void *addr[MAXTAGS];
  char tag[MAXTAGS][50];
  int errorFlag = 0;

  /* read parameter file on all processes for simplicty */

  nt = 0;

  strcpy(tag[nt], "Omega");
  addr[nt] = &Omega;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "OmegaLambda");
  addr[nt] = &OmegaLambda;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "HubbleParam");
  addr[nt] = &HubbleParam;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Redshift");
  addr[nt] = &Redshift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "GravitionalConstant");
  addr[nt] = &GravitionalConstant;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "ICFile");
  addr[nt] = ICFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "ICFileType");
  addr[nt] = &ICFileType;
  id[nt++] = INT;

  strcpy(tag[nt], "OutputFile");
  addr[nt] = OutputFile;
  id[nt++] = STRING;

  strcpy(tag[nt], "FileBase");
  addr[nt] = FileBase;
  id[nt++] = STRING;

  strcpy(tag[nt], "UnitVelocity_in_cm_per_s");
  addr[nt] = &UnitVelocity_in_cm_per_s;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitLength_in_cm");
  addr[nt] = &UnitLength_in_cm;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "UnitMass_in_g");
  addr[nt] = &UnitMass_in_g;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "DeltaA");
  addr[nt] = &DeltaA;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "EndRedShift");
  addr[nt] = &EndRedShift;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "MaxOctreeLevel");
  addr[nt] = &MaxOctreeLevel;
  id[nt++] = INT;

  strcpy(tag[nt], "CellSizeParameter");
  addr[nt] = &CellSizeParameter;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "NumberOfNeighbors");
  addr[nt] = &NumberOfNeighbors;
  id[nt++] = INT;

  strcpy(tag[nt], "DensityEstimator");
  addr[nt] = &DensityEstimator;
  id[nt++] = INT;

  strcpy(tag[nt], "BoundingBox_min_x");
  addr[nt] = &BoundingBox_min_x;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoundingBox_max_x");
  addr[nt] = &BoundingBox_max_x;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoundingBox_min_y");
  addr[nt] = &BoundingBox_min_y;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoundingBox_max_y");
  addr[nt] = &BoundingBox_max_y;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoundingBox_min_z");
  addr[nt] = &BoundingBox_min_z;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "BoundingBox_max_z");
  addr[nt] = &BoundingBox_max_z;
  id[nt++] = FLOAT;

  strcpy(tag[nt], "Periodic");
  addr[nt] = &Periodic;
  id[nt++] = INT;

  strcpy(tag[nt], "MassPerParticle");
  addr[nt] = &MassPerParticle;
  id[nt++] = FLOAT;

  if((fd = fopen(fname, "r")))
    {
      while(!feof(fd))
	{
	  buf[0] = 0;
	  fgets(buf, 200, fd);

	  if(sscanf(buf, "%s%s%s", buf1, buf2, buf3) < 2)
	    continue;

	  if(buf1[0] == '%')
	    continue;

	  for(i = 0, j = -1; i < nt; i++)
	    if(strcmp(buf1, tag[i]) == 0)
	      {
		j = i;
		tag[i][0] = 0;
		break;
	      }

	  if(j >= 0)
	    {
	      switch (id[j])
		{
		case FLOAT:
		  *((double *) addr[j]) = atof(buf2);
		  break;
		case STRING:
		  strcpy((char *)addr[j], buf2);
		  break;
		case INT:
		  *((int *) addr[j]) = atoi(buf2);
		  break;
		}
	    }
	  else
	    {
		fprintf(stdout, "Error in file %s:   Tag '%s' not allowed or multiple defined.\n", fname,
			buf1);
	      errorFlag = 1;
	    }
	}
      fclose(fd);

    }
  else
    {
	fprintf(stdout, "Parameter file %s not found.\n", fname);
      errorFlag = 1;
    }


  for(i = 0; i < nt; i++)
    {
      if(*tag[i])
	{
	    fprintf(stdout, "Error. I miss a value for tag '%s' in parameter file '%s'.\n", tag[i], fname);
	  errorFlag = 1;
	}
    }

  if(errorFlag)
    {
      exit(0);
    }

#undef FLOAT
#undef STRING
#undef INT
#undef MAXTAGS
}
