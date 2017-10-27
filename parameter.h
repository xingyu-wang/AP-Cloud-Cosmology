#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

void read_parameterfile(char *fname);

extern double Omega;
extern double OmegaLambda;
extern double Redshift;
extern double HubbleParam;
extern double GravitionalConstant;
extern double UnitLength_in_cm;
extern double UnitMass_in_g;
extern double UnitVelocity_in_cm_per_s;
extern char ICFile[200], FileBase[200], OutputFile[200];
extern int ICFileType;
extern double DeltaA;
extern double EndRedShift;
extern int MaxOctreeLevel;
extern double CellSizeParameter;
extern int NumberOfNeighbors;
extern int DensityEstimator;
extern double BoundingBox_min_x;
extern double BoundingBox_max_x;
extern double BoundingBox_min_y;
extern double BoundingBox_max_y;
extern double BoundingBox_min_z;
extern double BoundingBox_max_z;
extern int Periodic;
extern double MassPerParticle;
