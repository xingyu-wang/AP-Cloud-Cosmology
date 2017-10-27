/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \author Gaurish Telang <gaurish108@gaurish108> 
 * \date   Mon Jan. 30 2012
 * 
 * \brief  Contains the Octree data structure.  
 */

#ifndef __OCTREE_H__
#define __OCTREE_H__


#ifdef _MSC_VER  
typedef __int32 int32_t; 
typedef unsigned __int32 uint32_t; 
typedef __int64 int64_t; 
typedef unsigned __int64 uint64_t;  
#else 
#include <stdint.h> 
#endif 

#include <deque>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <numeric> 
#include <fstream>
#include <math.h>
#include "particle.h"

#include <omp.h>

/*!
\brief This struct contains the morton key and the index of a particle. 
       a vector of such structs is then sorted by the "key" field. 
 */
struct KeyIndex
{
//  int key;//Should be a 32/64 bit integer. With a 32 bit integer we can use Octrees of depth 10. With a 64 bit integer we can use octrees of depth 20.
  uint64_t key;
  int index;
  bool operator < (const KeyIndex& a) const
  {
    return (key < a.key);
  }
};


/*!
\brief This class does not have octree search algorithm. The search algorithm is implemented in DefaultParticles and its subclasses.
       The member fields starting from "m_v" must have same length with DefaultParticles::m_iNumOfParticles and m_iTotalNumberOfParticles. 
       That is, the length of these member field means the total number of particles.
       Each value of the array represent octree nodes.
       The location information can be fetched by decoding m_vKey value.
\see DefaultParticles
*/
class Octree {

private:

  /*!
  \brief array of x-coordinates of particles.
  */
//  const std::vector<double> m_vCoordX;
  /*!
  \brief array of y-coordinates of particles.
  */
//  const std::vector<double> m_vCoordY;
  /*!
  \brief array of z-coordinates of particles.
  */
//  const std::vector<double> m_vCoordZ;

	double *m_vCoordX,*m_vCoordY,*m_vCoordZ;

  /*!
  \brief Total number of particles. The value of this field must be same with DefaultParticles::m_iNumOfParticles.
  */
  int m_iTotalNumberOfParticles;
 
  /*!
  \brief Array of size m_iTotalNumberOfParticles which contains the Morton key and Index of each particle. 
  */
  KeyIndex *m_vParticleKeyIndex ;
                        
  /*!
  \brief The maximum depth of this octree
  \see m_vDepth
  */
  int m_iMaxDepth;

  /*!
  \brief This field is for search algorithm. Do not use this field except the implementor of octree search algorithm.
  */
  int m_iLinearSearchDepth;
  
  /*!
  \brief Minimum X-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_min_x ;

  /*!
  \brief Maximum X-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_max_x ;
 
  /*!
  \brief Minimum Y-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_min_y ;

  /*!
  \brief Maximum Y-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_max_y ;

  /*!
  \brief Minimum Z-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_min_z ;

  /*!
  \brief Maximum Z-coordinate of the bounding box in which we construct the octree
  */
  double m_dBoundingBox_max_z ;

  /*!
  \brief This is the length of the octree. Specifically it is the length of the arrays listed below
  */
  int m_iTreeLength;
  
  /*!
  \brief Morton Key. 
         This key has octal base number representation.
         So the total length having meaning is 3 * (m_iMaxDepth - 1) bits.
         The first 3 bits represent the first region divided by 2.
         For example,
         000 means x <= 1/2, y <= 1/2, z <= 1/2
         001 means x <= 1/2, y <= 1/2, z > 1/2
         010 means x <= 1/2, y > 1/2, z <= 1/2
         ...
         111 means x > 1/2, y > 1/2, z > 1/2
         The second successive three bits have the same meaning above in the region the first three bits represent.
         So the root node of octree must have 0 key value.
         In the second depth of octree, all bits execpt the first three bits have 0 values.
  */
//  int* m_vNodeKey;
  uint64_t* m_vNodeKey;
  /*!
  \brief Represent the depth of the node which is designated by the array index.
         The value of the root node is 0.
  */
  int* m_vDepth;

  /*!
  \brief The first index of particles of the nodes. In combination with m_vNumberOfParticle, we can obtain all indexes of particles in the node which is designated by the array index.
             In order to travel all particles in the node which is designated by the index, reference from m_vFirstParticleIndex to m_vFirstParticleIndex + m_vNumberOfParticle - 1.
  \see m_vNumberOfParticle
  */
  int* m_vFirstParticleIndex;
  
  /*!
  \brief The total number of particles including in the node which is designated by the array index.
  \see m_vFirstParticleIndex
  */
  int* m_vNumberOfContainedParticles;

  /*!
  \brief The first index of child node. In combination with m_vNumberOfChild, we can obtain all indexes of children nodes in the node which is designated by the array index.
         In order to travel all nodes in the node which is designated by the index, reference from m_vFirstChildIndex to m_vFirstChildIndex + m_vNumberOfChild - 1.
         If the node designated tye the array index is a leaf node, this field is set by -1.
  \see m_vNumberOfChild
  */
  int* m_vFirstChildIndex;
  
  /*!
  \brief The total number of children nodes including in the node which is designated by the array index.
         If the node designated tye the array index is a leaf node, this field is set by 0.
  \see m_vFirstChildIndex
  */
  int* m_vNumberOfChildren;

  /*!
  \brief The lower limit of x coordinates of the node which is designated by the array index.
  */
  double* m_vLowerLimitOfX;
  
  /*!
  \brief The lower limit of y coordinates of the node which is designated by the array index.
  */
  double* m_vLowerLimitOfY;
  
  /*!
  \brief The lower limit of z coordinates of the node which is designated by the array index.
  */
  double* m_vLowerLimitOfZ;

  /*!
  \brief The upper limit of x coordinates of the node which is designated by the array index.
  */
  double* m_vUpperLimitOfX;
  
  /*!
  \brief The upper limit of y coordinates of the node which is designated by the array index.
  */
  double* m_vUpperLimitOfY;
  
  /*!
  \brief The upper limit of z coordinates of the node which is designated by the array index.
  */
  double* m_vUpperLimitOfZ;

	int bc;

public:
/*  Octree(const std::vector<double> &xp,
         const std::vector<double> &yp,
         const std::vector<double> &zp,
         int treedepth, 
         int numOfParticles);*/
	Octree(Physical_particle *physical_particle, int treedepth, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max,int bc);
	~Octree();
	
  /*!
  \brief Method for building octree data structure from the positions of the particles.
  */
  int buildOctree();
  
  /*!
   \brief Given a search point, find all the neighbours within specified radius. Store all the indices of the
          particles on a deque.
  */
  int searchNeighbor(const double &search_x,
                     const double &search_y,
                     const double &search_z,
                     const double &radius,
		     std::vector<int> *result,
		     std::vector<double> *distance);
  int searchNeighbor(const double &search_x,
                     const double &search_y,
                     const double &search_z,
                     const double &radius1,
                     const double &radius2,
		     std::vector<int> *result,
		     std::vector<double> *distance);

int interpolation_from_node_to_particle_2to1_3d(Computational_particle *cp, Physical_particle *pp);
int select_computational_particles_2to1_3d(Computational_particle *computational_particle, double c, int order, int nmax, int mindepth,int maxdepth,int bc);
double periodic_wrap(double x, int dimension);

private:
    /*!
        \brief This method computes the Morton Key for a particle whose coordinates are x,y and z.
               Used for spatially sorting the particles.
     */
  inline uint32_t computeKey(const double& x, 
                             const double& y, 
                             const double& z) ;
  inline uint64_t computeKey64(const double& x, 
                             const double& y, 
                             const double& z) ;
};
#endif // __OCTREE_H__                                                
