/*!
 * \author Yu, Kwangmin <yukwangmin@gmail.com> 
 * \author Gaurish Telang <gaurish108@gaurish108> 
 * \date   Mon Jan. 30 2012
 * 
 * \brief Implemntations of the methods and helper functions mentioned in the Octree class.
 */
#include "octree.h"
#include "apcloud_macro.h"

//#define OUTPUT_TIME

inline bool is_node_intersect_search_region(const double min_x,
                                            const double max_x,
                                            const double min_y,
                                            const double max_y,
                                            const double min_z,
                                            const double max_z,
                                            const double& search_x, 
                                            const double& search_y, 
                                            const double& search_z, const double& radius) 
{
  //The following calculation has been fashioned such that:
  //if final value of squared_dmin == 0 then the point lies inside the node. 
  //if final value of squared_dmin !=0  then the point lies outside the node, AND 
  //             tells the SQUARE of the minimum distance of the point to the points on the node boundary(surface).
  float squared_dmin=0;
  float temp; //Used as a temporary variable to store some intermediate results.
  //Process the x cooridinates
  if( search_x < min_x ) 
    { 
       temp = search_x - min_x; 
       squared_dmin += temp*temp;
    }  
   else if( search_x > max_x )	
    { 
       temp         = search_x - max_x ; 
       squared_dmin += temp*temp                  ;
    }   
 

  //Process the Y-coorindtaes
  if( search_y < min_y )  
  { 
    temp = search_y - min_y ; 
    squared_dmin += temp*temp;
  } 
   else if( search_y > max_y ) 
   { 
    temp = search_y - max_y ; 
    squared_dmin += temp*temp;
   }   



  //Process the Z-coorindtaes
  if( search_z < min_z ) 
  { 
    temp          = search_z - min_z; 
    squared_dmin += temp*temp;
  } 
  
  else if( search_z > max_z ) 
   { 
    temp          = search_z - max_z; 
    squared_dmin += temp*temp;
   }   


  // Based on the calculated value give a  YES/NO answer to the question
  if (squared_dmin <= radius*radius) 
    {
       return true;
    }
  else 
    { 
       return false;
    }
}

Octree::~Octree() {
	delete[] m_vParticleKeyIndex;
	delete[] m_vNodeKey;
	delete[] m_vDepth;
	delete[] m_vFirstParticleIndex;
	delete[] m_vNumberOfContainedParticles;
	delete[] m_vFirstChildIndex;
	delete[] m_vNumberOfChildren;
	delete[] m_vLowerLimitOfX;
	delete[] m_vUpperLimitOfX;
	delete[] m_vLowerLimitOfY;
	delete[] m_vUpperLimitOfY;
	delete[] m_vLowerLimitOfZ;
	delete[] m_vUpperLimitOfZ;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
/*Octree::Octree(const std::vector<double> &xp,
               const std::vector<double> &yp,
               const std::vector<double> &zp,
               int treedepth, 
               int numOfParticles) :m_vCoordX(xp), m_vCoordY(yp), m_vCoordZ(zp), m_iTotalNumberOfParticles(numOfParticles), m_iMaxDepth(treedepth) {
	buildOctree();
}*/
Octree::Octree(Physical_particle *physical_particle, int treedepth, double x_min, double x_max, double y_min, double y_max, double z_min, double z_max, int BC) :m_vCoordX(physical_particle->x), m_vCoordY(physical_particle->y), m_vCoordZ(physical_particle->z), m_iTotalNumberOfParticles(physical_particle->size), m_iMaxDepth(treedepth), m_dBoundingBox_min_x(x_min), m_dBoundingBox_max_x(x_max), m_dBoundingBox_min_y(y_min), m_dBoundingBox_max_y(y_max), m_dBoundingBox_min_z(z_min), m_dBoundingBox_max_z(z_max), bc(BC){
	buildOctree();
}


int Octree::buildOctree() {

  //Calculate the bounding box limits
/*
   m_dBoundingBox_min_x = *std::min_element(m_vCoordX.begin(), m_vCoordX.end());
   m_dBoundingBox_max_x = *std::max_element(m_vCoordX.begin(), m_vCoordX.end());

   m_dBoundingBox_min_y = *std::min_element(m_vCoordY.begin(), m_vCoordY.end());
   m_dBoundingBox_max_y = *std::max_element(m_vCoordY.begin(), m_vCoordY.end());

   m_dBoundingBox_min_z = *std::min_element(m_vCoordZ.begin(), m_vCoordZ.end());
   m_dBoundingBox_max_z = *std::max_element(m_vCoordZ.begin(), m_vCoordZ.end());
*/
   //Create an array which will hold the morton key and index of each particle.
   //After we fill this array, it will be sorted by the "key" field in the "KeyIndex" struct.
   //The operator < has been overloaded for this purpose
   m_vParticleKeyIndex = new KeyIndex [ m_iTotalNumberOfParticles ];

#ifdef OUTPUT_TIME
	double time1;
  time1 = omp_get_wtime();
#endif
  #pragma omp parallel for
  for (int i = 0; i < m_iTotalNumberOfParticles; ++i)
    {
//      m_vParticleKeyIndex[i].key   = computeKey( m_vCoordX[i], m_vCoordY[i], m_vCoordZ[i] );//function call which computes the morton key from the coordinate arrays.defined below.
      m_vParticleKeyIndex[i].key   = computeKey64( m_vCoordX[i], m_vCoordY[i], m_vCoordZ[i] );//function call which computes the morton key from the coordinate arrays.defined below.
      m_vParticleKeyIndex[i].index = i;
    }
#ifdef OUTPUT_TIME
	time1= omp_get_wtime()-time1;
	std::cout << "Key Computation Time : " << time1 << std::endl;
#endif
  //Now sort the key-index array
#ifdef OUTPUT_TIME
  time1 = omp_get_wtime();
#endif
  std::sort(m_vParticleKeyIndex, m_vParticleKeyIndex + m_iTotalNumberOfParticles);
#ifdef OUTPUT_TIME
	time1= omp_get_wtime()-time1;
	std::cout << "Key Sorting Time : " << time1 << std::endl;
#endif
  //Now that we have the sorted key-index array we can begin the tree construction level by level
  //For that we need to calculate three helper vectors. bitmasks, baseaddresses and numnodes. 
  //numnodes[i]     is the number of nodes at level i of the octree. zero level correposnds to root nodes. hence numnodes[0] = 1
  //basaddresses[i] is the position in the octree arrays of first node at level i
  //bitmasks[i]     is the Bit-Mask for level i. used in the computation of the numnodes vector 
  std::vector<uint64_t> bitmasks(m_iMaxDepth+1);
  std::vector<int> baseaddresses(m_iMaxDepth+1);
  std::vector<int> numnodes(m_iMaxDepth+1);

  //CALCULATE THE BITMASKS VECTOR
   bitmasks[m_iMaxDepth]=1;
   bitmasks[m_iMaxDepth] = ( bitmasks[m_iMaxDepth] << (3 * m_iMaxDepth) ) - 1 ;//Deepest level bitmask 
   bitmasks[    0      ] = 0                       ;//Root    level bitmask 
   //Compute remaining level bitmasks
   for (int k = m_iMaxDepth -1 ; k >= 1 ; --k)
    {
      int shiftbits = 3 * ( m_iMaxDepth - k ) ;
      bitmasks[k] = ( bitmasks[ m_iMaxDepth ] >> shiftbits ) << shiftbits;  
    } 

#ifdef OUTPUT_TIME
  time1 = omp_get_wtime();
#endif
  //COMPUTE THE NUMNODES VECTOR
   numnodes[ 0 ] =  1 ;//Root level has only 1 node i.e. the root node itself.
  for (int k = 1; k < m_iMaxDepth + 1; ++k)
   {
    int count=1;
    
		#pragma omp parallel for default(shared) reduction(+:count)
     for ( int i = 1; i < m_iTotalNumberOfParticles ; ++i)
     {
       if ( (m_vParticleKeyIndex[i].key & bitmasks[k]) != (m_vParticleKeyIndex[i-1].key & bitmasks[k] ))
	 {
	   ++count;
	 }
     }
     numnodes[k]=count;   
   }
#ifdef OUTPUT_TIME
	time1= omp_get_wtime()-time1;
	std::cout << "NUMNODES Computation Time : " << time1 << std::endl;
#endif

  //COMPUTE THE BASEADDRESSES VECTOR.
  baseaddresses[0]=0; 
  baseaddresses[1]=1;
  for (unsigned int k = 2; k < baseaddresses.size(); ++k)
   {
     baseaddresses[k]=baseaddresses[k-1]+numnodes[k-1];
   }   

  //With the bitmasks, numnodes and baseaddresses vectors we can now begin the tree computation.
  
  //This is the length of the arrays used in descirbing octrees. 
 m_iTreeLength = std::accumulate(numnodes.begin(),numnodes.end(),0);

  //Now we allocate the arrays using the new operator and m_iTreeLength
  m_vNodeKey			= new uint64_t[m_iTreeLength];
  m_vDepth			= new int[m_iTreeLength];
  m_vFirstParticleIndex		= new int[m_iTreeLength];
  m_vNumberOfContainedParticles = new int[m_iTreeLength];
  m_vFirstChildIndex		= new int[m_iTreeLength];
  m_vNumberOfChildren		= new int[m_iTreeLength];  

  m_vLowerLimitOfX		= new double[m_iTreeLength];     
  m_vUpperLimitOfX		= new double[m_iTreeLength]; 

  m_vLowerLimitOfY		= new double[m_iTreeLength]; 
  m_vUpperLimitOfY		= new double[m_iTreeLength]; 

  m_vLowerLimitOfZ		= new double[m_iTreeLength]; 
  m_vUpperLimitOfZ		= new double[m_iTreeLength]; 
#ifdef OUTPUT_TIME
  time1 = omp_get_wtime();
#endif
  //Start building the last level
  
  //Construct the deepest level of the tree first.
  m_vNodeKey[baseaddresses[m_iMaxDepth]]   = m_vParticleKeyIndex[0].key;
  m_vDepth[baseaddresses[m_iMaxDepth]]     = m_iMaxDepth;

  m_vFirstParticleIndex[baseaddresses[m_iMaxDepth]]	= 0;
 //tree[baseaddressed[m_iMaxDepth]].pnum calculated in the for loop below.

  m_vNumberOfChildren[baseaddresses[m_iMaxDepth]]	= 0;
  m_vFirstChildIndex[baseaddresses[m_iMaxDepth]]	= (-1);
  

  int j   = 1; //This counter changes on encountering a new key.  
  int fix = 0; //Used in the pnum calculation. Records the index where the key change last took place.  
               //pnum = i-fix where i and fix are positions where successive key changes take place.  

   for (int i = 1; i < m_iTotalNumberOfParticles ; ++i)
    {
       if(m_vParticleKeyIndex[i].key!=m_vParticleKeyIndex[i-1].key)
  	{
	     m_vNodeKey[baseaddresses[m_iMaxDepth]+j]				=m_vParticleKeyIndex[i].key;  
     
	     m_vDepth[baseaddresses[m_iMaxDepth]+j]				=m_iMaxDepth; 

	     m_vFirstParticleIndex[baseaddresses[m_iMaxDepth]+j]		=i         ; 
  
       	     m_vNumberOfChildren[baseaddresses[m_iMaxDepth]+j]  		=0         ; //No child nodes. Deepest level.

	     m_vFirstChildIndex[baseaddresses[m_iMaxDepth]+j]   		=(-1);//Special sentinel value. No children since leaf node.    
                                                       
	     m_vNumberOfContainedParticles[baseaddresses[m_iMaxDepth]+j-1]	=i-fix; //Justified above!
             ++j;  

             fix								=i;//change origin of the measure tape. 
  	}//END computation on encountering a new label

    }
   m_vNumberOfContainedParticles[baseaddresses[m_iMaxDepth]+j-1]  = m_iTotalNumberOfParticles - fix;//This cannot be tackled within above for loop.
  //END DEEPEST LEVEL NODE CONSTRUCTION.
#ifdef OUTPUT_TIME
	time1= omp_get_wtime()-time1;
	std::cout << "Leaf Construction Time : " << time1 << std::endl;
  time1 = omp_get_wtime();
#endif
  //Now that the deepest level has been constructed lets us start filling in the rest of the octree nodes.
  //Fill the other levels.
   for (int level=m_iMaxDepth-1; level>=0; --level)
    {
      j=1;
      fix=0;    //used in cnum calculation. Records the index where key change LAST took place. 
                //cnum=i-fix. 'i' and 'fix' are positions where successive MASKED key changes take place
            
      //LEVEL CONSTRUCTION BEGINS. WE MUST SCAN LEVEL+1 SECTION TO INSERT APPROPRIATE VALUES IN THE LEVEL SECTION.
      m_vNodeKey[baseaddresses[level]]			= ( m_vNodeKey[baseaddresses[level+1]] & bitmasks[level])     ;
      m_vDepth[baseaddresses[level]]			= level;
      m_vFirstParticleIndex[baseaddresses[level]]	= 0                         ;
      m_vFirstChildIndex[baseaddresses[level]]		= baseaddresses[level+1]     ;
 
      int sum=m_vNumberOfContainedParticles[baseaddresses[level+1]];//This variable used in pnum calculations.
      m_vNumberOfContainedParticles[baseaddresses[level]]	= sum  ; 

      for (int i = 1; i < numnodes[level+1] ; ++i)
  	{
     	  if( ( m_vNodeKey[baseaddresses[level+1]+i] & bitmasks[level] ) != ( m_vNodeKey[baseaddresses[level+1] + i-1] & bitmasks[level] )  )//Relate i and i-1
  	    {
  	       m_vNodeKey[baseaddresses[level]+j]			= ( m_vNodeKey[baseaddresses[level+1]+i] & bitmasks[level] )  ; //Give it the bitmasked key.  
	       m_vDepth[baseaddresses[level]+j]				= level                       ;//Assign the depth
	       m_vFirstParticleIndex[baseaddresses[level]+j]		= m_vFirstParticleIndex[baseaddresses[level+1]+i]                              ;
	       m_vNumberOfContainedParticles[baseaddresses[level]+j-1]	= sum;//Calculate the difference where succesive changes take place in keys masked with a level mask.
	       m_vFirstChildIndex[baseaddresses[level]+j]		= baseaddresses[level+1]+i        ;//Since we encountered a new masked level key. we record the place where this is encountered/
	       m_vNumberOfChildren[baseaddresses[level]+j-1]		= i-fix                           ;//Will have to keep some kind of running sum
	       
  	      ++j  ;
  	      fix=i;
	      sum = m_vNumberOfContainedParticles[baseaddresses[level+1]+i];//initializing sum to pnum of the new key encountered.              
            }

      	  else sum += m_vNumberOfContainedParticles[baseaddresses[level+1]+i];//This is executed only if baseaddresses[level+1]+i and baseaddresses[level+1]+i-1 share the same parent . 
	                                                  //Hence we increment the sum counter which is keeping track of the number of particles within the parent node.
     	}//end inner for 

       m_vNumberOfChildren[baseaddresses[level]+j-1]                = numnodes[level+1]-fix   ;
       m_vNumberOfContainedParticles[baseaddresses[level]+j-1]      = sum                     ;
    }//end outer for

//END LEVEL CONSTRUCTION.
#ifdef OUTPUT_TIME
	time1= omp_get_wtime()-time1;
	std::cout << "Non-leaf Construction Time : " << time1 << std::endl;
#endif
  //Finally we fill in the vertex information required by all the nodes.
  //This is a separate computatin because of its length.
  //Loop over all the tree vector from left to right
  //for center & vertex computation. The only information necessary 
  //for the entire computation here is the DEPTH and the KEY of all the octree nodes.
  //These have already been computed above.

	//In the previous code, the vertex location of each node is computed as the limit of bounding box of the domain plus a series of displacement. Each level higher than the current level offers a displacement for each node, so the time complexity is O(m_iMaxDepth*m_iTreeLength). For some problem, the time to compute bounding box can take ~50% of the total time to build octree. In this code, the vertex location of each node is computed as the vertex location of its parent plus only one displacement given by the current level, so the new time complexity is O(m_iTreeLength). This modification reduces the time to compute vertex location by ~80%, compared to the previous code.
#ifdef OUTPUT_TIME
	time1= omp_get_wtime();
#endif
    double xwidth =  m_dBoundingBox_max_x - m_dBoundingBox_min_x   ;
    double ywidth =  m_dBoundingBox_max_y - m_dBoundingBox_min_y   ;
    double zwidth =  m_dBoundingBox_max_z - m_dBoundingBox_min_z   ;

    m_vLowerLimitOfX[0] = m_dBoundingBox_min_x; 
    m_vLowerLimitOfY[0] = m_dBoundingBox_min_y;
    m_vLowerLimitOfZ[0] = m_dBoundingBox_min_z; 

    m_vUpperLimitOfX[0] = m_dBoundingBox_max_x; 
    m_vUpperLimitOfY[0] = m_dBoundingBox_max_y;
    m_vUpperLimitOfZ[0] = m_dBoundingBox_max_z; 

	for(int level=0; level<m_iMaxDepth; ++level)
	{
		xwidth=0.5*xwidth;
		ywidth=0.5*ywidth;
		zwidth=0.5*zwidth;
    #pragma omp parallel for
		for(int i=0; i < numnodes[level] ; ++i)
		{
			for(int j=m_vFirstChildIndex[baseaddresses[level]+i]; j<m_vFirstChildIndex[baseaddresses[level]+i]+m_vNumberOfChildren[baseaddresses[level]+i]; ++j)
			{
				uint64_t bit_triplet_extracted = ( m_vNodeKey[j] >> (3*(m_iMaxDepth-level-1)) ) & 0x7;
      	//Case 000
        if (bit_triplet_extracted ==       0)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
    	  }
        //Case 001
        else if (bit_triplet_extracted  == 1)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
	  }
         //Case 010
        else if (bit_triplet_extracted  == 2)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
	  }
         //Case 011
        else if (bit_triplet_extracted  == 3)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i];
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i]-xwidth;
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
	  }
          //Case 100
        else if (bit_triplet_extracted ==  4)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
	  } 
        //Case 101
        else if (bit_triplet_extracted  == 5)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i];
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i]-ywidth;
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
	  }
          //Case 110 
        else if (bit_triplet_extracted  == 6)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i];
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i]-zwidth;
	  }
         //Case 111
        else if (bit_triplet_extracted  == 7)
	  {
            m_vLowerLimitOfX[j]=m_vLowerLimitOfX[baseaddresses[level]+i]+xwidth;
            m_vLowerLimitOfY[j]=m_vLowerLimitOfY[baseaddresses[level]+i]+ywidth;
            m_vLowerLimitOfZ[j]=m_vLowerLimitOfZ[baseaddresses[level]+i]+zwidth;
            m_vUpperLimitOfX[j]=m_vUpperLimitOfX[baseaddresses[level]+i];
            m_vUpperLimitOfY[j]=m_vUpperLimitOfY[baseaddresses[level]+i];
            m_vUpperLimitOfZ[j]=m_vUpperLimitOfZ[baseaddresses[level]+i];
	  }	
			
			}//j
		}//i
	}//level
#ifdef OUTPUT_TIME
	time1= omp_get_wtime()-time1;
	std::cout << "Bounding box Construction Time : " << time1 << std::endl;
#endif
  //Print the octree information to a file. 

  // std::ofstream octree_file;
  // octree_file.open("printed_octree_vector.txt")			;
  // octree_file << "BOUNDING_BOX INFORMATION :" << std::endl	;
  // octree_file << "Min_x:  "  << m_dBoundingBox_min_x   << std::endl			;
  // octree_file << "Max_x:  "  << m_dBoundingBox_max_x<< std::endl			;

  // octree_file << "\nMin_y:  " << m_dBoundingBox_min_y<< std::endl			;
  // octree_file << "Max_y:  "   << m_dBoundingBox_max_y<< std::endl			;

  // octree_file << "\nMin_z:  "  << m_dBoundingBox_min_z<< std::endl			;
  // octree_file << "Max_z:  "   << m_dBoundingBox_max_z<< std::endl			;
  
  // octree_file << "\nLEVEL\t"<<"BITMASK\t"<<"BASEADDRESSES\t"<<"NUMNODES"<<std::endl;
  // for (int i = 0; i < m_iMaxDepth+1; ++i)
  //   {
  //     octree_file << i << "\t" << bitmasks[i] << "\t" << baseaddresses[i] << "\t\t" << numnodes[i] << std::endl;
  //   }
  //  octree_file << "Length of the octree vector is m_iTreeLength = " << m_iTreeLength << std::endl;  
  // for (int level = 0; level <= m_iMaxDepth; ++level)
  // {
  //   octree_file<<"============================================================================================================================================================"<<std::endl;
  //   octree_file<<"RANK\t"<<"DEPTH\t"<< "KEY\t"<<"PIDX\t"<<"PNUM\t # "<<"CIDX\t"<<"CNUM\t"<<"GLOBAL POSN\t\t"<<"VMIN[0]\tVMIN[1]\tVMIN[2]\t\t\tVMAX[0]\tVMAX[1]\tVMAX[2]"<<std::endl;
  //   octree_file<<"============================================================================================================================================================"<<std::endl; 
  //    for (int i = 0; i < numnodes[level]; ++i)
  // 	{
  // 	  octree_file<<std::setprecision(5)<< i      <<"\t"
  //             << m_vDepth[baseaddresses[level]+i]    <<"\t"  
  // 	      << m_vNodeKey[baseaddresses[level]+i]      <<"\t"
  // 	      << m_vFirstParticleIndex[baseaddresses[level]+i] <<"\t"
  // 	      << m_vNumberOfContainedParticles[baseaddresses[level]+i] <<"\t # "
  // 	      << m_vFirstChildIndex[baseaddresses[level]+i] <<"\t"
  // 	      << m_vNumberOfChildren[baseaddresses[level]+i] <<"\t\t"
  // 		     <<        baseaddresses[level]+i     <<"\t\t" 

  // 		   << m_vLowerLimitOfX[baseaddresses[level]+i] << "\t"
  // 		   << m_vLowerLimitOfY[baseaddresses[level]+i] << "\t"
  // 		   << m_vLowerLimitOfZ[baseaddresses[level]+i] << "\t\t\t"

  // 		   << m_vUpperLimitOfX[baseaddresses[level]+i] << "\t"
  // 		   << m_vUpperLimitOfY[baseaddresses[level]+i] << "\t"
  // 		   << m_vUpperLimitOfZ[baseaddresses[level]+i] << "\t" <<std::endl;

  //      	}
  //     octree_file<<'\n';
  // }



  return 0;

}//End of the buildOctree method
 
////////////////////////////////////////////////////////////////////////////////////////////////////

inline uint32_t Octree::computeKey(const double& x, const double& y, const double& z)
{
  uint32_t    	key=0;//Key which will be returned. 
                      //Excluding this bit, the key can containg upto 30 bits 
                      //corresponding to a max-tree depth of 10. 

  double	left_x=m_dBoundingBox_min_x ,right_x=m_dBoundingBox_max_x ;
  double	left_y=m_dBoundingBox_min_y ,right_y=m_dBoundingBox_max_y ;
  double	left_z=m_dBoundingBox_min_z ,right_z=m_dBoundingBox_max_z ;

  //Midpoint of the various intervals in the following for loop;
  double	midpt_x;
  double        midpt_y;
  double        midpt_z;

   for (int i = 1; i <= m_iMaxDepth; ++i)
    {
      //Compute midpoints.
      midpt_x=(left_x+right_x)/2.0;
      midpt_y=(left_y+right_y)/2.0;
      midpt_z=(left_z+right_z)/2.0;

      //Now we consider 8 cases. EXACTLY ONE of the following will be true within this for loop. 
      //Hence we place it within an if-elsef-if condtruct which guarantees mutual exclusion
      
      //left_x, left_y et al. will get modified within these if-else constructs 
      //since the boundary of the containing interval is continously shrinking.

      //Case 1.
     if (x<midpt_x && \
         y<midpt_y && \
         z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc000 
 	 key=(key << 3);//insert three zeros at the end
 

         right_x=midpt_x;
       

         right_y=midpt_y;


         right_z=midpt_z;
       }
   

      //Case 2.
     else if (x>=midpt_x &&   \
              y<midpt_y &&        \
              z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc100 
 	 key=(key << 3) + (1<<2);
 
         left_x=midpt_x;

       

         right_y=midpt_y;


         right_z=midpt_z;
       }

     //Case 3
        else if (x<midpt_x &&   \
                 y>=midpt_y &&        \
                 z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc010 
 	 key=(key << 3) + (1<<1);
 

         right_x=midpt_x;
       
         left_y =midpt_y;



         right_z=midpt_z;
       }

     
     //Case 4
      else if ( x<midpt_x &&   \
                y<midpt_y &&        \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc001 
 	 key=(key << 3) + 1;
 

         right_x=midpt_x;
       

         right_y=midpt_y;

         left_z=midpt_z;

       }

     //Case 5
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc110 
 	 key=(key << 3) + (1<<2) + (1<<1) ;
 
         left_x=midpt_x;

       
         left_y=midpt_y;



         right_z=midpt_z;
       }

     //Case 6
      else if ( x>=midpt_x &&   \
                y<midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc101 
 	 key=(key << 3) + (1<<2) + 1 ;
 
         left_x =midpt_x;

       

         right_y=midpt_y;

         left_z=midpt_z;

       }

     //Case 7
      else if ( x<midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc011 
 	 key=(key << 3) + (1<<1) + 1 ;
 

         right_x=midpt_x;
       
         left_y =midpt_y;


         left_z =midpt_z;

       }

     //Case 8
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc111 
 	 key=(key << 3) + (1<<2) + (1<<1) +1 ;
 
         left_x=midpt_x;

       
         left_y =midpt_y;


         left_z=midpt_z;

       }

    }//End for loop
  return key ;
}

inline uint64_t Octree::computeKey64(const double& x, const double& y, const double& z)
{
  uint64_t    	key=0;//Key which will be returned. 
                      //Excluding this bit, the key can containg upto 30 bits 
                      //corresponding to a max-tree depth of 10. 

  double	left_x=m_dBoundingBox_min_x ,right_x=m_dBoundingBox_max_x ;
  double	left_y=m_dBoundingBox_min_y ,right_y=m_dBoundingBox_max_y ;
  double	left_z=m_dBoundingBox_min_z ,right_z=m_dBoundingBox_max_z ;

  //Midpoint of the various intervals in the following for loop;
  double	midpt_x;
  double        midpt_y;
  double        midpt_z;

   for (int i = 1; i <= m_iMaxDepth; ++i)
    {
      //Compute midpoints.
      midpt_x=(left_x+right_x)/2.0;
      midpt_y=(left_y+right_y)/2.0;
      midpt_z=(left_z+right_z)/2.0;

      //Now we consider 8 cases. EXACTLY ONE of the following will be true within this for loop. 
      //Hence we place it within an if-elsef-if condtruct which guarantees mutual exclusion
      
      //left_x, left_y et al. will get modified within these if-else constructs 
      //since the boundary of the containing interval is continously shrinking.

      //Case 1.
     if (x<midpt_x && \
         y<midpt_y && \
         z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc000 
 	 key=(key << 3);//insert three zeros at the end
 

         right_x=midpt_x;
       

         right_y=midpt_y;


         right_z=midpt_z;
       }
   

      //Case 2.
     else if (x>=midpt_x &&   \
              y<midpt_y &&        \
              z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc100 
 	 key=(key << 3) + (1<<2);
 
         left_x=midpt_x;

       

         right_y=midpt_y;


         right_z=midpt_z;
       }

     //Case 3
        else if (x<midpt_x &&   \
                 y>=midpt_y &&        \
                 z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc010 
 	 key=(key << 3) + (1<<1);
 

         right_x=midpt_x;
       
         left_y =midpt_y;



         right_z=midpt_z;
       }

     
     //Case 4
      else if ( x<midpt_x &&   \
                y<midpt_y &&        \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc001 
 	 key=(key << 3) + 1;
 

         right_x=midpt_x;
       

         right_y=midpt_y;

         left_z=midpt_z;

       }

     //Case 5
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z<midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc110 
 	 key=(key << 3) + (1<<2) + (1<<1) ;
 
         left_x=midpt_x;

       
         left_y=midpt_y;



         right_z=midpt_z;
       }

     //Case 6
      else if ( x>=midpt_x &&   \
                y<midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc101 
 	 key=(key << 3) + (1<<2) + 1 ;
 
         left_x =midpt_x;

       

         right_y=midpt_y;

         left_z=midpt_z;

       }

     //Case 7
      else if ( x<midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc011 
 	 key=(key << 3) + (1<<1) + 1 ;
 

         right_x=midpt_x;
       
         left_y =midpt_y;


         left_z =midpt_z;

       }

     //Case 8
      else if ( x>=midpt_x &&   \
                y>=midpt_y &&   \
                z>=midpt_z )
       {
         // -----------------------------abc becomes --------------------------abc111 
 	 key=(key << 3) + (1<<2) + (1<<1) +1 ;
 
         left_x=midpt_x;

       
         left_y =midpt_y;


         left_z=midpt_z;

       }

    }//End for loop
  return key ;
}

////////////////////////////////////////////////////////////////////////////////////////////////////
int Octree::searchNeighbor(const double &search_x,
                           const double &search_y,
                           const double &search_z,
                           const double &radius,
			   std::vector<int> *result,
			   std::vector<double> *distance) 
{
result->clear();
distance->clear();
/*printf("%d %d \n",m_iTotalNumberOfParticles, m_iMaxDepth);
printf("%f %f %f %f \n",search_x,search_y,search_z,radius);*/
   std::deque<int> nodelist;
   nodelist.push_back(0); // Push back the index of the root node, which is at position 0 in the octree vector.

   // In this while loop, we find all the leaf nodes of the octree which intersect the search region.  
   // If, at the end of the while loop, nodelist queue is empty, that means no leaf node intersects the search region and hence the search point has no particles in it neighbourhood.
   // In this we don't check whether an octree node is a subset of the search region. 
   // Paradoxically, this seems to be performing better than an algorithm which does. This needs some more investigating. 
   // Probably for deeper octrees and smaller search radii. that algorithm will endup outperforming this one. Also that means the subset check would need to be improved upon, since it could be the bottleneck. 
   while( 1 ) 
   {

     //std::cout << nodelist.size() << std::endl;
     // If the node-list is empty stop the execution of he current function, and return 0, since 0 neighbours have been founf. 
     // Also as a result, you see that nothing is in the deque pointed to be result.
    if( nodelist.empty() ) 
      { 
/*printf("Oh no..\n");*/
        return 0;
      }
    // Get the node number at the front of the queue.
    // The popping will be done at the end of the while loop.
    int inode       = nodelist[0];
    int first_index = m_vFirstChildIndex[inode];

    // If the head node of the deque is a leaf that means all the remaining nodes in the deque are also leaves. 
    // So we can now exit the while loop and process look through all the particles within the leaf nodes. 
    if( first_index == -1 ) 
      {
/*printf("1! ");*/
        break;
      }

    int last_index = first_index + m_vNumberOfChildren[ inode ];
    //  Find all the child nodes of the CURRENT node which intersect the sphere having radius "radius".
    //  After finding push_back their index onto the stack.
/*printf("%d ",last_index-first_index);*/
    for(int i=first_index ; i < last_index ; i++) 
     {
       if(   is_node_intersect_search_region( m_vLowerLimitOfX[i],\
                                              m_vUpperLimitOfX[i],\
                                              m_vLowerLimitOfY[i],\
                                              m_vUpperLimitOfY[i],\
                                              m_vLowerLimitOfZ[i],\
                                              m_vUpperLimitOfZ[i],\
                                              search_x,\
                                              search_y,\
                                              search_z, radius) )
	 {
             nodelist.push_back( i );
/*printf("2! ");*/
	 }
    }
    nodelist.pop_front( );
  }
  // When the  above while-loop ends normally, the queue consists only of leaf nodes intersecting the search region.
  // We now look at every particle in each leaf node looking for neighbouring particles. Probably as an optimization step,
  // it would be interesting to first mark all the leaf nodes which are a subset of the search region , so all the particles within
  // such leaf nodes can be directly added. 
  while( !nodelist.empty() ) 
   {
/*printf("Nodelist size=%lu \n",nodelist.size());    */
     // Number at the front of the queue
    int index = nodelist[0];
     // Positions of the first and the last particles within the node pointed to by computational_particle->IndexOfFirstPhysicalParticle. 
    int first = m_vFirstParticleIndex[index];
    int end   = first + m_vNumberOfContainedParticles[index];
    
    // Cycle through all the particles within the current leaf node and look for neighbours.
    for (int i = first; i < end; i++) 
     {
      double dx = search_x - m_vCoordX[ m_vParticleKeyIndex[i].index ];
      double dy = search_y - m_vCoordY[ m_vParticleKeyIndex[i].index ];
      double dz = search_z - m_vCoordZ[ m_vParticleKeyIndex[i].index ];
      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( squared_distance < radius*radius) 
        {
          result->push_back(m_vParticleKeyIndex[i].index);
	  distance->push_back(sqrt(squared_distance));
/*printf("A neighbor found! Index=%d \n",m_vParticleKeyIndex[i].index);*/
        }
    }
    // Since all the information needed from the queue head has been processes we have no need for it any more. So pop it.
    nodelist.pop_front();
  }
    // Return the number of neighbours found.
    // The actual neighbours are pointed to by the result pointer. 
    return result->size();
}                

int Octree::searchNeighbor(const double &search_x,
                           const double &search_y,
                           const double &search_z,
                           const double &radius1,
			   const double &radius2,
			   std::vector<int> *result,
			   std::vector<double> *distance) 
{

   std::deque<int> nodelist;
   nodelist.push_back(0); // Push back the index of the root node, which is at position 0 in the octree vector.

   // In this while loop, we find all the leaf nodes of the octree which intersect the search region.  
   // If, at the end of the while loop, nodelist queue is empty, that means no leaf node intersects the search region and hence the search point has no particles in it neighbourhood.
   // In this we don't check whether an octree node is a subset of the search region. 
   // Paradoxically, this seems to be performing better than an algorithm which does. This needs some more investigating. 
   // Probably for deeper octrees and smaller search radii. that algorithm will endup outperforming this one. Also that means the subset check would need to be improved upon, since it could be the bottleneck. 
   while( 1 ) 
   {

     //std::cout << nodelist.size() << std::endl;
     // If the node-list is empty stop the execution of he current function, and return 0, since 0 neighbours have been founf. 
     // Also as a result, you see that nothing is in the deque pointed to be result.
    if( nodelist.empty() ) 
      { 
        return 0;
      }
    // Get the node number at the front of the queue.
    // The popping will be done at the end of the while loop.
    int inode       = nodelist[0];
    int first_index = m_vFirstChildIndex[inode];

    // If the head node of the deque is a leaf that means all the remaining nodes in the deque are also leaves. 
    // So we can now exit the while loop and process look through all the particles within the leaf nodes. 
    if( first_index == -1 ) 
      {
        break;
      }

    int last_index = first_index + m_vNumberOfChildren[ inode ];
    //  Find all the child nodes of the CURRENT node which intersect the sphere having radius "radius".
    //  After finding push_back their index onto the stack.
    for(int i=first_index ; i < last_index ; i++) 
     {
       if(   is_node_intersect_search_region( m_vLowerLimitOfX[i],\
                                              m_vUpperLimitOfX[i],\
                                              m_vLowerLimitOfY[i],\
                                              m_vUpperLimitOfY[i],\
                                              m_vLowerLimitOfZ[i],\
                                              m_vUpperLimitOfZ[i],\
                                              search_x,\
                                              search_y,\
                                              search_z, radius2) )
	 {
             nodelist.push_back( i );
	 }
    }
    nodelist.pop_front( );
  }
  // When the  above while-loop ends normally, the queue consists only of leaf nodes intersecting the search region.
  // We now look at every particle in each leaf node looking for neighbouring particles. Probably as an optimization step,
  // it would be interesting to first mark all the leaf nodes which are a subset of the search region , so all the particles within
  // such leaf nodes can be directly added. 
  while( !nodelist.empty() ) 
   {
    
     // Number at the front of the queue
    int index = nodelist[0];
     // Positions of the first and the last particles within the node pointed to by computational_particle->IndexOfFirstPhysicalParticle. 
    int first = m_vFirstParticleIndex[index];
    int end   = first + m_vNumberOfContainedParticles[index];
    
    // Cycle through all the particles within the current leaf node and look for neighbours.
    for (int i = first; i < end; i++) 
     {
      double dx = search_x - m_vCoordX[ m_vParticleKeyIndex[i].index ];
      double dy = search_y - m_vCoordY[ m_vParticleKeyIndex[i].index ];
      double dz = search_z - m_vCoordZ[ m_vParticleKeyIndex[i].index ];
      double squared_distance = dx*dx + dy*dy + dz*dz;
      if( squared_distance < radius2*radius2 && squared_distance > radius1*radius1) 
        {
          result->push_back(m_vParticleKeyIndex[i].index);
	  distance->push_back(sqrt(squared_distance));
        }
    }
    // Since all the information needed from the queue head has been processes we have no need for it any more. So pop it.
    nodelist.pop_front();
  }
    // Return the number of neighbours found.
    // The actual neighbours are pointed to by the result pointer. 
    return result->size();
}                


//int Octree::select_computational_particles_2to1_3d(std::vector<double> &cx, std::vector<double> &cy, std::vector<double> &cz, std::vector<double> &ch, std::vector<int> &ni, std::vector<int> &node_index, std::vector<int> &neigh24, double c, int order, int nmax, std::vector<int>& number_of_particles, int mindepth,int maxdepth,int bc) 
int Octree::select_computational_particles_2to1_3d(Computational_particle *computational_particle, double c, int order, int nmax, int mindepth,int maxdepth,int bc) 
/* Given an quadtree of physical particles, select some nodes to be temporary particle by the error criterion

	h^(2order-2)n>c,
where h is the diamater of corresponding cell, n is the number of physical particles in the cell, c and order are given parameter.

After the top->down sweep, some other temporary particles are added to enforce 2 to 1 property in a bottom->up sweep. That is, the level of two adjacent temporary nodes must be adjacent. If not, then the coarse node must be opened, and its children are temporary particles.

In the end, if a temporary particle has no child selected as temporary particle, then the particle is a computational particle (if it is selected in top->down sweep) or vacuum particle (if it is selected in bottom->up sweep). If non-periodic boundary is used, then some boundary particles are added in this step.

The error criterion is deduced from error balance of Monte Carlo noise and discritezation error. Theoritically this mesh size should give optimal mesh for density estimation.

The main motivation to enforce 2 to 1 property is to make sure the mesh size changes smoothly and avoid imbalanced GFD stencil. Moreover, 2 to 1 property gives well defined neighbours (Neigh24) which is a constant running time neighbour search data structure, much more effieicnt than octree search.

Input:
c 	Parameter in criterion, need to be optimized experimentally. Usually 10^(-6)<c<10^(-2).
order 	Parameter in criterion, order of GFD. order=2 is used in current code.

nmax	Maximum number of computational+vacuum+boundary particles. (size of given memory)
mindepth	Minimum level of node that are possible to be choosen as computational particle. 
maxdepth	Maximum level of node that are possible to be choosen as computational particle. 


Return value:
total	Total number of computational+vacuum+boundary particles, or -1 if there is not enough memory (nmax too small).

Output:

computational_particles->size Total number of computational+vacuum+boundary particles.
computational_particles->nc Number of computational particles.
computational_particles->nv Number of vacuum particles.
computational_particles->nb Number of boundary particles.

computational_particles->x	vector (size=total) of x coordinate of computational/vacuum/boundary particles.
computational_particles->y	vector (size=total) of y coordinate of computational/vacuum/boundary particles.
computational_particles->z	vector (size=total) of z coordinate of computational/vacuum/boundary particles.

computational_particles->CellLength	vector (size=total) of diamater of cell of computational/vacuum/boundary particles.

computational_particle->NumberOfPhysicalParticles	vector (size=total) of number of physical particles in cell of computational particles.
computational_particle->IndexOfFirstPhysicalParticle	vector (size=total) of index of first (sorted by Morton key) physical particles in cell of computational particles.
computational_particle->Type vector of type of computational(=0), vacuum(=-1), boundary(=1) particles.
computational_particle->FirstNeigh vector of location of first neighbour in Neigh24.
computational_particle->NumberNeigh	vector of number of neighbours.
computational_particle->Neigh24	vector (size<=24*size) of index of neighbours. For example, the first neighbour of ith particle is Neigh24[FirstNeigh[i]], and its last neighbour is Neigh24[FirstNeigh[i]+NumberNeigh[i]-1]. It is called Neigh24 since theoretically the maximum number of neighbour of each particle is 24, which can be easily shown from 2 to 1 property.

All memory for computational particles are allocated in this function, both these mentioned above and those not initialized in this function, i.e., density and potential. Previous allocated memory for computational particles is freed.

Xingyu Wang, Jul 20, 2015.

*/
{
//Data structure for temporary particles. It should be better to use std::vector with dynamic memory allocation to storage these data since we do not know the size of temporary particles in advanve, but I have not get the time to modify the code. The current walk around is to set an estimated size nmax in advance, call the function, and double nmax and call the function again if nmax is smaller than the needed size. 
//tx ty tz: center position
//leaf: type of node: if it is opened (not a leaf), it will not be selected as a future computational particle, and leaf[i]=0. If it is not opened, then it is a potential computational particle, leaf[i]=1 if it is given in top->down sweep (potential computational particle), and leaf[i]=-1 if it is given in bottom->up sweep (potential vacuum particle)
//n_p: Number of physical particles in the cell. We do not need to count the number here; they are already given in Octree::build.
//old_index: Index of temporaty node in the octree.
//particle_index: Index of first physical particle in the cell, sorted by Morton key.
//parent: the index of parent in tx/ty/tz. Set but not used in current code.
//level: level
//child[8*i+j]: the index of jth child in tx/ty/tz of ith node if ith node is opened (i.e., is not a leaf). 
//neigh[6*i+j]: the index of jth neighbour of ith node in tx/ty/tz. Since octree is a hierachical data structure and node from different level can all be selected as temporary node, so ith node can has more than one "neighbour" in jth direction in the usual sense. What we call neighbour is the one has the same level as ith node (if there exists one), or the one that has the lowest level.
double tx[nmax],ty[nmax],tz[nmax];
int leaf[nmax],n_p[nmax],old_index[nmax],particle_index[nmax],parent[nmax],level[nmax],child[8*nmax],neigh[6*nmax];
double tempx,tempy,tempz,h,h0,temp,dis,mindis,hx,hy,hz,hx0,hy0,hz0;
int i,j,k,l,inode,i1,i2,ml,np0,np1,np2,np3,np4,np5,np6,np7,i0,di,i3,i4,ln;
i=0;
j=1;
inode=0;
tx[0]=(m_vUpperLimitOfX[inode]+m_vLowerLimitOfX[inode])/2;
ty[0]=(m_vUpperLimitOfY[inode]+m_vLowerLimitOfY[inode])/2;
tz[0]=(m_vUpperLimitOfZ[inode]+m_vLowerLimitOfZ[inode])/2;
n_p[0]=m_vNumberOfContainedParticles[inode];
old_index[0]=0;
particle_index[0]=m_vFirstParticleIndex[inode];
parent[0]=-1;//no parent for root node
level[0]=m_vDepth[inode];
h0=(m_vUpperLimitOfX[inode]-m_vLowerLimitOfX[inode])*(m_vUpperLimitOfY[inode]-m_vLowerLimitOfY[inode])*(m_vUpperLimitOfZ[inode]-m_vLowerLimitOfZ[inode]);
h0=pow(h0,1.0/3);
hx0=(m_vUpperLimitOfX[inode]-m_vLowerLimitOfX[inode]);
hy0=(m_vUpperLimitOfY[inode]-m_vLowerLimitOfY[inode]);
hz0=(m_vUpperLimitOfZ[inode]-m_vLowerLimitOfZ[inode]);

//bc=0: non periodic boundary condition. -1 means no neighbour
neigh[0]=-1;//x+
neigh[1]=-1;//x-
neigh[2]=-1;//y- It should be y+ for consistency, but this order does not change the final result anyway
neigh[3]=-1;//y+
neigh[4]=-1;//z+
neigh[5]=-1;//z-

if(bc==1)//periodic boundary condition. All neighbours of root node are itself. This is the only modification we need for periodic domain (besides the periodic wrap function for distance calculation), since all other neighbours will be computed from here recursively.
{
neigh[0]=0;//x+
neigh[1]=0;//x-
neigh[2]=0;//y-
neigh[3]=0;//y+
neigh[4]=0;//z+
neigh[5]=0;//z-
//Need to set the children of root node here
child[0]=1;
child[1]=2;
child[2]=3;
child[3]=4;
child[4]=5;
child[5]=6;
child[6]=7;
child[7]=8;
}

//Construct from top to bottom
while((j-i)>0)
{
	inode=old_index[i];
	temp=1.0*n_p[i];
	h=h0;
	hx=hx0;
	hy=hy0;
	hz=hz0;
	for(i1=0;i1<level[i];i1++)
	{
		h=0.5*h;
		hx=0.5*hx;
		hy=0.5*hy;
		hz=0.5*hz;
	}
	for(i1=0;i1<(2*order-2);i1++)
		temp=temp*h;

/*	double sigma=0.1,a1=1.0;
	temp=a1*exp(-((tx[i]-m_vCoordX[0])*(tx[i]-m_vCoordX[0])+(ty[i]-m_vCoordY[0])*(ty[i]-m_vCoordY[0]))/sigma/sigma/2);
	for(i1=0;i1<(2*k);i1++)
		temp=temp*h;*/
	if(((temp>c)||(level[i]<mindepth))&&(level[i]<maxdepth))// need refinement: open ith cell
	{
		leaf[i]=0;
		l=0;
		if (inode>=0)
			l=m_vNumberOfChildren[inode];
		np0=0;
		np1=0;
		np2=0;
		np3=0;
		np4=0;
		np5=0;
		np6=0;
		np7=0;
		if (!l)// i is a leaf in octree or l is not a node of octree
		{
			for(i1=particle_index[i];i1<particle_index[i]+n_p[i];i1++)
			{
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
					np0=np0+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
					np1=np1+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
					np2=np2+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
					np3=np3+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
					np4=np4+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
					np5=np5+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
					np6=np6+1;
				if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
					np7=np7+1;
			}
		}

			//first child: ---
			i2=-1;
			k=j;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx<tx[i])&&(tempy<ty[i])&&(tempz<tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// first child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// first child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i];
				n_p[k]=np0;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i]=k;
			tx[k]=tx[i]-0.25*hx;
			ty[k]=ty[i]-0.25*hy;
			tz[k]=tz[i]-0.25*hz;
			neigh[6*k]=j+4;
			neigh[6*k+3]=j+2;
			neigh[6*k+4]=j+1;
			i1=neigh[6*i+1];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+4];
				neigh[6*k+1]=i1;
				neigh[6*i1]=k;
			}
			else
			{
				neigh[6*k+1]=i1;
			}
			i1=neigh[6*i+2];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+2];
				neigh[6*k+2]=i1;
				neigh[6*i1+3]=k;
			}
			else
			{
				neigh[6*k+2]=i1;
			}
			i1=neigh[6*i+5];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+1];
				neigh[6*k+5]=i1;
				neigh[6*i1+4]=k;
			}
			else
			{
				neigh[6*k+5]=i1;
			}

			//second child: --+
			i2=-1;
			k=j+1;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx<tx[i])&&(tempy<ty[i])&&(tempz>tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// second child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// second child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0;
				n_p[k]=np1;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+1]=k;
			tx[k]=tx[i]-0.25*hx;
			ty[k]=ty[i]-0.25*hy;
			tz[k]=tz[i]+0.25*hz;
			neigh[6*k]=j+5;
			neigh[6*k+3]=j+3;
			neigh[6*k+5]=j;
			i1=neigh[6*i+1];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+5];
				neigh[6*k+1]=i1;
				neigh[6*i1]=k;
			}
			else
			{
				neigh[6*k+1]=i1;
			}
			i1=neigh[6*i+2];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+3];
				neigh[6*k+2]=i1;
				neigh[6*i1+3]=k;
			}
			else
			{
				neigh[6*k+2]=i1;
			}
			i1=neigh[6*i+4];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1];
				neigh[6*k+4]=i1;
				neigh[6*i1+5]=k;
			}
			else
			{
				neigh[6*k+4]=i1;
			}

			//third child: -+-
			i2=-1;
			k=j+2;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx<tx[i])&&(tempy>ty[i])&&(tempz<tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// third child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// third child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0+np1;
				n_p[k]=np2;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+2]=k;
			tx[k]=tx[i]-0.25*hx;
			ty[k]=ty[i]+0.25*hy;
			tz[k]=tz[i]-0.25*hz;
			neigh[6*k]=j+6;
			neigh[6*k+2]=j;
			neigh[6*k+4]=j+3;
			i1=neigh[6*i+1];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+6];
				neigh[6*k+1]=i1;
				neigh[6*i1]=k;
			}
			else
			{
				neigh[6*k+1]=i1;
			}
			i1=neigh[6*i+3];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1];
				neigh[6*k+3]=i1;
				neigh[6*i1+2]=k;
			}
			else
			{
				neigh[6*k+3]=i1;
			}
			i1=neigh[6*i+5];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+3];
				neigh[6*k+5]=i1;
				neigh[6*i1+4]=k;
			}
			else
			{
				neigh[6*k+5]=i1;
			}

			//fourth child: -++
			i2=-1;
			k=j+3;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx<tx[i])&&(tempy>ty[i])&&(tempz>tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// fourth child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// fourth child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2;
				n_p[k]=np3;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+3]=k;
			tx[k]=tx[i]-0.25*hx;
			ty[k]=ty[i]+0.25*hy;
			tz[k]=tz[i]+0.25*hz;
			neigh[6*k]=j+7;
			neigh[6*k+2]=j+1;
			neigh[6*k+5]=j+2;
			i1=neigh[6*i+1];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+7];
				neigh[6*k+1]=i1;
				neigh[6*i1]=k;
			}
			else
			{
				neigh[6*k+1]=i1;
			}
			i1=neigh[6*i+3];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+1];
				neigh[6*k+3]=i1;
				neigh[6*i1+2]=k;
			}
			else
			{
				neigh[6*k+3]=i1;
			}
			i1=neigh[6*i+4];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+2];
				neigh[6*k+4]=i1;
				neigh[6*i1+5]=k;
			}
			else
			{
				neigh[6*k+4]=i1;
			}

			//fifth child: +--
			i2=-1;
			k=j+4;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx>tx[i])&&(tempy<ty[i])&&(tempz<tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// fifth child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// fifth child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3;
				n_p[k]=np4;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+4]=k;
			tx[k]=tx[i]+0.25*hx;
			ty[k]=ty[i]-0.25*hy;
			tz[k]=tz[i]-0.25*hz;
			neigh[6*k+1]=j+0;
			neigh[6*k+3]=j+6;
			neigh[6*k+4]=j+5;
			i1=neigh[6*i];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1];
				neigh[6*k]=i1;
				neigh[6*i1+1]=k;
			}
			else
			{
				neigh[6*k]=i1;
			}
			i1=neigh[6*i+2];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+6];
				neigh[6*k+2]=i1;
				neigh[6*i1+3]=k;
			}
			else
			{
				neigh[6*k+2]=i1;
			}
			i1=neigh[6*i+5];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+5];
				neigh[6*k+5]=i1;
				neigh[6*i1+4]=k;
			}
			else
			{
				neigh[6*k+5]=i1;
			}

			//sixth child: +-+
			i2=-1;
			k=j+5;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx>tx[i])&&(tempy<ty[i])&&(tempz>tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// sixth child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// sixth child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3+np4;
				n_p[k]=np5;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+5]=k;
			tx[k]=tx[i]+0.25*hx;
			ty[k]=ty[i]-0.25*hy;
			tz[k]=tz[i]+0.25*hz;
			neigh[6*k+1]=j+1;
			neigh[6*k+3]=j+7;
			neigh[6*k+5]=j+4;
			i1=neigh[6*i];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+1];
				neigh[6*k]=i1;
				neigh[6*i1+1]=k;
			}
			else
			{
				neigh[6*k]=i1;
			}
			i1=neigh[6*i+2];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+7];
				neigh[6*k+2]=i1;
				neigh[6*i1+3]=k;
			}
			else
			{
				neigh[6*k+2]=i1;
			}
			i1=neigh[6*i+4];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+4];
				neigh[6*k+4]=i1;
				neigh[6*i1+5]=k;
			}
			else
			{
				neigh[6*k+4]=i1;
			}

			//seventh child: ++-
			i2=-1;
			k=j+6;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx>tx[i])&&(tempy>ty[i])&&(tempz<tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// seventh child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// seventh child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3+np4+np5;
				n_p[k]=np6;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+6]=k;
			tx[k]=tx[i]+0.25*hx;
			ty[k]=ty[i]+0.25*hy;
			tz[k]=tz[i]-0.25*hz;
			neigh[6*k+1]=j+2;
			neigh[6*k+2]=j+4;
			neigh[6*k+4]=j+7;
			i1=neigh[6*i];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+2];
				neigh[6*k]=i1;
				neigh[6*i1+1]=k;
			}
			else
			{
				neigh[6*k]=i1;
			}
			i1=neigh[6*i+3];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+4];
				neigh[6*k+3]=i1;
				neigh[6*i1+2]=k;
			}
			else
			{
				neigh[6*k+3]=i1;
			}
			i1=neigh[6*i+5];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+7];
				neigh[6*k+5]=i1;
				neigh[6*i1+4]=k;
			}
			else
			{
				neigh[6*k+5]=i1;
			}

			//eighth child: +++
			i2=-1;
			k=j+7;
			for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
			{
				tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
				tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
				tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
				if ((tempx>tx[i])&&(tempy>ty[i])&&(tempz>tz[i]))
				{	
					i2=i1;
					break;
				}
			}
			if(i2>0)// eighth child exists in octree
			{
				old_index[k]=i2;
				particle_index[k]=m_vFirstParticleIndex[i2];
				n_p[k]=m_vNumberOfContainedParticles[i2];
			}
			else// eighth child does not exist in octree
			{
				old_index[k]=-1;
				particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3+np4+np5+np6;
				n_p[k]=np7;
			}			
			parent[k]=i;
			level[k]=level[i]+1;
			leaf[k]=-1;
			child[8*i+7]=k;
			tx[k]=tx[i]+0.25*hx;
			ty[k]=ty[i]+0.25*hy;
			tz[k]=tz[i]+0.25*hz;
			neigh[6*k+1]=j+3;
			neigh[6*k+2]=j+5;
			neigh[6*k+5]=j+6;
			i1=neigh[6*i];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+3];
				neigh[6*k]=i1;
				neigh[6*i1+1]=k;
			}
			else
			{
				neigh[6*k]=i1;
			}
			i1=neigh[6*i+3];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+5];
				neigh[6*k+3]=i1;
				neigh[6*i1+2]=k;
			}
			else
			{
				neigh[6*k+3]=i1;
			}
			i1=neigh[6*i+4];
			if((i1>=0)&&(leaf[i1]==0))
			{
				i1=child[8*i1+6];
				neigh[6*k+4]=i1;
				neigh[6*i1+5]=k;
			}
			else
			{
				neigh[6*k+4]=i1;
			}


			j=j+8;
			if((j+8)>nmax)
				return -1;
	}
	else
	{
		leaf[i]=1;
	}
	i=i+1;
}

//Refine from bottom to top
ml=level[i-1];
for(ln=ml;ln>1;ln--)
{
	i0=0;
	while(j>i0)
	{
		if(level[i0]==ln)
		{
			for(di=0;di<6;di++)
			{
				i=neigh[6*i0+di];
				while (((i>0)&&(leaf[i]==0))&&(level[i]<ln))
				{
//bc
					mindis=periodic_wrap(tx[i0]-tx[child[8*i]],0)*periodic_wrap(tx[i0]-tx[child[8*i]],0)+periodic_wrap(ty[i0]-ty[child[8*i]],1)*periodic_wrap(ty[i0]-ty[child[8*i]],1)+periodic_wrap(tz[i0]-tz[child[8*i]],2)*periodic_wrap(tz[i0]-tz[child[8*i]],2);						
					
					i3=0;
					for(i4=0;i4<8;i4++)
					{
						dis=periodic_wrap(tx[i0]-tx[child[8*i+i4]],0)*periodic_wrap(tx[i0]-tx[child[8*i+i4]],0)+periodic_wrap(ty[i0]-ty[child[8*i+i4]],1)*periodic_wrap(ty[i0]-ty[child[8*i+i4]],1)+periodic_wrap(tz[i0]-tz[child[8*i+i4]],2)*periodic_wrap(tz[i0]-tz[child[8*i+i4]],2);					
						if(dis<mindis)
						{
							mindis=dis;
							i3=i4;
						}
					}
					neigh[6*i0+di]=child[8*i+i3];
					i=child[8*i+i3];
				}

				while((i>0)&&(level[i]<(ln-1)))//2to1 not satisfied, refine i
				{


					inode=old_index[i];
					leaf[i]=0;
					l=0;
					if (inode>=0)
						l=m_vNumberOfChildren[inode];
					hx=hx0;
					hy=hy0;
					hz=hz0;
					for(i1=0;i1<level[i];i1++)
					{
						hx=0.5*hx;
						hy=0.5*hy;
						hz=0.5*hz;
					}
					np0=0;
					np1=0;
					np2=0;
					np3=0;
					np4=0;
					np5=0;
					np6=0;
					np7=0;
					if (!l)// i is a leaf in octree or l is not a node of octree
					{
						for(i1=particle_index[i];i1<particle_index[i]+n_p[i];i1++)
						{
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
								np0=np0+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
								np1=np1+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
								np2=np2+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]<tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
								np3=np3+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
								np4=np4+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]<ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
								np5=np5+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]<tz[i]))
								np6=np6+1;
							if ((m_vCoordX[m_vParticleKeyIndex[i1].index]>=tx[i])&&(m_vCoordY[m_vParticleKeyIndex[i1].index]>=ty[i])&&(m_vCoordZ[m_vParticleKeyIndex[i1].index]>=tz[i]))
								np7=np7+1;
						}
					}

						//first child: ---
						i2=-1;
						k=j;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx<tx[i])&&(tempy<ty[i])&&(tempz<tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// first child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// first child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i];
							n_p[k]=np0;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i]=k;
						tx[k]=tx[i]-0.25*hx;
						ty[k]=ty[i]-0.25*hy;
						tz[k]=tz[i]-0.25*hz;
						neigh[6*k]=j+4;
						neigh[6*k+3]=j+2;
						neigh[6*k+4]=j+1;

						i1=neigh[6*i+1];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+1]=i1;
						if(level[i1]==level[k])
							neigh[6*i1]=k;					

						i1=neigh[6*i+2];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+2]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+3]=k;


						i1=neigh[6*i+5];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+5]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+4]=k;


						//second child: --+
						i2=-1;
						k=j+1;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx<tx[i])&&(tempy<ty[i])&&(tempz>tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// second child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// second child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0;
							n_p[k]=np1;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+1]=k;
						tx[k]=tx[i]-0.25*hx;
						ty[k]=ty[i]-0.25*hy;
						tz[k]=tz[i]+0.25*hz;
						neigh[6*k]=j+5;
						neigh[6*k+3]=j+3;
						neigh[6*k+5]=j;

						i1=neigh[6*i+1];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+1]=i1;
						if(level[i1]==level[k])
							neigh[6*i1]=k;					

						i1=neigh[6*i+2];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+2]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+3]=k;


						i1=neigh[6*i+4];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+4]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+5]=k;


						//third child: -+-
						i2=-1;
						k=j+2;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx<tx[i])&&(tempy>ty[i])&&(tempz<tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// third child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// third child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0+np1;
							n_p[k]=np2;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+2]=k;
						tx[k]=tx[i]-0.25*hx;
						ty[k]=ty[i]+0.25*hy;
						tz[k]=tz[i]-0.25*hz;
						neigh[6*k]=j+6;
						neigh[6*k+2]=j;
						neigh[6*k+4]=j+3;

						i1=neigh[6*i+1];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+1]=i1;
						if(level[i1]==level[k])
							neigh[6*i1]=k;					

						i1=neigh[6*i+3];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+3]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+2]=k;


						i1=neigh[6*i+5];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+5]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+4]=k;


						//fourth child: -++
						i2=-1;
						k=j+3;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx<tx[i])&&(tempy>ty[i])&&(tempz>tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// fourth child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// fourth child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2;
							n_p[k]=np3;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+3]=k;
						tx[k]=tx[i]-0.25*hx;
						ty[k]=ty[i]+0.25*hy;
						tz[k]=tz[i]+0.25*hz;
						neigh[6*k]=j+7;
						neigh[6*k+2]=j+1;
						neigh[6*k+5]=j+2;

						i1=neigh[6*i+1];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+1]=i1;
						if(level[i1]==level[k])
							neigh[6*i1]=k;					

						i1=neigh[6*i+3];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+3]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+2]=k;


						i1=neigh[6*i+4];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+4]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+5]=k;


						//fifth child: +--
						i2=-1;
						k=j+4;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx>tx[i])&&(tempy<ty[i])&&(tempz<tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// fifth child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// fifth child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3;
							n_p[k]=np4;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+4]=k;
						tx[k]=tx[i]+0.25*hx;
						ty[k]=ty[i]-0.25*hy;
						tz[k]=tz[i]-0.25*hz;
						neigh[6*k+1]=j+0;
						neigh[6*k+3]=j+6;
						neigh[6*k+4]=j+5;

						i1=neigh[6*i];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+1]=k;					

						i1=neigh[6*i+2];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+2]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+3]=k;


						i1=neigh[6*i+5];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+5]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+4]=k;


						//sixth child: +-+
						i2=-1;
						k=j+5;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx>tx[i])&&(tempy<ty[i])&&(tempz>tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// sixth child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// sixth child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3+np4;
							n_p[k]=np5;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+5]=k;
						tx[k]=tx[i]+0.25*hx;
						ty[k]=ty[i]-0.25*hy;
						tz[k]=tz[i]+0.25*hz;
						neigh[6*k+1]=j+1;
						neigh[6*k+3]=j+7;
						neigh[6*k+5]=j+4;

						i1=neigh[6*i];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+1]=k;					

						i1=neigh[6*i+2];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+2]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+3]=k;


						i1=neigh[6*i+4];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+4]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+5]=k;

						//seventh child: ++-
						i2=-1;
						k=j+6;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx>tx[i])&&(tempy>ty[i])&&(tempz<tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// seventh child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// seventh child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3+np4+np5;
							n_p[k]=np6;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+6]=k;
						tx[k]=tx[i]+0.25*hx;
						ty[k]=ty[i]+0.25*hy;
						tz[k]=tz[i]-0.25*hz;
						neigh[6*k+1]=j+2;
						neigh[6*k+2]=j+4;
						neigh[6*k+4]=j+7;

						i1=neigh[6*i];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+1]=k;					

						i1=neigh[6*i+3];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+3]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+2]=k;


						i1=neigh[6*i+5];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+5]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+4]=k;
						

						//eighth child: +++
						i2=-1;
						k=j+7;
						for(i1=m_vFirstChildIndex[inode];i1<m_vFirstChildIndex[inode]+l;i1++)
						{
							tempx=(m_vUpperLimitOfX[i1]+m_vLowerLimitOfX[i1])/2;
							tempy=(m_vUpperLimitOfY[i1]+m_vLowerLimitOfY[i1])/2;
							tempz=(m_vUpperLimitOfZ[i1]+m_vLowerLimitOfZ[i1])/2;
							if ((tempx>tx[i])&&(tempy>ty[i])&&(tempz>tz[i]))
							{	
								i2=i1;
								break;
							}
						}
						if(i2>0)// eighth child exists in octree
						{
							old_index[k]=i2;
							particle_index[k]=m_vFirstParticleIndex[i2];
							n_p[k]=m_vNumberOfContainedParticles[i2];
						}
						else// eighth child does not exist in octree
						{
							old_index[k]=-1;
							particle_index[k]=m_vFirstParticleIndex[i]+np0+np1+np2+np3+np4+np5+np6;
							n_p[k]=np7;
						}			
						parent[k]=i;
						level[k]=level[i]+1;
						leaf[k]=-1;
						child[8*i+7]=k;
						tx[k]=tx[i]+0.25*hx;
						ty[k]=ty[i]+0.25*hy;
						tz[k]=tz[i]+0.25*hz;
						neigh[6*k+1]=j+3;
						neigh[6*k+2]=j+5;
						neigh[6*k+5]=j+6;

						i1=neigh[6*i];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+1]=k;					

						i1=neigh[6*i+3];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+3]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+2]=k;


						i1=neigh[6*i+4];
						while (((i1>0)&&(leaf[i1]==0))&&(level[i1]<level[k]))
						{
							mindis=periodic_wrap(tx[k]-tx[child[8*i1]],0)*periodic_wrap(tx[k]-tx[child[8*i1]],0)+periodic_wrap(ty[k]-ty[child[8*i1]],1)*periodic_wrap(ty[k]-ty[child[8*i1]],1)+periodic_wrap(tz[k]-tz[child[8*i1]],2)*periodic_wrap(tz[k]-tz[child[8*i1]],2);						
							i3=0;
							for(i4=0;i4<8;i4++)
							{
								dis=periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)*periodic_wrap(tx[k]-tx[child[8*i1+i4]],0)+periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)*periodic_wrap(ty[k]-ty[child[8*i1+i4]],1)+periodic_wrap(tz[k]-tz[child[8*i1+i4]],2)*periodic_wrap(tz[k]-tz[child[8*i1+i4]],2);					
								if(dis<mindis)
								{
									mindis=dis;
									i3=i4;
								}
							}
							i1=child[8*i1+i3];
						}
						neigh[6*k+4]=i1;
						if(level[i1]==level[k])
							neigh[6*i1+5]=k;


						j=j+8;	
						if((j+8)>nmax)
							return -1;


					i=neigh[6*i0+di];
					while (((i>0)&&(leaf[i]==0))&&(level[i]<ln))
					{
						mindis=periodic_wrap(tx[i0]-tx[child[8*i]],0)*periodic_wrap(tx[i0]-tx[child[8*i]],0)+periodic_wrap(ty[i0]-ty[child[8*i]],1)*periodic_wrap(ty[i0]-ty[child[8*i]],1)+periodic_wrap(tz[i0]-tz[child[8*i]],2)*periodic_wrap(tz[i0]-tz[child[8*i]],2);
						i3=0;
						for(i4=0;i4<8;i4++)
						{
							dis=periodic_wrap(tx[i0]-tx[child[8*i+i4]],0)*periodic_wrap(tx[i0]-tx[child[8*i+i4]],0)+periodic_wrap(ty[i0]-ty[child[8*i+i4]],1)*periodic_wrap(ty[i0]-ty[child[8*i+i4]],1)+periodic_wrap(tz[i0]-tz[child[8*i+i4]],2)*periodic_wrap(tz[i0]-tz[child[8*i+i4]],2);
							if(dis<mindis)
							{
								mindis=dis;
								i3=i4;
							}
						}
						neigh[6*i0+di]=child[8*i+i3];
						i=child[8*i+i3];
					}					
				}
			}
		}
	i0=i0+1;
	}
	i0=0;
	while(j>i0)
	{
		if(level[i0]==ln)
		{
			for(di=0;di<6;di++)
			{
				i=neigh[6*i0+di];
				while (((i>0)&&(leaf[i]==0))&&(level[i]<ln))
				{
					mindis=periodic_wrap(tx[i0]-tx[child[8*i]],0)*periodic_wrap(tx[i0]-tx[child[8*i]],0)+periodic_wrap(ty[i0]-ty[child[8*i]],1)*periodic_wrap(ty[i0]-ty[child[8*i]],1)+periodic_wrap(tz[i0]-tz[child[8*i]],2)*periodic_wrap(tz[i0]-tz[child[8*i]],2);
					i3=0;
					for(i4=0;i4<8;i4++)
					{
						dis=periodic_wrap(tx[i0]-tx[child[8*i+i4]],0)*periodic_wrap(tx[i0]-tx[child[8*i+i4]],0)+periodic_wrap(ty[i0]-ty[child[8*i+i4]],1)*periodic_wrap(ty[i0]-ty[child[8*i+i4]],1)+periodic_wrap(tz[i0]-tz[child[8*i+i4]],2)*periodic_wrap(tz[i0]-tz[child[8*i+i4]],2);		
						if(dis<mindis)
						{
							mindis=dis;
							i3=i4;
						}
					}
					neigh[6*i0+di]=child[8*i+i3];
					i=child[8*i+i3];
				}
			}
		}
	i0=i0+1;
	}			
}

//Now with all temporary particles found, we can select computational particles now.

//free the previous allocated memory for computational particles
clear(computational_particle);

int nc, nv, nb;
nc=0;
nv=0;
nb=0;
double temp_length;
for(i=0;i<j;i++)//particles given by first sweep
	if(leaf[i]==1)
	{
		computational_particle->x.push_back (tx[i]);
		computational_particle->y.push_back (ty[i]);
		computational_particle->z.push_back (tz[i]);
		temp_length=h0;
		for(i1=0;i1<level[i];i1++)
			temp_length=0.5*temp_length;
		computational_particle->CellLength.push_back (temp_length);
		computational_particle->NumberOfPhysicalParticles.push_back (n_p[i]);
		computational_particle->IndexOfFirstPhysicalParticle.push_back (particle_index[i]);
		computational_particle->Type.push_back (0);
		particle_index[i]=nc;//particle_index becomes the new index in output
		nc=nc+1;
	}
for(i=0;i<j;i++)//particles given by second sweep
	if(leaf[i]==-1)
	{
		computational_particle->x.push_back (tx[i]);
		computational_particle->y.push_back (ty[i]);
		computational_particle->z.push_back (tz[i]);
		temp_length=h0;
		for(i1=0;i1<level[i];i1++)
			temp_length=0.5*temp_length;
		computational_particle->CellLength.push_back (temp_length);
		computational_particle->NumberOfPhysicalParticles.push_back (n_p[i]);
		computational_particle->IndexOfFirstPhysicalParticle.push_back (particle_index[i]);
		computational_particle->Type.push_back (-1);
		particle_index[i]=nc+nv;//particle_index becomes the new index in output
		nv=nv+1;
	}
i0=0;
int total_number_of_neighbors=0;
for(i=0;i<j;i++)//boundary particles and neigh24
	if(leaf[i])
	{
		computational_particle->FirstNeigh.push_back(total_number_of_neighbors);
		i2=0;

		if(neigh[6*i]==-1)//xmax boundary
		{
			computational_particle->x.push_back (1);
			computational_particle->y.push_back (ty[i]);
			computational_particle->z.push_back (tz[i]);
			temp_length=h0;
			for(i1=0;i1<level[i];i1++)
				temp_length=0.5*temp_length;
			computational_particle->CellLength.push_back (temp_length);
			computational_particle->NumberOfPhysicalParticles.push_back (0);
			computational_particle->IndexOfFirstPhysicalParticle.push_back (-1);
			computational_particle->Neigh24.push_back(nc+nv+nb);
			i2=i2+1;
			nb=nb+1;
		}
		else
		{
			i3=neigh[6*i];
			if(leaf[i3])//i3 is a node in output
			{
				computational_particle->Neigh24.push_back(particle_index[i3]);		
				i2=i2+1;
			}
			else//i3 is not a node in output, its children are in output because of 2to1 balance
			{
				computational_particle->Neigh24.push_back(particle_index[child[8*i3]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+1]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+2]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+3]]);
				i2=i2+4;
			}
		}

		if(neigh[6*i+1]==-1)//xmin boundary
		{
			computational_particle->x.push_back (-1);
			computational_particle->y.push_back (ty[i]);
			computational_particle->z.push_back (tz[i]);
			temp_length=h0;
			for(i1=0;i1<level[i];i1++)
				temp_length=0.5*temp_length;
			computational_particle->CellLength.push_back (temp_length);
			computational_particle->NumberOfPhysicalParticles.push_back (0);
			computational_particle->IndexOfFirstPhysicalParticle.push_back (-1);
			computational_particle->Neigh24.push_back(nc+nv+nb);
			i2=i2+1;
			nb=nb+1;
		}
		else
		{
			i3=neigh[6*i+1];
			if(leaf[i3])//i3 is a node in output
			{
				computational_particle->Neigh24.push_back(particle_index[i3]);	
				i2=i2+1;
			}
			else//i3 is not a node in output, its children are in output because of 2to1 balance
			{
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+4]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+5]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+6]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+7]]);
				i2=i2+4;
			}
		}

		if(neigh[6*i+3]==-1)//ymax boundary
		{

			computational_particle->x.push_back (tx[i]);
			computational_particle->y.push_back (1);
			computational_particle->z.push_back (tz[i]);
			temp_length=h0;
			for(i1=0;i1<level[i];i1++)
				temp_length=0.5*temp_length;
			computational_particle->CellLength.push_back (temp_length);
			computational_particle->NumberOfPhysicalParticles.push_back (0);
			computational_particle->IndexOfFirstPhysicalParticle.push_back (-1);
			computational_particle->Neigh24.push_back(nc+nv+nb);
			i2=i2+1;
			nb=nb+1;
		}
		else
		{
			i3=neigh[6*i+3];
			if(leaf[i3])//i3 is a node in output
			{
				computational_particle->Neigh24.push_back(particle_index[i3]);	
				i2=i2+1;
			}
			else//i3 is not a node in output, its children are in output because of 2to1 balance
			{
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+0]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+1]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+4]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+5]]);
				i2=i2+4;
			}
		}

		if(neigh[6*i+2]==-1)//ymin boundary
		{

			computational_particle->x.push_back (tx[i]);
			computational_particle->y.push_back (-1);
			computational_particle->z.push_back (tz[i]);
			temp_length=h0;
			for(i1=0;i1<level[i];i1++)
				temp_length=0.5*temp_length;
			computational_particle->CellLength.push_back (temp_length);
			computational_particle->NumberOfPhysicalParticles.push_back (0);
			computational_particle->IndexOfFirstPhysicalParticle.push_back (-1);
			computational_particle->Neigh24.push_back(nc+nv+nb);
			i2=i2+1;
			nb=nb+1;
		}
		else
		{
			i3=neigh[6*i+2];
			if(leaf[i3])//i3 is a node in output
			{
				computational_particle->Neigh24.push_back(particle_index[i3]);
				i2=i2+1;
			}
			else//i3 is not a node in output, its children are in output because of 2to1 balance
			{
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+2]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+3]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+6]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+7]]);
				i2=i2+4;
			}
		}

		if(neigh[6*i+4]==-1)//zmax boundary
		{

			computational_particle->x.push_back (tx[i]);
			computational_particle->y.push_back (ty[i]);
			computational_particle->z.push_back (1);
			temp_length=h0;
			for(i1=0;i1<level[i];i1++)
				temp_length=0.5*temp_length;
			computational_particle->CellLength.push_back (temp_length);
			computational_particle->NumberOfPhysicalParticles.push_back (0);
			computational_particle->IndexOfFirstPhysicalParticle.push_back (-1);
			computational_particle->Neigh24.push_back(nc+nv+nb);
			i2=i2+1;
			nb=nb+1;
		}
		else
		{
			i3=neigh[6*i+4];
			if(leaf[i3])//i3 is a node in output
			{
				computational_particle->Neigh24.push_back(particle_index[i3]);
				i2=i2+1;
			}
			else//i3 is not a node in output, its children are in output because of 2to1 balance
			{
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+0]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+2]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+4]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+6]]);
				i2=i2+4;
			}
		}

		if(neigh[6*i+5]==-1)//zmin boundary
		{

			computational_particle->x.push_back (tx[i]);
			computational_particle->y.push_back (ty[i]);
			computational_particle->z.push_back (-1);
			temp_length=h0;
			for(i1=0;i1<level[i];i1++)
				temp_length=0.5*temp_length;
			computational_particle->CellLength.push_back (temp_length);
			computational_particle->NumberOfPhysicalParticles.push_back (0);
			computational_particle->IndexOfFirstPhysicalParticle.push_back (-1);
			computational_particle->Neigh24.push_back(nc+nv+nb);
			i2=i2+1;
			nb=nb+1;
		}
		else
		{
			i3=neigh[6*i+5];
			if(leaf[i3])//i3 is a node in output
			{
				computational_particle->Neigh24.push_back(particle_index[i3]);
				i2=i2+1;
			}
			else//i3 is not a node in output, its children are in output because of 2to1 balance
			{
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+1]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+3]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+5]]);
				computational_particle->Neigh24.push_back(particle_index[child[8*i3+7]]);
				i2=i2+4;
			}
		}
		total_number_of_neighbors=total_number_of_neighbors+i2;
		computational_particle->NumberNeigh.push_back(i2);
		i0=i0+1;
	}
for(i=nc+nv;i<(nc+nv+nb);i++)
{
	computational_particle->NumberNeigh.push_back(0);
	computational_particle->FirstNeigh.push_back(total_number_of_neighbors);
	computational_particle->Type.push_back (1);
}
//Initialize other variable in computational particle.
computational_particle->nc=nc;
computational_particle->nv=nv;
computational_particle->nb=nb;
computational_particle->size=nc+nv+nb;

computational_particle->Density.resize(computational_particle->size);
computational_particle->Potential.resize(computational_particle->size);
computational_particle->Density.resize(computational_particle->size);
computational_particle->p_x.resize(computational_particle->size);
computational_particle->p_y.resize(computational_particle->size);
computational_particle->p_z.resize(computational_particle->size);
computational_particle->p_xx.resize(computational_particle->size);
computational_particle->p_yy.resize(computational_particle->size);
computational_particle->p_zz.resize(computational_particle->size);
computational_particle->p_xy.resize(computational_particle->size);
computational_particle->p_xz.resize(computational_particle->size);
computational_particle->p_yz.resize(computational_particle->size);

return (nc+nv+nb);
}

int Octree::interpolation_from_node_to_particle_2to1_3d(Computational_particle *cp, Physical_particle *pp)
{
	int i;
	#pragma omp parallel
	{
		int j,k,l;
		double dx,dy,dz;
	 	#pragma omp for
		for(i=0;i<(cp->nc+cp->nv);i++)
		{
			j=cp->IndexOfFirstPhysicalParticle[i];
			for(k=j;k<j+cp->NumberOfPhysicalParticles[i];k++)
			{
				l=m_vParticleKeyIndex[k].index;
				dx=m_vCoordX[l]-cp->x[i];
				dy=m_vCoordY[l]-cp->y[i];
				dz=m_vCoordZ[l]-cp->z[i];
				pp->a_x[l]=cp->p_x[i]+dx*cp->p_xx[i]+dy*cp->p_xy[i]+dz*cp->p_xz[i];
				pp->a_y[l]=cp->p_y[i]+dx*cp->p_xy[i]+dy*cp->p_yy[i]+dz*cp->p_yz[i];
				pp->a_z[l]=cp->p_z[i]+dx*cp->p_xz[i]+dy*cp->p_yz[i]+dz*cp->p_zz[i];

//				pp->a_x[l]=cp->p_x[i];
//				pp->a_y[l]=cp->p_y[i];
//				pp->a_z[l]=cp->p_z[i];
			}
		}
	}
	return 0;
}


double Octree::periodic_wrap(double x, int dimension)
{
double length;
  if(bc==1)
  {
		if(dimension==0)//x
			length=m_dBoundingBox_max_x-m_dBoundingBox_min_x;
		if(dimension==1)//y
			length=m_dBoundingBox_max_y-m_dBoundingBox_min_y;
		if(dimension==2)//z
			length=m_dBoundingBox_max_z-m_dBoundingBox_min_z;

   	while(x >= (0.5*length))
     	x -= length;

   	while(x < -(0.5*length))
     	x += length;
  }
  return x;
}










