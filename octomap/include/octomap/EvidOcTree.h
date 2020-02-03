#ifndef OCTOMAP_EVID_OCTREE_H
#define OCTOMAP_EVID_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <iostream>

namespace octomap {

	// forward declaration for "friend"
	class EvidOcTree;

	// node definition
	class EvidOcTreeNode : public OcTreeNode {
	public:
		friend class EvidOcTree;  
		
		struct EvidState {
			EvidState() : f(0.), o(0.), ign(1.0), conf(0.) {};
			EvidState(float _f, float _o, float _ign, float _conf) 
				: f(_f), o(_o), ign(_ign), conf(_conf) {};

			float f, o, ign, conf;  // power set of states 
		};

		EvidOcTreeNode() : OcTreeNode() {}

		EvidOcTreeNode(const EvidOcTreeNode& rhs) : OcTreeNode(rhs), evid(rhs.evid) {}	


	protected:
		EvidState evid;

	}





} // end namespace




#endif