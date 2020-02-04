#ifndef OCTOMAP_EVID_OCTREE_H
#define OCTOMAP_EVID_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <iostream>

template<typename T>
bool isNull(T x){return fabs(x) < 1e-5;}

namespace octomap {

	// forward declaration for "friend"
	class EvidOcTree;

	// node definition
	class EvidOcTreeNode : public OcTreeNode {
	public:
		friend class EvidOcTree;  
		
		class EvidMass {
		public:	
			EvidMass() : conf(0.), free(0.), occu(0.), igno(1.0) {};
			EvidMass(float _c, float _f, float _o, float _i) 
				: conf(_c), free(_f), occu(_o), igno(_i) {};

			void operator= (const EvidMass& rhs){
				conf = rhs.conf;
				free = rhs.free; 
				occu = rhs.occu; 
				igno = rhs.igno;}

			float conf_() const {return conf;}
			float free_() const {return free;}
			float occu_() const {return occu;}
			float igno_() const {return igno;}
		
		protected:
			float free, occu, igno, conf;  // power set of states 
		};

		EvidOcTreeNode() : OcTreeNode() {}
		EvidOcTreeNode(const EvidOcTreeNode& rhs) : OcTreeNode(rhs), mass(rhs.mass) {}	
		
		inline EvidMass getMass() const {return mass;}
		inline void setMass(float _c, float _f, float _o, float _i) {mass = EvidMass(_c, _f, _o, _i);}
		inline void setMass(const EvidMass& _mass) {mass = _mass;}

		inline bool isMassSet() const {return isNull(mass.igno_()-1.0);}

		void copyData(const EvidOcTreeNode& from){
			OcTreeNode::copyData(from);
			this->mass = from.getMass();
		}

		/// update this node's evidential mass according to its children's evidential mass
		void updateMassChildren();

		EvidOcTreeNode::EvidMass getAverageChildMass() const;

	protected:
		EvidMass mass;
		bool occ_before;  // whether occupied in previous pointcloud
	};

} // end namespace

#endif