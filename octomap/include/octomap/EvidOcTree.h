#ifndef OCTOMAP_EVID_OCTREE_H
#define OCTOMAP_EVID_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <iostream>

namespace octomap {

	// forward declaration for "friend"
	class EvidOcTree;

	// node definition ---------------------------------------------
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

		EvidOcTreeNode() : OcTreeNode(), igno_before(true) {}
		EvidOcTreeNode(const EvidOcTreeNode& rhs) : 
			OcTreeNode(rhs), mass(rhs.mass), occ_before(rhs.occ_before), igno_before(rhs.igno_before) {}	
		
		inline EvidMass getMass() const {return mass;}
		inline void setMass(float _c, float _f, float _o, float _i) {mass = EvidMass(_c, _f, _o, _i);}
		inline void setMass(const EvidMass& _mass) {mass = _mass;}

		inline bool isMassSet() const {return igno_before;}

		void copyData(const EvidOcTreeNode& from){
			OcTreeNode::copyData(from);
			this->mass = from.getMass();
			occ_before = from.occ_before;
			igno_before = from.igno_before;
		}

		/// update this node's evidential mass according to its children's evidential mass
		void updateMassChildren();

		EvidOcTreeNode::EvidMass getAverageChildMass() const;

	protected:
		EvidMass mass;
		bool occ_before;  // whether occupied in previous pointcloud. Only have meaning if igno_before = false
		bool igno_before; // whether is uknown in prev pointlcoud
	};

	// tree definition ---------------------------------------------
	class EvidOcTree : public OccupancyOcTreeBase <EvidOcTreeNode> {
	
	public:
		/// Default constructor, set resolutions of leafs
		EvidOcTree(double resolution);

		/// virtual constructor: create a new object of the same type
		EvidOcTree* create() const {return new EvidOcTree(resolution);}

		std::string getTreeType() const {return "EvidOcTree";}

		/**
		 * Prune a node when it is collapsible. This overloaded version only 
		 * considers the node occupancy for pruning while ingoring the evid mass
		 */
		bool pruneNode(EvidOcTreeNode* node);

		bool isNodeCollapsible(const EvidOcTreeNode* node) const;
		
		/**
		 * Basic Belief Assignment for assigning evidential mass to a cell
		 * given LiDAR measurement. 
		 */ 
		struct BasicBeliefAssignment{
			BasicBeliefAssignment(float _c, float _f, float _o, float _i, bool _isOcc) 
					: mc(_c), mf(_f), mo(_o), mi(_i), isOcc(_isOcc) {}
			BasicBeliefAssignment(const BasicBeliefAssignment& rhs) : mc(rhs.mc), mf(rhs.mf), mo(rhs.mo), mi(rhs.mi), isOcc(rhs.isOcc) {}

			float mc, mf, mo, mi;  // evidential mass
			bool isOcc;
		};

		EvidOcTreeNode* updateNode(const OcTreeKey& key, BasicBeliefAssignment bba_m2, bool lazy_eval = false);
		
		EvidOcTreeNode* updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                           unsigned int depth, const BasicBeliefAssignment& bba_m2, bool lazy_eval = false);
		
		/**
		 * Counter part of updateNodeLogOdds in Evidential Grid
		 * This function perform 2-step evidential fusion (conjunctive & dempster normalization)
		 * Node's log odds is preserved for pruning & occupancy queries but no longer directly
		 * updated. Instead, node's occupancy probability is retrieved by the Pignistic Transformation
		 */
		void upadteNodeEvidMass(EvidOcTreeNode* evidNode, const BasicBeliefAssignment& bba_m2) const;

	protected:
		/**
     * Static member object which ensures that this OcTree's prototype
     * ends up in the classIDMapping only once. You need this as a 
     * static member in any derived octree class in order to read .ot
     * files through the AbstractOcTree factory. You should also call
     * ensureLinking() once from the constructor.
     */
    class StaticMemberInitializer{
       public:
         StaticMemberInitializer() {
           EvidOcTree* tree = new EvidOcTree(0.1);
           tree->clearKeyRays();
           AbstractOcTree::registerTreeType(tree);
         }

         /**
         * Dummy function to ensure that MSVC does not drop the
         * StaticMemberInitializer, causing this tree failing to register.
         * Needs to be called from the constructor of this octree.
         */
         void ensureLinking() {};
    };
    /// static member to ensure static initialization (only once)
    static StaticMemberInitializer evidOcTreeMemberInit;

	};

} // end namespace

#endif