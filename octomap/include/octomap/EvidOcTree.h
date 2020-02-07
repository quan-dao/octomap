#ifndef OCTOMAP_EVID_OCTREE_H
#define OCTOMAP_EVID_OCTREE_H

#include <octomap/OcTreeNode.h>
#include <octomap/OccupancyOcTreeBase.h>
#include <iostream>
#include <vector>

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
			float conf, free, occu, igno;  // power set of states 
		};

		EvidOcTreeNode() : OcTreeNode() {}
		EvidOcTreeNode(const EvidOcTreeNode& rhs) : OcTreeNode(rhs), mass(rhs.mass) {}	
		
		inline EvidMass getMass() const {return mass;}
		inline void setMass(float _c, float _f, float _o, float _i) {mass = EvidMass(_c, _f, _o, _i);}
		inline void setMass(const EvidMass& _mass) {mass = _mass;}

		void copyData(const EvidOcTreeNode& from){
			OcTreeNode::copyData(from);
			this->mass = from.getMass();
		}

    inline void printMass() {
      std::cout << mass.conf_() << "\t" << mass.free_() << "\t" << mass.occu_() << "\t" << mass.igno_() << "\n";
    }

		/// update this node's evidential mass according to its children's evidential mass
		void updateMassChildren();

		/// update this node's occupancy probability according to its evidential mass computed from its children's evidential mass
		void updateOccupancyChildren();

		/// convert evidential mass to occupancy probability
		float evidMassToLogOdds() const;

		EvidOcTreeNode::EvidMass getAverageChildMass() const;

	protected:
		EvidMass mass;
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
			BasicBeliefAssignment() : mc(0.0), mf(0.0), mo(0.0), mi(0.0) {}
			
			float mc, mf, mo, mi;  // evidential mass
		};

    EvidOcTreeNode* updateNode(const point3d& value, bool occupied, bool lazy_eval = false);

		EvidOcTreeNode* updateNode(const OcTreeKey& key, bool occupied, bool lazy_eval = false);

		EvidOcTreeNode* updateNode(const OcTreeKey& key, BasicBeliefAssignment bba_m2, bool lazy_eval = false);
		
		EvidOcTreeNode* updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                           unsigned int depth, const BasicBeliefAssignment& bba_m2, bool lazy_eval = false);
		
		/**
		 * Counter part of updateNodeLogOdds in Evidential Grid
		 * This function perform 2-step evidential fusion (conjunctive & dempster normalization)
		 * Node's log odds is preserved for pruning & occupancy queries but no longer directly
		 * updated. Instead, node's occupancy probability is retrieved by the Pignistic Transformation
		 */
		void upadteNodeEvidMass(const OcTreeKey& key, EvidOcTreeNode* evidNode, const BasicBeliefAssignment& bba_m2);

		/**
		 * Function for updating occupancy of inner nodes after integrating measurement with lazy_eval on
		 * This function needs to be called after updating node with lazy_evale enabled to
		 * ensure multi resolution behavior.
		 */
		void updateInnerOccupancy();

	protected:
		void updateInnerOccupancyRecurs(EvidOcTreeNode* node, unsigned int depth);

		// initial evidential mass corresponding to whether a cell is occupied or free
		const float bba_mo = 0.7f, bba_mf = 0.85f;

		// threshold of conflict mass for detecting moving objects
		const float thres_conf = 0.35f;  //TODO: tune this

    std::vector<OcTreeKey> conf_keys;


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