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
			EvidMass() : conf(0.), free(0.), occu(0.), igno(0.99) {};
			EvidMass(float _c, float _f, float _o, float _i) 
				: conf(_c), free(_f), occu(_o), igno(_i) {};

			void operator= (const EvidMass& rhs){
				conf = rhs.conf;
				free = rhs.free; 
				occu = rhs.occu; 
				igno = rhs.igno;}

			bool operator== (const EvidMass& rhs){
				return (conf==rhs.conf) && (free == rhs.free) && (occu == rhs.occu) && (igno == rhs.igno);
			}

			float conf_() const {return conf;}
			float free_() const {return free;}
			float occu_() const {return occu;}
			float igno_() const {return igno;}

		protected:
			float conf, free, occu, igno;  // power set of states 
		};

		EvidOcTreeNode() : OcTreeNode(), timestamp(0) {}
		EvidOcTreeNode(const EvidOcTreeNode& rhs) : OcTreeNode(rhs), mass(rhs.mass), timestamp(rhs.timestamp) {}	
		
		inline EvidMass getMass() const {return mass;}
		inline void setMass(float _c, float _f, float _o, float _i) {mass = EvidMass(_c, _f, _o, _i);}
		inline void setMass(const EvidMass& _mass) {mass = _mass;}
		inline void decayMass(const float alp){
			EvidMass _mass = EvidMass(alp*mass.conf_(),
															alp*mass.free_(),
															alp*mass.occu_(),
															alp*mass.igno_() + 1.0f-alp);
			setMass(_mass);
		}

		void copyData(const EvidOcTreeNode& from){
			OcTreeNode::copyData(from);
			mass = from.getMass();
			setTimestamp(from.getTimestamp());
		}

		/// update this node's evidential mass according to its children's evidential mass
		void updateMassChildren();

		/// update this node's occupancy probability according to its evidential mass computed from its children's evidential mass
		void updateOccupancyChildrenStamped(uint32_t timestamp);

		/// convert evidential mass to occupancy probability
		float evidMassToLogOdds() const;

		/// Average mass of a node's children
		EvidOcTreeNode::EvidMass getAverageChildMass() const;

		/**
		 * Helper function for handling timestamp.
		 * Only provide function for setting timestamp equal to an input so that
		 * timestamp can be set according to ROS timestamp
		 */
		inline uint32_t getTimestamp() const { return timestamp; }
    inline void setTimestamp(uint32_t t) {timestamp = t; }

	protected:
		EvidMass mass;
		uint32_t timestamp;
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

    EvidOcTreeNode* updateNode(const point3d& value, bool occupied, uint32_t timestamp, bool lazy_eval = false);

		EvidOcTreeNode* updateNode(const OcTreeKey& key, bool occupied, uint32_t timestamp, bool lazy_eval = false);

		EvidOcTreeNode* updateNode(const OcTreeKey& key, BasicBeliefAssignment bba_m2, uint32_t timestamp, bool lazy_eval = false);
		
		EvidOcTreeNode* updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                           unsigned int depth, const BasicBeliefAssignment& bba_m2, uint32_t timestamp, bool lazy_eval = false);
		
		/**
		 * Counter part of updateNodeLogOdds in Evidential Grid
		 * This function perform 2-step evidential fusion (conjunctive & dempster normalization)
		 * Node's log odds is preserved for pruning & occupancy queries but no longer directly
		 * updated. Instead, node's occupancy probability is retrieved by the Pignistic Transformation
		 */
		void upadteNodeEvidMass(const OcTreeKey& key, EvidOcTreeNode* evidNode, const BasicBeliefAssignment& bba_m2, uint32_t timestamp);

		/**
		 * Function for updating occupancy of inner nodes after integrating measurement with lazy_eval on
		 * This function needs to be called after updating node with lazy_evale enabled to
		 * ensure multi resolution behavior.
		 */
		void updateInnerOccupancy(uint32_t timestamp);

		/**
		 * Function for publishing coordinates of cells having high conflict mass.
		 * To get coordinates of these cells, their keys are first converted "point3d"
		 * After publishing this function wipe out vector "conf_keys" for updating with incoming measurement 
		 */ 
		void publishMovingCells();

	protected:
		void updateInnerOccupancyRecurs(EvidOcTreeNode* node, unsigned int depth, uint32_t timestamp);

		// initial evidential mass corresponding to whether a cell is occupied or free
		const float bba_mo = 0.7f, bba_mf = 0.7f; // 0.8
		
		// time constance
		const float tau = 1.3f;

		// threshold of conflict mass for detecting moving objects
		const float thres_conf = 0.005f;  //TODO: tune this

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

	std::ostream& operator<<(std::ostream& out, const EvidOcTreeNode::EvidMass& mass) {
    return out << '(' << mass.conf_() << ' ' << mass.free_() << ' ' << mass.occu_() << ' ' << mass.igno_() << ')';
  }

} // end namespace

#endif