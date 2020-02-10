#include <octomap/EvidOcTree.h>
#include <octomap/octomap_utils.h>
#include <cmath>


namespace octomap {

	// node implementation ------------------------------------------------------
	EvidOcTreeNode::EvidMass EvidOcTreeNode::getAverageChildMass() const {
		float _mc = 0.0, _mf = 0.0, _mo = 0.0, _mi = 0.0;
		int num_child = 0;

		if (children != NULL){
			for (int i=0; i<8; i++){
				EvidOcTreeNode* child = static_cast<EvidOcTreeNode*>(children[i]);
				if (child != NULL){
					_mc += child->mass.conf_();
					_mf += child->mass.free_();
					_mo += child->mass.occu_();
					_mi += child->mass.igno_();
					num_child ++;
				}
			}
		}

		if (num_child > 0) {
			_mc /= num_child;
			_mf /= num_child;
			_mo /= num_child;
			_mi /= num_child;
			return EvidMass(_mc, _mf, _mo, _mi);
		} else {
			return EvidMass();
		}
	}

	void EvidOcTreeNode::updateMassChildren() {
		setMass(getAverageChildMass());
	}

	void EvidOcTreeNode::updateOccupancyChildrenStamped(uint32_t timestamp) {
		this->updateMassChildren();
		this->setLogOdds(this->evidMassToLogOdds());
		this->setTimestamp(timestamp);
	}

	float EvidOcTreeNode::evidMassToLogOdds() const {
		// use Pignistic Transformation to convert evidential mass to occupancy probability
		return logodds((double) (mass.occu_() + mass.igno_() * 0.5f));
	}

	// tree implementation ------------------------------------------------------
	EvidOcTree::EvidOcTree(double in_resolution)
	: OccupancyOcTreeBase<EvidOcTreeNode>(in_resolution) {
		evidOcTreeMemberInit.ensureLinking();
	};

	bool EvidOcTree::pruneNode(EvidOcTreeNode* node) {
		if (!isNodeCollapsible(node))
			return false;
		// std::cout << "standing in the wrong place\n";
		// copy occupancy probability from one of its children
		node->value = getNodeChild(node, 0)->value;
		// set evid mass to average children's values
		node->updateMassChildren();
		
		// delete child
		for (unsigned int i=0;i<8;i++) {
      deleteNodeChild(node, i);
    }
    delete[] node->children;
    node->children = NULL;

    return true;
	}

	bool EvidOcTree::isNodeCollapsible(const EvidOcTreeNode* node) const{
    // all children must exist, must not have children of
    // their own and have the same occupancy probability
    if (!nodeChildExists(node, 0))
      return false;

    const EvidOcTreeNode* firstChild = getNodeChild(node, 0);
    if (nodeHasChildren(firstChild))
      return false;

    for (unsigned int i = 1; i<8; i++) {
      // compare nodes only using their occupancy, ignoring color for pruning
      if (!nodeChildExists(node, i) || nodeHasChildren(getNodeChild(node, i)) || !(getNodeChild(node, i)->getValue() == firstChild->getValue()))
        return false;
    }

    return true;
  }

	void EvidOcTree::insertPointCloud(const Pointcloud& scan, const octomap::point3d& sensor_origin, uint32_t timestamp,
                   double maxrange, bool lazy_eval, bool discretize) {
		KeySet free_cells, occupied_cells;
    if (discretize)
      computeDiscreteUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);
    else
      computeUpdate(scan, sensor_origin, free_cells, occupied_cells, maxrange);

    // insert data into tree  -----------------------
    for (KeySet::iterator it = free_cells.begin(); it != free_cells.end(); ++it) {
      updateNode(*it, false, timestamp, lazy_eval);
    }
    for (KeySet::iterator it = occupied_cells.begin(); it != occupied_cells.end(); ++it) {
      updateNode(*it, true, timestamp, lazy_eval);
    }

		/// with this evidential grid implementation lazy_eval is always on, need to call updateInnerOccupancy() 
		updateInnerOccupancy(timestamp);
	}

	EvidOcTreeNode* EvidOcTree::updateNode(const point3d& value, bool occupied, uint32_t timestamp, bool lazy_eval) {
		OcTreeKey key;
    if (!this->coordToKeyChecked(value, key))
      return NULL;
		// std::cout << "Key found. Proceed to updateNode(key, occupied, lazy_eval)\n";
    return updateNode(key, occupied, timestamp, lazy_eval);
	}

	EvidOcTreeNode* EvidOcTree::updateNode(const OcTreeKey& key, bool occupied, uint32_t timestamp, bool lazy_eval) {
		BasicBeliefAssignment bba_m2;
		if(occupied) {
			bba_m2.mo = bba_mo;
			bba_m2.mi = 1.0f - bba_mo;
		} else {
			// free
			bba_m2.mf = bba_mf;
			bba_m2.mi = 1.0f - bba_mf;
		}
		// std::cout << "bba_m2 is created. Proceed to updateNode(key, bba_m2, lazy_eval)\n";
		return updateNode(key, bba_m2, timestamp, lazy_eval);
	}

	EvidOcTreeNode* EvidOcTree::updateNode(const OcTreeKey& key, BasicBeliefAssignment bba_m2, uint32_t timestamp, bool lazy_eval) {
		bool createdRoot = false;
    if (this->root == NULL){
      this->root = new EvidOcTreeNode();
      this->tree_size++;
      createdRoot = true;
    }
		// std::cout << "root is created. Proceed to recursive\n";
    return updateNodeRecurs(this->root, createdRoot, key, 0, bba_m2, timestamp, lazy_eval);
	}

	EvidOcTreeNode* EvidOcTree::updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                           unsigned int depth, const BasicBeliefAssignment& bba_m2, uint32_t timestamp, bool lazy_eval) {
		bool created_node = false;
		assert(node);

    // follow down to last level
    if (depth < this->tree_depth) {
      unsigned int pos = computeChildIdx(key, this->tree_depth -1 - depth);
      if (!this->nodeChildExists(node, pos)) {
        // child does not exist, but maybe it's a pruned node?
        if (!this->nodeHasChildren(node) && !node_just_created ) {
          // current node does not have children AND it is not a new node
          // -> expand pruned node
					this->expandNode(node);
				} else {
          // not a pruned node, create requested child
					this->createNodeChild(node, pos);
					created_node = true;
        }
      }

      if (lazy_eval)
        return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, bba_m2, timestamp, lazy_eval);
      else {
        EvidOcTreeNode* retval = updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, bba_m2, timestamp, lazy_eval);
				// prune node if possible, otherwise set own probability
        // note: combining both did not lead to a speedup!
        // if (this->pruneNode(node)){
        //   // return pointer to current parent (pruned), the just updated node no longer exists
        //   return node;  // ptr to parent
        // } else {
				// 	std::cout << "Preparing updateOccupancyChildren()\n";
        //   node->updateOccupancyChildren();
				// 	std::cout << "Finishing updateOccupancyChildren()\n";
				// 	std::cout << "Standing here !!!\n";
				// 	return retval;
        // }
				return retval;
      }
    }
    // at last level, update node, end of recursion
    else {
			upadteNodeEvidMass(key, node, bba_m2, timestamp);
			EvidOcTreeNode* leaf = node;
      return leaf;  // ptr to leaf
    }
	}

	void EvidOcTree::upadteNodeEvidMass(const OcTreeKey& key, EvidOcTreeNode* evidNode, const BasicBeliefAssignment &bba_m2, uint32_t timestamp) {
		// compute decay factor
		uint32_t delta_t = timestamp - evidNode->getTimestamp();  // in millisecond 
		assert(delta_t >= 0);
		
		float alpha = exp(-((float) delta_t /1000.0f) / tau);
		assert((alpha >= 0.0f) && (alpha <=1.0f));

		// decay mass of current node
		evidNode->decayMass(alpha);

		// conjuctive fusion
		float m12_c = evidNode->mass.free_() * bba_m2.mo + evidNode->mass.occu_() * bba_m2.mf;
		float m12_f = evidNode->mass.free_() * bba_m2.mf + evidNode->mass.free_() * bba_m2.mi + evidNode->mass.igno_() * bba_m2.mf;
		float m12_o = evidNode->mass.occu_() * bba_m2.mo + evidNode->mass.occu_() * bba_m2.mi + evidNode->mass.igno_() * bba_m2.mo;
		float m12_i = evidNode->mass.igno_() * bba_m2.mi;
		assert(evidNode->mass.isMassValid());

		// detect cell with high conflict mass
		if(m12_c > this->thres_conf)
			this->conf_keys.push_back(key);
		// else
		// 	std::cout << "conflict mass = " << m12_c << "\n";

		// dempster normalization
		float k = 1.0f / (1.0f - m12_c);
		evidNode->setMass(0.0f, k*m12_f, k*m12_o, k*m12_i);

		// update node's logodds according to its evidential mass
		evidNode->setLogOdds(evidNode->evidMassToLogOdds());

		// clamping logodds
		if (evidNode->getLogOdds() < this->clamping_thres_min) {
			evidNode->setLogOdds(this->clamping_thres_min);
			return;
		}
		if (evidNode->getLogOdds() > this->clamping_thres_max) {
			evidNode->setLogOdds(this->clamping_thres_max);
		}
	}

	void EvidOcTree::updateInnerOccupancy(uint32_t timestamp) {
		this->updateInnerOccupancyRecurs(this->root, 0, timestamp);
	}

	void EvidOcTree::updateInnerOccupancyRecurs(EvidOcTreeNode* node, unsigned int depth, uint32_t timestamp){
		// only update inner node
		if(this->nodeHasChildren(node)) {
			// return early to last level
			if(depth < this->tree_depth) {
				for(unsigned int i=0; i<8; i++) {
					if(this->nodeChildExists(node, i))
						updateInnerOccupancyRecurs(getNodeChild(node, i), depth+1, timestamp);
				}
			}
			// end of recursion, reach leaf node
			node->updateOccupancyChildrenStamped(timestamp);
			// try to prune node 
			this->pruneNode(node);
		}
	}
	
	void EvidOcTree::publishMovingCells(uint32_t timestamp) {
		// iterate conf_keys, convert key to coordinate & print out coordinate
		std::cout << "[TIME = " << timestamp <<"] Moving cells cooridnates:\n";
		for(unsigned  i = 0; i < conf_keys.size(); i++) {
			point3d p_ = keyToCoord(conf_keys[i]);
			std::cout << "( " << p_.x() << ", " << p_.y() << ", " << p_.z() << " )\n"; 
		}
		// clear conf_keys to prepare for next measurement
		conf_keys.clear();
	}

	EvidOcTree::StaticMemberInitializer EvidOcTree::evidOcTreeMemberInit;

}//end namespace