#include <octomap/EvidOcTree.h>
#include <octomap/octomap_utils.h>

template<typename T>
bool isNull(T x){return fabs(x) < 1e-5;}

template<typename T>
bool isEqual(T x, T y){return isNull(x-y);}

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

	void EvidOcTreeNode::updateMassChildren(){
		setMass(getAverageChildMass());
	}

	// tree implementation ------------------------------------------------------
	EvidOcTree::EvidOcTree(double in_resolution)
	: OccupancyOcTreeBase<EvidOcTreeNode>(in_resolution) {
		evidOcTreeMemberInit.ensureLinking();
	};

	bool EvidOcTree::pruneNode(EvidOcTreeNode* node) {
		if (!isNodeCollapsible(node))
			return false;

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

	bool EvidOcTree::isNodeCollapsible(const EvidOcTreeNode* node) const {
		// all children must exist & not have children of their own
		// while having the same occupancy probability
		if(!nodeHasChildren(node))
			return false;
		
		const EvidOcTreeNode* firstChild = getNodeChild(node, 0);
		if(nodeHasChildren(firstChild))
			return false;

		for(unsigned int i=1; i<8; i++){
			if(!nodeChildExists(node, i))
				return false;
			const EvidOcTreeNode* child = getNodeChild(node, i);
			// compare node using their occupancy only while ignoring evid mass
			if(nodeHasChildren(child) || !isEqual(child->getValue(), firstChild->getValue()))
				return false;
		}

		return true;
	}

	EvidOcTreeNode* EvidOcTree::updateNode(const OcTreeKey& key, bool occupied, bool lazy_eval) {
		BasicBeliefAssignment bba_m2;
		if(occupied) {
			bba_m2.mo = bba_mo;
			bba_m2.mi = 1.0f - bba_mo;
		} else {
			// free
			bba_m2.mf = bba_mf;
			bba_m2.mi = 1.0f - bba_mf;
		}
		return updateNode(key, bba_m2, lazy_eval);
	}

	EvidOcTreeNode* EvidOcTree::updateNode(const OcTreeKey& key, BasicBeliefAssignment bba_m2, bool lazy_eval) {
		bool createdRoot = false;
    if (this->root == NULL){
      this->root = new EvidOcTreeNode();
      this->tree_size++;
      createdRoot = true;
    }

    return updateNodeRecurs(this->root, createdRoot, key, 0, bba_m2, lazy_eval);
	}

	EvidOcTreeNode* EvidOcTree::updateNodeRecurs(EvidOcTreeNode* node, bool node_just_created, const OcTreeKey& key,
                           unsigned int depth, const BasicBeliefAssignment& bba_m2, bool lazy_eval) {
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
        }
        else {
          // not a pruned node, create requested child
          this->createNodeChild(node, pos);
          created_node = true;
        }
      }

      if (lazy_eval)
        return updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, bba_m2, lazy_eval);
      else {
        EvidOcTreeNode* retval = updateNodeRecurs(this->getNodeChild(node, pos), created_node, key, depth+1, bba_m2, lazy_eval);
        // prune node if possible, otherwise set own probability
        // note: combining both did not lead to a speedup!
        if (this->pruneNode(node)){
          // return pointer to current parent (pruned), the just updated node no longer exists
          retval = node;
        } else{
          node->updateOccupancyChildren();
        }

        return retval;
      }
    }
    // at last level, update node, end of recursion
    else {
      upadteNodeEvidMass(key, node, bba_m2);
      return node;
    }
	}

	void EvidOcTree::upadteNodeEvidMass(const OcTreeKey& key, EvidOcTreeNode* evidNode, const BasicBeliefAssignment &bba_m2) {
		// conjuctive fusion
		float m12_c = evidNode->mass.free_() * bba_m2.mo + evidNode->mass.occu_() * bba_m2.mf;
		float m12_f = evidNode->mass.free_() * bba_m2.mf + evidNode->mass.free_() * bba_m2.mi + evidNode->mass.igno_() * bba_m2.mf;
		float m12_o = evidNode->mass.occu_() * bba_m2.mo + evidNode->mass.occu_() * bba_m2.mi + evidNode->mass.igno_() * bba_m2.mo;
		float m12_i = evidNode->mass.igno_() * bba_m2.mi;
		
		// detect cell with high conflict mass
		if(m12_c > this->thres_conf)
			this->conf_keys.push_back(key);

		// dempster normalization
		float k = 1.0f / (1.0f - m12_c);
		evidNode->setMass(0.0f, k*m12_f, k*m12_o, k*m12_i);

		// update node's logodds value
		float occuProb = k*(m12_o + m12_i*0.5f);
		evidNode->setLogOdds(logodds((double) occuProb));
		if (evidNode->getLogOdds() < this->clamping_thres_min) {
			evidNode->setLogOdds(this->clamping_thres_min);
			return;
		}
		if (evidNode->getLogOdds() > this->clamping_thres_max) {
			evidNode->setLogOdds(this->clamping_thres_max);
		}
	}

	//TODO: A function to search for node with high conflict

}//end namespace