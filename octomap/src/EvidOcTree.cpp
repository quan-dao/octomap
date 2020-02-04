#include <octomap/EvidOcTree.h>

namespace octomap {

	// node implementation ------------------------------------------------------
	EvidOcTreeNode::EvidMass EvidOcTreeNode::getAverageChildMass() const {
		float _mc = 0, _mf = 0, _mo = 0, _mi = 0;
		int num_child = 0;

		if (children != NULL){
			for (int i=0; i<8; i++){
				EvidOcTreeNode* child = static_cast<EvidOcTreeNode*>(children[i]);
				if (child != NULL && child->isMassSet()){
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


}//end namespace