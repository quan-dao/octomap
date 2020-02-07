#include <octomap/octomap.h>
#include <octomap/EvidOcTree.h>

using namespace std;
using namespace octomap;

void print_query_info(point3d query, EvidOcTreeNode* node) {
  if (node != NULL) {
		
    cout << "Infor at " << query << ":\n";
		cout << "Evidential mass: "; 
		node->printMass(); 
    cout << "occupancy probability :\t " << node->getOccupancy() << endl;
  }
  else 
    cout << "occupancy probability at " << query << ":\t is unknown" << endl;   
  cout << "---------------------------\n"; 
}

int main(int argc, char** argv) {
	cout << endl;
  cout << "generating example map" << endl;

	EvidOcTree tree(0.1);  // create an empty tree with resolution 0.1;
	cout << "Tree is initialized\n";

	// // test 1 point
	// float x=-20*0.05, y=-20*0.05f, z=-20*0.05f;
	// point3d endpoint(x, y, z);
	// cout << "End point is created\n";
	// tree.updateNode(endpoint, true, false);
	// cout << "Finish updating tree\n";

  // insert some measurements of occupied cells
  for (int x=-20; x<20; x++) {
    for (int y=-20; y<20; y++) {
      for (int z=-20; z<20; z++) {
        point3d endpoint ((float) x*0.05f, (float) y*0.05f, (float) z*0.05f);
        tree.updateNode(endpoint, true); // integrate 'occupied' measurement
      }
    }
  }
	cout << "Occupied cells are integrated \n";

	// insert some measurements of free cells
  for (int x=-30; x<30; x++) {
    for (int y=-30; y<30; y++) {
      for (int z=-30; z<30; z++) {
        point3d endpoint ((float) x*0.02f-1.0f, (float) y*0.02f-1.0f, (float) z*0.02f-1.0f);
        tree.updateNode(endpoint, false);  // integrate 'free' measurement
      }
    }
  }

  cout<<"Recursively update inner node & do pruning\n";
  tree.updateInnerOccupancy();
  cout<<"finish recursively update inner node & do pruning\n";

	cout << endl;
  cout << "performing some queries:" << endl;
  
  point3d query (0., 0., 0.);
  EvidOcTreeNode* result = tree.search (query);
  print_query_info(query, result);

  query = point3d(-0.5,-0.5,-0.5);
  result = tree.search (query);
  print_query_info(query, result);

  query = point3d(-1.,-1.,-1.);
  result = tree.search (query);
  print_query_info(query, result);

  query = point3d(0.5,0.5,0.5);
  result = tree.search (query);
  print_query_info(query, result);

  query = point3d(1.,1.,1.);
  result = tree.search (query);
  print_query_info(query, result);


  cout << endl;
  tree.writeBinary("simple_tree_evid.bt");
  cout << "wrote example file simple_tree_evid.bt" << endl << endl;
  cout << "now you can use octovis to visualize: octovis simple_tree_evid.bt"  << endl;
  cout << "Hint: hit 'F'-key in viewer to see the freespace" << endl  << endl;  

	return 0;
}
