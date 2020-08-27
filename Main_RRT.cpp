#include "RRT_Prim.h"
#include "common.h"
#include <iostream>
  
const int max_tree_size = 100000; // Max tree size (to keep memory usage low)

bool path_found = false; // Whether a feasible path has been found yet
world w; // Environment information
float start_coord[3]; // Starting position

void get_user_input(int &nw);
node* RRT_loop(node* root); 

int main()
{
	// Initial pose (faux sensor measurements)
	float hdg_init = PI/2.0f;
	float p_init[3];
	p_init[0] = 0.0f;
	p_init[1] = 0.0f;
	p_init[2] = 5.0f;

	// User selected
	int nw;

	// Get user input for world selection
	get_user_input(nw); 

	// Initialize world
	node* root = initialize_world(nw, p_init, hdg_init);

	// RRT algorithm
	root = RRT_loop(root);

	// Garbage collection (tree and world memory)
	cleanup_tree_and_world(&root);
}

void get_user_input(int &nw)
{
	// World number
	printf("Which world? ");
	std::cin >> nw;
}

node* RRT_loop(node* root)
{
	clock_t t1, t2;
	bool done = false;
	float elapsed;
	int iter = 0;
	int gn = 0;

	node* prim_root;
	node* goal = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);
	goal->coord[0] = w.goals[gn][0]; goal->coord[1] = w.goals[gn][1]; goal->coord[2] = w.goals[gn][2];
	node* randn = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);
	srand(static_cast<int>(time(NULL)));

	while (!done){

		// Start timer
		t1 = clock();
		bool done_iter = false;

		while (!done_iter){

			// Keep track of real-time computation interval
			t2 = clock();
			elapsed = (float)(t2 - t1) / static_cast<float>(CLOCKS_PER_SEC);
			if (path_found && elapsed >= t_int){
				done_iter = true;
				break;
			}

			if (tree_size(root) < max_tree_size){

				// Generate a random node (change coordinates of random node)
				rand_node_coord(randn, goal, iter);

				// Extend tree
				prim_root = extend_tree(root, randn, goal, gn);

				// Check if current goal has been reached
				if (goal_reached(prim_root, goal, gn)){
					// If intermediate goal, update goal node
					if (gn < (w.n_goals - 1)){
						gn++;
						goal->coord[0] = w.goals[gn][0];
						goal->coord[1] = w.goals[gn][1];
						goal->coord[2] = w.goals[gn][2];
						//printf("Intermediate goal %d of %d reached!\n", gn, w.n_goals - 1);
					}
					// If final goal, feasible path has been found
					else if (!path_found){
						path_found = true;
						printf("Feasible path found in %.5f seconds!\n", elapsed);
					}
				}
			}

			iter++;
		}

		// Print number of iterations and size of tree
		printf("Number of iterations: %d,  Size of tree: %d\n", iter, tree_size(root));

		// Update tree in accordance with aircraft's real-time motion
		done = update_tree(&root, goal);
	}

	// Garbage collection
	free(randn);
	free(goal);

	return root;
}