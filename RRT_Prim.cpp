#include "mex.h"

#include <vector>
#include <stack>
#include <queue>
#include <math.h>
#include <time.h>
#include <algorithm>
#include <stdio.h>
#include <string.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <queue>

bool path_found;
world w; 
float start_coord[3];
std::queue<std::vector<float> > new_obs;
std::vector<std::vector<float> > found_obs; // Unknown obstacles found via laser range finder 

#define PI 3.1415f

const float V = 7.0f; // Flight velocity
const float bf = 3.0f; // Buffer around obstacles
const int max_tr_deg = 80; // Max turning rate
const int max_zr = 2; // Max climb rate
const float coll_dist = 0.2f; // Gap between collision checks
const float ATA_dist = 7.0f; // Distance required to do ATA
const float C2H_dist = 6.0f; // Distance required to do C2H
const float t_td = 0.23f; // Time-delay constant
const float t_int = 0.5f; // Planning horizon interval
const float ts = 0.5f; // Maximum time between nodes
const float ts_max = 1.0f; // Maximum time for a primitive
const int near_max = 10; // Maximum number of next near attempts
const int bias_freq = 40; // Frequency of selecting goal as "random" node
const float hdg_del_dome = 2.25147f; // Heading of x-axis of Stinger Dome relative to NED frame
const float ATA_t = 2.3462f; // ATA end time
const float C2H_x = 4.9336f; // C2H end x
const float C2H_z = -1.9029f; // C2H end z (the same sign as in CSV)
const float C2H_t = 2.04f; // C2H end time
const float H2C_x = 5.0f; // H2C end x
const float H2C_z = 0.52175f; // H2C end z (the same sign as in CSV)
const float H2C_t = 1.4177f; // H2C end time

// Tree
struct node{
	float coord[3];
	float hdg;
	int tr_deg;
	int zr;
	int type;
	float t;
	float cost;
	struct node* next;
	struct node* child;
	struct node* parent;
};

// Displacement
struct disp{
	float length;
	float cost;
};

// List for nearest sorting
struct list{
	node* n;
	float proximity;
	float nearness;
	float smart_nearness;
};

typedef float(*ptr_to_DCM)[3][3];

// World
struct world{
	float lims[3];
	int n_obs, n_goals;
	float** obs;
	float** goals;
	node* start;
	float offset[4];
	ptr_to_DCM DCM;
	bool two_D;
};

enum agile_man_t { ATA, C2H, H2C };

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, 
  const mxArray *prhs[])
{    int        i;
       
    /* Examine input (right-hand-side) arguments. */
    mexPrintf("\nThere are %d right-hand-side argument(s).", nrhs);
    for (i=0; i<nrhs; i++)  {
        mexPrintf("\n\tInput Arg %i is of type:\t%s ",i,mxGetClassName(prhs[i]));
    }
    
    /* Examine output (left-hand-side) arguments. */
    mexPrintf("\n\nThere are %d left-hand-side argument(s).\n", nlhs);
    if (nlhs > nrhs)
      mexErrMsgIdAndTxt( "MATLAB:mexfunction:inputOutputMismatch",
              "Cannot specify more outputs than inputs.\n");
    
    for (i=0; i<nlhs; i++)  {
        plhs[i]=mxCreateDoubleMatrix(1,1,mxREAL);
        *mxGetPr(plhs[i])=(double)mxGetNumberOfElements(prhs[i]);
    }
}

// Allocate dynamic memory and initialize new node
node* new_node(const float coord_x, const float coord_y, const float coord_z, const float hdg, const int tr_deg, const int zr,
	const int type, const float t, const float cost, node* const parent)
{
	node* new_node = (node*)malloc(sizeof(node));

	if (new_node) {
		new_node->coord[0] = coord_x;
		new_node->coord[1] = coord_y;
		new_node->coord[2] = coord_z;
		new_node->hdg = hdg;
		new_node->tr_deg = tr_deg;
		new_node->zr = zr;
		new_node->type = type;
		new_node->t = t;
		new_node->cost = cost;
		new_node->next = NULL;
		new_node->child = NULL;
		new_node->parent = parent;
	}

	return new_node;
}

// Generate random node (revisit because this will create garbage)
void rand_node_coord(node* rand_node, const node* const goal_node, const int iter)
{
	// Every bias_freq iterations, return goal node
	if (iter % bias_freq == 0) {
		rand_node->coord[0] = goal_node->coord[0];
		rand_node->coord[1] = goal_node->coord[1];
		rand_node->coord[2] = goal_node->coord[2];
	}
	else{
		rand_node->coord[0] = static_cast <float> (rand() % 30000) / 30000.0f * w.lims[0];
		rand_node->coord[1] = static_cast <float> (rand() % 30000) / 30000.0f * w.lims[1];
		rand_node->coord[2] = static_cast <float> (rand() % 30000) / 30000.0f * w.lims[2];
	}
}

// Add node to tree
node* add_child(node* existing_parent, node* new_child)
{
	if (existing_parent == NULL) return NULL;

	if (existing_parent->child) return add_sibling(existing_parent->child, new_child);

	else{
		new_child->parent = existing_parent;
		return existing_parent->child = new_child;
	}
}

// Add node as sibling if parent node already has a child
node* add_sibling(node* existing_child, node* new_sibling)
{
	if (existing_child == NULL) return NULL;

	while (existing_child->next) existing_child = existing_child->next;

	new_sibling->parent = existing_child->parent;
	return existing_child->next = new_sibling;
}

// Calculate trim states (coordinates and heading) at end of primitive
void trim_end_states(node* current, const node* const from, const int tr_deg, const int zr, const float dt)
{
	float tr_rad = tr_deg*PI / 180.0f;

	// If straight 
	if (tr_deg == 0){
		current->coord[0] = from->coord[0] + V*dt*cosf(from->hdg)*cosf(asinf(zr / V));
		current->coord[1] = from->coord[1] + V*dt*sinf(from->hdg)*cosf(asinf(zr / V));
		current->coord[2] = from->coord[2] + zr*dt;
	}
	// If turning
	else{
		current->coord[0] = from->coord[0] + (V / tr_rad*sinf(from->hdg + tr_rad*dt) - V / tr_rad*sinf(from->hdg))*cosf(asinf(zr / V));
		current->coord[1] = from->coord[1] + (V / -tr_rad*cosf(from->hdg + tr_rad*dt) + V / tr_rad*cosf(from->hdg))*cosf(asinf(zr / V));
		current->coord[2] = from->coord[2] + zr*dt;
	}

	// Calculate heading and keep between -PI and PI
	current->hdg = from->hdg + tr_rad*dt;
	limit_angle(current->hdg);
}

// Extend tree
node* extend_tree(node* root, node* const rand_node, node* const goal, const int gn)
{
	std::vector<list> near_vec;
	node* prim_root = NULL;
	node* prim_end = NULL;

	// Nearness criteria depends on whether feasible path has been found
	int n;
	if (!path_found) n = 0;
	else n = 0;

	// Get vector of nodes orderered by nearness
	near_vec = near(root, rand_node, n);

	int i = 0;
	int m = static_cast<int>(near_vec.size());

	// Look for a desirable primitive over a maximum number of nearest nodes
	while (i < std::min(near_max, m)){

		// Do not extend tree from a hover
		if (near_vec[i].n->type == 3){
			i++;
			continue;
		}

		// If the nearest node is within range of the final goal, only try C2H maneuver
		float D2G = sqrtf(powf(w.goals[w.n_goals - 1][0] - near_vec[i].n->coord[0], 2) + powf(w.goals[w.n_goals - 1][1] - near_vec[i].n->coord[1], 2));
		float C2H_range = C2H_x + V*t_td + w.goals[w.n_goals - 1][4];

		if (D2G <= C2H_range){
			prim_root = steer_agile(near_vec[i].n, C2H);
			// Add C2H maneuver if there is no collision and it ends in goal region
			if (w.two_D){
				D2G = sqrtf(powf(w.goals[w.n_goals - 1][0] - prim_root->child->coord[0], 2) + powf(w.goals[w.n_goals - 1][1] - prim_root->child->coord[1], 2));
			}
			else{
				D2G = sqrtf(powf(w.goals[w.n_goals - 1][0] - prim_root->child->coord[0], 2) + powf(w.goals[w.n_goals - 1][1] - prim_root->child->coord[1], 2) + powf(w.goals[w.n_goals - 1][2] - prim_root->child->coord[2], 2));
			}
			if (!collision(prim_root) && D2G <= w.goals[gn][4]){
				add_child(near_vec[i].n, prim_root);
				break;
			}
			else free_tree(&prim_root);
		}
		// Else, try normal primitive
		else{
			prim_root = steer_an(near_vec[i].n, rand_node);
			// If there is a collision on the first try, attempt ATA
			//if (collision(prim_root) && i == 0){
			if (collision(prim_root)){
				free_tree(&prim_root);
				prim_root = steer_agile(near_vec[i].n, ATA);
			}

			//If no collision, extend tree with primitive (unless it ends too close to goal region)
			prim_end = prim_root;
			while (prim_end->child != NULL){
				prim_end = prim_end->child;
			}

			D2G = sqrtf(powf(w.goals[w.n_goals - 1][0] - prim_end->coord[0], 2) + powf(w.goals[w.n_goals - 1][1] - prim_end->coord[1], 2));
			float ngr_2 = C2H_x + V*t_td - w.goals[w.n_goals - 1][4];

			if (!collision(prim_root) && D2G > ngr_2){
				add_child(near_vec[i].n, prim_root);
				break;
			}
			else free_tree(&prim_root);
		}

		i++;
	}

	return prim_root;
}

// Generate maneuver primitive steering towards node
node* steer_an(node* const from, node* const towards)
{
	int tr_deg, zr;
	float L;
	node* temp = NULL;

	// First node of primtive is time-delayed node
	node* td = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 1, from->t + t_td, 0.0f, from);
	trim_end_states(td, from, from->tr_deg, from->zr, t_td);
	disp d = disp_info(from, td);
	td->cost = from->cost + d.cost;

	// Find best trim primitive and time spent along it
	float r_G = norm(td, towards);
	float theta_G = atan2f(towards->coord[1] - td->coord[1], towards->coord[0] - td->coord[0]) - td->hdg;

	limit_angle(theta_G);

	float r_c = sqrtf(powf(towards->coord[0] - td->coord[0], 2) + powf(towards->coord[1] - td->coord[1], 2)) / (2 * sinf(theta_G));

	if (fabsf(theta_G) <= 0.01f) L = r_G; // if theta_G == 0
	else L = r_G*theta_G / sinf(theta_G); // if theta_G != 0

	float dt_max = std::min((L / V), ts_max); // Max time on primitive must not be greater than const ts_max
	float dt = std::min(ts, dt_max); // Time step must not be greater than dt_max;

	tr_deg = int(V / r_c * 180.0f / PI);
	tr_deg = ((tr_deg + 10 / 2) / 10) * 10; // Round turning rate to nearest 10 deg
	zr = int((towards->coord[2] - td->coord[2]) / dt_max);

	// Keep turning and climb/descent rate within allowable bounds
	if (tr_deg > max_tr_deg) tr_deg = max_tr_deg;
	else if (tr_deg < -max_tr_deg) tr_deg = -max_tr_deg;
	if (zr > max_zr) zr = max_zr;
	else if (zr < -max_zr) zr = -max_zr;

	// Assign primitive's turning and climb/descent rate to time-delayed node
	td->tr_deg = tr_deg;
	td->zr = zr;

	// Add intermediate nodes along trim primitive to tree
	node* last = td;
	while (dt < dt_max){
		temp = new_node(0.0f, 0.0f, 0.0f, 0.0f, tr_deg, zr, 0, td->t + dt, 0.0f, last);
		trim_end_states(temp, td, tr_deg, zr, dt);
		temp->cost = last->cost + disp_info(last, temp).cost;
		add_child(last, temp);

		last = temp;
		dt += ts;
	}

	// For dt = dt_max
	temp = new_node(0.0f, 0.0f, 0.0f, 0.0f, tr_deg, zr, 0, td->t + dt_max, 0.0f, last);
	trim_end_states(temp, td, tr_deg, zr, dt_max);
	temp->cost = last->cost + disp_info(last, temp).cost;
	add_child(last, temp);


	return td; // Time-delayed node is root of primitive
}

// Generate agile maneuver primitive 
node* steer_agile(node* const from, const agile_man_t agile_man)
{
	float dx_end{ 0 }, dy_end{ 0 }, dz_end{ 0 }, t_end{ 0 };
	int m_type = 0;

	switch (agile_man) {
	case ATA:
		m_type = 2;
		dx_end = 0.0f;
		dy_end = 0.0f;
		dz_end = 0.0f;
		t_end = ATA_t;
		break;
	case C2H:
		m_type = 3;
		dx_end = C2H_x;
		dy_end = 0.0f;
		dz_end = -C2H_z;
		t_end = C2H_t;
		break;
	case H2C:
		m_type = 4;
		dx_end = H2C_x;
		dy_end = 0.0f;
		dz_end = -H2C_z;
		//dz_end = 0.0f; //TEMP: for one plane flight
		t_end = H2C_t;
		break;
	default:
		printf("Unrecognized agile maneuver type.\n");
	}

	// Rotate to align with heading at from node
	ptr_to_DCM DCM = create_DCM(0.0f, 0.0f, -from->hdg);
	float temp[3] = { dx_end, dy_end, dz_end };
	rotate_arrvec(DCM, temp);
	dx_end = temp[0];
	dy_end = temp[1];
	dz_end = temp[2];
	free(DCM);
	DCM = NULL;

	// If H2C (i.e. first primitive)
	if (agile_man == H2C) {
		node* H2C_end = new_node(from->coord[0] + dx_end, from->coord[1] + dy_end, from->coord[2] + dz_end, from->hdg, 0, 0, m_type, from->t + t_end, 0.0f, from);
		return H2C_end;
	}

	// Root of primtive is time-delayed node
	node* td = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 1, from->t + t_td, 0.0f, from);
	trim_end_states(td, from, from->tr_deg, from->zr, t_td);
	td->cost = from->cost + disp_info(from, td).cost;

	// End node
	node* end = new_node(td->coord[0] + dx_end, td->coord[1] + dy_end, td->coord[2] + dz_end, 0.0f, 0, 0, m_type, td->t + t_end, 0.0f, td);
	end->cost = td->cost + 50.0f;
	if (agile_man == ATA){
		end->hdg = td->hdg + PI;
	}
	else {
		end->hdg = td->hdg;
	}

	limit_angle(end->hdg); // Make sure heading lies between -PI and PI

	add_child(td, end); // Append end node to time-delayed node

	return td; // Time-delayed node is root of primitive
}

// Update tree in accordance with aircraft's real-time motion
bool update_tree(node** root, node* const goal)
{
	static int count = 1;
	bool alg_end = false;
	bool large_pos_err = false;
	if(count % 5 == 0 && count <= 10){
		printf("large pos error\n");
		large_pos_err = true;
	}
	std::queue<node*> comm_q;

	update_tree_for_new_obstacles(root);

	// Get node nearest goal
	int n;
	if (!path_found) n = 0;
	else n = 1;
	std::vector<list> nearest_vec = near(*root, goal, n);
	node* end = nearest_vec[0].n;

	// Get stack of nodes after current root to that nearest node
	std::stack<node*> r2e_list = root_to_end(*root, end);

	// Get node at or just past time interval on way to nearest node
	node* comm_end = *root;
	while ((comm_end->t - (*root)->t) < t_int*2.0f){

		// Return if goal node has been reached
		float D2G;
		if (w.two_D){
			D2G = sqrtf(powf(w.goals[w.n_goals - 1][0] - comm_end->coord[0], 2) + powf(w.goals[w.n_goals - 1][1] - comm_end->coord[1], 2));
		}
		else
		{
			D2G = sqrtf(powf(w.goals[w.n_goals - 1][0] - comm_end->coord[0], 2) + powf(w.goals[w.n_goals - 1][1] - comm_end->coord[1], 2) + powf(w.goals[w.n_goals - 1][2] - comm_end->coord[2], 2));
		}
		if (D2G <= w.goals[w.n_goals - 1][4] && comm_end->type == 3){
			printf("Goal node reached :)\n");
			alg_end = true;
			break;
		}

		// If the current optimal node is reached, continue along children
		if (r2e_list.empty()){
			if (comm_end->child != NULL){
				comm_end = comm_end->child;
			}

			// If no children, algorithm has failed, go into hover
			else{
				printf("Caught up to trajectory end :(\n");

				// End with C2H
				node* C2H_root = steer_agile(comm_end, C2H);
				comm_q.push(C2H_root);
				comm_q.push(C2H_root->child);
				free_tree(&C2H_root);

				alg_end = true;
				break;
			}
		}
		else{
			comm_end = r2e_list.top();
			r2e_list.pop();
		}

		comm_q.push(comm_end);

		if (comm_end->type == 2) printf("ATA at [%.1f, %.1f, %.1f]\n", (double)comm_end->coord[0], (double)comm_end->coord[1], (double)comm_end->coord[2]);
		else if (comm_end->type == 3) printf("C2H at [%.1f, %.1f, %.1f]\n", (double)comm_end->coord[0], (double)comm_end->coord[1], (double)comm_end->coord[2]);

	}

	// Send commited nodes
	// Re-planning if position tracking error is large
	if (large_pos_err){
		// Add commited nodes starting from aircraft's actual position
		node* last = *root; // Note: should be node with actual position
		last->coord[0] += 1;
		last->coord[1] += 1;
		last->coord[2] += 0;
		node* next;
		while (!comm_q.empty()){
			next = comm_q.front();
			trim_end_states(next, last, next->tr_deg, next->zr, (next->t - last->t));
			add_to_commit_2(next);
			comm_q.pop();
			last = next;
		}
		// Free tree except for root node (node at end of commit)
		*root = next;
		(*root)->next = NULL;
		(*root)->child = NULL;
		(*root)->parent = NULL;
		free_tree_2(root, *root);
	}
	else {
		while (!comm_q.empty()){
			add_to_commit_2(comm_q.front());
			comm_q.pop();
		}
		// Prune tree
		prune_tree(root, comm_end);
		*root = comm_end;
		(*root)->parent = NULL;
	}

	++count;
	return alg_end;
}

// Displacement information for trajectory (length and cost)
disp disp_info(node* const from, node* const to)
{
	float L;

	float r = norm(from, to);
	float theta = atan2f(to->coord[1] - from->coord[1], to->coord[0] - from->coord[0]) - from->hdg;

	limit_angle(theta);

	if (fabsf(theta) <= 0.01f) L = r; // if theta == 0	
	else L = r*theta / sinf(theta); // if theta != 0

	//float xy_norm = sqrtf(powf(to->coord[0] - from->coord[0], 2) + powf(to->coord[1] - from->coord[1], 2));
	//float r_c = xy_norm / 2 * sinf(theta);
	//float fpa = atan2f(to->coord[2] - from->coord[2], xy_norm);

	disp d = { L, L }; // smoothness: 1/r_c + fabsf(fpa)

	return d;
}

// Check for collision
bool collision(node* const root)
{
	node* from = root;
	node* to = NULL;
	node* temp = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 1.0f, 0.0f, NULL);
	node* last = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 1.0f, 0.0f, NULL);
	float dt, dt_end;

	while (from->child != NULL){

		to = from->child;

		// Check if node is in world
		if (out_of_world(to)) {
			free(temp);
			free(last);
			return true;
		}

		dt = coll_dist / V; // Check for collision every coll_dist
		if (to->type < 2) dt_end = (to->t - from->t); // Non-agile
		else if (to->type == 2) dt_end = ATA_dist / V; // ATA
		else dt_end = C2H_dist / V; // C2H

		memcpy(last, from, sizeof(node));
		while (dt < dt_end){

			// Check if node collides with known obstacles or path to node intersects new found obstacles
			trim_end_states(temp, from, to->tr_deg, to->zr, dt);
			if (inside_object(temp) || intersects_new_found_obs(last, temp)){
				free(temp);
				free(last);
				return true;
			}

			memcpy(last, temp, sizeof(node));
			dt += coll_dist / V;
		}

		// Check end of path for collisions with known or new found obstacles
		trim_end_states(temp, from, to->tr_deg, to->zr, dt_end);
		if (inside_object(temp) || out_of_world(temp) || intersects_new_found_obs(last, temp)){
			free(temp);
			free(last);
			return true;
		}

		from = to;
	}

	free(temp);
	free(last);
	return false;
}

// Check if node is outside of world
bool out_of_world(node* const n)
{
	// For one-plane motion
	if (w.two_D){
		if (n->coord[0] < bf || n->coord[0] > (w.lims[0] - bf) ||
			n->coord[1] < bf || n->coord[1] > (w.lims[1] - bf)){
			return true;
		}
	}
	// For 3D motion
	else {
		if (n->coord[0] < bf || n->coord[0] > (w.lims[0] - bf) ||
			n->coord[1] < bf || n->coord[1] > (w.lims[1] - bf) ||
			n->coord[2] < bf || n->coord[2] > (w.lims[2] - bf)){
			return true;
		}
	}
	
	return false;
}

// Check if node is inside any of world's known objects
bool inside_object(node* const n)
{
	for (int i = 0; i < w.n_obs; i++){
		// For one-plane motion
		if (w.two_D){
			if (n->coord[0] >= (w.obs[i][0] - bf) &&
				n->coord[0] <= (w.obs[i][0] + w.obs[i][3] + bf) &&
				n->coord[1] >= (w.obs[i][1] - bf) &&
				n->coord[1] <= (w.obs[i][1] + w.obs[i][4] + bf)){
				return true;
			}
		}
		// For 3D motion
		else {
			if (n->coord[0] >= (w.obs[i][0] - bf) &&
				n->coord[0] <= (w.obs[i][0] + w.obs[i][3] + bf) &&
				n->coord[1] >= (w.obs[i][1] - bf) &&
				n->coord[1] <= (w.obs[i][1] + w.obs[i][4] + bf) &&
				n->coord[2] >= (w.obs[i][2] - bf) &&
				n->coord[2] <= (w.obs[i][2] + w.obs[i][5] + bf)){
				return true;
			}
		}
	}

	return false;
}

// Check if path between from and to nodes collides with new objects detected
bool intersects_new_found_obs(node* const from, node* const to)
{
	if (found_obs.empty()){
		return false;
	}

	for (int i = 0; i < found_obs.size(); ++i){
		if (intersection(from->coord[0], from->coord[1], to->coord[0], to->coord[1],
			found_obs[i][0], found_obs[i][1], found_obs[i][2], found_obs[i][3])){
			return true;
		}
	}

	return false;
}

// Deallocate memory of tree
void free_tree(node** n)
{
	if (*n == NULL)
		return;

	free_tree(&((*n)->child));
	free_tree(&((*n)->next));

	free(*n);
	*n = NULL;
}

// Deallocate memory of tree, except root
void free_tree_2(node** n, node* const root)
{
	if (*n == NULL)
		return;

	free_tree_2(&((*n)->child), root);
	free_tree_2(&((*n)->next), root);

	if (*n != root){
		free(*n);
		*n = NULL;
	}
}

// Prune tree (free memory of nodes that aren't in new_root's tree)
void prune_tree(node** root, node* const new_root)
{
	bool tbd = true;

	if (*root == NULL)
		return;

	if (*root == new_root) {
		tbd = false;
	}

	if (tbd){
		prune_tree(&(*root)->child, new_root);
		prune_tree(&(*root)->next, new_root);
		free(*root);
		*root = NULL;
	}
	else{
		prune_tree(&(*root)->next, new_root);
	}
}

// Create list of nodes in order of some nearness criteria
std::vector<list> near(node* const root, node* const to, const int near_type)
{
	std::vector<list> near_vec; // Vector of nodes with nearness information
	add_to_near_vec(root, to, &near_vec); // Add all nodes in tree, with nearness information, to vector
	if (near_type == 0) std::sort(near_vec.begin(), near_vec.end(), comp_near); // Sort via nearness criteria
	else if (near_type == 1) std::sort(near_vec.begin(), near_vec.end(), comp_prox); // Sort via proximity
	else std::sort(near_vec.begin(), near_vec.end(), comp_snear); // Sort via smart nearness criteria

	return near_vec;
}

// Euclidean distance comparison
bool comp_prox(const list &a, const list &b)
{
	return a.proximity < b.proximity;
}

// Nearness comparison
bool comp_near(const list &a, const list &b)
{
	return a.nearness < b.nearness;
}

// Smart nearness comparison
bool comp_snear(const list &a, const list &b)
{
	return a.smart_nearness < b.smart_nearness;
}

// Get size of tree
int tree_size(node* const root)
{
	if (root == NULL)
		return 0;
	else{
		return (tree_size(root->child) + 1 + tree_size(root->next));
	}
}

// Add nodes in tree, with nearness information, to vector
void add_to_near_vec(node* const n, node* const to, std::vector<list>* near_vec)
{
	if (n == NULL) return;

	// Nearness info is based on distance from time-delayed node to to node
	node* td = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 1, n->t + t_td, 0.0f, NULL);
	trim_end_states(td, n, n->tr_deg, n->zr, t_td);
	td->cost = n->cost + disp_info(n, td).cost;

	float D;
	if (w.two_D) {
		D = sqrtf(powf(n->coord[0] - start_coord[0], 2) + powf(n->coord[1] - start_coord[1], 2));
	}
	else{
		D = sqrtf(powf(n->coord[0] - start_coord[0], 2) + powf(n->coord[1] - start_coord[1], 2) + powf(n->coord[2] - start_coord[2], 2));
	}
	float L = n->cost;
	float D2G = norm(n, to);

	// Node address and nearness info
	list l;
	l.n = n;
	if (n->type == 3) { // C2H means node is in goal region
		l.proximity = 0.0f;
		l.nearness = 0.0f;
		l.smart_nearness = 0.0f;
	}
	else {
		l.proximity = D2G;
		l.nearness = disp_info(td, to).length;
		l.smart_nearness = (L + D2G) / D;
	}

	(*near_vec).push_back(l); // Add to vector

	free(td);
	td = NULL;

	// Recursively iterate through function for each node
	if (n->child != NULL) add_to_near_vec(n->child, to, near_vec);
	if (n->next != NULL) add_to_near_vec(n->next, to, near_vec);

	return;
}

// Create stack of nodes starting right after root and going to end
std::stack<node*> root_to_end(node* const root, node* const end)
{
	node* temp = end;
	std::stack<node*> r2e_dq;
	r2e_dq.push(temp); // Initialize deque with end node
	while (temp->parent != root && temp->parent != NULL){
		temp = temp->parent;
		r2e_dq.push(temp); // Add to front of deque
	}

	return r2e_dq;
}

// Add node to vector of committed nodes
void add_to_commit(const node* const n)
{
	/*

	float rel_pos[3] = { n->coord[0] - start_coord[0],
	n->coord[1] - start_coord[1],
	n->coord[2] - start_coord[2] };
	rotate_arrvec(w.DCM, rel_pos);
	float pos[3] = { rel_pos[0] + w.offset[0] + start_coord[0],
	rel_pos[1] + w.offset[1] + start_coord[1],
	rel_pos[2] + w.offset[2] + start_coord[2] };
	float hdg = n->hdg + w.offset[3];
	limit_angle(hdg);

	*/
}

// Add node to vector of committed nodes, converting from world frame to Stinger Dome
void add_to_commit_2(const node* const n)
{
	float pos[3] {n->coord[0], -n->coord[1], -n->coord[2]};

	rotate_arrvec(w.DCM, pos);

	float hdg = -w.offset[3] - n->hdg;

	limit_angle(hdg);

	// Write to CSV
	static bool first = true;
	static std::ofstream myCSV;

	if (first){
		myCSV.open("RRT_CSV.csv");
		first = false;
	}
	else {
		myCSV.open("RRT_CSV.csv", std::ofstream::out | std::ofstream::app);
	}

	myCSV << n->coord[0] << "," << n->coord[1] << "," << n->coord[2] << "," <<
		n->hdg << "," << n->tr_deg << "," << n->zr << "," << n->t << "," <<
		n->type << "\n";
	myCSV.close();

}

// Create Direction Cosine Matrix
ptr_to_DCM create_DCM(const float phi, const float the, const float psi)
{
	ptr_to_DCM DCM = (ptr_to_DCM)malloc(sizeof(*DCM));

	(*DCM)[0][0] = cosf(the)*cosf(psi);
	(*DCM)[0][1] = cosf(the)*sinf(psi);
	(*DCM)[0][2] = -sinf(the);

	(*DCM)[1][0] = -cosf(phi)*sinf(psi) + sinf(phi)*sinf(the)*cosf(psi);
	(*DCM)[1][1] = cosf(phi)*cosf(psi) + sinf(phi)*sinf(the)*sinf(psi);
	(*DCM)[1][2] = sinf(phi)*cosf(the);

	(*DCM)[2][0] = sinf(phi)*sinf(psi) + cosf(phi)*sinf(the)*cosf(psi);
	(*DCM)[2][1] = -sinf(phi)*cosf(psi) + cosf(phi)*sinf(the)*sinf(psi);
	(*DCM)[2][2] = cosf(phi)*cosf(the);

	return DCM;
}

// Calculate L2-norm between two nodes
float norm(node* const n1, node* const n2)
{
	return sqrtf(powf(n2->coord[0] - n1->coord[0], 2) + powf(n2->coord[1] - n1->coord[1], 2) + powf(n2->coord[2] - n1->coord[2], 2));
}

// Calculate L2-norm between two sets of 2D coordinates
float norm(const float p1_x, const float p1_y, const float p2_x, const float p2_y)
{
	return sqrtf(powf(p2_x - p1_x, 2) + powf(p2_y - p1_y, 2));
}

// Calculate L2-norm between two sets of positions
float norm(const float p1_x, const float p1_y, const float p1_z, const float p2_x, const float p2_y, const float p2_z)
{
	return sqrtf(powf(p2_x - p1_x, 2) + powf(p2_y - p1_y, 2) + powf(p2_z - p1_z, 2));
}

// Create world
void create_world(const int n)
{
	// Two corridor passages
	if (n == 0){
		
		w.two_D = false;

		// Starting node
		w.start = new_node(15.0f, 4.0f, 5.0f, PI, 0, 0, 0, 0.0f, 0.0f, NULL);

		// Environment boundaries
		w.lims[0] = 50.0f;
		w.lims[1] = 25.0f;
		w.lims[2] = 10.0f;

		// Obstacles
		w.n_obs = 2;
		float obs_arr[2][6] = { { 0.0f, 8.0f, 0.0f, 40.0f, 2.0f, 10.0f },
		{ 10.0f, 16.0f, 0.0f, 40.0f, 2.0f, 10.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		// Goals (intermediate and final)
		w.n_goals = 9;
		float goal_arr[9][6] = { { 45.0f, 8.5f, 5.0f, PI / 2.0f, 1.5f, PI },
		{ 30.0f, 13.0f, 5.0f, PI, 1.5f, PI },
		{ 20.0f, 13.0f, 5.0f, PI, 1.5f, PI },
		{ 10.0f, 13.0f, 5.0f, PI, 1.5f, PI },
		{ 5.0f, 16.5f, 5.0f, PI / 2.0f, 1.5f, PI },
		{ 10.0f, 21.0f, 5.0f, 0.0f, 1.5f, PI },
		{ 20.0f, 21.0f, 5.0f, 0.0f, 1.5f, PI },
		{ 30.0f, 21.0f, 5.0f, 0.0f, 1.5f, PI },
		{ 48.0f, 21.0f, 5.0f, 0.0f, 2.0f, PI } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// Two column obstacles
	else if (n == 1){

		w.two_D = false;

		w.start = new_node(3.0f, 3.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 50.0f;
		w.lims[1] = 25.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 2;
		float obs_arr[2][6] = { { 10.0f, 10.0f, 0.0f, 10.0f, 10.0f, 10.0f },
		{ 30.0f, 13.0f, 0.0f, 10.0f, 10.0f, 10.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 1;
		float goal_arr[1][6] = { 48.0f, 23.0f, 5.0f, PI / 2.0f, 2.0f, PI };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// Wall with passage
	else if (n == 2) {

		w.two_D = false;

		w.start = new_node(3.0f, 3.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 50.0f;
		w.lims[1] = 25.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 2;
		float obs_arr[2][6] = { { 24.0f, 0.0f, 0.0f, 2.0f, 10.0f, 10.0f },
		{ 24.0f, 15.0f, 0.0f, 2.0f, 10.0f, 10.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 2;
		float goal_arr[2][6] = { { 25.0f, 12.5f, 5.0f, 0.0f, 1.0f, PI },
		{ 48.0f, 23.0f, 5.0f, 0.0f, 2.0f, PI } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// S walls
	else if (n == 3) {

		w.two_D = false;

		w.start = new_node(3.0f, 3.0f, 5.0f, PI / 4, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 50.0f;
		w.lims[1] = 25.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 2;
		float obs_arr[2][6] = { { 15.0f, 0.0f, 0.0f, 5.0f, 20.0f, 10.0f },
		{ 30.0f, 5.0f, 0.0f, 5.0f, 20.0f, 10.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 3;
		float goal_arr[3][6] = { { 17.5f, 22.5f, 5.0f, 0.0f, 2.5f, PI },
		{ 32.5f, 2.5f, 5.0f, 0.0f, 2.5f, PI },
		{ 48.0f, 23.0f, 5.0f, 0.0f, 2.0f, PI } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// S walls in Dome
	else if (n == 4) {

		w.two_D = true;

		w.start = new_node(7.0f, 15.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 60.0f;
		w.lims[1] = 30.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 2;
		float obs_arr[2][6] = { { 22.0f, 15.0f, 0.0f, 2.0f, 15.0f, 10.0f },
		{ 38.0f, 0.0f, 0.0f, 2.0f, 15.0f, 10.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 1;
		float goal_arr[1][6] = { 50.0f, 10.0f, 5.0f, 0.0f, 2.5f, PI };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// U in Dome
	else if (n == 5) {

		w.two_D = true;

		w.start = new_node(5.0f, 5.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 60.0f;
		w.lims[1] = 30.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 1;
		float obs_arr[1][6] = { 0.0f, 15.0f, 0.0f, 30.0f, 1.0f, 10.0f };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 2;
		float goal_arr[2][6] = { { 40.0f, 15.0f, 5.0f, PI / 2.0f, 2.5f, PI },
		{ 7.0f, 21.0f, 5.0f, PI, 2.5f, 45.0f*PI / 180.0f } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// U in Dome that should force ATA
	else if (n == 6) {

		w.two_D = true;

		w.start = new_node(20.0f, 5.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 60.0f;
		w.lims[1] = 30.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 1;
		float obs_arr[1][6] = { 20.0f, 15.0f, 0.0f, 40.0f, 1.0f, 10.0f };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 2;
		float goal_arr[2][6] = { {5.0f, 5.0f, 5.0f, PI/2.0f, 2.5f, PI },
		{ 45.0f, 21.0f, 5.0f, 0.0f, 2.5f, 180.0f*PI / 180.0f } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// Big U corridors
	else if (n == 7) {

		w.two_D = true;

		w.start = new_node(6.0f, 50.0f, 5.0f, -PI / 2.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 100.0f;
		w.lims[1] = 100.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 1;
		float obs_arr[1][6] = { 12.0f, 0.0f, 0.0f, 76.0f, 88.0f, 10.0f };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 3;
		float goal_arr[3][6] = { { 6.0f, 95.0f, 5.0f, PI / 4.0f, 4.0f, PI },
								 { 94.0f, 95.0f, 5.0f, -PI / 4.0f, 4.0f, PI },
								 { 94.0f, 5.0f, 5.0f, -PI / 2.0f, 4.0f, PI } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// 3D Gaps
	else if (n == 8) {

		w.two_D = false;

		w.start = new_node(10.0f, 10.0f, 10.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 100.0f;
		w.lims[1] = 100.0f;
		w.lims[2] = 30.0f;

		w.n_obs = 6;
		float obs_arr[6][6] = { { 0.0f, 30.0f, 0.0f, 100.0f, 2.0f, 15.0f },
								{ 0.0f, 30.0f, 15.0f, 65.0f, 2.0f, 15.0f },
								{ 80.0f, 30.0f, 15.0f, 20.0f, 2.0f, 15.0f },
								{ 0.0f, 70.0f, 15.0f, 100.0f, 2.0f, 15.0f },
								{ 0.0f, 70.0f, 0.0f, 20.0f, 2.0f, 15.0f },
								{ 35.0f, 70.0f, 0.0f, 70.0f, 2.0f, 15.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 3;
		float goal_arr[3][6] = { { 72.5f, 31.0f, 22.5f, PI / 2.0f, 5.0f, PI },
								 { 27.5f, 71.0f, 7.5f, PI / 2.0f, 5.0f, PI },
								 { 25.0f, 95.0f, 15.0f, PI / 2.0f, 5.0f, PI } };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// Half of dome
	else if (n == 9) {

		w.two_D = false;

		w.start = new_node(5.0f, 5.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 60.0f;
		w.lims[1] = 60.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 2;
		float obs_arr[2][6] = { { 0.0f, 29.0f, 0.0f, 60.0f, 2.0f, 5.0f },
		{ 29.0f, 29.0f, 0.0f, 2.0f, 2.0f, 10.0f } };

		w.obs = (float**)malloc(w.n_obs * sizeof(float**));
		for (int i = 0; i < w.n_obs; i++){
			w.obs[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_obs; ++i){
			for (int j = 0; j < 6; ++j){
				w.obs[i][j] = obs_arr[i][j];
			}
		}

		w.n_goals = 1;
		float goal_arr[1][6] = { 55.0f, 55.0f, 5.0f, 0.0f, 3.5f, PI };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}

	// No obstacles
	else {

		w.two_D = true;

		w.start = new_node(5.0f, 15.0f, 5.0f, 0.0f, 0, 0, 0, 0.0f, 0.0f, NULL);

		w.lims[0] = 60.0f;
		w.lims[1] = 30.0f;
		w.lims[2] = 10.0f;

		w.n_obs = 0;

		w.n_goals = 1;
		float goal_arr[1][6] = { 40.0f, 15.0f, 5.0f, 0.0f, 2.5f, PI };

		w.goals = (float**)malloc(w.n_goals * sizeof(float**));
		for (int i = 0; i < w.n_goals; i++){
			w.goals[i] = (float*)malloc(6 * sizeof(float*));
		}

		for (int i = 0; i < w.n_goals; i++){
			for (int j = 0; j < 6; j++){
				w.goals[i][j] = goal_arr[i][j];
			}
		}
	}
}

// Check if goal region has been reached
bool goal_reached(node* const n, node* const goal, const int gn)
{
	if (n == NULL) return false;
	float d_hdg = n->hdg - w.goals[gn][3];
	// For one-plane motion
	if (w.two_D){
		if (gn == (w.n_goals - 1)){
			if (n->type == 3 && norm(n->coord[0], n->coord[1], goal->coord[0], goal->coord[1]) <= w.goals[gn][4] && fabsf(limit_angle(d_hdg)) <= w.goals[gn][5]) { // If in region
				return true;
			}
		}
		else {
			if (norm(n->coord[0], n->coord[1], goal->coord[0], goal->coord[1]) <= w.goals[gn][4] && fabsf(limit_angle(d_hdg)) <= w.goals[gn][5]) { // If in region
				return true;
			}
		}
	}
	// For 3D motion
	else {
		if (gn == (w.n_goals - 1)){
			if (n->type == 3 && norm(n, goal) <= w.goals[gn][4] && fabsf(limit_angle(d_hdg)) <= w.goals[gn][5]) { // If in region
				return true;
			}
		}
		else {
			if (norm(n, goal) <= w.goals[gn][4] && fabsf(limit_angle(d_hdg)) <= w.goals[gn][5]) { // If in region
				return true;
			}
		}
	}
	if (goal_reached(n->next, goal, gn)) return true; // Recursively iterate
	return goal_reached(n->child, goal, gn);
}

// Rotate array by DCM
void rotate_arrvec(ptr_to_DCM DCM, float arr_vec[])
{
	float tmp[3] = { arr_vec[0], arr_vec[1], arr_vec[2] };

	arr_vec[0] = (*DCM)[0][0] * tmp[0] + (*DCM)[0][1] * tmp[1] + (*DCM)[0][2] * tmp[2];
	arr_vec[1] = (*DCM)[1][0] * tmp[0] + (*DCM)[1][1] * tmp[1] + (*DCM)[1][2] * tmp[2];
	arr_vec[2] = (*DCM)[2][0] * tmp[0] + (*DCM)[2][1] * tmp[1] + (*DCM)[2][2] * tmp[2];
}

// Limit angle to -PI to PI range
float limit_angle(float &angle)
{
	angle = fmod(angle, 2 * PI);
	if (angle > PI) angle = -2 * PI + angle;
	else if (angle < -PI) angle = 2 * PI + angle;

	return angle;
}

// Get heading from euler array when aircraft is in hover
float hover_hdg(const float phi, const float the, const float psi)
{
	// Just return psi if not near hover
	if (fabsf(the - PI / 2.0f) > 0.35f) return psi;

	// Create DCM prime
	ptr_to_DCM DCM_p = (ptr_to_DCM)malloc(sizeof(*DCM_p));

	(*DCM_p)[0][0] = cosf(the)*cosf(psi);
	(*DCM_p)[1][0] = cosf(the)*sinf(psi);
	(*DCM_p)[2][0] = -sinf(the);

	(*DCM_p)[0][1] = -cosf(phi)*sinf(psi) + sinf(phi)*sinf(the)*cosf(psi);
	(*DCM_p)[1][1] = cosf(phi)*cosf(psi) + sinf(phi)*sinf(the)*sinf(psi);
	(*DCM_p)[2][1] = sinf(phi)*cosf(the);

	(*DCM_p)[0][2] = sinf(phi)*sinf(psi) + cosf(phi)*sinf(the)*cosf(psi);
	(*DCM_p)[1][2] = -sinf(phi)*cosf(psi) + cosf(phi)*sinf(the)*sinf(psi);
	(*DCM_p)[2][2] = cosf(phi)*cosf(the);

	float temp[3] = { 0.0f, 0.0f, 1.0f };
	rotate_arrvec(DCM_p, temp);

	float hdg = atan2f(temp[1], temp[0]);

	free(DCM_p);
	DCM_p = NULL;

	return hdg;
}

// Check if there is an intersection between two lines (p0--p1, p2--p3)
bool intersection(const float p0_x, const float p0_y, const float p1_x, const float p1_y,
	const float p2_x, const float p2_y, const float p3_x, const float p3_y)
{
	float d = norm(p2_x, p2_y, p3_x, p3_y);

	// Add buffer
	float p2_x_bf = p3_x + (p2_x - p3_x) + bf * (p2_x - p3_x) / d;
	float p2_y_bf = p3_y + (p2_y - p3_y) + bf * (p2_y - p3_y) / d;
	float p3_x_bf = p2_x + (p3_x - p2_x) + bf * (p3_x - p2_x) / d;
	float p3_y_bf = p2_y + (p3_y - p2_y) + bf * (p3_y - p2_y) / d;

	float s1_x, s1_y, s2_x, s2_y;
	s1_x = p1_x - p0_x;
	s1_y = p1_y - p0_y;
	s2_x = p3_x_bf - p2_x_bf;
	s2_y = p3_y_bf - p2_y_bf;

	float den = (-s2_x * s1_y + s1_x * s2_y);

	if (fabsf(den) < 0.0001f) return false;

	float s, t;
	s = (-s1_y * (p0_x - p2_x_bf) + s1_x * (p0_y - p2_y_bf)) / den;
	t = (s2_x * (p0_y - p2_y_bf) - s2_y * (p0_x - p2_x_bf)) / den;

	if (s >= 0 && s <= 1 && t >= 0 && t <= 1){
		return true; // Collision
	}

	return false; // No collision
}

// Update tree to take into account all new obstacles detected in real-time
void update_tree_for_new_obstacles(node** root)
{
	while (!new_obs.empty()){
		float d = norm((*root)->coord[0], (*root)->coord[1], new_obs.front()[0], new_obs.front()[1]) + 3.0f;
		float d2 = norm((*root)->coord[0], (*root)->coord[1], new_obs.front()[2], new_obs.front()[3]) + 3.0f;
		if (d2 > d) {
			d = d2;
		}

		prune_new_obs_collisions(root, root, d);
		found_obs.push_back(new_obs.front());
		new_obs.pop();
	}
}

// Delete tree nodes that collide with new obstacles
void prune_new_obs_collisions(node** n, node** root, const float d)
{
	if (*n == NULL){
		return;
	}

	// Do not keep searching node and children if it is further from root than furthest obstacle corner
	float dn = norm((*root)->coord[0], (*root)->coord[1], (*n)->coord[0], (*n)->coord[1]);
	bool past_obs = (dn > d) ? true : false;

	if (!past_obs && (*n)->child != NULL){
		if (intersection((*n)->coord[0], (*n)->coord[1], (*n)->child->coord[0], (*n)->child->coord[1],
			new_obs.front()[0], new_obs.front()[1], new_obs.front()[2], new_obs.front()[3])) {
			free_tree(n);
			return;
		}
	}

	if (!past_obs){
		prune_new_obs_collisions(&((*n)->child), root, d);
	}
	prune_new_obs_collisions(&((*n)->next), root, d);
}

// Create world and initialize tree
node* initialize_world(const int nw, const float p_init[3], const float hdg_init)
{
	// Create world
	create_world(nw);

	w.offset[3] = -hdg_del_dome;
	w.DCM = create_DCM(0.0f, 0.0f, w.offset[3]);

	start_coord[0] = w.start->coord[0];
	start_coord[1] = w.start->coord[1];
	start_coord[2] = w.start->coord[2];

	// Initialize tree with two nodes
	node* init = w.start;
	node* second = steer_agile(init, H2C); // start with H2C
	//node* second = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 1.0f, 0.0f, NULL);
	//trim_end_states(second, root, 0, 0, 0.5f);

	add_to_commit_2(init);
	add_to_commit_2(second);

	return second;
}

// Create world and initialize tree when in Dome
node* initialize_world_2(const int nw, const float p_init[3], const float hdg_init)
{
	// Create world
	create_world(nw);

	ptr_to_DCM DCM_init = create_DCM(0.0f, 0.0f, hdg_del_dome);
	float temp[3] {p_init[0], p_init[1], p_init[2]};
	rotate_arrvec(DCM_init, temp);
	free(DCM_init);

	w.start->coord[0] = temp[0];
	w.start->coord[1] = -temp[1];
	w.start->coord[2] = -temp[2];
	w.start->hdg = hdg_del_dome - hdg_init;
	w.offset[3] = -hdg_del_dome;
	w.DCM = create_DCM(0.0f, 0.0f, w.offset[3]);

	start_coord[0] = w.start->coord[0];
	start_coord[1] = w.start->coord[1];
	start_coord[2] = w.start->coord[2];

	printf("p_init in map: [%.1f, %.1f, %.1f]\n", start_coord[0], start_coord[1], start_coord[2]);

	// Initialize tree with two nodes
	node* init = w.start;
	node* second = steer_agile(init, H2C); // start with H2C
	//node* second = new_node(0.0f, 0.0f, 0.0f, 0.0f, 0, 0, 0, 1.0f, 0.0f, NULL);
	//trim_end_states(second, root, 0, 0, 0.5f);

	add_to_commit_2(init);
	add_to_commit_2(second);

	return second;
}

// Garbage cleanup for tree and world dynamic memory
void cleanup_tree_and_world(node** root)
{
	// Free tree
	free_tree(root);

	// Free world
	for (int i = 0; i < w.n_obs; i++){
		free(w.obs[i]);
	}

	free(w.obs);
	for (int i = 0; i < w.n_goals; i++){
		free(w.goals[i]);
	}

	free(w.goals);
	free(w.DCM);
}

void quat2eul(const float q[4], float eul[3])
{
	eul[0] = atan2f(2 * (q[0] * q[1] + q[2] * q[3]), 1 - 2 * (q[1] * q[1] + q[2] * q[2]));
	eul[1] = asinf(2 * (q[0] * q[2] - q[3] * q[1]));
	eul[2] = atan2f(2 * (q[0] * q[3] + q[1] * q[2]), 1 - 2 * (q[2] * q[2] + q[3] * q[3]));
}

