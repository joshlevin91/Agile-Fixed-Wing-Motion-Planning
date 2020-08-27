#ifndef TYPES_H
#define TYPES_H

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

#endif
