#pragma once

#include "Triangle.h"

FB_DECLARE_STRUCT(assignment, Ray)

FB_PACKAGE1(assignment)

enum node_axis
{
	axis_x,
	axis_y,
	axis_z,
	axis_none,
};

class KDNode
{
public:
	struct t_box3d
	{
		math::VC3 size;
		math::VC3 xyz_min;
		math::VC3 xyz_max;
	};
	t_box3d bounding_box;
	node_axis axis = axis_none;
	PodVector<Triangle> triangles;
	PodVector<Triangle> left_tris;
	PodVector<Triangle> right_tris;
	KDNode* left = nullptr;
	KDNode* right = nullptr;
	math::VC3 center;
	math::VC3 mid_point;
	uint32_t depth = 0;
	void boundingBoxSet();
	void initBoundingBox();
	node_axis boundingBoxLongestAxis();
	void boundingBoxSizeSet();
	void splitTriangles();
	void triangleVecMidpoint();
	bool boundingBoxRayHit(Ray& ray, math::VC3 ray_dir_inverse) const;
	~KDNode();
};

class KDTree
{
public:
	void createTree(PodVector<Triangle> &polygons);
	KDNode* createTreeRecursive(PodVector<Triangle> &triangles, uint32_t depth);
	KDNode* createNode(PodVector<Triangle> &triangles);
	bool treeRayHits(Ray &ray) const;
	bool treeRayHitsRecursive(KDNode* node, Ray& ray, math::VC3 ray_dir_inverse) const;
	bool findIntersection(KDNode* node, Ray& ray) const;
	uint32_t min_kd_node_triangles = 4;
	uint32_t max_kd_tree_depth = 16;
	KDNode* root = nullptr;
	~KDTree();
};

FB_END_PACKAGE1()
