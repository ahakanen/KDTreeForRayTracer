#include "Precompiled.h"
#include "KDTree.h"
#include "Assignment/Ray.h"
#include "fb/lang/NumericLimits.h"

FB_PACKAGE1(assignment)

/* Create the kd tree recursively */

void KDTree::createTree(PodVector<Triangle> &polygons)
{
	root = createTreeRecursive(polygons, 0);
}

KDNode *KDTree::createTreeRecursive(PodVector<Triangle> &triangles, uint32_t depth)
{
	KDNode *node;

	node = createNode(triangles);
	node->depth = depth;
	if (triangles.getSize() < min_kd_node_triangles || depth >= max_kd_tree_depth)
		return (node);
	node->axis = node->boundingBoxLongestAxis();
	node->splitTriangles();
	node->left = createTreeRecursive(node->left_tris, depth + 1);
	node->right = createTreeRecursive(node->right_tris, depth + 1);
	return (node);
}

/* Create a single kd tree node */

KDNode *KDTree::createNode(PodVector<Triangle> &triangles)
{
	KDNode* node;

	node = new KDNode();
	node->triangles = triangles;
	node->boundingBoxSet();
	return (node);
}

/* Precalculate ray direction inverse for later use and start searching the tree recursively */

bool KDTree::treeRayHits(Ray &ray) const
{
	bool result;
	math::VC3 ray_dir_inverse;

	ray_dir_inverse.x = 1.0f / ray.direction.getNormalized().x;
	ray_dir_inverse.y = 1.0f / ray.direction.getNormalized().y;
	ray_dir_inverse.z = 1.0f / ray.direction.getNormalized().z;
	result = treeRayHitsRecursive(root, ray, ray_dir_inverse);
	return (result);
}

/* Search through the tree, then check the triangles in the last node for intersections */

bool KDTree::treeRayHitsRecursive(KDNode* node, Ray &ray, math::VC3 ray_dir_inverse) const
{
	bool	hits_right;
	bool	hits_left;

	if (node && node->boundingBoxRayHit(ray, ray_dir_inverse))
	{
		if (node->left || node->right)
		{
			hits_left = treeRayHitsRecursive(node->left, ray, ray_dir_inverse);
			hits_right = treeRayHitsRecursive(node->right, ray, ray_dir_inverse);
			return (hits_left || hits_right);
		}
		return (findIntersection(node, ray));
	}
	return (false);
}

/* Search for triangle intersections in a similar manner to the brute force method */

bool KDTree::findIntersection(KDNode *node, Ray &ray) const
{
	bool somethingFound;

	somethingFound = false;
	for (SizeType i = 0; i < node->triangles.getSize(); ++i)
	{
		const Triangle& tri = node->triangles[i];
		if (ray.startingPoint != &tri && tri.intersects(ray) == 1)
		{
			/* We have a hit */
			somethingFound = true;
			if (ray.intersectionType == IntersectionTypeAny)
				break;
		}
	}
	return (somethingFound);
}

/* Define the bounding box of node, which will used when searching for intersections */

void KDNode::boundingBoxSet()
{
	float min;
	float max;

	initBoundingBox();
	for (SizeType i = 0; i < triangles.getSize(); i++)
	{
		for (SizeType j = 0; j < 3; j++)
		{
			min = std::fmin(triangles[i].vertices[0]->v[j], std::fmin(triangles[i].vertices[1]->v[j], triangles[i].vertices[2]->v[j]));
			bounding_box.xyz_min.v[j] = std::fmin(bounding_box.xyz_min.v[j], min);
			max = std::fmax(triangles[i].vertices[0]->v[j], std::fmax(triangles[i].vertices[1]->v[j], triangles[i].vertices[2]->v[j]));
			bounding_box.xyz_max.v[j] = std::fmax(bounding_box.xyz_max.v[j], max);
		}
	}
	boundingBoxSizeSet();
}

/* Initialize the bounding box */

void KDNode::initBoundingBox()
{
	bounding_box.xyz_max.x = lang::NumericLimits<float>::getLowest();
	bounding_box.xyz_max.y = lang::NumericLimits<float>::getLowest();
	bounding_box.xyz_max.z = lang::NumericLimits<float>::getLowest();
	bounding_box.xyz_min.x = lang::NumericLimits<float>::getHighest();
	bounding_box.xyz_min.y = lang::NumericLimits<float>::getHighest();
	bounding_box.xyz_min.z = lang::NumericLimits<float>::getHighest();
	if (triangles.getSize() == 0)
	{
		bounding_box.xyz_min.x = 0;
		bounding_box.xyz_min.y = 0;
		bounding_box.xyz_min.z = 0;
		bounding_box.xyz_max.x = 0;
		bounding_box.xyz_max.y = 0;
		bounding_box.xyz_max.z = 0;
		bounding_box.size.x = 0;
		bounding_box.size.y = 0;
		bounding_box.size.z = 0;
	}
}

/* Define the size of a bounding box */

void KDNode::boundingBoxSizeSet()
{
	bounding_box.size.x = bounding_box.xyz_max.x - bounding_box.xyz_min.x;
	bounding_box.size.y = bounding_box.xyz_max.y - bounding_box.xyz_min.y;
	bounding_box.size.z = bounding_box.xyz_max.z - bounding_box.xyz_min.z;
}

/* Find the longest axis of a bounding box */

node_axis KDNode::boundingBoxLongestAxis()
{
	float longest;

	longest = std::fmax(bounding_box.size.x, std::fmax(bounding_box.size.y, bounding_box.size.z));
	if (bounding_box.size.x == longest)
		return (axis_x);
	else if (bounding_box.size.y == longest)
		return (axis_y);
	else
		return (axis_z);
}

/* Split triangles in two based on if they're on left or right side of node axis */

void KDNode::splitTriangles()
{
	float triCenter = 0;

	triangleVecMidpoint();
	for (SizeType i = 0; i < triangles.getSize(); i++)
	{
		triCenter = (triangles[i].vertices[0]->v[axis] + triangles[i].vertices[1]->v[axis] + triangles[i].vertices[2]->v[axis]) / 3.0f;
		if (mid_point.v[axis] >= triCenter)
			right_tris.pushBack(triangles[i]);
		else
			left_tris.pushBack(triangles[i]);
	}
}

/* Find the "center of mass" of the triangles in the node */

void KDNode::triangleVecMidpoint()
{
	center.x = 0;
	center.y = 0;
	center.z = 0;
	mid_point.x = 0;
	mid_point.y = 0;
	mid_point.z = 0;
	for (SizeType i = 0; i < triangles.getSize(); i++)
	{
		center.x = (triangles[i].vertices[0]->x + triangles[i].vertices[1]->x + triangles[i].vertices[2]->x) / 3.0f;
		center.y = (triangles[i].vertices[0]->y + triangles[i].vertices[1]->y + triangles[i].vertices[2]->y) / 3.0f;
		center.z = (triangles[i].vertices[0]->z + triangles[i].vertices[1]->z + triangles[i].vertices[2]->z) / 3.0f;
		mid_point.x += center.x / (float)triangles.getSize();
		mid_point.y += center.y / (float)triangles.getSize();
		mid_point.z += center.z / (float)triangles.getSize();
	}
}

/* Determine if the ray hits a bounding box 
   https://tavianator.com/2011/ray_box.html
*/

bool KDNode::boundingBoxRayHit(Ray& ray, math::VC3 ray_dir_inverse) const
{
	float	t[8];

	t[0] = (bounding_box.xyz_min.x - ray.origin.x) * ray_dir_inverse.x;
	t[1] = (bounding_box.xyz_max.x - ray.origin.x) * ray_dir_inverse.x;
	t[2] = (bounding_box.xyz_min.y - ray.origin.y) * ray_dir_inverse.y;
	t[3] = (bounding_box.xyz_max.y - ray.origin.y) * ray_dir_inverse.y;
	t[4] = (bounding_box.xyz_min.z - ray.origin.z) * ray_dir_inverse.z;
	t[5] = (bounding_box.xyz_max.z - ray.origin.z) * ray_dir_inverse.z;
	t[6] = std::fmax(std::fmax(std::fmin(t[0], t[1]), std::fmin(t[2], t[3])), std::fmin(t[4], t[5]));
	t[7] = std::fmin(std::fmin(std::fmax(t[0], t[1]), std::fmax(t[2], t[3])), std::fmax(t[4], t[5]));
	if (t[7] < 0 || t[6] > t[7])
		return (false);
	return (true);
}

KDTree::~KDTree()
{
	if (root)
		delete root;
}

KDNode::~KDNode()
{
	if (left)
		delete left;
	if (right)
		delete right;
}

FB_END_PACKAGE1()
