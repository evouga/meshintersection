#ifndef APPROACHDISTANCE_H
#define APPROACHDISTANCE_H

#include <Eigen/Core>
#include <vector>

/*
 * Computes whether or not two meshes intersect (using a kDOP BVH, with 
 * Moller's triangle-triangle intersection test as the narrow phase).
 * If the meshes intersect, face1 and face2 are set a pair of triangles 
 * that intersect. (If multiple triangles intersect, one pair is chosen
 * arbitrarily.)
 */

bool meshIntersection(
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2,
    const Eigen::MatrixXi& F2,
    int& face1,
    int& face2);

#endif