# meshintersection

Computes whether or not two meshes intersect (using a kDOP BVH, with Moller's triangle-triangle intersection test as the narrow phase).
If the meshes intersect, face1 and face2 are set a pair of triangles that intersect. (If multiple triangles intersect, one pair is chosen arbitrarily.)

meshIntersection() in MeshIntersection.h is the function. See main.cpp for a usage example.
