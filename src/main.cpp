#include "MeshIntersection.h"
#include "igl/readOBJ.h"
#include <Eigen/Core>
#include <iostream>

int main()
{
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;

    if (!igl::readOBJ("bunny.obj", V, F))
    {
        if (!igl::readOBJ("../bunny.obj", V, F))
        {
            std::cerr << "Can't load bunny" << std::endl;
            return - 1;
        }
    }

    
    Eigen::MatrixXd V2 = V;
    Eigen::MatrixXi F2 = F;
    for (int i = 0; i < V2.rows(); i++)
    {
        V2(i, 0) += 1e-1;
    }

    int f1, f2;
    bool result;
    result = meshIntersection(V, F, V2, F, f1, f2);
    if (result)
        std::cout << "Faces " << f1 << " and " << f2 << " intersect" << std::endl;
    else
        std::cout << "No intersection" << std::endl;
}