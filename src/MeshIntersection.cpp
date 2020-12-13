#include "MeshIntersection.h"
#include "TriTriIntersection.h"
#include <iostream>
#include <algorithm>

struct DOPNode
{
    DOPNode* left, * right;    
    int face;
    double lo[13];
    double hi[13];    
    const static double DOPvecs[];

    DOPNode() : left(NULL), right(NULL) {}
    ~DOPNode() { delete left; delete right; }
};

const double DOPNode::DOPvecs[] = { 
    1, 0, 0,
    0, 1, 0,
    0, 0, 1,
    1, 1, 0,
    1, -1, 0,
    1, 0, 1,
    1, 0, -1,
    0, 1, 1,
    0, 1, -1,
    1, 1, 1,
    1, 1, -1,
    1, -1, 1,
    1, -1, -1 
};

DOPNode* fitFace(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, int face)
{
    DOPNode* result = new DOPNode;    
    result->face = face;    
    for (int i = 0; i < 13; i++)
    {
        result->lo[i] = std::numeric_limits<double>::infinity();
        result->hi[i] = -std::numeric_limits<double>::infinity();
    }

    for (int i = 0; i < 3; i++)
    {
        Eigen::Vector3d pt = V.row(F(face, i)).transpose();
        for (int j = 0; j < 13; j++)
        {
            Eigen::Vector3d DOPvec;
            for (int k = 0; k < 3; k++)
                DOPvec[k] = DOPNode::DOPvecs[3 * j + k];
            double coord = pt.dot(DOPvec) / DOPvec.norm();
            result->lo[j] = std::min(result->lo[j], coord);
            result->hi[j] = std::max(result->hi[j], coord);            
        }
    }
    return result;
}

DOPNode* buildDOP(std::vector<DOPNode*> &leaves, int first, int last)
{
    if (last - first == 1)
        return leaves[first];
    DOPNode* result = new DOPNode;
    result->face = -1;
    for (int i = 0; i < 13; i++)
    {
        result->lo[i] = std::numeric_limits<double>::infinity();
        result->hi[i] = -std::numeric_limits<double>::infinity();
    }
    for (int i = first; i < last; i++)
    {
        for (int j = 0; j < 13; j++)
        {
            result->lo[j] = std::min(result->lo[j], leaves[i]->lo[j]);
            result->hi[j] = std::max(result->hi[j], leaves[i]->hi[j]);
        }
    }
    double bestlen = 0;
    int bestaxis = 0;
    for (int j = 0; j < 13; j++)
    {
        double len = result->hi[j] - result->lo[j];
        if (len > bestlen)
        {
            bestlen = len;
            bestaxis = j;
        }
    }
    std::sort(leaves.begin() + first, leaves.begin() + last, 
        [bestaxis](DOPNode* a, DOPNode* b) -> bool
    {
        return a->lo[bestaxis] + a->hi[bestaxis] < b->lo[bestaxis] + b->hi[bestaxis];
    }
    );
    int mid = first + (last - first) / 2;
    result->left = buildDOP(leaves, first, mid);
    result->right = buildDOP(leaves, mid, last);
    return result;
}

static double triangleTriangleDist(const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1, int v1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2, int v2)
{    
    double result = std::numeric_limits<double>::infinity();

    Eigen::Vector3d p[3];
    Eigen::Vector3d q[3];
    for (int i = 0; i < 3; i++)
    {
        p[i] = V1.row(F1(v1, i)).transpose();
        p[i] = V2.row(F2(v2, i)).transpose();
    }
    return NoDivTriTriIsect(p[0].data(), p[1].data(), p[2].data(), q[0].data(), q[1].data(), q[2].data());
}

static bool boxesIntersect(DOPNode* tree1, DOPNode* tree2)
{
    for (int i = 0; i < 13; i++)
    {
        if (tree1->hi[i] < tree2->lo[i] || tree1->lo[i] > tree2->hi[i])
            return false;
    }
    return true;
}

static bool treeIntersection(DOPNode* tree1, DOPNode* tree2, const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2, const Eigen::MatrixXi& F2, int& face1, int& face2)
{
    if (!boxesIntersect(tree1, tree2))
        return false;

    if (tree1->left)
    {
        if (treeIntersection(tree1->left, tree2, V1, F1, V2, F2, face1, face2))
            return true;
        else if (treeIntersection(tree1->right, tree2, V1, F1, V2, F2, face1, face2))
            return true;
        else return false;
    }
    else if (tree2->left)
    {
        if (treeIntersection(tree1, tree2->left, V1, F1, V2, F2, face1, face2))
            return true;
        else if (treeIntersection(tree1, tree2->right, V1, F1, V2, F2, face1, face2))
            return true;
        else return false;
    }
    else
    {
        if (triangleTriangleDist(V1, F1, tree1->face, V2, F2, tree2->face) == 0.0)
        {
            face1 = tree1->face;
            face2 = tree2->face;
            return true;
        }
    }

    return false;
}

bool meshIntersection(
    const Eigen::MatrixXd& V1,
    const Eigen::MatrixXi& F1,
    const Eigen::MatrixXd& V2,
    const Eigen::MatrixXi& F2,
    int &face1,
    int &face2
)
{
    int nfaces1 = F1.rows();
    int nfaces2 = F2.rows();

    // quick reject

    double hi1[13];
    double lo1[13];
    double hi2[13];
    double lo2[13];
    for (int i = 0; i < 13; i++)
    {
        lo1[i] = std::numeric_limits<double>::infinity();
        lo2[i] = std::numeric_limits<double>::infinity();
        hi1[i] = -std::numeric_limits<double>::infinity();
        hi2[i] = -std::numeric_limits<double>::infinity();
    }

    int nverts1 = V1.rows();
    int nverts2 = V2.rows();
    for (int i = 0; i < nverts1; i++)
    {
        Eigen::Vector3d pt = V1.row(i).transpose();
        for (int j = 0; j < 13; j++)
        {
            Eigen::Vector3d DOPvec;
            for (int k = 0; k < 3; k++)
                DOPvec[k] = DOPNode::DOPvecs[3 * j + k];
            double coord = pt.dot(DOPvec) / DOPvec.norm();
            lo1[j] = std::min(lo1[j], coord);
            hi1[j] = std::max(hi1[j], coord);            
        }
    }
    for (int i = 0; i < nverts2; i++)
    {
        Eigen::Vector3d pt = V2.row(i).transpose();
        for (int j = 0; j < 13; j++)
        {
            Eigen::Vector3d DOPvec;
            for (int k = 0; k < 3; k++)
                DOPvec[k] = DOPNode::DOPvecs[3 * j + k];
            double coord = pt.dot(DOPvec) / DOPvec.norm();
            lo2[j] = std::min(lo2[j], coord);
            hi2[j] = std::max(hi2[j], coord);            
        }
    }

    for (int i = 0; i < 13; i++)
    {
        if (hi1[i] < lo2[i] || hi2[i] < lo1[i])
            return false;
    }
    
    std::vector<DOPNode*> leaves1(nfaces1);
    for (int i = 0; i < nfaces1; i++)
    {
        leaves1[i] = fitFace(V1, F1, i);
    }
    DOPNode* tree1 = buildDOP(leaves1, 0, leaves1.size());

    std::vector<DOPNode*> leaves2(nfaces2);
    for (int i = 0; i < nfaces2; i++)
    {
        leaves2[i] = fitFace(V2, F2, i);
    }
    DOPNode* tree2 = buildDOP(leaves2, 0, leaves2.size());
    
    bool result = treeIntersection(tree1, tree2, V1, F1, V2, F2, face1, face2);

    delete tree1;
    delete tree2;
    return result;
}