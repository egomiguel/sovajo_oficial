#ifndef BINARY_TREE_H
#define BINARY_TREE_H

#include "vtkSmartPointer.h"
#include "vtkPoints.h"
#include "vtkTriangle.h"
#include "vtkPolyData.h"

struct BinaryTree
{
    BinaryTree(BinaryTree * father, BinaryTree * right, BinaryTree * left, int downPoint, int upPoint, double value, bool finalNode);
    BinaryTree * right;
    BinaryTree * left;
    BinaryTree * father;
    int downPoint, upPoint;
    double value;
    bool isLeaf;
    bool finalNode;
    int extremeAverage;
    int downAmount;
    int upAmount;

    int GetBalance();
    BinaryTree * GenerateRight(int maxUp, int maxDown, const vtkSmartPointer<vtkPoints> upPoints, const vtkSmartPointer<vtkPoints> downPoints, bool closed);
    BinaryTree * GenerateLeft(int maxUp, int maxDown, const vtkSmartPointer<vtkPoints> upPoints, const vtkSmartPointer<vtkPoints> downPoints, bool closed);
    
    ~BinaryTree();
};

class GenerateTree
{
private:
    vtkSmartPointer<vtkPoints> downPoints, upPoints;
    bool closed;
    BinaryTree * myTree;
    BinaryTree * finalNode;
    int maxDown, maxUp;
    std::vector<BinaryTree *> nodes;
    void GeneratePath(std::vector<BinaryTree *> childs);
    void GetPath(const BinaryTree * pFinalNode, std::vector<std::pair<int, int>>& path);
public:
    ~GenerateTree();
    GenerateTree(const vtkSmartPointer<vtkPoints> downPoints, const vtkSmartPointer<vtkPoints> upPoints, bool closed = true);
    std::vector<std::pair<int, int>> GetPath();
};


#endif