#include <opencv2/calib3d/calib3d.hpp>
#include "BinaryTree.hpp"
#include "vtkNew.h"
#include "SPlane.hpp"

using namespace TKA::SEGMENTATION;

double CalculateValue(const double P1[3], const double P2[3], const double P3[3])
{
    cv::Point3d point1 = { P1[0], P1[1], P1[2] };
    cv::Point3d point2 = { P2[0], P2[1], P2[2] };
    cv::Point3d point3 = { P3[0], P3[1], P3[2] };
   
    cv::Point3d v1 = point1 - point3;
    cv::Point3d v2 = point2 - point3;
    cv::Point3d crossP = v1.cross(v2);

    return (sqrt(crossP.dot(crossP)));
}

void showNode(BinaryTree * node, std::string msg)
{

    if (node)
    {
        msg = msg + ": (";
        std::cout << msg << node->downPoint << "; " << node->upPoint << " ).  Value: " << node->value << " average Exteme: " << node->extremeAverage << std::endl;
    }
    else
    {
        msg = msg + ". Empty.";
        std::cout << msg << std::endl;
    }
}

void showPoint(double a[3])
{
    std::cout << a[0] << "; " << a[1] << "; " << a[2] << std::endl;
}


BinaryTree::BinaryTree(BinaryTree * father, BinaryTree * right, BinaryTree * left, int downPoint, int upPoint, double value, bool finalNode)
{
    this->father = father;
    this->right = right;
    this->left = left;
    this->downPoint = downPoint;
    this->upPoint = upPoint;
    this->value = value;
    this->finalNode = finalNode;
    this->downAmount = 0;
    this->upAmount = 0;
    this->isLeaf = false;
    if (father)
    {
        this->extremeAverage = father->extremeAverage;
        if (downPoint == father->downPoint)
        {
            downAmount = father->downAmount + 1;
        }

        if (upPoint == father->upPoint)
        {
            upAmount = father->upAmount + 1;
        }
    }
    else
    {
        this->extremeAverage = 0;
    }

}

int BinaryTree::GetBalance()
{
    return abs(downAmount - upAmount);
}

BinaryTree * BinaryTree::GenerateRight(int maxUp, int maxDown, const vtkSmartPointer<vtkPoints> upPoints, const vtkSmartPointer<vtkPoints> downPoints, bool closed)
{
    BinaryTree * result = nullptr;
    if (closed == true)
    {
        if (upPoint >= maxUp)
        {
            return result;
        }
    }
    else
    {
        if (upPoint >= maxUp - 1)
        {
            return result;
        }
    }

    double a[3];
    double b[3];
    double c[3];

    int p1, p2, p3;
    p1 = upPoint;
    p2 = upPoint + 1;
    p3 = downPoint;

    if (p2 == maxUp)
    {
        p2 = 0;
    }
    if (p3 == maxDown)
    {
        p3 = 0;
    }

    upPoints->GetPoint(p1, a);
    upPoints->GetPoint(p2, b);
    downPoints->GetPoint(p3, c);

    double newValue = CalculateValue(a, b, c) + this->value;

    bool tFinalNode;
    if (closed == true)
    {
        if (downPoint >= maxDown && upPoint + 1 >= maxUp)
        {
            tFinalNode = true;
        }
        else
        {
            tFinalNode = false;
        }
    }
    else
    {
        if (downPoint >= maxDown - 1 && upPoint + 1 >= maxUp - 1)
        {
            tFinalNode = true;
        }
        else
        {
            tFinalNode = false;
        }
    }

    result = new BinaryTree(this, nullptr, nullptr, downPoint, upPoint + 1, newValue, tFinalNode);

    if (closed == true)
    {
        if (result->downPoint == maxDown || result->upPoint == maxUp || result->downPoint == 0 || result->upPoint == 0)
        {
            result->extremeAverage = 1 + this->extremeAverage;
        }
    }
    else
    {
        if (result->downPoint == maxDown - 1 || result->upPoint == maxUp - 1 || result->downPoint == 0 || result->upPoint == 0)
        {
            result->extremeAverage = 1 + this->extremeAverage;
        }
    }

    return result;
}

BinaryTree * BinaryTree::GenerateLeft(int maxUp, int maxDown, const vtkSmartPointer<vtkPoints> upPoints, const vtkSmartPointer<vtkPoints> downPoints, bool closed)
{
    BinaryTree * result = nullptr;
    if (closed == true)
    {
        if (downPoint >= maxDown)
        {
            return result;
        }
    }
    else
    {
        if (downPoint >= maxDown - 1)
        {
            return result;
        }
    }

    double a[3];
    double b[3];
    double c[3];

    int p1, p2, p3;
    p1 = downPoint;
    p2 = downPoint + 1;
    p3 = upPoint;

    if (p2 == maxDown)
    {
        p2 = 0;
    }
    if (p3 == maxUp)
    {
        p3 = 0;
    }

    upPoints->GetPoint(p3, a);
    downPoints->GetPoint(p1, b);
    downPoints->GetPoint(p2, c);

    double newValue = CalculateValue(b, c, a) + this->value;

    bool tFinalNode;

    if (closed == true)
    {
        if (downPoint + 1 >= maxDown && upPoint >= maxUp)
        {
            tFinalNode = true;
        }
        else
        {
            tFinalNode = false;
        }
    }
    else
    {
        if (downPoint + 1 >= maxDown - 1 && upPoint >= maxUp - 1)
        {
            tFinalNode = true;
        }
        else
        {
            tFinalNode = false;
        }
    }

    result = new BinaryTree(this, nullptr, nullptr, downPoint + 1, upPoint, newValue, tFinalNode);

    if (closed == true)
    {
        if (result->downPoint == maxDown || result->upPoint == maxUp || result->downPoint == 0 || result->upPoint == 0)
        {
            result->extremeAverage = 1 + this->extremeAverage;
        }
    }
    else
    {
        if (result->downPoint == maxDown - 1 || result->upPoint == maxUp - 1 || result->downPoint == 0 || result->upPoint == 0)
        {
            result->extremeAverage = 1 + this->extremeAverage;
        }
    }

    return result;
}

BinaryTree::~BinaryTree()
{
    right = NULL;
    left = NULL;
    father = NULL;
    delete right;
    delete left;
    delete father;
}

GenerateTree::GenerateTree(const vtkSmartPointer<vtkPoints> downPoints, const vtkSmartPointer<vtkPoints> upPoints, bool closed)
{
    this->downPoints = downPoints;
    this->upPoints = upPoints;
    this->closed = closed;
    maxDown = downPoints->GetNumberOfPoints();
    maxUp = upPoints->GetNumberOfPoints();
    myTree = new BinaryTree(nullptr, nullptr, nullptr, 0, 0, 0, false);
    finalNode = new BinaryTree(nullptr, nullptr, nullptr, 0, 0, 0, false);
    BinaryTree * tRight = myTree->GenerateRight(maxUp, maxDown, upPoints, downPoints, closed);
    BinaryTree * tLeft = myTree->GenerateLeft(maxUp, maxDown, upPoints, downPoints, closed);
    myTree->right = tRight;
    myTree->left = tLeft;
    nodes.push_back(tRight);
    nodes.push_back(tLeft);
}

GenerateTree::~GenerateTree()
{
    myTree = NULL;
    finalNode = NULL;
    delete myTree;
    delete finalNode;
    for (auto& it : nodes)
    {
        delete it;
    }
}

void GenerateTree::GetPath(const BinaryTree * pFinalNode, std::vector<std::pair<int, int>>& path)
{
    if (pFinalNode)
    {
        int a = pFinalNode->downPoint;
        int b = pFinalNode->upPoint;
        if (a == maxDown)
        {
            a = 0;
        }
        if (b == maxUp)
        {
            b = 0;
        }
        path.push_back(std::make_pair(a, b));
        GetPath(pFinalNode->father, path);
    }
    else
    {
        return;
    }
}

std::vector<std::pair<int, int>> GenerateTree::GetPath()
{
    GeneratePath(nodes);

    std::vector<std::pair<int, int>> path;
    GetPath(finalNode, path);
    return path;
}

void GenerateTree::GeneratePath(std::vector<BinaryTree *> childs)
{
    rsize_t size = childs.size();
    if (size == 0)
    {
        return;
    }
    std::vector<BinaryTree *> newChilds(size + 1, nullptr);

    std::vector<BinaryTree *>::iterator it1 = childs.begin();
    std::vector<BinaryTree *>::iterator it2 = childs.end();
    bool finish = true;
    int cont = 0;

    for (; it1 != it2; ++it1)
    {
        if ((*it1) == nullptr)
        {
            continue;
        }
        if ((*it1)->finalNode == true)
        {
            if ((*it1)->value > finalNode->value || finalNode->value == 0)
            {
                finalNode = (*it1);
            }
            continue;
        }
        if ((*it1)->isLeaf == true)
        {
            continue;
        }
        cont++;
        finish = false;
        BinaryTree * tRight = (*it1)->GenerateRight(maxUp, maxDown, upPoints, downPoints, closed);
        BinaryTree * tLeft = (*it1)->GenerateLeft(maxUp, maxDown, upPoints, downPoints, closed);

        (*it1)->right = tRight;
        (*it1)->left = tLeft;

        if (tRight)
        {
            if (newChilds[tRight->downPoint] == nullptr)
            {
                newChilds[tRight->downPoint] = tRight;
            }
            else
            {
                if (tRight->value < newChilds[tRight->downPoint]->value)
                {
                    newChilds[tRight->downPoint]->isLeaf = true;
                    newChilds[tRight->downPoint] = tRight;
                }
                else if (tRight->value == newChilds[tRight->downPoint]->value)
                {
                    if (tRight->extremeAverage >= newChilds[tRight->downPoint]->extremeAverage || tRight->GetBalance() > newChilds[tRight->downPoint]->GetBalance())
                    {
                        tRight->isLeaf = true;
                    }
                    else
                    {
                        newChilds[tRight->downPoint]->isLeaf = true;
                        newChilds[tRight->downPoint] = tRight;
                    }
                }

            }
        }

        if (tLeft)
        {
            if (newChilds[tLeft->downPoint] == nullptr)
            {
                newChilds[tLeft->downPoint] = tLeft;
            }
            else
            {
                if (tLeft->value < newChilds[tLeft->downPoint]->value)
                {
                    newChilds[tLeft->downPoint]->isLeaf = true;
                    newChilds[tLeft->downPoint] = tLeft;
                }
                else if (tLeft->value == newChilds[tLeft->downPoint]->value)
                {
                    if (tLeft->extremeAverage >= newChilds[tLeft->downPoint]->extremeAverage || tLeft->GetBalance() > newChilds[tLeft->downPoint]->GetBalance())
                    {
                        tLeft->isLeaf = true;
                    }
                    else
                    {
                        newChilds[tLeft->downPoint]->isLeaf = true;
                        newChilds[tLeft->downPoint] = tLeft;
                    }
                }
            }
        }
    }

    if (finish == true)
    {
        return;
    }

    GeneratePath(newChilds);
}

