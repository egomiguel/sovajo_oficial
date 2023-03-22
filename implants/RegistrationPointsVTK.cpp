//#include "RegistrationPointsVTK.hpp"
//#include "ImplantTools.hpp"
//#include "vtkPlane.h"
//#include "vtkPlaneCollection.h"
//#include "vtkClipClosedSurface.h"
//#include "opencv2/imgproc.hpp"
//#include "vtkImplicitPolyDataDistance.h"
//
//
//RegistrationPointsVTK::RegistrationPointsVTK(const Knee& pKnee)
//{
//    mKnee = pKnee;
//
//    makeAxis();
//
//    mGoodLateralRef = false;
//    mGoodMediaRef = false;
//
//    Point vectorTEA = -mVectorTibiaTEALat;
//
//    Plane coronal;
//
//    mTibiaTransverse.reverseByPoint(mKnee.getAnkleCenter());
//    
//    coronal.init(mVectorTibiaAPFront, mKnee.getMedialPlateau());
//    coronal.reverseByPoint(mKnee.getTibiaTubercle());
//
//    myTuber = ImplantTools::GetFarestPoint(mKnee.GetTibiaPoly(), coronal, mTibiaTransverse);
//    Line helpLine(mVectorTibiaAPFront, (myTuber - vectorTEA));
//    Point myTuberIn = helpLine.getProjectPoint(mCenterEllipse);
//
//    double distanceRef = abs(mTibiaTransverse.eval(myTuberIn));
//
//    mTransverseRef = myTuber + 0.95 * distanceRef * mVectorTibiaAxisUp;
//    mObliqueRef = myTuber + 0.5 * distanceRef * mVectorTibiaAxisUp;
//    mTuberUpRef = myTuber + 0.5 * distanceRef * mVectorTibiaAxisUp;
//    mTuberDownRef = myTuber - 0.1 * distanceRef * mVectorTibiaAxisUp;
//    mTuberSideRef = myTuber + (0.6 * distanceRef) * vectorTEA;
//
//    Point a1 = myTuberIn + 0.95 * distanceRef * mVectorTibiaAxisUp;
//    Point b1 = myTuberIn + 0.5 * distanceRef * mVectorTibiaAxisUp;
//    Point c1 = myTuberIn + 0.5 * distanceRef * mVectorTibiaAxisUp;
//    Point d1 = myTuberIn - 0.1 * distanceRef * mVectorTibiaAxisUp;
//    Point e1 = myTuberIn + (0.6 * distanceRef) * vectorTEA;
//
//    Point a2 = a1 + 100.0 * mVectorTibiaAPFront;
//    Point b2 = b1 + 100.0 * mVectorTibiaAPFront;
//    Point c2 = c1 + 100.0 * mVectorTibiaAPFront;
//    Point d2 = d1 + 100.0 * mVectorTibiaAPFront;
//    Point e2 = e1 + 100.0 * mVectorTibiaAPFront;
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), a1, a2, mTransverseRef);
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), b1, b2, mObliqueRef);
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), c1, c2, mTuberUpRef);
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), d1, d2, mTuberDownRef);
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), e1, e2, mTuberSideRef);
//
//    Plane sagital;
//    sagital.init(mVectorFemurTEALat, mKnee.getFemurKneeCenter());
//    sagital.reverseByPoint(mKnee.getMedialEpicondyle());
//
//    Point planeNormal, planePoint;
//    vtkNew<vtkPlane> vtkPlaneMedial, vtkPlanelateral;
//
//    planeNormal = sagital.getNormalVector();
//    planePoint = sagital.getPoint();
//
//    vtkPlaneMedial->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
//    vtkPlaneMedial->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);
//
//    vtkNew<vtkPlaneCollection> medialPlanes, lateralPlanes;
//    medialPlanes->AddItem(vtkPlaneMedial);
//
//    vtkNew<vtkClipClosedSurface> medialClipper, lateralClipper;
//    medialClipper->SetInputData(mKnee.GetFemurPoly());
//    medialClipper->SetClippingPlanes(medialPlanes);
//    medialClipper->Update();
//
//    sagital.reverse();
//    planeNormal = sagital.getNormalVector();
//    planePoint = sagital.getPoint();
//
//    vtkPlanelateral->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
//    vtkPlanelateral->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);
//
//    lateralPlanes->AddItem(vtkPlanelateral);
//
//    lateralClipper->SetInputData(mKnee.GetFemurPoly());
//    lateralClipper->SetClippingPlanes(lateralPlanes);
//    lateralClipper->Update();
//
//    medialSideFemur = medialClipper->GetOutput();
//    lateralSideFemur = lateralClipper->GetOutput();
//}
//
//void RegistrationPointsVTK::makeAxis()
//{
//    mVectorFemurAxisUp = mKnee.getDirectVectorFemurAxis();
//    mVectorFemurAPFront = mKnee.getFemurDirectVectorAP();
//    mVectorFemurTEALat = mKnee.getFemurVectorLateralTEA();
//
//    std::pair<double, double> ellipseSize = mKnee.getTibiaAutomaticAxis(mVectorTibiaAPFront, mVectorTibiaTEALat, mCenterEllipse);
//
//    mVectorTibiaAxisUp = mVectorTibiaAPFront.cross(mVectorTibiaTEALat);
//
//    if (mKnee.getIsRight() == true)
//    {
//        mVectorTibiaAxisUp = (-1.0) * mVectorTibiaAxisUp;
//    }
//    mVectorTibiaAxisUp.normalice();
//
//    Point tempCenterLat = mCenterEllipse + (ellipseSize.second / 4.0) * mVectorTibiaTEALat;
//    Point tempCenterMed = mCenterEllipse - (ellipseSize.second / 4.0) * mVectorTibiaTEALat;
//
//    Point a = tempCenterLat + 100.0 * mVectorTibiaAxisUp;
//    Point b = tempCenterLat - 10.0 * mVectorTibiaAxisUp;
//
//    bool result1 = ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), a, b, mBigCenterLat);
//
//    if (result1 == false)
//    {
//        Line myLine = Line::makeLineWithPoints(a, b);
//        mBigCenterLat = myLine.getProjectPoint(mKnee.getLateralPlateau());
//    }
//
//    a = tempCenterMed + 100.0 * mVectorTibiaAxisUp;
//    b = tempCenterMed - 10.0 * mVectorTibiaAxisUp;
//
//    bool result2 = ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), a, b, mBigCenterMed);
//
//    if (result2 == false)
//    {
//        Line myLine = Line::makeLineWithPoints(a, b);
//        mBigCenterLat = myLine.getProjectPoint(mKnee.getMedialPlateau());
//    }
//    double transverseMargin = 4.0;
//    a = mBigCenterLat - transverseMargin * mVectorTibiaAxisUp;
//    b = mBigCenterMed - transverseMargin * mVectorTibiaAxisUp;
//    Point c = (mKnee.getTibiaKneeCenter() + 10.0 * mVectorTibiaAPFront) - (transverseMargin - 1.0) * mVectorTibiaAxisUp;
//    mTibiaTransverse.init(a, b, c);
//}
//
//RegistrationPointsVTK::Sphere::Sphere(const Point& pCenter, double pRadius)
//{
//    center = pCenter;
//    radius = pRadius;
//    real = true;
//}
//
//RegistrationPointsVTK::Sphere::Sphere()
//{
//    real = false;
//}
//
//bool RegistrationPointsVTK::Sphere::isPointInside(const Point& pnt) const
//{
//    if (real == false)
//    {
//        return false;
//    }
//
//    Point diff = pnt - center;
//    double temp = sqrt(diff.dot(diff));
//    if (temp <= radius)
//    {
//        return true;
//    }
//    return false;
//}
//
//Point RegistrationPointsVTK::getNearPointToLine(const Line& pLine, const vtkSmartPointer<vtkPoints> Ppoints, const Sphere& pExclude)
//{
//    int pointSize = Ppoints->GetNumberOfPoints();
//    double dist = -1;
//    Point nearPoint;
//    for (int j = 0; j < pointSize; j++)
//    {
//        double pnt[3];
//        Ppoints->GetPoint(j, pnt);
//        Point pntCV(pnt[0], pnt[1], pnt[2]);
//
//        if (pExclude.isPointInside(pntCV) == false)
//        {
//            double distTemp = pLine.getDistanceFromPoint(pntCV);
//
//            if (distTemp < dist || dist < 0)
//            {
//                dist = distTemp;
//                nearPoint = pntCV;
//            }
//        }
//    }
//    return nearPoint;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetLateralMainPointsFemur(int amount, int pTest)
//{
//    std::vector<Point> result;
//
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//    Point vectorTEA = mKnee.getLateralEpicondyle() - mKnee.getMedialEpicondylePerp();
//    Point vectorAP = vectorAxis.cross(vectorTEA);
//
//    Point a = (mKnee.getTopPointOnPatellaPath() + mKnee.getLateralCoronalDistalPoint()) / 2.0;
//    Point b = mKnee.getLateralCondyle();
//
//    if (pTest == 1)
//    {
//        a = mKnee.getLateralCoronalDistalPoint() + 0.3 * (mKnee.getTopPointOnPatellaPath() - mKnee.getLateralCoronalDistalPoint());
//    }
//
//    Plane transverse;
//    transverse.init(vectorAxis, mKnee.getFemurKneeCenter());
//    Plane onPlane = transverse.getPerpendicularPlane(a, b);
//
//    vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(mKnee.GetFemurPoly(), onPlane.getNormalVector(), onPlane.getPoint());
//
//    vtkIdType nearA, nearB;
//
//    nearA = ImplantTools::GetNearestPoints(mainContour, a);
//    nearB = ImplantTools::GetNearestPoints(mainContour, b);
//
//    Plane transverseCondyle;
//    transverseCondyle.init(vectorAxis, b);
//
//    if (transverseCondyle.eval(mKnee.getHipCenter()) > 0)
//    {
//        transverseCondyle.reverse();
//    }
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(mainContour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearB)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearB && pos == -1)
//    {
//        pos = tSize - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return result;
//    }
//
//    std::vector<Point> tempList;
//
//    if (transverseCondyle.eval(mainContour->GetPoint(allPoints[5])) < 0)
//    {
//        for (int i = tSize - 1; i > 0; i--)
//        {
//            if (allPoints[i] == nearA)
//            {
//                break;
//            }
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//        }
//    }
//    else
//    {
//        for (int i = 1; i < tSize; i++)
//        {
//            if (allPoints[i] == nearA)
//            {
//                break;
//            }
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//        }
//    }
//
//    if (amount <= 0)
//    {
//        return tempList;
//    }
//
//    std::reverse(tempList.begin(), tempList.end());
//
//    result = reduceSortPointsByAmount(tempList, amount + 2);
//
//    if (result.size() > amount)
//    {
//        return std::vector<Point>(result.begin(), result.begin() + amount);
//    }
//    else
//    {
//        return result;
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::GetMedialMainPointsFemur(int amount, int pTest)
//{
//    std::vector<Point> result;
//
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//    Point vectorTEA = mKnee.getLateralEpicondyle() - mKnee.getMedialEpicondylePerp();
//    Point vectorAP = vectorAxis.cross(vectorTEA);
//
//    Plane transverseCondyle;
//    
//    Point a = (mKnee.getTopPointOnPatellaPath() + mKnee.getMedialCoronalDistalPoint()) / 2.0;
//    Point b = mKnee.getMedialCondyle();
//
//    if (pTest == 1)
//    {
//        a = mKnee.getMedialInferiorFemurPoint();
//    }
//
//    if (pTest == 2)
//    {
//        a = mKnee.getMedialCoronalDistalPoint() + 0.2 * (mKnee.getTopPointOnPatellaPath() - mKnee.getMedialCoronalDistalPoint());
//        b = mKnee.getMedialInferiorFemurPoint();
//    }
//
//    if (pTest == 2)
//    {
//        transverseCondyle.init(vectorAP, b);
//
//        if (transverseCondyle.eval(mKnee.getMedialCondyle()) > 0)
//        {
//            transverseCondyle.reverse();
//        }
//    }
//    else
//    {
//        transverseCondyle.init(vectorAxis, b);
//
//        if (transverseCondyle.eval(mKnee.getHipCenter()) > 0)
//        {
//            transverseCondyle.reverse();
//        }
//    }
//
//    Plane transverse;
//    transverse.init(vectorAxis, mKnee.getFemurKneeCenter());
//    Plane onPlane = transverse.getPerpendicularPlane(a, b);
//
//    vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(mKnee.GetFemurPoly(), onPlane.getNormalVector(), onPlane.getPoint());
//
//    vtkIdType nearA, nearB;
//
//    nearA = ImplantTools::GetNearestPoints(mainContour, a);
//    nearB = ImplantTools::GetNearestPoints(mainContour, b);
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(mainContour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearB)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearB && pos == -1)
//    {
//        pos = tSize - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return result;
//    }
//
//    std::vector<Point> tempList;
//
//    if (transverseCondyle.eval(mainContour->GetPoint(allPoints[5])) < 0)
//    {
//        for (int i = tSize - 1; i > 0; i--)
//        {
//            if (allPoints[i] == nearA)
//            {
//                break;
//            }
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//        }
//    }
//    else
//    {
//        for (int i = 1; i < tSize; i++)
//        {
//            if (allPoints[i] == nearA)
//            {
//                break;
//            }
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//        }
//    }
//
//    if (amount <= 0)
//    {
//        return tempList;
//    }
//
//    std::reverse(tempList.begin(), tempList.end());
//
//    result = reduceSortPointsByAmount(tempList, amount + 2);
//
//    if (result.size() > amount)
//    {
//        return std::vector<Point>(result.begin(), result.begin() + amount);
//    }
//    else
//    {
//        return result;
//    }
//}
//
//void RegistrationPointsVTK::reduceSortPointsByRadius(const std::vector<Point>& sortPoints, double radius, std::vector<Point>& result)
//{
//    int tSize = sortPoints.size();
//
//    if (tSize < 3)
//    {
//        return;
//    }
//
//    result.push_back(sortPoints[0]);
//
//    for (int i = 1; i < tSize; i++)
//    {
//        if (ImplantTools::isPointInsideSphere(result[result.size() - 1], radius, sortPoints[i]) == false)
//        {
//            Line myLine = Line::makeLineWithPoints(sortPoints[i], sortPoints[i - 1]);
//            std::pair<Point, Point> intercep;
//
//            myLine.getInterceptionSphere(result[result.size() - 1], radius, intercep);
//
//            if (ImplantTools::getDistanceBetweenPoints(sortPoints[i], intercep.first, true) < ImplantTools::getDistanceBetweenPoints(sortPoints[i], intercep.second, true))
//            {
//                result.push_back(intercep.first);
//            }
//            else
//            {
//                result.push_back(intercep.second);
//            }
//        }
//
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::reduceSortPointsByRadius(const std::vector<Point>& sortPoints, double radius)
//{
//    std::vector<Point> result;
//
//    int tSize = sortPoints.size();
//
//    if (tSize < 3)
//    {
//        return sortPoints;
//    }
//
//    result.push_back(sortPoints[0]);
//
//    for (int i = 1; i < tSize; i++)
//    {
//        if (ImplantTools::isPointInsideSphere(result[result.size() - 1], radius, sortPoints[i]) == false)
//        {
//            Line myLine = Line::makeLineWithPoints(sortPoints[i], sortPoints[i - 1]);
//            std::pair<Point, Point> intercep;
//
//            myLine.getInterceptionSphere(result[result.size() - 1], radius, intercep);
//
//            if (ImplantTools::getDistanceBetweenPoints(sortPoints[i], intercep.first, true) < ImplantTools::getDistanceBetweenPoints(sortPoints[i], intercep.second, true))
//            {
//                result.push_back(intercep.first);
//            }
//            else
//            {
//                result.push_back(intercep.second);
//            }
//        }
//
//    }
//
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::reduceSortPointsByAmount(const std::vector<Point>& sortPoints, int pointsAmount)
//{
//    int tSize = sortPoints.size();
//
//    if (tSize < 3)
//    {
//        return sortPoints;
//    }
//
//    double dist = ImplantTools::getDistanceBetweenPoints(sortPoints[0], sortPoints[tSize - 1]);
//    double radius = (dist / pointsAmount);
//
//    std::vector<Point> result;
//    result = reduceSortPointsByRadius(sortPoints, radius);
//
//    while (result.size() > pointsAmount)
//    {
//        radius++;
//        result.clear();
//        result = reduceSortPointsByRadius(sortPoints, radius);
//    }
//
//    if (result.size() == pointsAmount)
//    {
//        return result;
//    }
//    else
//    {
//        radius--;
//        return reduceSortPointsByRadius(sortPoints, radius);
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::GetLateralBorderPointsFemur(int amount)
//{
//    std::vector<Point> result, resultFinal;
//    std::vector<Point> mainPoints = GetLateralMainPointsFemur(amount);
//   
//    Point lateralEpi = mKnee.getLateralEpicondyle();
//    Point vectorTEA = mKnee.getMedialEpicondylePerp() - mKnee.getLateralEpicondyle();
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//
//    int tSize = mainPoints.size();
//
//    if (tSize == 0)
//    {
//        return resultFinal;
//    }
//
//    double radius = (ImplantTools::getDistanceBetweenPoints(lateralEpi, mainPoints[tSize - 1])) * 0.6;
//
//    Plane sagital;
//    sagital.init(vectorTEA, mKnee.getFemurKneeCenter());
//
//    Sphere mySphere(lateralEpi, radius);
//
//    std::vector<vtkSmartPointer<vtkPolyData>> contoursList;
//
//    for (int i = 0; i < tSize; i++)
//    {
//        Plane newPlane = sagital.getPerpendicularPlane(mainPoints[i], lateralEpi);
//        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(lateralSideFemur, newPlane.getNormalVector(), newPlane.getPoint());
//        vtkSmartPointer<vtkPoints> tPoints = contour->GetPoints();
//        contoursList.push_back(contour);
//        int pointSize = tPoints->GetNumberOfPoints();
//
//        if (pointSize < 20)
//        {
//            return resultFinal;
//        }
//
//        Line line1(sagital.getNormalVector(), mainPoints[i]);
//        Plane onEpicondyle;
//        onEpicondyle.init(sagital.getNormalVector(), lateralEpi);
//
//        Point inter = onEpicondyle.getInterceptionLinePoint(line1);
//
//        Point vector1 = inter - mainPoints[i];
//        Point vector2 = inter - lateralEpi;
//        vector1.normalice();
//        vector2.normalice();
//
//        Point epiPoint, basePoint;
//        epiPoint = inter + 1000.0 * vector1;
//        basePoint = inter + 1000.0 * vector2;
//
//        Point topBase = getNearPointToLine(Line(vector1, basePoint), tPoints, mySphere);
//        Point topEpi = getNearPointToLine(Line(vector2, epiPoint), tPoints, mySphere);
//
//        onEpicondyle.movePlane(topEpi);
//        inter = onEpicondyle.getInterceptionLinePoint(Line(vector1, topBase));
//
//        vector1 = inter - topBase;
//        vector2 = inter - topEpi;
//        vector1.normalice();
//        vector2.normalice();
//
//        Point a, b;
//        a = inter + 10.0 * vector1;
//        b = inter + 10.0 * vector2;
//
//        Line obliqueLine = Line::makeLineWithPoints(a, b);
//
//        Point nearPoint = getNearPointToLine(obliqueLine, tPoints, mySphere);
//
//        result.push_back(nearPoint);
//    }
//
//    double radius2 = 5.0;
//
//    for (int i = 0; i < tSize; i++)
//    {
//        vtkIdType nearEpi;
//        nearEpi = ImplantTools::GetNearestPoints(contoursList[i], lateralEpi);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(contoursList[i], lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//        it1 = lines.begin();
//        it2 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//
//        int pos = 0;
//        int cont = 0;
//
//        for (; it1 != it2; ++it1)
//        {
//            allPoints.push_back(it1->first);
//            if (it1->first == nearEpi)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        if (lines.back().second == nearEpi && pos == -1)
//        {
//            pos = cont - 1;
//        }
//
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//
//        Point decisionVector = lateralEpi - result[i];
//        Plane decisionPlane;
//        decisionPlane.init(decisionVector, lateralEpi);
//        decisionPlane.reverseByPoint(result[i]);
//
//        double decisionPnt[3];
//        contoursList[i]->GetPoint(allPoints[5], decisionPnt);
//
//        if (decisionPlane.eval(decisionPnt) < 0)
//        {
//            std::reverse(allPoints.begin(), allPoints.end());
//        }
//
//        Sphere mySphereBorder(result[i], radius2);
//
//        for (int j = 0; j < allPoints.size() - 1; j++)
//        {
//            double pnt1[3];
//            contoursList[i]->GetPoint(allPoints[j], pnt1);
//            Point checkPoint1(pnt1[0], pnt1[1], pnt1[2]);
//
//            double pnt2[3];
//            contoursList[i]->GetPoint(allPoints[j + 1], pnt2);
//            Point checkPoint2(pnt2[0], pnt2[1], pnt2[2]);
//
//            if ((mySphereBorder.isPointInside(checkPoint1) == false && mySphereBorder.isPointInside(checkPoint2) == true) ||
//                (mySphereBorder.isPointInside(checkPoint1) == true && mySphereBorder.isPointInside(checkPoint2) == false))
//            {
//                Line myLine = Line::makeLineWithPoints(checkPoint1, checkPoint2);
//                std::pair<Point, Point> intercep;
//
//                myLine.getInterceptionSphere(mySphereBorder.center, mySphereBorder.radius, intercep);
//
//                Point outSise;
//                if (mySphereBorder.isPointInside(checkPoint1) == false)
//                {
//                    outSise = checkPoint1;
//                }
//                else
//                {
//                    outSise = checkPoint2;
//                }
//
//                if (ImplantTools::getDistanceBetweenPoints(outSise, intercep.first, true) < ImplantTools::getDistanceBetweenPoints(outSise, intercep.second, true))
//                {
//                    resultFinal.push_back(intercep.first);
//                }
//                else
//                {
//                    resultFinal.push_back(intercep.second);
//                }
//            }
//        }
//
//    }
//    //std::cout << rightPos << std::endl;
//    return resultFinal;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetLateralObliquePointsFemur(int amount)
//{
//    std::vector<Point> result;
//    Point vectorTEA = mKnee.getMedialEpicondylePerp() - mKnee.getLateralEpicondyle();
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//    Point vectorAP = vectorTEA.cross(vectorAxis);
//    Plane coronal, sagital;
//    coronal.init(vectorAP, mKnee.getFemurKneeCenter());
//    sagital.init(vectorTEA, mKnee.getFemurKneeCenter());
//
//    if (sagital.eval(mKnee.getLateralEpicondyle()) < 0)
//    {
//        sagital.reverse();
//    }
//
//    Plane obliquePlane = coronal.getPerpendicularPlane(lastMedialBorderFemur, mKnee.getLateralCoronalDistalPoint());
//    vtkSmartPointer<vtkPolyData> obliqueContour = ImplantTools::getMaxContour(mKnee.GetFemurPoly(), obliquePlane.getNormalVector(), obliquePlane.getPoint());
//
//    vtkIdType nearA, nearB;
//
//    nearA = ImplantTools::GetNearestPoints(obliqueContour, mKnee.getLateralCoronalDistalPoint());
//    nearB = ImplantTools::GetNearestPoints(obliqueContour, lastMedialBorderFemur);
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(obliqueContour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearA)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearA && pos == -1)
//    {
//        pos = cont - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return result;
//    }
//
//    std::vector<Point> tempList1, tempList2;
//
//    for (int i = tSize - 1; i > 0; i--)
//    {
//        double pnt[3];
//        obliqueContour->GetPoint(allPoints[i], pnt);
//
//        if (sagital.eval(pnt) >= 0)
//        {
//            tempList1.push_back(Point(pnt[0], pnt[1], pnt[2]));
//        }
//        else
//        {
//            break;
//        }
//    }
//
//    for (int i = 0; i < tSize; i++)
//    {
//        double pnt[3];
//        obliqueContour->GetPoint(allPoints[i], pnt);
//
//        if (sagital.eval(pnt) >= 0)
//        {
//            tempList2.push_back(Point(pnt[0], pnt[1], pnt[2]));
//        }
//        else
//        {
//            break;
//        }
//    }
//
//    if (tempList1.size() < tempList2.size())
//    {
//        std::reverse(tempList1.begin(), tempList1.end());
//        result = reduceSortPointsByAmount(tempList1, amount + 1);
//    }
//    else
//    {
//        std::reverse(tempList2.begin(), tempList2.end());
//        result = reduceSortPointsByAmount(tempList2, amount + 1);
//    }
//
//    if (result.size() > amount)
//    {
//        return std::vector<Point>(result.begin(), result.begin() + amount);
//    }
//    else
//    {
//        return result;
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::GetSagitalPointsFemur(int amount)
//{
//    std::vector<Point> result;
//    Point vectorTEA = mKnee.getMedialEpicondylePerp() - mKnee.getLateralEpicondyle();
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//
//    Plane sagital, transversal1, transversal2;
//
//    sagital.init(vectorTEA, mKnee.getTopPointOnPatellaPath());
//    transversal1.init(vectorAxis, mKnee.getTopPointOnPatellaPath());
//    transversal2.init(vectorAxis, mKnee.getAnteriorCortex());
//
//    if (transversal1.eval(mKnee.getHipCenter()) < 0)
//    {
//        transversal1.reverse();
//    }
//
//    if (transversal2.eval(mKnee.getHipCenter()) > 0)
//    {
//        transversal2.reverse();
//    }
//
//    vtkSmartPointer<vtkPolyData> sagitalContour = ImplantTools::getMaxContour(mKnee.GetFemurPoly(), sagital.getNormalVector(), sagital.getPoint());
//
//    vtkIdType nearA;
//
//    nearA = ImplantTools::GetNearestPoints(sagitalContour, mKnee.getTopPointOnPatellaPath());
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(sagitalContour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearA)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearA && pos == -1)
//    {
//        pos = cont - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return result;
//    }
//
//    std::vector<Point> tempList;
//
//    if (transversal1.eval(sagitalContour->GetPoint(allPoints[5])) < 0)
//    {
//        for (int i = tSize - 1; i > 0; i--)
//        {
//            double pnt[3];
//            sagitalContour->GetPoint(allPoints[i], pnt);
//
//            if (transversal2.eval(pnt) > 0)
//            {
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else
//            {
//                break;
//            }
//        }
//    }
//    else
//    {
//        for (int i = 1; i < tSize; i++)
//        {
//            double pnt[3];
//            sagitalContour->GetPoint(allPoints[i], pnt);
//
//            if (transversal2.eval(pnt) > 0)
//            {
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else
//            {
//                break;
//            }
//        }
//    }
//
//    result = reduceSortPointsByAmount(tempList, amount + 2);
//
//    if (result.size() > amount)
//    {
//        return std::vector<Point>(result.begin() + 1, result.begin() + amount + 1);
//    }
//    else
//    {
//        return result;
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::GetMedialBorderPointsFemur(int amount)
//{
//    std::vector<Point> result, resultFinal;
//
//    std::vector<Point> mainPoints = GetMedialMainPointsFemur(amount);
//
//    Point medialEpi = mKnee.getMedialEpicondyle();
//    Point vectorTEA = mKnee.getMedialEpicondylePerp() - mKnee.getLateralEpicondyle();
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//    Point vectorAP = vectorTEA.cross(vectorAxis);
//
//    int tSize = mainPoints.size();
//
//    if (tSize == 0)
//    {
//        return resultFinal;
//    }
//
//    double radius = (ImplantTools::getDistanceBetweenPoints(medialEpi, mainPoints[tSize - 1])) * 0.6;
//
//    Plane sagital;
//    sagital.init(vectorTEA, mKnee.getFemurKneeCenter());
//
//    Sphere mySphere(medialEpi, radius);
//
//    Point initVector, tempVector;
//
//    int rightPos = -1;
//
//    double lastAngle = -1;
//
//    std::vector<vtkSmartPointer<vtkPolyData>> contoursList;
//
//    for (int i = 0; i < tSize; i++)
//    {
//        Plane newPlane = sagital.getPerpendicularPlane(mainPoints[i], medialEpi);
//        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(medialSideFemur, newPlane.getNormalVector(), newPlane.getPoint());
//        vtkSmartPointer<vtkPoints> tPoints = contour->GetPoints();
//        contoursList.push_back(contour);
//        int pointSize = tPoints->GetNumberOfPoints();
//
//        if (pointSize < 20)
//        {
//            return resultFinal;
//        }
//
//        Line line1(sagital.getNormalVector(), mainPoints[i]);
//        Plane onEpicondyle;
//        onEpicondyle.init(sagital.getNormalVector(), medialEpi);
//
//        Point inter = onEpicondyle.getInterceptionLinePoint(line1);
//
//        Point vector1 = inter - mainPoints[i];
//        Point vector2 = inter - medialEpi;
//        vector1.normalice();
//        vector2.normalice();
//
//        Point epiPoint, basePoint;
//        epiPoint = inter + 1000.0 * vector1;
//        basePoint = inter + 1000.0 * vector2;
//
//        Point topBase = getNearPointToLine(Line(vector1, basePoint), tPoints, mySphere);
//        Point topEpi = getNearPointToLine(Line(vector2, epiPoint), tPoints, mySphere);
//
//        onEpicondyle.movePlane(topEpi);
//        inter = onEpicondyle.getInterceptionLinePoint(Line(vector1, topBase));
//
//        vector1 = inter - topBase;
//        vector2 = inter - topEpi;
//        vector1.normalice();
//        vector2.normalice();
//
//        Point a, b;
//        a = inter + 10.0 * vector1;
//        b = inter + 10.0 * vector2;
//
//        Line obliqueLine = Line::makeLineWithPoints(a, b);
//
//        Point nearPoint = getNearPointToLine(obliqueLine, tPoints, mySphere);
//
//        result.push_back(nearPoint);
//        if (i == 0)
//        {
//            initVector = nearPoint - mKnee.getMedialCoronalDistalPoint();
//        }
//        else if (rightPos < 0)
//        {
//            tempVector = nearPoint - mKnee.getMedialCoronalDistalPoint();
//            double temp = ImplantTools::getAngleBetweenVectorsDegree(initVector, tempVector);
//
//            if (lastAngle < 0)
//            {
//                lastAngle = temp;
//            }
//            else
//            {
//                double diff = abs(temp - lastAngle);
//                lastAngle = temp;
//                if (diff < 3)
//                {
//                    rightPos = i - 1;
//                }
//            }
//        }
//    }
//
//    lastMedialBorderFemur = result[tSize - 1];
//
//    if (rightPos > 0)
//    {
//        if (rightPos > 3)
//        {
//            rightPos = 3;
//        }
//
//        Line fixLine = Line::makeLineWithPoints(result[rightPos], mKnee.getMedialCoronalDistalPoint());
//        for (int i = 0; i < rightPos; i++)
//        {
//            Point nearPoint = getNearPointToLine(fixLine, contoursList[i]->GetPoints(), mySphere);
//            result[i] = nearPoint;
//        }
//    }
//
//    double radius2 = 5.0;
//
//    for (int i = 0; i < tSize; i++)
//    {
//        vtkIdType nearEpi;
//        nearEpi = ImplantTools::GetNearestPoints(contoursList[i], medialEpi);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(contoursList[i], lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//        it1 = lines.begin();
//        it2 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//
//        int pos = 0;
//        int cont = 0;
//
//        for (; it1 != it2; ++it1)
//        {
//            allPoints.push_back(it1->first);
//            if (it1->first == nearEpi)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        if (lines.back().second == nearEpi && pos == -1)
//        {
//            pos = cont - 1;
//        }
//
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//
//        Point decisionVector = medialEpi - result[i];
//        Plane decisionPlane;
//        decisionPlane.init(decisionVector, medialEpi);
//        decisionPlane.reverseByPoint(result[i]);
//
//        double decisionPnt[3];
//        contoursList[i]->GetPoint(allPoints[5], decisionPnt);
//
//        if (decisionPlane.eval(decisionPnt) < 0)
//        {
//            std::reverse(allPoints.begin(), allPoints.end());
//        }
//
//        Sphere mySphereBorder(result[i], radius2);
//
//        for (int j = 0; j < allPoints.size() - 1; j++)
//        {
//            double pnt1[3];
//            contoursList[i]->GetPoint(allPoints[j], pnt1);
//            Point checkPoint1(pnt1[0], pnt1[1], pnt1[2]);
//
//            double pnt2[3];
//            contoursList[i]->GetPoint(allPoints[j + 1], pnt2);
//            Point checkPoint2(pnt2[0], pnt2[1], pnt2[2]);
//
//            if ((mySphereBorder.isPointInside(checkPoint1) == false && mySphereBorder.isPointInside(checkPoint2) == true) ||
//                (mySphereBorder.isPointInside(checkPoint1) == true && mySphereBorder.isPointInside(checkPoint2) == false))
//            {
//                Line myLine = Line::makeLineWithPoints(checkPoint1, checkPoint2);
//                std::pair<Point, Point> intercep;
//
//                myLine.getInterceptionSphere(mySphereBorder.center, mySphereBorder.radius, intercep);
//
//                Point outSise;
//                if (mySphereBorder.isPointInside(checkPoint1) == false)
//                {
//                    outSise = checkPoint1;
//                }
//                else
//                {
//                    outSise = checkPoint2;
//                }
//
//                if (ImplantTools::getDistanceBetweenPoints(outSise, intercep.first, true) < ImplantTools::getDistanceBetweenPoints(outSise, intercep.second, true))
//                {
//                    resultFinal.push_back(intercep.first);
//                }
//                else
//                {
//                    resultFinal.push_back(intercep.second);
//                }
//            }
//
//        }
//
//    }
//
//    return resultFinal;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTibiaMedialObliquePoints(int amount)
//{
//    std::vector<Point> result;
//    if (mGoodMediaRef == false)
//    {
//        return result;
//    }
//    result = GetObliquePointsTibia(mMedialTibiaRef, mObliqueRef, amount);
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTibiaLateralObliquePoints(int amount)
//{
//    std::vector<Point> result;
//    if (mGoodLateralRef == false)
//    {
//        return result;
//    }
//
//    result = GetObliquePointsTibia(mLateralTibiaRef, mObliqueRef, amount);
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTibiaOnTuberAxisPoints(int amount)
//{
//    std::vector<Point> result;
//    result = GetObliquePointsTibia(mTuberUpRef, mTuberDownRef, amount);
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTibiaOnTuberSidePoints(int amount)
//{
//    std::vector<Point> result;
//    result = GetTransversalPointsTibia(mTuberDownRef, mTuberSideRef, amount);
//    return result;
//}
//
//Point RegistrationPointsVTK::GetTibiaCenterEllipse()
//{
//    return mCenterEllipse;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTibiaTemplatePointsRight()
//{
//    std::vector<Point> result;
//
//    Plane transverseEllipse, transverseOnTuber;
//    Plane coronal;
//    Plane sagitalLat, sagitalMed, sagitalCenter;
//
//    transverseEllipse.init(mVectorTibiaAxisUp, mCenterEllipse);
//    transverseOnTuber.init(mVectorTibiaAxisUp, (myTuber + 2.0 * mVectorTibiaAxisUp));
//    coronal.init(mVectorTibiaAPFront, mCenterEllipse);
//    sagitalLat.init(mVectorTibiaTEALat, mBigCenterLat);
//    sagitalMed.init(-mVectorTibiaTEALat, mBigCenterMed);
//    sagitalCenter.init(mVectorTibiaTEALat, (myTuber + 3.0 * mVectorTibiaTEALat));
//
//    vtkSmartPointer<vtkPolyData> contourTransverse1 = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), transverseEllipse.getNormalVector(), transverseEllipse.getPoint());
//    vtkSmartPointer<vtkPolyData> contourTransverse2 = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), transverseOnTuber.getNormalVector(), transverseOnTuber.getPoint());
//
//    vtkSmartPointer<vtkPolyData> contourCoronal = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), coronal.getNormalVector(), coronal.getPoint());
//
//    vtkSmartPointer<vtkPolyData> contourSagitalLat = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), sagitalLat.getNormalVector(), sagitalLat.getPoint());
//    vtkSmartPointer<vtkPolyData> contourSagitalMed = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), sagitalMed.getNormalVector(), sagitalMed.getPoint());
//
//    vtkSmartPointer<vtkPolyData> contourSagitalCenter = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), sagitalCenter.getNormalVector(), sagitalCenter.getPoint());
//
//    std::vector<Plane> conditions;
//    conditions.push_back(coronal);
//
//    std::vector<Point> temp;
//
//    ImplantTools::GetPointsOnContourSort(contourTransverse1, conditions, coronal, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    ImplantTools::GetPointsOnContourSort(contourTransverse2, conditions, coronal, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    conditions.clear();
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(sagitalLat);
//
//    ImplantTools::GetPointsOnContourSort(contourCoronal, conditions, mBigCenterLat, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//    
//    conditions.clear();
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(sagitalMed);
//
//    ImplantTools::GetPointsOnContourSort(contourCoronal, conditions, mBigCenterMed, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//    
//    conditions.clear();
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(coronal);
//
//    ImplantTools::GetPointsOnContourSort(contourSagitalLat, conditions, mBigCenterLat, temp);
//    std::vector<Point> resultTemp2;
//    reduceSortPointsByRadius(temp, 5, resultTemp2);
//    temp.clear();
//    for (int i = 0; i < 4; i++)
//    {
//        result.push_back(resultTemp2[i]);
//    }
//    
//    ImplantTools::GetPointsOnContourSort(contourSagitalMed, conditions, mBigCenterMed, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    
//    conditions.clear();
//    coronal.movePlaneOnNormal(10);
//    transverseOnTuber.movePlaneOnNormal(-10);
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(coronal);
//
//    ImplantTools::GetPointsOnContourSort(contourSagitalCenter, conditions, coronal, temp);
//    std::vector<Point> resultTemp1;
//    reduceSortPointsByRadius(temp, 5, resultTemp1);
//    temp.clear();
//    for (int i = 1; i < resultTemp1.size(); i++)
//    {
//        result.push_back(resultTemp1[i]);
//    }
//
//    Point a(17.58, -49.17, 743.10);
//    Point b(15.58, -52.17, 745.10);
//
//    Point c(28.28, -51.17, 747.10);
//    Point d(31.28, -53.67, 746.10);
//    Point e(34.28, -56.17, 746.10);
//    Point f(37.28, -58.67, 745.10);
//
//    Point myPointTemp;
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (a + 5.0 * mVectorTibiaAxisUp), (a - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (b + 5.0 * mVectorTibiaAxisUp), (b - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//    
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (c + 5.0 * mVectorTibiaAxisUp), (c - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (d + 5.0 * mVectorTibiaAxisUp), (d - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (e + 5.0 * mVectorTibiaAxisUp), (e - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (f + 5.0 * mVectorTibiaAxisUp), (f - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    return result;
//}
//
//void RegistrationPointsVTK::GetFemurObliqueLateralEpi(const std::vector<Point>& pRef, std::vector<Point>& result, double pRadius)
//{
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//    Point vectorTEA = mKnee.getLateralEpicondyle() - mKnee.getMedialEpicondylePerp();
//    Point vectorAP = vectorAxis.cross(vectorTEA);
//
//    Plane sagital;
//    sagital.init(vectorTEA, mKnee.getFemurKneeCenter());
//
//    for (int i = 0; i < pRef.size(); i++)
//    {
//        Plane myPlane = sagital.getPerpendicularPlane(pRef[i], mKnee.getLateralEpicondyle());
//        vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(lateralSideFemur, myPlane.getNormalVector(), myPlane.getPoint());
//
//        vtkIdType nearA, nearB;
//
//        nearA = ImplantTools::GetNearestPoints(mainContour, mKnee.getLateralEpicondyle());
//        nearB = ImplantTools::GetNearestPoints(mainContour, pRef[i]);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(mainContour, lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//        it1 = lines.begin();
//        it2 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//        int cont = 0;
//        int pos = -1;
//
//        for (; it1 != it2; ++it1)
//        {
//            allPoints.push_back(it1->first);
//            if (it1->first == nearB)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        int tSize = cont;
//
//        if (lines.back().second == nearB && pos == -1)
//        {
//            pos = tSize - 1;
//        }
//
//        if (pos > 0)
//        {
//            std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//        }
//
//        Plane transversePoint;
//        transversePoint.init(vectorTEA, pRef[i]);
//        if (transversePoint.eval(mKnee.getMedialEpicondyle()) > 0)
//        {
//            transversePoint.reverse();
//        }
//
//        std::vector<Point> tempList;
//
//        if (transversePoint.eval(mainContour->GetPoint(allPoints[2])) < 0)
//        {
//            for (int j = tSize - 1; j > 0; j--)
//            {
//                if (allPoints[j] == nearA)
//                {
//                    break;
//                }
//                double pnt[3];
//                mainContour->GetPoint(allPoints[j], pnt);
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//        else
//        {
//            for (int j = 1; j < tSize; j++)
//            {
//                if (allPoints[j] == nearA)
//                {
//                    break;
//                }
//                double pnt[3];
//                mainContour->GetPoint(allPoints[j], pnt);
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//
//        std::vector<Point> resultTemp;
//        reduceSortPointsByRadius(tempList, pRadius, resultTemp);
//        for (int j = 0; j < resultTemp.size(); j++)
//        {
//            double myDistance = ImplantTools::getDistanceBetweenPoints(mKnee.getLateralEpicondyle(), resultTemp[j]);
//            if (myDistance > 10)
//            {
//                result.push_back(resultTemp[j]);
//            }
//        }
//    }
//}
//
//void RegistrationPointsVTK::GetFemurObliqueMedialEpi(const std::vector<Point>& pRef, std::vector<Point>& result, double pRadius)
//{
//    Point vectorAxis = mKnee.getHipCenter() - mKnee.getFemurKneeCenter();
//    Point vectorTEA = mKnee.getLateralEpicondyle() - mKnee.getMedialEpicondylePerp();
//    Point vectorAP = vectorAxis.cross(vectorTEA);
//
//    Plane sagital;
//    sagital.init(vectorTEA, mKnee.getFemurKneeCenter());
//
//    for (int i = 0; i < pRef.size(); i++)
//    {
//        Plane myPlane = sagital.getPerpendicularPlane(pRef[i], mKnee.getMedialEpicondyle());
//        vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(medialSideFemur, myPlane.getNormalVector(), myPlane.getPoint());
//
//        vtkIdType nearA, nearB;
//
//        nearA = ImplantTools::GetNearestPoints(mainContour, mKnee.getMedialEpicondyle());
//        nearB = ImplantTools::GetNearestPoints(mainContour, pRef[i]);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(mainContour, lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//        it1 = lines.begin();
//        it2 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//        int cont = 0;
//        int pos = -1;
//
//        for (; it1 != it2; ++it1)
//        {
//            allPoints.push_back(it1->first);
//            if (it1->first == nearB)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        int tSize = cont;
//
//        if (lines.back().second == nearB && pos == -1)
//        {
//            pos = tSize - 1;
//        }
//
//        if (pos > 0)
//        {
//            std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//        }
//
//        Plane transversePoint;
//        transversePoint.init(vectorTEA, pRef[i]);
//        if (transversePoint.eval(mKnee.getLateralEpicondyle()) > 0)
//        {
//            transversePoint.reverse();
//        }
//
//        std::vector<Point> tempList;
//
//        if (transversePoint.eval(mainContour->GetPoint(allPoints[2])) < 0)
//        {
//            for (int j = tSize - 1; j > 0; j--)
//            {
//                if (allPoints[j] == nearA)
//                {
//                    break;
//                }
//                double pnt[3];
//                mainContour->GetPoint(allPoints[j], pnt);
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//        else
//        {
//            for (int j = 1; j < tSize; j++)
//            {
//                if (allPoints[j] == nearA)
//                {
//                    break;
//                }
//                double pnt[3];
//                mainContour->GetPoint(allPoints[j], pnt);
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//        std::reverse(tempList.begin(), tempList.end());
//        std::vector<Point> resultTemp;
//        reduceSortPointsByRadius(tempList, pRadius, resultTemp);
//        for (int j = 0; j < resultTemp.size(); j++)
//        {
//            double myDistance1 = ImplantTools::getDistanceBetweenPoints(pRef[i], resultTemp[j]);
//            double myDistance2 = ImplantTools::getDistanceBetweenPoints(mKnee.getMedialEpicondyle(), resultTemp[j]);
//            if (myDistance1 > 10 && myDistance2 > 10)
//            {
//                result.push_back(resultTemp[j]);
//            }
//        }
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::GetFemurTemplatePoints(bool isLeft)
//{
//    GetLateralBorderPointsFemur();
//    GetMedialBorderPointsFemur();
//
//    std::vector<Point> result, temp, temp2;
//
//    temp = GetLateralMainPointsFemur(0, 1);
//    reduceSortPointsByRadius(temp, 5, temp2);
//
//    for (int i = 0; i < temp2.size(); i++)
//    {
//        result.push_back(temp2[i]);
//    }
//
//    if (isLeft == true)
//    {
//        GetFemurObliqueLateralEpi(temp2, result, 15);
//    }
//    else
//    {
//        GetFemurObliqueLateralEpi(temp2, result, 15);
//    }
//    
//    temp.clear();
//    temp2.clear();
//
//    temp = GetMedialMainPointsFemur(0, 1);
//    reduceSortPointsByRadius(temp, 5, temp2);
//
//    for (int i = 0; i < temp2.size(); i++)
//    {
//        result.push_back(temp2[i]);
//    }
//
//    if (isLeft == true)
//    {
//        GetFemurObliqueMedialEpi(temp2, result, 15);
//    }
//    else
//    {
//        GetFemurObliqueMedialEpi(temp2, result, 20);
//    }
//
//    temp.clear();
//    temp2.clear();
//    
//
//    temp = GetMedialMainPointsFemur(0, 2);
//    reduceSortPointsByRadius(temp, 5, temp2);
//
//    for (int i = 0; i < temp2.size(); i++)
//    {
//        result.push_back(temp2[i]);
//    }
//
//    if (isLeft == true)
//    {
//        GetFemurObliqueMedialEpi(temp2, result, 15);
//    }
//    else
//    {
//        GetFemurObliqueMedialEpi(temp2, result, 20);
//    }
//
//    temp.clear();
//    temp2.clear();
//
//    temp = GetLateralObliquePointsFemur();
//    for (int i = 1; i < temp.size(); i++)
//    {
//        result.push_back(temp[i]);
//    }
//    temp.clear();
//
//    vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
//    implicitPolyDataDistance->SetInput(mKnee.GetFemurPoly());
//
//    if (isLeft == true)
//    {
//        temp.push_back(Point(9.06, 30.28, -257.41));
//        temp.push_back(Point(9.06, 32.28, -262.41));
//        temp.push_back(Point(9.06, 34.28, -267.41));
//
//        temp.push_back(Point(-15.94, 34.28, -254.41));
//        temp.push_back(Point(-15.94, 33.28, -259.41));
//
//        temp.push_back(Point(-29.94, 54.28, -256.41));
//        temp.push_back(Point(44.06, 60.28, -257.41));
//
//        for (int i = 0; i < temp.size(); i++)
//        {
//            double myClosest[3];
//            double pnt[3] = { temp[i].x, temp[i].y, temp[i].z };
//            implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
//            result.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
//        }
//    }
//    else
//    {
//        temp.push_back(Point(24.48, -75.11, 779.70));
//        temp.push_back(Point(24.48, -72.11, 774.70));
//        temp.push_back(Point(23.48, -66.11, 769.70));
//
//        temp.push_back(Point(51.48, -71.11, 785.70));
//        temp.push_back(Point(51.48, -72.11, 780.70));
//        temp.push_back(Point(51.48, -72.11, 775.70));
//        temp.push_back(Point(51.48, -71.11, 770.70));
//
//        temp.push_back(Point(67.23, -44.11, 782.70));
//        temp.push_back(Point(-10.52, -46.11, 782.70));
//
//        for (int i = 0; i < temp.size(); i++)
//        {
//            double myClosest[3];
//            double pnt[3] = { temp[i].x, temp[i].y, temp[i].z };
//            implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
//            result.push_back(Point(myClosest[0], myClosest[1], myClosest[2]));
//        }
//    }
//
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTibiaTemplatePoints()
//{
//    std::vector<Point> result;
//
//    Plane transverseEllipse, transverseOnTuber;
//    Plane coronal;
//    Plane sagitalLat, sagitalMed, sagitalCenter;
//
//    transverseEllipse.init(mVectorTibiaAxisUp, mCenterEllipse);
//    transverseOnTuber.init(mVectorTibiaAxisUp, (myTuber - 2.0 * mVectorTibiaAxisUp));
//    coronal.init(mVectorTibiaAPFront, mCenterEllipse);
//    sagitalLat.init(mVectorTibiaTEALat, mBigCenterLat);
//    sagitalMed.init(-mVectorTibiaTEALat, mBigCenterMed);
//    sagitalCenter.init(mVectorTibiaTEALat, (myTuber + 3.0 * mVectorTibiaTEALat));
//
//    vtkSmartPointer<vtkPolyData> contourTransverse1 = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), transverseEllipse.getNormalVector(), transverseEllipse.getPoint());
//    vtkSmartPointer<vtkPolyData> contourTransverse2 = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), transverseOnTuber.getNormalVector(), transverseOnTuber.getPoint());
//
//    vtkSmartPointer<vtkPolyData> contourCoronal = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), coronal.getNormalVector(), coronal.getPoint());
//    
//    vtkSmartPointer<vtkPolyData> contourSagitalLat = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), sagitalLat.getNormalVector(), sagitalLat.getPoint());
//    vtkSmartPointer<vtkPolyData> contourSagitalMed = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), sagitalMed.getNormalVector(), sagitalMed.getPoint());
//
//    vtkSmartPointer<vtkPolyData> contourSagitalCenter = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), sagitalCenter.getNormalVector(), sagitalCenter.getPoint());
//
//    std::vector<Plane> conditions;
//    conditions.push_back(coronal);
//
//    std::vector<Point> temp;
//
//    ImplantTools::GetPointsOnContourSort(contourTransverse1, conditions, coronal, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    ImplantTools::GetPointsOnContourSort(contourTransverse2, conditions, coronal, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    conditions.clear();
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(sagitalLat);
//
//    ImplantTools::GetPointsOnContourSort(contourCoronal, conditions, mBigCenterLat, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    conditions.clear();
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(sagitalMed);
//
//    ImplantTools::GetPointsOnContourSort(contourCoronal, conditions, mBigCenterMed, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//    conditions.clear();
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(coronal);
//
//    ImplantTools::GetPointsOnContourSort(contourSagitalLat, conditions, mBigCenterLat, temp);
//
//    std::vector<Point> resultTemp2;
//    reduceSortPointsByRadius(temp, 5, resultTemp2);
//    temp.clear();
//    for (int i = 0; i < 4; i++)
//    {
//        result.push_back(resultTemp2[i]);
//    }
//
//    ImplantTools::GetPointsOnContourSort(contourSagitalMed, conditions, mBigCenterMed, temp);
//    reduceSortPointsByRadius(temp, 5, result);
//    temp.clear();
//
//
//    conditions.clear();
//    coronal.movePlaneOnNormal(10);
//    transverseOnTuber.movePlaneOnNormal(-10);
//    conditions.push_back(transverseOnTuber);
//    conditions.push_back(coronal);
//
//    ImplantTools::GetPointsOnContourSort(contourSagitalCenter, conditions, coronal, temp);
//    std::vector<Point> resultTemp1;
//    reduceSortPointsByRadius(temp, 5, resultTemp1);
//    temp.clear();
//    for (int i = 1; i < resultTemp1.size(); i++)
//    {
//        result.push_back(resultTemp1[i]);
//    }
//
//    Point a(6.42, 48.98, -288.42);
//    Point b(8.42, 43.98, -288.42);
//
//    Point c(-3.08, 44.98, -286.42);
//    Point d(-5.58, 43.98, -286.42);
//    Point e(-8.08, 42.98, -286.42);
//    Point f(-10.58, 41.98, -287.42);
//
//    Point myPointTemp;
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (a + 5.0 * mVectorTibiaAxisUp), (a - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (b + 5.0 * mVectorTibiaAxisUp), (b - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (c + 5.0 * mVectorTibiaAxisUp), (c - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (d + 5.0 * mVectorTibiaAxisUp), (d - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (e + 5.0 * mVectorTibiaAxisUp), (e - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), (f + 5.0 * mVectorTibiaAxisUp), (f - 5.0 * mVectorTibiaAxisUp), myPointTemp);
//    result.push_back(myPointTemp);
//
//    return result;
//}
//
//std::pair<Point, double> RegistrationPointsVTK::getLateralPlateauCirclePoints(std::vector<Point>& circle)
//{
//    mGoodLateralRef = false;
//    double moveFar = 1000.0;
//    Plane transverse;
//    Point vectorAxis = mVectorTibiaAxisUp;
//    //transverse.init(vectorAxis, mCenterEllipse);
//    transverse.init(mTibiaTransverse.getNormalVector(), mTibiaTransverse.getPoint());
//    transverse.reverseByPoint(mKnee.getAnkleCenter(), false);
//
//    Point vectorAP = mVectorTibiaAPFront;
//    Point vectorTEA = mVectorTibiaTEALat;
//
//    double angleSign = 1.0;
//
//    if (mKnee.getIsRight() == true)
//    {
//        angleSign = -1.0;
//    }
//
//    std::vector<double> angles = { -10, 40, 90 };
//    std::vector<Point> result;
//
//    for (const auto myAngle : angles)
//    {
//        double rad = angleSign * (myAngle / 180.0) * PI;
//        cv::Mat rotMat = ImplantTools::getRotateMatrix(vectorAxis, rad);
//        cv::Mat teaMat = rotMat * vectorTEA.ToMatPoint();
//        Point tempTEA = Point(teaMat);
//        tempTEA.normalice();
//
//        vectorAP = vectorAxis.cross(tempTEA);
//        vectorAP.normalice();
//
//        Point refSagital = mBigCenterLat + moveFar * tempTEA;
//        Point refSagital2 = mBigCenterLat + 0.7 * moveFar * vectorAxis;
//
//        Plane coronal, sagital;
//        coronal.init(vectorAP, mBigCenterLat);
//        sagital.init(tempTEA, mBigCenterLat);
//        sagital.reverseByPoint(refSagital);
//
//        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), coronal.getNormalVector(), coronal.getPoint());
//        vtkIdType nearPoint = ImplantTools::GetNearestPoints(contour, Line::makeLineWithPoints(refSagital, refSagital2), transverse, sagital);
//        double pnt[3];
//        contour->GetPoint(nearPoint, pnt);
//        result.push_back(Point(pnt[0], pnt[1], pnt[2]));
//    }
//
//    Plane circlePlane;
//    circlePlane.init(vectorAxis, mBigCenterLat);
//
//    cv::Mat rotZ = ImplantTools::GetRotateZ(circlePlane.getNormalVector());
//    auto it1 = result.begin();
//    auto it2 = result.end();
//    double Z = 0;
//    std::vector<cv::Point2f> coplanar;
//
//    for (; it1 != it2; ++it1)
//    {
//        Point proj = circlePlane.getProjectionPoint(*it1);
//        cv::Mat pointMat = rotZ * proj.ToMatPoint();
//        Point temp = Point(pointMat);
//        coplanar.push_back(cv::Point2f(temp.x, temp.y));
//        Z = temp.z;
//    }
//
//    cv::Point2d center;
//    double radius;
//    std::pair<cv::Point2d, double> myCircle = ImplantTools::findCircle(coplanar[0], coplanar[1], coplanar[2]);
//    center = myCircle.first;
//    cv::Mat pointMat = rotZ.inv() * Point(center.x, center.y, Z).ToMatPoint();
//
//    Point finalCenter = Point(pointMat);
//    /*Point finalCenterOnSurface;
//    Point a = finalCenter + 100.0 * vectorAxis;
//    Point b = finalCenter - 10.0 * vectorAxis;
//    bool interceptionResult = ImplantTools::GetInterceptionWithSegment(mKnee.GetTibiaPoly(), a, b, finalCenterOnSurface);
//    if (interceptionResult == false)
//    {
//        Line tempLine = Line::makeLineWithPoints(a, b);
//        finalCenterOnSurface = tempLine.getProjectPoint(finalCenter);
//    }*/
//
//    radius = 0.8 * myCircle.second;
//
//    Point newAP = circlePlane.getProjectionVector(mVectorTibiaAPFront);
//    Point newTEA = vectorTEA;
//    newAP.normalice();
//    Point frontA = finalCenter + radius * newAP;
//    Point sideB = finalCenter + radius * newTEA;
//    double step = 1.0 / 5.0;
//
//    transverse.movePlane(mTransverseRef);
//
//    for (double i = 0; i < 5; i++)
//    {
//        Point temp = frontA + i * step * (sideB - frontA);
//        Point vectorTemp = temp - finalCenter;
//        vectorTemp.normalice();
//        Point pointOnCircle = finalCenter + radius * vectorTemp;
//
//        Plane oblique = circlePlane.getPerpendicularPlane(pointOnCircle, finalCenter);
//        Plane obliqueRef;
//        obliqueRef.init(vectorTemp, finalCenter);
//        obliqueRef.reverseByPoint(pointOnCircle);
//
//        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), oblique.getNormalVector(), oblique.getPoint());
//
//        vtkIdType nearID = ImplantTools::GetNearestPoints(contour, finalCenter);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(contour, lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it11, it22;
//        it11 = lines.begin();
//        it22 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//        int cont = 0;
//        int pos = -1;
//
//        for (; it11 != it22; ++it11)
//        {
//            allPoints.push_back(it11->first);
//            if (it11->first == nearID)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        int tSize = cont;
//
//        if (lines.back().second == nearID && pos == -1)
//        {
//            pos = cont - 1;
//        }
//
//        if (pos == -1)
//        {
//            continue;
//        }
//
//        if (pos > 0)
//        {
//            std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//        }
//
//        if (obliqueRef.eval(contour->GetPoint(allPoints[10])) < 0)
//        {
//            std::reverse(allPoints.begin(), allPoints.end());
//        }
//
//        bool finishSphere = false;
//
//        if (i == 0)
//        {
//            //finishSphere = true;
//        }
//
//        for (int j = 0; j < tSize - 1; j++)
//        {
//            double pnt1[3];
//            double pnt2[3];
//            contour->GetPoint(allPoints[j], pnt1);
//            contour->GetPoint(allPoints[j + 1], pnt2);
//
//            Point a = Point(pnt1[0], pnt1[1], pnt1[2]);
//            Point b = Point(pnt2[0], pnt2[1], pnt2[2]);
//
//            if ((ImplantTools::isPointInsideSphere(finalCenter, radius, a) == true && ImplantTools::isPointInsideSphere(finalCenter, radius, b) == false ||
//                ImplantTools::isPointInsideSphere(finalCenter, radius, a) == false && ImplantTools::isPointInsideSphere(finalCenter, radius, b) == true) &&
//                finishSphere == false)
//            {
//                Point exterior;
//
//                if (ImplantTools::isPointInsideSphere(finalCenter, radius, a) == false)
//                {
//                    exterior = a;
//                }
//                else
//                {
//                    exterior = b;
//                }
//
//                std::pair<Point, Point> resultIntercep;
//                Line intercep = Line::makeLineWithPoints(a, b);
//                intercep.getInterceptionSphere(finalCenter, radius, resultIntercep);
//
//                if (ImplantTools::getDistanceBetweenPoints(resultIntercep.first, exterior) < ImplantTools::getDistanceBetweenPoints(resultIntercep.second, exterior))
//                {
//                    circle.push_back(resultIntercep.first);
//                }
//                else
//                {
//                    circle.push_back(resultIntercep.second);
//                }
//                
//                finishSphere = true;
//            }
//
//            if (finishSphere == true)
//            {
//                if ((transverse.eval(a) > 0 && transverse.eval(b) < 0) ||
//                    (transverse.eval(a) < 0 && transverse.eval(b) > 0) ||
//                    transverse.eval(a) == 0 || transverse.eval(b) == 0)
//                {
//                    Line intercep = Line::makeLineWithPoints(a, b);
//
//                    if (i == 0)
//                    {
//                        mLateralTibiaRef = transverse.getInterceptionLinePoint(intercep);
//                        mGoodLateralRef = true;
//                    }
//                    
//                    break;
//                }
//            }
//
//        }
//
//    }
//
//    return std::make_pair(finalCenter, radius);
//}
//
//std::pair<Point, double> RegistrationPointsVTK::getMedialPlateauCirclePoints(std::vector<Point>& circle, std::vector<Point>& transversePoints)
//{
//    mGoodMediaRef = false;
//    double moveFar = 1000.0;
//    Plane transverse;
//    Point vectorAxis = mVectorTibiaAxisUp;
//    //transverse.init(vectorAxis, mCenterEllipse);
//    transverse.init(mTibiaTransverse.getNormalVector(), mTibiaTransverse.getPoint());
//    transverse.reverseByPoint(mKnee.getAnkleCenter(), false);
//
//    Point vectorAP = mVectorTibiaAPFront;
//    Point vectorTEA = -mVectorTibiaTEALat;
//
//    double angleSign = 1.0;
//
//    if (mKnee.getIsRight() == false)
//    {
//        angleSign = -1.0;
//    }
//
//    std::vector<double> angles = { -10, 40, 90 };
//    std::vector<Point> result;
//
//    for (const auto myAngle : angles)
//    {
//        double rad = angleSign * (myAngle / 180.0) * PI;
//        cv::Mat rotMat = ImplantTools::getRotateMatrix(vectorAxis, rad);
//        cv::Mat teaMat = rotMat * vectorTEA.ToMatPoint();
//        Point tempTEA = Point(teaMat);
//        tempTEA.normalice();
//
//        vectorAP = vectorAxis.cross(tempTEA);
//        vectorAP.normalice();
//
//        Point refSagital = mBigCenterMed + moveFar * tempTEA;
//        Point refSagital2 = mBigCenterMed + 0.7 * moveFar * vectorAxis;
//
//        Plane coronal, sagital;
//        coronal.init(vectorAP, mBigCenterMed);
//        sagital.init(tempTEA, mBigCenterMed);
//        sagital.reverseByPoint(refSagital);
//
//        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), coronal.getNormalVector(), coronal.getPoint());
//        vtkIdType nearPoint = ImplantTools::GetNearestPoints(contour, Line::makeLineWithPoints(refSagital, refSagital2), transverse, sagital);
//        double pnt[3];
//        contour->GetPoint(nearPoint, pnt);
//        result.push_back(Point(pnt[0], pnt[1], pnt[2]));
//    }
//
//    Plane circlePlane;
//    circlePlane.init(vectorAxis, mBigCenterMed);
//
//    cv::Mat rotZ = ImplantTools::GetRotateZ(circlePlane.getNormalVector());
//    auto it1 = result.begin();
//    auto it2 = result.end();
//    double Z = 0;
//    std::vector<cv::Point2f> coplanar;
//
//    for (; it1 != it2; ++it1)
//    {
//        Point proj = circlePlane.getProjectionPoint(*it1);
//        cv::Mat pointMat = rotZ * proj.ToMatPoint();
//        Point temp = Point(pointMat);
//        coplanar.push_back(cv::Point2f(temp.x, temp.y));
//        Z = temp.z;
//    }
//
//    cv::Point2d center;
//    double radius;
//    std::pair<cv::Point2d, double> myCircle = ImplantTools::findCircle(coplanar[0], coplanar[1], coplanar[2]);
//    center = myCircle.first;
//    cv::Mat pointMat = rotZ.inv() * Point(center.x, center.y, Z).ToMatPoint();
//
//    Point finalCenter = Point(pointMat);
//    radius = 0.8 * myCircle.second;
//
//    Point newAP = circlePlane.getProjectionVector(mVectorTibiaAPFront);
//    Point newTEA = vectorTEA;
//    newAP.normalice();
//    Point frontA = finalCenter + radius * newAP;
//    Point sideB = finalCenter + radius * newTEA;
//    double step = 1.0 / 6.0;
//
//    transverse.movePlane(mTransverseRef);
//
//    for (double i = 0; i < 6; i++)
//    {
//        Point temp = frontA + i * step * (sideB - frontA);
//        Point vectorTemp = temp - finalCenter;
//        vectorTemp.normalice();
//        Point pointOnCircle = finalCenter + radius * vectorTemp;
//
//        Plane oblique = circlePlane.getPerpendicularPlane(pointOnCircle, finalCenter);
//        Plane obliqueRef;
//        obliqueRef.init(vectorTemp, finalCenter);
//        obliqueRef.reverseByPoint(pointOnCircle);
//
//        vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), oblique.getNormalVector(), oblique.getPoint());
//
//        vtkIdType nearID = ImplantTools::GetNearestPoints(contour, finalCenter);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(contour, lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it11, it22;
//        it11 = lines.begin();
//        it22 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//        int cont = 0;
//        int pos = -1;
//
//        for (; it11 != it22; ++it11)
//        {
//            allPoints.push_back(it11->first);
//            if (it11->first == nearID)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        int tSize = cont;
//
//        if (lines.back().second == nearID && pos == -1)
//        {
//            pos = cont - 1;
//        }
//
//        if (pos == -1)
//        {
//            continue;
//        }
//
//        if (pos > 0)
//        {
//            std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//        }
//
//        if (obliqueRef.eval(contour->GetPoint(allPoints[10])) < 0)
//        {
//            std::reverse(allPoints.begin(), allPoints.end());
//        }
//
//        bool finishSphere = false;
//
//        if (i == 0)
//        {
//            finishSphere = true;
//        }
//
//        for (int j = 0; j < (tSize - 1); j++)
//        {
//            double pnt1[3];
//            double pnt2[3];
//            contour->GetPoint(allPoints[j], pnt1);
//            contour->GetPoint(allPoints[j + 1], pnt2);
//
//            Point a = Point(pnt1[0], pnt1[1], pnt1[2]);
//            Point b = Point(pnt2[0], pnt2[1], pnt2[2]);
//
//            if ((ImplantTools::isPointInsideSphere(finalCenter, radius, a) == true && ImplantTools::isPointInsideSphere(finalCenter, radius, b) == false ||
//                ImplantTools::isPointInsideSphere(finalCenter, radius, a) == false && ImplantTools::isPointInsideSphere(finalCenter, radius, b) == true) && 
//                finishSphere == false)
//            {
//                Point exterior;
//
//                if (ImplantTools::isPointInsideSphere(finalCenter, radius, a) == false)
//                {
//                    exterior = a;
//                }
//                else
//                {
//                    exterior = b;
//                }
//
//                std::pair<Point, Point> resultIntercep;
//                Line intercep = Line::makeLineWithPoints(a, b);
//                intercep.getInterceptionSphere(finalCenter, radius, resultIntercep);
//
//                if (ImplantTools::getDistanceBetweenPoints(resultIntercep.first, exterior) < ImplantTools::getDistanceBetweenPoints(resultIntercep.second, exterior))
//                {
//                    circle.push_back(resultIntercep.first);
//                }
//                else
//                {
//                    circle.push_back(resultIntercep.second);
//                }
//                finishSphere = true;
//
//            }
//
//            if (transversePoints.size() == 4)
//            {
//                continue;
//            }
//
//            if (finishSphere == true)
//            {
//                if ((transverse.eval(a) > 0 && transverse.eval(b) < 0) ||
//                    (transverse.eval(a) < 0 && transverse.eval(b) > 0) || 
//                    transverse.eval(a) == 0 || transverse.eval(b) == 0)
//                {
//                    Line intercep = Line::makeLineWithPoints(a, b);
//                    if (i == 0)
//                    {
//                        mMedialTibiaRef = transverse.getInterceptionLinePoint(intercep);
//                        mGoodMediaRef = true;
//                    }
//                    else
//                    {
//                        transversePoints.push_back(transverse.getInterceptionLinePoint(intercep));
//                    }
//
//                    break;
//                }
//            }
//
//        }
//    }
//
//    //for (double i = 0; i < 2.0 * PI; i += 0.3)
//    //{
//    //    double X = center.x + radius * cos(i);
//    //    double Y = center.y + radius * sin(i);
//    //    Point temp = Point(X, Y, Z);
//    //    pointMat = rotZ.inv() * temp.ToMatPoint();
//    //    circle.push_back(Point(pointMat));
//    //}
//
//    return std::make_pair(finalCenter, radius);
//}
//
//std::vector<Point> RegistrationPointsVTK::GetTransversalPointsTibia(const Point& fixPoint, const Point& movePoint, int amount)
//{
//    std::vector<Point> result;
//
//    Plane coronalPlane, sagitalFix, sagitalMove;
//    coronalPlane.init(mVectorTibiaAPFront, (mKnee.getLateralPlateau() + mKnee.getMedialPlateau()) / 2.0);
//
//    sagitalFix.init(mVectorTibiaTEALat, fixPoint);
//    sagitalFix.reverseByPoint(movePoint);
//
//    sagitalMove.init(mVectorTibiaTEALat, movePoint);
//    sagitalMove.reverseByPoint(fixPoint);
//
//    Plane oblique = coronalPlane.getPerpendicularPlane(fixPoint, movePoint);
//
//    vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), oblique.getNormalVector(), oblique.getPoint());
//
//    vtkIdType nearA;
//
//    nearA = ImplantTools::GetNearestPoints(mainContour, fixPoint);
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(mainContour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearA)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearA && pos == -1)
//    {
//        pos = cont - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return result;
//    }
//
//    std::vector<Point> tempList;
//
//    if (sagitalFix.eval(mainContour->GetPoint(allPoints[3])) < 0)
//    {
//        for (int i = tSize - 1; i > 0; i--)
//        {
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            if (sagitalMove.eval(pnt) > 0)
//            {
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else
//            {
//                break;
//            }
//        }
//    }
//    else
//    {
//        for (int i = 1; i < tSize; i++)
//        {
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            if (sagitalMove.eval(pnt) > 0)
//            {
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else
//            {
//                break;
//            }
//        }
//    }
//
//    result = reduceSortPointsByAmount(tempList, amount + 2);
//
//    if (result.size() > amount)
//    {
//        return std::vector<Point>(result.begin() + 1, result.begin() + 1 + amount);
//    }
//    else
//    {
//        return result;
//    }
//}
//
//std::vector<Point> RegistrationPointsVTK::GetObliquePointsTibia(const Point& fixPoint, const Point& movePoint, int amount)
//{
//    std::vector<Point> result;
//
//    Plane coronalPlane, transverseFix, transverseMove;
//    coronalPlane.init(mVectorTibiaAPFront, (mKnee.getLateralPlateau() + mKnee.getMedialPlateau()) / 2.0);
//    transverseFix.init(mVectorTibiaAxisUp, fixPoint);
//    transverseFix.reverseByPoint(movePoint);
//
//    transverseMove.init(mVectorTibiaAxisUp, movePoint);
//    transverseMove.reverseByPoint(fixPoint);
//
//    Plane oblique = coronalPlane.getPerpendicularPlane(fixPoint, movePoint);
//
//    vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), oblique.getNormalVector(), oblique.getPoint());
//
//    vtkIdType nearA;
//
//    nearA = ImplantTools::GetNearestPoints(mainContour, fixPoint);
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(mainContour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearA)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearA && pos == -1)
//    {
//        pos = cont - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return result;
//    }
//
//    std::vector<Point> tempList;
//
//    if (transverseFix.eval(mainContour->GetPoint(allPoints[3])) < 0)
//    {
//        for (int i = tSize - 1; i > 0; i--)
//        {
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            if (transverseMove.eval(pnt) > 0)
//            {
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else
//            {
//                break;
//            }
//        }
//    }
//    else
//    {
//        for (int i = 1; i < tSize; i++)
//        {
//            double pnt[3];
//            mainContour->GetPoint(allPoints[i], pnt);
//            if (transverseMove.eval(pnt) > 0)
//            {
//                tempList.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else
//            {
//                break;
//            }
//        }
//    }
//
//    result = reduceSortPointsByAmount(tempList, amount + 2);
//
//    if (result.size() > amount)
//    {
//        return std::vector<Point>(result.begin() + 1, result.begin() + 1 + amount);
//    }
//    else
//    {
//        return result;
//    }
//}
//
///*
//std::pair<std::vector<Point>, std::vector<Point>> RegistrationPointsVTK::GetTransversePointsTibia(int amount)
//{
//    std::pair<std::vector<Point>, std::vector<Point>> resultPair;
//
//    if (amount < 2)
//    {
//        return resultPair;
//    }
//
//    std::vector<Point> resultMed, resultLat;
//    Plane transverse, planeTEAMed, planeTEALat, planeAP;
//    transverse.init(mKnee.getNormalVectorTibiaPlane(), mKnee.getTibiaTubercle());
//
//    Point myVector = mKnee.getMedialPlateau() - transverse.getProjectionPoint(mKnee.getMedialPlateau());
//
//    Point transversePnt = mKnee.getTibiaTubercle() + (myVector) * 0.7;
//
//    transverse.movePlane(transversePnt);
//
//    Point tibiaVectorTEA = mKnee.getNormalVectorTibiaPlane().cross(mVectorTibiaAPFront);
//
//    vtkSmartPointer<vtkPolyData> contour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), transverse.getNormalVector(), transverse.getPoint());
//
//    planeTEAMed.init(mVectorTibiaAPFront, mKnee.getMedialPlateau());
//    planeTEAMed.reverseByPoint(mKnee.getTibiaTubercle());
//
//    planeTEALat.init(mVectorTibiaAPFront, mKnee.getLateralPlateau());
//    planeTEALat.reverseByPoint(mKnee.getTibiaTubercle());
//
//    if (planeTEAMed.eval(mKnee.getTibiaCenterPointFreeAP()) > 0)
//    {
//        planeTEAMed.movePlane(mKnee.getTibiaCenterPointFreeAP());
//    }
//
//    if (planeTEALat.eval(mKnee.getTibiaCenterPointFreeAP()) > 0)
//    {
//        planeTEALat.movePlane(mKnee.getTibiaCenterPointFreeAP());
//    }
//
//    planeAP.init(tibiaVectorTEA, mKnee.getLateralPlateau());
//    planeAP.reverseByPoint(mKnee.getMedialPlateau());
//    planeAP.movePlane(mKnee.getTibiaKneeCenter());
//
//    double medialDistance = 99999999.0, lateralDistance = 99999999.0;
//    double medialDistTemp, letaralDistTemp;
//    vtkIdType nearA, nearB;
//
//    for (int i = 0; i < contour->GetPoints()->GetNumberOfPoints(); i++)
//    {
//        double pnt[3];
//        contour->GetPoint(i, pnt);
//
//        if (planeAP.eval(pnt) > 0)
//        {
//            medialDistTemp = abs(planeTEAMed.eval(pnt));
//            if (medialDistance > medialDistTemp)
//            {
//                medialDistance = medialDistTemp;
//                nearA = i;
//            }
//        }
//        else
//        {
//            letaralDistTemp = abs(planeTEALat.eval(pnt));
//            if (lateralDistance > letaralDistTemp)
//            {
//                lateralDistance = letaralDistTemp;
//                nearB = i;
//            }
//        }
//    }
//
//    std::list<std::pair<vtkIdType, vtkIdType>> lines;
//    ImplantTools::ExtractSortLines(contour, lines);
//
//    std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//    it1 = lines.begin();
//    it2 = lines.end();
//
//    std::vector<vtkIdType> allPoints;
//    int cont = 0;
//    int pos = -1;
//
//    for (; it1 != it2; ++it1)
//    {
//        allPoints.push_back(it1->first);
//        if (it1->first == nearA)
//        {
//            pos = cont;
//        }
//        cont++;
//    }
//
//    int tSize = cont;
//
//    if (lines.back().second == nearA && pos == -1)
//    {
//        pos = cont - 1;
//    }
//
//    if (pos > 0)
//    {
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//    }
//    else if (pos == -1 || tSize < 20)
//    {
//        return resultPair;
//    }
//
//    std::vector<Point> tempListMedial, tempListLateral;
//
//    if (planeTEAMed.eval(contour->GetPoint(allPoints[5])) < 0)
//    {
//        for (int i = tSize - 1; i > 0; i--)
//        {
//            if (allPoints[i] == nearB)
//            {
//                break;
//            }
//            double pnt[3];
//            contour->GetPoint(allPoints[i], pnt);
//            if (planeAP.eval(pnt) > 5)
//            {
//                tempListMedial.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else if (planeAP.eval(pnt) < -5)
//            {
//                tempListLateral.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//    }
//    else
//    {
//        for (int i = 1; i < tSize; i++)
//        {
//            if (allPoints[i] == nearB)
//            {
//                break;
//            }
//            double pnt[3];
//            contour->GetPoint(allPoints[i], pnt);
//            if (planeAP.eval(pnt) > 5)
//            {
//                tempListMedial.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//            else if (planeAP.eval(pnt) < -5)
//            {
//                tempListLateral.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//    }
//
//    std::reverse(tempListMedial.begin(), tempListMedial.end());
//    resultMed = reduceSortPointsByAmount(tempListMedial, amount + 1);
//    //std::reverse(tempListMedial.begin(), tempListMedial.end());
//
//    //std::reverse(tempListLateral.begin(), tempListLateral.end());
//    resultLat = reduceSortPointsByAmount(tempListLateral, amount + 1);
//
//
//    if (resultMed.size() < amount || resultLat.size() < amount)
//    {
//        goodTibiaPoints = false;
//        return resultPair;
//    }
//
//    if (resultMed.size() > amount)
//    {
//        resultPair.second = std::vector<Point>(resultMed.begin(), resultMed.begin() + amount);
//    }
//    else
//    {
//        resultPair.second = resultMed;
//    }
//
//    if (resultLat.size() > amount)
//    {
//        resultPair.first = std::vector<Point>(resultLat.begin(), resultLat.begin() + amount);
//    }
//    else
//    {
//        resultPair.first = resultLat;
//    }
//
//    lateralTibiaPoint1 = resultPair.first[resultPair.first.size() - 1];
//    medialTibiaPoint1 = resultPair.second[resultPair.second.size() - 1];
//
//    lateralTibiaPoint2 = resultPair.first[resultPair.first.size() - 2];
//    medialTibiaPoint2 = resultPair.second[resultPair.second.size() - 2];
//
//    goodTibiaPoints = true;
//
//    return resultPair;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetUpPointsTibia(const std::vector<Point>& pPoints, const Point& pPlateau)
//{
//    std::vector<Point> result, resultFinal;
//
//    Point axisVector = mKnee.getNormalVectorTibiaPlane();
//    axisVector.normalice();
//
//    Plane transverse;
//    transverse.init(axisVector, mKnee.getTibiaKneeCenter());
//
//    std::vector<vtkSmartPointer<vtkPolyData>> allContours;
//
//    for (int i = 0; i < pPoints.size(); i++)
//    {
//        Line lineAxis(axisVector, pPoints[i]);
//        Line tempLine = lineAxis.getPerpendicularLine(pPlateau);
//
//        Plane tempPlane;
//        tempPlane.init(tempLine.getDirectVector(), pPoints[i]);
//        Point intercep = tempPlane.getInterceptionLinePoint(tempLine);
//
//        Point perpVector = tempPlane.getProjectionPoint(pPlateau) - pPlateau;
//        perpVector.normalice();
//        Point a = intercep + 1000.0 * axisVector;
//        Point b = intercep + 1000.0 * perpVector;
//
//        Line obliqueLine = Line::makeLineWithPoints(a, b);
//
//        Plane cutPlane = transverse.getPerpendicularPlane(pPoints[i], pPlateau);
//
//        vtkSmartPointer<vtkPolyData> mainContour = ImplantTools::getMaxContour(mKnee.GetTibiaPoly(), cutPlane.getNormalVector(), cutPlane.getPoint());
//        allContours.push_back(mainContour);
//
//        Sphere mySphere;
//        Point nearPoint = getNearPointToLine(obliqueLine, mainContour->GetPoints(), mySphere);
//
//        result.push_back(nearPoint);
//    }
//
//    double radius = 9999999.0;
//
//    for (int i = 0; i < result.size(); i++)
//    {
//        double temp = ImplantTools::getDistanceBetweenPoints(result[i], pPlateau);
//        if (radius > temp)
//        {
//            radius = temp;
//        }
//    }
//
//    radius = radius * 0.9;
//
//    Sphere centerSphere(pPlateau, radius);
//
//    for (int i = 0; i < allContours.size(); i++)
//    {
//        vtkSmartPointer<vtkPolyData> contour = allContours[i];
//        vtkIdType nearA, nearB;
//        nearA = ImplantTools::GetNearestPoints(contour, pPlateau);
//        nearB = ImplantTools::GetNearestPoints(contour, pPoints[i]);
//
//        std::list<std::pair<vtkIdType, vtkIdType>> lines;
//        ImplantTools::ExtractSortLines(contour, lines);
//
//        std::list<std::pair<vtkIdType, vtkIdType>>::iterator it1, it2;
//        it1 = lines.begin();
//        it2 = lines.end();
//
//        std::vector<vtkIdType> allPoints;
//
//        int pos = 0;
//        int cont = 0;
//
//        for (; it1 != it2; ++it1)
//        {
//            allPoints.push_back(it1->first);
//            if (it1->first == nearA)
//            {
//                pos = cont;
//            }
//            cont++;
//        }
//
//        int tSize = cont;
//
//        if (lines.back().second == nearA && pos == -1)
//        {
//            pos = cont - 1;
//        }
//
//        std::rotate(allPoints.begin(), allPoints.begin() + pos, allPoints.end());
//
//        Plane decidePlane;
//        decidePlane.init(mVectorTibiaAPFront, pPlateau);
//        decidePlane.reverseByPoint(pPoints[i]);
//
//        std::vector<Point> tempListPoints;
//
//        if (decidePlane.eval(contour->GetPoint(allPoints[5])) < 0)
//        {
//            for (int j = tSize - 1; j > 0; j--)
//            {
//                if (allPoints[j] == nearB)
//                {
//                    break;
//                }
//                double pnt[3];
//                contour->GetPoint(allPoints[j], pnt);
//                tempListPoints.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//        else
//        {
//            for (int j = 1; j < tSize; j++)
//            {
//                if (allPoints[j] == nearB)
//                {
//                    break;
//                }
//                double pnt[3];
//                contour->GetPoint(allPoints[j], pnt);
//                tempListPoints.push_back(Point(pnt[0], pnt[1], pnt[2]));
//            }
//        }
//
//        for (int j = 0; j < tempListPoints.size() - 1; j++)
//        {
//            if ((centerSphere.isPointInside(tempListPoints[j]) == false && centerSphere.isPointInside(tempListPoints[j + 1]) == true) ||
//                (centerSphere.isPointInside(tempListPoints[j]) == true && centerSphere.isPointInside(tempListPoints[j + 1]) == false))
//            {
//                Line myLine = Line::makeLineWithPoints(tempListPoints[j], tempListPoints[j + 1]);
//                std::pair<Point, Point> intercep;
//
//                myLine.getInterceptionSphere(centerSphere.center, centerSphere.radius, intercep);
//
//                Point outSise;
//                if (centerSphere.isPointInside(tempListPoints[j]) == false)
//                {
//                    outSise = tempListPoints[j];
//                }
//                else
//                {
//                    outSise = tempListPoints[j + 1];
//                }
//
//                if (ImplantTools::getDistanceBetweenPoints(outSise, intercep.first, true) < ImplantTools::getDistanceBetweenPoints(outSise, intercep.second, true))
//                {
//                    resultFinal.push_back(intercep.first);
//                }
//                else
//                {
//                    resultFinal.push_back(intercep.second);
//                }
//            }
//
//        }
//    }
//
//    return resultFinal;
//}
//
//std::pair<std::vector<Point>, std::vector<Point>> RegistrationPointsVTK::GetUpsPointsTibia(int amount)
//{
//    std::pair<std::vector<Point>, std::vector<Point>> mainPoints = GetTransversePointsTibia(amount);
//
//    std::vector<Point> medPoints, latPoints;
//
//    latPoints = GetUpPointsTibia(mainPoints.first, mKnee.getLateralPlateau());
//
//    medPoints = GetUpPointsTibia(mainPoints.second, mKnee.getMedialPlateau());
//
//    std::pair<std::vector<Point>, std::vector<Point>> result;
//    result.first = latPoints;
//    result.second = medPoints;
//
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetMedialObliqueOut(int amount)
//{
//    std::vector<Point> result;
//    if (goodTibiaPoints == false)
//    {
//        return result;
//    }
//
//    result = GetObliquePointsTibia(medialTibiaPoint2, mKnee.getTibiaTubercle(), amount);
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetMedialObliqueIn(int amount)
//{
//    Point movePoint = mKnee.getTibiaTubercle() + (mKnee.getTibiaKneeCenter() - mKnee.getTibiaTubercle()) * 0.3;
//
//    std::vector<Point> result;
//    if (goodTibiaPoints == false)
//    {
//        return result;
//    }
//
//    result = GetObliquePointsTibia(medialTibiaPoint1, movePoint, amount);
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetLateralObliqueOut(int amount)
//{
//    std::vector<Point> result;
//    if (goodTibiaPoints == false)
//    {
//        return result;
//    }
//
//    result = GetObliquePointsTibia(lateralTibiaPoint2, mKnee.getTibiaTubercle(), amount);
//    return result;
//}
//
//std::vector<Point> RegistrationPointsVTK::GetLateralObliqueIn(int amount)
//{
//    Point movePoint = mKnee.getTibiaTubercle() + (mKnee.getTibiaKneeCenter() - mKnee.getTibiaTubercle()) * 0.3;
//
//    std::vector<Point> result;
//    if (goodTibiaPoints == false)
//    {
//        return result;
//    }
//
//    result = GetObliquePointsTibia(lateralTibiaPoint1, movePoint, amount);
//    return result;
//}
//*/