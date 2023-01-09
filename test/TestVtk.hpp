#include <vtkNew.h>
#include <fstream>
#include <opencv2/calib3d/calib3d.hpp>
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "vtkFlyingEdges3D.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkNamedColors.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
#include "itkImageToVTKImageFilter.h"
#include "itkVTKImageToImageFilter.h"
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkProperty.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
#include "vtkCellData.h"
#include <vtkSTLReader.h>

#include "vtkClipPolyData.h"
#include "vtkInteractorObserver.h"
#include "vtkInteractorStyleSwitch.h"
#include "vtkPolyDataMapper.h"
#include "vtkRenderer.h"
#include "vtkSmartPointer.h"
#include "vtkSphereSource.h"
#include <vtkStripper.h>
#include "vtkFeatureEdges.h"
#include "vtkClipClosedSurface.h"
#include "vtkAppendPolyData.h"
#include "vtkPlaneCollection.h"
#include "vtkImplicitPolyDataDistance.h"
#include "vtkCellLocator.h"
#include "vtkUnsignedCharArray.h"

#include "TestFunction.hpp"
#include "itkImage.h"
#include "vtkPolyLine.h"
#include "vtkLine.h"
#include "vtkPolyDataReader.h"
#include "vtkCleanPolyData.h"
#include "vtkImageData.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkImageStencil.h"
#include "vtkTransformPolyDataFilter.h"
#include "itkImage.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "itkCastImageFilter.h"
#include "vtkXMLPolyDataReader.h"

#include <vtkCellArray.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkNamedColors.h>
#include <vtkObjectFactory.h>
#include <vtkPlaneSource.h>
#include <vtkPropPicker.h>
#include <vtkProperty.h>
#include <vtkRendererCollection.h>

using RegistrationPixelType = int16_t;
using RegistrationImageType = itk::Image<RegistrationPixelType, 3>;


namespace TestVTK
{
    using SegmentPixelType = uint8_t;
    using SegmentImageType = itk::Image<SegmentPixelType, 3>;
    using PointTypeITK = itk::Point<double, 3>;
    double medialPoint[3] = { 42.46, -102.62, -632.71 };
    double medialPos[3] = { 75.0, -102.62, -632.71 };

    double lateralPoint[3] = { 115.36, -86.94, -632.71 };
    double lateralPos[3] = { 120.0, -86.94, -632.71 };

    using RegistrationPixelType = int16_t;
    using RegistrationImageType = itk::Image<RegistrationPixelType, 3>;


    std::string femurStr = "D:\\3D_DICOM\\right_femur_segment.nrrd";


    // Handle mouse events
    class MouseInteractorStyle2 : public vtkInteractorStyleTrackballCamera
    {
    public:
        static MouseInteractorStyle2* New();
        vtkTypeMacro(MouseInteractorStyle2, vtkInteractorStyleTrackballCamera);
        vtkNew<vtkNamedColors> colors;

        virtual void OnRightButtonDown() override
        {
            int* clickPos = this->GetInteractor()->GetEventPosition();

            // Pick from this location.
            vtkNew<vtkPropPicker> picker;
            picker->Pick(clickPos[0], clickPos[1], 0, this->GetDefaultRenderer());

            double* pos = picker->GetPickPosition();
            std::cout << "Pick position (world coordinates) is: " << pos[0] << " "
                << pos[1] << " " << pos[2] << std::endl;

            auto pickedActor = picker->GetActor();
            if (pickedActor == nullptr)
            {
                std::cout << "No actor picked." << std::endl;
            }
            else
            {
                std::cout << "Picked actor: " << picker->GetActor() << std::endl;
                // Create a sphere
                // Create a sphere
                vtkNew<vtkSphereSource> sphereSource;
                sphereSource->SetCenter(pos[0], pos[1], pos[2]);
                sphereSource->SetRadius(1);

                // Create a mapper and actor
                vtkNew<vtkPolyDataMapper> mapper;
                mapper->SetInputConnection(sphereSource->GetOutputPort());

                vtkNew<vtkActor> actor;
                actor->SetMapper(mapper);
                actor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

                this->GetDefaultRenderer()->AddActor(actor);
                // Forward events
                vtkInteractorStyleTrackballCamera::OnRightButtonDown();
            }
        }

    private:
    };

    vtkStandardNewMacro(MouseInteractorStyle2);

    void picking(const vtkSmartPointer<vtkPolyData>& polydata)
    {
        vtkNew<vtkNamedColors> colors;
        // Create a mapper
        vtkNew<vtkPolyDataMapper> mapper;
        mapper->SetInputData(polydata);
        mapper->ScalarVisibilityOff();

        // Create an actor
        vtkNew<vtkActor> actor;
        actor->SetMapper(mapper);
        actor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

        std::cout << "Actor address: " << actor << std::endl;

        // A renderer and render window
        vtkNew<vtkRenderer> renderer;
        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->AddRenderer(renderer);
        renderWindow->SetWindowName("Picking");

        // An interactor
        vtkNew<vtkRenderWindowInteractor> renderWindowInteractor;
        renderWindowInteractor->SetRenderWindow(renderWindow);

        // Set the custom stype to use for interaction.
        vtkNew<MouseInteractorStyle2> style;
        style->SetDefaultRenderer(renderer);

        renderWindowInteractor->SetInteractorStyle(style);

        // Add the actors to the scene
        renderer->AddActor(actor);
        renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

        // Render and interact
        renderWindow->Render();
        renderWindowInteractor->Initialize();
        renderWindowInteractor->Start();
    }

    void show(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2)
    {
        vtkNew<vtkNamedColors> colors;

        vtkNew<vtkPolyDataMapper> contoursMapper;
        contoursMapper->SetInputData(poly1);
        contoursMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> contoursActor;
        contoursActor->SetMapper(contoursMapper);
        contoursActor->GetProperty()->SetRepresentationToWireframe();
        contoursActor->GetProperty()->ShadingOff();
        contoursActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

        vtkNew<vtkPolyDataMapper> surfaceMapper;
        surfaceMapper->SetInputData(poly2);
        surfaceMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> surfaceActor;
        surfaceActor->SetMapper(surfaceMapper);
        surfaceActor->GetProperty()->SetRepresentationToWireframe();
        surfaceActor->GetProperty()->ShadingOff();
        surfaceActor->GetProperty()->SetColor(
            colors->GetColor3d("MistyRose").GetData());

        // Create two renderers side by side to show the contours and the surface
        // separately
        //
        std::cout << "Press 't' for trackball interaction" << std::endl;
        std::cout << "Press 'r' to reset the camera" << std::endl;
        std::cout << "Press 'w' for wireframe representation" << std::endl;
        std::cout << "Press 's' for surface representation" << std::endl;

        vtkNew<vtkRenderer> renderer1;
        renderer1->SetViewport(0., 0., 0.5, 1.);
        renderer1->SetBackground(colors->GetColor3d("CadetBlue").GetData());

        vtkNew<vtkRenderer> renderer2;
        renderer2->SetViewport(0.5, 0., 1., 1.);
        renderer2->SetBackground(colors->GetColor3d("BurlyWood").GetData());
        //renderer2->SetBackground(colors->GetColor3d("SlateGray").GetData());

        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->SetSize(800, 400);
        renderWindow->SetWindowName("ContoursToSurface");

        renderWindow->AddRenderer(renderer1);
        renderWindow->AddRenderer(renderer2);

        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(renderWindow);

        renderer1->AddViewProp(contoursActor);
        renderer2->AddViewProp(surfaceActor);
        renderWindow->Render();

        interactor->Start();
    }

    void show(vtkSmartPointer<vtkPolyData> poly, const std::vector<cv::Point3d>& points, bool makePolyLine = false)
    {
        vtkNew<vtkNamedColors> colors;

        vtkNew<vtkPolyDataMapper> contoursMapper;
        contoursMapper->SetInputData(poly);
        contoursMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> polyActor;
        polyActor->SetMapper(contoursMapper);
        polyActor->GetProperty()->SetRepresentationToWireframe();
        polyActor->GetProperty()->ShadingOff();
        polyActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

        std::vector<vtkSmartPointer<vtkActor>> pointsActor;

        if (makePolyLine == true)
        {
            vtkNew<vtkPoints> tPoints;
            for (int i = 0; i < points.size(); i++)
            {
                double pnt[3];
                pnt[0] = points[i].x;
                pnt[1] = points[i].y;
                pnt[2] = points[i].z;

                tPoints->InsertNextPoint(pnt);
            }

            vtkNew<vtkCellArray> cells;

            for (unsigned int i = 0; i < tPoints->GetNumberOfPoints() - 1; i++)
            {
                vtkNew<vtkLine> myLine;

                myLine->GetPointIds()->SetId(0, i);
                myLine->GetPointIds()->SetId(1, i + 1);

                cells->InsertNextCell(myLine);
            }

            vtkNew<vtkPolyData> polyLine;
            polyLine->SetPoints(tPoints);
            polyLine->SetLines(cells);

            vtkNew<vtkPolyDataMapper> polyMapper;
            polyMapper->SetInputData(polyLine);
            polyMapper->ScalarVisibilityOff();

            vtkNew<vtkActor> polyActor;
            polyActor->SetMapper(polyMapper);
            polyActor->GetProperty()->SetRepresentationToWireframe();
            polyActor->GetProperty()->ShadingOff();
            polyActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

            pointsActor.push_back(polyActor);
        }
        else
        {
            for (int i = 0; i < points.size(); i++)
            {
                double pnt[3];
                pnt[0] = points[i].x;
                pnt[1] = points[i].y;
                pnt[2] = points[i].z;

                vtkNew<vtkSphereSource> sphere;
                sphere->SetCenter(pnt);
                sphere->SetRadius(1);
                sphere->Update();

                vtkNew<vtkPolyDataMapper> sphereMapper;
                sphereMapper->SetInputData(sphere->GetOutput());
                sphereMapper->ScalarVisibilityOff();

                vtkNew<vtkActor> sphereActor;
                sphereActor->SetMapper(sphereMapper);
                sphereActor->GetProperty()->SetRepresentationToWireframe();
                sphereActor->GetProperty()->ShadingOff();
                sphereActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

                pointsActor.push_back(sphereActor);
            }
        }

        vtkNew<vtkRenderer> renderer;
        //renderer->SetViewport(0., 0., 0.5, 1.);
        renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->SetSize(800, 400);
        renderWindow->SetWindowName("Surface");

        renderWindow->AddRenderer(renderer);

        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(renderWindow);

        renderer->AddActor(polyActor);

        for (int i = 0; i < pointsActor.size(); i++)
        {
            renderer->AddActor(pointsActor[i]);
        }

        renderWindow->Render();

        interactor->Start();
    }

    void show_join(vtkSmartPointer<vtkPolyData> poly1, vtkSmartPointer<vtkPolyData> poly2)
    {
        vtkNew<vtkNamedColors> colors;

        vtkNew<vtkPolyDataMapper> contoursMapper;
        contoursMapper->SetInputData(poly1);
        contoursMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> polyActor;
        polyActor->SetMapper(contoursMapper);
        polyActor->GetProperty()->SetRepresentationToWireframe();
        polyActor->GetProperty()->ShadingOff();
        polyActor->GetProperty()->SetColor(colors->GetColor3d("red").GetData());

        std::vector<vtkSmartPointer<vtkActor>> pointsActor;

        vtkNew<vtkPolyDataMapper> sphereMapper;
        sphereMapper->SetInputData(poly2);
        sphereMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> sphereActor;
        sphereActor->SetMapper(sphereMapper);
        sphereActor->GetProperty()->SetRepresentationToWireframe();
        sphereActor->GetProperty()->ShadingOff();
        sphereActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

        pointsActor.push_back(sphereActor);

        vtkNew<vtkRenderer> renderer;
        //renderer->SetViewport(0., 0., 0.5, 1.);
        renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->SetSize(800, 400);
        renderWindow->SetWindowName("Surface");

        renderWindow->AddRenderer(renderer);

        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(renderWindow);

        renderer->AddActor(polyActor);

        for (int i = 0; i < pointsActor.size(); i++)
        {
            renderer->AddActor(pointsActor[i]);
        }

        renderWindow->Render();

        interactor->Start();
    }

    void show(vtkSmartPointer<vtkPolyData> poly, std::vector<vtkSmartPointer<vtkPolyData>> polyList = {})
    {
        /*vtkNew<vtkNamedColors> colors;

        vtkNew<vtkPolyDataMapper> contoursMapper;
        contoursMapper->SetInputData(poly);
        contoursMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> polyActor;
        polyActor->SetMapper(contoursMapper);
        polyActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

        vtkNew<vtkRenderer> renderer;
        renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->SetSize(800, 400);
        renderWindow->SetWindowName("Surface");

        renderWindow->AddRenderer(renderer);

        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(renderWindow);

        renderer->AddActor(polyActor);

        renderWindow->Render();

        interactor->Start();*/

        vtkNew<vtkNamedColors> colors;

        vtkNew<vtkPolyDataMapper> contoursMapper;
        contoursMapper->SetInputData(poly);
        contoursMapper->ScalarVisibilityOff();

        vtkNew<vtkActor> polyActor;
        polyActor->SetMapper(contoursMapper);
        polyActor->GetProperty()->SetRepresentationToWireframe();
        polyActor->GetProperty()->ShadingOff();
        polyActor->GetProperty()->SetColor(colors->GetColor3d("MistyRose").GetData());

        std::vector<vtkSmartPointer<vtkActor>> pointsActor;

        for (int i = 0; i < polyList.size(); i++)
        {
            vtkNew<vtkPolyDataMapper> polyMapper;
            polyMapper->SetInputData(polyList[i]);
            polyMapper->ScalarVisibilityOff();

            vtkNew<vtkActor> sphereActor;
            sphereActor->SetMapper(polyMapper);
            sphereActor->GetProperty()->SetRepresentationToWireframe();
            sphereActor->GetProperty()->ShadingOff();
            sphereActor->GetProperty()->SetColor(colors->GetColor3d("blue").GetData());

            pointsActor.push_back(sphereActor);
        }

        vtkNew<vtkRenderer> renderer;
        //renderer->SetViewport(0., 0., 0.5, 1.);
        renderer->SetBackground(colors->GetColor3d("CadetBlue").GetData());

        vtkNew<vtkRenderWindow> renderWindow;
        renderWindow->SetSize(800, 400);
        renderWindow->SetWindowName("Surface");

        renderWindow->AddRenderer(renderer);

        vtkNew<vtkRenderWindowInteractor> interactor;
        interactor->SetRenderWindow(renderWindow);

        renderer->AddActor(polyActor);

        for (int i = 0; i < pointsActor.size(); i++)
        {
            renderer->AddActor(pointsActor[i]);
        }

        renderWindow->Render();

        interactor->Start();
    }


    vtkSmartPointer<vtkPolyData> itkImageToSurface3D(const std::string& pathImage)
    {
        SegmentImageType::Pointer imgOut;
        Test::readImage<SegmentImageType>(pathImage, imgOut);

        using FilterType = itk::ImageToVTKImageFilter<SegmentImageType>;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(imgOut);

        try
        {
            filter->Update();
        }
        catch (itk::ExceptionObject & error)
        {
            std::cout << "Error in conversion itk to vtk image: " << error << std::endl;
        }

        auto surface = vtkSmartPointer<vtkFlyingEdges3D>::New();

        surface->SetInputData(filter->GetOutput());

        surface->SetNumberOfContours(1);
        surface->SetValue(0, 1);
        surface->Update();

        vtkSmartPointer<vtkPolyData> poly = surface->GetOutput();
        return poly;
    }

    vtkSmartPointer<vtkPolyData> ReadPolyData(std::string name)
    {
        vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
        reader->SetFileName(name.c_str());
        reader->Update();
        return reader->GetOutput();
    }

    vtkSmartPointer<vtkPolyData> ReadPolyDataSTL(std::string name)
    {
        vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
        reader->SetFileName(name.c_str());
        reader->Update();
        return reader->GetOutput();
    }

    vtkSmartPointer<vtkPolyData> ReadPolyDataVTP(std::string name)
    {
        vtkSmartPointer<vtkXMLPolyDataReader> reader = vtkSmartPointer<vtkXMLPolyDataReader>::New();
        reader->SetFileName(name.c_str());
        reader->Update();
        return reader->GetOutput();
    }

    void SavePolyData(vtkSmartPointer<vtkPolyData> poly1, std::string name)
    {
        vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
        writer->SetFileName(name.c_str());
        writer->SetInputData(poly1);
        writer->Write();
        std::cout << "Save image" << std::endl;
    }

    vtkSmartPointer<vtkPolyData> GetCenterSlice(vtkSmartPointer<vtkPolyData> polyData, int axis)
    {
        vtkNew<vtkPlane> plane;
        vtkNew<vtkCutter> cutter;

        if (axis == 1)
        {
            plane->SetNormal(1, 0, 0);
        }
        else if (axis == 2)
        {
            plane->SetNormal(0, 1, 0);
        }
        else
        {
            plane->SetNormal(0, 0, 1);
        }

        plane->SetOrigin(polyData->GetCenter());
        //plane->SetOrigin(99.0, -100.0, -612.0);
        cutter->SetCutFunction(plane);
        cutter->SetInputData(polyData);
        cutter->Update();

        return cutter->GetOutput();
    }


    void split()
    {
        double medialProj[3];
        double lateralProj[3];
        vtkSmartPointer<vtkPolyData> femur = ReadPolyData(femurStr);
        //vtkSmartPointer<vtkClipPolyData> clip = vtkSmartPointer<vtkClipPolyData>::New();
        //clip->SetValue(0);
        ////clip->SetInsideOut(true);
        //clip->GenerateClippedOutputOn();
        //clip->SetInputData(femur);
        //vtkSmartPointer<vtkPlane> plane = vtkSmartPointer<vtkPlane>::New();
        //plane->SetOrigin(32.46, -102.62, -632.71);
        //plane->SetNormal(-1.0, 0.0, 0.0);
        //plane->ProjectPoint(medialPoint, medialProj);
        //plane->ProjectPoint(lateralPoint, lateralProj);
        //double medialNormal[3] = { medialPoint[0] - medialProj[0], medialPoint[1] - medialProj[1], medialPoint[2] - medialProj[2] };
        //double lateralNormal[3] = { lateralPoint[0] - lateralProj[0], lateralPoint[1] - lateralProj[1], lateralPoint[2] - lateralProj[2] };

        //plane->SetNormal(lateralNormal);

        //clip->SetClipFunction(plane);
        //clip->Update();
        //auto surface = clip->GetOutput();
        /////////////////////////////////////////////


        ////////////////////////////////////////////

        vtkNew<vtkFeatureEdges> boundaryEdges;
        boundaryEdges->SetInputData(femur);
        boundaryEdges->BoundaryEdgesOn();
        boundaryEdges->FeatureEdgesOff();
        boundaryEdges->NonManifoldEdgesOff();
        boundaryEdges->ManifoldEdgesOff();
        boundaryEdges->Update();

        auto border = boundaryEdges->GetOutput();

        std::cout << "First numeber cell: " << boundaryEdges->GetOutput()->GetNumberOfCells() << std::endl;

        /* vtkNew<vtkStripper> boundaryStrips;
         boundaryStrips->SetInputData(boundaryEdges->GetOutput());
         boundaryStrips->Update();

         vtkNew<vtkPolyData> boundaryPoly;
         boundaryPoly->SetPoints(boundaryStrips->GetOutput()->GetPoints());
         boundaryPoly->SetPolys(boundaryStrips->GetOutput()->GetLines());*/

        vtkNew<vtkAppendPolyData> filter;

        filter->AddInputData(femur);
        filter->AddInputData(border);
        filter->Update();
        auto surface2 = filter->GetOutput();

        auto center = GetCenterSlice(surface2, 1);


        /////////////////////

        vtkNew<vtkFeatureEdges> boundaryEdges2;
        boundaryEdges2->SetInputData(surface2);
        boundaryEdges2->BoundaryEdgesOn();
        boundaryEdges2->FeatureEdgesOff();
        boundaryEdges2->NonManifoldEdgesOn();
        boundaryEdges2->ManifoldEdgesOff();
        boundaryEdges2->Update();

        std::cout << "Second numeber cell: " << boundaryEdges2->GetOutput()->GetNumberOfCells() << std::endl;

        /*vtkNew<vtkStripper> boundaryStrips2;
        boundaryStrips2->SetInputConnection(boundaryEdges2->GetOutputPort());
        boundaryStrips2->Update();
        vtkNew<vtkPolyData> boundaryPoly2;
        boundaryPoly2->SetPoints(boundaryStrips2->GetOutput()->GetPoints());
        boundaryPoly2->SetPolys(boundaryStrips2->GetOutput()->GetLines());*/

        ///////////////////////
        show(surface2, center);
    }


    void splitClose()
    {
        vtkSmartPointer<vtkPolyData> femur = ReadPolyData(femurStr);

        double medialNormal[3] = { medialPoint[0] - medialPos[0], medialPoint[1] - medialPos[1], medialPoint[2] - medialPos[2] };
        double lateralNormal[3] = { -lateralPoint[0] + lateralPos[0], -lateralPoint[1] + lateralPos[1], -lateralPoint[2] + lateralPos[2] };

        vtkNew<vtkPlane> plane1;
        plane1->SetOrigin(medialPos);
        plane1->SetNormal(lateralNormal);
        vtkNew<vtkPlane> plane2;
        plane2->SetOrigin(lateralPos);
        plane2->SetNormal(medialNormal);

        vtkNew<vtkPlaneCollection> planes;
        planes->AddItem(plane1);
        planes->AddItem(plane2);

        vtkNew<vtkClipClosedSurface> clipper;
        clipper->SetInputData(femur);
        clipper->SetClippingPlanes(planes);
        /*clipper->SetActivePlaneId(1);
        clipper->SetScalarModeToColors();
        clipper->SetClipColor(colors->GetColor3d("Banana").GetData());
        clipper->SetBaseColor(colors->GetColor3d("Tomato").GetData());
        clipper->SetActivePlaneColor(colors->GetColor3d("SandyBrown").GetData());*/
        clipper->Update();

        auto surface = clipper->GetOutput();

        auto center = GetCenterSlice(surface, 3);

        show(surface, center);
    }

    void calculateDistance()
    {
        vtkSmartPointer<vtkPolyData> femur = ReadPolyData(femurStr);

        double p1[3] = { 96.71, -60.715, -633.41 };

        vtkNew<vtkCellLocator> cellLocator;
        cellLocator->SetDataSet(femur);
        cellLocator->BuildLocator();

        // Find the closest points to TestPoint
        double closestPoint[3];   // the coordinates of the closest point will be
                                  // returned here
        double closestPointDist2; // the squared distance to the closest point will be
                                  // returned here
        vtkIdType cellId; // the cell id of the cell containing the closest point will
                          // be returned here
        int subId;        // this is rarely used (in triangle strips only, I believe)
        cellLocator->FindClosestPoint(p1, closestPoint, cellId, subId, closestPointDist2);

        cv::Point3d diff = cv::Point3d(closestPoint[0] - p1[0], closestPoint[1] - p1[1], closestPoint[2] - p1[2]);
        double myDistancePointToPoint = sqrt(diff.dot(diff));
        cv::Point3d resultNear = cv::Point3d(closestPoint[0], closestPoint[1], closestPoint[2]);

        std::cout << "Distance result: " << sqrt(closestPointDist2) << " distance to point: " << myDistancePointToPoint << std::endl;
        std::cout << "Point result: " << resultNear << std::endl;

        vtkNew<vtkIdList> idList;
        femur->GetCellPoints(cellId, idList);

        for (int i = 0; i < idList->GetNumberOfIds(); i++)
        {
            double pnt[3];
            femur->GetPoint(idList->GetId(i), pnt);
            cv::Point3d cvPoint = cv::Point3d(pnt[0], pnt[1], pnt[2]);
            std::cout << "Points: " << cvPoint << std::endl;
        }


    }

    vtkSmartPointer<vtkPolyData> createOne(double z)
    {
        double p0[3] = { 2.0, 3.0, z };
        double p1[3] = { 3.0, 3.0, z };
        double p2[3] = { 3.0, 4.0, z };
        double p3[3] = { 2.0, 4.0, z };



        // Create a vtkPoints object and store the points in it
        vtkNew<vtkPoints> points;
        points->InsertNextPoint(p0);
        points->InsertNextPoint(p1);
        points->InsertNextPoint(p2);
        points->InsertNextPoint(p3);

        vtkNew<vtkCellArray> cells;
        for (int i = 0; i < points->GetNumberOfPoints(); i++)
        {
            vtkSmartPointer<vtkLine> line = vtkSmartPointer<vtkLine>::New();
            if (i == points->GetNumberOfPoints() - 1)
            {
                line->GetPointIds()->SetId(0, i);
                line->GetPointIds()->SetId(1, 0);
                cells->InsertNextCell(line);
            }
            else
            {
                line->GetPointIds()->SetId(0, i);
                line->GetPointIds()->SetId(1, i + 1);
                cells->InsertNextCell(line);
            }
        }
        vtkNew<vtkPolyData> polyData;
        polyData->SetPoints(points);
        polyData->SetLines(cells);
        return polyData;
    }

    vtkSmartPointer<vtkPolyData> createSlices(int amount)
    {
        vtkNew<vtkAppendPolyData> filter;
        for (int i = 0; i < amount; i++)
        {
            filter->AddInputData(createOne(double(i)));
        }
        filter->Update();
        return filter->GetOutput();
    }

    vtkSmartPointer<vtkPolyData> CreatePolyLine(vtkSmartPointer<vtkPoints> points)
    {
        vtkNew<vtkCellArray> cells;

        for (unsigned int i = 0; i < points->GetNumberOfPoints() - 1; i++)
        {
            vtkNew<vtkLine> myLine;

            myLine->GetPointIds()->SetId(0, i);
            myLine->GetPointIds()->SetId(1, i + 1);

            cells->InsertNextCell(myLine);
        }

        vtkNew<vtkPolyData> polyData;
        polyData->SetPoints(points);
        polyData->SetLines(cells);
        return polyData;
    }

    vtkSmartPointer<vtkPolyData> CreateClosePolyLine(vtkSmartPointer<vtkPoints> points)
    {
        vtkNew<vtkCellArray> cells;

        for (unsigned int i = 0; i < points->GetNumberOfPoints(); i++)
        {
            cells->InsertNextCell(2);
            cells->InsertCellPoint(i);
            cells->InsertCellPoint((i + 1) == points->GetNumberOfPoints() ? 0: i + 1);
        }

        vtkNew<vtkPolyData> polyData;
        polyData->SetLines(cells);
        polyData->SetPoints(points);
        return polyData;
    }

    vtkSmartPointer<vtkPolyData> CreatePolyLine(std::vector<PointTypeITK> pPoints)
    {
        vtkNew<vtkPoints> points;
        for (int i = 0; i < pPoints.size(); i++)
        {
            double pnt[3];
            pnt[0] = pPoints[i][0];
            pnt[1] = pPoints[i][1];
            pnt[2] = pPoints[i][2];

            points->InsertNextPoint(pnt);
        }

        return CreatePolyLine(points);
    }

    void ChangePolyColor(vtkSmartPointer<vtkPolyData>& poly, double rgb[3])
    {
        vtkNew<vtkUnsignedCharArray> cellData;
        cellData->SetNumberOfComponents(3);
        cellData->SetNumberOfTuples(poly->GetNumberOfCells());

        for (int i = 0; i < poly->GetNumberOfCells(); i++)
        {
            //double rgb[3] = {255, 0, 0};

            cellData->InsertTuple(i, rgb);
        }

        poly->GetCellData()->SetScalars(cellData);
    }

    void MakeSlicesX(const vtkSmartPointer<vtkPolyData> polyData, double distance)
    {
        double bounds[6];
        polyData->GetBounds(bounds);
        double min, max;
        vtkNew<vtkPlane> plane;
        vtkNew<vtkCutter> cutter;
        cutter->SetInputData(polyData);
        std::vector<vtkSmartPointer<vtkPolyData>> allSurface;
        vtkNew<vtkAppendPolyData> filter;
        if (distance <= 1)
        {
            distance = 1.0;
        }

        min = bounds[0];
        max = bounds[1];

        plane->SetNormal(1, 0, 0);

        for (double i = min + 0.2; i < max; i += distance)
        {
            plane->SetOrigin(i, bounds[2], bounds[4]);
            cutter->SetCutFunction(plane);
            cutter->Update();

            auto contour = cutter->GetOutput();

            if (contour->GetNumberOfPoints() > 0)
            {
                vtkNew<vtkPolyDataConnectivityFilter> connectivityFilter;
                connectivityFilter->SetInputData(contour);
                connectivityFilter->SetExtractionModeToLargestRegion();
                connectivityFilter->Update();
                filter->AddInputData(connectivityFilter->GetOutput());
            }
        }
        filter->Update();
        show(filter->GetOutput(), polyData);
    }

    vtkSmartPointer<vtkPolyData> getContours(const vtkSmartPointer<vtkPolyData> polyData, const cv::Point3d& pNormal, const cv::Point3d& pPoint)
    {
        vtkNew<vtkPlane> plane;
        vtkNew<vtkCutter> cutter;
        cutter->SetInputData(polyData);

        plane->SetNormal(pNormal.x, pNormal.y, pNormal.z);
        plane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
        cutter->SetCutFunction(plane);
        cutter->Update();

        auto contour = cutter->GetOutput();

        return contour;
    }

    vtkSmartPointer<vtkPolyData> CreateSphereTest(const cv::Point3d& pPoint)
    {
        double pnt[3];
        pnt[0] = pPoint.x;
        pnt[1] = pPoint.y;
        pnt[2] = pPoint.z;

        vtkNew<vtkSphereSource> sphere;
        sphere->SetCenter(pnt);
        sphere->SetRadius(10);
        sphere->Update();

        return sphere->GetOutput();
    }

    vtkSmartPointer<vtkPolyData> CleanPoly(const vtkSmartPointer<vtkPolyData> poly)
    {
        vtkNew<vtkPolyDataConnectivityFilter> femurConnectivityFilter;
        femurConnectivityFilter->SetInputData(poly);
        femurConnectivityFilter->SetExtractionModeToLargestRegion();
        femurConnectivityFilter->Update();

        vtkNew<vtkCleanPolyData> clean;
        clean->SetInputData(femurConnectivityFilter->GetOutput());
        clean->Update();
        return clean->GetOutput();
    }

    vtkSmartPointer<vtkPolyData> MergePolyWithSphere(const vtkSmartPointer<vtkPolyData> poly, const std::vector<cv::Point3d>& pPoints)
    {
        vtkNew<vtkAppendPolyData> filter;
        filter->AddInputData(poly);

        for (int i = 0; i < pPoints.size(); i++)
        {
            vtkSmartPointer<vtkPolyData> mySphere = TestVTK::CreateSphereTest(pPoints[i]);
            filter->AddInputData(mySphere);
        }

        filter->Update();
        return filter->GetOutput();
    }

    vtkSmartPointer<vtkPolyData> CreateSeveralSpheres(const std::vector<PointTypeITK>& pPoints)
    {
        vtkNew<vtkAppendPolyData> filter;

        for (int i = 0; i < pPoints.size(); i++)
        {
            cv::Point3d myPoint = { pPoints[i][0] , pPoints[i][1] , pPoints[i][2] };
            vtkSmartPointer<vtkPolyData> mySphere = TestVTK::CreateSphereTest(myPoint);
            filter->AddInputData(mySphere);
        }

        filter->Update();
        return filter->GetOutput();
    }

    RegistrationImageType::Pointer polyDataToImageData(vtkSmartPointer<vtkPolyData> polydata, double* spacing)
    {
        vtkSmartPointer<vtkImageData> imageData = vtkSmartPointer<vtkImageData>::New();
        double bounds[6];
        polydata->GetBounds(bounds);
        imageData->SetSpacing(spacing);

        //compute dimensions
        int dim[3];
        for (int i = 0; i < 3; ++i) {
            dim[i] = static_cast<int>(ceil((bounds[i * 2 + 1] - bounds[i * 2]) / spacing[i]));
        }
        imageData->SetDimensions(dim);
        imageData->SetExtent(0, dim[0] - 1, 0, dim[1] - 1, 0, dim[2] - 1);

        double origin[3];
        origin[0] = bounds[0] + spacing[0] / 2;
        origin[1] = bounds[2] + spacing[1] / 2;
        origin[2] = bounds[4] + spacing[2] / 2;
        imageData->SetOrigin(origin);

        imageData->AllocateScalars(VTK_UNSIGNED_CHAR, 1);

        //fill the image with foreground voxels:
        unsigned char inval = 255;
        unsigned char outval = 0;
        vtkIdType count = imageData->GetNumberOfPoints();
        for (vtkIdType i = 0; i < count; ++i) {
            imageData->GetPointData()->GetScalars()->SetTuple1(i, inval);
        }

        //polygonal data --> image stencil:
        vtkSmartPointer<vtkPolyDataToImageStencil> pdtoImageStencil = vtkSmartPointer<vtkPolyDataToImageStencil>::New();
        pdtoImageStencil->SetInputData(polydata);
        pdtoImageStencil->SetOutputOrigin(origin);
        pdtoImageStencil->SetOutputSpacing(spacing);
        pdtoImageStencil->SetOutputWholeExtent(imageData->GetExtent());
        pdtoImageStencil->Update();

        //cut the corresponding white image and set the background:
        vtkSmartPointer<vtkImageStencil> imageStencil = vtkSmartPointer<vtkImageStencil>::New();
        imageStencil->SetInputData(imageData);
        imageStencil->SetStencilConnection(pdtoImageStencil->GetOutputPort());
        imageStencil->ReverseStencilOff();
        imageStencil->SetBackgroundValue(outval);
        imageStencil->Update();

        imageData->DeepCopy(imageStencil->GetOutput());


        using FilterType = itk::VTKImageToImageFilter<SegmentImageType>;
        FilterType::Pointer filter = FilterType::New();
        filter->SetInput(imageData);

        try
        {
            filter->Update();
        }
        catch (itk::ExceptionObject & error)
        {
            std::cout << "**Error converting: " << error << std::endl;
        }

        auto caster = itk::CastImageFilter<SegmentImageType, RegistrationImageType>::New();
        caster->SetInput(filter->GetOutput());
        caster->Update();

        return caster->GetOutput();

        //return imageData;
    }

    vtkSmartPointer<vtkPolyData> TransformPoly(vtkSmartPointer<vtkPolyData> poly, itk::Matrix< double, 3, 3 > rotate, itk::Vector< double, 3 > translate)
    {
        vtkNew<vtkTransform> vtkTransform;

        vtkNew<vtkMatrix4x4> m;

        m->SetElement(0, 0, rotate[0][0]);
        m->SetElement(1, 0, rotate[1][0]);
        m->SetElement(2, 0, rotate[2][0]);
        m->SetElement(3, 0, 0);

        m->SetElement(0, 1, rotate[0][1]);
        m->SetElement(1, 1, rotate[1][1]);
        m->SetElement(2, 1, rotate[2][1]);
        m->SetElement(3, 1, 0);

        m->SetElement(0, 2, rotate[0][2]);
        m->SetElement(1, 2, rotate[1][2]);
        m->SetElement(2, 2, rotate[2][2]);
        m->SetElement(3, 2, 0);

        m->SetElement(0, 3, translate[0]);
        m->SetElement(1, 3, translate[1]);
        m->SetElement(2, 3, translate[2]);
        m->SetElement(3, 3, 1);

        vtkTransform->SetMatrix(m);

        vtkNew<vtkTransformPolyDataFilter> transformFilter;
        transformFilter->SetInputData(poly);
        transformFilter->SetTransform(vtkTransform);
        transformFilter->Update();

        auto resultTransform = transformFilter->GetOutput();

        //SavePolyData(resultTransform, "transform_poly.vtk");

        return resultTransform;
    }

    cv::Point3d getNearPoint(const vtkSmartPointer<vtkPolyData>& poly, const cv::Point3d& point)
    {
        vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
        implicitPolyDataDistance->SetInput(poly);
        double myClosest[3];
        double pnt[3] = { point.x, point.y, point.z };
        implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(pnt, myClosest);
        return cv::Point3d(myClosest[0], myClosest[1], myClosest[2]);
    }

    std::pair<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPolyData>> SplitPoly(vtkSmartPointer<vtkPolyData> poly, cv::Point3d a, cv::Point3d b, cv::Point3d c)
    {

        cv::Point3d directVector1 = b - a;
        cv::Point3d directVector2 = c - a;

        cv::Point3d mPoint = a;
        cv::Point3d normalVector = directVector1.cross(directVector2);
        normalVector = normalVector / sqrt(normalVector.dot(normalVector));

        vtkNew<vtkPlane> vtkPlaneMedial, vtkPlanelateral;

        vtkPlaneMedial->SetOrigin(mPoint.x, mPoint.y, mPoint.z);
        vtkPlaneMedial->SetNormal(normalVector.x, normalVector.y, normalVector.z);

        vtkNew<vtkPlaneCollection> medialPlanes, lateralPlanes;
        medialPlanes->AddItem(vtkPlaneMedial);

        vtkNew<vtkClipClosedSurface> medialClipper, lateralClipper;
        medialClipper->SetInputData(poly);
        medialClipper->SetClippingPlanes(medialPlanes);
        medialClipper->Update();

        normalVector = -normalVector;
        vtkPlanelateral->SetOrigin(mPoint.x, mPoint.y, mPoint.z);
        vtkPlanelateral->SetNormal(normalVector.x, normalVector.y, normalVector.z);

        lateralPlanes->AddItem(vtkPlanelateral);

        lateralClipper->SetInputData(poly);
        lateralClipper->SetClippingPlanes(lateralPlanes);
        lateralClipper->Update();

        std::pair<vtkSmartPointer<vtkPolyData>, vtkSmartPointer<vtkPolyData>> result;

        result.first = medialClipper->GetOutput();
        result.second = lateralClipper->GetOutput();

        return result;
    }

}