#ifndef TEST_PKA_H
#define TEST_PKA_H

#include <vtkNew.h>
#include <fstream>
#include "vtkCutter.h"
#include "vtkPlane.h"
#include "vtkFlyingEdges3D.h"
#include "vtkPolyDataConnectivityFilter.h"
#include "vtkNamedColors.h"
#include "vtkXMLPolyDataWriter.h"
#include "vtkPolyDataWriter.h"
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

#include <itkRigid3DTransform.h>
#include <itkVersorRigid3DTransform.h>
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"

#include <qjsonArray.h>
#include <qfile.h>
#include <qdir.h>
#include <qstring.h>
#include <qjsondocument.h>
#include <qjsonobject.h>
#include <qjsonvalue.h>

#include "uka_implants/Knee.hpp"
#include "uka_implants/FemurImplantMatch.hpp"
#include "uka_implants/TibiaImplantMatch.hpp"
#include "uka_implants/Point.hpp"
#include "uka_implants/TibiaSpacerImplant.hpp"
#include "uka_implants/TibiaSpacerImplantMatch.hpp"
#include "uka_implants/ImplantsMatchFinalInfo.hpp"
#include "uka_implants/FemurImplantThreePlane.hpp"
#include "uka_implants/FemurImplantOnePlane.hpp"
#include <vector>
#include <string>
#include <sstream>
#include <iostream>

namespace TEST_PKA
{
	vtkSmartPointer<vtkPolyData> ReadPolyData(const std::string& name)
	{
		vtkSmartPointer<vtkPolyDataReader> reader = vtkSmartPointer<vtkPolyDataReader>::New();
		reader->SetFileName(name.c_str());
		reader->Update();
		return reader->GetOutput();
	}

	vtkSmartPointer<vtkPolyData> ReadPolyDataSTL(const std::string& name)
	{
		vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
		reader->SetFileName(name.c_str());
		reader->Update();
		return reader->GetOutput();
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

		return resultTransform;
	}

	void show(vtkSmartPointer<vtkPolyData> poly, std::vector<vtkSmartPointer<vtkPolyData>> polyList = {})
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

	UKA::IMPLANTS::Point QJsonArrayToPoint(const QJsonArray& pArray)
	{
		if (pArray.size() == 3)
		{
			QJsonValue a = pArray[0];
			QJsonValue b = pArray[1];
			QJsonValue c = pArray[2];
			UKA::IMPLANTS::Point result = UKA::IMPLANTS::Point(a.toDouble(), b.toDouble(), c.toDouble());
			//std::cout << result << std::endl;
			return result;
		}
		else
		{
			return UKA::IMPLANTS::Point();
		}
	}


	UKA::IMPLANTS::Knee CreateKneeFromFile_NumbersPKA(const std::string& sourcePath, UKA::IMPLANTS::KneeSideEnum pSide = UKA::IMPLANTS::KneeSideEnum::KRight, UKA::IMPLANTS::SurgerySideEnum pSurgery = UKA::IMPLANTS::SurgerySideEnum::KLateral)
	{
		QDir directory(QString::fromStdString(sourcePath));
		QStringList files = directory.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
		if (files.size() == 0)
		{
			std::cout << "Error there are not json file" << std::endl;
			throw ("Error there are not json file");
		}

		QString jsonPath = directory.filePath(files[0]);
		QString var;
		QFile file;

		file.setFileName(jsonPath);
		file.open(QIODevice::ReadOnly | QIODevice::Text);
		var = file.readAll();
		file.close();
		QJsonDocument jsonDoc = QJsonDocument::fromJson(var.toUtf8());
		QJsonObject jsonObject = jsonDoc.object();

		UKA::IMPLANTS::Point hipCenter = QJsonArrayToPoint(jsonObject.value(QString("0")).toArray());
		UKA::IMPLANTS::Point anteriorCortex = QJsonArrayToPoint(jsonObject.value(QString("4")).toArray());
		UKA::IMPLANTS::Point femurKneeCenter = QJsonArrayToPoint(jsonObject.value(QString("3")).toArray());
		UKA::IMPLANTS::Point lateralEpicondyle = QJsonArrayToPoint(jsonObject.value(QString("2")).toArray());
		UKA::IMPLANTS::Point medialEpicondyle = QJsonArrayToPoint(jsonObject.value(QString("1")).toArray());

		UKA::IMPLANTS::Point tibiaKneeCenter = QJsonArrayToPoint(jsonObject.value(QString("5")).toArray());
		UKA::IMPLANTS::Point tibiaTubercle = QJsonArrayToPoint(jsonObject.value(QString("6")).toArray());
		UKA::IMPLANTS::Point tibiaPCL = QJsonArrayToPoint(jsonObject.value(QString("7")).toArray());
		UKA::IMPLANTS::Point lateralAnkle = QJsonArrayToPoint(jsonObject.value(QString("9")).toArray());
		UKA::IMPLANTS::Point medialAnkle = QJsonArrayToPoint(jsonObject.value(QString("8")).toArray());


		QString FemurPolyStr = jsonObject.value(QString("femur_poly")).toString();
		QString TibiaPolyStr = jsonObject.value(QString("tibia_poly")).toString();
		//QString PatellaPolyStr = jsonObject.value(QString("patella_poly")).toString();

		if (!(directory.exists(FemurPolyStr) && directory.exists(TibiaPolyStr)))
		{
			std::cout << "Error there are not vtk file" << std::endl;
			throw ("Error there are not vtk file");
		}

		QString femurPolyPath = directory.filePath(FemurPolyStr);
		QString tibiaPolyPath = directory.filePath(TibiaPolyStr);
		//QString patellaPolyPath = directory.filePath(PatellaPolyStr);

		vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly;

		femurPoly = ReadPolyData(femurPolyPath.toStdString());
		tibiaPoly = ReadPolyData(tibiaPolyPath.toStdString());

		UKA::IMPLANTS::Point ankleCenter = UKA::IMPLANTS::Knee::getComputeAnkleCenter(lateralAnkle, medialAnkle);

		//Patella myPatella;
		//myPatella.init(kneeCap, patellaLat, patellaMed, patellaInf, patellaPoly);

		UKA::IMPLANTS::Knee knee;

		knee.init(hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondyle, tibiaKneeCenter,
			tibiaTubercle, tibiaPCL, ankleCenter, femurPoly, tibiaPoly, pSide, pSurgery);

		return knee;
	}


	void MatchPKA()
	{
		std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data"; // right
		UKA::IMPLANTS::Knee myKnee = CreateKneeFromFile_NumbersPKA(folder, UKA::IMPLANTS::KRight, UKA::IMPLANTS::KMedial);
		std::cout << "111111111111111111111111111111" << std::endl;
		UKA::IMPLANTS::FemurImplant* femurImplant = new UKA::IMPLANTS::FemurImplantOnePlane();
		std::cout << "2222222222222222222222222222222222222" << std::endl;
		UKA::IMPLANTS::TibiaImplant tibiaImplant;

		//////////////////////////////////////////////Implants

		UKA::IMPLANTS::Plane pPosterior;
		pPosterior.init(UKA::IMPLANTS::Point(-1.98705693113, 16.0609375982, 8.76559347481), UKA::IMPLANTS::Point(-2.50396500712, 16.0613293561, -8.94437729719), UKA::IMPLANTS::Point(-7.92687551969, 16.0599612205, 8.9119337399));

		UKA::IMPLANTS::Point pRodBasePoint(17.1500611935, -3.64306703195, 0.0115669254999);
		UKA::IMPLANTS::Point pRodTopPoint(-1.57067067764, 0.0963311269987, -0.0533355484788);
		std::vector<UKA::IMPLANTS::Point> pSideBorder1 = { UKA::IMPLANTS::Point(8.48264081644,13.1976585493,-8.74732637542), UKA::IMPLANTS::Point(10.9804957096,10.5705588957,-8.91229039377),
														   UKA::IMPLANTS::Point(12.8703327684,7.76764520155,-8.80442944992), UKA::IMPLANTS::Point(14.1367283503,4.61611911302,-8.77980126905) };

		std::vector<UKA::IMPLANTS::Point> pSideBorder2 = { UKA::IMPLANTS::Point(8.46616178139,13.336942077,8.47133670984), UKA::IMPLANTS::Point(10.8359483977,10.9111407919,8.6214219954),
														   UKA::IMPLANTS::Point(12.9502176686,7.6594488483,8.67962746283), UKA::IMPLANTS::Point(14.174456106,4.73284273416,8.55774650189) };


		auto femurModel = ReadPolyDataSTL("D:\\sovajo\\Cases_Plan_PKA\\UKA-data\\implants\\femur.STL");
		auto tibiaModel = ReadPolyDataSTL("D:\\sovajo\\Cases_Plan_PKA\\UKA-data\\implants\\tibia.STL");

		UKA::IMPLANTS::FemurImplantInfo femurInfo;
		femurInfo.femurDistalThickness = 4.0;
		femurInfo.femurPosteriorThickness = 4.0;

		((UKA::IMPLANTS::FemurImplantOnePlane*)femurImplant)->init(pPosterior, pRodBasePoint, pRodTopPoint, pSideBorder1, pSideBorder2, femurModel, femurInfo);
		std::cout << "333333333333333333333333333333333333333" << std::endl;
		UKA::IMPLANTS::Point apLinePclPoint(10.35, 11.28, 1.15);
		UKA::IMPLANTS::Point apLineTuberPoint(10.01, -39.601, 1.06);
		UKA::IMPLANTS::Point sidePointUp(-15.2245, -10.73, 1.15);
		UKA::IMPLANTS::Point exteriorPointDown(-14.5152, -10.8294, 4.41175);
		UKA::IMPLANTS::Point planeSidePoint(12.10, -12.08, 2.59);
		UKA::IMPLANTS::TibiaImplantInfo tibiaInfo;
		tibiaInfo.tibiaThickness = 4.0;
		tibiaInfo.tibiaSpacer = 4.0;

		tibiaImplant.init(apLinePclPoint, apLineTuberPoint, sidePointUp, exteriorPointDown, planeSidePoint, tibiaInfo);

		////////////////////////////////////////////////////////

		UKA::IMPLANTS::FemurImplantMatch femurImplantMatch;
		UKA::IMPLANTS::TibiaImplantMatch tibiaImplantMatch;

		tibiaImplantMatch.init(tibiaImplant, myKnee);

		itk::Rigid3DTransform<double>::Pointer transformTibia = itk::VersorRigid3DTransform<double>::New();
		transformTibia->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
		transformTibia->SetOffset(tibiaImplantMatch.GetTranslationMatrix());
		femurImplantMatch.init(femurImplant, myKnee);
		std::cout << "4444444444444444444444444444444444444444" << std::endl;
		itk::Rigid3DTransform<double>::Pointer transformFemur = itk::VersorRigid3DTransform<double>::New();

		transformFemur->SetMatrix(femurImplantMatch.GetRotationMatrix());
		transformFemur->SetOffset(femurImplantMatch.GetTranslationMatrix());


		vtkSmartPointer<vtkPolyData> newImplantFemur = TransformPoly(femurModel, femurImplantMatch.GetRotationMatrix(), femurImplantMatch.GetTranslationMatrix());

		vtkSmartPointer<vtkPolyData> newImplantTibia = TransformPoly(tibiaModel, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());


		////////////////////// First Match
		std::vector<vtkSmartPointer<vtkPolyData>> polyList1, polyList2;
		//polyList.push_back(newImplantFemur);
		polyList1.push_back(newImplantTibia);
		show(myKnee.GetTibiaPoly(), polyList1);

		polyList2.push_back(newImplantFemur);
		show(myKnee.GetFemurPoly(), polyList2);

		///////////////////////////////////////////////////// Second Match
		/*
		UKA::IMPLANTS::ImplantsMatchFinalInfo info(&myKnee, femurImplant, tibiaImplant, transformFemur, transformTibia);
		//auto transformFemurFinal = info.FemurImplantToTibiaImplant();
		auto transformTibiaFinal = info.TibiaImplantToFemurImplant();

		//vtkSmartPointer<vtkPolyData> newImplantFemurFinal = TestVTK::TransformPoly(femurModel, transformFemurFinal->GetMatrix(), transformFemurFinal->GetOffset());
		vtkSmartPointer<vtkPolyData> newImplantTibiaFinal = TransformPoly(tibiaModel, transformTibiaFinal->GetMatrix(), transformTibiaFinal->GetOffset());

		//std::vector<vtkSmartPointer<vtkPolyData>> polyListFinal;
		//polyListFinal.push_back(newImplantFemur);
		//polyListFinal.push_back(newImplantTibiaFinal);
		//show(myKnee.GetTibiaPoly(), polyListFinal);
		*/
		itk::Rigid3DTransform<double>::Pointer tibiaTransformIn = itk::VersorRigid3DTransform<double>::New();
		itk::Rigid3DTransform<double>::Pointer tibiaTransformOut = itk::VersorRigid3DTransform<double>::New();

		itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
		itk::Rigid3DTransform<double>::Pointer femurTransformOut = itk::VersorRigid3DTransform<double>::New();

		tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
		tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

		femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
		femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrix());

		std::vector<PointTypeITK> hullFemur = femurImplantMatch.GetHullPointsOnePlane(femurTransformIn, femurTransformOut, UKA::IMPLANTS::FemurImplantMatch::KOnePlaneAnteriorAndDistalCurve);
		
		std::vector<PointTypeITK> hull = tibiaImplantMatch.GetHullPoints(tibiaTransformIn, tibiaTransformOut);

		
		std::cout << "Hull size: " << hull.size() << std::endl;

		std::vector<cv::Point3d> tPoints;

		for (int i = 0; i < hull.size(); i++)
		{
			cv::Point3d myPoint(hull[i][0], hull[i][1], hull[i][2]);

			tPoints.push_back(myPoint);
		}

		//show(myKnee.GetTibiaPoly(), tPoints, true);
		//show(newImplantTibia, tPoints, true);

		UKA::IMPLANTS::ImplantsMatchFinalInfo matchFinalInfo(&myKnee, femurImplant, tibiaImplant, femurTransformIn, tibiaTransformIn);
		std::cout << "5555555555555555555555555555555555555555" << std::endl;
		matchFinalInfo.test();
		std::cout << "666666666666666666666666666666666" << std::endl;
		matchFinalInfo.SetTibiaProtrudes(-5);
		matchFinalInfo.SetFemurProtrudesAxial(-5);
		matchFinalInfo.SetFemurProtrudesCoronal(-5);

		std::cout << "777777777777777777777777777777777777777777" << std::endl;

		matchFinalInfo.test();

		std::cout << "888888888888888888888888888888888888888888" << std::endl;

		vtkSmartPointer<vtkPolyData> newImplantTibia3 = TransformPoly(tibiaModel, matchFinalInfo.getITKTibiaTransform()->GetMatrix(), matchFinalInfo.getITKTibiaTransform()->GetTranslation());
		std::vector<vtkSmartPointer<vtkPolyData>> polyList3;
		//polyList.push_back(newImplantFemur);
		polyList3.push_back(newImplantTibia3);
		show(myKnee.GetTibiaPoly(), polyList3);
		std::cout << "9999999999999999999999999999999999999999999999999999" << std::endl;
		delete femurImplant;
		femurImplant = NULL;
	}

	void MatchPKAThreePlanes()
	{
		std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data_2"; // left
		UKA::IMPLANTS::Knee myKnee = CreateKneeFromFile_NumbersPKA(folder, UKA::IMPLANTS::KLeft, UKA::IMPLANTS::KMedial);

		UKA::IMPLANTS::FemurImplant* femurImplant = new UKA::IMPLANTS::FemurImplantThreePlane();
		UKA::IMPLANTS::TibiaImplant tibiaImplant;

		//////////////////////////////////////////////Implants

		UKA::IMPLANTS::Plane pPosterior, pCenter, pAnterior;
		UKA::IMPLANTS::Point pRodTopPoint = UKA::IMPLANTS::Point(-29.55, -6.25, 21.69);
		UKA::IMPLANTS::Point pRodBaseExtremeSide1 = UKA::IMPLANTS::Point(-38.03, -8.55, 5.91);
		UKA::IMPLANTS::Point pRodBaseExtremeSide2 = UKA::IMPLANTS::Point(-20.95, -8.52, 5.96);

		pPosterior.init(UKA::IMPLANTS::Point(-29.2892, -18.0309, 27.4683), UKA::IMPLANTS::Point(-23.6379, -18.0308, 19.2064), UKA::IMPLANTS::Point(-34.6305, -18.0303, 20.6091));
		pCenter.init(UKA::IMPLANTS::Point(-24.5741, -15.8175, 12.3777), UKA::IMPLANTS::Point(-33.6575, -15.7641, 12.3235), UKA::IMPLANTS::Point(-29.2341, -13.4633, 10.0236));
		pAnterior.init(UKA::IMPLANTS::Point(-24.2632, 1.76137, 5.37675), UKA::IMPLANTS::Point(-30.5815, 2.39802, 5.3786), UKA::IMPLANTS::Point(-28.4761, -3.36245, 5.38769));

		auto femurModel = ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\3_planes\\femur_LL_RM.stl");
		auto tibiaModel = ReadPolyDataSTL("D:\\sovajo\\Cases_Plan_PKA\\UKA-data\\implants\\tibia.STL");

		UKA::IMPLANTS::FemurImplantInfo femurInfo;
		femurInfo.femurDistalThickness = 4.0;
		femurInfo.femurPosteriorThickness = 4.0;
		std::cout << "aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa" << std::endl;
		((UKA::IMPLANTS::FemurImplantThreePlane*)femurImplant)->init(pPosterior, pCenter, pAnterior, pRodTopPoint, pRodBaseExtremeSide1, pRodBaseExtremeSide2, femurModel, femurInfo);
		std::cout << "bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb" << std::endl;
		UKA::IMPLANTS::Point apLinePclPoint(10.35, 11.28, 1.15);
		UKA::IMPLANTS::Point apLineTuberPoint(10.01, -39.601, 1.06);
		UKA::IMPLANTS::Point sidePointUp(-15.2245, -10.73, 1.15);
		UKA::IMPLANTS::Point exteriorPointDown(-14.5152, -10.8294, 4.41175);
		UKA::IMPLANTS::Point planeSidePoint(12.10, -12.08, 2.59);
		UKA::IMPLANTS::TibiaImplantInfo tibiaInfo;
		tibiaInfo.tibiaThickness = 4.0;
		tibiaInfo.tibiaSpacer = 4.0;

		tibiaImplant.init(apLinePclPoint, apLineTuberPoint, sidePointUp, exteriorPointDown, planeSidePoint, tibiaInfo);
		std::cout << "cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc" << std::endl;
		////////////////////////////////////////////////////////

		UKA::IMPLANTS::FemurImplantMatch femurImplantMatch;
		UKA::IMPLANTS::TibiaImplantMatch tibiaImplantMatch;

		tibiaImplantMatch.init(tibiaImplant, myKnee);
		std::cout << "ddddddddddddddddddddddddddddddddddddddddddddddddddddddddddd" << std::endl;
		itk::Rigid3DTransform<double>::Pointer transformTibia = itk::VersorRigid3DTransform<double>::New();
		transformTibia->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
		transformTibia->SetOffset(tibiaImplantMatch.GetTranslationMatrix());
		femurImplantMatch.init(femurImplant, myKnee);
		std::cout << "eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee" << std::endl;
		itk::Rigid3DTransform<double>::Pointer transformFemur = itk::VersorRigid3DTransform<double>::New();

		transformFemur->SetMatrix(femurImplantMatch.GetRotationMatrix());
		transformFemur->SetOffset(femurImplantMatch.GetTranslationMatrix());
		std::cout << "ffffffffffffffffffffffffffffffffffffffffffff" << std::endl;
		vtkSmartPointer<vtkPolyData> newImplantFemur = TransformPoly(femurModel, femurImplantMatch.GetRotationMatrix(), femurImplantMatch.GetTranslationMatrix());

		vtkSmartPointer<vtkPolyData> newImplantTibia = TransformPoly(tibiaModel, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());
		std::cout << "ggggggggggggggggggggggggggggggggggggggggggggggggggg" << std::endl;
		////////////////////// First Match
		std::vector<vtkSmartPointer<vtkPolyData>> polyList1, polyList2;
		//polyList.push_back(newImplantFemur);
		polyList1.push_back(newImplantTibia);
		show(myKnee.GetTibiaPoly(), polyList1);

		polyList2.push_back(newImplantFemur);
		show(myKnee.GetFemurPoly(), polyList2);
		std::cout << "hhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh" << std::endl;
		///////////////////////////////////////////////////// Second Match
		/*
		UKA::IMPLANTS::ImplantsMatchFinalInfo info(&myKnee, femurImplant, tibiaImplant, transformFemur, transformTibia);
		//auto transformFemurFinal = info.FemurImplantToTibiaImplant();
		auto transformTibiaFinal = info.TibiaImplantToFemurImplant();

		//vtkSmartPointer<vtkPolyData> newImplantFemurFinal = TestVTK::TransformPoly(femurModel, transformFemurFinal->GetMatrix(), transformFemurFinal->GetOffset());
		vtkSmartPointer<vtkPolyData> newImplantTibiaFinal = TransformPoly(tibiaModel, transformTibiaFinal->GetMatrix(), transformTibiaFinal->GetOffset());

		//std::vector<vtkSmartPointer<vtkPolyData>> polyListFinal;
		//polyListFinal.push_back(newImplantFemur);
		//polyListFinal.push_back(newImplantTibiaFinal);
		//show(myKnee.GetTibiaPoly(), polyListFinal);
		*/


		itk::Rigid3DTransform<double>::Pointer tibiaTransformIn = itk::VersorRigid3DTransform<double>::New();
		itk::Rigid3DTransform<double>::Pointer tibiaTransformOut = itk::VersorRigid3DTransform<double>::New();

		itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
		itk::Rigid3DTransform<double>::Pointer femurTransformOut = itk::VersorRigid3DTransform<double>::New();

		tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
		tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

		femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
		femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrix());

		std::cout << "111111111111111111111111111111111111111111" << std::endl;
		
		std::vector<PointTypeITK> hullFemur = femurImplantMatch.GetHullPointsThreePlanes(femurTransformIn, femurTransformOut, UKA::IMPLANTS::FemurImplantMatch::KThreePlaneAnterior, 1, 1);
		std::cout << "22222222222222222222222222222222222222222" << std::endl;
		std::vector<PointTypeITK> hullTest = femurImplantMatch.GetHullPointsOnePlane(femurTransformIn, femurTransformOut, UKA::IMPLANTS::FemurImplantMatch::KOnePlaneAnteriorAndDistalCurve);

		std::cout << "33333333333333333333333333333333333333333" << std::endl;
		std::vector<PointTypeITK> hull = hullFemur;
		std::cout << "Hull size: " << hull.size() << std::endl;
		std::cout << "Hull test size: " << hullTest.size() << std::endl;

		std::vector<cv::Point3d> tPoints;

		for (int i = 0; i < hull.size(); i++)
		{
			cv::Point3d myPoint(hull[i][0], hull[i][1], hull[i][2]);

			tPoints.push_back(myPoint);
		}

		show(myKnee.GetFemurPoly(), tPoints, true);
		show(newImplantFemur, tPoints, true);
		
		
		UKA::IMPLANTS::ImplantsMatchFinalInfo matchFinalInfo(&myKnee, femurImplant, tibiaImplant, femurTransformIn, tibiaTransformIn);
		std::cout << "5555555555555555555555555555555555555555" << std::endl;
		matchFinalInfo.test();
		std::cout << "666666666666666666666666666666666" << std::endl;
		matchFinalInfo.SetTibiaProtrudes(-5);
		matchFinalInfo.SetFemurProtrudesAxial(-5);
		matchFinalInfo.SetFemurProtrudesCoronal(-5);

		std::cout << "777777777777777777777777777777777777777777" << std::endl;

		matchFinalInfo.test();

		std::cout << "888888888888888888888888888888888888888888" << std::endl;

		vtkSmartPointer<vtkPolyData> newImplantTibia3 = TransformPoly(tibiaModel, matchFinalInfo.getITKTibiaTransform()->GetMatrix(), matchFinalInfo.getITKTibiaTransform()->GetTranslation());
		std::vector<vtkSmartPointer<vtkPolyData>> polyList3;
		//polyList.push_back(newImplantFemur);
		polyList3.push_back(newImplantTibia3);
		show(myKnee.GetTibiaPoly(), polyList3);
		std::cout << "9999999999999999999999999999999999999999999999999999" << std::endl;
		delete femurImplant;
		femurImplant = NULL;
		
		
	}
}


#endif