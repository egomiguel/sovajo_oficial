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
#include "vtkCamera.h"
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
#include "vtkAxesActor.h"
#include "vtkCaptionActor2D.h"
#include "vtkTextActor.h"
#include "vtkTextProperty.h"
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
	enum  LandmarkType
	{
		kHipCenter = 0,
		kMedialEpicondyle,
		kLateralEpicondyle,
		kFemurKneeCenter,

		// Tibia

		kTibiaKneeCenter,
		kTibiaTuberosity,
		kPCLInsertionPoint,
		kMedialMalleolus,
		kLateralMalleolus,

		//cutPoints
		kFemurDistalLateral,
		kFemurDistalMedial,
		kFemurMedialPosteriorCondyle,
		kFemurLateralPosteriorCondyle,
		kTibiaMedialPlatformPoint,
		kTibiaLateralPlatformPoint,

		kNoneLandmark

	};

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
		itk::Rigid3DTransform<double>::Pointer tibiaSideTransformOut = itk::VersorRigid3DTransform<double>::New();

		itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
		itk::Rigid3DTransform<double>::Pointer femurTransformOut = itk::VersorRigid3DTransform<double>::New();

		tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
		tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

		femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
		femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrix());

		std::vector<PointTypeITK> hullFemur = femurImplantMatch.GetHullPointsOnePlane(femurTransformIn, femurTransformOut, UKA::IMPLANTS::FemurImplantMatch::KOnePlanePosterior, 1, 1, 0);
		
		UKA::IMPLANTS::TibiaImplantMatch::HullPoints hullTemp = tibiaImplantMatch.GetHullPoints(tibiaTransformIn, tibiaTransformOut, tibiaSideTransformOut);
		std::vector<PointTypeITK> hull = hullFemur;
		
		std::cout << "Hull size: " << hull.size() << std::endl;

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
		matchFinalInfo.setFemurPCAAngle(0);

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


	UKA::IMPLANTS::Point toPoint(itk::Point<double> &p)
	{
		return UKA::IMPLANTS::Point(p[0], p[1], p[2]);
	}

	UKA::IMPLANTS::Point toPoint(const QJsonValue& val)
	{
		auto arr = val.toArray();

		itk::Point<double> p;
		for (int i = 0; i < 3; i++) {
			p[i] = arr[i].toDouble();
		}
		return UKA::IMPLANTS::Point(p[0], p[1], p[2]);;
	}

	UKA::IMPLANTS::Plane toPlane(const QJsonValue &val)
	{
		auto arr = val.toArray();
		std::vector<itk::Point<double>> points;
		for (auto i = 0; i < arr.size(); ++i)
		{
			points.push_back(toPoint(arr[i]).ToITKPoint());
		}
		double center[3];
		for (int i = 0; i < 3; ++i)
		{
			center[i] = (points[0][i] + points[1][i] + points[2][i] + points[3][i]) / 4;
		}
		auto p1 = points[1] - points[0];
		p1.Normalize();
		auto p2 = points[1] - points[2];
		p2.Normalize();
		itk::Vector<double> normal;
		vtkMath::Cross(p1.GetDataPointer(), p2.GetDataPointer(), normal.GetDataPointer());
		UKA::IMPLANTS::Plane plane;
		plane.init(UKA::IMPLANTS::Point(normal[0], normal[1], normal[2]), UKA::IMPLANTS::Point(center[0], center[1], center[2]));
		return plane;
	}

	std::shared_ptr<UKA::IMPLANTS::FemurImplant> createFemurThreeImplant(vtkSmartPointer<vtkPolyData> femurImplantData, const QString &femurImplantPath)
	{
		QFile file(femurImplantPath);
		if (!file.open(QIODevice::ReadOnly))
		{
			return nullptr;
		}
		QJsonParseError parseError;
		auto obj = QJsonDocument::fromJson(file.readAll()).object();


		UKA::IMPLANTS::FemurImplantInfo info;
		info.femurDistalThickness = obj["distal_thickness"].toDouble();
		info.femurPosteriorThickness = obj["posterior_thickness"].toDouble();
		auto femurImplant = std::make_shared<UKA::IMPLANTS::FemurImplantThreePlane>();
		femurImplant->init(toPlane(obj["posterior_plane"]), toPlane(obj["center_plane"]),
			toPlane(obj["anterior_plane"]), toPoint(obj["p1"]), toPoint(obj["p2"]), toPoint(obj["p3"]), femurImplantData, info);
		return femurImplant;
	}

	std::map<int, itk::Point<double>> readLandmarks(const QString &filePath)
	{
		std::map<int, itk::Point<double>> points;
		QFile file(filePath);
		if (!file.open(QIODevice::ReadOnly))
		{
			return points;
		}
		QJsonParseError err;
		auto obj = QJsonDocument::fromJson(file.readAll(), &err).object();
		if (err.error != QJsonParseError::NoError)
		{
			return points;
		}
		auto keys = obj.keys();

		for (int i = 0; i < keys.size(); ++i)
		{

			auto value = obj[keys[i]].toArray();
			itk::Point<double> point;
			point[0] = value[0].toDouble();
			point[1] = value[1].toDouble();
			point[2] = value[2].toDouble();

			points[keys[i].toInt()] = point;

		}
		return points;
	}

	vtkSmartPointer<vtkTransform> getImplantToFemur(const QString &jsonPath)
	{
		vtkNew<vtkTransform> transform;
		QFile file(jsonPath);
		if (file.open(QIODevice::ReadOnly))
		{
			auto obj = QJsonDocument::fromJson(file.readAll()).object();
			if (obj.contains("femur_matrix"))
			{
				auto array = obj["femur_matrix"].toArray();
				vtkNew<vtkMatrix4x4> matrix;
				for (int i = 0; i < 4; ++i)
				{
					for (int j = 0; j < 4; ++j)
					{
						matrix->SetElement(i, j, array[i * 4 + j].toDouble());
					}
				}
				vtkNew<vtkTransform> tranform;
				tranform->SetMatrix(matrix);
				return tranform;
			}
		}
		return vtkSmartPointer<vtkTransform>::New();
	}

	itk::Rigid3DTransform<>::Pointer toItkTransform(vtkSmartPointer<vtkTransform> trans)
	{
		itk::Vector<double> translation;
		itk::Matrix<double, 3, 3> rotation;
		auto matrix = trans->GetMatrix();
		for (int i = 0; i < 3; ++i)
		{
			for (int j = 0; j < 3; ++j)
			{
				rotation(i, j) = matrix->GetElement(i, j);
			}
			translation[i] = matrix->GetElement(i, 3);
		}
		auto transform = itk::VersorRigid3DTransform<>::New();
		transform->SetMatrix(rotation);
		transform->SetTranslation(translation);
		return transform;
	}

	itk::Rigid3DTransform<>::Pointer toItkTransform(const itk::Matrix<double, 3, 3> &rotate, const itk::Vector<double, 3> &translate)
	{
		auto trans = itk::VersorRigid3DTransform<>::New();
		trans->SetMatrix(rotate);
		trans->SetTranslation(translate);
		return trans;
	}

	vtkSmartPointer<vtkTransform> toVtkTransform(const itk::Matrix<double, 3, 3> &rotate, const itk::Vector<double, 3> &translate)
	{
		vtkNew<vtkMatrix4x4> matrix;
		matrix->Identity();

		for (int row = 0; row < 3; row++) {
			for (int col = 0; col < 3; col++) {
				matrix->SetElement(row, col, rotate(row, col));
			}
			matrix->SetElement(row, 3, translate[row]);
		}
		vtkNew<vtkTransform> transform;
		transform->SetMatrix(matrix);
		return transform;
	}

	vtkSmartPointer<vtkAxesActor> createAxesActor()
	{
		vtkNew<vtkAxesActor> axesActor;
		axesActor->GetXAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
		axesActor->GetYAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
		axesActor->GetZAxisCaptionActor2D()->GetTextActor()->SetTextScaleModeToNone();
		axesActor->GetXAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(16);
		axesActor->GetYAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(16);
		axesActor->GetZAxisCaptionActor2D()->GetCaptionTextProperty()->SetFontSize(16);
		return axesActor;
	}

	std::shared_ptr<UKA::IMPLANTS::TibiaImplant> createTibiaImplant(const QString &tibiaImplantFilePath)
	{
		QFile file(tibiaImplantFilePath);
		if (!file.open(QIODevice::ReadOnly))
		{
			return nullptr;
		}
		QJsonParseError parseError;
		auto obj = QJsonDocument::fromJson(file.readAll()).object();

		auto tibiaImplant = std::make_shared<UKA::IMPLANTS::TibiaImplant>();
		UKA::IMPLANTS::TibiaImplantInfo implantInfo;
		implantInfo.tibiaThickness = obj["thickness"].toDouble();
		implantInfo.tibiaSpacer = 0;
		tibiaImplant->init(toPoint(obj["pcl_point"]),
			toPoint(obj["tuber_point"]),
			toPoint(obj["ref_point"]),
			toPoint(obj["exterior_point"]),
			toPoint(obj["side_point"]),
			implantInfo);
		return tibiaImplant;
	}

	vtkSmartPointer<vtkTransform> getImplantToTibia(const QString &jsonPath)
	{
		vtkNew<vtkTransform> transform;
		QFile file(jsonPath);
		if (file.open(QIODevice::ReadOnly))
		{
			auto obj = QJsonDocument::fromJson(file.readAll()).object();
			if (obj.contains("tibia_matrix"))
			{
				auto array = obj["tibia_matrix"].toArray();
				vtkNew<vtkMatrix4x4> matrix;
				for (int i = 0; i < 4; ++i)
				{
					for (int j = 0; j < 4; ++j)
					{
						matrix->SetElement(i, j, array[i * 4 + j].toDouble());
					}
				}
				vtkNew<vtkTransform> tranform;
				tranform->SetMatrix(matrix);
				return tranform;
			}
		}
		return vtkSmartPointer<vtkTransform>::New();
	}

	void MatchPKAThreePlanes()
	{
		//std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data_2"; // left
		std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data"; // right
		UKA::IMPLANTS::Knee myKnee = CreateKneeFromFile_NumbersPKA(folder, UKA::IMPLANTS::KRight, UKA::IMPLANTS::KMedial);

		UKA::IMPLANTS::FemurImplant* femurImplant = new UKA::IMPLANTS::FemurImplantThreePlane();
		UKA::IMPLANTS::TibiaImplant tibiaImplant;

		//////////////////////////////////////////////Implants

		UKA::IMPLANTS::Plane pPosterior, pCenter, pAnterior;
		/*UKA::IMPLANTS::Point pRodTopPoint = UKA::IMPLANTS::Point(-29.55, -6.25, 21.69);
		UKA::IMPLANTS::Point pRodBaseExtremeSide1 = UKA::IMPLANTS::Point(-38.03, -8.55, 5.91);
		UKA::IMPLANTS::Point pRodBaseExtremeSide2 = UKA::IMPLANTS::Point(-20.95, -8.52, 5.96);

		pPosterior.init(UKA::IMPLANTS::Point(-29.2892, -18.0309, 27.4683), UKA::IMPLANTS::Point(-23.6379, -18.0308, 19.2064), UKA::IMPLANTS::Point(-34.6305, -18.0303, 20.6091));
		pCenter.init(UKA::IMPLANTS::Point(-24.5741, -15.8175, 12.3777), UKA::IMPLANTS::Point(-33.6575, -15.7641, 12.3235), UKA::IMPLANTS::Point(-29.2341, -13.4633, 10.0236));
		pAnterior.init(UKA::IMPLANTS::Point(-24.2632, 1.76137, 5.37675), UKA::IMPLANTS::Point(-30.5815, 2.39802, 5.3786), UKA::IMPLANTS::Point(-28.4761, -3.36245, 5.38769));
*/
		UKA::IMPLANTS::Point pRodTopPoint = UKA::IMPLANTS::Point(-34.7, 21.7, -6.7);
		UKA::IMPLANTS::Point pRodBaseExtremeSide1 = UKA::IMPLANTS::Point(-42.9, 5.9, -8.5);
		UKA::IMPLANTS::Point pRodBaseExtremeSide2 = UKA::IMPLANTS::Point(-26.0, 6.0, -8.6);

		pPosterior.init(UKA::IMPLANTS::Point(-39.3325, 24.2145, -18.0487), UKA::IMPLANTS::Point(-29.68, 23.2605, -18.0496), UKA::IMPLANTS::Point(-39.1445, 17.3601, -18.0491));
		pCenter.init(UKA::IMPLANTS::Point(-38.1176, 12.5853, -16.0614), UKA::IMPLANTS::Point(-31.2094, 12.7796, -16.2579), UKA::IMPLANTS::Point(-37.8839, 10.2226, -13.6989));
		pAnterior.init(UKA::IMPLANTS::Point(-34.1404, 5.39091, -3.37154), UKA::IMPLANTS::Point(-30.4084, 5.39085, 2.24808), UKA::IMPLANTS::Point(-36.0448, 5.39103, 2.7041));


		//auto femurModel = ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\3_planes\\femur_LL_RM.stl");
		auto femurModel = ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\testImplants\\femur_LM_RL_SZ1.stl");
		auto tibiaModel = ReadPolyDataSTL("D:\\sovajo\\Cases_Plan_PKA\\UKA-data\\implants\\tibia.STL");

		//QString dir = "D:\\sovajo\\Implants\\pka\\testImplants";
		//auto femurImplant = createFemurThreeImplant(femurModel, QString("%1\\femur_LM_RL_SZ1_data.json").arg(dir));

		UKA::IMPLANTS::FemurImplantInfo femurInfo;
		femurInfo.femurDistalThickness = 4.0;
		femurInfo.femurPosteriorThickness = 4.0;

		((UKA::IMPLANTS::FemurImplantThreePlane*)femurImplant)->init(pPosterior, pCenter, pAnterior, pRodTopPoint, pRodBaseExtremeSide1, pRodBaseExtremeSide2, femurModel, femurInfo);

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

		itk::Rigid3DTransform<double>::Pointer transformFemur = itk::VersorRigid3DTransform<double>::New();

		transformFemur->SetMatrix(femurImplantMatch.GetRotationMatrix());
		transformFemur->SetOffset(femurImplantMatch.GetTranslationMatrix());

		vtkSmartPointer<vtkPolyData> newImplantFemur = TransformPoly(femurModel, femurImplantMatch.GetRotationMatrix(), femurImplantMatch.GetTranslationMatrix());

		vtkSmartPointer<vtkPolyData> newImplantTibia = TransformPoly(tibiaModel, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());

		////////////////////// First Match
		std::vector<vtkSmartPointer<vtkPolyData>> polyList1, polyList2;
		//polyList.push_back(newImplantFemur);
		//polyList1.push_back(newImplantTibia);
		//show(myKnee.GetTibiaPoly(), polyList1);

		polyList2.push_back(newImplantTibia);
		show(myKnee.GetTibiaPoly(), polyList2);

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
		itk::Rigid3DTransform<double>::Pointer tibiaTransformOutSide = itk::VersorRigid3DTransform<double>::New();

		itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
		itk::Rigid3DTransform<double>::Pointer femurTransformOut = itk::VersorRigid3DTransform<double>::New();

		tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
		tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

		femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
		femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrix());
		
		std::vector<PointTypeITK> hullFemur = femurImplantMatch.GetHullPointsThreePlanes(femurTransformIn, femurTransformOut, UKA::IMPLANTS::FemurImplantMatch::KThreePlaneAnterior, 1, 1, 0);
		std::vector<PointTypeITK> hullTibia = (tibiaImplantMatch.GetHullPoints(tibiaTransformIn, tibiaTransformOut, tibiaTransformOutSide, 1, 1, 0, 5, 0)).implantPoints;

		std::vector<PointTypeITK> hull = hullTibia;
		std::cout << "Hull femur size: " << hullFemur.size() << std::endl;
		std::cout << "Hull tibia size: " << hullTibia.size() << std::endl;

		std::vector<cv::Point3d> tPoints;

		for (int i = 0; i < hull.size(); i++)
		{
			cv::Point3d myPoint(hull[i][0], hull[i][1], hull[i][2]);

			tPoints.push_back(myPoint);
		}

		show(myKnee.GetTibiaPoly(), tPoints, true);
		show(newImplantTibia, tPoints, true);

		std::vector<vtkSmartPointer<vtkPolyData>> polyList3;
		polyList3.push_back(newImplantTibia);
		show(myKnee.GetTibiaPoly(), polyList3);
		
		
		UKA::IMPLANTS::ImplantsMatchFinalInfo matchFinalInfo(&myKnee, femurImplant, tibiaImplant, femurTransformIn, tibiaTransformIn);
		matchFinalInfo.test();
		/*matchFinalInfo.SetTibiaProtrudes(-5);
		matchFinalInfo.SetFemurProtrudesAxial(-5);
		matchFinalInfo.SetFemurProtrudesCoronal(-5);
		matchFinalInfo.test();*/

		vtkSmartPointer<vtkPolyData> newImplantTibia3 = TransformPoly(tibiaModel, matchFinalInfo.getITKTibiaTransform()->GetMatrix(), matchFinalInfo.getITKTibiaTransform()->GetTranslation());
		//std::vector<vtkSmartPointer<vtkPolyData>> polyList3;
		//polyList.push_back(newImplantFemur);
		//polyList3.push_back(newImplantTibia3);
		//show(myKnee.GetTibiaPoly(), polyList3);
		//delete femurImplant;
		//femurImplant = NULL;
		
		
	}

	void testTibiaBounary()
	{
		QString dir = "D:\\sovajo\\Test_Cases\\octubre_2025\\utka_test_251015\\data";

		auto side = UKA::IMPLANTS::KneeSideEnum::KLeft;
		auto surgerySide = UKA::IMPLANTS::SurgerySideEnum::KMedial;

		auto femurData = ReadPolyData(QString("%1\\femur.vtk").arg(dir).toStdString());
		auto tibiaData = ReadPolyData(QString("%1\\tibia.vtk").arg(dir).toStdString());

		auto landmarks = readLandmarks(QString("%1\\landmark.json").arg(dir));
		auto ankleCenter = landmarks.at(kMedialMalleolus) + (landmarks.at(kLateralMalleolus) - landmarks.at(kMedialMalleolus))*0.45;;
		UKA::IMPLANTS::Knee knee;
		knee.init(toPoint(landmarks[LandmarkType::kHipCenter]), toPoint(landmarks[LandmarkType::kFemurKneeCenter]),
			toPoint(landmarks[LandmarkType::kLateralEpicondyle]), toPoint(landmarks[LandmarkType::kMedialEpicondyle]),
			toPoint(landmarks[LandmarkType::kTibiaKneeCenter]), toPoint(landmarks[LandmarkType::kTibiaTuberosity]),
			toPoint(landmarks[LandmarkType::kPCLInsertionPoint]), toPoint(ankleCenter), femurData, tibiaData, side, surgerySide, false);
		knee.setLateralAndMedialInferiorFemurPoints(toPoint(landmarks[LandmarkType::kFemurDistalLateral]), toPoint(landmarks[LandmarkType::kFemurDistalMedial]));
		knee.setLateralAndMedialPosteriorFemurPoints(toPoint(landmarks[LandmarkType::kFemurLateralPosteriorCondyle]), toPoint(landmarks[LandmarkType::kFemurMedialPosteriorCondyle]));
		knee.setLateralAndMedialPlateauPoints(toPoint(landmarks[LandmarkType::kTibiaLateralPlatformPoint]), toPoint(landmarks[LandmarkType::kTibiaMedialPlatformPoint]));

		auto tibiaImplantData = ReadPolyDataSTL(QString("%1/tibia_LM_RL_A+_A#8mm.stl").arg(dir).toStdString());
		auto tibiaImplant = createTibiaImplant(QString("%1/tibia_LM_RL_A+_A#8mm_data.json").arg(dir));
		auto implantToTibiaTrans = getImplantToTibia(QString("%1/plan.json").arg(dir));

		UKA::IMPLANTS::TibiaImplantMatch tibiaImplantMatch;
		tibiaImplantMatch.init(*tibiaImplant, knee);
		itk::Rigid3DTransform<>::Pointer boneToPlane = itk::VersorRigid3DTransform<>::New();
		itk::Rigid3DTransform<>::Pointer sideToPlane = itk::VersorRigid3DTransform<>::New();
		auto pointsInBone = tibiaImplantMatch.GetHullPoints(toItkTransform(implantToTibiaTrans), boneToPlane, sideToPlane,
			1, 2, 1,
			500).implantPoints;

		std::vector<cv::Point3d> tPoints;

		for (int i = 0; i < pointsInBone.size(); i++)
		{
			cv::Point3d myPoint(pointsInBone[i][0], pointsInBone[i][1], pointsInBone[i][2]);

			tPoints.push_back(myPoint);
		}

		show(knee.GetTibiaPoly(), tPoints, true);

		vtkSmartPointer<vtkPolyData> newImplantTibia = TransformPoly(tibiaImplantData, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());
		show(newImplantTibia, tPoints, true);

		/*
		vtkNew<vtkPoints> points;
		for (auto& p : pointsInBone)
		{
			points->InsertNextPoint(p.GetDataPointer());
		}
		vtkNew<vtkCellArray> lines;
		for (size_t i = 1; i < pointsInBone.size(); i++)
		{
			vtkNew<vtkLine> line;
			line->GetPointIds()->SetId(0, i);
			line->GetPointIds()->SetId(1, i - 1);
			lines->InsertNextCell(line);

		}
		vtkNew<vtkPolyData> borderdata;
		borderdata->SetPoints(points);
		borderdata->SetLines(lines);

		vtkNew<vtkActor> borderActor;
		borderActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		borderActor->GetMapper()->SetInputDataObject(borderdata);
		borderActor->GetProperty()->SetColor(1, 0, 0);
		borderActor->GetProperty()->SetLineWidth(2);

		vtkNew<vtkActor> implantActor;
		implantActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		implantActor->GetMapper()->SetInputDataObject(tibiaImplantData);
		implantActor->GetProperty()->SetColor(0, 1, 0);
		implantActor->SetUserTransform(implantToTibiaTrans);

		vtkNew<vtkActor> tibiaActor;
		tibiaActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		tibiaActor->GetMapper()->SetInputDataObject(tibiaData);
		tibiaActor->GetProperty()->SetColor(1, 1, 1);
		tibiaActor->GetProperty()->SetOpacity(0.7);

		vtkNew<vtkRenderer> render;
		render->AddActor(borderActor);
		render->AddActor(implantActor);
		render->AddActor(tibiaActor);

		vtkNew<vtkRenderWindow> renderWindow;
		renderWindow->AddRenderer(render);
		auto camera = render->GetActiveCamera();
		camera->SetFocalPoint(0, 0, 0);
		camera->SetPosition(0, -1, 0);
		camera->SetViewUp(0, 0, 1);
		render->ResetCamera();
		renderWindow->Render();
		vtkNew<vtkRenderWindowInteractor> interactor;
		renderWindow->SetInteractor(interactor);
		interactor->Start();
		*/
	}

	void testFemurAnteriodPlane()
	{
		QString dir = "D:\\sovajo\\Test_Cases\\agosto_2025\\UKA_Test";
		auto side = UKA::IMPLANTS::KneeSideEnum::KRight;
		auto surgerySide = UKA::IMPLANTS::SurgerySideEnum::KMedial;
		auto femurData = ReadPolyData(QString("%1\\femur.vtk").arg(dir).toStdString());
		auto tibiaData = ReadPolyData(QString("%1\\tibia.vtk").arg(dir).toStdString());
		auto femurImplantData = ReadPolyDataSTL(QString("%1\\femur_LL_RM_SZ4.stl").arg(dir).toStdString());
		auto femurThreePlaneImplant = createFemurThreeImplant(femurImplantData, QString("%1\\femur_LL_RM_SZ4_data.json").arg(dir));

		auto landmarks = readLandmarks(QString("%1\\landmark.json").arg(dir));
		auto ankleCenter = landmarks.at(kMedialMalleolus) + (landmarks.at(kLateralMalleolus) - landmarks.at(kMedialMalleolus))*0.45;;
		UKA::IMPLANTS::Knee knee;
		knee.init(toPoint(landmarks[LandmarkType::kHipCenter]), toPoint(landmarks[LandmarkType::kFemurKneeCenter]),
			toPoint(landmarks[LandmarkType::kLateralEpicondyle]), toPoint(landmarks[LandmarkType::kMedialEpicondyle]),
			toPoint(landmarks[LandmarkType::kTibiaKneeCenter]), toPoint(landmarks[LandmarkType::kTibiaTuberosity]),
			toPoint(landmarks[LandmarkType::kPCLInsertionPoint]), toPoint(ankleCenter), femurData, tibiaData, side, surgerySide, false);
		knee.setLateralAndMedialInferiorFemurPoints(toPoint(landmarks[LandmarkType::kFemurDistalLateral]), toPoint(landmarks[LandmarkType::kFemurDistalMedial]));
		knee.setLateralAndMedialPosteriorFemurPoints(toPoint(landmarks[LandmarkType::kFemurLateralPosteriorCondyle]), toPoint(landmarks[LandmarkType::kFemurMedialPosteriorCondyle]));
		knee.setLateralAndMedialPlateauPoints(toPoint(landmarks[LandmarkType::kTibiaLateralPlatformPoint]), toPoint(landmarks[LandmarkType::kTibiaMedialPlatformPoint]));
		UKA::IMPLANTS::FemurImplantMatch femurMatch;
		femurMatch.init(femurThreePlaneImplant.get(), knee);


		vtkNew<vtkRenderer> render1;
		vtkNew<vtkRenderer> render2;
		{
			auto implantToFemur = getImplantToFemur(QString("%1\\plan_1.json").arg(dir));;
			itk::Rigid3DTransform<>::Pointer boneToPlane = itk::VersorRigid3DTransform<>::New();
			femurMatch.GetHullPointsThreePlanes(toItkTransform(implantToFemur), boneToPlane, UKA::IMPLANTS::FemurImplantMatch::BoneAreaThreePlanes::KThreePlanePosterior);

			auto planeToBone = toVtkTransform(boneToPlane->GetMatrix(), boneToPlane->GetTranslation());
			planeToBone->Inverse();
			auto axesActor = createAxesActor();
			axesActor->SetTotalLength(30, 30, 30);
			axesActor->SetUserTransform(planeToBone);

			vtkNew<vtkActor> femurActor;
			femurActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
			femurActor->GetMapper()->SetInputDataObject(femurData);
			femurActor->GetProperty()->SetColor(1, 1, 1);
			femurActor->GetProperty()->SetOpacity(0.6);
			vtkNew<vtkActor> femurImplantActor;
			femurImplantActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
			femurImplantActor->GetMapper()->SetInputDataObject(femurImplantData);
			femurImplantActor->GetProperty()->SetColor(0, 1, 0);
			femurImplantActor->SetUserTransform(implantToFemur);
			femurImplantActor->GetProperty()->SetOpacity(0.8);

			render1->AddActor(femurActor);
			render1->AddActor(axesActor);
			render1->AddActor(femurImplantActor);
		}
		//2
		{
			auto implantToFemur = getImplantToFemur(QString("%1\\plan_2.json").arg(dir));

			itk::Rigid3DTransform<>::Pointer boneToPlane = itk::VersorRigid3DTransform<>::New();
			femurMatch.GetHullPointsThreePlanes(toItkTransform(implantToFemur), boneToPlane, UKA::IMPLANTS::FemurImplantMatch::BoneAreaThreePlanes::KThreePlanePosterior);
			auto planeToBone = toVtkTransform(boneToPlane->GetMatrix(), boneToPlane->GetTranslation());
			planeToBone->Inverse();

			auto axesActor = createAxesActor();
			axesActor->SetTotalLength(30, 30, 30);
			axesActor->SetUserTransform(planeToBone);


			vtkNew<vtkActor> femurActor;
			femurActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
			femurActor->GetMapper()->SetInputDataObject(femurData);
			femurActor->GetProperty()->SetColor(1, 1, 1);
			femurActor->GetProperty()->SetOpacity(0.6);
			vtkNew<vtkActor> femurImplantActor;
			femurImplantActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
			femurImplantActor->GetMapper()->SetInputDataObject(femurImplantData);
			femurImplantActor->GetProperty()->SetColor(0, 1, 0);
			femurImplantActor->SetUserTransform(implantToFemur);
			femurImplantActor->GetProperty()->SetOpacity(0.8);

			render2->AddActor(femurActor);
			render2->AddActor(axesActor);
			render2->AddActor(femurImplantActor);
		}

		vtkNew<vtkRenderWindow> renderWindow;
		render1->SetViewport(0, 0, 0.5, 1);
		render2->SetViewport(0.5, 0, 1.0, 1);
		renderWindow->AddRenderer(render1);
		renderWindow->AddRenderer(render2);

		auto camera = render1->GetActiveCamera();
		render2->SetActiveCamera(camera);
		camera->SetFocalPoint(0, 0, 0);
		camera->SetPosition(0, 0, -1);
		camera->SetViewUp(0, -1, 0);
		render1->ResetCamera();
		render2->ResetCamera();
		renderWindow->Render();
		vtkNew<vtkRenderWindowInteractor> interactor;
		renderWindow->SetInteractor(interactor);
		interactor->Start();
	}

	void testFemurAnteriodPlaneOctobre()
	{
		QString dir = "D:\\sovajo\\Test_Cases\\octubre_2025\\utka_test_251014\\data";


		auto side = UKA::IMPLANTS::KneeSideEnum::KRight;
		auto surgerySide = UKA::IMPLANTS::SurgerySideEnum::KMedial;


		auto femurData = ReadPolyData(QString("%1\\femur.vtk").arg(dir).toStdString());
		auto tibiaData = ReadPolyData(QString("%1\\tibia.vtk").arg(dir).toStdString());
		auto femurImplantData = ReadPolyDataSTL(QString("%1\\femur_LL_RM_SZ5.stl").arg(dir).toStdString());
		auto femurThreePlaneImplant = createFemurThreeImplant(femurImplantData, QString("%1\\femur_LL_RM_SZ5_data.json").arg(dir));

		auto landmarks = readLandmarks(QString("%1\\landmark.json").arg(dir));
		auto ankleCenter = landmarks.at(kMedialMalleolus) + (landmarks.at(kLateralMalleolus) - landmarks.at(kMedialMalleolus))*0.45;;
		UKA::IMPLANTS::Knee knee;
		knee.init(toPoint(landmarks[LandmarkType::kHipCenter]), toPoint(landmarks[LandmarkType::kFemurKneeCenter]),
			toPoint(landmarks[LandmarkType::kLateralEpicondyle]), toPoint(landmarks[LandmarkType::kMedialEpicondyle]),
			toPoint(landmarks[LandmarkType::kTibiaKneeCenter]), toPoint(landmarks[LandmarkType::kTibiaTuberosity]),
			toPoint(landmarks[LandmarkType::kPCLInsertionPoint]), toPoint(ankleCenter), femurData, tibiaData, side, surgerySide, false);
		knee.setLateralAndMedialInferiorFemurPoints(toPoint(landmarks[LandmarkType::kFemurDistalLateral]), toPoint(landmarks[LandmarkType::kFemurDistalMedial]));
		knee.setLateralAndMedialPosteriorFemurPoints(toPoint(landmarks[LandmarkType::kFemurLateralPosteriorCondyle]), toPoint(landmarks[LandmarkType::kFemurMedialPosteriorCondyle]));
		knee.setLateralAndMedialPlateauPoints(toPoint(landmarks[LandmarkType::kTibiaLateralPlatformPoint]), toPoint(landmarks[LandmarkType::kTibiaMedialPlatformPoint]));
		UKA::IMPLANTS::FemurImplantMatch femurMatch;
		femurMatch.init(femurThreePlaneImplant.get(), knee);
		auto implantToFemur = getImplantToFemur(QString("%1\\plan.json").arg(dir));;
		itk::Rigid3DTransform<>::Pointer boneToPlane = itk::VersorRigid3DTransform<>::New();
		auto pointsInBone = femurMatch.GetHullPointsThreePlanes(toItkTransform(implantToFemur), boneToPlane, UKA::IMPLANTS::FemurImplantMatch::BoneAreaThreePlanes::KThreePlaneAnterior);

		vtkNew<vtkPoints> points;
		for (auto& p : pointsInBone)
		{
			points->InsertNextPoint(p.GetDataPointer());
		}
		vtkNew<vtkCellArray> lines;
		for (size_t i = 1; i < pointsInBone.size(); i++)
		{
			vtkNew<vtkLine> line;
			line->GetPointIds()->SetId(0, i);
			line->GetPointIds()->SetId(1, i - 1);
			lines->InsertNextCell(line);

		}
		vtkNew<vtkPolyData> borderdata;
		borderdata->SetPoints(points);
		borderdata->SetLines(lines);

		vtkNew<vtkActor> borderActor;
		borderActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		borderActor->GetMapper()->SetInputDataObject(borderdata);
		borderActor->GetProperty()->SetColor(1, 0, 0);
		borderActor->GetProperty()->SetLineWidth(2);


		vtkNew<vtkActor> femurActor;
		femurActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		femurActor->GetMapper()->SetInputDataObject(femurData);
		femurActor->GetProperty()->SetColor(1, 1, 1);
		femurActor->GetProperty()->SetOpacity(0.6);

		vtkNew<vtkActor> femurImplantActor;
		femurImplantActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		femurImplantActor->GetMapper()->SetInputDataObject(femurImplantData);
		femurImplantActor->GetProperty()->SetColor(0, 1, 0);
		femurImplantActor->SetUserTransform(implantToFemur);
		femurImplantActor->GetProperty()->SetOpacity(0.8);


		vtkNew<vtkRenderWindow> renderWindow;
		vtkNew<vtkRenderer> render;
		render->AddActor(borderActor);
		render->AddActor(femurImplantActor);
		render->AddActor(femurActor);
		renderWindow->AddRenderer(render);

		auto camera = render->GetActiveCamera();
		camera->SetFocalPoint(0, 0, 0);
		camera->SetPosition(0, 0, -1);
		camera->SetViewUp(0, -1, 0);
		render->ResetCamera();
		renderWindow->Render();
		vtkNew<vtkRenderWindowInteractor> interactor;
		renderWindow->SetInteractor(interactor);
		interactor->Start();
	}


}


#endif