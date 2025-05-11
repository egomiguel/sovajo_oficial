#ifndef TEST_PKA_SUEN_NEW_H
#define TEST_PKA_SUEN_NEW_H

#include <iostream>
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include <QJsonDocument>
#include <QJsonObject>
#include <QJsonArray>
#include "itkPoint.h"
#include <QFile>
#include "uka_implants/FemurImplant.hpp"
#include "uka_implants/TibiaImplant.hpp"
#include "uka_implants/FemurImplantMatch.hpp"
#include "uka_implants/TibiaImplantMatch.hpp"
#include "uka_implants/TibiaSpacerImplant.hpp"
#include "uka_implants/ImplantsMatchFinalInfo.hpp"
#include "vtkPolyDataReader.h"
#include "vtkTransform.h"
#include "vtkMatrix4x4.h"
#include "vtkActor.h"
#include "vtkCamera.h"
#include "vtkPolyDataMapper.h"
#include <QFileInfo>
#include "vtkProperty.h"
#include "vtkSTLReader.h"
#include "vtkRenderWindowInteractor.h"
#include "uka_implants/FemurImplantThreePlane.hpp"
#include "itkVersorRigid3DTransform.h"

namespace TEST_PKA_SUEN_NEW
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
	vtkSmartPointer<vtkPolyData> readSTLData(QString &filePath)
	{
		vtkNew<vtkSTLReader> reader;
		reader->SetFileName(filePath.toStdString().c_str());
		reader->Update();
		return reader->GetOutput();
	}
	vtkSmartPointer<vtkPolyData> readVTKData(QString &filePath)
	{
		vtkNew<vtkPolyDataReader> reader;
		reader->SetFileName(filePath.toStdString().c_str());
		reader->Update();
		return reader->GetOutput();
	}
	itk::VersorRigid3DTransform<>::Pointer toItkTransform(const itk::Matrix<double, 3, 3> &rotate, const itk::Vector<double, 3> &translate)
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

	void testMatch()
	{
		QString dir = "D:\\sovajo\\Implants\\pka\\testImplants";
		auto landmarks = readLandmarks(QString("%1\\landmark.json").arg(dir));

		auto side = UKA::IMPLANTS::KneeSideEnum::KRight;
		auto surgerySide = UKA::IMPLANTS::SurgerySideEnum::KLateral;
		UKA::IMPLANTS::Knee knee;
		auto ankleCenter = landmarks.at(kMedialMalleolus) + (landmarks.at(kLateralMalleolus) - landmarks.at(kMedialMalleolus))*0.45;;

		auto femurData = readVTKData(QString("%1\\femur.vtk").arg(dir));
		auto tibiaData = readVTKData(QString("%1\\tibia.vtk").arg(dir));
		auto femurImplantData = readSTLData(QString("%1\\femur_LM_RL_SZ1.stl").arg(dir));
		auto tibiaImplant = createTibiaImplant(QString("%1\\tibia_LM_RL_A+_A#8mm_data.json").arg(dir));
		auto femurThreePlaneImplant = createFemurThreeImplant(femurImplantData, QString("%1\\femur_LM_RL_SZ1_data.json").arg(dir));
		
		knee.init(toPoint(landmarks[LandmarkType::kHipCenter]), toPoint(landmarks[LandmarkType::kFemurKneeCenter]),
			toPoint(landmarks[LandmarkType::kLateralEpicondyle]), toPoint(landmarks[LandmarkType::kMedialEpicondyle]),
			toPoint(landmarks[LandmarkType::kTibiaKneeCenter]), toPoint(landmarks[LandmarkType::kTibiaTuberosity]),
			toPoint(landmarks[LandmarkType::kPCLInsertionPoint]), toPoint(ankleCenter), 
			femurData, tibiaData, side, surgerySide, false);

		knee.setLateralAndMedialInferiorFemurPoints(toPoint(landmarks[LandmarkType::kFemurDistalLateral]), toPoint(landmarks[LandmarkType::kFemurDistalMedial]));
		knee.setLateralAndMedialPosteriorFemurPoints(toPoint(landmarks[LandmarkType::kFemurLateralPosteriorCondyle]), toPoint(landmarks[LandmarkType::kFemurMedialPosteriorCondyle]));
		knee.setLateralAndMedialPlateauPoints(toPoint(landmarks[LandmarkType::kTibiaLateralPlatformPoint]), toPoint(landmarks[LandmarkType::kTibiaMedialPlatformPoint]));

		UKA::IMPLANTS::TibiaImplantMatch tibiaMatch;
		tibiaMatch.init(*tibiaImplant, knee);
		auto implantToTibia = toItkTransform(tibiaMatch.GetRotationMatrix(), tibiaMatch.GetTranslationMatrix());

		UKA::IMPLANTS::FemurImplantMatch femurMatch;
		femurMatch.init(femurThreePlaneImplant.get(), knee);
		auto implantToFemur = toItkTransform(femurMatch.GetRotationMatrix(), femurMatch.GetTranslationMatrix());

		UKA::IMPLANTS::ImplantsMatchFinalInfo finalMatch(&knee, femurThreePlaneImplant.get(), *tibiaImplant, implantToFemur.GetPointer(), implantToTibia.GetPointer());

		std::cout << "FemurVarusAngle:" << finalMatch.GetFemurVarusAngle() << std::endl;
		std::cout << "FemurFlexionAngle:" << finalMatch.GetFemurFlexionAngle() << std::endl;
		//ERROR:this two angle is so big!
		std::cout << "FemurPCAAngle:" << finalMatch.GetFemurImplantPCAAngle() << std::endl;
		std::cout << "FemurTEAAngle:" << finalMatch.GetFemurImplantTEAAngle() << std::endl;
		//
		//finalMatch.setFemurVarusAngle(0);
		//finalMatch.setFemurFlexionAngle(0);
		//ERROR
		finalMatch.setFemurPCAAngle(0);
		auto implantToFemurNew = finalMatch.getITKFemurTransform();
		vtkNew<vtkActor> femurActor;
		femurActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		femurActor->GetMapper()->SetInputDataObject(femurData);
		femurActor->GetProperty()->SetColor(1, 1, 1);
		femurActor->GetProperty()->SetOpacity(0.7);

		vtkNew<vtkActor> femurImplantActor;
		femurImplantActor->SetMapper(vtkSmartPointer<vtkPolyDataMapper>::New());
		femurImplantActor->GetMapper()->SetInputDataObject(femurImplantData);
		femurImplantActor->GetProperty()->SetColor(0, 1, 0);
		femurImplantActor->SetUserTransform(toVtkTransform(implantToFemurNew->GetMatrix(), implantToFemurNew->GetTranslation()));

		vtkNew<vtkRenderer> render;
		render->AddActor(femurActor);
		render->AddActor(femurImplantActor);

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

	}
}


#endif