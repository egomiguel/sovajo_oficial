#include <iostream>
#include <QJsonDocument>
#include <QFile>
#include "vtkPolyData.h"
#include "vtkPolyDataReader.h"
#include <QJsonObject>
#include <QJsonArray>
#include "vtkstlreader.h"
#include "itkPoint.h"
#include <QDate>
#include <QDateTime>
#include "vtkFlyingEdges3D.h"
#include "vtkPolyDataWriter.h"
#include "itkImageFileReader.h"
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
//#include "itkImageToVtkImageFilter.h"
#include "vtkMarchingCubes.h"
#include <QFileInfo>
#include <QJsonArray>
#include <QDir>
#include "itkVersorRigid3DTransform.h"
#include "implants/Knee.hpp"
#include "implants/FemurImplant.hpp"
#include "implants/FemurImplantMatch.hpp"
#include "implants/TibiaImplant.hpp"
#include "implants/TibiaImplantMatch.hpp"
#include "implants/ImplantsMatchFinalInfo.hpp"
#include "implants/Point.hpp"
#include "implants/Utils.hpp"
#include "implants/ImplantTools.hpp"

using namespace TKA::IMPLANTS;

namespace TEST_IMPLANTS_TIBIA
{
	enum LandmarkType {
		// Femur
		kHipCenter = 0,
		kMedialEpicondyle,
		kLateralEpicondyle,
		kFemurKneeCenter,
		kFemurAnteriorCortex,//kAnteriorFemurPoint,
		kFemurMedialPosteriorCondyle,
		kFemurLateralPosteriorCondyle,

		// Tibia
		kTibiaMedialPlatformPoint,
		kTibiaLateralPlatformPoint,
		kTibiaKneeCenter,
		kTibiaTuberosity,
		kPCLInsertionPoint,
		kMedialMalleolus,
		kLateralMalleolus,

		// KneeCap

		kKneeCapLateralSide,
		kKneeCapMedialSide,
		kKneeCapInferiorSide,
		kKneeCapRotatePoint,

		//auto 
		kFemurDistalLateral,
		kFemurDistalMedial,
		kNoneLandmark

	};
	vtkSmartPointer<vtkPolyData> readSTL(const char *stlPath)
	{
		vtkNew<vtkSTLReader> reader;
		reader->SetFileName(stlPath);
		reader->Update();
		return reader->GetOutput();
	}
	vtkSmartPointer<vtkPolyData> readVTK(const char *stlPath)
	{
		vtkNew<vtkPolyDataReader> reader;
		reader->SetFileName(stlPath);
		reader->Update();
		return reader->GetOutput();
	}
	TKA::IMPLANTS::Point readPoint(const QJsonValue& val)
	{
		auto arr = val.toArray();
		return TKA::IMPLANTS::Point(arr[0].toDouble(), arr[1].toDouble(), arr[2].toDouble());
	}
	TKA::IMPLANTS::Plane readPlane(const QJsonValue& val) {
		std::array<double, 6> res;

		auto arr = val.toArray();

		for (int i = 0; i < 6; i++) {
			res[i] = arr[i].toDouble();
		}
		TKA::IMPLANTS::Plane plane;
		plane.init(TKA::IMPLANTS::Point(res[0], res[1], res[2]),
			TKA::IMPLANTS::Point(res[3], res[4], res[5]));
		return plane;

	}

	std::shared_ptr<TKA::IMPLANTS::TibiaImplant> createTibiaImplant(const char *jsonPath, const char *stlPath)
	{
		QFile file(jsonPath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();

		TKA::IMPLANTS::TibiaImplantInfo tibiaInfo;
		tibiaInfo.tibiaLateralThickness = obj["thickness_left"].toDouble();;
		tibiaInfo.tibiaMedialThickness = obj["thickness_right"].toDouble();;

		auto tibiaImplant = std::make_shared<TKA::IMPLANTS::TibiaImplant>();
		tibiaImplant->init(readPoint(obj["point_pcl1"]), readPoint(obj["point_pcl2"]),
			readPoint(obj["point_front"]), readPoint(obj["point_external"]), tibiaInfo);
		return tibiaImplant;
	}

	std::shared_ptr<TKA::IMPLANTS::FemurImplant> createFemurImplant(const char *jsonPath, const char *stlPath)
	{
		QFile file(jsonPath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();

		auto femurImplant = std::make_shared<TKA::IMPLANTS::FemurImplant>();

		TKA::IMPLANTS::FemurImplantInfo  femurInfo;
		femurInfo.femurPosteriorLateralThickness = obj["posterior_lateral_thickness"].toDouble();
		femurInfo.femurDistalLateralThickness = obj["distal_lateral_thickness"].toDouble();
		femurInfo.femurPosteriorMedialThickness = obj["posterior_medial_thickness"].toDouble();
		femurInfo.femurDistalMedialThickness = obj["distal_medial_thickness"].toDouble();

		std::vector<itk::Point<double>> patellaPoints;

		femurImplant->init(
			readPlane(obj["face_a"]),
			readPlane(obj["face_b"]),
			readPlane(obj["face_c"]),
			readPlane(obj["face_d"]),
			readPlane(obj["face_e"]),
			readPoint(obj["p1"]), readPoint(obj["p2"]), readPoint(obj["p3"]),
			readSTL(stlPath),
			patellaPoints, femurInfo
		);

		return femurImplant;
	}
	std::map<LandmarkType, TKA::IMPLANTS::Point> readLandmarks(const char *filePath)
	{
		std::map<LandmarkType, TKA::IMPLANTS::Point> items;
		QFile file(filePath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read landmarks failed" << std::endl;
			return items;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();
		for (auto key : obj.keys())
		{

			items[(LandmarkType)key.toInt()] = readPoint(obj[key]);
		}
		return items;
	}
	std::shared_ptr<TKA::IMPLANTS::Knee> createKnee(std::map<LandmarkType, TKA::IMPLANTS::Point> &landmarks,
		const char *femurVtkPath, const char *tibiaVtkPath, bool left)
	{
		auto knee = std::make_shared<TKA::IMPLANTS::Knee>();
		auto ankleCenter = landmarks.at(kMedialMalleolus).ToITKPoint() +
			(landmarks.at(kLateralMalleolus).ToITKPoint() - landmarks.at(kMedialMalleolus).ToITKPoint())*0.45;
		TKA::IMPLANTS::Patella patella;//
		knee->init(landmarks[kHipCenter],
			landmarks[kFemurAnteriorCortex],
			landmarks[kFemurKneeCenter],
			landmarks[kLateralEpicondyle],
			landmarks[kMedialEpicondyle],
			landmarks[kTibiaKneeCenter],
			landmarks[kTibiaTuberosity],
			landmarks[kPCLInsertionPoint], TKA::IMPLANTS::Point(ankleCenter[0], ankleCenter[1], ankleCenter[2]), patella,
			readVTK(femurVtkPath), readVTK(tibiaVtkPath), left ? TKA::IMPLANTS::KneeSideEnum::KLeft : TKA::IMPLANTS::KneeSideEnum::KRight,
			true);
		//knee->setLateralAndMedialInferiorFemurPoints(landmarks[kFemurDistalLateral], landmarks[kFemurDistalMedial]);
		//knee->setLateralAndMedialPosteriorFemurPoints(landmarks[kFemurLateralPosteriorCondyle], landmarks[kFemurMedialPosteriorCondyle]);
		//knee->setLateralAndMedialPlateauPoints(landmarks[kTibiaLateralPlatformPoint], landmarks[kTibiaMedialPlatformPoint]);
		return knee;
	}


	itk::VersorRigid3DTransform<double>::Pointer toTransformInverse(double ele[16])
	{
		itk::Rigid3DTransform<>::MatrixType m;
		itk::Rigid3DTransform<>::OffsetType offset;
		for (int row = 0; row < 3; row++)
		{
			for (int col = 0; col < 3; col++)
			{
				m(row, col) = ele[row * 4 + col];
			}
			offset[row] = ele[row * 4 + 3];
		}
		auto trans = itk::VersorRigid3DTransform<>::New();
		trans->SetMatrix(m);
		trans->SetOffset(offset);
		auto trans2 = itk::VersorRigid3DTransform<>::New();
		trans->GetInverse(trans2);
		return trans2;
	}

	void printInfo(const TKA::IMPLANTS::ImplantsMatchFinalInfo &info)
	{
		std::cout << "FemurVarusAngle:" << info.GetFemurVarusAngle() << std::endl;
		std::cout << "FemurImplantTEAAngle:" << info.GetFemurImplantTEAAngle() << std::endl;
		std::cout << "FemurImplantPCAAngle:" << info.GetFemurImplantPCAAngle() << std::endl;
		std::cout << "FemurFlexionAngle:" << info.GetFemurFlexionAngle() << std::endl;
		auto axial = info.GetFemurResectionAxial();
		std::cout << "FemurResectionAxial:" << axial.lateral << "," << axial.medial << std::endl;
		auto coronal = info.GetFemurResectionCoronal();
		std::cout << "FemurResectionCoronal:" << coronal.lateral << "," << coronal.medial << std::endl;

		std::cout << "****************************************" << std::endl;
	}
	void printTibia(const TKA::IMPLANTS::ImplantsMatchFinalInfo &info)
	{
		std::cout << "TibiaImplantRotationAngle:" << info.GetTibiaImplantRotationAngle() << std::endl;
		std::cout << "TibiaVarusAngle:" << info.GetTibiaVarusAngle() << std::endl;
		std::cout << "TibiaImplantSlopeAngle:" << info.GetTibiaImplantSlopeAngle() << std::endl;
		auto resection = info.GetTibiaResection();
		std::cout << "TibiaResection:" << resection.lateral << "," << resection.medial << std::endl;
		std::cout << "****************************************" << std::endl;
	}
	bool equal(double a, double b, double tolerance = 1e-6)
	{
		return std::abs(a - b) < tolerance;
	}

	bool equal(itk::VersorRigid3DTransform<>::Pointer p1, itk::VersorRigid3DTransform<>::Pointer p2)
	{
		auto t1 = p1->GetTranslation();
		auto t2 = p1->GetTranslation();
		for (auto i = 0; i < 3; ++i)
		{
			if (!equal(t1[i], t2[i]))
			{
				std::cout << "i:" << t1[i] << "," << t2[i] << std::endl;
				return false;
			}
		}
		auto v1 = p1->GetVersor();
		auto v2 = p2->GetVersor();
		return equal(v1.GetW(), v2.GetW()) && equal(v1.GetX(), v2.GetX())
			&& equal(v1.GetY(), v2.GetY()) && equal(v1.GetZ(), v2.GetZ());
	}

	itk::VersorRigid3DTransform<>::Pointer matchTibia(std::shared_ptr<TKA::IMPLANTS::Knee> knee, std::shared_ptr<TKA::IMPLANTS::TibiaImplant> tibiaImplant)
	{
		auto tibiaMatch = std::make_shared<TKA::IMPLANTS::TibiaImplantMatch>();
		tibiaMatch->init(*tibiaImplant, *knee);

		auto trans = itk::VersorRigid3DTransform<>::New();
		trans->SetMatrix(tibiaMatch->GetRotationMatrix());
		trans->SetTranslation(tibiaMatch->GetTranslationMatrix());
		return trans;
	}
	itk::VersorRigid3DTransform<>::Pointer matchFemur(std::shared_ptr<TKA::IMPLANTS::Knee> knee, std::shared_ptr<TKA::IMPLANTS::FemurImplant> femurImplant)
	{
		auto femurMatch = std::make_shared<TKA::IMPLANTS::FemurImplantMatch>();
		femurMatch->init(*femurImplant, *knee);

		auto trans = itk::VersorRigid3DTransform<>::New();
		trans->SetMatrix(femurMatch->GetRotationMatrix());
		trans->SetTranslation(femurMatch->GetTranslationMatrix());
		return trans;



	}
	void testImplants()
	{
		const char *dirPath = "D:\sovajo\Test_Cases\TKA_Test";

		auto femur10Implants = createFemurImplant(QString("%1/femur_right_10_data.json").arg(dirPath).toStdString().c_str(),
			QString("%1/femur_right_10.stl").arg(dirPath).toStdString().c_str());
		auto femur11Implants = createFemurImplant(QString("%1/femur_right_11_data.json").arg(dirPath).toStdString().c_str(),
			QString("%1/femur_right_11.stl").arg(dirPath).toStdString().c_str());
		auto tibia17_6Implants = createTibiaImplant(QString("%1/tibia_17_6mm_data.json").arg(dirPath).toStdString().c_str(),
			QString("%1/tibia_17_6mm.stl").arg(dirPath).toStdString().c_str());

		auto landmarks = readLandmarks(QString("%1/landmark.json").arg(dirPath).toStdString().c_str());
		//Right Leg
		auto knee = createKnee(landmarks, QString("%1/femur.vtk").arg(dirPath).toStdString().c_str(),
			QString("%1/tibia.vtk").arg(dirPath).toStdString().c_str(), false);
		//femur to 11 femurimplant

		auto implantToFemur = matchFemur(knee, femur11Implants);
		auto implantToTibia = matchTibia(knee, tibia17_6Implants);
		TKA::IMPLANTS::ImplantsMatchFinalInfo info(knee.get(), *femur11Implants, *tibia17_6Implants,
			implantToFemur, implantToTibia);
		//it is fine
		printInfo(info);

		//change femur implant from 11 to 10,keep angle and distance
		auto implantToFemur2 = matchFemur(knee, femur10Implants);
		//Code inside braces will cause problems
		{

			//the Two Inverse operations will cause NAN
			//VersorRigid3DTransform->vtkMatrix->Invert->vtkTransform->Inverse->VersorRigid3DTransform
			auto rotation = implantToFemur2->GetMatrix();
			auto translation = implantToFemur2->GetTranslation();
			vtkNew<vtkMatrix4x4> matrix;
			matrix->Identity();
			for (int row = 0; row < 3; row++) {
				for (int col = 0; col < 3; col++) {
					matrix->SetElement(row, col, rotation(row, col));
				}
				matrix->SetElement(row, 3, translation[row]);
			}
			//Mark: the Two Inverse operations will cause NAN
			matrix->Invert();

			auto boneToImplant = vtkSmartPointer<vtkTransform>::New();
			boneToImplant->SetMatrix(matrix);
			boneToImplant->Inverse();


			vtkNew<vtkMatrix4x4> matrix2;
			boneToImplant->GetMatrix(matrix2);
			itk::Rigid3DTransform<>::MatrixType m;
			itk::Rigid3DTransform<>::OffsetType offset;
			for (int row = 0; row < 3; row++)
			{
				for (int col = 0; col < 3; col++)
				{
					m(row, col) = matrix2->GetElement(row, col);
				}
				offset[row] = matrix2->GetElement(row, 3);
			}
			auto trans2 = itk::VersorRigid3DTransform<>::New();
			trans2->SetMatrix(m);
			trans2->SetOffset(offset);

			std::cout << "******equal:" << equal(implantToFemur2, trans2) << std::endl;
			implantToFemur2 = trans2;
		}
		//
		TKA::IMPLANTS::ImplantsMatchFinalInfo info2(knee.get(), *femur10Implants, *tibia17_6Implants,
			implantToFemur2, implantToTibia);
		//it is fine
		printInfo(info2);
		//keep angle and distance
		//The errors start here
		info2.setFemurVarusAngle(info.GetFemurVarusAngle());
		info2.setFemurPCAAngle(info.GetFemurImplantPCAAngle());
		info2.setFemurTEAAngle(info.GetFemurImplantTEAAngle());
		info2.setFemurFlexionAngle(info.GetFemurFlexionAngle());
		auto axialDistance = info.GetFemurResectionAxial();
		if (axialDistance.lateral < axialDistance.medial)
		{
			info2.SetThicknessFemurAxialMedial(axialDistance.medial);
		}
		else
		{
			info2.SetThicknessFemurAxialLateral(axialDistance.lateral);
		}
		auto cornalDistance = info.GetFemurResectionCoronal();
		if (cornalDistance.lateral < cornalDistance.medial)
		{
			info2.SetThicknessFemurCoronalMedial(cornalDistance.medial);
		}
		else
		{
			info2.SetThicknessFemurCoronalLateral(cornalDistance.lateral);
		}
		printInfo(info2);
	}

	
	void testTibiaImplant2()
	{
		const char *dirPath = "D:\\sovajo\\Test_Cases\\TKA_Test";

		auto femur4Implants = createFemurImplant(QString("%1/femur_right_4_data.json").arg(dirPath).toStdString().c_str(),
			QString("%1/femur_right_4.stl").arg(dirPath).toStdString().c_str());
		auto tibia3_10Implants = createTibiaImplant(QString("%1/tibia_right_3_10+_data.json").arg(dirPath).toStdString().c_str(),
			QString("%1/tibia_right_3_10+.stl").arg(dirPath).toStdString().c_str());
		auto tibia3_12Implants = createTibiaImplant(QString("%1/tibia_right_4_12+_data.json").arg(dirPath).toStdString().c_str(),
			QString("%1/tibia_right_3_12+.stl").arg(dirPath).toStdString().c_str());
		auto landmarks = readLandmarks(QString("%1/landmark.json").arg(dirPath).toStdString().c_str());
		//Right Leg
		auto knee = createKnee(landmarks, QString("%1/femur.vtk").arg(dirPath).toStdString().c_str(),
			QString("%1/tibia.vtk").arg(dirPath).toStdString().c_str(), false);
		auto implantToFemur = matchFemur(knee, femur4Implants);
		auto implantToTibia = matchTibia(knee, tibia3_10Implants);

		//1
		implantToTibia->Print(std::cout);

		TKA::IMPLANTS::ImplantsMatchFinalInfo info(&knee, *femur4Implants, *tibia3_10Implants,
			implantToFemur, implantToTibia);
		info.setTibiaVarusAngle(3.1562700920978872);
		info.setTibiaRotationAngle(-30.746164266024422);
		info.setTibiaSlopeAngle(-0.40464877902511232);
		auto newImplantToTibia = info.getITKTibiaTransform();
		//2
		newImplantToTibia->Print(std::cout);

	}
}
