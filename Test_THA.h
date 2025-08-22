#ifndef TEST_THA_SUEN_NEW_H
#define TEST_THA_SUEN_NEW_H

#include <QJsonObject>
#include <QJsonValue>
#include <iostream>
#include <QStringList>
#include "vtkMatrix4x4.h"
#include "vtkTransform.h"
#include "itkVersorRigid3DTransform.h"
#include "tha_implants/HipFemurOppside.hpp"
#include "tha_implants/HipPelvis.hpp"
#include "tha_implants/HipFemurStemImplant.hpp"
#include "tha_implants/HipFemurStemHeadImplant.hpp"
#include "tha_implants/HipPelvisCupImplant.hpp"
#include "tha_implants/HipFemurStemImplantMatch.hpp"
#include "tha_implants/HipPelvisImplantsMatchInfo.hpp"
#include "tha_implants/HipFemurStemImplantMatch.hpp"
#include "tha_implants/HipFemurStemHeadImplantMatch.hpp"
#include "tha_implants/HipPelvisCupImplantMatch.hpp"
#include "tha_implants/Point.hpp"
#include "vtkSTLReader.h"
#include "vtkPolyDataReader.h"
#include <QJsonDocument>
#include <QFile>
#include <QJsonArray>

namespace TEST_THA_SUEN
{
	enum  LandmarkType
	{
		//acetabulum
		kASISSurgicalSide = 0,
		kASISNonSurgicalSide,
		kAcetabulumFeaturePoint1,
		kAcetabulumFeaturePoint2,
		kAcetabulumFeaturePoint3,
		kSurgicalPubicTubercle,
		kNonSurgicalPubicTubercle,
		kAcetabulumRotateCenter = 99,//the landmarktype need 10 points to generate

		//surgical femur
		kFemurSurgicalFeaturePoint1 = 100,
		kFemurSurgicalLesserTrochanter,
		kFemurSurgicalHeadCenter,
		kFemurSurgicalFeaturePoint2,
		kFemurSurgicalFeaturePoint3,
		kMedialEpicondyle,
		kLaterialEpicondyle,
		kFemurSurgicalKneeCenter,
		KFemurSurgicalDistalCanalCenter,
		kFemurSurgicalNeckCenter,
		kFemurSurgicalProximalCanalCenter,
		//None-surgical Femur

		kFemurNonSurgicalLesserTrochanter,
		kFemurNonSurigicalHeadCenter,
		kFemurNonSurgicalKneeCenter,
		KFemurNonSurgicalDistalCanalCenter,
		kFemurNonSurgicalProximalCanalCenter,
		kNone

	};

	vtkSmartPointer<vtkPolyData> readSTL(const QString & stlPath)
	{
		vtkNew<vtkSTLReader> reader;
		reader->SetFileName(stlPath.toStdString().c_str());
		reader->Update();
		return reader->GetOutput();
	}
	vtkSmartPointer<vtkPolyData> readVTK(const QString & stlPath)
	{
		vtkNew<vtkPolyDataReader> reader;
		reader->SetFileName(stlPath.toStdString().c_str());
		reader->Update();
		return reader->GetOutput();
	}

	THA::IMPLANTS::Point readPoint(const QJsonValue& val)
	{
		auto arr = val.toArray();
		return THA::IMPLANTS::Point(arr[0].toDouble(), arr[1].toDouble(), arr[2].toDouble());
	}

	std::shared_ptr<THA::IMPLANTS::HipPelvisCupImplant> createCupImplant(const QString & jsonPath)
	{
		QFile file(jsonPath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();

		auto cupImplant = std::make_shared<THA::IMPLANTS::HipPelvisCupImplant>();

		cupImplant->init(readPoint(obj["top"]), readPoint(obj["p1"]), readPoint(obj["p2"]), readPoint(obj["p3"]),
			readPoint(obj["center"]), obj["diameter"].toDouble()*0.5,
			10.0 / VTK_INT_MAX, obj["thickness"].toDouble());
		return cupImplant;
	}

	std::shared_ptr<THA::IMPLANTS::HipFemurStemImplant> createHipFemurStemImplant(const QString & jsonPath)
	{
		QFile file(jsonPath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();
		auto surfaceArr = obj["head_plane"].toArray();
		std::vector<THA::IMPLANTS::Point> headPlanePoints;
		for (auto p : surfaceArr)
		{
			headPlanePoints.push_back(readPoint(p));
		}
		auto stemImplant = std::make_shared<THA::IMPLANTS::HipFemurStemImplant>();
		stemImplant->init(readPoint(obj["rod_center"]), readPoint(obj["rod_base"]),
			readPoint(obj["head_center"]), headPlanePoints);

		return stemImplant;
	}

	std::shared_ptr<THA::IMPLANTS::HipFemurStemHeadImplant> createHipFemurHeadImplant(const QString & jsonPath)
	{
		QFile file(jsonPath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();
		auto headImplant = std::make_shared<THA::IMPLANTS::HipFemurStemHeadImplant>();
		headImplant->init(readPoint(obj["p1"]), readPoint(obj["p2"]),
			readPoint(obj["p3"]), readPoint(obj["top"]),
			readPoint(obj["center"]),
			obj["diameter"].toDouble());
		return headImplant;
	}
	std::map<LandmarkType, THA::IMPLANTS::Point> readLandmarks(const QString & filePath)
	{
		std::map<LandmarkType, THA::IMPLANTS::Point> items;
		QFile file(filePath);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read landmarks failed" << std::endl;
			return items;
		}
		auto obj = QJsonDocument::fromJson(file.readAll()).object();
		for (auto key : obj.keys())
		{
			auto landmarkType = (LandmarkType)key.toInt();
			if (landmarkType == kAcetabulumRotateCenter)
			{
				std::vector<THA::IMPLANTS::Point> points;
				for (auto item : obj[key].toArray())
				{
					points.push_back(readPoint(item));
				}
				auto center = THA::IMPLANTS::HipPelvis::getNativeCenterOfRotation(points);
				items[landmarkType] = THA::IMPLANTS::Point(center.first.x, center.first.y, center.first.z);
			}
			else
			{
				items[landmarkType] = readPoint(obj[key]);
			}
		}
		return items;
	}

	std::shared_ptr<THA::IMPLANTS::HipPelvis> createHipPelvis(const QString &dir)
	{

		auto pelvisSide = THA::IMPLANTS::PelvisSide::LEFT_SIDE;
		auto landmarks = readLandmarks(QString("%1/landmark.json").arg(dir));
		auto leftFemurData = readVTK(QString("%1/left_femur.vtk").arg(dir));
		auto rightFemurData = readVTK(QString("%1/right_femur.vtk").arg(dir));
		auto pevlisData = readVTK(QString("%1/pelvis.vtk").arg(dir));
		THA::IMPLANTS::HipFemur  hipFemur;

		hipFemur.init(landmarks[LandmarkType::kFemurSurgicalHeadCenter], landmarks[LandmarkType::kFemurSurgicalNeckCenter],
			landmarks[LandmarkType::KFemurSurgicalDistalCanalCenter],
			landmarks[LandmarkType::kFemurSurgicalProximalCanalCenter], landmarks[LandmarkType::kFemurSurgicalLesserTrochanter],
			landmarks[LandmarkType::kMedialEpicondyle], landmarks[LandmarkType::kLaterialEpicondyle],
			landmarks[LandmarkType::kFemurSurgicalKneeCenter],
			leftFemurData);



		THA::IMPLANTS::HipFemurOppside hipFemurOppside;
		hipFemurOppside.init(landmarks[LandmarkType::kFemurNonSurigicalHeadCenter], landmarks[LandmarkType::kFemurNonSurgicalProximalCanalCenter],
			landmarks[LandmarkType::kFemurNonSurgicalLesserTrochanter], landmarks[LandmarkType::kFemurNonSurgicalKneeCenter]);

		auto hipPelvis = std::make_shared<THA::IMPLANTS::HipPelvis>();

		THA::IMPLANTS::Plane coronalCT;
		coronalCT.init(THA::IMPLANTS::Point(0, 1, 0), THA::IMPLANTS::Point(0, 0, 0));

		hipPelvis->init(landmarks[LandmarkType::kASISSurgicalSide], landmarks[LandmarkType::kASISNonSurgicalSide],
			landmarks[LandmarkType::kSurgicalPubicTubercle], landmarks[LandmarkType::kNonSurgicalPubicTubercle],
			pevlisData, hipFemur,
			hipFemurOppside,
			pelvisSide, landmarks[kAcetabulumRotateCenter], coronalCT);

		return hipPelvis;
	}

	itk::Rigid3DTransform<>::Pointer matchCupToAcetabulum(std::shared_ptr<THA::IMPLANTS::HipPelvisCupImplant> cupImplant,
		std::shared_ptr<THA::IMPLANTS::HipPelvis> hipPelvis, double pAbductionAngle, double pAnteversionAngle)
	{

		THA::IMPLANTS::HipPelvisCupImplantMatch cupMatch;
		cupMatch.init(*hipPelvis, *cupImplant);
		auto matrix = cupMatch.getTransform(pAbductionAngle, pAnteversionAngle, 0, 0, 0);
		return matrix;
	}

	itk::Rigid3DTransform<>::Pointer matchStemToFemur(std::shared_ptr<THA::IMPLANTS::HipFemurStemImplant> stemImplant,
		std::shared_ptr<THA::IMPLANTS::HipPelvis> hipPelvis)
	{

		THA::IMPLANTS::HipFemurStemImplantMatch stemMatch;
		stemMatch.init(*hipPelvis, *stemImplant);
		auto trans = itk::VersorRigid3DTransform<>::New();
		trans->SetMatrix(stemMatch.GetRotationMatrix());
		trans->SetOffset(stemMatch.GetTranslationMatrix());
		return trans;

	}

	itk::Rigid3DTransform<>::Pointer matchHeadToStem(std::shared_ptr<THA::IMPLANTS::HipFemurStemImplant> stemImplant,
		std::shared_ptr<THA::IMPLANTS::HipFemurStemHeadImplant> headImplant)
	{

		THA::IMPLANTS::HipFemurStemHeadImplantMatch headMatch;
		headMatch.init(*headImplant, *stemImplant);
		return headMatch.getStemHeadTransform();
	}

	void testImplant()
	{
		QString dataDir = "E:/data/tha/implant_test";
		auto hipPelvis = createHipPelvis(dataDir);
		auto cupImplant = createCupImplant(QString("%1/cup48mm.json").arg(dataDir));
		auto stemImplant = createHipFemurStemImplant(QString("%1/stem3#.json").arg(dataDir));
		auto headImplant = createHipFemurHeadImplant(QString("%1/CEH32mm_M.json").arg(dataDir));
		auto cupToAcetabulumTransform = matchCupToAcetabulum(cupImplant, hipPelvis, 40, 20);
		auto stemToFemurTransform = matchStemToFemur(stemImplant, hipPelvis);
		auto headToStemTransform = matchHeadToStem(stemImplant, headImplant);
		THA::IMPLANTS::HipPelvisImplantsMatchInfo matchInfo(*hipPelvis,
			*cupImplant, *stemImplant, *headImplant,
			cupToAcetabulumTransform, stemToFemurTransform, headToStemTransform);
		auto defaultTiltAngle = hipPelvis->getCoronalTiltAngleDegree();
		std::cout << "pelvis tilt angle:" << defaultTiltAngle << std::endl;
		matchInfo.setPelvisTiltAngleDegree(defaultTiltAngle + 1);
		auto newCupToAcetabulumTransform = matchInfo.setCupAngles(40, 20);

	}
	void testImplant1()
	{
		QString dataDir = "D:\\sovajo\\Test_Cases\\octubre_2025\\Pelvis_Test_1";
		auto hipPelvis = createHipPelvis(dataDir);
		auto hipLenth = hipPelvis->getHipLengthDistance();
		auto oppsiteLength = hipPelvis->getHipLengthDistanceOppsite();
		auto offset = hipPelvis->getCombinedOffsetDistance();
		auto oppsiteOffset = hipPelvis->getCombinedOffsetDistanceOppsite();
		std::cout << hipLenth << "," << oppsiteLength << std::endl;
		std::cout << offset << "," << oppsiteOffset << std::endl;
	}

	void testImplant2()
	{
		QString dataDir = "D:\\sovajo\\Test_Cases\\octubre_2025\\pelvis_test_2";
		auto hipPelvis = createHipPelvis(dataDir);
		auto hipLenth = hipPelvis->getHipLengthDistance();
		auto oppsiteLength = hipPelvis->getHipLengthDistanceOppsite();
		auto offset = hipPelvis->getCombinedOffsetDistance();
		auto oppsiteOffset = hipPelvis->getCombinedOffsetDistanceOppsite();
		//!!!
		//hip length and combined offset are wrong
		std::cout << hipLenth << "," << oppsiteLength << std::endl;
		std::cout << offset << "," << oppsiteOffset << std::endl;
		//!!! 
		//This pelvis has a positive inclination angle.but the angle from getCoronalTiltAngleDegree is Negative 
		std::cout << hipPelvis->getCoronalTiltAngleDegree() << std::endl;

		//
		auto cupImplant = createCupImplant(QString("%1/DTUCA56mm.json").arg(dataDir));
		auto stemImplant = createHipFemurStemImplant(QString("%1/DTUCS8#.json").arg(dataDir));
		auto headImplant = createHipFemurHeadImplant(QString("%1/CEH36mm_S.json").arg(dataDir));
		auto cupToAcetabulumTransform = matchCupToAcetabulum(cupImplant, hipPelvis, 40, 20);
		auto stemToFemurTransform = matchStemToFemur(stemImplant, hipPelvis);
		auto headToStemTransform = matchHeadToStem(stemImplant, headImplant);
		THA::IMPLANTS::HipPelvisImplantsMatchInfo matchInfo(*hipPelvis,
			*cupImplant, *stemImplant, *headImplant,
			cupToAcetabulumTransform, stemToFemurTransform, headToStemTransform);

		//!!! we need radiographic version,not true anatomic version  measured relative to the pelvic frontal plane
		//!!! we need radiographic inclination,not true anatomic version measured relative to the pelvic frontal plane
		std::cout << "cupVersion:" << matchInfo.getCupAntversion() << std::endl;//20
		std::cout << "cupInclination:" << matchInfo.getCupInclination() << std::endl;//40

		THA::IMPLANTS::Plane sagital, coronal;
		sagital.init(THA::IMPLANTS::Point(1, 0, 0), THA::IMPLANTS::Point(0, 0, 0));
		coronal.init(THA::IMPLANTS::Point(0, 1, 0), THA::IMPLANTS::Point(0, 0, 0));

		std::cout << "cupVersion CT:" << matchInfo.getCupAntversion(sagital, coronal) << std::endl;//20
		std::cout << "cupInclination CT:" << matchInfo.getCupInclination(sagital, coronal) << std::endl;//40

		{
			//!The anterverison and inclination angle i calculated are completely different from the angles from HipPelvisImplantsMatchInfo.
			auto cupPlaneNormalVector = cupToAcetabulumTransform->TransformPoint(cupImplant->getVectorZ().ToITKPoint());
			double orgin[3] = { 0,0,0 };
			double saggital[3] = { 1,0,0 };
			double coronal[3] = { 0,1,0 };
			double z[3] = { 0,0,1 };
			double xproj[3];

			// project the normal vector of the cup plane  onto the sagittal plane, and calculate the cup anterversion angle 
			vtkPlane::ProjectVector(cupPlaneNormalVector.GetDataPointer(), orgin, saggital, xproj);
			auto cupAnteverion = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(xproj, z));
			vtkPlane::ProjectVector(cupPlaneNormalVector.GetDataPointer(), orgin, coronal, xproj);
			// project the normal vector of the cup plane  onto the cornal plane, and calculate the  cup inclination angle 
			auto cupInclination = vtkMath::DegreesFromRadians(vtkMath::AngleBetweenVectors(xproj, z));
			std::cout << "cupVersion2:" << cupAnteverion << std::endl;//0.2
			std::cout << "cupInclination2:" << cupInclination << std::endl;//10.9

		}
	}
}

#endif