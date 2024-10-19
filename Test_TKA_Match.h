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

namespace TEST_IMPLANTS
{
	std::string casePlanPathFemur = "D:\\sovajo\\Errores\\Error9\\test_implants\\fine_case_plan\\femur.vtk"; 
	std::string casePlanPathTibia = "D:\\sovajo\\Errores\\Error9\\test_implants\\fine_case_plan\\tibia.vtk";
	std::string casePlanPathLandmark = "D:\\sovajo\\Errores\\Error9\\test_implants\\fine_case_plan\\landmark.json";

	std::string implantsPathTibia10 = "D:\\sovajo\\Errores\\Error9\\test_implants\\implants\\tibia_10_6_data.json";
	std::string implantsPathFemur10 = "D:\\sovajo\\Errores\\Error9\\test_implants\\implants\\femur_left_10_data.json";
	std::string implantsPathFemur18 = "D:\\sovajo\\Errores\\Error9\\test_implants\\implants\\femur_left_18_data.json";

	enum Landmarks {
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

	std::map<Landmarks, itk::Point<double>> readCtLandmarks(const char *filePath)
	{
		std::map<Landmarks, itk::Point<double>> ctLanmarks;
		QFile file(filePath);
		file.open(QIODevice::ReadOnly);
		auto obj = QJsonDocument::fromJson(file.readAll()).object();
		auto keys = obj.keys();
		for (int i = 0; i < keys.size(); ++i)
		{
			auto val = obj[keys[i]].toArray();
			itk::Point<double> p;
			p[0] = val[0].toDouble();
			p[1] = val[1].toDouble();
			p[2] = val[2].toDouble();
			ctLanmarks[(Landmarks)keys[i].toInt()] = p;
		}
		return ctLanmarks;
	}

	vtkSmartPointer<vtkPolyData> readVtkFile(const char *filePath)
	{
		vtkNew<vtkPolyDataReader> reader;
		reader->SetFileName(filePath);
		reader->Update();
		return reader->GetOutput();
	}
	vtkSmartPointer<vtkPolyData> readSTLFile(const char *filePath)
	{
		vtkNew<vtkSTLReader> reader;
		reader->SetFileName(filePath);
		reader->Update();
		return reader->GetOutput();
	}

	class LandmarksFind
	{
	public:
		LandmarksFind(std::map<Landmarks, itk::Point<double>> &landmarks) :m_landmarks(landmarks)
		{
		}
		Point operator[](Landmarks type)
		{
			auto p = m_landmarks[type];
			return Point(p[0], p[1], p[2]);
		}
	private:
		std::map<Landmarks, itk::Point<double>> m_landmarks;
	};

	enum BoundingPlaneType {
		kPlaneFemurPosterior,         // A
		kPlaneFemurPosteriorOblique,  // B
		kPlaneFemurDistal,            // C
		kPlaneFemurAnteriorOblique,   // D
		kPlaneFemurAnterior,          // E
		kPlaneTibia,
		kPlanePatella,
	};


	std::shared_ptr<TibiaImplant> createTibiaImplant(const char *jsonFile)
	{
		QFile file(jsonFile);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed:" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll().toLower()).object();

		std::map<QString, Point> tibiaKeyPoints;
		for (auto key : { "point_pcl1", "point_pcl2", "point_front", "point_external", "point_plateau_left", "point_plateau_right" }) {
			auto val = obj[key].toArray();
			tibiaKeyPoints[key] = Point(val[0].toDouble(), val[1].toDouble(), val[2].toDouble());
		}
		TibiaImplantInfo tibiaInfo;
		tibiaInfo.tibiaLateralThickness = obj["thickness_left"].toDouble();
		tibiaInfo.tibiaMedialThickness = obj["thickness_right"].toDouble();
		auto tibiaImplant = std::make_shared<TibiaImplant>();
		tibiaImplant->init(tibiaKeyPoints["point_pcl1"], tibiaKeyPoints["point_pcl2"], tibiaKeyPoints["point_front"], tibiaKeyPoints["point_external"], tibiaInfo);
		return tibiaImplant;
	}

	std::shared_ptr<FemurImplant> createFemurImplant(const char * jsonFile)
	{
		QFile file(jsonFile);
		if (!file.open(QIODevice::ReadOnly))
		{
			std::cout << "read failed:" << std::endl;
			return nullptr;
		}
		auto obj = QJsonDocument::fromJson(file.readAll().toLower()).object();
		FemurImplantInfo  femurInfo;
		femurInfo.femurPosteriorLateralThickness = obj["posterior_lateral_thickness"].toDouble();;
		femurInfo.femurDistalLateralThickness = obj["distal_lateral_thickness"].toDouble();
		femurInfo.femurPosteriorMedialThickness = obj["posterior_medial_thickness"].toDouble();
		femurInfo.femurDistalMedialThickness = obj["distal_medial_thickness"].toDouble();

		std::map<BoundingPlaneType, Plane> planes;
		for (int faceIdx = 0; faceIdx < 5; faceIdx++) {
			auto val = obj[QString("face_") + static_cast<char>('a' + faceIdx)];
			std::array<double, 6> res;
			auto arr = val.toArray();
			for (int i = 0; i < 6; i++) {
				res[i] = arr[i].toDouble();
			}
			Point origP(res[3], res[4], res[5]);
			Point normP(res[0], res[1], res[2]);
			Plane plane;
			plane.init(normP, origP);
			planes[(BoundingPlaneType)faceIdx] = plane;
		}
		std::map<QString, Point> P;
		for (auto key : { "p1", "p2", "p3" }) {
			auto arr = obj[key].toArray();
			P[key] = Point(arr[0].toDouble(), arr[1].toDouble(), arr[2].toDouble());

		}
		QFileInfo fileInfo(jsonFile);
		auto stlPath = fileInfo.absoluteDir().absolutePath() + "/" + obj["filename"].toString();
		auto implantData = readSTLFile(stlPath.toStdString().c_str());
		auto femurImplant = std::make_shared<FemurImplant>();
		std::vector<itk::Point<double>> res;

		{
			auto arrList = obj["patella"].toArray();
			for (auto point : arrList) {
				auto arr = point.toArray();
				itk::Point<double> p;
				p[0] = arr[1].toDouble();
				p[1] = arr[2].toDouble();
				p[2] = arr[3].toDouble();
				res.push_back(p);
			}
		}
		femurImplant->init(
			planes[kPlaneFemurPosterior],
			planes[kPlaneFemurPosteriorOblique],
			planes[kPlaneFemurDistal],
			planes[kPlaneFemurAnteriorOblique],
			planes[kPlaneFemurAnterior],
			P["p1"], P["p2"], P["p3"],
			implantData,
			res, femurInfo
		);


		return femurImplant;
	}


	itk::Rigid3DTransform<>::Pointer getTransform(const itk::Matrix<double>& rotation, const itk::Vector<double>& translation)
	{
		auto res = itk::VersorRigid3DTransform<>::New();
		res->SetMatrix(rotation);
		res->SetTranslation(translation);
		return res.GetPointer();


	}

	void testImplant()
	{
		auto femurData = readVtkFile(casePlanPathFemur.c_str());
		auto tibiaData = readVtkFile(casePlanPathTibia.c_str());
		auto landmarks = readCtLandmarks(casePlanPathLandmark.c_str());
		LandmarksFind anno(landmarks);
		auto ankleCenter = Knee::getComputeAnkleCenter(anno[Landmarks::kLateralMalleolus], anno[Landmarks::kMedialMalleolus]);
		Knee knee;
		Patella patella;

		knee.init(anno[Landmarks::kHipCenter],
			anno[Landmarks::kFemurAnteriorCortex],
			anno[Landmarks::kFemurKneeCenter],
			anno[Landmarks::kLateralEpicondyle],
			anno[Landmarks::kMedialEpicondyle],
			anno[Landmarks::kTibiaKneeCenter],
			anno[Landmarks::kTibiaTuberosity],
			anno[Landmarks::kPCLInsertionPoint], ankleCenter, patella,
			femurData, tibiaData, KneeSideEnum::KRight
		);
		knee.setLateralAndMedialInferiorFemurPoints(anno[Landmarks::kFemurDistalLateral], anno[Landmarks::kFemurDistalMedial]);
		knee.setLateralAndMedialPosteriorFemurPoints(anno[Landmarks::kFemurLateralPosteriorCondyle], anno[Landmarks::kFemurMedialPosteriorCondyle]);
		knee.setLateralAndMedialPlateauPoints(anno[Landmarks::kTibiaLateralPlatformPoint], anno[Landmarks::kTibiaMedialPlatformPoint]);

		//
		auto tibiaImplant10 = createTibiaImplant(implantsPathTibia10.c_str());
		auto femurImplant10 = createFemurImplant(implantsPathFemur10.c_str());
		auto femurImplant18 = createFemurImplant(implantsPathFemur18.c_str());

		FemurImplantMatch matchFemur10;
		matchFemur10.init(*femurImplant10, knee);
		auto implant10Tofemur = getTransform(matchFemur10.GetRotationMatrix(), matchFemur10.GetTranslationMatrix());

		TibiaImplantMatch matchTibai10;
		matchTibai10.init(*tibiaImplant10, knee);
		auto implant10Totibia = getTransform(matchTibai10.GetRotationMatrix(), matchTibai10.GetTranslationMatrix());


		FemurImplantMatch matchFemur18;
		matchFemur18.init(*femurImplant18, knee);
		auto implant18Tofemur = getTransform(matchFemur18.GetRotationMatrix(), matchFemur18.GetTranslationMatrix());


		//implantFinalInfo with femur implant 10
		ImplantsMatchFinalInfo implantFinalInfo10(&knee, *femurImplant10, *tibiaImplant10, implant10Tofemur, implant10Totibia);

		auto femurVarus = implantFinalInfo10.GetFemurVarusAngle();
		auto pcaAngle = implantFinalInfo10.GetFemurImplantPCAAngle();
		auto teaAngle = implantFinalInfo10.GetFemurImplantTEAAngle();
		auto flexionAngle = implantFinalInfo10.GetFemurFlexionAngle();
		auto femurResectionAxial = implantFinalInfo10.GetFemurResectionAxial();
		auto femurResectionCoronal = implantFinalInfo10.GetFemurResectionCoronal();
		std::cout << "----------------test0-----------------------" << std::endl;
		std::cout << "femurVarus=" << femurVarus << std::endl;
		std::cout << "pcaAngle=" << pcaAngle << std::endl;
		std::cout << "teaAngle=" << teaAngle << std::endl;
		std::cout << "flexionAngle=" << flexionAngle << std::endl;
		std::cout << "femurResectionAxial=" << femurResectionAxial.lateral << "," << femurResectionAxial.medial << std::endl;
		std::cout << "femurResectionCoronal=" << femurResectionCoronal.lateral << "," << femurResectionCoronal.medial << std::endl;
		std::cout << "Tibia varus = " << implantFinalInfo10.GetTibiaVarusAngle() << std::endl;
		std::cout << "Tibia slope = " << implantFinalInfo10.GetTibiaImplantSlopeAngle() << std::endl;
		std::cout << "Tibia rotations = " << implantFinalInfo10.GetTibiaImplantRotationAngle() << std::endl;
		
		std::cout << "-----------------------------------------" << std::endl;

		//set angles
		/*implantFinalInfo10.setFemurFlexionAngle(5);
		implantFinalInfo10.setFemurTEAAngle(15);
		implantFinalInfo10.setFemurVarusAngle(-6);*/
		//implantFinalInfo10.setFemurVarusAngle(6);
		//implantFinalInfo10.setFemurTEAAngle(-4);
		//implantFinalInfo10.setFemurTEAAngle(-4);
		
		implantFinalInfo10.setFemurVarusAngle(-6);
		implantFinalInfo10.setFemurFlexionAngle(5);

		//implantFinalInfo10.setFemurPCAAngle(-9);
		//implantFinalInfo10.setFemurTEAAngle(5);

		//Point cross = (-knee.getFemurVectorLateralTEA()).cross(knee.getDirectVectorFemurAxis());
		//cross.normalice();
		//Point myCenter = knee.getFemurKneeCenter() + 40 * cross;
		//std::vector<Point> lista = { myCenter };
		//ImplantTools::show(knee.GetFemurPoly(), lista);

		//implantFinalInfo10.setTibiaSlopeAngle(4);
		//implantFinalInfo10.setTibiaRotationAngle(6);
		implantFinalInfo10.setTibiaSlopeAngle(-3);

		implantFinalInfo10.setTibiaVarusAngle(4);
		//implantFinalInfo10.setTibiaRotationAngle(6);
		//implantFinalInfo10.setTibiaVarusAngle(-5);
		

		femurVarus = implantFinalInfo10.GetFemurVarusAngle();
		pcaAngle = implantFinalInfo10.GetFemurImplantPCAAngle();
		teaAngle = implantFinalInfo10.GetFemurImplantTEAAngle();
		flexionAngle = implantFinalInfo10.GetFemurFlexionAngle();
		femurResectionAxial = implantFinalInfo10.GetFemurResectionAxial();
		femurResectionCoronal = implantFinalInfo10.GetFemurResectionCoronal();
		std::cout << "----------------test1-----------------------" << std::endl;
		std::cout << "femurVarus=" << femurVarus << std::endl;
		std::cout << "pcaAngle=" << pcaAngle << std::endl;
		std::cout << "teaAngle=" << teaAngle << std::endl;
		std::cout << "flexionAngle=" << flexionAngle << std::endl;
		std::cout << "femurResectionAxial=" << femurResectionAxial.lateral << "," << femurResectionAxial.medial << std::endl;
		std::cout << "femurResectionCoronal=" << femurResectionCoronal.lateral << "," << femurResectionCoronal.medial << std::endl;
		std::cout << "Tibia varus = " << implantFinalInfo10.GetTibiaVarusAngle() << std::endl;
		std::cout << "Tibia slope = " << implantFinalInfo10.GetTibiaImplantSlopeAngle() << std::endl;
		std::cout << "Tibia rotations = " << implantFinalInfo10.GetTibiaImplantRotationAngle() << std::endl;
		
		std::cout << "-----------------------------------------" << std::endl;

		//implantFinalInfo with femur implant 18

		ImplantsMatchFinalInfo implantFinalInfo18(&knee, *femurImplant18, *tibiaImplant10, implant18Tofemur, implant10Totibia);
		implantFinalInfo18.setFemurVarusAngle(femurVarus);
		implantFinalInfo18.setFemurPCAAngle(pcaAngle);
		implantFinalInfo18.setFemurTEAAngle(teaAngle);
		implantFinalInfo18.setFemurFlexionAngle(flexionAngle);
		implantFinalInfo18.SetThicknessFemurAxialMedial(femurResectionAxial.medial);
		implantFinalInfo18.SetThicknessFemurAxialLateral(femurResectionAxial.lateral);
		implantFinalInfo18.SetThicknessFemurCoronalMedial(femurResectionCoronal.medial);
		implantFinalInfo18.SetThicknessFemurCoronalLateral(femurResectionCoronal.lateral);

		femurVarus = implantFinalInfo18.GetFemurVarusAngle();
		pcaAngle = implantFinalInfo18.GetFemurImplantPCAAngle();
		teaAngle = implantFinalInfo18.GetFemurImplantTEAAngle();
		flexionAngle = implantFinalInfo18.GetFemurFlexionAngle();
		femurResectionAxial = implantFinalInfo18.GetFemurResectionAxial();
		femurResectionCoronal = implantFinalInfo18.GetFemurResectionCoronal();

		std::cout << "----------------------test2----------------------" << std::endl;
		std::cout << "femurVarus=" << femurVarus << std::endl;
		std::cout << "pcaAngle=" << pcaAngle << std::endl;
		std::cout << "teaAngle=" << teaAngle << std::endl;
		std::cout << "flexionAngle=" << flexionAngle << std::endl;
		std::cout << "femurResectionAxial=" << femurResectionAxial.lateral << "," << femurResectionAxial.medial << std::endl;
		std::cout << "femurResectionCoronal=" << femurResectionCoronal.lateral << "," << femurResectionCoronal.medial << std::endl;
		std::cout << "Tibia varus = " << implantFinalInfo18.GetTibiaVarusAngle() << std::endl;
		std::cout << "Tibia slope = " << implantFinalInfo18.GetTibiaImplantSlopeAngle() << std::endl;
		std::cout << "Tibia rotations = " << implantFinalInfo18.GetTibiaImplantRotationAngle() << std::endl;
		
		std::cout << "----------------------------------------------" << std::endl;


		//get new transform
		auto newImplant18ToFemur = implantFinalInfo18.getITKFemurTransform();

		ImplantsMatchFinalInfo implantFinalInfo18_2(&knee, *femurImplant18, *tibiaImplant10, newImplant18ToFemur, implant10Totibia);
		femurVarus = implantFinalInfo18_2.GetFemurVarusAngle();
		pcaAngle = implantFinalInfo18_2.GetFemurImplantPCAAngle();
		teaAngle = implantFinalInfo18_2.GetFemurImplantTEAAngle();
		flexionAngle = implantFinalInfo18_2.GetFemurFlexionAngle();
		femurResectionAxial = implantFinalInfo18_2.GetFemurResectionAxial();
		femurResectionCoronal = implantFinalInfo18_2.GetFemurResectionCoronal();
		std::cout << "---------------------------test3-------------------" << std::endl;

		std::cout << "femurVarus=" << femurVarus << std::endl;
		std::cout << "pcaAngle=" << pcaAngle << std::endl;
		std::cout << "teaAngle=" << teaAngle << std::endl;
		std::cout << "flexionAngle=" << flexionAngle << std::endl;
		std::cout << "femurResectionAxial=" << femurResectionAxial.lateral << "," << femurResectionAxial.medial << std::endl;
		std::cout << "femurResectionCoronal=" << femurResectionCoronal.lateral << "," << femurResectionCoronal.medial << std::endl;
		std::cout << "Tibia varus = " << implantFinalInfo18_2.GetTibiaVarusAngle() << std::endl;
		std::cout << "Tibia slope = " << implantFinalInfo18_2.GetTibiaImplantSlopeAngle() << std::endl;
		std::cout << "Tibia rotations = " << implantFinalInfo18_2.GetTibiaImplantRotationAngle() << std::endl;
		
		std::cout << "----------------------------------------------" << std::endl;

	}
}


namespace TEST_PKA_SUEN
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
		auto arr = obj["patella"].toArray();
		for (auto item : arr)
		{
			itk::Point<double> p;
			p[0] = item.toArray()[0].toDouble();
			p[1] = item.toArray()[0].toDouble();
			p[2] = item.toArray()[0].toDouble();
			patellaPoints.push_back(p);
		}
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

	itk::Rigid3DTransform<>::Pointer matchTibia(std::shared_ptr<TKA::IMPLANTS::Knee> knee, std::shared_ptr<TKA::IMPLANTS::TibiaImplant> tibiaImplant)
	{
		auto tibiaMatch = std::make_shared<TKA::IMPLANTS::TibiaImplantMatch>();
		tibiaMatch->init(*tibiaImplant, *knee);

		auto trans = itk::VersorRigid3DTransform<>::New();
		trans->SetMatrix(tibiaMatch->GetRotationMatrix());
		trans->SetTranslation(tibiaMatch->GetTranslationMatrix());
		return trans;
	}

	itk::Rigid3DTransform<>::Pointer matchFemur(std::shared_ptr<TKA::IMPLANTS::Knee> knee, std::shared_ptr<TKA::IMPLANTS::FemurImplant> femurImplant)
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
		const char *dirPath = "D:/sovajo/Test_Cases/implants_change";

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

		{
			//Code inside braces will cause problems
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

			//std::cout << "******equal:" << equal(implantToFemur2, trans2) << std::endl;
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

}
