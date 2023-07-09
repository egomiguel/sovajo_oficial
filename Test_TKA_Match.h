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
