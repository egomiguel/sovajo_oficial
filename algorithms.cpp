#include <iostream>
//#include "test/TestFunction.hpp"
//#include "Test_TKA_Tibia_Rotation.h"
#include "test/TestVtk.hpp"
#include "Test_PKA.h"
#include "Test_PKA_SUEN.h"
#include "test/segmentationTest.cpp"

#include "implants/FemurImplantMatch.hpp"
#include "implants/TibiaImplantMatch.hpp"
#include "implants/PatellaImplantMatch.hpp"
#include "implants/LegAngle.hpp"
#include "implants/Balance.hpp"
#include "implants/ImplantsException.hpp"
//#include "implants/HipPoints.hpp"
#include "ImplantsMatchFinalInfo.hpp"
#include "implants/ImplantTools.hpp"
#include "implants/Patella.hpp"
#include "implants/PatellaImplantMatchInfo.hpp"

#include "tha_implants/HipPelvis.hpp"
#include "tha_implants/HipPelvisCupImplant.hpp"
#include "tha_implants/HipPelvisCupImplantMatch.hpp"

#include "tha_implants/HipFemurStemImplant.hpp"
#include "tha_implants/HipFemurStemImplantMatch.hpp"
#include "tha_implants/HipPelvisImplantsMatchInfo.hpp"
#include "tha_implants/HipFemurStemHeadImplantMatch.hpp"

#include "tha_implants/HipFemur.hpp"
#include "tha_implants/HipFemurOppside.hpp"
#include "tha_implants/Point.hpp"

#include "uka_implants/Knee.hpp"
#include "uka_implants/FemurImplantMatch.hpp"
#include "uka_implants/TibiaImplantMatch.hpp"
#include "uka_implants/Point.hpp"
#include "uka_implants/TibiaSpacerImplant.hpp"
#include "uka_implants/TibiaSpacerImplantMatch.hpp"
#include "uka_implants/ImplantsMatchFinalInfo.hpp"

#include "hip/HipCenter.hpp"
#include "segmentation/AutomaticSegmentation.hpp"
#include "segmentation/ManualSegmentation.hpp"
#include "segmentation/SliceBorder.hpp"
#include "segmentation/SPlane.hpp"
#include "segmentation/MakeSurfaceTest.cpp"
#include "registration/Registration.hpp"
#include "registration/FemurRegistration.hpp"
#include "registration/TibiaRegistration.hpp"
#include "registration/RPlane.hpp"
#include "registration/RegistrationException.hpp"
#include "registration/KneeCapRegistration.hpp"
#include "general_registration/GeneralRegistrationPointsToPolydata.hpp"
#include "general_registration/Types.hpp"
//#include "registration/HipPelvisRegistration.hpp"
//#include "registration/HipFemurRegistration.hpp"
#include "registration/LeastSquaresICP.hpp"
#include "registration/GeneralRegistrationPointsToPolydata.hpp"
//#include <itkViewImage.h>
#include "implants/Point.hpp"
#include "vtkAppendPolyData.h"
#include "vtkRuledSurfaceFilter.h"
#include "vtkContourTriangulator.h"
#include <cmath>
#include <vtkVersion.h>

//#include "slicer/vtkPlanarContourToClosedSurfaceConversionRule.h"
#include "segmentation/surface/BinaryTree.hpp"
#include "segmentation/surface/Tiling.hpp"
#include "vtkPointsProjectedHull.h"
#include "vtkKdTree.h"
#include "vtkConnectivityFilter.h"
#include "vtkElevationFilter.h"
#include "vtkWindowedSincPolyDataFilter.h"
#include "vtkSmoothPolyDataFilter.h"
#include "registration/LeastSquaresICP.hpp"
//#include "registration/CoherentPoint/CoherentPointDrift.hpp"
#include "implants/BoneRotulaPath.hpp"
#include "itkResampleImageFilter.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkTranslationTransform.h"
#include "itkAffineTransform.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkChangeInformationImageFilter.h"
#include "uka_implants/FindRegistrationPoints.hpp"
//#include "segmentation/coherent/rigid.hpp"
#include "load_files/LoadFiles.hpp"

//#include "HipPelvis.hpp"
//#include "HipPelvisCupImplant.hpp"
//#include "HipPelvisCupImplantMatch.hpp"
//#include "HipFemurStemImplant.hpp"
//#include "HipFemurStemImplantMatch.hpp"
#include "vtkCardinalSpline.h"
#include "vtkSplineFilter.h"

#include <qfile.h>
#include <qdir.h>
#include <qstring.h>
#include <qjsondocument.h>
#include <qjsonobject.h>
#include <qjsonvalue.h>
#include <qjsonArray.h>

#include <vector>
#include <string>
#include <sstream>

#include "itkSmoothingRecursiveGaussianImageFilter.h"
#include "test/TestFusionImagen.hpp"

#include "CheckImageRod/CheckImageRod.hpp"

#include <vtkPlane.h>
#include <vtkDistancePolyDataFilter.h>

#include "spine_segmentation/SpineSegmentation.hpp"
#include "spine_segmentation/BallsManualSegmentation.hpp"


std::string dicom = "C:\\DICOM\\dcm2";
std::string wFemoralDir = "C:\\DICOM\\femur";
std::string wTibiaDir = "C:\\DICOM\\tibia";
std::string dicom3d = "C:\\DICOM\\volume\\femur_kevin.nrrd";
std::string dicom_femur = "D:\\3D_DICOM\\femur_seg.nrrd";
std::string dicom_tibia = "D:\\3D_DICOM\\tibia_seg.nrrd";
std::string vtk_kneeCap_right = "D:\\3D_DICOM\\KneeCap_seg_right_surf.vtk";
std::string dicom_knee_cap = "D:\\3D_DICOM\\kneeCap_seg_right.nrrd";
std::string dicom_femur_vtk = "D:\\3D_DICOM\\VTK_BONE\\femur_right_vtk.vtk";
std::string dicom_tibia_vtk = "D:\\3D_DICOM\\VTK_BONE\\tibia_right_vtk.vtk";
std::string dicom_femur_left = "D:\\3D_DICOM\\femur_left.nrrd";
std::string dicom_tibia_left = "D:\\3D_DICOM\\tibia_left.nrrd";
std::string dicom_tibia_left2 = "D:\\3D_DICOM\\tibia_picked.nrrd";
std::string hip3d = "D:\\3D_DICOM\\hip.nrrd";
std::string path = "D:\\Registro\\Points\\Left\\hip_left.txt";
std::string hip = "D:\\workspace\\hip.txt";

std::string dicom_femur_kevin = "D:\\3D_DICOM\\_femur_seg_kevin_last.nrrd";
std::string dicom_tibia_kevin = "D:\\3D_DICOM\\_tibia_seg_kevin_last.nrrd";

std::string femur_path = "D:\\Left_Model\\femurSeg.nrrd";
std::string tibia_path = "D:\\Left_Model\\tibiaSeg.nrrd";
std::string knee_cap = "D:\\Left_Model\\KneeCap.nrrd";

std::string femur_path_poly = "D:\\Left_Model\\femurSegSurf.vtk";
std::string tibia_path_poly = "D:\\Left_Model\\tibiaSegSurf.vtk";

std::string femur_implant_poly = "D:\\Left_Model\\femur_implant.vtk";
std::string knee_cap_poly = "D:\\Left_Model\\KneeCapSurf.vtk";

std::string femur_right_modo_nrrd = "D:\\3D_DICOM\\Virtual_Modo\\Right_Modo\\Femur.nrrd";
std::string tibia_right_modo_nrrd = "D:\\3D_DICOM\\Virtual_Modo\\Right_Modo\\Tibia.nrrd";

using namespace TKA::SEGMENTATION;
using namespace TKA::HIP;
using namespace TKA::IMPLANTS;
using namespace TKA::REGISTRATION;
using namespace IMAGE_ROD;

itk::Point<double, 3> makeItkPoint(double x, double y, double z)
{
	itk::Point<double, 3> point;
	point[0] = x;
	point[1] = y;
	point[2] = z;
	return point;
}


vtkSmartPointer<vtkPolyData> SplitPoly(vtkSmartPointer<vtkPolyData> poly, Plane pPlane, Point ref)
{
	pPlane.reverseByPoint(ref);
	vtkNew<vtkPlane> vtkPlaneMedial;

	Point mPoint = pPlane.getPoint();
	Point normalVector = pPlane.getNormalVector();

	vtkPlaneMedial->SetOrigin(mPoint.x, mPoint.y, mPoint.z);
	vtkPlaneMedial->SetNormal(normalVector.x, normalVector.y, normalVector.z);

	vtkNew<vtkPlaneCollection> medialPlanes;
	medialPlanes->AddItem(vtkPlaneMedial);

	vtkNew<vtkClipClosedSurface> medialClipper;
	medialClipper->SetInputData(poly);
	medialClipper->SetClippingPlanes(medialPlanes);
	medialClipper->Update();

	return medialClipper->GetOutput();
}

void writeListPoints(const std::vector<cv::Point3d>& pPoints, const std::string& pName)
{
	std::ofstream myfile;
	myfile.open(pName);
	for (int i = 0; i < pPoints.size(); i++)
	{
		myfile << "cv::Point3d(" << pPoints[i].x << ", " << pPoints[i].y << ", " << pPoints[i].z << ")," << "\n";

	}
	myfile.close();
}

Point QJsonArrayToPoint(const QJsonArray& pArray)
{
	if (pArray.size() == 3)
	{
		QJsonValue a = pArray[0];
		QJsonValue b = pArray[1];
		QJsonValue c = pArray[2];
		Point result = Point(a.toDouble(), b.toDouble(), c.toDouble());
		//std::cout << result << std::endl;
		return result;
	}
	else
	{
		return Point();
	}
}

Knee CreateKneeFromFile(const std::string& sourcePath)
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

	Point hipCenter = QJsonArrayToPoint(jsonObject.value(QString("hip_center")).toArray());
	Point anteriorCortex = QJsonArrayToPoint(jsonObject.value(QString("anterior_cortex")).toArray());
	Point femurKneeCenter = QJsonArrayToPoint(jsonObject.value(QString("femur_knee_center")).toArray());
	Point lateralEpicondyle = QJsonArrayToPoint(jsonObject.value(QString("lateral_epicondyle")).toArray());
	Point medialEpicondyle = QJsonArrayToPoint(jsonObject.value(QString("medial_epicondyle")).toArray());
	Point lateralCondyle = QJsonArrayToPoint(jsonObject.value(QString("lateral_condyle")).toArray());
	Point medialCondyle = QJsonArrayToPoint(jsonObject.value(QString("medial_condyle")).toArray());
	Point lateralPlateau = QJsonArrayToPoint(jsonObject.value(QString("lateral_plateau")).toArray());
	Point medialPlateau = QJsonArrayToPoint(jsonObject.value(QString("medial_plateau")).toArray());
	Point tibiaKneeCenter = QJsonArrayToPoint(jsonObject.value(QString("tibia_knee_center")).toArray());
	Point tibiaTubercle = QJsonArrayToPoint(jsonObject.value(QString("tibia_tuber")).toArray());
	Point tibiaPCL = QJsonArrayToPoint(jsonObject.value(QString("tibia_pcl")).toArray());
	Point lateralAnkle = QJsonArrayToPoint(jsonObject.value(QString("lateral_ankle")).toArray());
	Point medialAnkle = QJsonArrayToPoint(jsonObject.value(QString("medial_ankle")).toArray());
	Point kneeCap = QJsonArrayToPoint(jsonObject.value(QString("knee_cap")).toArray());

	Point patellaLat = QJsonArrayToPoint(jsonObject.value(QString("patella_lateral")).toArray());
	Point patellaMed = QJsonArrayToPoint(jsonObject.value(QString("patella_medial")).toArray());
	Point patellaInf = QJsonArrayToPoint(jsonObject.value(QString("patella_inf")).toArray());

	QString FemurPolyStr = jsonObject.value(QString("femur_poly")).toString();
	QString TibiaPolyStr = jsonObject.value(QString("tibia_poly")).toString();
	QString PatellaPolyStr = jsonObject.value(QString("patella_poly")).toString();

	if (!(directory.exists(FemurPolyStr) && directory.exists(TibiaPolyStr) && directory.exists(PatellaPolyStr)))
	{
		std::cout << "Error there are not vtk file" << std::endl;
		throw ("Error there are not vtk file");
	}

	QString femurPolyPath = directory.filePath(FemurPolyStr);
	QString tibiaPolyPath = directory.filePath(TibiaPolyStr);
	QString patellaPolyPath = directory.filePath(PatellaPolyStr);

	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly, patellaPoly;

	femurPoly = TestVTK::ReadPolyData(femurPolyPath.toStdString());
	tibiaPoly = TestVTK::ReadPolyData(tibiaPolyPath.toStdString());
	patellaPoly = TestVTK::ReadPolyData(patellaPolyPath.toStdString());

	Point ankleCenter = Knee::getComputeAnkleCenter(lateralAnkle, medialAnkle);

	Patella myPatella;
	myPatella.init(kneeCap, patellaLat, patellaMed, patellaInf, patellaPoly);

	Knee knee;
	knee.init(hipCenter, anteriorCortex, femurKneeCenter, lateralEpicondyle, medialEpicondyle, tibiaKneeCenter,
		tibiaTubercle, tibiaPCL, ankleCenter, myPatella, femurPoly, tibiaPoly, KneeSideEnum::KRight);

	return knee;
}

Knee CreateKneeFromFile_Numbers(const std::string& sourcePath, KneeSideEnum pSide = KneeSideEnum::KRight)
{
	QDir directory(QString::fromStdString(sourcePath));
	QStringList files = directory.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
	if (files.size() == 0)
	{
		std::cout << "Error there are not json file" << std::endl;
		throw ("Error there are not json file");
	}

	QString jsonPath;

	if (files.size() == 1)
	{
		jsonPath = directory.filePath(files[0]);
	}
	else
	{
		QString kneeLandmark = "landmark.json";
		jsonPath = directory.filePath(kneeLandmark);
	}

	QString var;
	QFile file;

	file.setFileName(jsonPath);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	var = file.readAll();
	file.close();
	QJsonDocument jsonDoc = QJsonDocument::fromJson(var.toUtf8());
	QJsonObject jsonObject = jsonDoc.object();

	Point hipCenter = QJsonArrayToPoint(jsonObject.value(QString("0")).toArray());
	Point anteriorCortex = QJsonArrayToPoint(jsonObject.value(QString("4")).toArray());
	Point femurKneeCenter = QJsonArrayToPoint(jsonObject.value(QString("3")).toArray());
	Point lateralEpicondyle = QJsonArrayToPoint(jsonObject.value(QString("2")).toArray());
	Point medialEpicondyle = QJsonArrayToPoint(jsonObject.value(QString("1")).toArray());
	Point lateralCondyle = QJsonArrayToPoint(jsonObject.value(QString("6")).toArray());
	Point medialCondyle = QJsonArrayToPoint(jsonObject.value(QString("5")).toArray());
	Point lateralPlateau = QJsonArrayToPoint(jsonObject.value(QString("8")).toArray());
	Point medialPlateau = QJsonArrayToPoint(jsonObject.value(QString("7")).toArray());
	Point tibiaKneeCenter = QJsonArrayToPoint(jsonObject.value(QString("9")).toArray());
	Point tibiaTubercle = QJsonArrayToPoint(jsonObject.value(QString("10")).toArray());
	Point tibiaPCL = QJsonArrayToPoint(jsonObject.value(QString("11")).toArray());
	Point lateralAnkle = QJsonArrayToPoint(jsonObject.value(QString("13")).toArray());
	Point medialAnkle = QJsonArrayToPoint(jsonObject.value(QString("12")).toArray());
	Point kneeCap = QJsonArrayToPoint(jsonObject.value(QString("17")).toArray());

	Point patellaLat = QJsonArrayToPoint(jsonObject.value(QString("14")).toArray());
	Point patellaMed = QJsonArrayToPoint(jsonObject.value(QString("15")).toArray());
	Point patellaInf = QJsonArrayToPoint(jsonObject.value(QString("16")).toArray());

	QString FemurPolyStr = jsonObject.value(QString("femur_poly")).toString();
	QString TibiaPolyStr = jsonObject.value(QString("tibia_poly")).toString();
	QString PatellaPolyStr = jsonObject.value(QString("patella_poly")).toString();

	if (!(directory.exists(FemurPolyStr) && directory.exists(TibiaPolyStr) && directory.exists(PatellaPolyStr)))
	{
		std::cout << "Error there are not vtk file" << std::endl;
		throw ("Error there are not vtk file");
	}

	QString femurPolyPath = directory.filePath(FemurPolyStr);
	QString tibiaPolyPath = directory.filePath(TibiaPolyStr);
	QString patellaPolyPath = directory.filePath(PatellaPolyStr);

	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly, patellaPoly;

	femurPoly = TestVTK::ReadPolyData(femurPolyPath.toStdString());
	tibiaPoly = TestVTK::ReadPolyData(tibiaPolyPath.toStdString());
	patellaPoly = TestVTK::ReadPolyData(patellaPolyPath.toStdString());

	Point ankleCenter = Knee::getComputeAnkleCenter(lateralAnkle, medialAnkle);

	Patella myPatella;
	//myPatella.init(kneeCap, patellaLat, patellaMed, patellaInf, patellaPoly);

	Knee knee;
	/*knee.init(hipCenter, anteriorCortex, femurKneeCenter, lateralEpicondyle, medialEpicondyle, lateralCondyle, medialCondyle,
		lateralPlateau, medialPlateau, tibiaKneeCenter, tibiaTubercle, tibiaPCL, ankleCenter, myPatella, femurPoly, tibiaPoly);*/

	knee.init(hipCenter, anteriorCortex, femurKneeCenter, lateralEpicondyle, medialEpicondyle, /*lateralPlateau,
		medialPlateau,*/ tibiaKneeCenter, tibiaTubercle, tibiaPCL, ankleCenter, myPatella, femurPoly, tibiaPoly, pSide);

	return knee;
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

	femurPoly = TestVTK::ReadPolyData(femurPolyPath.toStdString());
	tibiaPoly = TestVTK::ReadPolyData(tibiaPolyPath.toStdString());

	UKA::IMPLANTS::Point ankleCenter = Knee::getComputeAnkleCenter(lateralAnkle, medialAnkle);

	//Patella myPatella;
	//myPatella.init(kneeCap, patellaLat, patellaMed, patellaInf, patellaPoly);

	UKA::IMPLANTS::Knee knee;
	
	knee.init(hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondyle, tibiaKneeCenter, 
			  tibiaTubercle, tibiaPCL, ankleCenter, femurPoly, tibiaPoly, pSide, pSurgery);

	return knee;
}

QJsonObject getJSONFromPath(const QString& path, const QDir& directory)
{
	QString jsonPath = directory.filePath(path);
	QString var;
	QFile file;

	file.setFileName(jsonPath);
	file.open(QIODevice::ReadOnly | QIODevice::Text);
	var = file.readAll();
	file.close();

	QJsonDocument jsonDoc = QJsonDocument::fromJson(var.toUtf8());
	QJsonObject jsonObject = jsonDoc.object();

	return jsonObject;
}

double FemurRegistrationFromNumbersTKA(const Knee& myKnee, const std::string& sourcePath, bool& result)
{
	QDir directory(QString::fromStdString(sourcePath));
	QStringList files = directory.entryList(QStringList() << "*.json" << "*.JSON", QDir::Files);
	if (files.size() == 0)
	{
		std::cout << "Error there are not json file" << std::endl;
		throw ("Error there are not json file");
	}
	
	QString femurPoints = "femur_reg_data";
	QString femurLandmark = "femur_reg_landmarks";

	bool var1 = false, var2 = false;

	for (int i = 0; i < files.size(); i++)
	{
		if (files[i].contains(femurPoints) && var1 == false)
		{
			femurPoints = files[i];
			var1 = true;
		}
		
		if (files[i].contains(femurLandmark) && var2 == false)
		{
			femurLandmark = files[i];
			var2 = true;
		}

	}

	if (var1 == false || var2 == false)
	{
		std::cout << "Error there are not femur json files" << std::endl;
		throw ("Error there are not femur json files");
	}

	QJsonObject jsonObjectLandmark = getJSONFromPath(femurLandmark, directory);
	QJsonObject jsonObjectPoints = getJSONFromPath(femurPoints, directory);

	auto femurLandmarkPoints = jsonObjectLandmark.value(QString("points")).toObject();
	auto femurRegistrationPoints = jsonObjectPoints.value(QString("points")).toObject();
	
	Point hipCenter = QJsonArrayToPoint(femurLandmarkPoints.value(QString("0")).toArray());
	Point kneeCenter = QJsonArrayToPoint(femurLandmarkPoints.value(QString("1")).toArray());
	Point medialEpicondyle = QJsonArrayToPoint(femurLandmarkPoints.value(QString("3")).toArray());

	std::vector<itk::Point<double, 3>> pBonePoints;

	for (int i = 0; i < 7; i++)
	{
		auto obj = femurRegistrationPoints.value(QString::number(i)).toObject();

		for (auto it = obj.begin(); it != obj.end(); ++it)
		{
			Point temp = QJsonArrayToPoint(it->toArray());
			itk::Point<double, 3> pointITK;
			pointITK[0] = temp.x;
			pointITK[1] = temp.y;
			pointITK[2] = temp.z;

			pBonePoints.push_back(pointITK);
		}
	}

	FemurRegistration registration(myKnee.GetFemurPoly(), myKnee.getHipCenter().ToITKPoint(),
									myKnee.getFemurKneeCenter().ToITKPoint(), 
									myKnee.getMedialEpicondyle().ToITKPoint());

	result = registration.MakeRegistration(pBonePoints, Registration::makeItkPoint(hipCenter.x, hipCenter.y, hipCenter.z),
									Registration::makeItkPoint(kneeCenter.x, kneeCenter.y, kneeCenter.z),
									Registration::makeItkPoint(medialEpicondyle.x, medialEpicondyle.y, medialEpicondyle.z), true);
	
	return registration.getError();
}

FemurImplant CreateFemurImplantFromFile(const std::string& sourcePath)
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

	QJsonArray planeA = jsonObject.value(QString("face_a")).toArray();
	QJsonArray planeB = jsonObject.value(QString("face_b")).toArray();
	QJsonArray planeC = jsonObject.value(QString("face_c")).toArray();
	QJsonArray planeD = jsonObject.value(QString("face_d")).toArray();
	QJsonArray planeE = jsonObject.value(QString("face_e")).toArray();

	Plane A, B, C, D, E;

	A.init(Point(planeA[0].toDouble(), planeA[1].toDouble(), planeA[2].toDouble()), Point(planeA[3].toDouble(), planeA[4].toDouble(), planeA[5].toDouble()));
	B.init(Point(planeB[0].toDouble(), planeB[1].toDouble(), planeB[2].toDouble()), Point(planeB[3].toDouble(), planeB[4].toDouble(), planeB[5].toDouble()));
	C.init(Point(planeC[0].toDouble(), planeC[1].toDouble(), planeC[2].toDouble()), Point(planeC[3].toDouble(), planeC[4].toDouble(), planeC[5].toDouble()));
	D.init(Point(planeD[0].toDouble(), planeD[1].toDouble(), planeD[2].toDouble()), Point(planeD[3].toDouble(), planeD[4].toDouble(), planeD[5].toDouble()));
	E.init(Point(planeE[0].toDouble(), planeE[1].toDouble(), planeE[2].toDouble()), Point(planeE[3].toDouble(), planeE[4].toDouble(), planeE[5].toDouble()));

	Point P1 = QJsonArrayToPoint(jsonObject.value(QString("p1")).toArray());
	Point P2 = QJsonArrayToPoint(jsonObject.value(QString("p2")).toArray());
	Point P3 = QJsonArrayToPoint(jsonObject.value(QString("p3")).toArray());

	FemurImplantInfo implantFemur;
	implantFemur.femurPosteriorLateralThickness = jsonObject.value(QString("posterior_lateral_thickness")).toDouble();
	implantFemur.femurDistalLateralThickness = jsonObject.value(QString("distal_lateral_thickness")).toDouble();
	implantFemur.femurPosteriorMedialThickness = jsonObject.value(QString("posterior_medial_thickness")).toDouble();
	implantFemur.femurDistalMedialThickness = jsonObject.value(QString("distal_medial_thickness")).toDouble();

	std::vector<PointTypeITK> kneeCapPath;

	QString FemurPolyStr = jsonObject.value(QString("filename")).toString();

	if (!(directory.exists(FemurPolyStr)))
	{
		std::cout << "Error there are not vtk file" << std::endl;
		throw ("Error there are not vtk file");
	}

	QString femurPolyPath = directory.filePath(FemurPolyStr);

	vtkSmartPointer<vtkPolyData> femurPoly;

	femurPoly = TestVTK::ReadPolyDataSTL(femurPolyPath.toStdString());


	FemurImplant femurImplant;
	femurImplant.init(A, B, C, D, E, P1, P2, P3, femurPoly, kneeCapPath, implantFemur);

	return femurImplant;
}

TibiaImplant CreateTibiaImplantFromFile(const std::string& sourcePath, vtkSmartPointer<vtkPolyData>& pTibiaPoly)
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

	Point pcl1 = QJsonArrayToPoint(jsonObject.value(QString("point_pcl1")).toArray());
	Point pcl2 = QJsonArrayToPoint(jsonObject.value(QString("point_pcl2")).toArray());
	Point front = QJsonArrayToPoint(jsonObject.value(QString("point_front")).toArray());
	Point exterior = QJsonArrayToPoint(jsonObject.value(QString("point_external")).toArray());

	TibiaImplantInfo implantTibia;
	implantTibia.tibiaLateralThickness = jsonObject.value(QString("thickness_left")).toDouble();
	implantTibia.tibiaMedialThickness = jsonObject.value(QString("thickness_right")).toDouble();


	QString TibiaPolyStr = jsonObject.value(QString("filename")).toString();

	if (!(directory.exists(TibiaPolyStr)))
	{
		std::cout << "Error there are not vtk file" << std::endl;
		throw ("Error there are not vtk file");
	}

	QString tibiaPolyPath = directory.filePath(TibiaPolyStr);

	pTibiaPoly = TestVTK::ReadPolyDataSTL(tibiaPolyPath.toStdString());

	TibiaImplant tibiaImplant;
	tibiaImplant.init(pcl1, pcl2, front, exterior, implantTibia);

	return tibiaImplant;
}

PatellaImplant CreatePatellaImplantFromFile(const std::string& sourcePath, vtkSmartPointer<vtkPolyData>& pPatellaPoly)
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

	Point point1 = QJsonArrayToPoint(jsonObject.value(QString("base_point1")).toArray());
	Point point2 = QJsonArrayToPoint(jsonObject.value(QString("base_point2")).toArray());
	Point point3 = QJsonArrayToPoint(jsonObject.value(QString("base_point3")).toArray());
	Point exterior = QJsonArrayToPoint(jsonObject.value(QString("top_point")).toArray());

	QString PolyStr = jsonObject.value(QString("filename")).toString();

	if (!(directory.exists(PolyStr)))
	{
		std::cout << "Error there are not vtk file" << std::endl;
		throw ("Error there are not vtk file");
	}

	QString PolyPath = directory.filePath(PolyStr);

	pPatellaPoly = TestVTK::ReadPolyDataSTL(PolyPath.toStdString());

	PatellaImplant Implant;
	Implant.init(point1, point2, point3, exterior);

	return Implant;
}

std::vector<cv::Point3d> GetSpineLineFromFile(const std::string& sourcePath)
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

	auto pointList = jsonObject.value(QString("centerPhysicalPoints")).toArray();
	std::vector<cv::Point3d> result;

	for (auto& item : pointList)
	{
		Point point1 = QJsonArrayToPoint(item.toArray());
		result.push_back(point1);
	}
	return result;
}

std::vector<Point> ITKVectorToCV(const std::vector<PointTypeITK>& pData)
{
	std::vector<Point> result;
	auto it1 = pData.begin();
	auto it2 = pData.end();

	for (; it1 != it2; ++it1)
	{
		Point temp = { (*it1)[0], (*it1)[1], (*it1)[2] };
		result.push_back(temp);
	}
	return result;
}

void EraseImagePartUsingPlane(ImageType::Pointer pImg, Plane& pPlane)
{
	ImageType::PixelType* fullBuffer = pImg->GetBufferPointer();
	ImageType::SizeType dims = pImg->GetLargestPossibleRegion().GetSize();

	for (int k = 0; k < dims[2]; k++)
	{
		for (int j = 0; j < dims[1]; j++)
		{
			for (int i = 0; i < dims[0]; i++)
			{
				int index = k * dims[1] * dims[0] + j * dims[0] + i;
				if (pPlane.eval(Point(i, j, k)) < 0)
				{
					fullBuffer[index] = 0;
				}
			}
		}
	}
}

std::vector<cv::Point3d> PointToPoint3D(const std::vector<Point>& data)
{
	std::vector<cv::Point3d> result;
	for (int i = 0; i < data.size(); i++)
	{
		result.push_back(data[i]);
	}
	return result;
}

void executeHip()
{
	//Help for hip center lib

	/*
	Constructor Parameters :
		std::vector<HipPoint> pPointList : A vector of the structure (HipPoint). It is a vectors of points.

	Methods for use:
		getHipCenterBySphere(double pointOut[3]) : Return Hip center point as array[3] parameter.
		getHipCenterByThreeSphere(double pointOut[3]) : Return Hip center point as array[3] parameter.
	*/

	//Example reading data from a file for testing.

	double hipCenter_1[3], hipCenter_2[3];
	//std::string path = "bone.txt";

	std::ifstream infile(path);
	std::vector<HipPoint> points;
	double a, b, c;
	HipPoint myPoint;
	while (infile >> a >> b >> c)
	{
		myPoint.x = a;
		myPoint.y = b;
		myPoint.z = c;
		points.push_back(myPoint);
	}

	HipCenter hip(points);
	hip.GetHipCenterBySphere(hipCenter_1);
	hip.GetHipCenterByThreeSphere(hipCenter_2);
	std::cout << "Hip Center By Method 1: " << hipCenter_1[0] << " " << hipCenter_1[1] << " " << hipCenter_1[2] << std::endl;
	std::cout << "Hip Center By Method 2: " << hipCenter_2[0] << " " << hipCenter_2[1] << " " << hipCenter_2[2] << std::endl;
}

void executeSegmentation(const std::string& kneePath)
{
	//Help for segmentation lib

	/*
	Constructor Parameters :
		const std::string& dirName : Path of the folder with CT images.

	Methods for use:
		- ExecuteSegmentation(const char* wFemoralDir, const char* wTibiaDir, const std::string& seriesName = "") :
		This method execute the segmentation and receives as parameters the path of the folders to save the segmentation of the femur and tibia.

		- vis3d(SegmentImageType* image) : This visualize a segmented image
	*/

	// Example of use

	//std::string right_knee = "C:\\DICOM\\dcm2";
	/*std::string dicomLeft = "C:\\DICOM\\left_leg\\knee";
	std::string dicomLeftAnkle = "C:\\DICOM\\left_leg\\ankle";
	std::string dicomLeftHip = "C:\\DICOM\\left_leg\\hip";
	std::string pelvis = "C:\\DICOM\\Pelvis\\bone";*/

	//ImageType::Pointer inputImg = Test::ReadSeries<ImageType>(right_knee);
	//Test::SaveImage<ImageType>(inputImg, "original_right_knee");

	ImageType::Pointer inputImg;

	Test::readImage<ImageType>(kneePath, inputImg);

	ManualSegmentation segmentation(inputImg);

	auto start = std::chrono::system_clock::now();

	short min = 300;
	short max = 3036;

	ImageType::Pointer resultSeg = segmentation.ApplyThresholds(min, max);

	itk::Vector<double, 3> target;
	target[0] = -1;
	target[1] = 0;
	target[2] = 1;

	itk::Matrix< double, 3, 3 > pRotationOut;

	ImageType::Pointer resultImg = segmentation.MakeRotation(inputImg, ManualSegmentation::kSagital_X, target, pRotationOut);
	resultSeg = segmentation.MakeRotation(resultSeg, ManualSegmentation::kSagital_X, target, pRotationOut, 0);

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	//auto restore = segmentation.FillAndRestoreImageSegmentation(resultSeg, pRotationOut);

	Test::SaveImage<ImageType>(resultImg, "Rotate_Image");
	Test::SaveImage<ImageType>(resultSeg, "Rotate_Image_Seg");
	//Test::SaveImage<SegmentImageType>(restore, "Rotate_Image_Seg2");

	/*
	auto start = std::chrono::system_clock::now();
	// Some computation here
	AutomaticSegmentation segmentation(inputImg);
	segmentation.ExecuteSegmentation();

	SegmentImageType::Pointer femur = segmentation.GetFemoralSegment();
	SegmentImageType::Pointer tibia = segmentation.GetTibiaSegment();
	SegmentImageType::Pointer kneeCap = segmentation.GetKneeCapSegment();
	SegmentImageType::Pointer fibula = segmentation.GetFibulaSegment();

	auto end = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = end - start;

	std::cout << "elapsed time: " << elapsed_seconds.count() << "s\n";

	Test::SaveImage<SegmentImageType>(tibia, "segment_tibia_auto");
	Test::SaveImage<SegmentImageType>(femur, "segment_femur_auto");
	Test::SaveImage<SegmentImageType>(kneeCap, "segment_kneecap_auto");
	Test::SaveImage<SegmentImageType>(fibula, "segment_fibula_auto");
	*/
}

void executeImplantsMatch()
{
	//Help for implants lib

	/*
		In order to work with this library you must create the objects FemurImplantMatch and TibiaImplantMatch.
		These are the steps:
		1-) You should have the following points get from the CT image (These points are getting  by the doctor during the pre-planning operation.).
			hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondyle, medialCondyle,
			lateralCondyle, anteriorCortex, ankleCenter, lateralPlateau, medialPlateau,
			tibiaKneeCenter, tibiaTubercle, pclCenter.
			With each of these points you must create an object of class Point(x, y, z);

		2-) You should have a pointer to a 3d volume of the femur and another of the tibia.

		3-) You should have 3 points for each face of the femur knee implant (A, B, C, D, E).

		4-) You must have 7 specific points of the femur knee implant (P1, P2, P3, P4, P5, P6, P7)

		5-) You should have 3 points for the tibia knee implants face (tibiaP1, tibiaP2, tibiaP3).

		6-) You must have 1 specific point of the tibia knee implant (tibiaExternal)

		** With all this data you can start creating the objects how is explained.

		0-) Create ImplantInfo structure object. This structure stores the information
			on the dimensions of the implant.

		1-) Using the plane class you must create a Plane object for eah femur implant face
			using the face points.Exmple: Plane.init (PointFaceA1, PointFaceA2, PointFaceA3).

		2-) Now you must create the objects (FemurImplant, TibiaImplant, Knee)
			- FemurImplant : Receive the 5 planes of the implant and the its 7 points.
				femurImplant.init(A, B, C, D, E, P1, P2, P3, P4, P5, P6, P7);

			- TibiaImplant : Receive the 4 points of the tibia implant.
				tibiaImplant.init(tibiaP1, tibiaP2, tibiaP3, tibiaExternal);

			- Knee : Receive the points obtained from the CT images,
					the characteristics of the implant in the ImplantInfo structure and
					the pointers to the 3d volumes of the femur and tibia, as shown in the example below.

		3-) Finally, you can create the objects FemurImplantMatch and TibiaImplantMatch
			- FemurImplantMatch : Receive FemurImplant objects and Knee objects.
			- TibiaImplantMatch : Receive TibiaImplant objects and Knee objects.

		Note: The objective of these objects is to provide the equation of the planes
			of each of the faces of the implants when they adhere to the bone. With these planes
			you can know the position of the implant in the bone and color the area of the bone
			where the implant is placed.

			The implant can be positioned in the bone with reference to the condyles or the cortex.
			If you want to get the position according to the condyles you must pass the getPlane
			method parameter in true. If you use false the positioning will be referenced to the cortex.

	*/

	/*ImplantImageType::Pointer femur, tibia;
	Test::readImage<ImplantImageType>(femur_path, femur);
	Test::readImage<ImplantImageType>(tibia_path, tibia);*/

	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly, femurImplantPoly, patellaPoly;

	femurPoly = TestVTK::ReadPolyData("D:\\Excel\\femurSegSurf.vtk");//TestVTK::ReadPolyData("D:\\3D_DICOM\\Real_Person\\segmentation\\femur_poly.vtk"); //TestVTK::ReadPolyData(femur_path_poly);
	tibiaPoly = TestVTK::ReadPolyData("D:\\Excel\\tibiaSegSurf.vtk");//TestVTK::ReadPolyData("D:\\3D_DICOM\\Real_Person\\segmentation\\tibia_poly.vtk"); //TestVTK::ReadPolyData(tibia_path_poly);

	femurImplantPoly = TestVTK::ReadPolyData("D:\\Excel\\femur_implant.vtk"); //TestVTK::ReadPolyData(femur_implant_poly);

	FemurImplantInfo implantFemur;
	implantFemur.femurPosteriorLateralThickness = 9;
	implantFemur.femurDistalLateralThickness = 9;
	implantFemur.femurPosteriorMedialThickness = 9;
	implantFemur.femurDistalMedialThickness = 9;

	TibiaImplantInfo implantTibia;
	implantTibia.tibiaLateralThickness = 11;
	implantTibia.tibiaMedialThickness = 11;

	/* left model
	Point hipCenter(28.8433, -57.3193, 498.549);
	Point femurKneeCenter(-13.5082, -7.24707, 92.3906);
	Point lateralEpicondyle(-45.6781, -17.499, 106.142);
	Point medialEpicondyle(28.9137, -13.2568, 106.36);
	Point anteriorCortex(-20.225, 13.6104, 138.595);
	Point medialCondyle(13.359, -39.0635, 98.9912);
	Point lateralCondyle(-31.891, -38.3564, 101.191);

	Point lateralPlateau(14.7731, -11.4893, 74.239);
	Point medialPlateau(-18.4574, -6.54004, 75.3391);
	Point tibiaPCL(-6.79137, -29.5186, 67.0884);
	Point tibiaKneeCenter(-6.79137, 2.29785, 73.689);
	Point tibiaTubercle(1.693, 15.7314, 50.0369);
	//Point lateralAnkle(-19.30, -74.89, 384.90);
	//Point medialAnkle(40.89, -83.05, 396.90);
	Point ankleCenter(2.21859, -7.30133, -291.76);*/

	/* Real bone
	Point hipCenter(55.13, -118.85, -283.40);
	Point femurKneeCenter(82.72, -101.77, -648.11);
	Point lateralEpicondyle(113.24, -85.67, -633.41);
	Point medialEpicondyle(42.88, -93.72, -632.01);
	Point anteriorCortex(93.32, -121.27, -599.81);
	Point medialCondyle(58.56, -76.77, -638.31);
	Point lateralCondyle(97.56, -69.56, -636.21);

	Point lateralPlateau(101.79, -89.48, -658.61);
	Point medialPlateau(67.89, -97.11, -660.71);
	Point tibiaPCL(78.06, -75.50, -667.01);
	Point tibiaKneeCenter(91.20, -103.47, -658.61);
	Point tibiaTubercle(95.01, -121.27, -690.11);
	Point lateralAnkle(105.71, -70.22, -969.80);
	Point medialAnkle(55.42, -84.32, -961.40);*/


	Point hipCenter(52.5264, -122.331, -282.699);
	Point femurKneeCenter(82.2979, -105.165, -648.108);
	Point lateralEpicondyle(113.237, -85.6689, -635.509);
	Point medialEpicondyle(42.2708, -98.3838, -630.02);
	Point anteriorCortex(93.7412, -121.271, -591.412);
	Point medialCondyle(55.5967, -76.3447, -639.009);
	Point lateralCondyle(98.8271, -69.1396, -636.909);

	Point lateralPlateau(100.946, -90.7549, -658.607);
	Point medialPlateau(66.1924, -97.5361, -660.707);
	Point tibiaPCL(77.2119, -80.1592, -659.307);
	Point tibiaKneeCenter(82.7217, -111.522, -657.207);
	Point tibiaTubercle(94.165, -121.694, -694.305);
	Point ankleCenter(77.622, -79.5502, -969.079);

	Point planeA1(-24.6383, 55.6821, 34.4108);
	Point planeA2(0, 55.8612, 24.148);
	Point planeA3(24.6383, 55.6821, 34.4108);

	Point planeB1(-31.33, 56.294, 15.783);
	Point planeB2(-31.33, 63.8159, 8.39129);
	Point planeB3(31.33, 63.8159, 8.39129);

	Point planeC1(-31.2913, 64.5168, 8.10454);
	Point planeC2(-29.421, 83.7645, 8.10454);
	Point planeC3(28.2063, 83.7645, 8.10454);

	Point planeD1(-29.3421, 84.4655, 8.39129);
	Point planeD2(-27.0605, 98.6397, 22.3203);
	Point planeD3(22.6364, 98.6397, 22.3203);

	Point planeE1(-26.967, 98.935, 22.9464);
	Point planeE2(22.404, 98.935, 22.9464);
	Point planeE3(-14.8076, 101.72, 54.7848);

	Point P1(-31.33, 56.0022, 16.0697);
	Point P2(31.33, 56.0022, 16.0697);
	Point P3(-31.33, 64.1077, 8.10454);
	Point P4(31.33, 64.1077, 8.10454);
	Point P5(-29.3813, 84.1737, 8.10454);
	Point P6(28.1412, 84.1737, 8.10454);
	Point P7(-14.8076, 101.72, 54.7848);

	Point implantPCLP1(9.45, -15.18, -3.02);
	Point implantPCLP2(-9.45, -15.18, -3.02);
	Point implantFront(0, 28.35, 3.48);

	Point firstPlateau(-22.31, -1.25, 7.0);
	Point secondPlateau(22.37, -1.24, 7.0);

	Point kneeCap(0, 0, 0);

	Plane A, B, C, D, E;
	A.init(planeA1, planeA2, planeA3);
	B.init(planeB1, planeB2, planeB3);
	C.init(planeC1, planeC2, planeC3);
	D.init(planeD1, planeD2, planeD3);
	E.init(planeE1, planeE2, planeE3);

	std::cout << "**************Planes******************" << std::endl;

	auto startTime = std::chrono::system_clock::now();

	std::vector<PointTypeITK> kneeCapPath;
	FemurImplant femurImplant;
	TibiaImplant tibiaImplant;
	femurImplant.init(A, B, C, D, E, P4, P3, P7, femurImplantPoly, kneeCapPath, implantFemur);
	tibiaImplant.init(implantPCLP1, implantPCLP2, implantFront, firstPlateau, implantTibia);

	std::cout << "**************Transformation Planes******************" << std::endl;

	Patella myPatella;
	myPatella.init(kneeCap, lateralEpicondyle, medialEpicondyle, femurKneeCenter, patellaPoly); // Just for avoid problem with patella

	Knee knee;
	knee.init(hipCenter, anteriorCortex, femurKneeCenter, lateralEpicondyle, medialEpicondyle, tibiaKneeCenter, 
		tibiaTubercle, tibiaPCL, ankleCenter, myPatella, femurPoly, tibiaPoly, KneeSideEnum::KRight);
	FemurImplantMatch femurImplantMatch;
	TibiaImplantMatch tibiaImplantMatch;

	try
	{
		femurImplantMatch.init(femurImplant, knee);
		tibiaImplantMatch.init(tibiaImplant, knee);
	}
	catch (const ImplantsException& e)
	{
		std::cout << e.what() << std::endl;
	}

	femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneA).show();
	femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneB).show();
	femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneC).show();
	femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneD).show();
	femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneE).show();
	femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneMid).show();
	tibiaImplantMatch.getTibiaPlane().show();

	auto endTime = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = endTime - startTime;

	std::cout << "****elapsed time genaral implant match: " << elapsed_seconds.count() << std::endl;

	////////////////////////////Getting Hull Points//////////////////

	itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer femurTransformOut = itk::VersorRigid3DTransform<double>::New();

	itk::Rigid3DTransform<double>::Pointer tibiaTransformIn = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer tibiaTransformOut = itk::VersorRigid3DTransform<double>::New();

	femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
	femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrix());

	try
	{
		std::vector<PointTypeITK> hull = femurImplantMatch.GetHullPoints(femurTransformIn, femurTransformOut, FemurImplantMatch::kPlaneA, 0);
		std::cout << "femurTransformOut: " << femurTransformOut << std::endl;
		femurImplantMatch.GetHullPoints(femurTransformIn, femurTransformOut, FemurImplantMatch::kPlaneB, 0);
		femurImplantMatch.GetHullPoints(femurTransformIn, femurTransformOut, FemurImplantMatch::kPlaneC, 0);
		femurImplantMatch.GetHullPoints(femurTransformIn, femurTransformOut, FemurImplantMatch::kPlaneD, 0);
		femurImplantMatch.GetHullPoints(femurTransformIn, femurTransformOut, FemurImplantMatch::kPlaneE, 0);
		std::cout << "Hull femur size: " << hull.size() << std::endl;
	}
	catch (const ImplantsException& e)
	{
		Test::myPrint("Error femur get hull");
		std::cout << e.what() << std::endl;
	}

	tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
	tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());
	try
	{
		std::vector<PointTypeITK> hull = tibiaImplantMatch.GetHullPoints(tibiaTransformIn, tibiaTransformOut, 0, 0);
		std::cout << "Hull tibia size: " << hull.size() << std::endl;
	}
	catch (const ImplantsException& e)
	{
		Test::myPrint("Error tibia get hull");
		std::cout << e.what() << std::endl;
	}

	std::cout << "*************************-Femur-*****************************************" << std::endl;

	TestVTK::TransformPoly(femurImplantPoly, femurImplantMatch.GetRotationMatrix(), femurImplantMatch.GetTranslationMatrix());

	Plane test = femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneC);
	test.normalizeNormalVector();

	Point vectorAxis = hipCenter - femurKneeCenter;
	vectorAxis = vectorAxis / sqrt(vectorAxis.dot(vectorAxis));

	Point midTemp = (P3 + P4) / 2.0;
	cv::Mat midTempMat(3, 1, CV_64F);
	midTempMat.at<double>(0, 0) = midTemp.x;
	midTempMat.at<double>(1, 0) = midTemp.y;
	midTempMat.at<double>(2, 0) = midTemp.z;

	cv::Mat RotationMat = Test::Rigid3DTransformToCVRotation(femurTransformIn);
	cv::Mat translation = Test::Rigid3DTransformToCVTranslation(femurTransformIn);

	cv::Mat midMat = RotationMat * midTempMat + translation;
	Point mid = Point(midMat);
	double evalMid = femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneMid).eval(mid);
	double evalKneeCenter = femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneMid).eval(femurKneeCenter);

	double cortexEval = femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneE).eval(knee.getAnteriorCortex());

	double cortexKneeCenterOnC = femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneC).eval(femurKneeCenter);

	std::cout << "********************* Implant: " << test.getNormalVector() << std::endl;
	std::cout << "************************ Axis: " << vectorAxis << std::endl;
	std::cout << "********************mid point: " << mid << std::endl;
	std::cout << "***************eval mid point: " << evalMid << std::endl;
	std::cout << "*************eval knee center: " << evalKneeCenter << std::endl;
	std::cout << "*************eval cortex: " << cortexEval << std::endl;
	std::cout << "*************eval knee center on C: " << cortexKneeCenterOnC << std::endl;

	std::cout << "*************************-Tibia-*****************************************" << std::endl;

	Point tibiaAxis = tibiaKneeCenter - ankleCenter;
	double angle = Line::getAngleBetweenVectors(tibiaAxis, tibiaImplantMatch.getTibiaPlane().getNormalVector());

	Point centroidTibiaImplant = (implantPCLP1 + implantPCLP2 + implantFront) / 3.0;
	cv::Mat centroidTibiaMat(3, 1, CV_64F);
	centroidTibiaMat.at<double>(0, 0) = centroidTibiaImplant.x;
	centroidTibiaMat.at<double>(1, 0) = centroidTibiaImplant.y;
	centroidTibiaMat.at<double>(2, 0) = centroidTibiaImplant.z;

	cv::Mat rotationTibia = Test::Rigid3DTransformToCVRotation(tibiaTransformIn);
	cv::Mat translationTibia = Test::Rigid3DTransformToCVTranslation(tibiaTransformIn);

	cv::Mat centroidTibiaTransformMat = rotationTibia * centroidTibiaMat + translationTibia;


	Point centroidTibia = Point(centroidTibiaTransformMat);

	double evalTibia = tibiaImplantMatch.getTibiaPlane().eval(centroidTibia);

	std::cout << "******************** angle tibia: " << angle * 180 / 3.14 << std::endl;
	std::cout << "**********eval implant mid point: " << evalTibia << std::endl;
	//std::cout << "**********eval tibia knee center: " << tibiaImplantMatch.getTibiaPlane().eval(newTibiaCenter) << std::endl;*/

}

void MatchEasy()
{
	/*Knee kneeLeftModo = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Modo\\Right_Modo");
	Knee kneeLeft = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Right");

	Knee myKnee = kneeLeftModo;

	FemurImplantMatch femurImplantMatch;*/

	std::string femurImplantStr = "D:\\sovajo\\Errores\\Error6\\Femur_Implant";
	std::string tibiaImplantStr = "D:\\sovajo\\Errores\\Error6\\Tibia_Implant";
	//std::string patellaImplantStr = "D:\\3D_DICOM_DATA\\patella_Implant";

	/*Knee kneeLeftModo = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Modo\\Right_Modo");
	Knee kneeLeft = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Right");*/

	Knee myKnee = CreateKneeFromFile_Numbers("D:\\sovajo\\Errores\\Error6", KneeSideEnum::KRight);

	vtkSmartPointer<vtkPolyData> polyTibiaImplant, polyPatellaImplant;

	FemurImplant femurImplant = CreateFemurImplantFromFile(femurImplantStr);

	TibiaImplant tibiaImplant = CreateTibiaImplantFromFile(tibiaImplantStr, polyTibiaImplant);

	//PatellaImplant patellaImplant = CreatePatellaImplantFromFile(patellaImplantStr, polyPatellaImplant);

	FemurImplantMatch femurImplantMatch;
	TibiaImplantMatch tibiaImplantMatch;
	//PatellaImplantMatch patellaImplantMatch;

	femurImplantMatch.init(femurImplant, myKnee);

	tibiaImplantMatch.init(tibiaImplant, myKnee);


	//patellaImplantMatch.init(patellaImplant, myKnee);
	/*vtkSmartPointer<vtkPolyData> newImplant = TestVTK::TransformPoly(polyPatellaImplant, patellaImplantMatch.GetRotationMatrix(), patellaImplantMatch.GetTranslationMatrix());
	std::vector<vtkSmartPointer<vtkPolyData>> polyList;
	polyList.push_back(newImplant);
	TestVTK::show(patellaImplantMatch.GetCuttingPatella(), polyList);*/

	std::vector<double> my_Transform = { -0.8510644835740226, 0.4940021459999185, 0.17790762925282486, 107.34124402717282,
										  0.49079754865801617, 0.8688639863636932, -0.06475445491607344, 129.99908421613884,
										-0.1865663716487183, 0.032206411592580066, -0.981914322139241, -993.3514342643502 };

	itk::Rigid3DTransform<double>::Pointer transformIn = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer transformTibia = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer transformInTemp = Test::MakeTransformITK(my_Transform);

	transformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
	transformIn->SetOffset(femurImplantMatch.GetTranslationMatrixByCortex());

	transformTibia->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
	transformTibia->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

	/*itk::Rigid3DTransform<double>::Pointer composedTransform = itk::VersorRigid3DTransform<double>::New();
	composedTransform->Compose(transformIn);
	composedTransform->Compose(transformInTemp);*/

	std::vector<PointTypeITK> hull1;

	itk::Rigid3DTransform<double>::Pointer transformOut = itk::VersorRigid3DTransform<double>::New();

	std::vector<cv::Point3d> tPoints, tPoints2;

	try
	{

		//hull1 = femurImplantMatch.GetHullPoints(transformIn, transformOut, FemurImplantMatch::kPlaneB, 1, 2, 10, 15, 200);
		hull1 = tibiaImplantMatch.GetHullPoints(transformTibia, transformOut, 1, 10);
		std::cout << transformIn << std::endl;

		for (int i = 0; i < hull1.size(); i++)
		{
			Point myPoint(hull1[i][0], hull1[i][1], hull1[i][2]);

			tPoints.push_back(myPoint);
		}

		std::cout << " ********* Hull size: " << hull1.size() << std::endl;
	}
	catch (const ImplantsException& e)
	{
		Test::myPrint("Error get hull");
		std::cout << e.what() << std::endl;
	}
	//std::cout << std::endl;
	//std::cout << std::endl;

	//ImplantsMatchFinalInfo matchInfo(&myKnee, femurImplant, tibiaImplant, transformIn, transformTibia);

	/*std::cout << "Femur varus angle: " << matchInfo.GetFemurVarusAngle() << std::endl;
	std::cout << "Femur implant PCA angle: " << matchInfo.GetFemurImplantPCAAngle() << std::endl;
	std::cout << "Femur implant TEA angle " << matchInfo.GetFemurImplantTEAAngle() << std::endl;
	std::cout << "Femur flexion angle: " << matchInfo.GetFemurFlexionAngle() << std::endl;

	std::cout << std::endl;*/
	/*
	matchInfo.setFemurVarusAngle(30);
	matchInfo.setFemurPCAAngle(30);
	matchInfo.setFemurTEAAngle(30);
	matchInfo.setFemurFlexionAngle(45);

	matchInfo.SetThicknessFemurAxialLateral(2);
	matchInfo.SetThicknessFemurCoronalLateral(3);

	std::cout << "Femur varus angle: " << matchInfo.GetFemurVarusAngle() << std::endl;
	std::cout << "Femur implant PCA angle: " << matchInfo.GetFemurImplantPCAAngle() << std::endl;
	std::cout << "Femur implant TEA angle: " << matchInfo.GetFemurImplantTEAAngle() << std::endl;
	std::cout << "Femur flexion angle: " << matchInfo.GetFemurFlexionAngle() << std::endl;

	std::cout << std::endl;
	std::cout << "***************** After change *********************" << std::endl;
	std::cout << std::endl;

	ImplantsMatchFinalInfo matchInfo1(&myKnee, femurImplant, tibiaImplant, matchInfo.getITKFemurTransform(), transformTibia);

	std::cout << "Femur varus angle: " << matchInfo1.GetFemurVarusAngle() << std::endl;
	std::cout << "Femur implant PCA angle: " << matchInfo1.GetFemurImplantPCAAngle() << std::endl;
	std::cout << "Femur implant TEA angle: " << matchInfo1.GetFemurImplantTEAAngle() << std::endl;
	std::cout << "Femur flexion angle: " << matchInfo.GetFemurFlexionAngle() << std::endl;

	std::cout << std::endl;
	*/


	//std::cout << "SI Angle: " << matchInfo.getAngleSI() << std::endl;
	//std::cout << "Cartilage before : " << myKnee.getFemurCartilage() << std::endl;

	//vtkSmartPointer<vtkPolyData> polyNew1 = TestVTK::CreatePolyLine(hull1);
	//TestVTK::show(myKnee.GetFemurPoly());
	TestVTK::show(myKnee.GetTibiaPoly(), tPoints, true);

	//std::vector<cv::Point3d> testRef = { myKnee.getMedialInferiorFemurPoint(), myKnee.getLateralInferiorFemurPoint(), myKnee.getMedialCondyle(), myKnee.getLateralCondyle() };

	/*tPoints2.push_back(tPoints[0]);
	tPoints2.push_back(tPoints[tPoints.size() - 1]);*/

	//TestVTK::show(myKnee.GetFemurPoly(), testRef);
}

void MatchTibiaPKA()
{
	std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data";
	/*
		Here I am creating a right knee from a JSON file. The surgical side is medial.
	*/
	UKA::IMPLANTS::Knee myKnee = CreateKneeFromFile_NumbersPKA(folder, UKA::IMPLANTS::KRight, UKA::IMPLANTS::KMedial);
	UKA::IMPLANTS::TibiaImplant tibiaImplant;

	auto tibiaModel = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\tibia.STL");
	UKA::IMPLANTS::Point apLinePclPoint(10.4801, 9.06646, 0.934849);
	UKA::IMPLANTS::Point apLineTuberPoint(10.2438, -37.5031, 0.951554);
	UKA::IMPLANTS::Point sidePoint(-14.1021, -11.1526, 1.1545);
	UKA::IMPLANTS::Point exteriorPoint(-14.5152, -10.8294, 4.41175);
	UKA::IMPLANTS::Point planeSidePoint(12.10, -12.08, 2.59);
	UKA::IMPLANTS::TibiaImplantInfo tibiaInfo;
	tibiaInfo.tibiaThickness = 2.0;

	tibiaImplant.init(apLinePclPoint, apLineTuberPoint, sidePoint, exteriorPoint, planeSidePoint, tibiaInfo);

	UKA::IMPLANTS::TibiaImplantMatch tibiaImplantMatch;
	tibiaImplantMatch.init(tibiaImplant, myKnee);

	vtkSmartPointer<vtkPolyData> newImplantTibia = TestVTK::TransformPoly(tibiaModel, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());
	
	std::vector<vtkSmartPointer<vtkPolyData>> polyList;
	polyList.push_back(newImplantTibia);
	TestVTK::show(myKnee.GetTibiaPoly(), polyList);
}

void MatchEasyPKA()
{
	// Error 4 is right
	std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data"; // right
	UKA::IMPLANTS::Knee myKnee = CreateKneeFromFile_NumbersPKA(folder, UKA::IMPLANTS::KRight, UKA::IMPLANTS::KMedial);

	//vtkSmartPointer<vtkPolyData> polyTibiaImplant, polyPatellaImplant;

	UKA::IMPLANTS::FemurImplant* femurImplant = new UKA::IMPLANTS::FemurImplantOnePlane();

	UKA::IMPLANTS::TibiaImplant tibiaImplant;

	//////////////////////////////////////////////Implants

	UKA::IMPLANTS::Plane pPosterior;
	/*
	"D:\\sovajo\\Implants\\pka\\femoral.STL"
	pPosterior.init(UKA::IMPLANTS::Point(-7.36697, 16.0254, -8.73743), UKA::IMPLANTS::Point(-7.85389, 16.0254, 8.8855), UKA::IMPLANTS::Point(-14.8693, 16.0254, 0.351634));
	UKA::IMPLANTS::Point pRodBasePoint(17.1906, -2.90042, -1.59991);
	UKA::IMPLANTS::Point pRodTopPoint(-1.58235, 0.0278676, -0.0651274);
	std::vector<UKA::IMPLANTS::Point> pSideBorder1 = { UKA::IMPLANTS::Point(6.1356, 14.2949, 9.61137), UKA::IMPLANTS::Point(14.2886, 1.29773, 9.60394),
													   UKA::IMPLANTS::Point(11.6988, -9.03007, 9.62781)};

	std::vector<UKA::IMPLANTS::Point> pSideBorder2 = { UKA::IMPLANTS::Point(6.1543, 14.2949, -9.63153), UKA::IMPLANTS::Point(14.3147, 1.53173, -9.65267),
													   UKA::IMPLANTS::Point(11.5081, -9.33871, -9.64367)};*/


	pPosterior.init(UKA::IMPLANTS::Point(-1.98705693113, 16.0609375982, 8.76559347481), UKA::IMPLANTS::Point(-2.50396500712, 16.0613293561, -8.94437729719), UKA::IMPLANTS::Point(-7.92687551969, 16.0599612205, 8.9119337399));
	
	UKA::IMPLANTS::Point pRodBasePoint(17.1500611935, -3.64306703195, 0.0115669254999);
	UKA::IMPLANTS::Point pRodTopPoint(-1.57067067764, 0.0963311269987, -0.0533355484788);
	std::vector<UKA::IMPLANTS::Point> pSideBorder1 = { UKA::IMPLANTS::Point(8.48264081644,13.1976585493,-8.74732637542), UKA::IMPLANTS::Point(10.9804957096,10.5705588957,-8.91229039377),
													   UKA::IMPLANTS::Point(12.8703327684,7.76764520155,-8.80442944992), UKA::IMPLANTS::Point(14.1367283503,4.61611911302,-8.77980126905) };

	std::vector<UKA::IMPLANTS::Point> pSideBorder2 = { UKA::IMPLANTS::Point(8.46616178139,13.336942077,8.47133670984), UKA::IMPLANTS::Point(10.8359483977,10.9111407919,8.6214219954),
													   UKA::IMPLANTS::Point(12.9502176686,7.6594488483,8.67962746283), UKA::IMPLANTS::Point(14.174456106,4.73284273416,8.55774650189) };



	auto femurModel = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Cases_Plan_PKA\\UKA-data\\implants\\femur.STL");
	auto tibiaModel = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\tibia.STL");

	UKA::IMPLANTS::FemurImplantInfo femurInfo;
	femurInfo.femurDistalThickness = 1.0;
	femurInfo.femurPosteriorThickness = 2.0;

	((UKA::IMPLANTS::FemurImplantOnePlane*)femurImplant)->init(pPosterior, pRodBasePoint, pRodTopPoint, pSideBorder1, pSideBorder2, femurModel, femurInfo);

	UKA::IMPLANTS::Point apLinePclPoint(10.35, 11.28, 1.15);
	UKA::IMPLANTS::Point apLineTuberPoint(10.01, -39.601, 1.06);
	UKA::IMPLANTS::Point sidePoint(-15.2245, -12.73, 1.15);
	UKA::IMPLANTS::Point exteriorPoint(-14.5152, -10.8294, 4.41175);
	UKA::IMPLANTS::Point planeSidePoint(12.10, -12.08, 2.59);
	UKA::IMPLANTS::TibiaImplantInfo tibiaInfo;
	tibiaInfo.tibiaThickness = 2.0;

	tibiaImplant.init(apLinePclPoint, apLineTuberPoint, sidePoint, exteriorPoint, planeSidePoint, tibiaInfo);

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

	std::vector<PointTypeITK> hull1;

	itk::Rigid3DTransform<double>::Pointer transformOut = itk::VersorRigid3DTransform<double>::New();

	std::vector<cv::Point3d> tPoints, tPoints2;

	std::cout << transformFemur << std::endl;

	std::cout << "base: " << transformFemur->GetMatrix().GetVnlMatrix().get_column(0) << std::endl;

	try
	{

		hull1 = femurImplantMatch.GetHullPointsOnePlane(transformFemur, transformOut, UKA::IMPLANTS::FemurImplantMatch::KOnePlaneAnteriorAndDistalCurve, 10, 2, 200);
		//hull1 = tibiaImplantMatch.GetHullPoints(transformTibia, transformOut, 1, 1, 1.0, 0.8);
		//std::cout << transformOut << std::endl;

		for (int i = 0; i < hull1.size(); i++)
		{
			Point myPoint(hull1[i][0], hull1[i][1], hull1[i][2]);

			tPoints.push_back(myPoint);
		}
		std::cout << " ********* Hull size: " << hull1.size() << std::endl;
	}
	catch (const ImplantsException& e)
	{
		Test::myPrint("Error get hull");
		std::cout << e.what() << std::endl;
	}

	vtkSmartPointer<vtkPolyData> newImplantFemur = TestVTK::TransformPoly(femurModel, femurImplantMatch.GetRotationMatrix(), femurImplantMatch.GetTranslationMatrix());

	vtkSmartPointer<vtkPolyData> newImplantTibia = TestVTK::TransformPoly(tibiaModel, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());

	std::vector<vtkSmartPointer<vtkPolyData>> polyList;
	polyList.push_back(newImplantFemur);
	polyList.push_back(newImplantTibia);
	TestVTK::show(myKnee.GetFemurPoly(), polyList);

	///////////////////////////////////////////////////// Implant Match Info

	UKA::IMPLANTS::ImplantsMatchFinalInfo info(&myKnee, femurImplant, tibiaImplant, transformFemur, transformTibia);
	//auto transformFemurFinal = info.FemurImplantToTibiaImplant();
	auto transformTibiaFinal = info.TibiaImplantToFemurImplant();

	//vtkSmartPointer<vtkPolyData> newImplantFemurFinal = TestVTK::TransformPoly(femurModel, transformFemurFinal->GetMatrix(), transformFemurFinal->GetOffset());
	vtkSmartPointer<vtkPolyData> newImplantTibiaFinal = TestVTK::TransformPoly(tibiaModel, transformTibiaFinal->GetMatrix(), transformTibiaFinal->GetOffset());

	std::vector<vtkSmartPointer<vtkPolyData>> polyListFinal;
	polyListFinal.push_back(newImplantFemur);
	polyListFinal.push_back(newImplantTibiaFinal);
	TestVTK::show(myKnee.GetFemurPoly(), polyListFinal);

	
	//TestVTK::show(myKnee.GetFemurPoly(), tPoints, true);

}

void videoToImage(const std::string& videoPathIn, const std::string& videoPathOut)
{
	cv::VideoCapture cap(videoPathIn);

	if (!cap.isOpened())
	{
		std::cout << "Error open video" << std::endl;
		return;
	}

	int frameNumber = 0;
	cv::Mat frame;

	typedef unsigned char PixelType;
	typedef itk::Image<PixelType, 2> ImageType;

	while (cap.read(frame)) {
		cvtColor(frame, frame, cv::COLOR_BGR2GRAY);

		int width = frame.cols;
		int height = frame.rows;

		ImageType::Pointer itkImage = ImageType::New();
		ImageType::RegionType region;
		ImageType::IndexType start;
		start.Fill(0);

		ImageType::SizeType size;
		size[0] = width;  // X dimension
		size[1] = height; // Y dimension

		region.SetSize(size);
		region.SetIndex(start);

		itkImage->SetRegions(region);
		itkImage->Allocate();

		PixelType* itkBufferPointer = itkImage->GetBufferPointer();
		std::memcpy(itkBufferPointer, frame.data, width * height * sizeof(PixelType));

		typedef itk::ImageFileWriter<ImageType> WriterType;
		WriterType::Pointer writer = WriterType::New();

		std::ostringstream filename;
		filename << videoPathOut << "\\frame_" << frameNumber << ".dcm";
		writer->SetFileName(filename.str());
		writer->SetInput(itkImage);

		try {
			writer->Update();
		}
		catch (itk::ExceptionObject &error) {
			std::cout << "Error writing the image: " << error << std::endl;
		}

		std::cout << "Frame " << frameNumber << " save as " << filename.str() << std::endl;
		frameNumber++;
	}

}


void executeBalance()
{
	/*
		In order to work with this library you must create the objects Knee and LegAngle.
		These are the steps:
		1-) You should have the following points get from the CT image (These points are getting  by the doctor during the pre-planning operation.).
			hipCenter, femurKneeCenter, lateralEpicondyle, medialEpicondyle, medialCondyle,
			lateralCondyle, anteriorCortex, ankleCenter, lateralPlateau, medialPlateau,
			tibiaKneeCenter, tibiaTubercle, pclCenter.
			With each of these points you must create an object of class Point(x, y, z);

		2-) You should have a pointer to a 3d volume of the femur and another of the tibia.

		** With all this data you can start creating the objects how is explained.

		1-) Now you must create the objects (Knee and LegAngle)

			- Knee : Receive the points obtained from the CT images,
					the characteristics of the implant in the ImplantInfo structure and
					the pointers to the 3d volumes of the femur and tibia, as shown in the example below.

			- LegAngle : It has an empty constructor, it does not receive parameters.

		3-) Finally, you can create the Balance Object.
			- Receive Knee objects and LegAngle objects.

		Note: It is so important update Balance object transform matrix before to use Balance methods.
			  Every time you use the methods of the Balance object, you must update the transformation
			  matrix with the methods setTransformFemurCtToCamera(double matrix[16]) and
			  setTransformTibiaCtToCamera(double matrix[16]);
	*/

	FemurImplantInfo implantFemur;
	implantFemur.femurPosteriorLateralThickness = 9;
	implantFemur.femurDistalLateralThickness = 9;
	implantFemur.femurPosteriorMedialThickness = 9;
	implantFemur.femurDistalMedialThickness = 9;

	TibiaImplantInfo implantTibia;
	implantTibia.tibiaLateralThickness = 11;
	implantTibia.tibiaMedialThickness = 11;

	ImplantImageType::Pointer femur, tibia;
	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly, femurImplantSTL;

	Test::readImage<ImplantImageType>(dicom_femur, femur);
	Test::readImage<ImplantImageType>(dicom_tibia, tibia);

	femurPoly = TestVTK::ReadPolyData(dicom_femur_vtk);
	tibiaPoly = TestVTK::ReadPolyData(dicom_tibia_vtk);
	femurImplantSTL = TestVTK::ReadPolyData(femur_implant_poly);

	Point hipCenter(57.2, -47.89, 1158.81);
	Point femurKneeCenter(23.86, -55.92, 771.39);
	Point lateralEpicondyle(-11.02, -48.33, 782.79);
	Point medialEpicondyle(67.73, -49.17, 777.39);
	Point anteriorCortex(17.11, -73.36, 816.39);
	Point medialCondyle(48.61, -20.77, 771.39);
	Point lateralCondyle(2.77, -23.3, 778.59);

	Point lateralPlateau(3.89, -41.30, 743.80);
	Point medialPlateau(41.86, -39.05, 740.80);
	Point tibiaPCL(18.80, -24.70, 733.0);
	Point tibiaKneeCenter(17.67, -53.11, 745.00);
	Point tibiaTubercle(16.27, -74.77, 700.60);
	Point lateralAnkle(-19.30, -74.89, 384.90);
	Point medialAnkle(40.89, -83.05, 396.90);
	Point ankleCenter = Knee::getComputeAnkleCenter(lateralAnkle, medialAnkle);


	Point planeA1(-24.6383, 55.6821, 34.4108);
	Point planeA2(0, 55.8612, 24.148);
	Point planeA3(24.6383, 55.6821, 34.4108);

	Point planeB1(-31.33, 56.294, 15.783);
	Point planeB2(-31.33, 63.8159, 8.39129);
	Point planeB3(31.33, 63.8159, 8.39129);

	Point planeC1(-31.2913, 64.5168, 8.10454);
	Point planeC2(-29.421, 83.7645, 8.10454);
	Point planeC3(28.2063, 83.7645, 8.10454);

	Point planeD1(-29.3421, 84.4655, 8.39129);
	Point planeD2(-27.0605, 98.6397, 22.3203);
	Point planeD3(22.6364, 98.6397, 22.3203);

	Point planeE1(-26.967, 98.935, 22.9464);
	Point planeE2(22.404, 98.935, 22.9464);
	Point planeE3(-14.8076, 101.72, 54.7848);

	Point P1(-31.33, 56.0022, 16.0697);
	Point P2(31.33, 56.0022, 16.0697);
	Point P3(-31.33, 64.1077, 8.10454);
	Point P4(31.33, 64.1077, 8.10454);
	Point P5(-29.3813, 84.1737, 8.10454);
	Point P6(28.1412, 84.1737, 8.10454);
	Point P7(-14.8076, 101.72, 54.7848);

	Point implantPCLP1(9.45, -15.18, -3.02);
	Point implantPCLP2(-9.45, -15.18, -3.02);
	Point implantFront(0, 28.35, 3.48);

	Point firstPlateau(-22.31, -1.25, 7.0);
	Point secondPlateau(22.37, -1.24, 7.0);

	Point kneeCap(0, 0, 0);

	Plane A, B, C, D, E;
	A.init(planeA1, planeA2, planeA3);
	B.init(planeB1, planeB2, planeB3);
	C.init(planeC1, planeC2, planeC3);
	D.init(planeD1, planeD2, planeD3);
	E.init(planeE1, planeE2, planeE3);


	///////////////////////////////////////////

	std::vector<PointTypeITK> kneeCapPath;

	std::vector<Point> tPoints = { Point(-5.83, 103.52, 52.29), Point(-4.83, 103.52, 48.29), Point(-3.83, 103.52, 45.29), Point(-3.83, 103.52, 36.29),
								   Point(-2.83, 103.52, 31.29), Point(-1.83, 102.52, 26.29), Point(-0.83, 100.52, 21.29), Point(-0.83, 98.52, 18.29),
								   Point(-0.83, 95.52, 15.29), Point(0.17, 93.52, 13.29) };

	for (int i = 0; i < tPoints.size(); i++)
	{
		PointTypeITK itkPoint;
		itkPoint[0] = tPoints[i].x;
		itkPoint[1] = tPoints[i].y;
		itkPoint[2] = tPoints[i].z;
		kneeCapPath.push_back(itkPoint);
	}

	/////////////////////////////////////////////

	FemurImplant femurImplant;
	TibiaImplant tibiaImplant;

	femurImplant.init(A, B, C, D, E, P3, P4, P7, femurImplantSTL, kneeCapPath, implantFemur);

	tibiaImplant.init(implantPCLP1, implantPCLP2, implantFront, firstPlateau, implantTibia);

	Patella myPatella;
	myPatella.init(kneeCap, lateralEpicondyle, medialEpicondyle, femurKneeCenter, femurPoly); // Just for avoid problem with patella

	Knee knee;
	knee.init(hipCenter, anteriorCortex, femurKneeCenter, lateralEpicondyle, medialEpicondyle, 
		tibiaKneeCenter, tibiaTubercle, tibiaPCL, ankleCenter, myPatella, femurPoly, tibiaPoly, KneeSideEnum::KRight);

	auto femurImplantMatch = new FemurImplantMatch();
	auto tibiaImplantMatch = new TibiaImplantMatch();
	try
	{
		femurImplantMatch->init(femurImplant, knee);
		tibiaImplantMatch->init(tibiaImplant, knee);
	}
	catch (const ImplantsException& e)
	{
		std::cout << e.what() << std::endl;
	}

	itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer tibiaTransformIn = itk::VersorRigid3DTransform<double>::New();

	femurTransformIn->SetMatrix(femurImplantMatch->GetRotationMatrix());
	femurTransformIn->SetOffset(femurImplantMatch->GetTranslationMatrix());

	tibiaTransformIn->SetMatrix(tibiaImplantMatch->GetRotationMatrix());
	tibiaTransformIn->SetOffset(tibiaImplantMatch->GetTranslationMatrix());

	Balance balance;
	balance.init(knee, femurImplant, tibiaImplant, femurTransformIn, tibiaTransformIn);

	//////////////////////////

	itk::Rigid3DTransform<double>::Pointer famurCtToMarker = Test::generateRigid3DTransform();
	itk::Rigid3DTransform<double>::Pointer tibiaCtToMarker = Test::generateRigid3DTransform();

	itk::Rigid3DTransform<double>::Pointer famurMarkerToCamera = Test::generateRigid3DTransform();
	itk::Rigid3DTransform<double>::Pointer tibiaMarkerToCamera = Test::generateRigid3DTransform();

	balance.setTransformFemurCtToMarker(famurCtToMarker);

	balance.setTransformTibiaCtToMarker(tibiaCtToMarker);

	balance.setTransformFemurMarkerToCamera(famurMarkerToCamera);
	balance.setTransformTibiaMarkerToCamera(tibiaMarkerToCamera);

	BalanceInfo info = balance.distanceByAngleBeforeResectionBone();

	BalanceInfo info2 = balance.distanceByAngleWithImplant(Registration::makeItkPoint(3.89, -41.30, 743.80), Registration::makeItkPoint(41.86, -39.05, 740.80));

	std::cout << "Lateral Distance:  " << info.distanceLateral << "  Medial Distance: " << info.distanceMedial << std::endl;

	std::cout << " Lateral Distance:   " << info2.distanceLateral << "  Medial Distance: " << info2.distanceMedial << std::endl;

	AnglesInfo angles = balance.anglesByMotion();
	std::cout << "Flexion: " << angles.flexion_angle << "  Varus: " << angles.varus_angle << " Rotation:" << angles.rotation_angle << std::endl;

	std::cout << "Patella:  " << balance.getKneeCapPoint() << std::endl;

	//TestVTK::show(balance.GetFemurPolyCut(), balance.GetTibiaPolyCut());


	std::cout << "-----------------------------------------Info ----------------------------: " << std::endl;

	//BoneRotulaPath obj(knee);

	//vtkSmartPointer<vtkPolyData> allSurface = obj.JoinSagitalPointTest();

	//SliceBorder::ShowPolyData(allSurface, allSurface);

	//return;


	ImplantsMatchFinalInfo myInfo(&knee, femurImplant, tibiaImplant, femurTransformIn, tibiaTransformIn);

	/* myInfo.SetThicknessFemurImplant(5.0, 8.0);
	 myInfo.SetThicknessTibiaImplant(5.0, 4.0);
	 myInfo.SetCartilageFemur(2.0);
	 myInfo.SetCartilageTibia(2.0);*/
	myInfo.test();
	std::cout << "*********************************" << std::endl;

	std::vector<PointTypeITK> bonePath, impantPath;

	bonePath = myInfo.getBoneKneeCapPath();

	impantPath = myInfo.getImplantKneeCapPath();

	vtkNew<vtkPoints> bonePolyPoints, implantPolyPoints;

	for (int i = 0; i < bonePath.size(); i++)
	{
		//std::cout << bonePath[i] << std::endl;
		double pnt[3];
		pnt[0] = bonePath[i][0];
		pnt[1] = bonePath[i][1];
		pnt[2] = bonePath[i][2];

		bonePolyPoints->InsertNextPoint(pnt);
	}
	//std::cout << "********************************************" << std::endl;
	for (int i = 0; i < impantPath.size(); i++)
	{
		//std::cout << impantPath[i] << std::endl;

		double pnt[3];
		pnt[0] = impantPath[i][0];
		pnt[1] = impantPath[i][1];
		pnt[2] = impantPath[i][2];

		implantPolyPoints->InsertNextPoint(pnt);
	}

	vtkSmartPointer<vtkPolyData> bonePolyLine = TestVTK::CreatePolyLine(bonePolyPoints);
	vtkSmartPointer<vtkPolyData> implantPolyLine = TestVTK::CreatePolyLine(implantPolyPoints);

	vtkNew<vtkAppendPolyData> filter;
	filter->AddInputData(bonePolyLine);
	filter->AddInputData(implantPolyLine);
	filter->Update();

	/*vtkSmartPointer<vtkPolyData> topPointPath = TestVTK::CreateSphereTest(knee.getTopPointOnPatellaPath());
	vtkNew<vtkAppendPolyData> filter;

	filter->AddInputData(knee.GetFemurPoly());
	filter->AddInputData(topPointPath);
	filter->Update();*/

	//////////////////////////////////////////////////////////////////////////descomentar
	BoneRotulaPath obj(knee);
	vtkSmartPointer<vtkPolyData> allSurface = obj.JoinSagitalPointTest();
	SliceBorder::ShowPolyData(allSurface, filter->GetOutput());
	///////////////////////////////////////////////////////////////////////////

	/* SegmentImageType::Pointer imgOut;
	 Test::readImage<SegmentImageType>(dicom_femur, imgOut);
	 SliceBorder obj;
	 vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData(dicom_femur_vtk);

	 Point inferiorLateral = knee.getLateralInferiorFemurPoint();
	 Point inferiorMedial = knee.getMedialInferiorFemurPoint();

	 double pnt[3] = { inferiorMedial.x, inferiorMedial.y, inferiorMedial.z};

	 vtkNew<vtkSphereSource> sphere;
	 sphere->SetCenter(pnt);
	 sphere->SetRadius(1);
	 sphere->Update();

	 vtkSmartPointer<vtkPolyData> poly;
	 poly = sphere->GetOutput();

	 vtkNew<vtkAppendPolyData> filter;
	 filter->AddInputData(BoneVTK);
	 filter->AddInputData(poly);
	 filter->Update();

	 SliceBorder::ShowPolyData(filter->GetOutput(), filter->GetOutput());
 */
}

void ExecuteKneeCap()
{
	ImplantImageType::Pointer femur, tibia;
	Test::readImage<ImplantImageType>(femur_path, femur);
	Test::readImage<ImplantImageType>(tibia_path, tibia);

	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly, femurImplantPoly, patellaPoly;

	femurPoly = TestVTK::ReadPolyData(femur_path_poly);
	tibiaPoly = TestVTK::ReadPolyData(tibia_path_poly);

	femurImplantPoly = TestVTK::ReadPolyData(femur_implant_poly);

	Point hipCenter(28.8433, -57.3193, 498.549);
	Point femurKneeCenter(-10.35, -8.79, 95.65);
	Point lateralEpicondyle(-45.6781, -17.499, 106.142);
	Point medialEpicondyle(28.9137, -13.2568, 106.36);
	Point anteriorCortex(-20.225, 13.6104, 138.595);
	Point medialCondyle(13.359, -39.0635, 98.9912);
	Point lateralCondyle(-31.891, -38.3564, 101.191);

	Point lateralPlateau(14.7731, -11.4893, 74.239);
	Point medialPlateau(-18.4574, -6.54004, 75.3391);
	Point tibiaPCL(-6.79137, -29.5186, 67.0884);
	Point tibiaKneeCenter(-6.79137, 2.29785, 73.689);
	Point tibiaTubercle(1.693, 15.7314, 50.0369);
	//Point lateralAnkle(-19.30, -74.89, 384.90);
	//Point medialAnkle(40.89, -83.05, 396.90);
	Point ankleCenter(2.21859, -7.30133, -291.76);
	Point kneeCap(0, 0, 0);

	auto startTime = std::chrono::system_clock::now();

	Patella myPatella;
	myPatella.init(kneeCap, lateralEpicondyle, medialEpicondyle, femurKneeCenter, femurPoly); // Just for avoid problem with patella

	Knee knee;
	knee.init(hipCenter, anteriorCortex, femurKneeCenter, lateralEpicondyle, medialEpicondyle, 
		tibiaKneeCenter, tibiaTubercle, tibiaPCL, ankleCenter, myPatella, femurPoly, tibiaPoly, KneeSideEnum::KRight);

	auto endTime = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = endTime - startTime;
	std::cout << "elapsed time genaral knee: " << elapsed_seconds.count() << std::endl;

	//vtkSmartPointer<vtkPolyData> topPointPath = TestVTK::CreateSphereTest(knee.getTopPointOnPatellaPath());
	vtkNew<vtkAppendPolyData> filter;

	filter->AddInputData(knee.GetFemurPoly());
	//filter->AddInputData(topPointPath);
	filter->Update();


	BoneRotulaPath obj(knee);
	vtkSmartPointer<vtkPolyData> surface = obj.JoinSagitalPointTest();
	SliceBorder::ShowPolyData(filter->GetOutput(), surface);

	//SliceBorder::ShowPolyData(filter->GetOutput(), surface);
}

void DrawBonePlanes()
{
	/*
	cv::Point3d normalA(3.89014417579818, 505.394071270987, -19.6718255985779);
	double biasA = 31632.6251725764;
	LitePlane myPlaneA(normalA, biasA);

	cv::Point3d normalB(-36.1460400919554, 453.304036093569, -479.444350325885);
	double biasB = 386381.737013022;
	LitePlane myPlaneB(normalB, biasB);

	cv::Point3d normalC(4.00148899628732, 0.963765946016413, 46.4984063269836);
	double biasC = -35644.7260922921;
	LitePlane myPlaneC(normalC, biasC);

	cv::Point3d normalD(-66.7441356127742, -706.582140677784, -686.781748446354);
	double biasD = 484410.628959288;
	LitePlane myPlaneD(normalD, biasD);

	cv::Point3d normalE(26.2320591123625, 1574.30241308204, 103.17106472463);
	double biasE = 37587.0940589373;
	LitePlane myPlaneE(normalE, biasE);

	cv::Point3d normalTibia(10.560898922732, 128.712684730322, 789.507132458888);
	double biasTibia = -576396.794215868;
	//test
	/*cv::Point3d normalTibia(79.673925739382, 117.453098827502, 787.30987237116);
	double biasTibia = - 576248.615990424;*/
	/*
	LitePlane myPlaneTibia(normalTibia, biasTibia);

	cv::Point3d normalMid(-62.4267030138798, 0.689064939177189, 5.35793992829788);
	double biasMid = -2605.02763597974;
	LitePlane myPlaneMid(normalMid, biasMid);

	myPlaneA.setName("points_A");
	myPlaneB.setName("points_B");
	myPlaneC.setName("points_C");
	myPlaneD.setName("points_D");
	myPlaneE.setName("points_E");
	myPlaneTibia.setName("tibia_plane");

	std::vector<RegistrationImageType::Pointer> imgPoints;
	RegistrationImageType::Pointer imgOut;
	readImage<RegistrationImageType>(dicom_femur, imgOut);
	imgPoints.push_back(imgOut);
	Registration regis(imgPoints);
	std::vector<PointTypeITK> PointsITK = regis.getPointsCT();

	pcl::PointCloud<pcl::PointXYZ>::Ptr Points = itkVectorToPCL(PointsITK);
	pcl::PointCloud<pcl::PointXYZ>::Ptr redPoints(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr redPoints2(new pcl::PointCloud<pcl::PointXYZ>);

	auto it1 = Points->points.begin();
	auto it2 = Points->points.end();

	cv::Point3d temp;
	LitePlane current = myPlaneD;
	for (; it1 != it2; ++it1)
	{
		temp = cv::Point3d((*it1).x, (*it1).y, (*it1).z);
		if (current.isPointNear(temp))
		{
			redPoints->points.push_back(pcl::PointXYZ(temp.x, temp.y, temp.z));
		}
	}

	//Registration::drawCloud(Points, redPoints, true);

	std::string myFile = current.name;
	std::ifstream infile(myFile);
	double a, b, c;
	while (infile >> a >> b >> c)
	{
		redPoints2->points.push_back(pcl::PointXYZ(a, b, c));
	}

	//Registration::drawCloud(redPoints, redPoints2, true);*/
}

void checkAngles()
{
	Point hipCenter(57.2, -47.89, 1158.81);
	Point femurKneeCenter(23.86, -55.92, 771.39);
	Point lateralEpicondyle(-11.02, -48.33, 782.79);
	Point medialEpicondyle(67.73, -49.17, 777.39);
	Point anteriorCortex(17.11, -73.36, 816.39);
	Point medialCondyle(48.61, -20.77, 771.39);
	Point lateralCondyle(2.77, -23.3, 778.59);

	Point lateralPlateau(3.89, -41.30, 743.80);
	Point medialPlateau(41.86, -39.05, 740.80);
	Point tibiaPCL(18.80, -24.70, 733.0);
	Point tibiaKneeCenter(17.67, -53.11, 745.00);
	Point tibiaTubercle(16.27, -74.77, 700.60);
	Point lateralAnkle(-19.30, -74.89, 384.90);
	Point medialAnkle(40.89, -83.05, 396.90);

	Point Vf = hipCenter - femurKneeCenter;
	Vf.normalice();

	Plane coronal, sagital, transverse;
	transverse.init(Vf, femurKneeCenter);

	Point P1, P2;
	P1 = transverse.getProjectionPoint(lateralEpicondyle);
	P2 = transverse.getProjectionPoint(medialEpicondyle);
	Point Vs = P1 - P2;
	Vs.normalice();

	sagital.init(Vs, femurKneeCenter);

	Point Vc = Vf.cross(Vs);
	Vc.normalice();
	coronal.init(Vc, femurKneeCenter);

	Point transProjCor = coronal.getProjectionVector(Vf);
	transProjCor.normalice();

	Point corProjTrans = transverse.getProjectionVector(Vc);
	corProjTrans.normalice();

	std::cout << "Coronal Normal: " << Vc << std::endl;
	std::cout << "**Proj Coronal: " << corProjTrans << std::endl;
	std::cout << "\n" << std::endl;
	std::cout << "Transverse Normal: " << Vf << std::endl;
	std::cout << "**Proj Transverse: " << transProjCor << std::endl;


}

void getSlicesAndSurfaceBone()
{

	cv::Point3d min(40.7627, -128.052, -656.507);
	cv::Point3d center(77.4238, -98.1719, -590.012);
	cv::Point3d mainMax(114.085, -68.292, -523.517);
	cv::Point3d maxX(114.085, -128.052, -656.507);

	cv::Point3d diff = center - min;
	double distance = sqrt(diff.dot(diff));
	cv::Point3d vector = diff / distance;
	cv::Point3d max = min + 2.0 * distance * vector;

	SPlane sp1, sp2;
	sp1.init(diff, min);
	sp2.init(diff, max);


	SegmentImageType::Pointer imgOut;
	std::string imgName = "D:\\3D_DICOM\\right_femur_segment.nrrd";
	Test::readImage<SegmentImageType>(imgName, imgOut);

	SliceBorder obj;
	vtkSmartPointer<vtkPolyData> image = TestVTK::ReadPolyData(dicom_femur_vtk);

	std::vector<vtkSmartPointer<vtkPolyData>> slices = obj.MakeSlicesFromPolyData(image, 1, 1.0);



	//std::vector<vtkSmartPointer<vtkPolyData>> slices = obj.MakeSlicesFromPolyData(image, sp1, sp2, 2.0);
	//vtkSmartPointer<vtkPolyData> surface = obj.GetSurface3D(slices);
	//vtkSmartPointer<vtkPolyData> plane = SliceBorder::GetCenterSlice(surface, 3);

	vtkNew<vtkAppendPolyData> filter;
	for (int i = 0; i < slices.size(); i++)
	{
		filter->AddInputData(slices[i]);
		std::cout << "Original line number: " << slices[i]->GetNumberOfPoints() << std::endl;
	}
	filter->Update();

	/* vtkSegmentationConverterRule * slicer = new vtkPlanarContourToClosedSurfaceConversionRule();

	 auto surface2 = filter->GetOutput();
	 vtkNew<vtkPolyData> closedSurfacePolyData;
	 bool result;
	 result = slicer->Convert(surface2, closedSurfacePolyData);

	if (result == true)
	 {
		 std::cout << "Todo bien*******" << std::endl;
	 }
	 else
	 {
		 std::cout << "Todo mallll*******" << std::endl;
	 }*/
	 /* vtkNew<vtkPolyData> closedSurfacePolyData;

	  if (MySurface::Convert(surface2, closedSurfacePolyData))
	  {
		  std::cout << "Fine!!!!!!!!!!!!!!!!!!!!" << std::endl;
	  }*/


	  //SliceBorder::ShowPolyData(surface2, closedSurfacePolyData);

	  //SliceBorder::ShowPolyData(surface, plane);
}

void oldTilling()
{
	SegmentImageType::Pointer imgOut;
	std::string imgName = "D:\\3D_DICOM\\right_femur_segment.nrrd";
	Test::readImage<SegmentImageType>(imgName, imgOut);

	SliceBorder obj;
	vtkSmartPointer<vtkPolyData> image = TestVTK::ReadPolyData(dicom_femur_vtk);

	std::vector<vtkSmartPointer<vtkPolyData>> slices = obj.MakeSlicesFromPolyData(image, 1, 10.0);

	vtkSmartPointer<vtkPoints> tempUp;

	for (int i = 0; i < slices.size() - 1; i++)
	{
		auto slice1 = slices[i];
		auto slice2 = slices[i + 1];

		std::cout << "Iteration: " << i << std::endl;

		vtkNew<vtkAppendPolyData> filter;
		filter->AddInputData(slice1);
		filter->AddInputData(slice2);
		filter->Update();
		auto margeSlices = filter->GetOutput();

		if (i == 0)
		{
			Tiling myTiling(slice1, slice2, 0, true);
			tempUp = myTiling.GetOriginalUp();
			auto surface = myTiling.MakeTiling();
			SliceBorder::ShowPolyData(margeSlices, surface);
		}
		else
		{
			Tiling myTiling(tempUp, slice2, 0, true);
			tempUp = myTiling.GetOriginalUp();

			auto surface = myTiling.MakeTiling();

			SliceBorder::ShowPolyData(margeSlices, surface);
		}
	}

}

void newTiling()
{
	/*SegmentImageType::Pointer imgOut;
	std::string imgName = "D:\\3D_DICOM\\right_femur_segment.nrrd";
	Test::readImage<SegmentImageType>(imgName, imgOut);*/
	/*SliceBorder obj;
	std::vector<vtkSmartPointer<vtkPolyData>> new_slices;
	std::string myPath = "D:\\crash\\Contour\\contour-data-";
	for (int i = 0; i < 74; i++)
	{
		std::string myPathTemp = myPath + std::to_string(i);
		myPathTemp = myPathTemp + ".vtp";

		vtkSmartPointer<vtkPolyData> image1 = TestVTK::ReadPolyDataVTP(myPathTemp);

		std::vector<PointTypeITK> mainPoints = obj.GetMainPoints(image1, 160);
		vtkNew<vtkPoints> myPoints;
		for (int j = 0; j < mainPoints.size(); j++)
		{
			double pnt[3] = { mainPoints[j][0], mainPoints[j][1] , mainPoints[j][2] };
			myPoints->InsertNextPoint(pnt);
		}

		auto spheres = TestVTK::CreatePolyLine(myPoints);
		SliceBorder::ShowPolyData(image1, spheres);
	}*/

	SliceBorder obj;
	//vtkSmartPointer<vtkPolyData> image = TestVTK::ReadPolyData("D:\\3D_DICOM\\Real_Person\\segmentation\\tibia_poly.vtk");//"D:\\Big_Lee_Data\\tibia_seg_d.vtk");
	vtkSmartPointer<vtkPolyData> image = TestVTK::itkImageToSurface3D("D:\\Algorith_Project\\build\\segment_femur_auto.nrrd");


	////obj.SavePolyData(image, "right_femur.vtk");

	//SPlane sp1, sp2;

	//cv::Point3d point1 = { 185.28808593749832, -26.668747139198167, -797.51324391911260 };
	//cv::Point3d point2 = { 185.28808593749832, -134.95683307669816, -557.42999196598760 };
	//cv::Point3d point3 = { 186.28808593749832, -80.812790107948160, -677.47161794255010 };

	//cv::Point3d point4 = { 185.28808593750000, -88.000000000000000, -763.59997558593750 };
	//cv::Point3d point5 = { 185.28808593750000, -196.28808593750000, -523.51672363281250 };
	//cv::Point3d point6 = { 186.28808593750000, -142.14404296875000, -643.55834960937500 };

	//sp1.init(point1, point2, point3);
	//sp2.init(point4, point5, point6);

	//

	std::vector<vtkSmartPointer<vtkPolyData>> slices = obj.MakeSlicesFromPolyData(image, 1, 2);
	//std::vector<vtkSmartPointer<vtkPolyData>> slices_test = obj.MakeSlicesFromPolyData(image, 1, 4.0);
	//bool result;
	//std::string msg;
	//auto startTime = std::chrono::system_clock::now();

	//obj.GetMainPoints(slices[10]);

	//auto endTime = std::chrono::system_clock::now();

	/////////////////////////////////////////////////////////////////


	std::vector<vtkSmartPointer<vtkPolyData>> slices_test;
	for (int i = 0; i < slices.size(); i++)
	{
		/*vtkSmartPointer<vtkPolyData> poly = slices[i];
		vtkNew<vtkCardinalSpline> cardinalSpline;
		vtkNew<vtkSplineFilter> spline;
		spline->SetSpline(cardinalSpline);
		spline->SetInputData(poly);
		spline->SetSubdivideToLength();
		spline->SetLength(1);
		spline->Update();*/

		slices_test.push_back(TestVTK::CreateClosePolyLine((slices[i])->GetPoints()));
	}

	//TestVTK::show(slices_test[0], slices_test[20]);
	////////////////////////////////////////////////////////////////


	bool result;
	std::string msg;
	double normal[3] = { 0, -0.707091, 0.707123 };
	auto surface = obj.GetSurface3D(slices, result, msg);
	//auto surface_spline = obj.GetSurface3D(slices_test, result, msg);
	auto smoothSurface = obj.SmoothSurface(surface);

	/*std::cout << "Slices number: " << slices.size() << " Result: " << result << " msg: " << msg << std::endl;
	std::chrono::duration<double> elapsed_seconds = endTime - startTime;
	std::cout << "elapsed time genaral: " << elapsed_seconds.count() << std::endl;*/

	//TestVTK::SavePolyData(surface, "surface_raw.vtk");
	//TestVTK::SavePolyData(smoothSurface, "surface_smooth.vtk");

	vtkNew<vtkAppendPolyData> filter;
	for (int i = 0; i < slices.size(); i++)
	{
		/*std::string slicesName = "slice_" + std::to_string(i) + ".vtk";
		TestVTK::SavePolyData(slices[i], slicesName);*/

		filter->AddInputData(slices[i]);
	}

	//filter->AddInputData(surface);

	filter->Update();

	//SliceBorder::ShowPolyData(filter->GetOutput(), smoothSurface);
	TestVTK::show(surface, smoothSurface);

}

void testTree(int n)
{
	vtkNew<vtkPoints> down, up;
	for (int i = 1; i <= n; i++)
	{
		double pntDown[3];
		double pntUp[3];
		pntDown[0] = i;
		pntDown[1] = 1;
		pntDown[2] = 0;

		pntUp[0] = i;
		pntUp[1] = 2;
		pntUp[2] = 0;

		down->InsertNextPoint(pntDown);
		up->InsertNextPoint(pntUp);
	}

	std::vector<std::pair<int, int>> path;
	GenerateTree myTree(down, up);
	path = myTree.GetPath();

	for (int i = 0; i < path.size(); i++)
	{
		std::cout << "( " << path[i].first << ", " << path[i].second << " )" << std::endl;
	}
}

void LS_Registration_Femur()
{
	SegmentImageType::Pointer imgOut;
	Test::readImage<SegmentImageType>(dicom_femur, imgOut);

	vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData(dicom_femur_vtk);

	PointTypeITK hipCenter, kneeCenter, lateralEpi, medialEpi, cortex, lateralCondyle, medialCondyle;
	PointTypeITK hipCenter2, kneeCenter2, lateralEpi2, medialEpi2, cortex2;

	hipCenter = Registration::makeItkPoint(55.0469, -47.3281, 1158.21);

	cortex = Registration::makeItkPoint(19.6406, -72.5156, 825.392);

	lateralEpi = Registration::makeItkPoint(-10.7344, -47.4844, 780.994);

	medialEpi = Registration::makeItkPoint(66.6094, -47.7656, 770.837);

	kneeCenter = Registration::makeItkPoint(23.0156, -56.4844, 767.795);

	std::vector<PointTypeITK> bones_points;

	std::string point_path = "D:\\Registro\\Points\\registro_femur.txt";
	std::ifstream infile(point_path);
	double a, b, c;

	while (infile >> a >> b >> c)
	{
		PointTypeITK itkPoint;
		itkPoint[0] = a;
		itkPoint[1] = b;
		itkPoint[2] = c;
		bones_points.push_back(itkPoint);
	}

	hipCenter2 = Registration::makeItkPoint(324.73, 302.025, -1559.25);

	cortex2 = Registration::makeItkPoint(212.75, 6.7, -1688.75);

	lateralEpi2 = Registration::makeItkPoint(226.73, -23.01, -1731.10);

	medialEpi2 = Registration::makeItkPoint(235.15, -52.78, -1658.5);

	kneeCenter2 = Registration::makeItkPoint(214.92, -47.01, -1702.20);

	std::cout << "Bones Size: " << bones_points.size() << std::endl;

	std::vector<RegistrationImageType::Pointer> imgPoints;
	RegistrationImageType::Pointer knee;
	Test::readImage<RegistrationImageType>(dicom_femur, knee);
	//Test::readImage<RegistrationImageType>("D:\\Kevin errors\\crash_femur.nrrd", knee);
	imgPoints.push_back(knee);

	FemurRegistration * regis = new FemurRegistration(BoneVTK, hipCenter, kneeCenter, medialEpi);

	bool result;
	double error;

	std::cout << "**************** Before land mark: " << std::endl;

	//regis->RegistrationLandmarks(hipCenter2, cortex2, kneeCenter2, lateralEpi2, medialEpi2, error);

	//std::cout << "**************** Land mark Error: " << error << std::endl;

	regis->MakeRegistration(bones_points, hipCenter2, kneeCenter2, medialEpi2, true);

	std::cout << "**************** Femur Error1: " << regis->getError() << std::endl;

	std::cout << "**************** Test new so so " << std::endl;

	//KneeRegistration * regis2 = new KneeRegistration(BoneVTK, BoneVTK, hipCenter, kneeCenter, lateralEpi, medialEpi);

	//double newError = regis2->MakeFemurRegistration(bones_points, hipCenter2, kneeCenter2, lateralEpi2, medialEpi2);

	//std::cout << "**************** Femur Error2: " << newError << std::endl;

	/*
	std::cout << "******************* Transform1: \n" << regis->getTransformMarkerToCt() << std::endl;

	std::cout << "******************* Random ***************************** \n" << std::endl;

	FemurRegistration * regis2 = new FemurRegistration(BoneVTK, hipCenter, cortex, kneeCenter, lateralEpi, medialEpi);
	regis2->MakeRegistration(bones_points, hipCenter2, cortex2, kneeCenter2, lateralEpi2, medialEpi2, true);

	std::cout << "**************** Femur Error2: " << regis2->getError() << std::endl;

	std::cout << "******************* Transform2: \n" << regis2->getTransformMarkerToCt() << std::endl;*/
}

void LS_Registration_Tibia()
{
	PointTypeITK tibiaKnee, tubercle, lateralAnkle, medialAnkle;
	PointTypeITK tibiaKnee2, tubercle2, lateralAnkle2, medialAnkle2;

	tibiaKnee = Registration::makeItkPoint(14.0156, -62.3906, 741.996);

	tubercle = Registration::makeItkPoint(14.5781, -73.6406, 720.996);

	lateralAnkle = Registration::makeItkPoint(-16.4844, -86.1406, 386.801);

	medialAnkle = Registration::makeItkPoint(40.8906, -83.8906, 396.801);

	std::vector<RegistrationImageType::Pointer> imgPoints;
	RegistrationImageType::Pointer knee;
	Test::readImage<RegistrationImageType>(dicom_tibia, knee);
	imgPoints.push_back(knee);

	vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData(dicom_tibia_vtk);

	//////////////////////////////////////////////

	/*Point pNormal = (Test::ITKToCV(lateralAnkle) + Test::ITKToCV(medialAnkle)) / 2.0;
	Point pPoint = Test::ITKToCV(tibiaKnee);

	vtkNew<vtkPlane> plane;
	vtkNew<vtkCutter> cutter;
	cutter->SetInputData(BoneVTK);

	plane->SetNormal(pNormal.x, pNormal.y, pNormal.z);
	plane->SetOrigin(pPoint.x, pPoint.y, pPoint.z);
	cutter->SetCutFunction(plane);
	cutter->Update();

	auto contour = cutter->GetOutput();

	TestVTK::show(contour, contour);

	return;*/

	/////////////////////////////////////////////////

	TibiaRegistration * regis = new TibiaRegistration(BoneVTK, tubercle, lateralAnkle, medialAnkle);

	std::vector<PointTypeITK> bones_points;

	std::string point_path = "D:\\Registro\\Big_Lee_data\\sample_data_on_bone_3\\tibia.ini";
	std::ifstream infile(point_path);
	double a, b, c;
	double scale = 2.0;
	while (infile >> a >> b >> c)
	{
		PointTypeITK itkPoint;
		itkPoint[0] = a * scale;
		itkPoint[1] = b * scale;
		itkPoint[2] = c * scale;
		bones_points.push_back(itkPoint);
	}

	std::cout << "Bones Size: " << bones_points.size() << std::endl;

	std::string land_mark = "D:\\Registro\\Big_Lee_data\\sample_data_on_bone_3\\landmark.ini";
	std::ifstream infile2(land_mark);
	std::vector<PointTypeITK> landMarks;

	while (infile2 >> a >> b >> c)
	{
		PointTypeITK itkPoint;
		itkPoint[0] = a;
		itkPoint[1] = b;
		itkPoint[2] = c;
		landMarks.push_back(itkPoint);
	}

	std::cout << "Land mark Size: " << landMarks.size() << std::endl;

	tibiaKnee2 = landMarks[5];

	tubercle2 = landMarks[6];

	lateralAnkle2 = landMarks[9];

	medialAnkle2 = landMarks[8];

	bool result;
	double error;

	//regis->RegistrationLandmarks(tibiaKnee2, tubercle2, lateralAnkle2, medialAnkle2, error);

	std::cout << "Land mark error: " << error << std::endl;

	regis->MakeRegistration(bones_points, tubercle2, lateralAnkle2, medialAnkle2);

	std::cout << "**************** Tibia Error: " << regis->getError() << std::endl;

	std::cout << "******************* Transform: \n" << regis->getTransformMarkerToCt() << std::endl;

	std::cout << "**************************************************************" << std::endl;

	/*TibiaRegistration * regis2 = new TibiaRegistration(BoneVTK, tibiaKnee, tubercle, lateralAnkle, medialAnkle);

	regis2->MakeRegistration(bones_points, tibiaKnee2, tubercle2, lateralAnkle2, medialAnkle2, true);

	std::cout << "**************** Tibia Error2: " << regis2->getError() << std::endl;

	std::cout << "******************* Transform2: \n" << regis2->getTransformMarkerToCt() << std::endl;*/
}

void LS_Rotula_Registration()
{
	RegistrationImageType::Pointer imgOut;
	Test::readImage<RegistrationImageType>(dicom_knee_cap, imgOut);

	vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData(vtk_kneeCap_right);

	PointTypeITK hipCenter, kneeCenter, lateralEpi, medialEpi;

	hipCenter = Registration::makeItkPoint(57.2, -47.89, 1158.81);

	kneeCenter = Registration::makeItkPoint(23.86, -55.92, 771.39);

	lateralEpi = Registration::makeItkPoint(-11.02, -48.33, 782.79);

	medialEpi = Registration::makeItkPoint(67.73, -49.17, 777.39);

	std::vector<RegistrationImageType::Pointer> imgPoints;
	imgPoints.push_back(imgOut);

	KneeCapRegistration* regis = new KneeCapRegistration(BoneVTK, hipCenter, kneeCenter, lateralEpi, medialEpi);

	/*PointTypeITK upPoint = regis->GetTopPoint();
	PointTypeITK downPoint = regis->GetDownPoint();
	PointTypeITK latPoint = regis->GetLateralPoint();
	PointTypeITK medPoint = regis->GetMedialPoint();

	Point upPointCV(upPoint[0], upPoint[1], upPoint[2]);
	Point downPointCV(downPoint[0], downPoint[1], downPoint[2]);
	Point latPointCV(latPoint[0], latPoint[1], latPoint[2]);
	Point medPointCV(medPoint[0], medPoint[1], medPoint[2]);

	Plane planeTest;
	planeTest.init(upPointCV, latPointCV, medPointCV);

	std::cout << "Distance: " << planeTest.eval(downPointCV) << std::endl;

	std::cout << upPointCV << std::endl;
	std::cout << downPointCV << std::endl;
	std::cout << latPointCV << std::endl;
	std::cout << medPointCV << std::endl;*/

	/*vtkSmartPointer<vtkPolyData> esfera1 = TestVTK::CreateSphereTest(upPointCV);
	vtkSmartPointer<vtkPolyData> esfera2 = TestVTK::CreateSphereTest(downPointCV);
	vtkSmartPointer<vtkPolyData> esfera3 = TestVTK::CreateSphereTest(latPointCV);
	vtkSmartPointer<vtkPolyData> esfera4 = TestVTK::CreateSphereTest(medPointCV);*/

	vtkNew<vtkAppendPolyData> filter;

	std::vector<PointTypeITK> pointList = regis->GetRegistrationPoints();

	std::vector<PointTypeITK> realPoints;

	cv::Mat rot = Test::GenerateRotationMatrix();
	cv::Mat trans = Test::GenerateTranslationMatrix();

	for (int i = 0; i < pointList.size(); i++)
	{
		cv::Mat pointMat = Test::ITKToMat(pointList[i]);
		cv::Mat newPoint = rot * pointMat + trans;

		PointTypeITK myITKPointTemp = Test::MatToITK(newPoint);
		PointTypeITK myITKPoint = Test::AddNoiseITK(myITKPointTemp, 2.0);

		realPoints.push_back(myITKPoint);
	}

	regis->MakeRegistration(realPoints);

	std::cout << "Error: " << regis->getError() << std::endl;

	//return;

	for (int i = 0; i < pointList.size(); i++)
	{
		Point cvPoint(pointList[i][0], pointList[i][1], pointList[i][2]);
		vtkSmartPointer<vtkPolyData> esfera = TestVTK::CreateSphereTest(cvPoint);
		filter->AddInputData(esfera);
	}

	/*filter->AddInputData(esfera1);
	filter->AddInputData(esfera2);
	filter->AddInputData(esfera3);
	filter->AddInputData(esfera4);*/

	filter->AddInputData(BoneVTK);
	filter->Update();


	///////////////////

	/*vtkNew<vtkPlane> vtkPlaneA;
	vtkPlaneA->SetOrigin(18.0132, -77.5969, 744.714);
	vtkPlaneA->SetNormal(-0.1836, -0.9735, -0.1363);
	vtkNew<vtkPlaneCollection> planes;
	planes->AddItem(vtkPlaneA);

	vtkNew<vtkClipClosedSurface> femurClipper;
	femurClipper->SetInputData(BoneVTK);
	femurClipper->SetClippingPlanes(planes);
	femurClipper->Update();

	auto cutRotula = femurClipper->GetOutput();*/

	/////////////

	SliceBorder::ShowPolyData(filter->GetOutput(), BoneVTK);
}

void ChangeCoordenate()
{
	std::string path = "D:\\3D_DICOM\\Left_Leg\\Knee_cap_good_seg_bad_coord.nrrd";
	std::string pathTest = "D:\\3D_DICOM\\Left_Leg\\new\\Femur_good_seg.nrrd";

	ImageType::Pointer imageIn, imageIn2;

	Test::readImage<ImageType>(path, imageIn);
	Test::readImage<ImageType>(pathTest, imageIn2);

	PointTypeITK boneLat, boneMed, boneTop, boneDown, kneeCapLat, kneeCapMed, kneeCapTop, kneeCapDown;

	boneLat = Registration::makeItkPoint(11.05, -3.54, -229.46);
	boneMed = Registration::makeItkPoint(-14.68, 1.00, -229.46);
	boneTop = Registration::makeItkPoint(-4.35, 2.43, -222.46);
	boneDown = Registration::makeItkPoint(-2.79, 8.68, -238.46);

	kneeCapLat = Registration::makeItkPoint(-27.30, 22.8, 110.0);
	kneeCapMed = Registration::makeItkPoint(2.05, 17.50, 108.0);
	kneeCapTop = Registration::makeItkPoint(-9.97, 14.32, 113.85);
	kneeCapDown = Registration::makeItkPoint(-9.97, 15.38, 92.95);

	std::vector<PointTypeITK> source = { kneeCapLat, kneeCapMed, kneeCapDown, kneeCapTop };

	std::vector<PointTypeITK> target = { boneLat, boneMed, boneDown, boneTop };

	itk::Rigid3DTransform<double>::Pointer transformMine = Registration::GetTransformBetweenPoints(source, target);

	itk::Matrix< double, 3, 3 > rotation = transformMine->GetMatrix();
	itk::Vector< double, 3 > translate = transformMine->GetTranslation();

	//std::cout << transform << std::endl;

	itk::ChangeInformationImageFilter<ImageType>::Pointer filter = itk::ChangeInformationImageFilter<ImageType>::New();
	filter->SetInput(imageIn);

	const ImageType::DirectionType direction = imageIn->GetDirection();
	const ImageType::DirectionType newDirection = direction * transformMine->GetMatrix();
	filter->SetOutputDirection(newDirection);
	filter->ChangeDirectionOn();


	filter->Update();


	try
	{
		filter->UpdateOutputInformation();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cout << "Error filter: " << error << std::endl;
	}
	std::cout << "**************************************" << std::endl;
	ImageType::Pointer outputImage = filter->GetOutput();



	using FilterType = itk::ResampleImageFilter<ImageType, ImageType>;
	FilterType::Pointer filter2 = FilterType::New();
	using TransformType = itk::AffineTransform< double, 3 >;

	TransformType::Pointer transform = TransformType::New();
	filter2->SetTransform(transform);

	using InterpolatorType = itk::NearestNeighborInterpolateImageFunction<ImageType, double >;

	InterpolatorType::Pointer interpolator = InterpolatorType::New();
	filter2->SetInterpolator(interpolator);
	filter2->SetDefaultPixelValue(0);

	const ImageType::SpacingType & spacing = outputImage->GetSpacing();
	const ImageType::PointType & origin = outputImage->GetOrigin();
	ImageType::SizeType size = outputImage->GetLargestPossibleRegion().GetSize();

	filter2->SetOutputOrigin(origin);
	filter2->SetOutputSpacing(spacing);
	filter2->SetSize(size);

	filter2->SetOutputDirection(outputImage->GetDirection());

	TransformType::OutputVectorType translation;
	translation[0] = translate[0];  // X translation in millimeters
	translation[1] = translate[1];  // Y translation in millimeters
	translation[2] = translate[2];
	transform->Translate(translation);

	std::cout << translate[0] << ", " << translate[1] << ", " << translate[2] << std::endl;

	filter2->SetInput(outputImage);
	filter2->Update();

	using WriterType = itk::ImageFileWriter<ImageType>;
	WriterType::Pointer writer = WriterType::New();
	writer->SetFileName("ZimageOut.nrrd");
	writer->SetInput(filter2->GetOutput());
	try
	{
		writer->Update();
	}
	catch (itk::ExceptionObject & error)
	{
		std::cout << "Error: " << error << std::endl;
	}
}

void Test30PointsVTK()
{
	/*
	Knee myKnee = CreateKneeFromFile_Numbers("D:\\sovajo\\Errores\\Error3", KneeSideEnum::KLeft);

	double tError;
	std::vector<Point> pCheckPointsFemur, pCheckPointsTibia;
	FindRegistrationPoints myPointsObj(myKnee);

	std::vector<RegistrationPoints> femurPoints = myPointsObj.GetRegistrationPointsFemur(pCheckPointsFemur, tError);
	std::cout << "Error femur: " << tError << std::endl;
	std::vector<RegistrationPoints> tibiaPoints = myPointsObj.GetRegistrationPointsTibia(pCheckPointsTibia, tError);
	std::cout << "Error tibia: " << tError << std::endl;

	std::vector<cv::Point3d> myPointsFemur, myPointsTibia;

	for (int i = 0; i < femurPoints.size(); i++)
	{
		RegistrationPoints obj = femurPoints[i];

		std::vector<Point> temp = ITKVectorToCV(obj.points);

		for (int j = 0; j < temp.size(); j++)
		{
			myPointsFemur.push_back(temp[j]);
		}
	}

	for (int i = 0; i < tibiaPoints.size(); i++)
	{
		RegistrationPoints obj = tibiaPoints[i];

		std::vector<Point> temp = ITKVectorToCV(obj.points);

		for (int j = 0; j < temp.size(); j++)
		{
			myPointsTibia.push_back(temp[j]);
		}
	}

	myPointsTibia.push_back(myKnee.getLateralPlateau());
	myPointsTibia.push_back(myKnee.getMedialPlateau());

	//auto surface1 = TestVTK::MergePolyWithSphere(myKnee.GetFemurPoly(), myPointsFemur);
	//auto surface2 = TestVTK::MergePolyWithSphere(myKnee.GetTibiaPoly(), myPointsTibia);

	TestVTK::show(myKnee.GetFemurPoly(), PointToPoint3D(pCheckPointsFemur));
	TestVTK::show(myKnee.GetFemurPoly(), myPointsFemur);
	TestVTK::show(myKnee.GetTibiaPoly(), PointToPoint3D(pCheckPointsTibia));
	TestVTK::show(myKnee.GetTibiaPoly(), myPointsTibia);

	//writeListPoints(myPointsFemur, "Template_Official_Femur_Left.json");

	//SliceBorder::ShowPolyData(myPointsObj.GetMedialSide(), myPointsObj.GetLateralSide());
	*/
}

void TestSlice()
{
	SliceBorder obj;

	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly;

	femurPoly = TestVTK::ReadPolyData("D:\\fei_error\\femur_fei.vtk");
	tibiaPoly = TestVTK::ReadPolyData("D:\\fei_error\\tibia_fei.vtk");

	SPlane ObliqueTibiaFei;
	ObliqueTibiaFei.init(cv::Point3d(0.78, -0.11, 0.61), cv::Point3d(32.9, 51.6, -289.1));
	cv::Point3d PointTibiaFei(-24.10, 64.27, -291.27);

	std::vector<vtkSmartPointer<vtkPolyData>> slices = obj.MakeSlicesFromPolyData(TestVTK::CleanPoly(femurPoly), 1, 4.0);
	bool result;
	std::string msg;

	vtkNew<vtkAppendPolyData> appendFilter;

	for (int i = 0; i < slices.size(); i++)
	{
		/*std::vector<PointTypeITK> mainPoints = obj.GetMainPoints(slices[i], 179);
		vtkNew<vtkPoints> myPoints;
		for (int j = 0; j < mainPoints.size(); j++)
		{
			double pnt[3] = { mainPoints[j][0], mainPoints[j][1] , mainPoints[j][2] };
			myPoints->InsertNextPoint(pnt);
		}

		auto spheres = TestVTK::CreatePolyLine(myPoints);
		SliceBorder::ShowPolyData(slices[i], spheres);*/

		appendFilter->AddInputData(slices[i]);
	}

	appendFilter->Update();

	auto surface = obj.GetSurface3D(slices, result, msg);

	auto smoothSurface = obj.SmoothSurface(surface);

	SliceBorder::ShowPolyData(appendFilter->GetOutput(), smoothSurface);

}

void getPelvisPoint()
{
	/*
	//vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData("D:\\3D_DICOM\\Pelvis\\pelvis_right_vtk.vtk");

	vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData("D:\\Mega_Trabajo\\Pelvis\\Segmentation.vtk");

	PointTypeITK anteriorAcetabulum, posteriorAcetabulum, superiorAcetabulum;

	////////////////////////////////Left
	anteriorAcetabulum = Registration::makeItkPoint(77.27, -7.34, 309.29);
	posteriorAcetabulum = Registration::makeItkPoint(77.27, 20.66, 303.29);
	superiorAcetabulum = Registration::makeItkPoint(105.27, 10.66, 329.29);

	//anteriorAcetabulum = Registration::makeItkPoint(-42.48, -3.12, 298.29);
	//posteriorAcetabulum = Registration::makeItkPoint(-50.48, 19.88, 299.29);
	//superiorAcetabulum = Registration::makeItkPoint(-73.48, -1.12, 334.29);

	std::vector<cv::Point3d> landMark = { Point(-42.48, -3.12, 298.29), Point(-50.48, 19.88, 299.29), Point(-73.48, -1.12, 334.29) };
	TestVTK::show(BoneVTK, landMark);


	PelvisRegistration* regis;

	try
	{
		regis = new PelvisRegistration(BoneVTK, anteriorAcetabulum, posteriorAcetabulum, superiorAcetabulum, RegisterSide::LEFT);
	}
	catch (const std::exception& e)
	{
		std::cout << e.what() << std::endl;
		return;
	}

	double error;

	std::vector<PointTypeITK> verificationPointsITK;
	std::vector<RegistrationPointsHip> points = regis->getRegistrationPointPelvis(verificationPointsITK, error);
	std::vector <cv::Point3d > pointsList, verificationPoints;

	for (int i = 0; i < points.size(); i++)
	{
		RegistrationPointsHip obj = points[i];

		std::vector<Point> temp = ITKVectorToCV(obj.points);

		for (int j = 0; j < temp.size(); j++)
		{
			pointsList.push_back(temp[j]);
		}
	}

	std::cout << "Error: " << error << std::endl;

	//writeListPoints(finalPoint, "Template_Pelvis_Left.json");

	auto verificationPointTemp = ITKVectorToCV(verificationPointsITK);
	for (int j = 0; j < verificationPointTemp.size(); j++)
	{
		verificationPoints.push_back(verificationPointTemp[j]);
	}

	TestVTK::show(BoneVTK, pointsList);
	TestVTK::show(BoneVTK, verificationPoints);
	*/
}

void RegistrationScale()
{
	//Knee kneeRight = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Right");
	//Knee kneeLeft = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Left");

	std::string folder = "D:\\sovajo\\Cases_Plan_PKA\\UKA-data"; // right
	UKA::IMPLANTS::Knee myKnee = CreateKneeFromFile_NumbersPKA(folder, UKA::IMPLANTS::KRight, UKA::IMPLANTS::KMedial);

	double error;

	UKA::IMPLANTS::FindRegistrationPoints myPointsObj(myKnee);
	std::vector<UKA::IMPLANTS::Point> pCheckPoints;
	std::vector<UKA::IMPLANTS::RegistrationPoints> points = myPointsObj.GetRegistrationPointsTibia(pCheckPoints, error);
	std::vector < cv::Point3d > pointsList;

	for (int i = 0; i < points.size(); i++)
	{
		UKA::IMPLANTS::RegistrationPoints obj = points[i];

		std::vector<Point> temp = ITKVectorToCV(obj.points);

		for (int j = 0; j < temp.size(); j++)
		{
			pointsList.push_back(temp[j]);
		}
	}

	/*cv::Mat mData(6, 1, CV_64F);
	mData.at<double>(0, 0) = 5;
	mData.at<double>(1, 0) = -5;
	mData.at<double>(2, 0) = 8;

	mData.at<double>(3, 0) = 0;
	mData.at<double>(4, 0) = 0;
	mData.at<double>(5, 0) = 0;

	LeastSquaresICP leastSquare(pointsList);

	error = leastSquare.LeastSquaresTest(myKnee.GetTibiaPoly(), mData);

	std::cout << "Error: " << error << std::endl;*/

	TestVTK::show(myKnee.GetTibiaPoly(), pointsList);

	/*auto surface = TestVTK::MergePolyWithSphere(myKnee.GetFemurPoly(), pointsList);

	SliceBorder::ShowPolyData(myKnee.GetFemurPoly(), surface);*/
}

void RegistrationScale2()
{

	Knee kneeRight = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Modo\\Right_Modo");
	Knee kneeLeft = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Modo\\Left_Modo");

	double error;

	Knee kneeTemplate = kneeLeft;

	Point diff = kneeTemplate.getMedialEpicondyle() - kneeTemplate.getLateralEpicondyle();

	std::cout << "AP: " << kneeTemplate.getFemurDirectVectorAP() << std::endl;
	std::cout << "TEA: " << kneeTemplate.getFemurVectorLateralTEA() << std::endl;
	std::cout << "Center: " << (kneeTemplate.getLateralEpicondyle() + kneeTemplate.getMedialEpicondyle()) / 2.0 << std::endl;
	std::cout << "Size: " << sqrt(diff.dot(diff)) << std::endl;

	std::vector<cv::Point3d> myPoints = { kneeTemplate.getLateralEpicondyle(), kneeTemplate.getMedialEpicondyle() };
	TestVTK::show(kneeTemplate.GetFemurPoly(), myPoints);

	/*FindRegistrationPoints tempObj(kneeTemplate);

	tempObj.RegisterFemurTest(tempObj.GetRegistrationPointsFemurTemplateLikeCV(true), kneeTemplate.getFemurVectorLateralTEA(),
		kneeTemplate.getFemurDirectVectorAP(), (kneeTemplate.getLateralEpicondyle() + kneeTemplate.getMedialEpicondyle()) / 2.0, 1.0);*/

		//writeListPoints(tempObj.GetRegistrationPointsFemurTemplateLikeCV(true), "Template_Femur_Left.json");

		/*auto surface = TestVTK::MergePolyWithSphere(kneeTemplate.GetFemurPoly(), tempObj.GetRegistrationPointsFemurTemplateLikeCV(true));

		SliceBorder::ShowPolyData(kneeTemplate.GetFemurPoly(), surface);*/
}

void TestPointsOnBone()
{
	vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly;

	femurPoly = TestVTK::ReadPolyData("D:\\fei_error\\femur_fei.vtk");
	tibiaPoly = TestVTK::ReadPolyData("D:\\fei_error\\tibia_fei.vtk");
	/*tibiaPoly = TestVTK::itkImageToSurface3D(tibia_right_modo_nrrd);
	femurPoly = TestVTK::itkImageToSurface3D(femur_right_modo_nrrd);*/

	std::vector<cv::Point3d> myPointsLeft = { Point(31.4973, 28.8212, -256.417), Point(33.473, 34.6371, -265.722), Point(34.1499, 41.6148, -271.478),
											Point(36.5492, 49.4876, -274.569), Point(39.10, 57.7607, -273.30), Point(39.10, 65.427, -272.496),

											Point(21.06, 23.28, -256.41), Point(24.06, 26.28, -267.41), Point(25.06, 35.28, -276.11),
											Point(27.06, 45.28, -280.11), Point(30.06, 56.28, -282.41), Point(31.06, 67.28, -281.31),

											Point(8.06, 32.78, -264.41), Point(12.06, 30.28, -261.41), Point(16.06, 26.28, -259.41),

											Point(7.93961, 29.7385, -250.75), Point(7.64742, 31.5059, -244.439), Point(7.39682, 33.2148, -237.953),

											Point(0.747355, 28.0879, -259.398), Point(-4.94, 29.78, -269.41), Point(-11.94, 37.28, -275.41),
											Point(-16.181, 46.3011, -280.567), Point(-17.2306, 55.7779, -283.318), Point(-19.0695, 65.1437, -282.218),

											Point(-11.64, 30.28, -259.52), Point(-17.94, 36.28, -265.41), Point(-20.8727, 41.3887, -270.579),
											Point(-24.169, 47.8868, -273.372), Point(-24.6325, 54.7879, -275.694), Point(-25.3161, 61.7255, -275.836) };

	std::vector<cv::Point3d> myPointsRight = { Point(1.10369, -76.6342, 784.868), Point(-0.609375, -73.3192, 776.212), Point(-2.52, -66.11, 771.70),
											   Point(-5.11817, -56.7823, 769.054), Point(-6.90556, -49.903, 766.008), Point(-8.20312, -42.7279, 767.164),

											   Point(9.48, -80.95, 784.78), Point(8.48, -78.41, 774.70), Point(8.48, -71.00, 766.70),
											   Point(5.48, -60.11, 763.00), Point(0.48, -51.11, 759.90), Point(-2.52, -40.11, 759.50),

											   Point(24.52, -71.2057, 773.8), Point(20.10, -75.00, 777.121), Point(15.5546, -78.5255, 780.304),

											   Point(25.7886, -76.7344, 787.879), Point(26.1991, -76.7396, 792.989), Point(26.6035, -76.1719, 798.068),

											   Point(31.48, -77.11, 780.25), Point(36.00, -75.11, 768.70), Point(43.48, -67.11, 761.20),
											   Point(49.48, -60.11, 755.70), Point(52.48, -50.11, 752.20), Point(52.48, -38.11, 752.30),

											   Point(43.0236, -76.7344, 779.76), Point(51.7267, -70.4775, 770.702), Point(56.4092, -63.872, 764.11),
											   Point(59.5057, -55.6463, 759.846), Point(60.6198, -47.4987, 758.022), Point(61.48, -39.11, 759.70) };


	std::vector<cv::Point3d> myPointsTibiaLeft = { Point(27.92, 49.98, -286.42), Point(25.42, 46.98, -285.42), Point(21.42, 45.98, -285.42), Point(17.42, 47.98, -286.42),
															Point(-4.08, 44.98, -286.42), Point(-6.58, 43.98, -286.42), Point(-9.08, 42.98, -286.42), Point(-11.58, 41.98, -287.42),
															Point(-15.58, 49.98, -286.42), Point(-19.58, 51.98, -287.42), Point(-22.58, 55.98, -289.42), Point(-24.58, 60.98, -290.42),
															Point(15.42, 34.98, -292.42), Point(13.42, 33.98, -295.42), Point(11.42, 32.98, -298.42), Point(9.42, 30.98, -301.42)/*, Point(7.42, 29.98, -304.42)*/,
															Point(-4.58, 35.98, -291.42), Point(-2.58, 35.98, -294.42), Point(-0.58, 33.98, -297.42), Point(1.42, 31.98, -300.42)/*, Point(3.42, 30.98, -303.42)*/,
															Point(-15.58, 41.98, -299.42), Point(-18.58, 44.98, -299.42), Point(-20.58, 47.98, -299.42), Point(-22.58, 50.98, -299.42)/*, Point(-24.58, 53.98, -299.42)*/,
															Point(5.42, 26.98, -308.42), Point(5.42, 25.98, -313.42), Point(5.42, 26.98, -318.42) };


	//auto surface = TestVTK::MergePolyWithSphere(femurPoly, myPointsLeft);

	TestVTK::show(TestVTK::CleanPoly(femurPoly), myPointsLeft);
}

//#include <pcl/recognition/ransac_based/trimmed_icp.h>
/*
void GetRaymondError()
{
	std::vector<Point> source = { Point(43.63,-43.26,59.36), Point(46.63,-37.26,-27.14), Point(-0.87,-38.26,58.36),
		Point(1.63,-38.16,-28.14), Point(-45.87,-39.76,58.11), Point(-43.37,-37.76,-28.89) };

	std::vector<Point> target = { Point(307.02, -323.32, -1860.49), Point(307.63, -238.79, -1891.96), Point(315.44, -307.97, -1818.28),
		Point(310.52, -223.58, -1849.89), Point(318.30, -293.04, -1776.88), Point(315.09, -208.23, -1807.85) };

	Eigen::Matrix4d transform;

	pcl::PointCloud<pcl::PointXYZ> cloud_in;
	pcl::PointCloud<pcl::PointXYZ> cloud_out;
	pcl::registration::TransformationEstimationSVD<pcl::PointXYZ, pcl::PointXYZ, double> TESVD;

	int tSize = source.size();

	for (int i = 0; i < tSize; i++)
	{
		cloud_in.points.push_back(pcl::PointXYZ(source[i].x, source[i].y, source[i].z));
		cloud_out.points.push_back(pcl::PointXYZ(target[i].x, target[i].y, target[i].z));
	}

	TESVD.estimateRigidTransformation(cloud_in, cloud_out, transform);

	double error = 0;

	for (int i = 0; i < tSize; i++)
	{
		Eigen::Vector4d myVector;
		myVector << source[i].x, source[i].y, source[i].z, 1.0;

		Eigen::Vector4d result = transform * myVector;

		Eigen::Vector4d myVector2;
		myVector2 << target[i].x, target[i].y, target[i].z, 1.0;

		Eigen::Vector4d diff = result - myVector2;
		double currentError = sqrt(diff.dot(diff));
		error = error + currentError;
		std::string msg = "Point " + std::to_string(i) + " error: ";
		std::cout << msg << currentError << std::endl;
	}

	std::cout << "Total error: " << error << std::endl;
}
*/

void GetHipPoints()
{
	vtkSmartPointer<vtkPolyData> hipRight, hipLeft;

	hipLeft = TestVTK::ReadPolyData("D:\\3D_DICOM_DATA\\Modo\\Left_Modo\\hip_left.vtk");
	hipRight = TestVTK::ReadPolyData("D:\\3D_DICOM_DATA\\Modo\\Right_Modo\\hip_right.vtk");

	Point hipTopRight(56.55, -48.89, 1180.50);
	Point hipTopLeft(-0.28, 57.54, 155.63);

	Plane obliqueRight, obliqueLeft;
	obliqueRight.init(Point(20.55, -48.89, 1153.50), Point(45.55, -54.89, 1132.50), Point(46.55, -33.89, 1132.50));
	obliqueLeft.init(Point(31.72, 56.54, 132.63), Point(8.72, 46.54, 110.63), Point(6.72, 63.54, 110.63));

	std::vector<cv::Point3d> myPoints;

	///////////////////////////////// Right
	/*HipPoints myHip(hipRight, hipTopRight, obliqueRight);
	myHip.GetPointOnTop(Point(20.55, -48.89, 1153.50), Point(56.55, -48.89, 1136.50), myPoints);
	myHip.GetPointOnTop(Point(46.55, -29.89, 1148.50), Point(43.55, -31.89, 1132.50), myPoints);
	myHip.GetTransversal(Point(43.55, -24.89, 1097.50), Point(30.55, -53.89, 1097.50), Point(14.55, -22.89, 1097.50), myPoints);
	myHip.GetSagitalDown(Point(37.55, -26.89, 1110.50), Point(42.55, -22.89, 1098.50), Point(33.55, -29.89, 1082.50), myPoints);*/

	/////////////////////////////// Lefta
	/*HipPoints myHip(hipLeft, hipTopLeft, obliqueLeft);
	myHip.GetPointOnTop(Point(29.72, 55.54, 132.63), Point(2.72, 54.54, 110.63), myPoints);
	myHip.GetPointOnTop(Point(14.72, 66.54, 111.63), Point(12.72, 45.54, 111.63), myPoints);
	myHip.GetTransversal(Point(9.72, 78.54, 86.63), Point(10.72, 60.54, 86.63), Point(24.72, 75.54, 86.63), myPoints);
	myHip.GetSagitalDown(Point(16.72, 70.54, 96.63), Point(9.72, 78.54, 85.63), Point(19.72, 68.54, 69.63), myPoints);*/

	////////////////////////////////// landMark Right

	/*myPoints.push_back(Point(31.55, -60.89, 1148.50));
	myPoints.push_back(Point(29.55, -58.89, 1127.50));
	myPoints.push_back(Point(-7.55, -42.89, 1137.50));*/

	//////////////////////////////////// landmark left

	/*myPoints.push_back(Point(23.72, 46.54, 127.63));
	myPoints.push_back(Point(20.72, 43.54, 98.63));
	myPoints.push_back(Point(61.72, 55.54, 112.63));*/


	TestVTK::show(hipLeft, myPoints);
}

void HipFemoralRegistration()
{
	/*
	vtkSmartPointer<vtkPolyData> hipRight, hipLeft;

	hipLeft = TestVTK::ReadPolyData("D:\\Mega_Trabajo\\Modo\\Left_Modo\\hip_left.vtk");
	hipRight = TestVTK::ReadPolyData("D:\\Mega_Trabajo\\Modo\\Right_Modo\\hip_right.vtk");

	PointTypeITK anteriorFemoralNeck, anteriorDistalTrochanter, lateralTrochanter;

	////////////////////////////////Left
	anteriorFemoralNeck = Registration::makeItkPoint(23.72, 47.1258, 127.53);
	anteriorDistalTrochanter = Registration::makeItkPoint(18.19, 44.49, 92.25);
	lateralTrochanter = Registration::makeItkPoint(61.9828, 55.54, 112.63);

	///////////////////////////////////Right
	//anteriorFemoralNeck = Registration::makeItkPoint(30.7806, -60.6593, 1149.39);
	//anteriorDistalTrochanter = Registration::makeItkPoint(33.015, -53.22, 1108.26);
	//lateralTrochanter = Registration::makeItkPoint(-7.67188, -42.89, 1137.5);

	HipRegistrationFemur regis(hipLeft, anteriorFemoralNeck, anteriorDistalTrochanter, lateralTrochanter, RegisterSide::LEFT);

	double error;
	std::vector<PointTypeITK> verificationPointsTemp;
	std::vector<RegistrationPointsHip> points = regis.getRegistrationPointAnterior(verificationPointsTemp, error);

	std::cout << "Error: " << error << std::endl;

	std::vector < cv::Point3d > pointsList, verificationPoint;

	for (int i = 0; i < points.size(); i++)
	{
		RegistrationPointsHip obj = points[i];

		std::vector<Point> temp = ITKVectorToCV(obj.points);

		for (int j = 0; j < temp.size(); j++)
		{
			pointsList.push_back(temp[j]);
		}
	}

	std::vector<Point> temp = ITKVectorToCV(verificationPointsTemp);

	for (int j = 0; j < temp.size(); j++)
	{
		verificationPoint.push_back(temp[j]);
	}

	std::vector < cv::Point3d > pp = {cv::Point3d(23.72, 47.1258, 127.53), cv::Point3d(18.19, 44.49, 92.25) , cv::Point3d(61.9828, 55.54, 112.63) };
	TestVTK::show(hipLeft, pointsList);
	TestVTK::show(hipLeft, verificationPoint);
	*/
}

int getPivot(std::vector<int>& data, int initPos, int endPos)
{
	if (endPos - initPos == 1)
	{
		return initPos;
	}

	int half = (initPos + endPos) / 2;

	int a = data[initPos];
	int b = data[half];
	int c = data[endPos - 1];
	int pos = initPos;
	int item;

	if (a >= b && b >= c)
	{
		item = b;
	}
	else if (b >= a && a >= c)
	{
		item = a;
	}
	else
	{
		item = c;
	}

	int pivot = -1;

	for (int i = initPos; i < endPos; i++)
	{
		if (data[i] <= item)
		{
			if (data[i] == item)
			{
				pivot = pos;
			}

			int temp = data[i];
			data[i] = data[pos];
			data[pos] = temp;

			pos++;
		}
	}

	data[pivot] = data[pos - 1];
	data[pos - 1] = item;

	return pos - 1;
}

void QuickSort(std::vector<int>& data, int initPos, int endPos)
{
	if (initPos < endPos)
	{
		int pivot = getPivot(data, initPos, endPos);
		QuickSort(data, initPos, pivot);
		QuickSort(data, pivot + 1, endPos);
	}
}

void QuickSort(std::vector<int>& data, int initPos, int endPos, int size)
{
	if (initPos < endPos)
	{
		int pivot = getPivot(data, initPos, endPos);

		if (size % 2 != 0 && pivot == (size / 2))
		{
			return;
		}

		if (pivot >= (size / 2))
		{
			QuickSort(data, initPos, pivot, size);
		}
		else
		{
			QuickSort(data, pivot + 1, endPos, size);
		}
	}
}

int readExample()
{
	int c;
	std::vector<int> media;
	std::cin >> c;

	while (c != 0)
	{
		std::cin.ignore();
		std::string S;
		std::getline(std::cin, S);
		std::istringstream T(S);
		int i;
		std::vector<int> V;
		while (T >> i)
		{
			V.push_back(i);
		}
		QuickSort(V, 0, c, c);
		if (c % 2 == 0)
		{
			media.push_back(V[c / 2] + V[((c / 2) - 1)]);
		}
		else
		{
			media.push_back(2 * (V[c / 2]));
		}
		std::cin >> c;
	}

	for (auto& elem : media)
	{
		std::cout << elem << std::endl;
	}
	return 0;
}

void PelvisImplantMatch()
{

	std::string pelvisPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\pelvis_full.vtk";
	std::string pelvisRealPath = "D:\\3D_DICOM_DATA\\Person_2\\Pelvis.vtk";
	//std::string implantFullPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\pelvis_implant.vtk";
	//std::string implantCupPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\acetabular_cup.vtk";
	//std::string implantStemPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\femoral_stem.vtk";

	auto pelvis3D_Real = TestVTK::ReadPolyData("D:\\Mega_Trabajo\\Person_2\\Pelvis.vtk");
	auto implant_cup_3d = TestVTK::ReadPolyDataSTL("D:\\Mega_Trabajo\\Person_2\\Implants\\acetabular_shell.stl");
	auto implant_stem_3d = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Implants\\THA\\femoral_stem.stl");
	auto implant_head_3d = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Implants\\THA\\femoral_head.stl");

	//auto pelvis3D_modo = TestVTK::ReadPolyData(pelvisPath);
	//auto pelvis3D_Real = TestVTK::ReadPolyData(pelvisRealPath);
	//auto implant3D = TestVTK::ReadPolyData(implantCupPath);
	//auto stem3D = TestVTK::ReadPolyData(implantStemPath);

	//////////////////// Modo

	//THA::IMPLANTS::Point centerOfRotation = Point(-63.83, 3.24, 305.26);
	//THA::IMPLANTS::Point pLeftASIS = Point(134.261, -43.3698, 338.323);
	//THA::IMPLANTS::Point pRightASIS = Point(-101.607, -51.7796, 338.064);
	//THA::IMPLANTS::Point pLeftPubicTubercle = Point(45.966, -30.0529, 279.477);
	//THA::IMPLANTS::Point pRightPubicTubercle = Point(-7.69984, -31.1355, 275.443);
	//THA::IMPLANTS::Point pLeftLesserTrochanter = Point(56.4254, 55.5531, 275.961); // not a good point (is not on the femur)
	//THA::IMPLANTS::Point pRightLesserTrochanter = Point(-30.7444, 53.9456, 272.75);// not a good point (is not on the femur)


	//////////////////Person_2
	THA::IMPLANTS::Point centerOfRotation = THA::IMPLANTS::Point(-72.30, 31.62, 1695.36);
	THA::IMPLANTS::Point pLeftASIS = THA::IMPLANTS::Point(122.481, -28.3346, 1758.14);
	THA::IMPLANTS::Point pRightASIS = THA::IMPLANTS::Point(-102.948, -23.1456, 1748.88);
	THA::IMPLANTS::Point pLeftPubicTubercle = THA::IMPLANTS::Point(47.1252, -19.7216, 1680.3);
	THA::IMPLANTS::Point pRightPubicTubercle = THA::IMPLANTS::Point(-11.6714, -18.741, 1682.92);
	THA::IMPLANTS::Point pLeftLesserTrochanter = THA::IMPLANTS::Point(103.508, 40.6431, 1650.47); // not a good point (is not on the femur)
	THA::IMPLANTS::Point pRightLesserTrochanter = THA::IMPLANTS::Point(-65.9917, 43.6621, 1643.16);// not a good point (is not on the femur)

	////////////////////////////////////////////////////////
	THA::IMPLANTS::Point implantTop = THA::IMPLANTS::Point(-32.16, 12.43, -1.95);
	THA::IMPLANTS::Point implantCup1 = THA::IMPLANTS::Point(-12.7411, -0.943873, -27.6937);
	THA::IMPLANTS::Point implantCup2 = THA::IMPLANTS::Point(-13.9138, -2.33964, 21.7003);
	THA::IMPLANTS::Point implantCup3 = THA::IMPLANTS::Point(2.47188, 17.1851, -0.863611);

	THA::IMPLANTS::Point stemTop = THA::IMPLANTS::Point(18.72, -22.79, -2.36);
	THA::IMPLANTS::Point stemBase = THA::IMPLANTS::Point(19.61, -162.36, -2.12);
	THA::IMPLANTS::Point stemHead = THA::IMPLANTS::Point(-12.13, -3.33, -2.35);

	THA::IMPLANTS::HipFemurOppside pFemurOppside;
	THA::IMPLANTS::HipFemur pFemur;
	pFemur.init(THA::IMPLANTS::Point(-66.68, 10.9341, 1705.03), THA::IMPLANTS::Point(-105.86, 33.7224, 1687.91), 
				THA::IMPLANTS::Point(-121.987, 54.1122, 1692.62), THA::IMPLANTS::Point(-77.6764, 42.9761, 1492.28),
				THA::IMPLANTS::Point(-65.2036, 46.8067, 1644.55),
				THA::IMPLANTS::Point(-67.9472, 29.0898, 1498.7), THA::IMPLANTS::Point(-93.7552, 37.0987, 1496.58), 
				THA::IMPLANTS::Point(-77.6764, 42.9761, 1492.28), pelvis3D_Real);

	pFemurOppside.init(THA::IMPLANTS::Point(117.477, 23.6717, 1723.36), THA::IMPLANTS::Point(156.129, 51.8456, 1698.98),
						THA::IMPLANTS::Point(102.612, 40.8368, 1653.34), THA::IMPLANTS::Point(115.105, 42.1657, 1491.97));

	THA::IMPLANTS::HipPelvis objPelvis;
	objPelvis.init(pLeftASIS, pRightASIS, pLeftPubicTubercle, pRightPubicTubercle, pelvis3D_Real, pFemur, pFemurOppside, THA::IMPLANTS::PelvisSide::RIGHT_SIDE, centerOfRotation, THA::IMPLANTS::Plane());

	THA::IMPLANTS::HipPelvisCupImplant objImplant;
	std::vector<THA::IMPLANTS::Point> pExternalPoints;
	objImplant.init(implantTop, implantCup1, implantCup2, implantCup3, pExternalPoints);

	THA::IMPLANTS::HipPelvisCupImplantMatch objMatch;
	objMatch.init(objPelvis, objImplant);

	itk::Rigid3DTransform<>::Pointer transformMatch = objMatch.getTransform(50, 30, 1, 1, 1);

	std::vector<THA::IMPLANTS::Point> stemHeadPlane = { THA::IMPLANTS::Point(-21.15, 3.66, -2.29),
														THA::IMPLANTS::Point(-20.5775, 4.32753, -5.33431),
														THA::IMPLANTS::Point(-21.8184, 2.84712, 0.554241),
														THA::IMPLANTS::Point(-19.1527, 6.02483, -0.912171),
														THA::IMPLANTS::Point(-23.4287, 0.92785, -3.51287) };

	THA::IMPLANTS::HipFemurStemImplant pImplantStem;
	pImplantStem.init(THA::IMPLANTS::Point(18.13, -22.78, -2.3), THA::IMPLANTS::Point(19.56, -162.34, -1.99), 
					THA::IMPLANTS::Point(-21.15, 3.66, -2.29), stemHeadPlane);

	std::vector<THA::IMPLANTS::Point> head_sphere = { THA::IMPLANTS::Point(-26.8174, 4.99074, -3.20225),
														THA::IMPLANTS::Point(-17.1012, 8.9965, 6.55039),
														THA::IMPLANTS::Point(-13.6094, 10.4183, -9.08246),
														THA::IMPLANTS::Point(-17.7142, 0.184736, -15.6298),
														THA::IMPLANTS::Point(-25.0233, -8.69219, -8.15378),
														THA::IMPLANTS::Point(-25.0151, -8.56904, 4.10177),
														THA::IMPLANTS::Point(-20.8054, -2.02465, 10.3861),
														THA::IMPLANTS::Point(-3.92857, 6.77136, 1.19525),
														THA::IMPLANTS::Point(-7.71749, -0.00787395, -14.2418),
														THA::IMPLANTS::Point(-18.6449, -14.4463, -6.11132) };
	THA::IMPLANTS::HipFemurStemHeadImplant pImplantHead;
	pImplantHead.init(THA::IMPLANTS::Point(-2.68735, -2.12955, 0.975769), THA::IMPLANTS::Point(-10.5758, -11.5306, 4.53979),
		THA::IMPLANTS::Point(-7.91941, -8.3649, -10.6989), THA::IMPLANTS::Point(-21.7549, 4.16821, -1.64382), head_sphere);
	
	THA::IMPLANTS::HipFemurStemImplantMatch stemMatch;
	stemMatch.init(objPelvis, pImplantStem);

	THA::IMPLANTS::HipFemurStemHeadImplantMatch stemHeadMatch;
	stemHeadMatch.init(pImplantHead, pImplantStem);

	auto pelvisInfo = THA::IMPLANTS::HipPelvisImplantsMatchInfo(objPelvis, objImplant, pImplantStem, pImplantHead,
																transformMatch, stemMatch.getStemTransform(), 
																stemHeadMatch.getStemHeadTransform());

	auto headToStemTransform = stemHeadMatch.getStemHeadTransform();
	auto transformImplantHeadToStemTemp = TestVTK::TransformPoly(implant_head_3d, headToStemTransform->GetMatrix(), headToStemTransform->GetTranslation());


	double antVersion = pelvisInfo.getCupAntversion();
	double inclination = pelvisInfo.getCupInclination();
	double shiftSuperior = pelvisInfo.getCupShiftSuperior();
	double shiftAnterior = pelvisInfo.getCupShiftAnterior();
	double shiftLateral = pelvisInfo.getCupShiftLateral();

	//std::cout << "Incination: " << inclination << ", Version: " << antVersion << std::endl;
	//std::cout << "Superior: " << shiftSuperior << ", Anterior: " << shiftAnterior << ", Lateral: " << shiftLateral << std::endl;

	//std::cout << "Offset: " << objPelvis.getCombinedOffsetDistance() << " Offset Opside: " << objPelvis.getCombinedOffsetDistanceOppsite() << std::endl;
	//std::cout << "Lenght: " << objPelvis.getHipLengthDistance() << " Offset Opside: " << objPelvis.getHipLengthDistanceOppsite() << std::endl;

	std::cout << "-------------------------------------------" << std::endl;
	/*
	double antVersionTemp = pelvisInfo.getCupAntversion();
	double inclinationTemp = pelvisInfo.getCupInclination();

	for (float i = 1; i < 10; i ++)
	{
		inclinationTemp += 5;
		antVersionTemp += 5;
		pelvisInfo.setCupAngles(inclinationTemp, antVersionTemp);

		antVersion = pelvisInfo.getCupAntversion();
		inclination = pelvisInfo.getCupInclination();
		shiftSuperior = pelvisInfo.getCupShiftSuperior();
		shiftAnterior = pelvisInfo.getCupShiftAnterior();
		shiftLateral = pelvisInfo.getCupShiftLateral();

		std::cout << "Incination: " << inclination << ", Version: " << antVersion << std::endl;
		std::cout << "Superior: " << shiftSuperior << ", Anterior: " << shiftAnterior << ", Lateral: " << shiftLateral << std::endl;

		std::cout << "-------------------------------------------" << std::endl;
	}
	*/
	
	std::cout << ": Current angle 1: " << pelvisInfo.getStemVersion() << std::endl;
	
	pelvisInfo.matchStemToCupRotationCenter();

	std::cout << ": Current angle 2: " << pelvisInfo.getStemVersion() << std::endl;

	for (float i = -20; i < 21; i += 5)
	{
		
		if (i == 20)
		{
			//std::cout << " ***************** Matching Stem with Hip Center: " <<std::endl;
			//pelvisInfo.matchStemToHipRotationCenter();
			//pelvisInfo.matchStemToCupRotationCenter();
		}
		else
		{
			//pelvisInfo.setStemHeadTranslation(i, i, i);
			//pelvisInfo.setStemVersionAngle(i);
		}

		//std::cout << "Setting stem version angle: " << i << ": Getting Current angle: " << pelvisInfo.getStemVersion() <<std::endl;

		//std::cout << "Hip Lenght: " << pelvisInfo.getHipLengthDistance() << ": Hip offset: " << pelvisInfo.getCombinedOffsetDistance()<< std::endl;
		//std::cout << "Set stem version: " << i << ": Current angle: " << pelvisInfo.getStemVersion() << std::endl;

		auto stemTransformMatch = pelvisInfo.getITKStemTransform();

		auto transformImplantCup = TestVTK::TransformPoly(implant_cup_3d, transformMatch->GetMatrix(), transformMatch->GetTranslation());
		auto transformImplantStem = TestVTK::TransformPoly(implant_stem_3d, stemTransformMatch->GetMatrix(), stemTransformMatch->GetTranslation());
		auto transformImplantHead = TestVTK::TransformPoly(transformImplantHeadToStemTemp, stemTransformMatch->GetMatrix(), stemTransformMatch->GetTranslation());

		//std::cout << "Stem superior distance: " << pelvisInfo.getStemAxisShiftSuperiorHip() << ", Stem Lateral Distance: " << pelvisInfo.getStemAxisShiftLateralHip() <<", Stem Anterior Distance: " << pelvisInfo.getStemAxisShiftAnteriorHip()<< std::endl;
		
		//pelvisInfo.setStemAxisShiftHip(10, 20, 30);

		//std::cout << "Stem superior distance: " << pelvisInfo.getStemAxisShiftSuperiorHip() << ", Stem Lateral Distance: " << pelvisInfo.getStemAxisShiftLateralHip() << ", Stem Anterior Distance: " << pelvisInfo.getStemAxisShiftAnteriorHip() << std::endl;


		//pelvisInfo.setStemTranslation(i, i, i);

		//std::cout << "Set stem version: " << i << " gest Stem version Angle: " << pelvisInfo.getStemVersion() << ", Stem Superior Distance: " << pelvisInfo.getStemAxisShiftSuperiorCup() << ", Stem Lateral Distance: " << pelvisInfo.getStemAxisShiftLateralCup() << ", Stem Anterior Distance: " << pelvisInfo.getStemAxisShiftAnteriorCup() << std::endl;

		double abductionCup = 40 + i / 20;
		double versionCup = 20 + i / 20;
		
		pelvisInfo.setPelvisTiltAngleDegree(i / 20 + 1);
		pelvisInfo.setCupAngles(abductionCup, versionCup);
		
		//pelvisInfo.setCupTranslation(0, 0, 10);

		
		std::cout << "Cup new abduction: " << abductionCup << ", Cup new version: " << versionCup << std::endl;
		std::cout << "Cup abduction: " << pelvisInfo.getCupInclination() << ", Cup version: " << pelvisInfo.getCupAntversion() << std::endl;
		std::cout << "Tilt angle: " << i/20 + 1 << std::endl;
		

		//std::cout << i << " Cup: " << pelvisInfo.getCupShiftAnterior() << ", Cup: " << pelvisInfo.getCupShiftLateral() << ", Cup: " << pelvisInfo.getCupShiftSuperior() <<std::endl;
		
		//std::cout << " Hip lenght: " << objPelvis.getHipLengthDistance() << ", new lenght: " << pelvisInfo.getHipLengthDistance() << ", Hip Offset: " << objPelvis.getCombinedOffsetDistance() << ", Hip New Offset: " << pelvisInfo.getCombinedOffsetDistance() << std::endl;

		
		//vtkNew<vtkAppendPolyData> appendFilter;
		//appendFilter->AddInputData(transformImplantCup);
		//appendFilter->AddInputData(transformImplantStem);
		//appendFilter->AddInputData(transformImplantHead);
		//appendFilter->AddInputData(pelvis3D_Real);
		//appendFilter->Update();
		std::vector<vtkSmartPointer<vtkPolyData>> polyList;
		polyList.push_back(transformImplantCup);

		TestVTK::show(pelvis3D_Real, polyList);

	}

	/*
	//pelvisInfo.matchStemToCupRotationCenter();
	std::cout << " ***************** Matching Stem with Hip Center: " << std::endl;
	std::cout << "Hip Lenght: " << pelvisInfo.getHipLengthDistance() << ": Hip offset: " << pelvisInfo.getCombinedOffsetDistance() << std::endl;

	pelvisInfo.setCupTranslation(10, 0, 0);
	std::cout << "Shift Cup Superior: " << 10 << std::endl;
	std::cout << "Hip Lenght: " << pelvisInfo.getHipLengthDistance() << ": Hip offset: " << pelvisInfo.getCombinedOffsetDistance() << std::endl;
	//pelvisInfo.matchStemToCupRotationCenter();
	pelvisInfo.setCupTranslation(0, 10, 0);
	std::cout << "Shift Cup Lateral: " << 10 << std::endl;
	std::cout << "Hip Lenght: " << pelvisInfo.getHipLengthDistance() << ": Hip offset: " << pelvisInfo.getCombinedOffsetDistance() << std::endl;
	//pelvisInfo.matchStemToCupRotationCenter();
	pelvisInfo.setCupTranslation(0, 0, 10);
	std::cout << "Shift Cup Anterior: " << 10 << std::endl;
	std::cout << "Hip Lenght: " << pelvisInfo.getHipLengthDistance() << ": Hip offset: " << pelvisInfo.getCombinedOffsetDistance() << std::endl;
	*/
	//std::vector<vtkSmartPointer<vtkPolyData>> polyList = { transformImplantCup, transformImplantStem };

	//std::vector<vtkSmartPointer<vtkPolyData>> polyList = { transformImplantCup };

	//cv::Point3d temp;
	//std::vector<cv::Point3d> axis;
	/* for (float i = 0; i < 100; i++)
	 {
		 temp = objPelvis.getPubicJoin() + i * objPelvis.getVectorAP();
		 axis.push_back(temp);
		 appendFilter->AddInputData(TestVTK::CreateSphereTest(temp));

		 temp = objPelvis.getPubicJoin() - i * objPelvis.getVectorAP();
		 axis.push_back(temp);
		 appendFilter->AddInputData(TestVTK::CreateSphereTest(temp));

		 temp = objPelvis.getPubicJoin() + i * objPelvis.getVectorASIS();
		 axis.push_back(temp);
		 appendFilter->AddInputData(TestVTK::CreateSphereTest(temp));

		 temp = objPelvis.getPubicJoin() - i * objPelvis.getVectorASIS();
		 axis.push_back(temp);
		 appendFilter->AddInputData(TestVTK::CreateSphereTest(temp));

		 temp = objPelvis.getPubicJoin() + i * objPelvis.getVectorInfSup();
		 axis.push_back(temp);
		 appendFilter->AddInputData(TestVTK::CreateSphereTest(temp));

		 temp = objPelvis.getPubicJoin() - i * objPelvis.getVectorInfSup();
		 axis.push_back(temp);
		 appendFilter->AddInputData(TestVTK::CreateSphereTest(temp));

	 }*/

	
	//TestVTK::show(objPelvis.getPelvisVTK(), polyList);

	//TestVTK::show(implant3D, polyList);

	/*auto startTime = std::chrono::system_clock::now();

	auto extractData = objPelvis.extractFemurAxisVector();

	auto endTime = std::chrono::system_clock::now();

	std::chrono::duration<double> elapsed_seconds = endTime - startTime;

	std::cout << "****elapsed time: " << elapsed_seconds.count() << std::endl;

	TestVTK::show(extractData);*/
}

void PelvisImplantMatch_2()
{
	/*
	std::string pelvisPath = "D:\\3D_DICOM_DATA\\Pelvis Fei\\pelvis_full.vtk";
	std::string implantFullPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\pelvis_implant.vtk";
	std::string implantCupPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\acetabular_cup.vtk";
	std::string implantStemPath = "D:\\3D_DICOM_DATA\\Modo_Pelvis\\femoral_stem.vtk";

	auto pelvis3D = TestVTK::ReadPolyData(pelvisPath);
	auto implant3D = TestVTK::ReadPolyData(implantCupPath);
	auto stem3D = TestVTK::ReadPolyData(implantStemPath);

	Point centerOfRotation = Point(-88.13, 12.56, 1579.16);

	Point pLeftASIS = Point(142.074, -24.3413, 1641.76);
	Point pRightASIS = Point(-127.583, -28.5375, 1640.24);
	Point pLeftPubicTubercle = Point(42.7281, -40.0887, 1549.96);
	Point pRightPubicTubercle = Point(-36.5689, -45.6939, 1559.42);
	Point pLeftLesserTrochanter = Point(96.4447, 35.576, 1514.4);
	Point pRightLesserTrochanter = Point(-88.0578, 31.6672, 1512.73);

	Point implantTop = Point(-32.16, 12.43, -1.95);
	Point implantCup1 = Point(-12.7411, -0.943873, -27.6937);
	Point implantCup2 = Point(-13.9138, -2.33964, 21.7003);
	Point implantCup3 = Point(2.47188, 17.1851, -0.863611);

	Point stemTop = Point(18.72, -22.79, -2.36);
	Point stemBase = Point(19.61, -162.36, -2.12);
	Point stemHead = Point(-12.13, -3.33, -2.35);

	HipPelvis objPelvis;
	objPelvis.init(pLeftASIS, pRightASIS, pLeftPubicTubercle, pRightPubicTubercle, pLeftLesserTrochanter, pRightLesserTrochanter, pelvis3D, PelvisSide::RIGHT_SIDE);

	HipPelvisCupImplant objImplant;
	objImplant.init(implantTop, implantCup1, implantCup2, implantCup3);

	HipPelvisCupImplantMatch objMatch;
	objMatch.init(objPelvis, objImplant, centerOfRotation);

	itk::Rigid3DTransform<>::Pointer transformMatch = objMatch.getTransform();

	HipFemurStemImplant objImplantStem;
	objImplantStem.init(stemTop, stemBase, stemHead);

	HipFemurStemImplantMatch objMatchStem;
	objMatchStem.init(objPelvis, objImplantStem, centerOfRotation);

	auto transformImplantCup = TestVTK::TransformPoly(implant3D, transformMatch->GetMatrix(), transformMatch->GetTranslation());

	auto transformImplantStem = TestVTK::TransformPoly(stem3D, objMatchStem.GetRotationMatrix(), objMatchStem.GetTranslationMatrix());

	std::vector<vtkSmartPointer<vtkPolyData>> polyList = { transformImplantCup, transformImplantStem };

	TestVTK::show(objPelvis.getPelvisVTK(), polyList);
	*/
}

void RotulaGroovePath()
{
	//Knee kneeRight = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Right");
	//Knee kneeLeft = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Left");
	//"D:\\sovajo\\Errores\\Error5", TKA::IMPLANTS::KRight
	//"D:\\sovajo\\Errores\\Error8", TKA::IMPLANTS::KLeft
	Knee newKnew = CreateKneeFromFile_Numbers("D:\\sovajo\\Errores\\Error3", TKA::IMPLANTS::KLeft);

	Knee myKnee = newKnew;
	myKnee.findFemurCondylesUsingDistalPoints();
	myKnee.makeKneeGroovePathNew();

	std::vector<Point> pointPath = myKnee.getKneeGroovePath();
	std::vector<Point> pointPathOutLiers; // = myKnee.getKneeGrooveOutLiers();
	std::vector<cv::Point3d> allPoints;

	for (int i = 0; i < pointPath.size(); i++)
	{
		allPoints.push_back(pointPath[i]);
	}

	for (int i = 0; i < pointPathOutLiers.size(); i++)
	{
		//allPoints.push_back(pointPathOutLiers[i]);
	}
	TestVTK::show(myKnee.GetFemurPoly(), myKnee.getFitGroovePointsTest());
	TestVTK::show(myKnee.GetFemurPoly(), myKnee.getFitVertexGroovePointsTest());
	TestVTK::show(myKnee.GetFemurPoly(), allPoints);
}

void TestHullPoints()
{
	std::string femurImplantStr = "D:\\fei_error\\problem2\\femur_implant";
	std::string tibiaImplantStr = "D:\\fei_error\\TibiaImplant_Problem\\tibia_implant";

	//Knee newKnew = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Modo\\Right_Modo");
	//Knee newKnew = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Right");

	Knee newKnew = CreateKneeFromFile_Numbers("D:\\fei_error\\problem2\\knee");

	Knee myKnee = newKnew;

	vtkSmartPointer<vtkPolyData> polyTibiaImplant;
	FemurImplant femurImplant = CreateFemurImplantFromFile(femurImplantStr);
	TibiaImplant tibiaImplant = CreateTibiaImplantFromFile(tibiaImplantStr, polyTibiaImplant);

	FemurImplantMatch femurImplantMatch;
	TibiaImplantMatch tibiaImplantMatch;

	femurImplantMatch.init(femurImplant, myKnee);
	tibiaImplantMatch.init(tibiaImplant, myKnee);

	itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer tibiaTransformIn = itk::VersorRigid3DTransform<double>::New();

	femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
	femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrixByCortex());

	tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
	tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

	itk::Rigid3DTransform<double>::Pointer femurTransformOut = itk::VersorRigid3DTransform<double>::New();

	std::vector<PointTypeITK> hullFemur, hullTibia;

	try
	{
		hullFemur = femurImplantMatch.GetHullPoints(femurTransformIn, femurTransformOut, FemurImplantMatch::kPlaneB);
		std::cout << "Hull femur size: " << hullFemur.size() << std::endl;
		std::cout << femurTransformOut << std::endl;
	}
	catch (const ImplantsException& e)
	{
		Test::myPrint("Error femur get hull");
		std::cout << e.what() << std::endl;
	}

	itk::Rigid3DTransform<double>::Pointer tibiaTransformOutNew = itk::VersorRigid3DTransform<double>::New();

	try
	{
		hullTibia = tibiaImplantMatch.GetHullPoints(tibiaTransformIn, tibiaTransformOutNew, 5, 10);
		std::cout << "Hull tibia new size: " << hullTibia.size() << std::endl;
		std::cout << tibiaTransformOutNew << std::endl;
	}
	catch (const ImplantsException& e)
	{
		Test::myPrint("Error tibia get hull new");
		std::cout << e.what() << std::endl;
	}

	std::vector<PointTypeITK> tFinalVector = hullFemur;

	vtkSmartPointer<vtkPolyData> polyNew = TestVTK::CreatePolyLine(tFinalVector);

	std::vector<vtkSmartPointer<vtkPolyData>> polyList;
	polyList.push_back(polyNew);

	std::vector<cv::Point3d> testPointLandMarkLat = { myKnee.getLateralCondyle(), myKnee.getLateralInferiorFemurPoint() };

	std::vector<cv::Point3d> testPointLandMarkMed = { myKnee.getMedialCondyle(), myKnee.getMedialInferiorFemurPoint() };

	std::vector<Point> testPoint = myKnee.getKneeGroovePath();
	std::vector<cv::Point3d> testPoint2;

	for (int i = 0; i < testPoint.size(); i++)
	{
		testPoint2.push_back(testPoint[i]);
	}

	TestVTK::show(myKnee.GetFemurPoly(), polyList);
}

void KneeimplantMatch()
{
	std::string femurImplantStr = "D:\\3D_DICOM_DATA\\Knee_Implant\\femur\\left";
	std::string tibiaImplantStr = "D:\\3D_DICOM_DATA\\Knee_Implant\\tibia";

	Knee kneeLeftModo = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Modo\\Left_Modo");
	Knee kneeLeft = CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Left");

	Knee myKnee = kneeLeft;

	vtkSmartPointer<vtkPolyData> polyTibiaImplant;

	FemurImplant femurImplant = CreateFemurImplantFromFile(femurImplantStr);
	TibiaImplant tibiaImplant = CreateTibiaImplantFromFile(tibiaImplantStr, polyTibiaImplant);

	FemurImplantMatch femurImplantMatch;
	TibiaImplantMatch tibiaImplantMatch;

	femurImplantMatch.init(femurImplant, myKnee);
	tibiaImplantMatch.init(tibiaImplant, myKnee);

	itk::Rigid3DTransform<double>::Pointer femurTransformIn = itk::VersorRigid3DTransform<double>::New();
	itk::Rigid3DTransform<double>::Pointer tibiaTransformIn = itk::VersorRigid3DTransform<double>::New();

	femurTransformIn->SetMatrix(femurImplantMatch.GetRotationMatrix());
	femurTransformIn->SetOffset(femurImplantMatch.GetTranslationMatrixByCortex());

	tibiaTransformIn->SetMatrix(tibiaImplantMatch.GetRotationMatrix());
	tibiaTransformIn->SetOffset(tibiaImplantMatch.GetTranslationMatrix());

	Balance balance;
	balance.init(myKnee, femurImplant, tibiaImplant, femurTransformIn, tibiaTransformIn);

	Plane planeE = femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneE);
	planeE.reverseByPoint(myKnee.getMedialCondyle());
	planeE.movePlaneOnNormal(-3.0);

	vtkSmartPointer<vtkPolyData> femurCut = SplitPoly(balance.getFemurPolyCut(), planeE, myKnee.getMedialCondyle());

	vtkSmartPointer<vtkPolyData> tibiaCut = SplitPoly(myKnee.GetTibiaPoly(), tibiaImplantMatch.getTibiaPlane(), myKnee.getAnkleCenter());

	auto transformImplantFemur = TestVTK::TransformPoly(femurImplant.GetImplantModel(), femurImplantMatch.GetRotationMatrix(), femurImplantMatch.GetTranslationMatrixByCortex());

	auto transformImplantTibia = TestVTK::TransformPoly(polyTibiaImplant, tibiaImplantMatch.GetRotationMatrix(), tibiaImplantMatch.GetTranslationMatrix());

	std::vector<vtkSmartPointer<vtkPolyData>> polyList = { transformImplantTibia };

	//TestVTK::show(tibiaCut, polyList);

	std::cout << "Plane C: " << femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneC).getNormalVector() << std::endl;
	std::cout << "Knee Axis: " << myKnee.getDirectVectorFemurAxis() << std::endl;

	std::cout << "Plane C: " << femurImplantMatch.GetPlane(FemurImplantMatch::kPlaneA).getNormalVector() << std::endl;
	std::cout << "Knee Axis: " << myKnee.getFemurDirectVectorAP() << std::endl;
}

std::vector<Point> getSpherePoints(std::string dir)
{
	std::ifstream infile(dir);
	std::vector<Point> points;
	double a, b, c;

	while (infile >> a >> b >> c)
	{
		Point myPoint;
		myPoint.x = a;
		myPoint.y = b;
		myPoint.z = c;
		points.push_back(myPoint);
	}

	return points;
}

void GenerateFilesJSON()
{
	std::string source, target;
	source = "D:\\Excel\\New_30_9_2022\\Raymond_Files\\Source";
	target = "D:\\Excel\\New_30_9_2022\\Raymond_Files\\Target";

	LoadFile myFiles("brand", "model", source, target);
	std::vector<std::string> errorList, overwriteWarning;
	myFiles.readFiles(errorList, overwriteWarning);

	for (int i = 0; i < errorList.size(); i++)
	{
		std::cout << errorList[i] << std::endl;
	}

	for (int i = 0; i < overwriteWarning.size(); i++)
	{
		std::cout << overwriteWarning[i] << std::endl;
	}

	errorList.clear();
	myFiles.writeFiles(errorList, true);

	for (int i = 0; i < errorList.size(); i++)
	{
		std::cout << errorList[i] << std::endl;
	}

	std::cout << "Finish!!!!" << std::endl;
}

enum kind_dataset { VALIDATION, TRAINING, TESTING };

void GenerateDataSet(kind_dataset dataset, int labels, bool gauss = false, int pInit = 0)
{
	std::string a = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\1_Bad\\";//******************762 begin bone
	std::string b = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\2_Good\\";
	std::string c = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\3_Good\\";
	std::string d = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\4_SoSo\\";//***************677
	std::string e = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\5_Bad\\";//*********************704 begin bone
	std::string f = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\6_Good\\";//0
	std::string g = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\7_SoSo_half\\";//*****************
	std::string h = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\8_SoSo_Half_good\\";//************** 0
	std::string m = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\9_SoSo_half\\";//*************
	std::string n = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\10_SoSo_Half_good\\";//************ 0

	std::vector<std::string> ct_train = { a, b, c, d, f, h, n };
	std::vector<int> slices_train = { 762, 719, 711, 667, 0, 0, 0 };

	std::vector<std::string> ct_val = { d };
	std::vector<int> slices_val = { 677 };

	std::vector<std::string> ct_test = { e };
	std::vector<int> slices_test = { 704 };

	std::string images_pathh, label_path;
	std::vector <std::string> ct;
	std::vector <int> slices_begin;

	if (dataset == kind_dataset::TRAINING)
	{
		images_pathh = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\Data\\Train\\images\\";
		label_path = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\Data\\Train\\labels\\";

		ct = ct_train;
		slices_begin = slices_train;
	}
	else if (dataset == kind_dataset::VALIDATION)
	{
		images_pathh = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\Data\\Val\\images\\";
		label_path = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\Data\\Val\\labels\\";

		ct = ct_val;
		slices_begin = slices_val;
	}
	else if (dataset == kind_dataset::TESTING)
	{
		images_pathh = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\Data\\Test\\images\\";
		label_path = "D:\\3D_DICOM_DATA\\pelvis_12\\xiazhi\\Data\\Test\\labels\\";

		ct = ct_test;
		slices_begin = slices_test;
	}

	int cont1 = pInit;
	int cont2 = pInit;

	for (int i = 0; i < ct.size(); i++)
	{
		ImageType::Pointer inputImg;

		if (labels == 3)
		{
			Test::readImage<ImageType>(ct[i] + "Seg.nrrd", inputImg);
		}
		else
		{
			Test::readImage<ImageType>(ct[i] + "Seg_New.nrrd", inputImg);
		}

		ImageType::Pointer inputImg2;

		if (gauss == true)
		{
			using FilterType = itk::SmoothingRecursiveGaussianImageFilter<ImageType, ImageType>;
			FilterType::Pointer smoothFilter = FilterType::New();

			ImageType::Pointer gaussImage;
			Test::readImage<ImageType>(ct[i] + "Raw_New.nrrd", gaussImage);

			smoothFilter->SetSigma(0.5);
			smoothFilter->SetInput(gaussImage);
			smoothFilter->Update();
			inputImg2 = smoothFilter->GetOutput();
		}
		else
		{
			Test::readImage<ImageType>(ct[i] + "Raw_New.nrrd", inputImg2);
		}

		std::vector<SegmentImageType::Pointer> images = Test::DICOM_TO_GrayScale_Slices_CV(inputImg2, false, slices_begin[i] - 1, labels);
		std::vector<SegmentImageType::Pointer> masks = Test::DICOM_TO_GrayScale_Slices_CV(inputImg, true, slices_begin[i] - 1, labels);

		for (int k = 0; k < images.size(); k++)
		{
			std::string var = images_pathh;
			std::string name = var + "image_" + std::to_string(cont1);
			Test::SaveImage<SegmentImageType>(images[k], name, false);
			cont1++;
		}

		for (int k = 0; k < masks.size(); k++)
		{
			std::string var = label_path;
			std::string name = var + "mask_" + std::to_string(cont2);
			Test::SaveImage<SegmentImageType>(masks[k], name, false);
			cont2++;
		}
		std::cout << ct[i] << std::endl;
	}

}

enum MyException
{
	except_a,
	except_b,
	except_c
};

void function_f() {
	throw except_b;
}

struct greater
{
	template<class T>
	bool operator()(T const &a, T const &b) const { return a > b; }
};

std::vector<Point> readPointsRegistrationTest(std::string path, std::vector<Point>& fill)
{
	std::ifstream infile(path);
	std::string line;

	double num;
	int cont = 0;
	std::vector<double> temp;
	std::vector<Point> points;
	while (infile.good())
	{
		infile >> num;
		cont++;
		temp.push_back(num);
		if (cont % 4 == 0)
		{
			points.push_back(Point(temp[1], temp[2], temp[3]));
			fill.push_back(Point(temp[1], temp[2], temp[3]));
			temp.clear();
		}
	}
	return points;
}

void getRandomPoints(std::vector<Point> source, int amount, std::vector<PointTypeITK>& target)
{
	int tSize = amount;
	if (amount > source.size())
	{
		tSize = source.size();
	}

	std::mt19937 generator(std::random_device{}());

	std::uniform_int_distribution<std::size_t> distribution(0, tSize - 1);

	for (std::size_t n = 0; n < amount; ++n)
	{
		std::size_t number = distribution(generator);
		Point temp = source[number];
		target.push_back(Registration::makeItkPoint(temp.x, temp.y, temp.z));
	}
}

void redPointsTestError()
{
	std::string path1 = "D:\\fei_error\\Registration_Error\\cross\\surface_verify_points_lateral_out.ini";
	std::string path2 = "D:\\fei_error\\Registration_Error\\cross\\surface_verify_points_lateral_condyle.ini";
	std::string path3 = "D:\\fei_error\\Registration_Error\\cross\\surface_verify_points_medial_condyle.ini";
	std::string path4 = "D:\\fei_error\\Registration_Error\\cross\\surface_verify_points_medial_out.ini";
	std::string path5 = "D:\\fei_error\\Registration_Error\\cross\\surface_verify_points_lateral_cross.ini";
	std::string path6 = "D:\\fei_error\\Registration_Error\\cross\\surface_verify_points_center_cortex.ini";

	std::vector<Point> points1, points2, points3, points4, points5, points6, points;
	std::vector<PointTypeITK> pointsITK;

	points1 = readPointsRegistrationTest(path1, points);
	points2 = readPointsRegistrationTest(path2, points);
	points3 = readPointsRegistrationTest(path3, points);
	points4 = readPointsRegistrationTest(path4, points);
	points5 = readPointsRegistrationTest(path5, points);
	points6 = readPointsRegistrationTest(path6, points);

	PointTypeITK realHip, realMedial, realKnee, virtualHip, virtualMedial, virtualKnee;
	realHip = Registration::makeItkPoint(92.3035, 82.9739, -109.181);
	realMedial = Registration::makeItkPoint(-90.9859, 122.137, -111.156);
	realKnee = Registration::makeItkPoint(-109.713, 121.427, -140.303);

	virtualHip = Registration::makeItkPoint(-26.44, 33.43, -92.89);
	virtualMedial = Registration::makeItkPoint(-32.30, 51.64, -298.00);
	virtualKnee = Registration::makeItkPoint(2.70, 40.64, -326.00);

	vtkSmartPointer<vtkPolyData> surface = TestVTK::ReadPolyData("D:\\3D_DICOM_DATA\\Pig_Cerdo\\Pig_left\\Pig_Segment_Knee_femur.vtk");
	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(surface);

	FemurRegistration femurRegis(surface, virtualHip, virtualKnee, virtualMedial);

	std::vector<double> result, goodPoints, goodPoints_2;
	float tGood = 0, tGood2 = 0;
	float resultGood = 0;
	for (int i = 0; i < 1000; i++)
	{
		getRandomPoints(points1, 6, pointsITK);
		getRandomPoints(points2, 6, pointsITK);
		getRandomPoints(points3, 6, pointsITK);
		getRandomPoints(points4, 6, pointsITK);
		getRandomPoints(points5, 3, pointsITK);
		getRandomPoints(points6, 6, pointsITK);

		femurRegis.MakeRegistration(pointsITK, realHip, realKnee, realMedial);
		double tError = femurRegis.getError();
		result.push_back(tError);

		if (tError <= 0.4)
		{
			resultGood++;
		}

		pointsITK.clear();

		/////////////////////////////////////

		auto transform = femurRegis.getTransformMarkerToCt();
		cv::Mat Rotation = Test::Rigid3DTransformToCVRotation(transform);
		cv::Mat Translation = Test::Rigid3DTransformToCVTranslation(transform);

		std::vector<cv::Point3d> points3d;
		tGood = 0;
		tGood2 = 0;
		for (int i = 0; i < points.size(); i++)
		{
			cv::Mat matPoint = Rotation * points[i].ToMatPoint() + Translation;
			Point tPoint = Point(matPoint);
			double myClosest[3];
			double arrayPoint[3] = { tPoint.x, tPoint.y, tPoint.z };
			double signedDistance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(arrayPoint, myClosest);
			if (abs(signedDistance) <= 0.4)
			{
				tGood++;
			}
			else if (abs(signedDistance) <= 0.6)
			{
				tGood2++;
			}
		}
		goodPoints.push_back(100.0 * tGood / float(points.size()));
		goodPoints_2.push_back(100.0 * tGood2 / float(points.size()));
	}

	std::sort(result.begin(), result.end(), greater());
	double suma = 0, sumaGood = 0, sumaGood2 = 0;
	for (int i = 0; i < result.size(); i++)
	{
		std::cout << "Error: " << result[i] << std::endl;
		suma = suma + result[i];
		sumaGood = sumaGood + goodPoints[i];
		sumaGood2 = sumaGood2 + goodPoints_2[i];
	}

	std::cout << "Average: " << suma / double(result.size()) << std::endl;
	std::cout << "Average good: " << 100.0 * resultGood / double(result.size()) << std::endl;
	std::cout << "Average Transform Points 0.4: " << sumaGood / double(goodPoints.size()) << std::endl;
	std::cout << "Average Transform Points 0.6: " << sumaGood2 / double(goodPoints_2.size()) << std::endl;

	/*
	Point translation(-885.604, -94.2699, 211.388);
	std::vector<double> rotationValues = { 0.410159, -0.20648, 0.888333, translation.x,
										   -0.473189, -0.880854, 0.0137382, translation.y,
											0.779655, -0.425984, -0.458995, translation.z,
										  0, 0, 0, 1};
	 cv::Mat rotation(4, 4, CV_64F, rotationValues.data());


	vtkNew<vtkImplicitPolyDataDistance> implicitPolyDataDistance;
	implicitPolyDataDistance->SetInput(surface);

	std::vector<cv::Point3d> points3d;
	int goodPoints = 0;
	for (int i = 0; i < points.size(); i++)
	{
		cv::Mat matPoint(4, 1, CV_64F);
		matPoint.at<double>(0, 0) = points[i].x;
		matPoint.at<double>(1, 0) = points[i].y;
		matPoint.at<double>(2, 0) = points[i].z;
		matPoint.at<double>(3, 0) = 1;

		matPoint = rotation.inv() * matPoint;
		Point tPoint = Point(matPoint.at<double>(0, 0), matPoint.at<double>(1, 0), matPoint.at<double>(2, 0));
		double myClosest[3];
		double arrayPoint[3] = {tPoint.x, tPoint.y, tPoint.z};
		double signedDistance = implicitPolyDataDistance->EvaluateFunctionAndGetClosestPoint(arrayPoint, myClosest);
		if (abs(signedDistance) < 0.4)
		{
			//std::cout << abs(signedDistance) << std::endl;
			points3d.push_back(tPoint.ToCVPoint());
			goodPoints++;
		}
		//points3d.push_back(tPoint.ToCVPoint());
	}

	LeastSquaresICP testICP(points3d);
	cv::Mat data(6, 1, CV_64F);

	data.at<double>(0, 0) = 0;
	data.at<double>(1, 0) = 0;
	data.at<double>(2, 0) = 0;
	data.at<double>(3, 0) = 0;
	data.at<double>(4, 0) = 0;
	data.at<double>(5, 0) = 0;

	//double error = testICP.LeastSquares(surface, data);

	//std::cout << "Least Square Error: " << error << std::endl;
	std::cout << "Good Per Cent: " << 100 * (float(goodPoints) / float(points.size())) << std::endl;

	TestVTK::show(surface, points3d);
	*/
}

double ex_euler(int n)
{
	if (n == 0)
	{
		return 0;
	}
	double h = 1. / n;
	int tCont = n;
	double errors = 0;
	double x = 0, y = 1., z = 1.;
	while (tCont > 0)
	{
		y = y + h * (2. - exp(-4. * x) - 2. * y);
		x = x + h;
		z = 1. + 0.5 * exp(-4. * x) - 0.5 * exp(-2. * x);

		errors += (std::abs(y - z) / z);
		tCont--;
	}

	errors = errors / (n + 1);
	errors = floor(errors * 1000000.0) / 1000000.0;

	return errors;
}

std::string printRoman(int num)
{
	int nums[13] = { 1, 4, 5, 9, 10, 40, 50, 90, 100, 400, 500, 900, 1000 };
	std::string roman[13] = { "I", "IV", "V", "IX", "X", "XL", "L", "XC", "C", "CD", "D", "CM", "M" };
	int tSize = 12;
	std::string result;
	while (num)
	{
		int div = num / nums[tSize];
		num %= nums[tSize];

		while (div)
		{
			result += roman[tSize];
			div--;
		}
		tSize--;
	}
	return result;
}

#include "Test_TKA_Match.h"
#include <vtkMassProperties.h>
#include <vtkIntersectionPolyDataFilter.h>
#include <vtkPolygon.h>
#include <vtkTriangle.h>


double getOverlappingArea(const vtkSmartPointer<vtkPolyData> containedPoly, const vtkSmartPointer<vtkPolyData> containerPoly, const vtkSmartPointer<vtkImplicitPolyDataDistance> implicitContainerDistance)
{
	double area = 0;

	for (vtkIdType cellId = 0; cellId < containedPoly->GetNumberOfCells(); ++cellId) {
		vtkCell* cell = containedPoly->GetCell(cellId);
		vtkPoints* points = cell->GetPoints();

		bool isContained = true;

		for (vtkIdType pointId = 0; pointId < points->GetNumberOfPoints(); ++pointId) {
			double point[3];
			points->GetPoint(pointId, point);

			if (implicitContainerDistance->FunctionValue(point) > 0)
			{
				isContained = false;
				break;
			}
		}

		if (isContained && points->GetNumberOfPoints() == 3) {

			double p0[3], p1[3], p2[3];
			points->GetPoint(0, p0);
			points->GetPoint(1, p1);
			points->GetPoint(2, p2);

			area += vtkTriangle::TriangleArea(p0, p1, p2);
		}
	}

	return area;
}

#include <vtkCubeSource.h>
#include <vtkTriangleFilter.h>

void PolydataInterception()
{
	vtkSmartPointer<vtkPolyData> sphere1 = TestVTK::CreateSphereTest(cv::Point3d(1,0,0), 2);
	vtkSmartPointer<vtkPolyData> sphere2 = TestVTK::CreateSphereTest(cv::Point3d(0, 0, 0), 2);

	vtkNew<vtkCubeSource> cubeSource;
	cubeSource->SetCenter(0, 0, 0);
	cubeSource->SetXLength(1);
	cubeSource->SetYLength(1);
	cubeSource->SetZLength(1);
	cubeSource->Update();

	auto pelvis = TestVTK::ReadPolyData("D:\\Mega_Trabajo\\Pelvis\\Segmentation.vtk");

	double pnt[3] = { 87, 6, 320 };
	double radius = 20;
	vtkNew<vtkSphereSource> sphere;
	sphere->SetCenter(pnt);
	sphere->SetRadius(radius);
	sphere->SetThetaResolution(sphere->GetThetaResolutionMaxValue() / 10);
	sphere->SetPhiResolution(sphere->GetPhiResolutionMaxValue() / 10);
	sphere->Update();

	/*vtkSmartPointer<vtkMassProperties> massProperties = vtkSmartPointer<vtkMassProperties>::New();
	massProperties->SetInputData(sphere->GetOutput());
	massProperties->Update();
	double area = massProperties->GetSurfaceArea();*/

	vtkNew<vtkImplicitPolyDataDistance> implicitContainerDistance;
	implicitContainerDistance->SetInput(pelvis);

	double intersection = getOverlappingArea(sphere->GetOutput(), pelvis, implicitContainerDistance);

	std::cout << "Area " << 0.5 * 4. * PI * radius * radius << std::endl;
	std::cout << "Subarea " << intersection << std::endl;

	//vtkSmartPointer<vtkTriangleFilter> triangleFilter = vtkSmartPointer<vtkTriangleFilter>::New();
	//triangleFilter->SetInputData(pelvis);
	//triangleFilter->Update();

	//vtkSmartPointer<vtkIntersectionPolyDataFilter> intersectionFilter = vtkSmartPointer<vtkIntersectionPolyDataFilter>::New();
	//intersectionFilter->SetInputData(0, sphere->GetOutput());
	//intersectionFilter->SetInputData(1, triangleFilter->GetOutput());
	//intersectionFilter->Update();
	//vtkSmartPointer<vtkPolyData> intersectionOutput = intersectionFilter->GetOutput();


	TestVTK::show_join(pelvis, sphere->GetOutput());
	//TestVTK::show(pelvis, intersectionOutput);

}

void Resgistration_General_test()
{
	
	Knee myKneeSource = CreateKneeFromFile_Numbers("D:\\sovajo\\Errores\\Error5", KneeSideEnum::KRight);

	Knee myKneeTarget = CreateKneeFromFile_Numbers("D:\\sovajo\\Errores\\Error2", KneeSideEnum::KLeft);


	std::vector<itk::Point<double, 3>> source = { myKneeSource.getHipCenter().ToITKPoint(),
												  myKneeSource.getLateralEpicondyle().ToITKPoint(),
												  myKneeSource.getMedialEpicondyle().ToITKPoint()};

	std::vector<itk::Point<double, 3>> target = { myKneeTarget.getHipCenter().ToITKPoint(), 
													myKneeTarget.getLateralEpicondyle().ToITKPoint(), 
													myKneeTarget.getMedialEpicondyle().ToITKPoint()};


	std::vector<itk::Point<double, 3>> rest_points = {myKneeSource.getLateralEpicondyle().ToITKPoint(),
													myKneeSource.getMedialEpicondyle().ToITKPoint(),
													myKneeSource.getMedialCondyle().ToITKPoint(), 
													myKneeSource.getLateralCondyle().ToITKPoint(), 
													myKneeSource.getFemurKneeCenter().ToITKPoint()};

	double error;
	auto regis = GENERAL::REGISTRATION::GeneralRegistrationPointsToPolydata(myKneeTarget.GetFemurPoly(), target, source);
	auto result = regis.MakeFinalAlignment(rest_points, error);
	std::cout << error << std::endl;
}

double getDistancePoints(const cv::Point3d& A, const cv::Point3d& B)
{
	auto diff = B - A;
	return sqrt(diff.dot(diff));
}

void PointsTo2D()
{
	std::vector<cv::Point3d> listaPoints = { cv::Point3d(-6.50229,76.922,-32.6216),
											cv::Point3d(-6.50229,76.4315,-31.1501),
											cv::Point3d(-6.50229,71.0362,-29.6786),
											cv::Point3d(-6.99277,66.6218,-28.2072),
											cv::Point3d(-6.50229,65.1503,-26.2452),
											cv::Point3d(-6.50229,63.6789,-24.2833),
											cv::Point3d(-6.0118,62.6979,-21.8309),
											cv::Point3d(-4.54035,61.2264,-19.8689),
											cv::Point3d(-4.54035,59.755,-17.907),
											cv::Point3d(-4.04986,58.774,-16.4355),
											cv::Point3d(-4.04986,57.793,-14.4736),
											cv::Point3d(-4.04986,56.3216,-11.5307),
											cv::Point3d(-4.04986,54.3596,-9.07825),
											cv::Point3d(-4.04986,53.8692,-7.60679),
											cv::Point3d(-4.04986,52.8882,-5.15436),
											cv::Point3d(-4.04986,52.8882,-3.19242),
											cv::Point3d(-4.04986,52.8882,-1.23048),
											cv::Point3d(-4.04986,52.8882,1.22195),
											cv::Point3d(-4.04986,52.8882,2.69341),
											cv::Point3d(-4.04986,52.8882,4.16486),
											cv::Point3d(-4.04986,52.8882,5.63632),
											cv::Point3d(-4.04986,52.8882,7.10778),
											cv::Point3d(-4.04986,52.8882,9.06972),
											cv::Point3d(-4.04986,52.8882,10.5412),
											cv::Point3d(-4.04986,52.8882,12.5031),
											cv::Point3d(-4.04986,52.8882,13.9746),
											cv::Point3d(-4.04986,52.8882,15.446),
											cv::Point3d(-4.04986,53.3787,17.8985),
											cv::Point3d(-4.54035,53.8692,19.3699),
											cv::Point3d(-6.0118,54.3596,21.8223),
											cv::Point3d(-6.0118,54.8501,23.2938),
											cv::Point3d(-6.50229,55.3406,25.2557),
											cv::Point3d(-6.99277,55.3406,26.7272),
											cv::Point3d(-6.99277,55.3406,28.1987),
											cv::Point3d(-6.99277,55.3406,30.1606),
											cv::Point3d(-7.97375,55.3406,31.6321),
											cv::Point3d(-7.97375,55.3406,33.1035),
											cv::Point3d(-7.97375,55.3406,35.5559),
											cv::Point3d(-9.93569,55.3406,38.4989),
											cv::Point3d(-9.93569,55.3406,39.9703),
											cv::Point3d(-9.93569,55.3406,41.4418),
											cv::Point3d(-10.9167,55.3406,42.9132),
											cv::Point3d(-11.4071,55.3406,45.3657),
											cv::Point3d(-11.4071,55.8311,47.3276),
											cv::Point3d(-11.4071,56.3216,48.7991),
											cv::Point3d(-11.8976,57.793,52.2325),
											cv::Point3d(-11.8976,57.793,53.7039),
											cv::Point3d(-12.8786,58.2835,55.1754),
											cv::Point3d(-12.8786,59.2645,57.1373),
											cv::Point3d(-13.3691,60.2455,59.0993),
											cv::Point3d(-13.3691,60.2455,61.0612),
											cv::Point3d(-13.3691,61.2264,64.0041),
											cv::Point3d(-13.3691,61.7169,65.966),
											cv::Point3d(-13.3691,62.2074,67.4375),
											cv::Point3d(-13.3691,62.2074,69.3994),
											cv::Point3d(-13.3691,62.6979,71.8519),
											cv::Point3d(-13.3691,63.6789,73.8138),
											cv::Point3d(-13.3691,64.1694,75.2853),
											cv::Point3d(-13.3691,64.6598,78.2282),
											cv::Point3d(-13.3691,65.1503,80.6806),
											cv::Point3d(-13.3691,65.1503,82.1521),
											cv::Point3d(-13.3691,65.6408,83.6235),
											cv::Point3d(-12.8786,65.6408,85.5855),
											cv::Point3d(-12.8786,65.6408,87.0569),
											cv::Point3d(-12.8786,67.1123,89.0189),
											cv::Point3d(-12.8786,67.1123,90.4903),
											cv::Point3d(-12.8786,68.0932,93.4332),
											cv::Point3d(-12.8786,68.5837,95.3952),
											cv::Point3d(-11.4071,69.0742,96.8666),
											cv::Point3d(-11.4071,70.0552,98.3381),
											cv::Point3d(-11.4071,70.0552,99.8096),
											cv::Point3d(-11.4071,71.5266,101.281),
											cv::Point3d(-11.4071,72.0171,102.752),
											cv::Point3d(-11.4071,72.5076,104.224),
											cv::Point3d(-11.4071,73.4886,106.186),
											cv::Point3d(-11.4071,73.4886,107.657),
											cv::Point3d(-11.4071,73.4886,109.129),
											cv::Point3d(-11.4071,73.4886,110.6),
											cv::Point3d(-11.4071,73.4886,112.072),
											cv::Point3d(-11.4071,72.5076,114.524),
											cv::Point3d(-11.4071,72.5076,116.486),
											cv::Point3d(-11.4071,72.0171,117.958),
											cv::Point3d(-11.4071,71.5266,119.429),
											cv::Point3d(-11.4071,71.5266,120.9),
											cv::Point3d(-11.4071,71.5266,122.372),
											cv::Point3d(-11.4071,71.5266,123.843),
											cv::Point3d(-11.4071,71.5266,125.315),
											cv::Point3d(-11.4071,71.5266,127.767),
											cv::Point3d(-11.4071,71.5266,129.239),
											cv::Point3d(-11.4071,71.5266,130.71),
											cv::Point3d(-11.4071,71.5266,132.182),
											cv::Point3d(-11.4071,71.5266,133.653),
											cv::Point3d(-11.4071,72.5076,135.125),
											cv::Point3d(-11.4071,72.5076,136.596),
											cv::Point3d(-11.4071,72.5076,138.067),
											cv::Point3d(-11.4071,73.4886,141.01),
											cv::Point3d(-11.4071,73.4886,142.972),
											cv::Point3d(-11.4071,74.4696,144.444),
											cv::Point3d(-11.4071,73.9791,145.915) };


}


void testSegSliceBordes()
{
	auto femurData = TestVTK::ReadPolyData("D:\\sovajo\\Test_Cases\\Segmentation_TKA\\test-segmentation\\femur.vtk");
	double normal[2][3] = { {1.0,0.0,0.0},{1.0,0.0,0.0} };
	double center[2][3] = { {37.93057714926232,198.0709686279297,-312.04998779296875},
							{126.2567716473318,198.0709686279297,-312.04998779296875} };
	TKA::SEGMENTATION::SPlane plane[2];
	plane[0].init(normal[0], center[0]);
	plane[1].init(normal[1], center[1]);
	TKA::SEGMENTATION::SliceBorder sliceBorder;
	auto slices = sliceBorder.MakeSlicesFromPolyData(femurData, plane[0], plane[1]);

	std::vector<std::vector<cv::Point3d>> lista;
	
	for (int i = 0; i < slices.size(); ++i)
	{
		auto points = sliceBorder.GetMainPoints(slices[i], 179);

		//std::cout << "Size: " << points.size() << std::endl;
		/*
		for (int j = 0; j < points.size(); j++)
		{
			if (std::isnan(points[j][0]))
			{

				std::cout << "i:NAN " << i << std::endl;
			}
			
		}
		*/

		std::vector<cv::Point3d> oneSlice;

		for (int j = 0; j < points.size(); j++)
		{
			oneSlice.push_back(cv::Point3d(points[j][0], points[j][1], points[j][2]));
		}
		lista.push_back(oneSlice);
	}
	

	//TestVTK::show(lista);
	for (int i = 9; i < 10; i++)
	{
		bool result;
		std::string msg;
		std::vector<vtkSmartPointer<vtkPolyData>> polySlices;
		polySlices.push_back(slices[i]);
		polySlices.push_back(slices[i+1]);
		auto surface = sliceBorder.GetSurface3D(polySlices, result, msg);
		//std::cout << "********************* " << i << std::endl;
		double pntDown[3];
		double pntUP[3];
		int pos1, pos2;
		if (i == 8)
		{
			pos1 = 170;
			pos2 = 170 + 251;
		}
		else
		{
			pos1 = 244;
			pos2 = 20 + 260;
		}
		//slices[i]->GetPoint(pos1, pntDown);
		//slices[i + 1]->GetPoint(pos2, pntUP);
		
		surface->GetPoint(pos1, pntDown);
		surface->GetPoint(pos2, pntUP);
		std::vector<cv::Point3d> nearPoints = {cv::Point3d(pntDown[0], pntDown[1], pntDown[2]), 
												cv::Point3d(pntUP[0], pntUP[1], pntUP[2])};

		std::cout << "Point down end: " << cv::Point3d(pntDown[0], pntDown[1], pntDown[2]) << ", Points up end: " << cv::Point3d(pntUP[0], pntUP[1], pntUP[2]) <<std::endl;
		
		TestVTK::show(surface, nearPoints);
		

		
	}
	//std::cout <<"********************* " << msg << std::endl;

}

#include "itkLabelImageToShapeLabelMapFilter.h"
#include "itkConnectedComponentImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkMaskImageFilter.h"

void segment_balls(int upperThreshold=100, int lowerThreshold=0, int minObjectSize=50)
{

	auto reader = itk::ImageFileReader<SPINE::SEGMENTATION::ImageType>::New();
	reader->SetFileName("D:\\sovajo\\Spine_Images\\avi_images\\DICOM_Format\\1_1.dcm");
	std::vector<itk::Point<double, 3>> centroid;
	SPINE::SEGMENTATION::BallsManualSegmentation balls(reader->GetOutput());
	SPINE::SEGMENTATION::SegmentImageType::Pointer result = balls.getSegmentBallsAndCentroids(lowerThreshold, upperThreshold, minObjectSize, centroid);

	/*
	using PixelType = unsigned short;
	using ImageType = itk::Image<PixelType, 3>;
	using LabelImageType = itk::Image<unsigned short, 3>;

	auto thresholdFilter = itk::BinaryThresholdImageFilter<ImageType, ImageType>::New();
	thresholdFilter->SetInput(reader->GetOutput());
	thresholdFilter->SetUpperThreshold(upperThreshold);
	thresholdFilter->SetLowerThreshold(lowerThreshold);
	thresholdFilter->SetInsideValue(1); 
	thresholdFilter->SetOutsideValue(0); 


	using LabelImageType = itk::Image<unsigned short, 3>;
	auto connectedFilter = itk::ConnectedComponentImageFilter<ImageType, LabelImageType>::New();
	connectedFilter->SetInput(thresholdFilter->GetOutput());

	auto relabelFilter = itk::RelabelComponentImageFilter<LabelImageType, LabelImageType>::New();
	relabelFilter->SetInput(connectedFilter->GetOutput());   
	relabelFilter->SetMinimumObjectSize(minObjectSize); 

	using LabelMapType = itk::LabelMap<itk::ShapeLabelObject<unsigned short, 3>>;
	auto labelMapFilter = itk::LabelImageToShapeLabelMapFilter<LabelImageType, LabelMapType>::New();
	labelMapFilter->SetInput(relabelFilter->GetOutput());
	labelMapFilter->Update();

	LabelMapType *labelMap = labelMapFilter->GetOutput();
	std::vector<unsigned long> sizes;
	for (unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); ++i) {
		auto labelObject = labelMap->GetNthLabelObject(i);
		sizes.push_back(labelObject->GetNumberOfPixels());
	}

	std::sort(sizes.begin(), sizes.end());
	unsigned long medianSize = sizes[sizes.size() / 2];
	//unsigned long rangeMin = static_cast<unsigned long>(medianSize * 0.8);
	unsigned long rangeMax = static_cast<unsigned long>(medianSize * 2);
	
	auto filteredLabels = LabelImageType::New();
	filteredLabels->SetRegions(relabelFilter->GetOutput()->GetLargestPossibleRegion());
	filteredLabels->Allocate();
	filteredLabels->FillBuffer(0);

	for (unsigned int i = 0; i < labelMap->GetNumberOfLabelObjects(); ++i) {
		auto labelObject = labelMap->GetNthLabelObject(i);
		unsigned long size = labelObject->GetNumberOfPixels();

		if (size <= rangeMax) {
			auto boundingBox = labelObject->GetBoundingBox();
			centroid.push_back(labelObject->GetCentroid());

			itk::ImageRegionIterator<LabelImageType> it(filteredLabels, boundingBox);
			for (it.GoToBegin(); !it.IsAtEnd(); ++it) {
				if (labelObject->HasIndex(it.GetIndex())) {
					it.Set(labelObject->GetLabel());
				}
			}
		}
	}

	///////////////////////////////////////////////////

	
	// Aplicar la máscara final
	auto maskFilter = itk::MaskImageFilter<ImageType, LabelImageType, ImageType>::New();
	maskFilter->SetInput(reader->GetOutput());
	maskFilter->SetMaskImage(filteredLabels);

	// Guardar la imagen resultante
	**/
	
	for (int i = 0; i < centroid.size(); i++)
	{
		std::cout << centroid[i] << std::endl;
	}

	std::cout << "Size: "<< centroid.size() << std::endl;
	
	auto writer = itk::ImageFileWriter<SPINE::SEGMENTATION::SegmentImageType>::New();
	writer->SetFileName("D:\\sovajo\\Spine_Images\\avi_images\\DICOM_Format\\output.dcm");
	writer->SetInput(result);

	try
	{
		writer->Update();
		std::cout << "Finishhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhhh" << std::endl;
	}
	catch (const itk::ExceptionObject & error)
	{
		std::cout << "Error: " << error.GetDescription() << std::endl;
	}

}


int main()
{
	testSegSliceBordes();
	//MatchEasy();
	//::testTibiaImplant2();
	//segment_balls();
	//RotulaGroovePath();
	std::cout << "VTK Version: " << vtkVersion::GetVTKVersion() << std::endl;
	/*
	std::string casePlan = "D:\\sovajo\\Cases_Plan_TKA\\case1_left";
	Knee myKnee;
	
	UKA::IMPLANTS::Point apLinePclPoint(10.4801, 9.06646, 0.934849);
	UKA::IMPLANTS::Point apLineTuberPoint(10.2438, -37.5031, 0.951554);
	UKA::IMPLANTS::Point sidePoint(-0.1021, -11.1526, 1.1545);
	UKA::IMPLANTS::Point exteriorPoint(-14.5152, -10.8294, 4.41175);
	UKA::IMPLANTS::TibiaImplantInfo tibiaInfo;
	tibiaInfo.tibiaThickness = 2.0;

	UKA::IMPLANTS::TibiaImplant tibiaImplant;
	tibiaImplant.init(apLinePclPoint, apLineTuberPoint, sidePoint, exteriorPoint, tibiaInfo);

	UKA::IMPLANTS::Point apLinePclPointSpacer(1.06782, -13.9567, -10.0878);
	UKA::IMPLANTS::Point apLineTuberPointSpacer(1.0678, 14.099, -10.1435);
	UKA::IMPLANTS::Point plateauRefPointUpSpacer(-2.50552, -0.14465, 0.412728);
	UKA::IMPLANTS::Point exteriorPointDownSpacer(1.06831, 0.435704, 0.00826925);

	UKA::IMPLANTS::TibiaSpacerImplant tibiaSpacer;
	tibiaSpacer.init(apLineTuberPointSpacer, apLinePclPointSpacer, plateauRefPointUpSpacer, exteriorPointDownSpacer);

	UKA::IMPLANTS::TibiaSpacerImplantMatch tibiaSpacerMatch;
	tibiaSpacerMatch.init(tibiaImplant, tibiaSpacer);

	auto tibiaImplant3D = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\tibia.STL");
	auto spacerImplant3D = TestVTK::ReadPolyDataSTL("D:\\sovajo\\Implants\\pka\\Tibial_spacer.STL");

	auto spacerTransform = TestVTK::TransformPoly(spacerImplant3D, tibiaSpacerMatch.GetRotationMatrix(), tibiaSpacerMatch.GetTranslationMatrix());
	*/
	/*
	vtkNew<vtkAppendPolyData> appendFilter;
	appendFilter->AddInputData(tibiaImplant3D);
	appendFilter->AddInputData(spacerTransform);
	appendFilter->Update();
	TestVTK::show(appendFilter->GetOutput());
	*/

	//bool result;

	
	SPINE::SEGMENTATION::ImageType::Pointer inputImg;

	std::string spine = "C:\\Users\\Miguel\\Desktop\\Algoritmos_Neuvos\\CT\\Spine_CT_17.nrrd";
	spine = "D:\\sovajo\\Case_Spine\\spine_ct.nrrd";
	Test::readImage<ImageType>(spine, inputImg);

	std::vector<cv::Point3d> listaPoints = GetSpineLineFromFile("D:\\sovajo\\Case_Spine");

	std::vector<ImageType::PointType> allPoints;

	for (auto& tempPoint: listaPoints)
	{
		ImageType::PointType physicalPoint;
		physicalPoint[0] = tempPoint.x;
		physicalPoint[1] = tempPoint.y;
		physicalPoint[2] = tempPoint.z;

		allPoints.push_back(physicalPoint);
	}

	SPINE::SEGMENTATION::SpineSegmentation segmentation(inputImg);
	ImageType::RegionType region;
	std::vector<SPINE::SEGMENTATION::SpineSegmentation::Plane> result = segmentation.getIntervertebralPlanes(allPoints, region);

	for (auto & item : result)
	{
		std::cout << "Center: " << item.center << std::endl;
	}
	
	
	//std::vector<SpineSegmentation::Plane> getIntervertebralPlanes(const std::vector<ImageType::PointType>& physicalPoints)

	std::string videoPathIn = "D:\\sovajo\\Spine_Images\\avi_images\\1_2.avi";
	std::string videoPathOut = "D:\\sovajo\\Spine_Images\\avi_images\\DICOM_Format";

	//videoToImage(videoPathIn, videoPathOut);

	//myKnee = CreateKneeFromFile_Numbers(casePlan, KneeSideEnum::KLeft);
	//double error = FemurRegistrationFromNumbersTKA(myKnee, casePlan, result);

	//std::cout <<"Result: "  << result << " Error: " << error << std::endl;

	//MatchEasyPKA();
	//RegistrationScale();
	//Resgistration_General_test();
	PelvisImplantMatch();
	//std::cout << "tttttttttttttttttt" << std::endl;
	//TestHullPoints();
	//Test30PointsVTK();

	//std::cout << Point(result) << "; " << proj << std::endl;

	//PolydataInterception();
	//TEST_PKA::MatchPKAThreePlanes();
	//TEST_PKA::MatchPKAThreePlanes();

	//double pnt[3] = { 0, 0, 0 };
	//Plane planeTemp;
	//planeTemp.init(Point(1, 0, 0), Point(1, 0, 0));

	//vtkNew<vtkSphereSource> sphere;
	//sphere->SetCenter(pnt);
	//sphere->SetRadius(2);
	//sphere->SetThetaResolution(sphere->GetThetaResolutionMaxValue() / 2);
	//sphere->SetPhiResolution(sphere->GetPhiResolutionMaxValue() / 2);
	//sphere->Update();

	//vtkNew<vtkPlane> vtkPlaneA;

	//auto planeNormal = planeTemp.getNormalVector();
	//auto planePoint = planeTemp.getPoint();
	//vtkPlaneA->SetOrigin(planePoint.x, planePoint.y, planePoint.z);
	//vtkPlaneA->SetNormal(planeNormal.x, planeNormal.y, planeNormal.z);

	//vtkNew<vtkPlaneCollection> cutPlanes;
	//cutPlanes->AddItem(vtkPlaneA);

	//vtkNew<vtkClipClosedSurface> Clipper;
	//Clipper->SetInputData(sphere->GetOutput());
	//Clipper->SetClippingPlanes(cutPlanes);
	//Clipper->Update();

	//auto hemiSphere = Clipper->GetOutput();

	//std::cout << sphere->GetThetaResolution() << ", " << sphere->GetThetaResolutionMaxValue() << ", "<< sphere->GetPhiResolution() << std::endl;

	//TestVTK::show(sphere->GetOutput(), hemiSphere);





	//TEST_IMPLANTS::testImplant();

	//Plane A, B;

	//A.init(Point(0, 1, 0), Point(0, 0, 0));
	//B.init(Point(3, 1, 2), Point(1, 2, 3));
	//Point vectorX(-1, 3, 10);
	//vectorX.normalice();

	//Point Aproj = A.getProjectionVector(vectorX);
	//Aproj.normalice();
	//
	//Point projTemp = ((Aproj.dot(B.getNormalVector()))/(B.getNormalVector().dot(B.getNormalVector()))) * B.getNormalVector();

	//Point vectorY = Aproj - projTemp;
	//vectorY.normalice();

	//Point Aproj2 = A.getProjectionVector(vectorY);
	//Aproj2.normalice();

	//std::cout << "pertence: " << B.getNormalVector().dot(vectorY) <<std::endl;
	//std::cout << "Aproj" << Aproj << std::endl;
	//std::cout << "vectorY" << vectorY << std::endl;
	//std::cout << "Aproj2" << Aproj2 << std::endl;

	std::string path, pig, pathHoles;
	path = "D:\\work\\Varilla\\split_image.nrrd";
	//path = "D:\\work\\Varilla\\Fei\\bone.nrrd";
	path = "D:\\fei_error\\problem4\\dead_knee.nrrd";
	//path = "D:\\3D_DICOM_DATA\\Person_1\\nrrd\\knee.nrrd";
	path = "C:\\Users\\miguel\\Downloads\\knee.nrrd";

	pathHoles = "C:\\Users\\miguel\\Downloads\\Holes_Test.nrrd";

	pig = "D:\\3D_DICOM_DATA\\Pig_Cerdo\\Pig_Right_2\\Tibia.nrrd";

	std::string imageHoleStr = "D:\\Mega_Trabajo\\Person_2\\Raw_Right_Femur_Seg.nrrd";
	std::string imageHoleStr2 = "D:\\Mega_Trabajo\\Person_2\\Femur_Right_very_smooth.nrrd";

	//SegmentImageType::Pointer imageHole;
	//Test::readImage<SegmentImageType>(imageHoleStr2, imageHole);

	////auto close1 = Test::CloseImage3(imageHole);
	//auto close1 = Test::FillHole<SegmentImageType>(imageHole);

	//Test::SaveImage<SegmentImageType>(close1, "D:\\Mega_Trabajo\\Person_2\\close2.nrrd");

	//TestVTK::SavePolyData(close1, "D:\\Mega_Trabajo\\Person_2\\close1.vtk");

	/*
	Test::SaveImage<SegmentImageType>(close1, "D:\\Mega_Trabajo\\Person_2\\close1.nrrd");

	Test::CloseImage2(imageHole);

	Test::SaveImage<SegmentImageType>(imageHole, "D:\\Mega_Trabajo\\Person_2\\close2.nrrd");*/


	/*SegmentImageType::Pointer myImage;

	RegistrationImageType::Pointer imageCT, imageMRI, imageOut;*/

	//PelvisImplantMatch();

	//getPelvisPoint();

	//HipFemoralRegistration();

	/*auto pelvis = TestVTK::itkImageToSurface3D("D:\\Mega_Trabajo\\Pelvis\\Segmentation_Left.nrrd");
	TestVTK::SavePolyData(pelvis, "D:\\Mega_Trabajo\\Pelvis\\Segmentation_Left.vtk");*/

	/*Test::readImage<RegistrationImageType>("D:\\work\\SPINE\\Spine_CT.nrrd", imageCT);
	Test::readImage<RegistrationImageType>("D:\\work\\SPINE\\Spine_MRI.nrrd", imageMRI);
	auto transform = GeneralRegistrationPointsToPolydata::MultiResImageRegistration(imageCT, imageMRI, imageOut);
	std::cout << transform << std::endl;
	Test::SaveImage<ImageType>(imageOut, "Registration_Out_IMG");*/

	//auto poly = ManualSegmentation::BinaryImageToPolyData(myImage);

	//TestVTK::SavePolyData(poly, "Suen_Poly");

	/*
	ImageType::Pointer myImage, outImage, outSlice;

	Test::readImage<ImageType>(pathHoles, myImage);

	Test::Get2DSlice<ImageType>(myImage, outSlice, 101);

	Test::SaveImage<ImageType>(outSlice, "slice_hole_fill_prev1");
	outImage = ManualSegmentation::FloodFillBinaryImageSlice(outSlice);
	Test::SaveImage<ImageType>(outImage, "slice_hole_fill_after1");*/

	//newTiling();

	/*
	try {
		function_f();
	}
	catch (MyException e)
	{
		std::cout << "Exepcion: " << e << std::endl;
	}*/

	/*
	Test::readImage<SegmentImageType>(path, myImage);
	SegmentImageType::Pointer fillImage = Test::FillHole<SegmentImageType>(myImage);
	Test::SaveImage<SegmentImageType>(fillImage, "test_bone_willian");*/

	//MatchEasy();

	//LS_Registration_Femur();

	//RegistrationScale();

	//executeSegmentation(path);

	//GenerateFilesJSON();



   // CheckImageRod checkRod(myImage);
   // if (checkRod.Execute())
   // {
   //     std::cout << "Fine!!!!";
   // }

	//Test::SaveImage<ImageType>(check.mTest, "Test_Rod_new");


	/*Plane splitPlane;
	splitPlane.init(Point(1, 0, 0), Point(250, 0, 225.4));
	EraseImagePartUsingPlane(myImage, splitPlane);
	Test::SaveImage<ImageType>(myImage, "split_image");*/

	//Fusion_Test::SaveImage<ImageType>(imageCT, "tttt");

	//Fusion_Test::fusion(imageCT, imageMRI);

	//Test::SaveImage<ImageType>(image, "Test_Series2");

	//MatchEasy();

	//executeSegmentation("");

	//executeImplantsMatch();

	//TestHullPoints();
	//GenerateFilesJSON();

	//GenerateDataSet(kind_dataset::TRAINING, 4, true, 3791);
	//GenerateDataSet(kind_dataset::TESTING, 4, true);

	/*std::vector<Point> a = getSpherePoints("D:\\workspace\\bone.txt");
	std::vector<Point> b = getSpherePoints("D:\\workspace\\hip.txt");
	std::vector<Point> c = getSpherePoints("D:\\workspace\\test.txt");
	std::vector<Point> d = getSpherePoints("D:\\workspace\\test2.txt");

	auto result1 = ImplantTools::fitSphere(a);
	auto result2 = ImplantTools::fitSphere(b);
	auto result3 = ImplantTools::fitSphere(c);
	auto result4 = ImplantTools::fitSphere(d);

	std::cout << result1.first << "    "<< result1.second <<  std::endl;
	std::cout << result2.first << "    " << result2.second << std::endl;
	std::cout << result3.first << "    " << result3.second << std::endl;
	std::cout << result4.first << "    " << result4.second << std::endl;*/

	//RegistrationScale();
	//PelvisImplantMatch();
	//KneeimplantMatch();
	//TestSlice();
	//PelvisImplantMatch_2();

	//segmentationDNN();

	//std::cout << CV_VERSION <<std::endl;

	/*std::vector<cv::Point3d> myPoints = { {-186.3414,-1501.3035,-80.8088},
{-183.7069,-1504.7179,-80.3374},
{-181.4981,-1506.221,-83.6954},
{-181.9024,-1504.2946,-87.595},
{-184.5381,-1500.9026,-88.0484},
{-186.7279,-1499.393,-84.7183},
{-186.3348,-1501.3098,-80.8225},
{-202.2921,-1519.2171,-207.4927},
{-183.2787,-1501.7539,-146.7955},
{-186.2205,-1501.3323,-81.0919},
{-210.6344,-1518.4702,-22.3104},
{-251.9301,-1549.8698,18.4318},
{-186.3029,-1501.3146,-80.8461},
{-288.2075,-1446.2173,-94.9103},
{-237.9902,-1460.6898,-86.2577},
{-186.2842,-1501.3197,-80.8458},
{-156.0011,-1559.8731,-81.9131},
{-152.4866,-1625.1453,-89.1679},
{-186.3243,-1501.3332,-80.8135},
{-110.5485,-1440.0069,-58.4176},
{-186.3156,-1501.3471,-80.8126},
{-156.9563,-1483.8357,-227.1448} };

	auto result = Test::fitSphere(myPoints);

	std::cout << "center: " << result.first << "  radius: " << result.second << std::endl;*/

	/*Point A(-232.9108, -1639.5415, -259.8053);
	Point B(-196.5311, -1660.1665, -195.3617);

	Line myLIne12 = Line::makeLineWithPoints(Point(-232.9495, -1638.4385, -260.2543), Point(-196.0186, -1660.1991, -196.4777));

	Point C(-253.732, -1683.9716, -247.7632);
	Point D(-232.0117, -1685.5656, -184.1435);

	Line myLIne34 = Line::makeLineWithPoints(Point(-253.9221, -1683.3173, -249.1907), Point(-231.8522, -1685.2375, -185.3496));

	double a = myLIne12.getDistanceFromPoint(A);
	double b = myLIne12.getDistanceFromPoint(B);
	double c = myLIne34.getDistanceFromPoint(C);
	double d = myLIne34.getDistanceFromPoint(D);

	std::cout << "A: " << a << "\n" << "B: " << b << "\n" << "C: " << c << "\n" << "D: " << d << "\n";*/

	/*Plane A, B, C, D, E, F;
	A.init(Point(-0.0831969, 0.0177094, 0.996376), Point(43.5555, -43.1478, 59.9849));
	B.init(Point(-0.0153994, -0.073256, 0.997194), Point(46.2422, -37.4291, -28.0257));

	C.init(Point(0.000162398, -0.000141977, 1), Point(-1.58163, -38.5354, 59.4116));
	D.init(Point(-0.105419, 0.088994, 0.990438), Point(1.62674, -38.378, -29.0893));

	E.init(Point(-0.00160365, 3.04789e-05, 0.999999), Point(-45.7437, -39.7215, 58.867));
	F.init(Point(2.01723e-05, -0.000134865, 1), Point(-43.3233, -37.9347, -29.6941));

	A.reverseByPoint(B.getPoint(), false);
	B.reverseByPoint(A.getPoint(), false);

	C.reverseByPoint(D.getPoint(), false);
	D.reverseByPoint(C.getPoint(), false);

	E.reverseByPoint(F.getPoint(), false);
	F.reverseByPoint(E.getPoint(), false);

	Point a = A.getProjectionPoint(Point(43.5863, -43.3087, 59.9543)) + 0.5 * A.getNormalVector();
	Point b = (Point(46.2644, -37.2907, -27.5038)) + B.getNormalVector();
	Point c = (Point(-1.59243, -38.3661, 58.9566)) + C.getNormalVector();
	Point d = (Point(1.43764, -38.6334, -28.5982)) + D.getNormalVector();
	Point e = (Point(-45.7503, -39.5615, 58.3166)) + E.getNormalVector();
	Point f = (Point(-43.3424, -37.9466, -29.2181)) + F.getNormalVector();

	std::cout << "A: " << a << "\n" << "B: " << b << "\n" << "C: " << c << "\n" << "D: " << d << "\n" <<"E: " << e << "\n" << "F: " << f << "\n";

*/
/* std::vector<int> data1, data2;
 srand(time(NULL));
 int tSize = 10000000;
 for (int i = 0; i < tSize; i++)
 {
	 int var = rand() % tSize;
	 data1.push_back(var);
	 data2.push_back(var);
 }

 auto start = std::chrono::system_clock::now();

 QuickSort(data1, 0, data1.size());

 auto end = std::chrono::system_clock::now();
 std::chrono::duration<double> elapsed_seconds = end - start;
 std::cout << "Mine time: " << elapsed_seconds.count() << "s\n";

 start = std::chrono::system_clock::now();
 std::sort(data2.begin(), data2.end());
 end = std::chrono::system_clock::now();

 elapsed_seconds = end - start;
 std::cout << "STD time: " << elapsed_seconds.count() << "s\n";

 for (int i = 0; i < tSize; i++)
 {
	 if (data1[i] != data2[i])
	 {
		 std::cout << "Esta mal " <<"s\n";
	 }
 }

 std::cout << "Finish " << "s\n";*/


 /*for (auto& elem : data)
 {
	 std::cout << elem << std::endl;
 }*/

 //readExample();
 //HipFemoralRegistration();
 //GetHipPoints();

 //CreateKneeFromFile("D:\\3D_DICOM_DATA\\Person_2\\Left");

 //TestPointsOnBone();
 //RegistrationScale();
 //executeSegmentation();
 //getSlicesAndSurfaceBone();
 //TestVTK::calculateDistance();

 //executeImplantsMatch();
 //executeBalance();

 //getPelvisPoint();

 /*Point vectorTest(2, 3, 1);
 Plane XZ;
 XZ.init(Point(0, 1, 0), Point(0, 0, 0));
 Point proj1 = XZ.getProjectionVector(vectorTest);

 cv::Mat rotation = Line::getRotateMatrix(Point(0, 1, 0), (30.0 * PI / 180.0));
 cv::Mat vectorTestRotate = rotation * vectorTest.ToMatPoint();
 Point test2 = XZ.getProjectionVector(Point(vectorTestRotate));

 cv::Mat test1Mat = rotation * proj1.ToMatPoint();
 Point test1 = (Point(test1Mat));

 test2.normalice();
 test1.normalice();

 std::cout << test1 << "  " << test2 << " "<< Line::getAngleBetweenVectors(vectorTest, Point(vectorTestRotate)) * 180.0 / PI<< std::endl;*/

 //getSlicesAndSurfaceBone();
 //tiling();
 //testTree(4);

 //newTiling();

 //oldTilling();

 /*cv::Point3d a = { 2, 1, 0 };
 cv::Point3d b = { -2, -1, 0 };*/
 //cv::Point3d c = { 5, 2, 0 };

 /*double sign = a.x -

  double area = (sqrt(V3.dot(V3))) / 2.0;*/
  //std::cout << "cross: " << a.cross(b) << std::endl;




  /*SegmentImageType::Pointer imgOut;
  std::string imgName = "D:\\3D_DICOM\\right_femur_segment.nrrd";

  Test::readImage<SegmentImageType>(imgName, imgOut);
  Segmentation::createBoneSurface(imgOut);*/

  //crossValidationTibiaRegistration();

  //femurRegistrationFast();
  //tibiaRegistrationFast();
  //Segmentation::ContourToSufarce();
  //Segmentation::GetContour();
  //checkAngles();
  //executeSegmentation();
  //executeHip();
  /*
  auto start = std::chrono::high_resolution_clock::now();
  double sumaTibia = 0.0;
  double sumaFemur = 0.0;
  for (int i = 0; i < 50; i++)
  {
	  sumaTibia = sumaTibia + tibiaRegistration();
  }

  for (int i = 0; i < 50; i++)
  {
	  sumaFemur = sumaFemur + femurRegistration();
  }
  auto finish = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> elapsed = finish - start;
  std::cout << "Time elapsed: " << elapsed.count() << std::endl;

  std::cout << "**** Average tibia: " << sumaTibia / 50.0 << std::endl;
  std::cout << "**** Average femur: " << sumaFemur / 50.0 << std::endl;*/

  //tibiaRegistration();
  //femurRegistration();
  //executeImplantsMatch();
  //DrawBonePlanes();
  //executeBalance();

 //LS_Rotula_Registration();

 //vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::itkImageToSurface3D(dicom_knee_cap);
 /*vtkSmartPointer<vtkPolyData> BoneVTK = TestVTK::ReadPolyData(vtk_kneeCap_right);
 vtkSmartPointer<vtkPolyData> clean = TestVTK::CleanPoly(BoneVTK);
 TestVTK::SavePolyData(clean, "rotula_right");*/


 //LS_Registration_Femur();
 //crossValidationTibiaRegistration();

 //LS_Registration_Femur_Random();

 //newTiling();

 //executeSegmentation();
 //std::cout << cv::getBuildInformation().c_str() << std::endl;

 //ExecuteKneeCap();

 //executeImplantsMatch();

 //BoneRotulaPath::FitPoints();

 //Coherent();

 //LS_Registration_Femur();

 //femurRegistrationFast();
 //tibiaRegistrationFast();

 //ImageType::Pointer myHip = Test::ReadSeries<ImageType>("D:\\new_ct_data\\person\\hip");
 //ImageType::Pointer myAnkle = Test::ReadSeries<ImageType>("D:\\new_ct_data\\person\\ankle");
 //ImageType::Pointer myKnee = Test::ReadSeries<ImageType>("D:\\new_ct_data\\person\\knee");


 //Test::SaveImage<ImageType>(myHip, "myHip");
 //Test::SaveImage<ImageType>(myAnkle, "myAnkle");
 //Test::SaveImage<ImageType>(myKnee, "myKnee");

 //vtkSmartPointer<vtkPolyData> bone = TestVTK::itkImageToSurface3D("D:\\3D_DICOM\\Left_Leg\\new\\Femur_good_seg.nrrd");
 //TestVTK::SavePolyData(bone, "My_famur.vtk") ;

 //ChangeCoordenate();

 //Test30PointsVTK();

 //TestSlice();

 //newTiling();

 //LS_Rotula_Registration();

 /*std::string source, target;
 source = "D:\\Excel\\Files";
 target = "D:\\Excel\\TargetFiles";

 LoadFile myFiles("brand", "model", source, target);
 std::vector<std::string> errorList, overwriteWarning;
 myFiles.readFiles(errorList, overwriteWarning);

 for (int i = 0; i < errorList.size(); i++)
 {
	 std::cout << errorList[i] << std::endl;
 }

 for (int i = 0; i < overwriteWarning.size(); i++)
 {
	 std::cout << overwriteWarning[i] << std::endl;
 }

 errorList.clear();
 myFiles.writeFiles(errorList, true);

 for (int i = 0; i < errorList.size(); i++)
 {
	 std::cout << errorList[i] << std::endl;
 }
*/
//executeSegmentation();

//getPelvisPoint();

/*std::string line;
ifstream myfile("D:\\Algorith_Project\\build\\femur_data.txt");
std::vector<std::string> data;
int cont = 0;
if (myfile.is_open())
{
	while (getline(myfile, line))
	{
		data.push_back(line);
	}
	if (data.size() % 1000000 == 0)
	{
		cont++;
		std::cout << cont << std::endl;
	}
	myfile.close();
}

for (int i = data.size() - 10; i < data.size(); i++)
{
	std::cout << data[i] << std::endl;
}*/

/*Point A(59.83, -76.34, -639.01);
Point B(97.56, -69.14, -636.91);
Point C(97.56, -76.34, -636.91);

Point latCondy(42.46, -93.72, -629.91);
Point medCondy(113.24, -85.67, -633.41);

Plane axial;
axial.init(A, B, C);

Point TEA = latCondy - medCondy;
TEA = axial.getProjectionVector(TEA);

Point PCA = A - B;

double angle = Line::getAngleBetweenVectors(TEA, PCA);
angle = angle * 180.0 / PI;

std::cout << "angle: " << angle << std::endl;*/

//vtkSmartPointer<vtkPolyData> femur = TestVTK::itkImageToSurface3D("D:\\3D_DICOM_DATA\\Modo\\Right_Modo\\Segmentation_nrrd\\Hip.nrrd");
//vtkSmartPointer<vtkPolyData> tibia = TestVTK::itkImageToSurface3D("D:\\3D_DICOM\\VTK_BONE\\tibia_right_vtk_.nrrd");

 //TestVTK::SavePolyData(femur, "hip_right.vtk");
 //TestVTK::SavePolyData(tibia, "tibia_new");

//std::cout << sin(0.78) * sin(0.78) + 0.5 << std::endl;

   /*cv::Mat a1 = cv::Mat::eye(3, 3, CV_64F);
   cv::Mat a2 = cv::Mat::eye(3, 3, CV_64F);
   cv::Mat a3 = cv::Mat::eye(3, 3, CV_64F);

   a1.at<double>(0, 0) = 5;
   a1.at<double>(1, 1) = 0;
   a1.at<double>(2, 2) = 0;

   a2.at<double>(1, 1) = 10;
   a3.at<double>(2, 2) = 20;

   std::cout << a1 * a2 * a3 << std::endl;*/

	std::cout << "Finish .... press some letter" << std::endl;
	char cc;
	std::cin >> cc;
	return EXIT_SUCCESS;

}