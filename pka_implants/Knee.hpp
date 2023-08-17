#ifndef IMPLANT_KNEE_H
#define IMPLANT_KNEE_H

#include "Plane.hpp"
#include "Utils.hpp"
#include <vector>
#include "Types.hpp"
#include "pka_implants_export.h"
#include "vtkPolyData.h"
#include "vtkSmartPointer.h"
//#include <fstream>

namespace PKA
{
	namespace IMPLANTS
	{
		struct TibiaImplantDimensions
		{
			double long_axis;
			double short_axis;
		};

		enum KneeSideEnum
		{
			KRight,
			KLeft
		};

		class PKA_IMPLANTS_EXPORT Knee
		{
		public:
			Knee();

			void init(const Point& hipCenter, const Point& anteriorCortex, const Point& femurKneeCenter, const Point& lateralEpicondyle,
				const Point& medialEpicondyle, const Point& tibiaKneeCenter, const Point& tibiaTubercle, const Point& pclCenter,
				const Point& ankleCenter, const Patella& pPatella, const vtkSmartPointer<vtkPolyData> femurPoly,
				const vtkSmartPointer<vtkPolyData> tibiaPoly, KneeSideEnum pSide, bool findRefPoints = true, double cartilage = 2.0, uint8_t imageValueMax = 1);

			///////////////////////////////////////////////////
			/*
				Be very careful with these functions, if the "findRefPoints" parameter of the init function is false, 
				then two of these methods, one for femur and one for tibia, must be called immediately after the init. 
				By default if init is True then it uses the function based on least squares for femur.
			*/

			void findFemurCondylesUsingDistalPoints();

			void findFemurCondylesUsingLeastSquare();

			void findAutomaticPlateaus();
			//////////////////////////////////////////////////
			
			Point getFemurVectorTEA() const;

			Point getTibiaVectorTEA() const;

			Point getFemurVectorLateralTEA() const;

			Point getTibiaVectorLateralTEA() const;

			Point getDirectVectorFemurAxis() const;

			Point getFemurDirectVectorAP() const;

			Point getTibiaDirectVectorAP() const;

			Point getNormalVectorTibiaPlane() const;

			Point getAnteriorCortex() const;

			Point getMoveCondyle(const FemurImplantInfo& pImplant) const;

			Point getFemurKneeCenter() const;

			Point getTibiaKneeCenter() const;

			Point getLateralCondyle() const;

			Point getMedialCondyle() const;

			Point getLateralPlateau() const;

			Point getMedialPlateau() const;

			Point getMovePlateau(const TibiaImplantInfo& pImplant) const;

			Point getAnkleCenter() const;

			Point getLateralEpicondyle() const;

			Point getMedialEpicondyle() const;

			Point getMedialEpicondylePerp() const;

			Point getHipCenter() const;

			Point getTibiaTubercle() const;

			Point getLateralInferiorFemurPoint() const;

			Point getMedialInferiorFemurPoint() const;

			void setLateralAndMedialInferiorFemurPoints(const Point& pLateral, const Point& pMedial);

			void setLateralAndMedialPosteriorFemurPoints(const Point& pLateral, const Point& pMedial);

			void setLateralAndMedialPlateauPoints(const Point& pLateral, const Point& pMedial);

			Point getInferiorMoveFemurPoint(const FemurImplantInfo& pImplant) const;

			Point getTibiaRotatePoint() const;

			Point getPclCenterPoint() const;

			cv::Mat getTibiaCenterPointApAutomatic() const;

			cv::Mat getTibiaCenterPointOnImplantAP(const TibiaImplantInfo& pImplant) const;

			Point getMoveTibiaKneeCenter(const TibiaImplantInfo& pImplant) const;

			Point getMoveFemurKneeCenter(const FemurImplantInfo& pImplant) const;

			Point getPatellaRotationPoint() const;

			//static void readImageTest(std::string path , ImplantImageType::Pointer& image);

			static Point getComputeAnkleCenter(const Point& lateralMalleolus, const Point& medialMalleolus);

			vtkSmartPointer<vtkPolyData> GetFemurPoly() const;

			vtkSmartPointer<vtkPolyData> GetTibiaPoly() const;

			vtkSmartPointer<vtkPolyData> GetPatellaPoly() const;

			double getEstimateFemurImplantsDimensions() const;

			TibiaImplantDimensions getEstimateTibiaImplantsDimensions() const;

			double getFemurCartilage() const;

			double getTibiaCartilage() const;

			bool getIsVarus() const;

			/*Point getLateralCoronalDistalPoint() const;

			Point getMedialCoronalDistalPoint() const;

			Point getTopPointOnPatellaPath() const;*/

			void setTibiaSlope(double angleInDegrees);

			std::vector<Point> getKneeGroovePath() const; //Now it is returning an empty list ************************************

			//std::vector<Point> getKneeGrooveOutLiers() const;

			bool getIsRight() const;

			Plane getEquisPlaneTibia() const;

			std::pair<double, double> getTibiaAutomaticAxis(Point& pVectorApFront, Point& pVectorTeaLat, Point& pCenterPoint) const;

			Plane getPatellaFrontPlane() const;

			Point getPatellaCenter() const;

			Point getPatellaDistalPosteriorPoint() const;

			double getPatellaThickness() const;

			double getPatellaDiameter() const;

			void setFemurCartilage(double pCartilage);

			void setTibiaCartilage(double pCartilage);

			Patella getPatella() const;

			static double getInitialAnglePCA(const Point& hipCenter, const Point& femurKneeCenter, const Point& lateralEpicondyle,
				const Point& medialEpicondyle, const Point& lateralCondyle, const Point& medialCondyle);
		private:
			Point lateralEpicondyle, medialEpicondyle, medialEpicondylePerp;
			Point hipCenter, femurKneeCenter, tibiaKneeCenter, tibiaRotatePoint;
			Point anteriorCortex;
			Point lateralInferiorFemurPoint, medialInferiorFemurPoint;
			Point lateralCondyle, medialCondyle;
			Point lateralPlateau, medialPlateau;
			Point tibiaTubercle;
			Point pclCenter;
			Point ankleCenter;
			Point femurDirectVectorAP;
			Point coronalDistalLat, coronalDistalMed;
			Point topPointOnPatellaPath;
			Side goodSide;
			double femurCartilage, tibiaCartilage; // , femurLatMedDiffDistal, femurLatMedDiffPosterior, tibiaLatMedDiff;
			Point tibiaNormalPlaneVector;
			Plane mPatellaPlane;
			Point mPatellaDistalPosteriorPoint;
			std::vector<Point> mFemur, mTibia, mKneeGroovePath, mKneeGrooveOutLiers;
			vtkSmartPointer<vtkPolyData> femurPoly, tibiaPoly;
			Patella mPatella;

			bool isInit;

			bool isVarus;

			bool isRight;

			void findTibiaPlaneNormalVector(double degreeAngle = 5.0);

			//void findTibiaDirectVectorAP();

			//void readImagePoints(const ImplantImageType::Pointer& image, uint8_t imageValueMax, std::vector<Point>& meshPoints) const;

			void FillFemurPoints();

			void FillTibiaPoints();

			void UpdateTopPointOnGroove();

			Point getPointAtSquareDistance(const Line& line, const Point& pPoint, const Point& nearReferencePoint, float distance, bool closest = true) const;

			Point getPointAtDistance(const Point& directVector, const Point& fixPoint, const Point& nearReferencePoint, float distance, bool closest = true) const;

			//void getDistalInferiorCondyleOld(Point& lateralPoint, Point& medialPoint) const;

			void getDistalInferiorCondyle(Point& lateralPoint, Point& medialPoint);

			//Point getDistalInferiorCondylePoint(const cv::RotatedRect& box, const cv::Mat& rotate, const Point& epicondyle, const double axisZ) const;

			//void get2DSlice(const ImplantImageType::Pointer& imageIn, ImplantImageType::Pointer& imageOut, const int sliceNumber) const;

			Side getGoodSide(const Point& hipCenter, const Point& kneeCenter, const Point& latEpi, const Point& medEpi, const Point& ankle) const;

			void makeKneeGroovePath();
		};
	}
}

#endif
