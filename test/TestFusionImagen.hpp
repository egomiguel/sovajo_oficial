#include "itkImageRegistrationMethodv4.h"
#include "itkMeanSquaresImageToImageMetricv4.h"

#include "itkVersorRigid3DTransform.h"
#include "itkCenteredTransformInitializer.h"

#include "itkRegularStepGradientDescentOptimizerv4.h"
// Software Guide : EndCodeSnippet

#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkResampleImageFilter.h"
#include "itkCastImageFilter.h"
#include "itkSubtractImageFilter.h"
#include "itkRescaleIntensityImageFilter.h"
#include "itkExtractImageFilter.h"

//  The following section of code implements a Command observer
//  that will monitor the evolution of the registration process.
//
#include "itkCommand.h"

namespace Fusion_Test
{

    class CommandIterationUpdate : public itk::Command
    {
    public:
        using Self = CommandIterationUpdate;
        using Superclass = itk::Command;
        using Pointer = itk::SmartPointer<Self>;
        itkNewMacro(Self);

    protected:
        CommandIterationUpdate() = default;

    public:
        using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
        using OptimizerPointer = const OptimizerType *;
        void
            Execute(itk::Object * caller, const itk::EventObject & event) override
        {
            Execute((const itk::Object *)caller, event);
        }
        void
            Execute(const itk::Object * object, const itk::EventObject & event) override
        {
            auto optimizer = static_cast<OptimizerPointer>(object);
            if (!itk::IterationEvent().CheckEvent(&event))
            {
                return;
            }
            std::cout << optimizer->GetCurrentIteration() << "   ";
            std::cout << optimizer->GetValue() << "   ";
            std::cout << optimizer->GetCurrentPosition() << std::endl;
        }
    };

    int fusion(std::string imageFix, std::string imageMove)
    {

        constexpr unsigned int Dimension = 3;
        using PixelType = float;

        using FixedImageType = itk::Image<PixelType, Dimension>;
        using MovingImageType = itk::Image<PixelType, Dimension>;

        using TransformType = itk::VersorRigid3DTransform<double>;
        // Software Guide : EndCodeSnippet

        using OptimizerType = itk::RegularStepGradientDescentOptimizerv4<double>;
        using MetricType = itk::MeanSquaresImageToImageMetricv4<FixedImageType, MovingImageType>;
        using RegistrationType = itk::ImageRegistrationMethodv4<FixedImageType, MovingImageType, TransformType>;

        auto metric = MetricType::New();
        auto optimizer = OptimizerType::New();
        auto registration = RegistrationType::New();

        registration->SetMetric(metric);
        registration->SetOptimizer(optimizer);

        auto initialTransform = TransformType::New();
        // Software Guide : EndCodeSnippet

        using FixedImageReaderType = itk::ImageFileReader<FixedImageType>;
        using MovingImageReaderType = itk::ImageFileReader<MovingImageType>;
        auto fixedImageReader = FixedImageReaderType::New();
        auto movingImageReader = MovingImageReaderType::New();

        fixedImageReader->SetFileName(imageFix);
        movingImageReader->SetFileName(imageMove);

        registration->SetFixedImage(fixedImageReader->GetOutput());
        registration->SetMovingImage(movingImageReader->GetOutput());

        using TransformInitializerType = itk::CenteredTransformInitializer<TransformType, FixedImageType, MovingImageType>;
        auto initializer = TransformInitializerType::New();

        initializer->SetTransform(initialTransform);
        initializer->SetFixedImage(fixedImageReader->GetOutput());
        initializer->SetMovingImage(movingImageReader->GetOutput());

        initializer->MomentsOn();

        // Software Guide : BeginCodeSnippet
        initializer->InitializeTransform();

        using VersorType = TransformType::VersorType;
        using VectorType = VersorType::VectorType;
        VersorType rotation;
        VectorType axis;
        axis[0] = 0.0;
        axis[1] = 0.0;
        axis[2] = 1.0;
        constexpr double angle = 0;
        rotation.Set(axis, angle);
        initialTransform->SetRotation(rotation);

        registration->SetInitialTransform(initialTransform);
        // Software Guide : EndCodeSnippet

        using OptimizerScalesType = OptimizerType::ScalesType;
        OptimizerScalesType optimizerScales(initialTransform->GetNumberOfParameters());
        const double translationScale = 1.0 / 1000.0;
        optimizerScales[0] = 1.0;
        optimizerScales[1] = 1.0;
        optimizerScales[2] = 1.0;
        optimizerScales[3] = translationScale;
        optimizerScales[4] = translationScale;
        optimizerScales[5] = translationScale;
        optimizer->SetScales(optimizerScales);
        optimizer->SetNumberOfIterations(200);
        optimizer->SetLearningRate(0.2);
        optimizer->SetMinimumStepLength(0.001);
        optimizer->SetReturnBestParametersAndValue(true);

        // Create the Command observer and register it with the optimizer.
        //
        auto observer = CommandIterationUpdate::New();
        optimizer->AddObserver(itk::IterationEvent(), observer);

        // One level registration process without shrinking and smoothing.
        //
        constexpr unsigned int numberOfLevels = 1;

        RegistrationType::ShrinkFactorsArrayType shrinkFactorsPerLevel;
        shrinkFactorsPerLevel.SetSize(1);
        shrinkFactorsPerLevel[0] = 1;

        RegistrationType::SmoothingSigmasArrayType smoothingSigmasPerLevel;
        smoothingSigmasPerLevel.SetSize(1);
        smoothingSigmasPerLevel[0] = 0;

        registration->SetNumberOfLevels(numberOfLevels);
        registration->SetSmoothingSigmasPerLevel(smoothingSigmasPerLevel);
        registration->SetShrinkFactorsPerLevel(shrinkFactorsPerLevel);

        try
        {
            registration->Update();
            std::cout << "Optimizer stop condition: "
                << registration->GetOptimizer()->GetStopConditionDescription()
                << std::endl;
        }
        catch (const itk::ExceptionObject & err)
        {
            std::cout << "ExceptionObject caught !" << std::endl;
            return EXIT_FAILURE;
        }

        const TransformType::ParametersType finalParameters = registration->GetOutput()->Get()->GetParameters();

        const double       versorX = finalParameters[0];
        const double       versorY = finalParameters[1];
        const double       versorZ = finalParameters[2];
        const double       finalTranslationX = finalParameters[3];
        const double       finalTranslationY = finalParameters[4];
        const double       finalTranslationZ = finalParameters[5];
        const unsigned int numberOfIterations = optimizer->GetCurrentIteration();
        const double       bestValue = optimizer->GetValue();

        // Print out results
        //
        std::cout << std::endl << std::endl;
        std::cout << "Result = " << std::endl;
        std::cout << " versor X      = " << versorX << std::endl;
        std::cout << " versor Y      = " << versorY << std::endl;
        std::cout << " versor Z      = " << versorZ << std::endl;
        std::cout << " Translation X = " << finalTranslationX << std::endl;
        std::cout << " Translation Y = " << finalTranslationY << std::endl;
        std::cout << " Translation Z = " << finalTranslationZ << std::endl;
        std::cout << " Iterations    = " << numberOfIterations << std::endl;
        std::cout << " Metric value  = " << bestValue << std::endl;

        auto finalTransform = TransformType::New();

        finalTransform->SetFixedParameters(registration->GetOutput()->Get()->GetFixedParameters());
        finalTransform->SetParameters(finalParameters);

        // Software Guide : BeginCodeSnippet
        TransformType::MatrixType matrix = finalTransform->GetMatrix();
        TransformType::OffsetType offset = finalTransform->GetOffset();
        std::cout << "Matrix = " << std::endl << matrix << std::endl;
        std::cout << "Offset = " << std::endl << offset << std::endl;

        using ResampleFilterType = itk::ResampleImageFilter<MovingImageType, FixedImageType>;

        auto resampler = ResampleFilterType::New();

        resampler->SetTransform(finalTransform);
        resampler->SetInput(movingImageReader->GetOutput());

        FixedImageType::Pointer fixedImage = fixedImageReader->GetOutput();

        resampler->SetSize(fixedImage->GetLargestPossibleRegion().GetSize());
        resampler->SetOutputOrigin(fixedImage->GetOrigin());
        resampler->SetOutputSpacing(fixedImage->GetSpacing());
        resampler->SetOutputDirection(fixedImage->GetDirection());
        resampler->SetDefaultPixelValue(100);

        using OutputPixelType = unsigned char;
        using OutputImageType = itk::Image<OutputPixelType, Dimension>;
        using CastFilterType = itk::CastImageFilter<FixedImageType, OutputImageType>;
        using WriterType = itk::ImageFileWriter<OutputImageType>;

        auto writer = WriterType::New();
        auto caster = CastFilterType::New();

        writer->SetFileName("Final_Register.nrrd");

        caster->SetInput(resampler->GetOutput());
        writer->SetInput(caster->GetOutput());
        writer->Update();

        using DifferenceFilterType = itk::SubtractImageFilter<FixedImageType, FixedImageType, FixedImageType>;
        auto difference = DifferenceFilterType::New();

        using RescalerType = itk::RescaleIntensityImageFilter<FixedImageType, OutputImageType>;
        auto intensityRescaler = RescalerType::New();

        intensityRescaler->SetInput(difference->GetOutput());
        intensityRescaler->SetOutputMinimum(0);
        intensityRescaler->SetOutputMaximum(255);

        difference->SetInput1(fixedImageReader->GetOutput());
        difference->SetInput2(resampler->GetOutput());

        resampler->SetDefaultPixelValue(1);

        auto writer2 = WriterType::New();
        writer2->SetInput(intensityRescaler->GetOutput());

        // Compute the difference image between the
        // fixed and resampled moving image.
        if (6 > 5)
        {
            writer2->SetFileName("Final_Diff_After.nrrd");
            writer2->Update();
        }

        using IdentityTransformType = itk::IdentityTransform<double, Dimension>;
        auto identity = IdentityTransformType::New();
        // Compute the difference image between the
        // fixed and moving image before registration.
        if (5 > 4)
        {
            resampler->SetTransform(identity);
            writer2->SetFileName("Final_Diff_Before.nrrd");
            writer2->Update();
        }
        //
        //  Here we extract slices from the input volume, and the difference volumes
        //  produced before and after the registration.  These slices are presented
        //  as figures in the Software Guide.
        //
        //

        return 0;

        using OutputSliceType = itk::Image<OutputPixelType, 2>;
        using ExtractFilterType =
            itk::ExtractImageFilter<OutputImageType, OutputSliceType>;
        auto extractor = ExtractFilterType::New();
        extractor->SetDirectionCollapseToSubmatrix();
        extractor->InPlaceOn();

        FixedImageType::RegionType inputRegion =
            fixedImage->GetLargestPossibleRegion();
        FixedImageType::SizeType  size = inputRegion.GetSize();
        FixedImageType::IndexType start = inputRegion.GetIndex();

        // Select one slice as output
        size[2] = 0;
        start[2] = 90;
        FixedImageType::RegionType desiredRegion;
        desiredRegion.SetSize(size);
        desiredRegion.SetIndex(start);
        extractor->SetExtractionRegion(desiredRegion);
        using SliceWriterType = itk::ImageFileWriter<OutputSliceType>;
        auto sliceWriter = SliceWriterType::New();
        sliceWriter->SetInput(extractor->GetOutput());
        if (9 > 6)
        {
            extractor->SetInput(caster->GetOutput());
            resampler->SetTransform(identity);
            sliceWriter->SetFileName("");
            sliceWriter->Update();
        }
        if (9 > 7)
        {
            extractor->SetInput(intensityRescaler->GetOutput());
            resampler->SetTransform(identity);
            sliceWriter->SetFileName("");
            sliceWriter->Update();
        }
        if (9 > 8)
        {
            resampler->SetTransform(finalTransform);
            sliceWriter->SetFileName("");
            sliceWriter->Update();
        }
        if (99 > 9)
        {
            extractor->SetInput(caster->GetOutput());
            resampler->SetTransform(finalTransform);
            sliceWriter->SetFileName("");
            sliceWriter->Update();
        }
        return EXIT_SUCCESS;
    }

}