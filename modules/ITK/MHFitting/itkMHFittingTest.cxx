#include "StatismoUI.h"
  #include <iostream>
#include <itkMesh.h>
#include <itkImage.h>
#include <itkImageFileReader.h>
#include "itkStatismoIO.h"
#include <itkStandardMeshRepresenter.h>
#include "ASMFitting.h"
#include <itkEuler3DTransform.h>
#include <itkCenteredTransformInitializer.h>
#include <itkCenteredVersorTransformInitializer.h>
#include <itkPosteriorModelBuilder.h>
#include <itkCastImageFilter.h>
#include "itkASMNormalDirectionPointSampler.h"
#include "itkASMNormalDirectionFeatureExtractor.h"
#include "itkASMGaussianGradientImagePreprocessor.h"
#include "itkActiveShapeModel.h"
#include "itkASMFitting.h"
#include "itkMeshFileWriter.h"
#include "itkTimeProbe.h"
#include "itkMeshFileReader.h"
#include "../cli/utils/statismo-fitting-utils.h"
#include "itkReducedVarianceModelBuilder.h"
//#include "itkRigidTransformModelBuilder.h"
#include "itkPosteriorModelBuilder.h"
#include "itkASMIdentityImagePreprocessor.h"
#include "itkMHFitting.h"

typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::Image<short, 3> ShortImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;
typedef statismo::ASMPreprocessedImage<MeshType> PreprocessedImageType;
typedef itk::ASMNormalDirectionPointSampler<MeshType, ImageType> SamplerType;

typedef itk::StandardMeshRepresenter<float, 3> RepresenterType;
typedef RepresenterType::RigidTransformType RigidTransformType;
typedef RepresenterType::PointType PointType;
typedef itk::StatisticalModel<MeshType> StatisticalModelType;
typedef itk::Euler3DTransform< float > TransformType;

typedef itk::MHFittingStepper<MeshType, ImageType> FittingStepType;
typedef itk::MHFittingResult<MeshType, ImageType> FittingResultType;
typedef itk::PosteriorModelBuilder<MeshType> PosteriorModelBuilderType;


// FIXME: these conversions have to go.
typedef vnl_vector<statismo::ScalarType> VnlVectorType;
VnlVectorType toVnlVector(const statismo::VectorType& v) {
    return VnlVectorType(v.data(), v.rows());
}
statismo::VectorType fromVnlVector(const VnlVectorType& v) {
    return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

}

PointType computeCenterOfMass(const MeshType* mesh) {

    PointType centerOfMass;
    centerOfMass.Fill(0);

    for (unsigned i = 0; i < mesh->GetNumberOfPoints(); ++i) {
        for (unsigned d = 0; d < 3; ++d) {
            centerOfMass.SetElement(d, centerOfMass.GetElement(d) + mesh->GetPoint(i).GetElement(d));
        }
    }
    for (unsigned d = 0; d < 3; ++d) {
        centerOfMass[d] /= mesh->GetNumberOfPoints();
    }



    return centerOfMass;
}

PointType computeCenterOfMass(const std::vector<PointType>& points) {

    PointType centerOfMass;
    centerOfMass.Fill(0);

    for (std::vector<PointType>::const_iterator it = points.begin(); it != points.end(); ++it) {
        for (unsigned d = 0; d < 3; ++d) {
            centerOfMass.SetElement(d, centerOfMass.GetElement(d) + it->GetElement(d));
        }
    }
    for (unsigned d = 0; d < 3; ++d) {
        centerOfMass[d] /= points.size();
    }



    return centerOfMass;
}


ActiveShapeModelType::Pointer loadASM(const std::string& filename) {
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    //statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMIdentityImagePreprocessorFactory<MeshType, ImageType>::GetInstance());
    ActiveShapeModelType::Pointer asmodel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    asmodel->Load(representer,  filename.c_str());

    return asmodel;
}

ImageType::Pointer readFloatImage(const std::string& filename) {
    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName(filename);
    reader->Update();
    ImageType::Pointer image = reader->GetOutput();
    return image;
}

ShortImageType::Pointer castImageToShort(ImageType* image) {
    itk::CastImageFilter<ImageType, ShortImageType>::Pointer caster = itk::CastImageFilter<ImageType, ShortImageType>::New();
    caster->SetInput(image);
    caster->Update();
    ShortImageType::Pointer castedImage = caster->GetOutput();
    return castedImage;
}



itk::VersorRigid3DTransform<float>::Pointer setupInitialTransform(const StatisticalModelType* model, const std::vector<PointType>& refPoints, const std::vector<PointType>& targetPoints) {

    itk::VersorRigid3DTransform<float>::Pointer transform = itk::VersorRigid3DTransform<float>::New();

    transform->SetIdentity();


    // compute the initial transformation

    PointType centerReference = computeCenterOfMass(refPoints);
    PointType centerTarget = computeCenterOfMass(targetPoints);
    itk::Vector<float> initialTranslation = transform->GetTranslation();
    for (unsigned d = 0; d < 3; ++d) {
        initialTranslation.SetElement(d, centerTarget.GetElement(d) - centerReference.GetElement(d));
    }


    PointType centerModelMean = computeCenterOfMass(model->DrawMean());

    transform->SetTranslation(initialTranslation);
    transform->SetCenter(centerModelMean);
    return transform;
}

FittingStepType::CorrespondencePoints  getCorrespondingPointsToReference( MeshType* reference,  const std::vector<PointType>& refPoints, const std::vector<PointType>& targetPoints) {

    typedef itk::PointsLocator< typename RepresenterType::MeshType::PointsContainer > PointsLocatorType;
    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(reference->GetPoints());
    ptLocator->Initialize();

    if (refPoints.size() != targetPoints.size()) {
        std::cout << "need the same number of reference and target points" << std::endl;
        exit(-1);
    }

    FittingStepType::CorrespondencePoints correspondingPoints;
    for (unsigned i = 0; i < refPoints.size(); ++i) {
        long ptId = ptLocator->FindClosestPoint(refPoints[i]);
        correspondingPoints.push_back(std::make_pair(ptId, targetPoints[i]));
    }
    return correspondingPoints;
}



FittingResultType::Pointer fitPose( StatismoUI::StatismoUI& ui,
              StatismoUI::ShapeModelTransformationView& ssmTransformationView,
             ImageType* image,
             PreprocessedImageType* preprocessedImage,
             const FittingStepType::CorrespondencePoints& correspondingPoints,
             const std::vector<PointType>& linePoints,
             ActiveShapeModelType* model,
             SamplerType* fitSampler,
                           mhfitting::MHFittingConfiguration  configuration,
             itk::VersorRigid3DTransform<float>* transform,
             statismo::VectorType coeffs) {

    itk::VersorRigid3DTransform<float>::Pointer currentRigidTransform = transform->Clone();

        // very ITK unlike, we use a init method instead of setting all fields manually.
    // This avoids 99% of all core dumps :-)
    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->init(image, preprocessedImage, correspondingPoints, linePoints, model, FittingStepType::SamplerPointerType(fitSampler), configuration, currentRigidTransform, coeffs);
    fittingStep->SetChainToPoseOnly(correspondingPoints, linePoints, transform, coeffs);
    std::cout << "Initialization done." << std::endl;

    vnl_vector<float> lastCoeffs;
    itk::OptimizerParameters<float> lastTransformParameters;
    for (int i =1; i <= 100; ++i) {
        FittingResultType::Pointer result;

        fittingStep->NextSample();
        result = fittingStep->GetOutput();


        if (lastTransformParameters != result->GetRigidTransformParameters() || lastCoeffs != result->GetCoefficients()) {
            ssmTransformationView.SetPoseTransformation(StatismoUI::PoseTransformation(currentRigidTransform)).SetShapeTransformation(result->GetCoefficients());
            ui.updateShapeModelTransformationView(ssmTransformationView);
            lastTransformParameters = result->GetRigidTransformParameters();
            lastCoeffs = result->GetCoefficients();
        }

        currentRigidTransform->SetParameters(result->GetRigidTransformParameters());

        std::cout << "Done with iteration " << i << std::endl;
    }


    std::cout << "current transform " << currentRigidTransform << std::endl;
    FittingResultType::Pointer result = fittingStep->GetOutput();
    return result;
}




FittingResultType::Pointer fitPoseAndShape( StatismoUI::StatismoUI& ui,
                                    StatismoUI::ShapeModelTransformationView& ssmTransformationView,
                                    ImageType* image,
                                    PreprocessedImageType* preprocessedImage,
                                    const FittingStepType::CorrespondencePoints& correspondingPoints,
                                    const std::vector<PointType>& linePoints,
                                    ActiveShapeModelType* model,
                                    SamplerType* fitSampler,
                                    mhfitting::MHFittingConfiguration  configuration,
                                    itk::VersorRigid3DTransform<float>* transform,
                                    statismo::VectorType coeffs) {

    itk::VersorRigid3DTransform<float>::Pointer currentRigidTransform = transform->Clone();

    // very ITK unlike, we use a init method instead of setting all fields manually.
    // This avoids 99% of all core dumps :-)
    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->init(image, preprocessedImage, correspondingPoints, linePoints, model, FittingStepType::SamplerPointerType(fitSampler), configuration, currentRigidTransform, coeffs);
    fittingStep->SetChainToPoseAndShape(correspondingPoints, linePoints, transform, coeffs);


    vnl_vector<float> lastCoeffs;
    itk::OptimizerParameters<float> lastTransformParameters;
    for (int i =1; i <= 100; ++i) {
        FittingResultType::Pointer result;

        fittingStep->NextSample();
        result = fittingStep->GetOutput();


        if (lastTransformParameters != result->GetRigidTransformParameters() || lastCoeffs != result->GetCoefficients()) {
            ssmTransformationView.SetPoseTransformation(StatismoUI::PoseTransformation(currentRigidTransform)).SetShapeTransformation(result->GetCoefficients());
            ui.updateShapeModelTransformationView(ssmTransformationView);
            lastTransformParameters = result->GetRigidTransformParameters();
            lastCoeffs = result->GetCoefficients();
        }

        currentRigidTransform->SetParameters(result->GetRigidTransformParameters());

        std::cout << "Done with iteration " << i << std::endl;
    }


    std::cout << "current transform " << currentRigidTransform << std::endl;
    FittingResultType::Pointer result = fittingStep->GetOutput();
    return result;
}


FittingResultType::Pointer fitWithHU( StatismoUI::StatismoUI& ui,
                                    StatismoUI::ShapeModelTransformationView& ssmTransformationView,
                                    ImageType* image,
                                    PreprocessedImageType* preprocessedImage,
                                    const FittingStepType::CorrespondencePoints& correspondingPoints,
                                    const std::vector<PointType>& linePoints,
                                    ActiveShapeModelType* model,
                                    SamplerType* fitSampler,
                                    mhfitting::MHFittingConfiguration  configuration,
                                    itk::VersorRigid3DTransform<float>* transform,
                                    vnl_vector<float> coeffs) {

    itk::VersorRigid3DTransform<float>::Pointer currentTransform = transform->Clone();

    FittingStepType::Pointer fittingStepFlexibleShapeAndHU = FittingStepType::New();
    fittingStepFlexibleShapeAndHU->init(image, preprocessedImage, correspondingPoints, linePoints, model, FittingStepType::SamplerPointerType(fitSampler), configuration, transform, fromVnlVector(coeffs));
    fittingStepFlexibleShapeAndHU->SetChainToLmAndHU(correspondingPoints, linePoints, transform, fromVnlVector(coeffs));


    vnl_vector<float> lastCoeffs;
    itk::OptimizerParameters<float> lastTransformParameters;

    for (int i =1; i <= 2000; ++i) {

        std::cout << "iteration: " << i << std::endl;
        FittingResultType::Pointer result;

        fittingStepFlexibleShapeAndHU->NextSample();
        result = fittingStepFlexibleShapeAndHU->GetOutput();
        currentTransform->SetParameters(result->GetRigidTransformParameters());

        if (lastTransformParameters != result->GetRigidTransformParameters() || lastCoeffs != result->GetCoefficients()) {

            ssmTransformationView.SetPoseTransformation(StatismoUI::PoseTransformation(currentTransform)).SetShapeTransformation(result->GetCoefficients());
            ui.updateShapeModelTransformationView(ssmTransformationView);

            lastTransformParameters = result->GetRigidTransformParameters();
            lastCoeffs = result->GetCoefficients();
        }



    }

    FittingResultType::Pointer result = fittingStepFlexibleShapeAndHU->GetOutput();
    return result;
}



void computePosteriorModel(
        ImageType* image,
        PreprocessedImageType* preprocessedImage,
        const FittingStepType::CorrespondencePoints& correspondingPoints,
        const std::vector<PointType>& linePoints,
        ActiveShapeModelType* poseModel,
        vnl_vector<float> lastModelCoefficientsPoseModel,
        ActiveShapeModelType* flexibleModel,
        SamplerType* fitSampler,
        mhfitting::MHFittingConfiguration  configuration,
        itk::VersorRigid3DTransform<float>* initialTransform,
        //statismo::VectorType coeffs,
        ActiveShapeModelType* resultingPosteriorModel,
        vnl_vector<float>& posteriorModelCoefficients) {

    // compute the uncertainty and set the model accordingly
    PosteriorModelBuilderType::Pointer posteriorModelBuilder = PosteriorModelBuilderType::New();
    PosteriorModelBuilderType::PointValueWithCovarianceListType constraints;

    itk::VersorRigid3DTransform<float>::Pointer currentTransform = initialTransform->Clone();

    FittingStepType::Pointer fittingStepFlexibleShapeUncertainty = FittingStepType::New();
    vnl_vector<float> coeffsForFlexibleModel  = flexibleModel->GetStatisticalModel()->ComputeCoefficients(poseModel->GetStatisticalModel()->DrawSample(lastModelCoefficientsPoseModel));
    fittingStepFlexibleShapeUncertainty->init(image, preprocessedImage, correspondingPoints, linePoints, flexibleModel, FittingStepType::SamplerPointerType(fitSampler), configuration, currentTransform, fromVnlVector(coeffsForFlexibleModel));

    typedef std::map<unsigned, statismo::MultiVariateNormalDistribution> UncertaintyMap;
    UncertaintyMap uncertaintyMap = fittingStepFlexibleShapeUncertainty->computePointUncertainty(correspondingPoints, linePoints);

    MeshType::Pointer ref = flexibleModel->GetStatisticalModel()->GetRepresenter()->GetReference();


    std::list<PointType> targetPointsForVisualization;
    for (unsigned i = 0; i < correspondingPoints.size(); ++i) {
        unsigned id = correspondingPoints[i].first;

        PointType refPt = ref->GetPoint(id);
        statismo::MatrixType uncertainty = uncertaintyMap.find(id)->second.covariance;
        vnl_matrix<double> cov(3, 3); cov.set_identity();
        for (unsigned k = 0; k < 3; ++k) {
            for (unsigned l = 0; l < 3; ++l) {
                cov(k,l) = uncertainty(k,l);
            }
        }
        std::ostringstream os;
        os << "landmark " << i;
//        ui.showLandmark(modelgroup, targetPoint, cov, os.str());

        vnl_matrix<float> Rfloat = initialTransform->GetMatrix().GetVnlMatrix();
        vnl_matrix<double> R(Rfloat.rows(), Rfloat.cols());
        for (unsigned k = 0; k < R.rows(); ++k) {
            for(unsigned l = 0; l < R.cols(); l++) {
                R(k,l) = Rfloat(k,l);
            }
        }
        vnl_matrix<double> covTargetPos = R * cov * R.transpose();

        StatisticalModelType::PointValuePairType pointValue(refPt , initialTransform->GetInverseTransform()->TransformPoint(correspondingPoints[i].second));
        // we use the original target point and not the one estimated.
        StatisticalModelType::PointValueWithCovariancePairType  pointValueCov(pointValue, uncertainty);
        constraints.push_back(pointValueCov);
    }


    //StatisticalModelType::Pointer posteriorModel = posteriorModelBuilder->BuildNewModelFromModel(pcaAsmModel->GetStatisticalModel(), constraints, false);
    StatisticalModelType::Pointer posteriorModel = posteriorModelBuilder->BuildNewModelFromModel(flexibleModel->GetStatisticalModel(), constraints, false);

    // we need to project the old solution in to the model


    posteriorModelCoefficients  = posteriorModel->ComputeCoefficients(poseModel->GetStatisticalModel()->DrawSample(lastModelCoefficientsPoseModel));


    resultingPosteriorModel->SetStatisticalModel(posteriorModel);
}


int main(int argc, char *argv[]) {

    // parameters, file names, etc.
    std::string flexibleModelName("/home/luetma00/workspaces/projects/varian/data/asm-pca-all-data-localized.h5");
    std::string pcaModelName("/home/luetma00/workspaces/projects/varian/data/asm-pca-all-data-localized.h5");
    std::string imageFilename("/export/skulls/data/shapes/esophagus/varian/raw/normalized/volume-ct/varian-0001.nii");
    std::string targetLinePointsFilename("/tmp/varian-0001-line-lms.csv");
    std::string referenceLandmarkPointsFilename("/tmp/ref-lms.csv");
    std::string targetLandmarkPointsFilename("/tmp/varian-0001-lm.csv");

    // configuring the asm and the mh fitting.
    SamplerType::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);
    statismo::ASMFittingConfiguration asmFitConfig(3,5,3);
    mhfitting::MHFittingConfiguration mhFitConfig(asmFitConfig);


    // set up the ui
    StatismoUI::StatismoUI ui;
    StatismoUI::Group modelgroup = ui.createGroup("model");
    StatismoUI::Group targetgroup = ui.createGroup("target");


    // handling the images
    ActiveShapeModelType::Pointer flexibleModel = loadASM(flexibleModelName);
    ActiveShapeModelType::Pointer pcaAsmModel = loadASM(pcaModelName);
    ImageType::Pointer image = readFloatImage(imageFilename);
    PreprocessedImageType *preprocessedImage = pcaAsmModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);


    // loading landmarks
    std::vector<PointType> refPoints = readLandmarksFile<MeshType>(referenceLandmarkPointsFilename);
    std::vector<PointType> targetPoints = readLandmarksFile<MeshType>(targetLandmarkPointsFilename);
    std::vector<PointType> linePoints = readLandmarksFile<MeshType>(targetLinePointsFilename);


    FittingStepType::CorrespondencePoints correspondingPoints = getCorrespondingPointsToReference(pcaAsmModel->GetStatisticalModel()->GetRepresenter()->GetReference(), refPoints, targetPoints);


    StatismoUI::ShapeModelView ssmView = ui.showStatisticalShapeModel(modelgroup, pcaAsmModel->GetStatisticalModel(), "a model");
    StatismoUI::ShapeModelTransformationView ssmTransformationView = ssmView.GetShapeModelTransformationView();
    ui.showImage(targetgroup, castImageToShort(image), "target image");


    // fitting pose
    itk::VersorRigid3DTransform<float>::Pointer initialTransform = setupInitialTransform(pcaAsmModel->GetStatisticalModel(), refPoints, targetPoints);
    statismo::VectorType coeffs = statismo::VectorType::Zero(pcaAsmModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());

    FittingResultType::Pointer fitPoseResult = fitPose(ui, ssmTransformationView, image, preprocessedImage, correspondingPoints, linePoints, pcaAsmModel, fitSampler, mhFitConfig, initialTransform, coeffs);
    FittingResultType::Pointer fitPoseAndShapeResult = fitPoseAndShape(ui, ssmTransformationView, image, preprocessedImage, correspondingPoints, linePoints, pcaAsmModel, fitSampler, mhFitConfig, fitPoseResult->GetRigidTransformation(), fromVnlVector(fitPoseResult->GetCoefficients()));


    // computing the posterior

    vnl_vector<float> coeffsPosteriorModel(pcaAsmModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());
    computePosteriorModel(image, preprocessedImage, correspondingPoints, linePoints, pcaAsmModel, fitPoseAndShapeResult->GetCoefficients(), flexibleModel, fitSampler, mhFitConfig, fitPoseAndShapeResult->GetRigidTransformation(), flexibleModel, coeffsPosteriorModel);


    StatismoUI::Group modelgroupPosterior = ui.createGroup("poster");
    StatismoUI::ShapeModelView shapeModelViewPosterior = ui.showStatisticalShapeModel(modelgroupPosterior, flexibleModel->GetStatisticalModel(), "a model");
    StatismoUI::ShapeModelTransformationView ssmTransformationViewPosterior = shapeModelViewPosterior.GetShapeModelTransformationView();
    ui.removeShapeModel(ssmView);


    // fitting with intensity
    FittingResultType::Pointer fitHUResult = fitWithHU(ui, ssmTransformationViewPosterior, image, preprocessedImage, correspondingPoints, linePoints, flexibleModel, fitSampler, mhFitConfig, fitPoseResult->GetRigidTransformation(), coeffsPosteriorModel);

    ui.showTriangleMesh(targetgroup, fitHUResult->GetMesh(), "final Result");





    return 0;
}

