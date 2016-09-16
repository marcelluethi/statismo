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
#include "itkASMNormalDirectionPointSampler.h"
#include "itkASMNormalDirectionFeatureExtractor.h"
#include "itkASMGaussianGradientImagePreprocessor.h"
#include "itkActiveShapeModel.h"
#include "itkASMFitting.h"
#include "MHFitting.h"
#include "itkMHFitting.h"
#include "itkMeshFileWriter.h"
#include "itkTimeProbe.h"
#include "itkMeshFileReader.h"
#include "../cli/utils/statismo-fitting-utils.h"
#include "itkReducedVarianceModelBuilder.h"
//#include "itkRigidTransformModelBuilder.h"
#include "itkPosteriorModelBuilder.h"


typedef itk::Mesh<float, 3> MeshType;
typedef itk::Image<float, 3> ImageType;
typedef itk::Image<short, 3> ShortImageType;
typedef itk::ActiveShapeModel<MeshType, ImageType> ActiveShapeModelType;
typedef itk::ImageFileReader<ImageType> ImageReaderType;
typedef itk::ImageFileReader<ShortImageType> ShortImageReaderType;

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


int main(int argc, char *argv[]) {

    StatismoUI::StatismoUI ui;
    StatismoUI::Group modelgroup = ui.createGroup("model");
    StatismoUI::Group targetgroup = ui.createGroup("target");



    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    //std::string pcaModelName("/export/skulls/data/shapes/ulna-right/aligned/registered-pami-ams/model-asm/asm-pca-3.h5");
    std::string flexibleModelName("/tmp/fancyasm5.h5");
    //std::string pcaModelName("//home/marcel/data/ulna-right/test/asm-pca-3.h5");(
    std::string pcaModelName("//tmp/asm-pca.h5");

    ActiveShapeModelType::Pointer pcaAsmModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    pcaAsmModel->Load(representer,  pcaModelName.c_str());

    ActiveShapeModelType::Pointer flexibleModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representerFlexibleModel = RepresenterType::New();
    flexibleModel->Load(representerFlexibleModel,  flexibleModelName.c_str());


    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);


    statismo::ASMFittingConfiguration asmFitConfig(3,5,3);
    statismo::MHFittingConfiguration mhFitConfig(asmFitConfig);

    // read and preprocess image
    ImageReaderType::Pointer reader = ImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    //reader->SetFileName("/export/skulls/data/shapes/ulna-right/aligned/initial/volume-ct/downsampled-2/vsd-0.nii");
    //reader->SetFileName("/export/skulls/data/shapes/esophagus/raw/normalized-varian/volume-ct/varian-0021.nii");
    //reader->SetFileName("/export/skulls/data/shapes/esophagus/raw/normalized-varian/volume-ct/varian-0001.nii");
    reader->SetFileName("/home/luetma00/Download/LUCRUSH02.vtk");
    //reader->SetFileName("/home/luetma00/Download/LUCRU2.vtk");
    //reader->SetFileName("/tmp/lucrush03-0.vtk");
////    reader->SetFileName("//home/marcel/data/ulna-right/test/image.nii");

    reader->Update();
    ImageType::Pointer image = reader->GetOutput();


    ShortImageReaderType::Pointer shortreader = ShortImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    //reader->SetFileName("/export/skulls/data/shapes/ulna-right/aligned/initial/volume-ct/downsampled-2/vsd-0.nii");
    //shortreader->SetFileName("/export/skulls/data/shapes/esophagus/raw/normalized-varian/volume-ct/varian-0021.nii");
    //shortreader->SetFileName("/export/skulls/data/shapes/esophagus/raw/normalized-varian/volume-ct/varian-0001.nii");
    shortreader->SetFileName("/home/luetma00/Download/LUCRUSH02.vtk");
    //shortreader->SetFileName("/tmp/lucrush03-0.nii");
    shortreader->Update();
    ShortImageType::Pointer shortImage = shortreader->GetOutput();





    statismo::ASMPreprocessedImage<MeshType> *pimage = pcaAsmModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    // just for testing
    itk::ReducedVarianceModelBuilder<MeshType>::Pointer redModelBuilder = itk::ReducedVarianceModelBuilder<MeshType>::New();
    

    std::vector<PointType> linePoints;
    //linePoints = readLandmarksFile<MeshType>(std::string("/tmp/lucrush2-lms.csv"));
    //linePoints = readLandmarksFile<MeshType>(std::string("/tmp/varian-0001-line-lms.csv"));
    //linePoints = readLandmarksFile<MeshType>(std::string("/tmp/lucrush3-line-lms.csv"));
    linePoints = readLandmarksFile<MeshType>(std::string("/home/luetma00/Download/LUCRUSH02-line-lms.csv"));
    //linePoints = readLandmarksFile<MeshType>(std::string("/tmp/0021lms-line.csv"));
    // Get poitn ids of reference and target points
    std::vector<PointType> refPoints;
    refPoints = readLandmarksFile<MeshType>(std::string("/tmp/fancylms.csv"));

    std::vector<PointType> targetPoints;
    //targetPoints = readLandmarksFile<MeshType>(std::string("/tmp/lucrush3-lms.csv"));
    targetPoints = readLandmarksFile<MeshType>(std::string("/home/luetma00/Download/LUCRUSH02-lms.csv"));
    //targetPoints = readLandmarksFile<MeshType>(std::string("//tmp/varian-0001-lm.csv"));
    //targetPoints = readLandmarksFile<MeshType>(std::string("/tmp/0021lms.csv"));





    // You should use here
    // > representer->ComputeRigidTransformFromLandmarks
    // if you have landmarks
    // Currently I only want an identity transform
    itk::VersorRigid3DTransform<float>::Pointer currentTransform = itk::VersorRigid3DTransform<float>::New();
    //k::Euler3DTransform<float>::Pointer currentTransform = itk::Euler3DTransform<float>::New();
    //RigidTransformType::Pointer currentTransform(versorTransform.GetPointer());
    currentTransform->SetIdentity();


    PointType centerModelMean = computeCenterOfMass(pcaAsmModel->GetStatisticalModel()->DrawMean());
    PointType centerReference = computeCenterOfMass(refPoints);
    PointType centerTarget =computeCenterOfMass(targetPoints);

    std::list<PointType> lr; lr.push_back(centerReference);
    ui.showPointCloud(modelgroup, lr, "refCenter");

    std::list<PointType> lt; lt.push_back(centerTarget);
    ui.showPointCloud(targetgroup, lt, "targetCenter");

    std::list<PointType> lm; lt.push_back(centerModelMean);
    ui.showPointCloud(targetgroup, lt, "meanCenter");


    itk::Vector<float> initialTranslation = currentTransform->GetTranslation();
    for (unsigned d = 0; d < 3; ++d) {
        initialTranslation.SetElement(d, centerTarget.GetElement(d) - centerReference.GetElement(d));
    }


    currentTransform->SetTranslation(initialTranslation);
    std::cout << "initial translation " << currentTransform->GetTranslation() << std::endl;
    currentTransform->SetCenter(centerModelMean);

    typedef itk::PointsLocator< typename RepresenterType::MeshType::PointsContainer > PointsLocatorType;
    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(pcaAsmModel->GetStatisticalModel()->GetRepresenter()->GetReference()->GetPoints());
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



    statismo::VectorType coeffs = statismo::VectorType::Zero(pcaAsmModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());

    // very ITK unlike, we use a init method instead of setting all fields manually.
    // This avoids 99% of all core dumps :-)
    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->init(image, pimage, correspondingPoints, linePoints, pcaAsmModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, coeffs);

    std::cout << "Initialization done." << std::endl;

	StatismoUI::ShapeModelView ssmView = ui.showStatisticalShapeModel(modelgroup, pcaAsmModel->GetStatisticalModel(), "a model");
    StatismoUI::ShapeModelTransformationView ssmTransformationView = ssmView.GetShapeModelTransformationView();


    ui.showImage(targetgroup, shortImage, "target image");

    vnl_vector<float> lastCoeffs;
    itk::OptimizerParameters<float> lastTransformParameters;
    for (int i =1; i <= 2000; ++i) {
        FittingResultType::Pointer result;

        fittingStep->NextSample();
        result = fittingStep->GetOutput();


       if (lastTransformParameters != result->GetRigidTransformParameters() || lastCoeffs != result->GetCoefficients()) {
           ssmTransformationView.SetPoseTransformation(StatismoUI::PoseTransformation(currentTransform)).SetShapeTransformation(result->GetCoefficients());
            ui.updateShapeModelTransformationView(ssmTransformationView);
            lastTransformParameters = result->GetRigidTransformParameters();
            lastCoeffs = result->GetCoefficients();
        }



        currentTransform->SetParameters(result->GetRigidTransformParameters());


        std::cout << "Writing result of iteration " << i << std::endl;



    }


    std::cout << "current transform " << currentTransform << std::endl;


    // compute the uncertainty and set the model accordingly
    PosteriorModelBuilderType::Pointer posteriorModelBuilder = PosteriorModelBuilderType::New();
    PosteriorModelBuilderType::PointValueWithCovarianceListType constraints;


    FittingStepType::Pointer fittingStepFlexibleShapeUncertainty = FittingStepType::New();
    vnl_vector<float> coeffsForFlexibleModel  = flexibleModel->GetStatisticalModel()->ComputeCoefficients(pcaAsmModel->GetStatisticalModel()->DrawSample(fittingStep->GetOutput()->GetCoefficients()));
    fittingStepFlexibleShapeUncertainty->init(image, pimage, correspondingPoints, linePoints, flexibleModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, fromVnlVector(coeffsForFlexibleModel));

        typedef std::map<unsigned, statismo::MultiVariateNormalDistribution> UncertaintyMap;
    UncertaintyMap uncertaintyMap = fittingStepFlexibleShapeUncertainty->computePointUncertainty(correspondingPoints, targetPoints);

    MeshType::Pointer ref = pcaAsmModel->GetStatisticalModel()->GetRepresenter()->GetReference();


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

        vnl_matrix<float> Rfloat = currentTransform->GetMatrix().GetVnlMatrix();
        vnl_matrix<double> R(Rfloat.rows(), Rfloat.cols());
        for (unsigned k = 0; k < R.rows(); ++k) {
            for(unsigned l = 0; l < R.cols(); l++) {
                R(k,l) = Rfloat(k,l);
            }
        }
        vnl_matrix<double> covTargetPos = R * cov * R.transpose();

        StatisticalModelType::PointValuePairType pointValue(refPt , currentTransform->GetInverseTransform()->TransformPoint(correspondingPoints[i].second));
        // we use the original target point and not the one estimated.
        StatisticalModelType::PointValueWithCovariancePairType  pointValueCov(pointValue, uncertainty);
        constraints.push_back(pointValueCov);
        //targetPointsForVisualization.push_back(currentTransform->TransformPoint(targetPoint));
        ui.showLandmark(targetgroup, correspondingPoints[i].second, covTargetPos, os.str());
    }

    StatismoUI::Group group = ui.createGroup("pts");
    ui.showPointCloud(group, targetPointsForVisualization, "target Points");


    //StatisticalModelType::Pointer posteriorModel = posteriorModelBuilder->BuildNewModelFromModel(pcaAsmModel->GetStatisticalModel(), constraints, false);
    StatisticalModelType::Pointer posteriorModel = posteriorModelBuilder->BuildNewModelFromModel(flexibleModel->GetStatisticalModel(), constraints, false);

    // we need to project the old solution in to the model


    vnl_vector<float> newCoeffs  = posteriorModel->ComputeCoefficients(pcaAsmModel->GetStatisticalModel()->DrawSample(fittingStep->GetOutput()->GetCoefficients()));


    pcaAsmModel->SetStatisticalModel(posteriorModel);

    FittingStepType::Pointer fittingStepFlexibleShapeAndHU = FittingStepType::New();
    fittingStepFlexibleShapeAndHU->init(image, pimage, correspondingPoints, linePoints, pcaAsmModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, fromVnlVector(newCoeffs));


    fittingStepFlexibleShapeAndHU->SetChainToLmAndHU(correspondingPoints, targetPoints, currentTransform, fromVnlVector(newCoeffs));

    ui.removeShapeModel(ssmView);

    StatismoUI::Group modelgroupPosterior = ui.createGroup("poster");
    StatismoUI::ShapeModelView shapeModelViewPosterior = ui.showStatisticalShapeModel(modelgroupPosterior, posteriorModel, "a model");
    StatismoUI::ShapeModelTransformationView ssmTransformationViewPosterior = shapeModelViewPosterior.GetShapeModelTransformationView();



    for (int i =1; i <= 2000; ++i) {

        std::cout << "iteration: " << i << std::endl;
        FittingResultType::Pointer result;

        fittingStepFlexibleShapeAndHU->NextSample();
        result = fittingStepFlexibleShapeAndHU->GetOutput();
        currentTransform->SetParameters(result->GetRigidTransformParameters());

        if (lastTransformParameters != result->GetRigidTransformParameters() || lastCoeffs != result->GetCoefficients()) {

            ssmTransformationViewPosterior.SetPoseTransformation(StatismoUI::PoseTransformation(currentTransform)).SetShapeTransformation(result->GetCoefficients());
            ui.updateShapeModelTransformationView(ssmTransformationViewPosterior);

            lastTransformParameters = result->GetRigidTransformParameters();
            lastCoeffs = result->GetCoefficients();
        }



    }


    ui.showTriangleMesh(targetgroup, fittingStepFlexibleShapeAndHU->GetOutput()->GetMesh(), "final Result");





    return 0;
}

