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

int main(int argc, char *argv[]) {

    StatismoUI::StatismoUI ui;

    std::cout << "Initializing..." << std::endl;
    // FIXME: these should go somewhere less "intrusive"
    statismo::ASMFeatureExtractorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMNormalDirectionFeatureExtractorFactory<MeshType, ImageType>::GetInstance());
    statismo::ASMImagePreprocessorFactory<MeshType, ImageType>::RegisterImplementation(itk::ASMGaussianGradientImagePreprocessorFactory<MeshType, ImageType>::GetInstance());

    //std::string modelname("/export/skulls/data/shapes/ulna-right/aligned/registered-pami-ams/model-asm/asm-pca-3.h5");
    std::string modelname("//tmp/fancyasm3"
                                  ".h5");
    //std::string modelname("//home/marcel/data/ulna-right/test/asm-pca-3.h5");(

    ActiveShapeModelType::Pointer aModel = ActiveShapeModelType::New();
    RepresenterType::Pointer representer = RepresenterType::New();
    aModel->Load(representer,  modelname.c_str());

    itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::Pointer fitSampler = itk::ASMNormalDirectionPointSampler<MeshType, ImageType>::New();
    fitSampler->SetNumberOfPoints(25);
    fitSampler->SetPointSpacing(1);


    statismo::ASMFittingConfiguration asmFitConfig(3,5,3);
    statismo::MHFittingConfiguration mhFitConfig(asmFitConfig);

    // read and preprocess image
    ImageReaderType::Pointer reader = ImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    //reader->SetFileName("/export/skulls/data/shapes/ulna-right/aligned/initial/volume-ct/downsampled-2/vsd-0.nii");
    reader->SetFileName("/export/skulls/data/shapes/esophagus/raw/normalized-varian/volume-ct/varian-0021.nii");
//    reader->SetFileName("/home/luetma00/Download/LUCRUSH02.vtk");
    //reader->SetFileName("/home/luetma00/Download/LUCRU2.vtk");
//    reader->SetFileName("/tmp/lucrush03-0.nii");
//    reader->SetFileName("//home/marcel/data/ulna-right/test/image.nii");

    reader->Update();
    ImageType::Pointer image = reader->GetOutput();


    ShortImageReaderType::Pointer shortreader = ShortImageReaderType::New();
    //reader->SetFileName("/export/skulls/data/shapes/submandibular_gland_l/aligned/initial/volume-ct/pddca-0522c0002.nii");
    //reader->SetFileName("/export/skulls/data/shapes/ulna-right/aligned/initial/volume-ct/downsampled-2/vsd-0.nii");
    shortreader->SetFileName("/export/skulls/data/shapes/esophagus/raw/normalized-varian/volume-ct/varian-0021.nii");
    //shortreader->SetFileName("/home/luetma00/Download/LUCRUSH02.vtk");
    //shortreader->SetFileName("/tmp/lucrush03-0.nii");
    shortreader->Update();
    ShortImageType::Pointer shortImage = shortreader->GetOutput();



    // You should use here
    // > representer->ComputeRigidTransformFromLandmarks
    // if you have landmarks
    // Currently I only want an identity transform
    itk::VersorRigid3DTransform<float>::Pointer currentTransform = itk::VersorRigid3DTransform<float>::New();
    //RigidTransformType::Pointer currentTransform(versorTransform.GetPointer());
    currentTransform->SetIdentity();
    itk::Point<float, 3> center;
    for (unsigned d =0; d < 3; ++d) {
        center[d] = image->GetOrigin()[d] + image->GetLargestPossibleRegion().GetSize()[d] * image->GetSpacing()[d] * 0.5;
    }
    currentTransform->SetCenter(center);


    statismo::ASMPreprocessedImage<MeshType> *pimage = aModel->GetstatismoImplObj()->GetImagePreprocessor()->Preprocess(image);

    // just for testing
    itk::ReducedVarianceModelBuilder<MeshType>::Pointer redModelBuilder = itk::ReducedVarianceModelBuilder<MeshType>::New();
    

    std::vector<PointType> linePoints;
    //linePoints = readLandmarksFile<MeshType>(std::string("/tmp/lucrush2-lms.csv"));
    //linePoints = readLandmarksFile<MeshType>(std::string("/tmp/lucrush3-line-lms.csv"));
    //linePoints = readLandmarksFile<MeshType>(std::string("/home/luetma00/Download/LUCRUSH02-line-lms.csv"));
    linePoints = readLandmarksFile<MeshType>(std::string("/tmp/0021lms-line.csv"));
    // Get poitn ids of reference and target points
    std::vector<PointType> refPoints;
    refPoints = readLandmarksFile<MeshType>(std::string("/tmp/fancylms.csv"));

    std::vector<PointType> targetPoints;
    //targetPoints = readLandmarksFile<MeshType>(std::string("/tmp/lucrush3-lms.csv"));
    //targetPoints = readLandmarksFile<MeshType>(std::string("/home/luetma00/Download/LUCRUSH02-lms.csv"));
    targetPoints = readLandmarksFile<MeshType>(std::string("/tmp/0021lms.csv"));


    typedef itk::PointsLocator< typename RepresenterType::MeshType::PointsContainer > PointsLocatorType;
    PointsLocatorType::Pointer ptLocator = PointsLocatorType::New();
    ptLocator->SetPoints(aModel->GetStatisticalModel()->GetRepresenter()->GetReference()->GetPoints());
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



    statismo::VectorType coeffs = statismo::VectorType::Zero(aModel->GetStatisticalModel()->GetNumberOfPrincipalComponents());

    // very ITK unlike, we use a init method instead of setting all fields manually.
    // This avoids 99% of all core dumps :-)
    FittingStepType::Pointer fittingStep = FittingStepType::New();
    fittingStep->init(image, pimage, correspondingPoints, linePoints, aModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, coeffs);

    std::cout << "Initialization done." << std::endl;

    StatismoUI::Group modelgroup = ui.createGroup("model");
	StatismoUI::ShapeModelView ssmView = ui.showStatisticalShapeModel(modelgroup, aModel->GetStatisticalModel(), "a model");
    StatismoUI::ShapeModelTransformationView ssmTransformationView = ssmView.GetShapeModelTransformationView();

	StatismoUI::Group targetgroup = ui.createGroup("target");
    ui.showImage(targetgroup, shortImage, "target image");

    for (int i =1; i <= 1000; ++i) {

        std::cout << "iteration: " << i << std::endl;
        FittingResultType::Pointer result;

        fittingStep->NextSample();
        result = fittingStep->GetOutput();
        ssmTransformationView.SetPoseTransformation(StatismoUI::PoseTransformation(currentTransform)).SetShapeTransformation(result->GetCoefficients());
        ui.updateShapeModelTransformationView(ssmTransformationView);

        currentTransform->SetParameters(result->GetRigidTransformParameters());

        std::cout << "rigid params " << result->GetRigidTransformParameters() << std::endl;

        std::cout << "Writing result of iteration " << i << std::endl;



    }


    std::cout << "current transform " << currentTransform << std::endl;


    // compute the uncertainty and set the model accordingly
    PosteriorModelBuilderType::Pointer posteriorModelBuilder = PosteriorModelBuilderType::New();
    PosteriorModelBuilderType::PointValueWithCovarianceListType constraints;

    typedef std::map<unsigned, statismo::MultiVariateNormalDistribution> UncertaintyMap;
    UncertaintyMap uncertaintyMap = fittingStep->computePointUncertainty(correspondingPoints, targetPoints);

    MeshType::Pointer ref = aModel->GetStatisticalModel()->GetRepresenter()->GetReference();
    //for (UncertaintyMap::const_iterator it = uncertaintyMap.begin(); it != uncertaintyMap.end(); ++it) {



    std::list<PointType> targetPointsForVisualization;
    for (unsigned i = 0; i < correspondingPoints.size(); ++i) {
        unsigned id = correspondingPoints[i].first;
        //PointType targetPoint = correspondingPoints[i].second;
        PointType targetPoint;
        targetPoint.SetElement(0, uncertaintyMap.find(id)->second.mean[0]);
        targetPoint.SetElement(1, uncertaintyMap.find(id)->second.mean[1]);
        targetPoint.SetElement(2, uncertaintyMap.find(id)->second.mean[2]);

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
        ui.showLandmark(modelgroup, targetPoint, cov, os.str());

        StatisticalModelType::PointValuePairType pointValue(refPt ,targetPoint);
        StatisticalModelType::PointValueWithCovariancePairType  pointValueCov(pointValue, uncertainty);
        constraints.push_back(pointValueCov);
        targetPointsForVisualization.push_back(currentTransform->TransformPoint(targetPoint));
    }

    StatismoUI::Group group = ui.createGroup("pts");
    ui.showPointCloud(group, targetPointsForVisualization, "target Points");


    StatisticalModelType::Pointer posteriorModel = posteriorModelBuilder->BuildNewModelFromModel(aModel->GetStatisticalModel(), constraints, false);


    // we need to project the old solution in to the model


    vnl_vector<float> newCoeffs  = posteriorModel->ComputeCoefficients(aModel->GetStatisticalModel()->DrawSample(fittingStep->GetOutput()->GetCoefficients()));


    aModel->SetStatisticalModel(posteriorModel);

    FittingStepType::Pointer fittingStep2 = FittingStepType::New();
    fittingStep2->init(image, pimage, correspondingPoints, linePoints, aModel, FittingStepType::SamplerPointerType(fitSampler.GetPointer()), mhFitConfig, currentTransform, fromVnlVector(newCoeffs));


    fittingStep2->SetChainToLmAndHU(correspondingPoints, targetPoints, currentTransform, fromVnlVector(newCoeffs));

    ui.removeShapeModel(ssmView);

    StatismoUI::Group modelgroupPosterior = ui.createGroup("poster");
    StatismoUI::ShapeModelView shapeModelViewPosterior = ui.showStatisticalShapeModel(modelgroupPosterior, posteriorModel, "a model");
    StatismoUI::ShapeModelTransformationView ssmTransformationViewPosterior = shapeModelViewPosterior.GetShapeModelTransformationView();

    for (int i =1; i <= 2000; ++i) {

        std::cout << "iteration: " << i << std::endl;
        FittingResultType::Pointer result;

        fittingStep2->NextSample();
        result = fittingStep2->GetOutput();
        currentTransform->SetParameters(result->GetRigidTransformParameters());
        ssmTransformationViewPosterior.SetPoseTransformation(StatismoUI::PoseTransformation(currentTransform)).SetShapeTransformation(result->GetCoefficients());
        ui.updateShapeModelTransformationView(ssmTransformationViewPosterior);



    }







    return 0;
}

