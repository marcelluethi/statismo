#include <itkImageFileReader.h>
#include "itkASMTraits.h"
#include "itkASMNormalDirectionGradientGaussianFeatureExtractor.h"

typedef itk::ASMTraits<float, float, 3> ASM;
typedef itk::ImageFileReader<ASM::ImageType> ImageReaderType;

int main(int argc, char *argv[]) {
    statismo::ASMFeatureExtractorFactory<ASM>::RegisterImplementation(itk::ASMNormalDirectionGradientGaussianFeatureExtractorFactory<ASM>::GetInstance());

    ImageReaderType::Pointer reader = ImageReaderType::New();
    reader->SetFileName("/home/langguth/workspaces.exported/stk.idea/bladderdemo/src/main/resources/bladder/volumes/24.vtk");
    reader->Update();
    ASM::ImagePointerType image = reader->GetOutput();

    ASM::ActiveShapeModelPointerType model = ASM::ActiveShapeModelType::New();

    std::string modelname("/home/langguth/workspaces.exported/stk.idea/bladderdemo/src/main/resources/bladder/asmModels/asmLevel-1.h5");
    model->Load(ASM::MeshRepresenterType::New(), modelname.c_str());

    ASM::MeshPointerType mesh = model->GetStatisticalModel()->DrawMean();

    itk::ASMNormalDirectionPointSampler<ASM>::Pointer sampler = itk::ASMNormalDirectionPointSampler<ASM>::New();
    //sampler->Init(90, (float)(2.0/3.0));
    sampler->Init(25, (float)(0.5));
    //sampler->Init(125, (float)(0.5));

    ASM::FitterConfigurationPointerType fitConfig = ASM::FitterConfigurationType::New();
    fitConfig->Init(5,3,3);
    //fitConfig->Init(5,5,5);

    ASM::FitterPointerType fitter = ASM::FitterType::New();
    fitter->Init(fitConfig, model, mesh, image, ASM::PointSamplerPointerType(sampler.GetPointer()));
    for (int i=0; i < 25; ++i) {
        std::cout << "iteration " << i << " start" << std::endl;
        ASM::FitterResultPointerType result = fitter->Fit();
        if (!result->IsValid()) {
            std::cout << "invalid result, aborting " <<std::endl;
            exit(0);
        }
        std::cout << "coeffs (adj)" << result->GetCoefficients() << std::endl;
        fitter = fitter->SetMesh(result->GetMesh());
    }
    return 0;
}


