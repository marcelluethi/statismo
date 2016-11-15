#ifndef STATISMO_ITKASMIDENTITYIMAGEPREPROCESSOR_H
#define STATISMO_ITKASMIDENTITYIMAGEPREPROCESSOR_H

#include "ASMImagePreprocessor.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "HDF5Utils.h"

namespace itk {

    template <typename TPointSet, typename TImage>
    class ASMIdentityPreprocessedImage: public statismo::ASMPreprocessedImage<TPointSet> {
    private:
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;
        //typedef typename InterpolatedImageType::CovariantVectorType CovariantVectorType;
        typedef CovariantVector<float, TImage::ImageDimension> CovariantVectorType;
        typedef vnl_vector<statismo::ScalarType> VectorType;


        const TImage* m_inputImage;
        const InterpolatedImagePointerType m_interpolatedImage;

        ASMIdentityPreprocessedImage(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage): m_inputImage(inputImage), m_interpolatedImage(interpolatedImage) {}

        static statismo::VectorType fromVnlVector(const VectorType& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());

        }
    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        virtual bool IsDefined(const PointType& point) const {
            typename TImage::IndexType index;
            // we don't care about the actual index, just whether it's present
            return m_inputImage->TransformPhysicalPointToIndex(point, index);
        }

        virtual statismo::VectorType Evaluate(const PointType& point) const {
            CovariantVectorType cv = m_interpolatedImage->EvaluateDerivative(point);
            return fromVnlVector(cv.GetVnlVector());
        }

        static ASMIdentityPreprocessedImage* Create(const TImage* inputImage, const InterpolatedImagePointerType interpolatedImage) {
            return new ASMIdentityPreprocessedImage(inputImage, interpolatedImage);
        }
    };

    template <typename TPointSet, typename TImage>
    class ASMIdentityImagePreprocessor: public statismo::ASMImagePreprocessor<TPointSet, TImage> {
    private:
        typedef BSplineInterpolateImageFunction<TImage, typename TImage::PixelType, typename TImage::PixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::Pointer InterpolatedImagePointerType;


        InterpolatedImagePointerType Interpolate(const TImage* const image) const {
            InterpolatedImagePointerType inter = InterpolatedImageType::New();
            inter->SetSplineOrder(1);

            inter->SetInputImage(image);
            return inter;
        }

    public:
        typedef ASMIdentityPreprocessedImage<TPointSet, TImage> PreprocessedImplType;

        ASMIdentityImagePreprocessor() {}

        virtual ASMIdentityImagePreprocessor<TPointSet, TImage>* Clone() const {
            return new ASMIdentityImagePreprocessor();
        };

        virtual PreprocessedImplType* Preprocess(const TImage* image) const {
            return PreprocessedImplType::Create(image, Interpolate(image));
        };
    };

    template<typename TPointSet, typename TImage>
    class ASMIdentityImagePreprocessorFactory : public statismo::ASMImagePreprocessorFactory<TPointSet, TImage> {

        typedef ASMIdentityImagePreprocessor<TPointSet, TImage> InstanceType;

    public:

        static const ASMIdentityImagePreprocessorFactory *GetInstance() {
            static ASMIdentityImagePreprocessorFactory *instance = new ASMIdentityImagePreprocessorFactory();
            return instance;
        }

        virtual std::string GetDescriptor() const {
            return "builtin::Identity";
        }

        virtual const statismo::ASMImagePreprocessor<TPointSet, TImage> *Instantiate(
                const H5::Group &h5Group) const {

            InstanceType* instance = new InstanceType();
            return instance;
        }
    };


}
#endif //STATISMO_ITKASMGAUSSIANGRADIENTIMAGEPREPROCESSOR_H
