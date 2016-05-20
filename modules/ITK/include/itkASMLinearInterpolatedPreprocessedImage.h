#ifndef STATISMO_ITKASMLinearInterpolateIMAGEPREPROCESSOR_H
#define STATISMO_ITKASMLinearInterpolateIMAGEPREPROCESSOR_H

#include "ASMImagePreprocessor.h"
#include "itkLinearInterpolateImageFunction.h"

namespace itk {

    template <typename TPointSet, typename TImage>
    class ASMLinearInterpolatedPreprocessedImage: public statismo::ASMPreprocessedImage<TPointSet> {
    private:
      typedef TImage ImageType;
      typedef LinearInterpolateImageFunction<ImageType> InterpolatorType;
      typedef typename InterpolatorType::Pointer InterpolatorPointer;

      const ImageType* m_inputImage;
      InterpolatorPointer m_interpolator;

      ASMLinearInterpolatedPreprocessedImage(const ImageType* inputImage) : m_inputImage(inputImage)
      {
        m_interpolator = InterpolatorType::New();
        m_interpolator->SetInputImage(m_inputImage);
      }


    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        virtual bool IsDefined(const PointType& point) const {

          return m_interpolator->IsInsideBuffer(point);

        }

        virtual statismo::VectorType Evaluate(const PointType& point) const {
          auto value = m_interpolator->Evaluate(point);
          statismo::VectorType outputVector(1);
          outputVector[0] = value;
          return outputVector;
        }

        static ASMLinearInterpolatedPreprocessedImage* Create(const TImage* inputImage) {
            return new ASMLinearInterpolatedPreprocessedImage(inputImage);
        }
    };


}
#endif //STATISMO_ITKASMLinearInterpolateIMAGEPREPROCESSOR_H
