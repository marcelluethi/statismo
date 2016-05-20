///*
// * This file is part of the statismo library.
// *
// * Author: Christoph Langguth (christoph.langguth@unibas.ch)
// *
// * Copyright (c) 2011-2015 University of Basel
// * All rights reserved.
// *
// * Redistribution and use in source and binary forms, with or without
// * modification, are permitted provided that the following conditions
// * are met:
// *
// * Redistributions of source code must retain the above copyright notice,
// * this list of conditions and the following disclaimer.
// *
// * Redistributions in binary form must reproduce the above copyright
// * notice, this list of conditions and the following disclaimer in the
// * documentation and/or other materials provided with the distribution.
// *
// * Neither the name of the project's author nor the names of its
// * contributors may be used to endorse or promote products derived from
// * this software without specific prior written permission.
// *
// * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
// * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
// * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// *
// */
//
#ifndef STATISMO_itkASMNormalDirectionImageValueExtractor_H
#define STATISMO_itkASMNormalDirectionImageValueExtractor_H

#include "ASMFeatureExtractor.h"
#include "itkASMFeatureExtractor.h"
#include "HDF5Utils.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkASMNormalDirectionPointSampler.h"
#include "ASMImagePreprocessor.h"


namespace itk {
    template<typename TPointSet, typename TImage>
    class ASMNormalDirectionImageValueExtractor : public statismo::ASMFeatureExtractor<TPointSet, TImage> {
    private:
        typedef TImage ImageType;
        typedef statismo::ASMPreprocessedImage<TPointSet> PreprocessedImageType;
        typedef LinearInterpolateImageFunction<ImageType> InterpolatorType;
        typedef statismo::ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
        typedef statismo::VectorType VectorType;
        typedef CovariantVector<float, TImage::ImageDimension> CovariantVectorType;
        typedef vnl_vector<statismo::ScalarType> VnlVectorType;
        typedef ASMNormalDirectionPointSampler<TPointSet, TImage> SamplerType;
        typedef typename ASMNormalDirectionPointSampler<TPointSet, TImage>::Pointer SamplerPointerType;


        const SamplerPointerType m_sampler;

        static statismo::VectorType fromVnlVector(const VnlVectorType& v) {
            return Eigen::Map<const statismo::VectorType>(v.data_block(), v.size());
        }

    public:
        typedef typename statismo::Representer<TPointSet>::PointType PointType;

        ASMNormalDirectionImageValueExtractor(SamplerPointerType sampler) : m_sampler(sampler) {
        }

        virtual ~ASMNormalDirectionImageValueExtractor() {
        }

        virtual void Delete() {
            m_sampler->Delete();
            delete this;
        }

        virtual ASMNormalDirectionImageValueExtractor* Clone() const {
            SamplerPointerType sampler = SamplerType::New();
            sampler->SetNumberOfPoints(m_sampler->GetNumberOfPoints());
            sampler->SetPointSpacing(m_sampler->GetPointSpacing());
            return new ASMNormalDirectionImageValueExtractor(sampler);
        }

        virtual ASMNormalDirectionImageValueExtractor* CloneForTarget(const ActiveShapeModelType* const model, const VectorType& coefficients, const RigidTransformPointerType transform) const {
            SamplerPointerType copySampler = m_sampler->CloneForTarget(model, coefficients, transform);
            return new ASMNormalDirectionImageValueExtractor(copySampler);
        }


        virtual bool ExtractFeatures(statismo::VectorType& output, const PreprocessedImageType* const image, const PointType& point) const {
          
          statismo::VectorType features(m_sampler->GetNumberOfPoints());

          int i = 0;
          float sum = 0;
           
          std::vector<PointType> samples = m_sampler->Sample(point);
          for (typename std::vector<PointType>::iterator it = samples.begin(); it != samples.end(); ++it) {
            PointType sample = *it;

            if (image->IsDefined(sample)) {
              VectorType value = image->Evaluate(sample);
              features[i] = value[0];
              //std::cout << i << ": " << samples[i] << " value: " << features[i] << std::endl;

            }
            else {
              // fail fast
              return false;
            }
            ++i;
          }



          output = features;
          return true;
        }
    };

    template<typename TPointSet, typename TImage>
    class ASMNormalDirectionImageValueExtractorFactory : public statismo::ASMFeatureExtractorFactory<TPointSet, TImage> {

        typedef ASMNormalDirectionPointSampler<TPointSet, TImage> SamplerType;
        typedef typename SamplerType::Pointer SamplerPointerType;
        typedef ASMNormalDirectionImageValueExtractor<TPointSet, TImage> InstanceType;

    public:

      static const ASMNormalDirectionImageValueExtractorFactory *GetInstance() {
        static ASMNormalDirectionImageValueExtractorFactory *instance = new ASMNormalDirectionImageValueExtractorFactory();
            return instance;
        }

        virtual std::string GetDescriptor() const {
            return "builtin::NormalDirection";
        }

        virtual const statismo::ASMFeatureExtractor<TPointSet, TImage> *Instantiate(
                const H5::Group &h5Group) const {

            unsigned int numberOfPoints = statismo::HDF5Utils::readInt(h5Group, "numberOfPoints");
            float spacing = statismo::HDF5Utils::readFloat(h5Group, "spacing");
            SamplerPointerType sampler = SamplerType::New();
            sampler->SetNumberOfPoints(numberOfPoints);
            sampler->SetPointSpacing(spacing);

            InstanceType* instance = new InstanceType(sampler);
            return instance;
        }

        virtual const statismo::ASMFeatureExtractor<TPointSet, TImage> *Instantiate( unsigned int numberOfPoints,
          float spacing) const {
          SamplerPointerType sampler = SamplerType::New();
          sampler->SetNumberOfPoints(numberOfPoints);
          sampler->SetPointSpacing(spacing);

          InstanceType* instance = new InstanceType(sampler);
          return instance;
        }

    };

}
#endif //STATISMO_itkASMNormalDirectionImageValueExtractor_H
