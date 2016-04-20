/*
 * This file is part of the statismo library.
 *
 * Author: Christoph Langguth (christoph.langguth@unibas.ch)
 *
 * Copyright (c) 2011-2015 University of Basel
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 * Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * Redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution.
 *
 * Neither the name of the project's author nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
 * TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 */

#ifndef STATISMO_itkAsmNormalDirectionGradientGaussianFeatureExtractor_H
#define STATISMO_itkAsmNormalDirectionGradientGaussianFeatureExtractor_H

#include "itkASMFeatureExtractor.h"
#include "HDF5Utils.h"
#include "itkBSplineInterpolateImageFunction.h"
#include "itkASMNormalDirectionPointSampler.h"
#include "itkDiscreteGaussianImageFilter.h"

namespace itk {
    template <typename ASM>
    class ASMNormalDirectionGradientGaussianFeatureExtractor : public itk::ASMFeatureExtractor<ASM> {

    private:
        typedef BSplineInterpolateImageFunction<typename ASM::ImageType, typename ASM::ImagePixelType, typename ASM::ImagePixelType> InterpolatedImageType;
        typedef typename InterpolatedImageType::CovariantVectorType InterpolatedGradientType;
        typedef DiscreteGaussianImageFilter<typename ASM::ImageType, typename ASM::ImageType> GaussianFilterType;
        float m_sigma;
        typename ASM::ImagePointerType m_image;
        typename InterpolatedImageType::Pointer m_interpolatedImage;
        typename ASMNormalDirectionPointSampler<ASM>::Pointer m_sampler;

        void Init(typename ASMNormalDirectionPointSampler<ASM>::Pointer sampler, typename ASM::ImagePointerType image, typename InterpolatedImageType::Pointer interpolatedImage, float sigma) {
            m_sampler = sampler;
            m_image = image;
            m_interpolatedImage = interpolatedImage;
            m_sigma = sigma;
        }

        typename InterpolatedImageType::Pointer Interpolate(typename ASM::ImagePointerType image) {
            if (m_sigma != 0) {
                typename GaussianFilterType::Pointer smooth = GaussianFilterType::New();
                smooth->SetVariance(m_sigma * m_sigma);
                // FIXME: smooth->SetMaximumKernelWidth ???
                smooth->SetInput(image);
                smooth->Update();
                image = smooth->GetOutput();
            }

            typename InterpolatedImageType::Pointer inter = InterpolatedImageType::New();
            inter->SetSplineOrder(1);
            inter->SetInputImage(image);
            return inter;
        }

    public:
        typedef ASMNormalDirectionGradientGaussianFeatureExtractor Self;
        typedef itk::ASMFeatureExtractor<ASM> Superclass;
        typedef itk::SmartPointer<Self>                Pointer;
        typedef itk::SmartPointer<const Self>          ConstPointer;
        itkSimpleNewMacro( Self );
        itkTypeMacro( ASMNormalDirectionGradientGaussianFeatureExtractor, itk::ASMFeatureExtractor<ASM> );

        ASMNormalDirectionGradientGaussianFeatureExtractor() {
            //std::cout << "constructor: itk::ASMNormalDirectionGradientGaussianFeatureExtractor " << this << std::endl;
        }

        virtual ~ASMNormalDirectionGradientGaussianFeatureExtractor() {
            //std::cout << "destructor: itk::ASMNormalDirectionGradientGaussianFeatureExtractor " << this << std::endl;
        }


        void Load(const H5::Group &h5Group) {
            unsigned numberOfPoints = statismo::HDF5Utils::readInt(h5Group, "numberOfPoints");
            float spacing = statismo::HDF5Utils::readFloat(h5Group, "spacing");
            m_sampler = ASMNormalDirectionPointSampler<ASM>::New();
            m_sampler->Init(numberOfPoints, spacing);
            m_sigma = statismo::HDF5Utils::readFloat(h5Group, "sigma");;
        }

        virtual typename ASM::FeatureExtractorPointerType SetImage(typename ASM::ImagePointerType image) {
            if (m_image == image) {
                return typename ASM::FeatureExtractorPointerType(this);
            }
            typename ASMNormalDirectionGradientGaussianFeatureExtractor::Pointer copy = New();
            copy->Init(m_sampler, image, Interpolate(image), m_sigma);
            return typename ASM::FeatureExtractorPointerType(copy.GetPointer());
        }

        virtual typename ASM::FeatureExtractorPointerType SetMesh(typename ASM::MeshPointerType mesh) {
            if (m_sampler->GetMesh() == mesh) {
                return typename ASM::FeatureExtractorPointerType(this);
            }

            typename ASMNormalDirectionGradientGaussianFeatureExtractor::Pointer copy = New();
            typename ASM::PointSamplerPointerType samplerCopy = m_sampler->SetMesh(mesh);
            // make sure that garbage collection does not kick in while we're casting the naked pointer back to its actual class
            samplerCopy->Register();
            typename ASMNormalDirectionPointSampler<ASM>::Pointer typedSamplerCopy = typename ASMNormalDirectionPointSampler<ASM>::Pointer(static_cast<ASMNormalDirectionPointSampler<ASM>*>(samplerCopy.GetPointer()));
            copy->Init(typedSamplerCopy, m_image, m_interpolatedImage, m_sigma);
            samplerCopy->UnRegister();
            return typename ASM::FeatureExtractorPointerType(copy.GetPointer());
        }

        virtual statismo::VectorType ExtractFeatures(const typename ASM::PointType& point) const {
            statismo::VectorType features (m_sampler->GetNumberOfPoints());

            int i = 0;
            float sum = 0;
            InterpolatedGradientType normal;
            normal.CastFrom(m_sampler->GetNormalForPoint(point));
            std::vector<typename ASM::PointType> samples = m_sampler->SampleAtPoint(point);

            for(typename std::vector<typename ASM::PointType>::iterator it = samples.begin(); it != samples.end(); ++it) {
                typename ASM::PointType sample = *it;

                typename ASM::ImageType::IndexType index;
                if (m_image->TransformPhysicalPointToIndex(sample, index)) {
                    
                    features[i] = m_interpolatedImage->EvaluateDerivative(sample) * normal;
                } else {
                    features[i] = 999;
                }
                // keep track of sum of absolute values
                if (features[i] < 0) {
                    sum -= features[i];
                } else {
                    sum += features[i];
                }
                ++i;
            }
            //normalize
            if (sum > 0) {
                for (i = 0; i < features.size(); ++i) {
                    features[i] /= sum;
                }
            }
            return features;
        }
    };

    template <typename ASM>
    class ASMNormalDirectionGradientGaussianFeatureExtractorFactory : public itk::ASMFeatureExtractorFactory<ASM> {
    public:
        typedef ASMNormalDirectionGradientGaussianFeatureExtractorFactory Self;
        typedef itk::ASMFeatureExtractorFactory<ASM> Superclass;
        typedef itk::SmartPointer<Self>                Pointer;
        typedef itk::SmartPointer<const Self>          ConstPointer;
        itkSimpleNewMacro( Self );
        itkTypeMacro( ASMNormalDirectionGradientGaussianFeatureExtractorFactory, itk::ASMFeatureExtractorFactory<ASM> );

        static typename ASM::FeatureExtractorFactoryPointerType GetInstance() {
            static typename ASM::FeatureExtractorFactoryPointerType instance = typename ASM::FeatureExtractorFactoryPointerType(ASMNormalDirectionGradientGaussianFeatureExtractorFactory<ASM>::New().GetPointer());
            return instance;
        }

        static const std::string Descriptor() {
            static const std::string s("builtin::NormalDirectionGradientGaussian");
            return s;
        }

        virtual ~ASMNormalDirectionGradientGaussianFeatureExtractorFactory() {
        }

        virtual std::string GetDescriptor() const {
            return Descriptor();
        }

        typename ASM::FeatureExtractorPointerType Instantiate(const H5::Group &h5Group) const {
            typename ASMNormalDirectionGradientGaussianFeatureExtractor<ASM>::Pointer instance = ASMNormalDirectionGradientGaussianFeatureExtractor<ASM>::New();
            instance->Load(h5Group);
            return typename ASM::FeatureExtractorPointerType(instance.GetPointer());
        }

    protected:
        ASMNormalDirectionGradientGaussianFeatureExtractorFactory() {
        }

    private:
        ASMNormalDirectionGradientGaussianFeatureExtractorFactory(const ASMNormalDirectionGradientGaussianFeatureExtractorFactory &o) { }

        ASMNormalDirectionGradientGaussianFeatureExtractorFactory &operator=(const ASMNormalDirectionGradientGaussianFeatureExtractorFactory &o) { }
    };
}


#endif //STATISMO_itkAsmNormalDirectionGradientGaussianFeatureExtractor_H
