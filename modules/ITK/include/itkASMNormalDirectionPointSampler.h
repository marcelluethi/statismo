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

#ifndef STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H
#define STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H

#include "itkPointsLocator.h"
#include "itkASMPointSampler.h"

namespace itk {
    template <typename ASM>
    class ASMNormalDirectionPointSampler : public ASMPointSampler<ASM> {
    public:

        typedef ASMNormalDirectionPointSampler Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(ASMNormalDirectionPointSampler, Object);

        virtual ~ASMNormalDirectionPointSampler() {
                //std::cout << "destructor: itk::ASMNormalDirectionPointSampler" << std::endl;
        };

        void Init(unsigned numberOfPoints, float spacing, typename ASM::MeshPointerType mesh = 0) {
            m_numberOfPoints = numberOfPoints;
            m_spacing = spacing;
            SetMeshInternal(mesh);
        }

        virtual typename ASM::PointSamplerPointerType SetMesh(typename ASM::MeshPointerType mesh) {
            if (mesh == m_mesh) {
                //std::cout << "PointSampler::SetMesh() returning self" << std::endl;
                return typename ASM::PointSamplerPointerType(this);
            }
            //std::cout << "PointSampler::SetMesh() returning new" << std::endl;
            typename ASMNormalDirectionPointSampler<ASM>::Pointer clone = New();
            clone->Init(m_numberOfPoints, m_spacing, mesh);
            return typename ASM::PointSamplerPointerType(clone.GetPointer());
        }

        unsigned GetNumberOfPoints() const {
            return m_numberOfPoints;
        }

        float GetSpacing() const {
            return m_spacing;
        }

        typename ASM::MeshPointerType GetMesh() const {
            return m_mesh;
        }

        virtual std::vector<typename ASM::PointType> SampleAtPoint(const typename ASM::PointType &targetPoint) {
            return Sample(targetPoint, FindClosestPointId(targetPoint));
        }

        virtual std::vector<typename ASM::PointType> SampleAtPointId(const unsigned &targetPointId) {
            typename ASM::PointType point = m_mesh->GetPoint(targetPointId);
            return Sample(point, targetPointId);
        }

        typename ASM::MeshAdapterType::PointNormalType GetNormalForPointId(const unsigned &targetPointId) const {
            typename ASM::MeshAdapterType::PointNormalType normal = m_normals->GetElement(targetPointId);
            return normal * (1.0 / normal.GetNorm());
        }

        typename ASM::MeshAdapterType::PointNormalType GetNormalForPoint(const typename ASM::PointType &targetPoint) {
            return GetNormalForPointId(FindClosestPointId(targetPoint));
        }

    private:
        typedef itk::PointsLocator<typename ASM::MeshType::PointsContainer> PointsLocatorType;
        unsigned m_numberOfPoints;
        float m_spacing;
        typename ASM::MeshPointerType m_mesh;
        typename ASM::MeshAdapterType::PointNormalsContainer::Pointer m_normals;
        typename PointsLocatorType::Pointer m_locator;
        bool m_locatorSet;

        std::vector<typename ASM::PointType> Sample(const typename ASM::PointType &targetPoint, const unsigned normalPointId) const {

            std::vector<typename ASM::PointType> samples;
            samples.reserve(m_numberOfPoints);

            int startInclusive = -(m_numberOfPoints / 2);
            int endExclusive = (m_numberOfPoints + 1) / 2;

            typename ASM::MeshAdapterType::PointNormalType normal = GetNormalForPointId(normalPointId);
            // Convert to an itk::(ContraVariant)Vector, because no operators are overloaded for adding the
            // covariant vector normal to a point.
            typename ASM::VectorType normalVector(normal.GetDataPointer());

            //std::cout << "sampler indexes: startInclusive="<<startInclusive<<" endExlusive=" << endExclusive << std::endl;
            for(int i = startInclusive; i < endExclusive; ++i) {
                typename ASM::PointType sample = targetPoint + normalVector * i * m_spacing;
                //std::cout << "sample(" << i <<") = " << sample << std::endl;
                samples.push_back(sample);
            }

            return samples;
        }

        void EnsureLocatorIsInitialized() {
            if (!m_locatorSet) {
                m_locator = PointsLocatorType::New();
                m_locator->SetPoints(m_mesh->GetPoints());
                m_locator->Initialize();
                m_locatorSet = true;
            }
        }

        const unsigned FindClosestPointId(const typename ASM::PointType &targetPoint) {
            EnsureLocatorIsInitialized();
            return m_locator->FindClosestPoint(targetPoint);
        }

        void SetMeshInternal(typename ASM::MeshPointerType mesh) {
            m_mesh = mesh;
            m_locatorSet = false;

            if (mesh) {
                typename ASM::MeshAdapterType::Pointer adapter = ASM::MeshAdapterType::New();
                adapter->SetMesh(mesh);
                m_normals = adapter->GetPointNormals();

            }
        }
    };


}
#endif //STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H
