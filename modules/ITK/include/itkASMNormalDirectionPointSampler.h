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

#include <itkObject.h>
#include <itkMacro.h>
#include <itkPointsLocator.h>
#include "itkASMPointSampler.h"
#include "itkTriangleMeshAdapter.h"
#include "itkCovariantVector.h"

#include <itkVTKPolyDataWriter.h>

namespace itk {
    template<typename TPointSet, typename TImage>
    class ASMNormalDirectionPointSampler : public ASMPointSampler<TPointSet, TImage> {
    private:
        typedef TPointSet MeshType;
        typedef typename MeshType::Pointer MeshPointerType;
        typedef PointsLocator<typename MeshType::PointsContainer> PointsLocatorType;
        typedef typename PointsLocatorType::Pointer PointsLocatorPointerType;
        typedef typename PointsLocatorType::PointsContainer PointsContainerType;
        typedef Vector<typename MeshType::PixelType> VectorType;
        typedef TriangleMeshAdapter<typename MeshType::PixelType> MeshAdapterType;
        typedef typename MeshAdapterType::Pointer MeshAdapterPointerType;
        typedef typename MeshAdapterType::PointNormalType PointNormalType;
        typedef typename MeshAdapterType::PointNormalsContainerPointer PointNormalsContainerPointer;

        unsigned int m_numberOfPoints;
        float m_pointSpacing;
        const MeshPointerType m_mesh;
        const PointsLocatorPointerType m_locator;
        const PointNormalsContainerPointer m_normals;

    public:

        /* standard ITK typedefs and macros */
        typedef ASMNormalDirectionPointSampler Self;
        typedef ASMPointSampler<TPointSet, TImage> Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);

        itkTypeMacro(ASMNormalDirectionPointSampler, Object);


        typedef statismo::ActiveShapeModel<TPointSet, TImage> ActiveShapeModelType;
        typedef typename ActiveShapeModelType::RepresenterType::RigidTransformPointerType RigidTransformPointerType;
        typedef typename statismo::Representer<TPointSet>::PointType PointType;
        typedef ASMNormalDirectionPointSampler<TPointSet, TImage> SelfType;

        ASMNormalDirectionPointSampler():
                m_mesh(0),
                m_locator(0),
                m_normals(0),
                m_numberOfPoints(0),
                m_pointSpacing(0)
        {}

        ~ASMNormalDirectionPointSampler() {
        }

        // FIXME: This should be either protected or private, and requires *all* fields
        ASMNormalDirectionPointSampler(MeshPointerType mesh, PointsLocatorPointerType locator, PointNormalsContainerPointer normals, unsigned int numberOfPoints, float pointSpacing):
                m_mesh (mesh),
                m_locator(locator),
                m_normals(normals),
                m_numberOfPoints(numberOfPoints),
                m_pointSpacing(pointSpacing) {

        }

        unsigned int GetNumberOfPoints() const {
            return m_numberOfPoints;
        }

        void SetNumberOfPoints(unsigned int numberOfPoints) {
            m_numberOfPoints = numberOfPoints;
        }

        float GetPointSpacing() const {
            return m_pointSpacing;
        }

        void SetPointSpacing(float pointSpacing) {
            m_pointSpacing = pointSpacing;
        }

        virtual SelfType *CloneForTarget(const ActiveShapeModelType *const model,
                                         const statismo::VectorType &coefficients, const RigidTransformPointerType transform) const {

            MeshPointerType mesh = model->GetStatisticalModel()->DrawSample(coefficients, false);
            if (transform) {
                mesh = model->GetRepresenter()->TransformMesh(mesh, transform);
            }

            PointsLocatorPointerType locator = PointsLocatorType::New();
            locator->SetPoints(const_cast<PointsContainerType *>(mesh->GetPoints()));
            locator->Initialize();

            MeshAdapterPointerType adapter = MeshAdapterType::New();
            adapter->SetMesh(mesh);
            PointNormalsContainerPointer normals = adapter->GetPointNormals();

            Pointer copy = new ASMNormalDirectionPointSampler(mesh, locator, normals, m_numberOfPoints, m_pointSpacing);
            //Pointer copy = ASMNormalDirectionPointSampler::Clone(); // this does NOT work: at runtime, it's resolved to the superclass....
            return copy.GetPointer();
        }

        virtual std::vector<PointType> Sample(const PointType &targetPoint) const {
            //FIXME: this only works for TriangleMeshes for now.

            //std::cout << "SAMPLING AT " << samplePoint << std::endl;
            std::vector<PointType> samples;
            samples.reserve(m_numberOfPoints);

            unsigned int normalPointId = m_locator->FindClosestPoint(targetPoint);

            // Convert to an itk::(ContraVariant)Vector, because no operators are overloaded for adding the
            // covariant vector normal to a point.
            VectorType normalVector(GetNormalForPointId(normalPointId).GetDataPointer());

            int startInclusive = -(m_numberOfPoints / 2);
            int endExclusive = (m_numberOfPoints + 1) / 2;

            for (int i = startInclusive; i < endExclusive; ++i) {
                PointType sample = targetPoint + normalVector * i * m_pointSpacing;
                samples.push_back(sample);
            }

            return samples;
        }


        PointNormalType GetNormalForPointId(const unsigned &targetPointId) const {
            PointNormalType normal = m_normals->GetElement(targetPointId);
            return normal * (1.0 / normal.GetNorm());
        }

        PointNormalType GetNormalForPoint(const PointType &targetPoint) const {
            return GetNormalForPointId(m_locator->FindClosestPoint(targetPoint));
        }


        void WriteMesh(const char* fileName)
        {
          

          auto polyDataWriter = VTKPolyDataWriter<MeshType>::New();
          polyDataWriter->SetFileName(fileName);
          polyDataWriter->SetInput(m_mesh);
          polyDataWriter->Update();


        }

    private:

    };

}
#endif //STATISMO_ASMNORMALDIRECTIONPOINTSAMPLER_H
