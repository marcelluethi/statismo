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

#ifndef STATISMO_ITKASMFITTER_H
#define STATISMO_ITKASMFITTER_H

#include <itkObject.h>
#include <itkMacro.h>
#include "ASMFitter.h"
#include "itkActiveShapeModel.h"
#include "itkASMPointSampler.h"
namespace itk {
    template<typename TPointSet, typename TImage>
    class ASMFitterStep : public Object {
    public:
        typedef ASMFitterStep Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;
        itkNewMacro( Self );
        itkTypeMacro( Self, Object);

        typedef typename ActiveShapeModel<TPointSet, TImage>::Pointer ModelPointerType;
        typedef typename TPointSet::Pointer PointSetPointerType;
        typedef typename TImage::Pointer ImagePointerType;
        typedef statismo::ASMFitterConfiguration ConfigurationType;
        typedef typename ASMPointSampler<TPointSet>::Pointer SamplerPointerType;
        typedef statismo::ASMFitterResult ResultType;
        typedef statismo::ASMFitterStep<TPointSet, TImage> ImplType;

        ASMFitterStep(): m_model(0), m_target(0), m_configuration(0,0,0) {}

        void SetModel(ModelPointerType model) {
            m_model = model;
        }

        void SetSource(statismo::VectorType source) {
            m_source = source;
        }

        void SetTarget(ImagePointerType target) {
            m_target = target;
        }

        void SetSampler(SamplerPointerType sampler) {
            m_sampler = sampler;
        }

        void SetConfiguration(const ConfigurationType& configuration) {
            m_configuration = configuration;
        }

        void Update() {
            ImplType* impl = ImplType::Create(m_configuration, m_model->GetstatismoImplObj(), m_source, m_target, m_sampler);
            m_result = impl->Perform();
            delete impl;
        }

        ResultType GetOutput() {
            return m_result;
        }

    private:
        ModelPointerType m_model;
        statismo::VectorType m_source;
        ImagePointerType m_target;
        SamplerPointerType m_sampler;
        ConfigurationType m_configuration;
        ResultType m_result;
        //statismo::VectorType m_coefficients; // FIXME: WTF?
    };
}
//
//
#endif //STATISMO_ITKASMFITTER_H