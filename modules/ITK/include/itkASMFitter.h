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

#include "ASMFitter.h"

namespace itk {

    class ASMFitterConfiguration : public Object, public statismo::ASMFitterConfiguration {
    public:

        typedef ASMFitterConfiguration Self;
        typedef Object Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(ASMFitterConfiguration, Object);
    };

    template<typename ASM>
    class ASMFitterResult : public Object, public statismo::ASMFitterResult<ASM> {
    public:

        typedef ASMFitterResult Self;
        typedef Object Superclass;
        typedef SmartPointer<Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro(Self);
        itkTypeMacro(ASMFitterResult, Object);
    };

    template<typename ASM>
    class ASMFitter : public Object, public statismo::ASMFitter<ASM> {
    public:

        typedef ASMFitter Self;
        typedef Object Superclass;
        typedef SmartPointer <Self> Pointer;
        typedef SmartPointer<const Self> ConstPointer;

        itkNewMacro( Self );
        itkTypeMacro( ASMFitter, Object);

        virtual typename ASM::FitterPointerType SetMesh(typename ASM::MeshPointerType mesh) {
            ASMFitter::Pointer copy = ASMFitter::New();
            copy->Init(this->m_configuration, this->m_model, mesh, this->m_targetImage, this->m_sampler->SetMesh(mesh), this->m_featureExtractor->SetMesh(mesh));
            return copy;
        }

    protected:
        virtual typename ASM::FitterResultPointerType InstantiateResult() {
            return ASM::FitterResultType::New();
        }
    };
}


#endif //STATISMO_ITKASMFITTER_H
