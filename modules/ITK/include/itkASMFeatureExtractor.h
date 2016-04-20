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

#ifndef STATISMO_ITKASMFEATUREEXTRACTOR_H
#define STATISMO_ITKASMFEATUREEXTRACTOR_H

#include <itkObjectFactory.h>
#include "ASMFeatureExtractor.h"

namespace itk {
    template <typename ASM>
    class ASMFeatureExtractor: public Object, public statismo::ASMFeatureExtractor<ASM> {
    public:
        typedef ASMFeatureExtractor Self;
        typedef Object Superclass;
        typedef itk::SmartPointer<Self>                Pointer;
        typedef itk::SmartPointer<const Self>          ConstPointer;

        itkSimpleNewMacro( Self );
        itkTypeMacro( ASMFeatureExtractor, Object );

        virtual typename ASM::FeatureExtractorPointerType SetImage(typename ASM::ImagePointerType) = 0;
        virtual typename ASM::FeatureExtractorPointerType SetMesh(typename ASM::MeshPointerType) = 0;
        virtual statismo::VectorType ExtractFeatures(const typename ASM::PointType& point) const = 0;
    };

    template <typename ASM>
    class ASMFeatureExtractorFactory: public Object, public statismo::ASMFeatureExtractorFactory<ASM> {
    public:
        typedef ASMFeatureExtractorFactory Self;
        typedef Object Superclass;
        typedef itk::SmartPointer<Self>                Pointer;
        typedef itk::SmartPointer<const Self>          ConstPointer;

        itkSimpleNewMacro( Self );
        itkTypeMacro( ASMFeatureExtractorFactory, Object );
    };

}
#endif //STATISMO_ITKASMFEATUREEXTRACTOR_H
