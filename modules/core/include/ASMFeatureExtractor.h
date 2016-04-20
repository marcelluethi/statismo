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

#ifndef STATISMO_ASMFEATUREEXTRACTOR_H
#define STATISMO_ASMFEATUREEXTRACTOR_H

#include "CommonTypes.h"

namespace H5 {
    class Group;
}

namespace statismo {
    template <typename ASM>
    class ASMFeatureExtractor {
    public:
        virtual typename ASM::FeatureExtractorPointerType SetImage(typename ASM::ImagePointerType) = 0;
        virtual typename ASM::FeatureExtractorPointerType SetMesh(typename ASM::MeshPointerType) = 0;
        virtual statismo::VectorType ExtractFeatures(const typename ASM::PointType& point) const = 0;
    };

    template<typename ASM>
    class ASMFeatureExtractorFactory {
    private:
        static std::vector<typename ASM::FeatureExtractorFactoryPointerType> *implementations() {
            static std::vector<typename ASM::FeatureExtractorFactoryPointerType> impls;
            return &impls;
        }
        ASMFeatureExtractorFactory(const ASMFeatureExtractorFactory& o) { }
        ASMFeatureExtractorFactory& operator=(const ASMFeatureExtractorFactory& o) {}

    protected:
        ASMFeatureExtractorFactory() {}

    public:
        virtual std::string GetDescriptor() const = 0;
        virtual typename ASM::FeatureExtractorPointerType Instantiate(const H5::Group& h5Group) const = 0;

        static std::vector<typename ASM::FeatureExtractorFactoryPointerType> GetImplementations() {
            return *implementations();
        }

        static void RegisterImplementation(typename ASM::FeatureExtractorFactoryPointerType impl) {
            if (GetImplementation(impl->GetDescriptor())) {
                //ignoring already-registered implementation
                return;
            }
            implementations()->push_back(impl);
        }

        static typename ASM::FeatureExtractorFactoryPointerType GetImplementation(std::string descriptor) {
            std::vector<typename ASM::FeatureExtractorFactoryPointerType> impls = GetImplementations();
            for (typename std::vector<typename ASM::FeatureExtractorFactoryPointerType>::iterator impl = impls.begin(); impl != impls.end(); ++impl) {
                if ((*impl)->GetDescriptor() == descriptor) {
                    return *impl;
                }
            }
            return 0;
        }
    };
}

#endif //STATISMO_ASMFEATUREEXTRACTOR_H
