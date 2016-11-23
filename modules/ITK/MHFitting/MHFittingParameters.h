//
// Created by luetma00 on 15.11.16.
//

#ifndef STATISMO_MHFITTINGPARAMETERS_H
#define STATISMO_MHFITTINGPARAMETERS_H


#include "CommonTypes.h"
#include "Types.h"
#include <boost/serialization/strong_typedef.hpp>

namespace mhfitting {


    class MHFittingParameters {


    public:

        struct Parameters {
            Parameters(const statismo::VectorType& v) : m_v(v) {}
            const statismo::VectorType& GetParameters() const { return m_v; }
        private:
            statismo::VectorType m_v;
        };

        struct RotationParameters : public Parameters {
            RotationParameters(const statismo::VectorType& v) : Parameters(v) {}
        };

        struct TranslationParameters : public Parameters {
            TranslationParameters(const statismo::VectorType& v) : Parameters(v) {}
        };

        struct ModelParameters : public Parameters {
            ModelParameters(const statismo::VectorType& v) : Parameters(v) {}
        };

        struct RotationCenter : public Parameters {
            RotationCenter(const statismo::VectorType& v) : Parameters(v) {}
        };


        MHFittingParameters() :
                m_coefficients(ModelParameters(statismo::VectorType::Zero(200))),
                m_rotationParameters(RotationParameters(statismo::VectorType::Zero(3))),
                m_translationParameters(TranslationParameters(statismo::VectorType::Zero(3))),
                m_rotationCenter(statismo::VectorType::Zero(3))
        {}


        MHFittingParameters(const ModelParameters& coefficients, const RotationParameters& rotationParameters, const TranslationParameters& translationParameters, const RotationCenter& rotationCenter) :
                m_coefficients(coefficients),
                m_rotationParameters(rotationParameters),
                m_translationParameters(translationParameters),
                m_rotationCenter(rotationCenter) {}

        const ModelParameters& GetCoefficients() const {
            return m_coefficients;
        }

       const RotationParameters& GetRotationParameters() const {
            return m_rotationParameters;
        }

        const TranslationParameters& GetTranslationParameters() const {
            return m_translationParameters;
        }


        const RotationCenter& GetRotationCenter() const {
            return m_rotationCenter;
        }


        int size() const {
            // TODO: if this is really used we should change to a different vectorized representation of all parameters
            return m_coefficients.GetParameters().size() + m_rotationParameters.GetParameters().size() + m_translationParameters.GetParameters().size();
        }

        float &operator[](int i) {
            // TODO: if this is really used we should change to a different vectorized representation of all parameters
            throw new std::runtime_error("operator[] not implemented for MHFittingResult");
            statismo::VectorType p = m_coefficients.GetParameters();
            return p(0);
        }

    private:
        ModelParameters m_coefficients;
        RotationParameters m_rotationParameters;
        TranslationParameters m_translationParameters;
        RotationCenter m_rotationCenter;
    };

}

#endif //STATISMO_MHFITTINGPARAMETERS_H
