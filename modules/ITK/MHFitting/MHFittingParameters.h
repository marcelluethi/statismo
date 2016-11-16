//
// Created by luetma00 on 15.11.16.
//

#ifndef STATISMO_MHFITTINGPARAMETERS_H
#define STATISMO_MHFITTINGPARAMETERS_H


#include "CommonTypes.h"
#include "Types.h"

namespace mhfitting {


    class MHFittingParameters {


    public:
        MHFittingParameters() {
        }

        MHFittingParameters(statismo::VectorType coefficients, statismo::VectorType rigidParameters, statismo::VectorType rotationCenter) :
                m_coefficients(coefficients),
                m_rigidParameters(rigidParameters),
                m_rotationCenter(rotationCenter) {}

        statismo::VectorType GetCoefficients() const {
            return m_coefficients;
        }

        statismo::VectorType GetRigidTransformParameters() const {
            return m_rigidParameters;
        }

        statismo::VectorType GetRotationCenter() const {
            return m_rotationCenter;
        }


        int size() {
            // TODO: if this is really used we should change to a different vectorized representation of all parameters
            return m_coefficients.size() + m_rigidParameters.size();
        }

        float &operator[](int i) {
            // TODO: if this is really used we should change to a different vectorized representation of all parameters
            throw new std::runtime_error("operator[] not implemented for MHFittingResult");
            return m_coefficients(0, 0);
        }

    private:
        statismo::VectorType m_coefficients;
        statismo::VectorType m_rigidParameters;
        statismo::VectorType m_rotationCenter;
    };

}

#endif //STATISMO_MHFITTINGPARAMETERS_H
