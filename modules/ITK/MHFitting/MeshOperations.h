//
// Created by luetma00 on 15.11.16.
//

#ifndef STATISMO_MESHOPERATIONS_H
#define STATISMO_MESHOPERATIONS_H

#include "Types.h"
#include "MHFittingParameters.h"

namespace mhfitting {

    template <class MeshType, class PointType>
    class MeshOperations {
    public:
        virtual std::pair<PointType, long> findClosestPoint(MeshType mesh, PointType pt) const = 0;
        virtual statismo::VectorType normalAtPoint(MeshType mesh, PointType pt) const = 0;
        virtual short huAtPoint(MeshType mesh, long ptId) const = 0;
        virtual MeshType transformMesh(const MHFittingParameters& fittingParameters) const = 0;
        virtual MeshType transformCenterLine(const MHFittingParameters& fittingParameters) const = 0;
        virtual PointType getPointWithId(MeshType mesh, unsigned id) const = 0;
        virtual PointType transformToModelSpace(const mhfitting::MHFittingParameters & fittingParameters,  PointType pt) const = 0;
        virtual MHFittingParameters::RigidParameters alignModelInstanceToTargetPoints(const MHFittingParameters &fittingParameters,
                                                                     const std::vector<PointType> &targetPoints) const = 0;
        virtual MHFittingParameters::ModelParameters projectModelInstanceToTargetPoints(const MHFittingParameters &fittingParameters,
                                                                       const CorrespondencePoints& correpsondencePoints,
                                                                        const std::vector<PointType> &targetPoints) const = 0;

    };
}

#endif //STATISMO_MESHOPERATIONS_H