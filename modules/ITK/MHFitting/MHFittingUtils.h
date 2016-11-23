//
// Created by luetma00 on 17.11.16.
//

#ifndef STATISMO_UTILS_H
#define STATISMO_UTILS_H

#include "Types.h"
#include <Eigen/SVD>
#include "../../core/include/CommonTypes.h"
#include <vector>
#include <stdexcept>

namespace mhfitting {



    class MHFittingUtils {

    public:

        struct ProcrustesResult {
            statismo::MatrixTypeDoublePrecision rotationMatrix;
            statismo::VectorTypeDoublePrecision translationVector;
        };

    static ProcrustesResult computeProcrustes(const std::vector<PointType>& fixedPoints, const std::vector<PointType>& movingPoints, const PointType& center) {

        if (fixedPoints.size() != movingPoints.size()) {
            throw std::logic_error("we need the same number of fixed and moving points");
        }

        unsigned n = fixedPoints.size();

        if (n < 3) {
            throw std::logic_error("we need at least 3 poitns ot compute the solution");
        }

        statismo::MatrixTypeDoublePrecision X = statismo::MatrixTypeDoublePrecision::Zero(n, 3);
        statismo::MatrixTypeDoublePrecision Y = statismo::MatrixTypeDoublePrecision::Zero(n, 3);

        for (unsigned i = 0; i < n; ++i) {
            for (unsigned d = 0; d < 3; ++d) {
                X(i, d) = fixedPoints[i].GetElement(d) - center.GetElement(d);
                Y(i, d) = movingPoints[i].GetElement(d) - center.GetElement(d);
            }
        }


        statismo::VectorTypeDoublePrecision meanX = X.colwise().mean();


        assert(meanX.rows() == 3);

        statismo::VectorTypeDoublePrecision meanY = Y.colwise().mean();
        assert(meanY.rows() == 3);

        double sigma2_x = 0.0;
        for (unsigned i = 0; i < X.cols(); ++i) {
            sigma2_x += (X.row(i).transpose() - meanX).dot(X.row(i).transpose() - meanX);
        }
        sigma2_x /= n;

        statismo::MatrixTypeDoublePrecision sigma_xy = statismo::MatrixTypeDoublePrecision::Zero(3,3);
        for (unsigned i = 0; i < n; ++i) {
            sigma_xy += (Y.row(i).transpose() - meanY) * (X.row(i).transpose() - meanX).transpose();
        }
        sigma_xy /= n;


        Eigen::JacobiSVD<statismo::MatrixTypeDoublePrecision> svd(sigma_xy,  Eigen::ComputeFullU | Eigen::ComputeFullV);

        statismo::MatrixTypeDoublePrecision S = statismo::MatrixTypeDoublePrecision::Identity(3,3);
        if (sigma_xy.determinant() < 0) {
            S(2, 2) = -1;
        }
        statismo::MatrixTypeDoublePrecision R = svd.matrixU() * S * svd.matrixV().transpose();


        double trDS = S.diagonal().dot(svd.singularValues());

        double c = (1.0 / (n * sigma2_x)) * trDS;

        statismo::VectorTypeDoublePrecision t = meanY - R * meanX;
        ProcrustesResult procResult;
        procResult.rotationMatrix = R;
        procResult.translationVector = t;


        return procResult;


    }
};
}

#endif //STATISMO_UTILS_H
