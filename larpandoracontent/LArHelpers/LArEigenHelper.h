/**
 *  @file   larpandoracontent/LArHelpers/LArEigenHelper.h
 *
 *  @brief  Header file for the principal curve analysis helper class.
 *
 *  $Log: $
 */
#ifndef LAR_EIGEN_HELPER_H
#define LAR_EIGEN_HELPER_H 1

#include "Objects/CartesianVector.h"

#include <Eigen/Dense>

#include <vector>

namespace lar_content
{

/**
 *  @brief  LArEigenHelper class
 */
class LArEigenHelper
{
public:
    /**
     *  @brief  Convert a container of calo hits into an Eigen matrix.
     *
     *  @param  caloHitContainer the calo hit list containing the hits from which to construct a maxtrix
     *  @param  hitMatrix the output Eigen matrix
     */
    template <class T>
    static void Vectorize(const T &caloHitContainer, Eigen::MatrixXf &hitMatrix);

    /**
     *  @brief  Convert a container of calo hits into a collection of Eigen matrices representing the centre, low and high coordinates of hits.
     *
     *  @param  caloHitContainer the calo hit list containing the hits from which to construct matrices
     *  @param  centre the output Eigen matrix of hit centres
     *  @param  low the output Eigen matrix of hit low edges (centre - width / 2)
     *  @param  high the output Eigen matrix of hit high edges (centre + width / 2)
     */
    template <class T>
    static void Vectorize(const T &caloHitContainer, Eigen::MatrixXf &centre, Eigen::MatrixXf &low, Eigen::MatrixXf &high);

    /**
     *  @brief  Convert a container of 3D calo hits into an Eigen matrix.
     *
     *  @param  caloHitContainer the calo hit list containing the hits from which to construct a maxtrix
     *  @param  hitMatrix the output Eigen matrix
     */
    template <class T>
    static void Vectorize3D(const T &caloHitContainer, Eigen::MatrixXf &hitMatrix);

    /**
     *  @brief  Retrieve the angle, coutner-clockwise relative to the x axis, between all hits in a matrix and a specified origin.
     *
     *  @param  hitMatrix the input collection of hits
     *  @param  origin the origin from which to measure hits
     *  @param  phis the output vector of angles in the range [0, 2pi]
     */
    static void GetAngles(const Eigen::MatrixXf &hitMatrix, const Eigen::RowVectorXf &origin, Eigen::RowVectorXf &phis);
};

} // namespace lar_content

#endif // #ifndef LAR_EIGEN_HELPER_H
