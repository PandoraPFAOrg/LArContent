/**
 *  @file   larpandoracontent/LArHelpers/LArEigenHelper.cc
 *
 *  @brief  Implementation of the principal curve analysis helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArHelpers/LArEigenHelper.h"
#include "larpandoracontent/LArHelpers/LArClusterHelper.h"
#include "larpandoracontent/LArHelpers/LArObjectHelper.h"

#include <Eigen/Dense>

#define _USE_MATH_DEFINES
#include <cmath>

using namespace pandora;

namespace lar_content
{

//------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
void LArEigenHelper::Vectorize(const T &caloHitContainer, Eigen::MatrixXf &hitMatrix)
{
    int i{0};
    for (const CaloHit *const pCaloHit : caloHitContainer)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        hitMatrix(i, 0) = pos.GetX();
        hitMatrix(i, 1) = pos.GetZ();
        ++i;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
void LArEigenHelper::Vectorize(const T &caloHitContainer, Eigen::MatrixXf &centre, Eigen::MatrixXf &low, Eigen::MatrixXf &high)
{
    LArEigenHelper::Vectorize(caloHitContainer, centre);
    int i{0};
    for (const CaloHit *const pCaloHit : caloHitContainer)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        low(i, 0) = pos.GetX() - pCaloHit->GetCellSize1() * 0.5f;
        low(i, 1) = pos.GetZ();
        high(i, 0) = pos.GetX() + pCaloHit->GetCellSize1() * 0.5f;
        high(i, 1) = pos.GetZ();
        ++i;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

template <class T>
void LArEigenHelper::Vectorize3D(const T &caloHitContainer, Eigen::MatrixXf &hitMatrix)
{
    int i{0};
    for (const CaloHit *const pCaloHit : caloHitContainer)
    {
        const CartesianVector &pos{pCaloHit->GetPositionVector()};
        hitMatrix(i, 0) = pos.GetX();
        hitMatrix(i, 1) = pos.GetY();
        hitMatrix(i, 2) = pos.GetZ();
        ++i;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArEigenHelper::GetAngles(const Eigen::MatrixXf &hitMatrix, const Eigen::RowVectorXf &origin, Eigen::RowVectorXf &phis)
{
    const float pi{static_cast<float>(M_PI)};
    Eigen::RowVectorXf piVec{Eigen::RowVectorXf::Constant(hitMatrix.rows(), pi)};
    Eigen::RowVectorXf zeroVec{Eigen::RowVectorXf::Zero(hitMatrix.rows())};
    Eigen::MatrixXf deltas(hitMatrix.rowwise() - origin);
    for (int i = 0; i < deltas.rows(); ++i)
        phis(i) = std::atan2(deltas(i, 1), deltas(i, 0));
    // Move from [-pi, +pi] to [0, 2pi)
    phis = (phis.array() < 0).select(piVec * 2 + phis, phis);
    phis = (phis.array() < 2 * pi).select(phis, zeroVec);
}

//------------------------------------------------------------------------------------------------------------------------------------------

template void LArEigenHelper::Vectorize(const CaloHitList &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize(const CaloHitVector &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize(const CaloHitSet &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize3D(const CaloHitList &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize3D(const CaloHitVector &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize3D(const CaloHitSet &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize(const CaloHitList &, Eigen::MatrixXf &, Eigen::MatrixXf &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize(const CaloHitVector &, Eigen::MatrixXf &, Eigen::MatrixXf &, Eigen::MatrixXf &);
template void LArEigenHelper::Vectorize(const CaloHitSet &, Eigen::MatrixXf &, Eigen::MatrixXf &, Eigen::MatrixXf &);

} // namespace lar_content
