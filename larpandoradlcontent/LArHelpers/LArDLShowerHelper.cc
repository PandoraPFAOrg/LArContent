/**
 *  @file   larpandoradlcontent/LArHelpers/LArDLShowerHelper.cc
 *
 *  @brief  Implementation of the lar deep learning helper helper class.
 *
 *  $Log: $
 */

#include "Objects/CaloHit.h"
#include "Objects/CartesianVector.h"

#include "larpandoradlcontent/LArHelpers/LArDLShowerHelper.h"

#include <set>

namespace lar_dl_content
{

using namespace pandora;

//-----------------------------------------------------------------------------------------------------------------------------------------

LArDLShowerHelper::HitFeatures::HitFeatures() :
    m_xRel{0.f},
    m_zRel{0.f},
    m_rRel{0.f},
    m_cosThetaRel{0.f},
    m_sinThetaRel{0.f},
    m_distToXGap{0.f},
    m_xWidth{0.f},
    m_energy{0.f}
{
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void LArDLShowerHelper::CalculateHitFeatures(const CaloHit *const pCaloHit, const std::set<float> &detXGaps, const CartesianVector &vtxPos, HitFeatures &hitFeatures)
{
    const double x(pCaloHit->GetPositionVector().GetX());
    const double xRel(x - vtxPos.GetX());
    const double z(pCaloHit->GetPositionVector().GetZ());
    const double zRel(z - vtxPos.GetZ());
    const double rRel(std::sqrt(std::pow(xRel, 2.) + std::pow(zRel, 2.)));
    const double cosThetaRel(rRel != 0. ? xRel / rRel : 0.);
    const double sinThetaRel(rRel != 0. ? zRel / rRel : 0.);

    double distToXGap(std::numeric_limits<double>::max());
    for (const double xGap : detXGaps)
    {
        const double dist(x - xGap);
        if (std::abs(dist) < std::abs(distToXGap))
            distToXGap = dist;
    }

    hitFeatures.m_xRel = static_cast<float>(xRel);
    hitFeatures.m_zRel = static_cast<float>(zRel);
    hitFeatures.m_rRel = static_cast<float>(rRel);
    hitFeatures.m_cosThetaRel = static_cast<float>(cosThetaRel);
    hitFeatures.m_sinThetaRel = static_cast<float>(sinThetaRel);
    hitFeatures.m_distToXGap = static_cast<float>(distToXGap);
    hitFeatures.m_xWidth = pCaloHit->GetCellSize1();
    hitFeatures.m_energy = pCaloHit->GetMipEquivalentEnergy();
}

} // namespace lar_dl_content
