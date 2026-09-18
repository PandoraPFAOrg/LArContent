/**
 *  @file   larpandoracontent/LArObjects/LArPlaneContextObject.cc
 *
 *  @brief  Implementation of the cluster helper class.
 *
 *  $Log: $
 */

#include "larpandoracontent/LArObjects/LArPlaneContextObject.h"

using namespace pandora;

namespace lar_content
{

bool LArPlaneContextObject::AddHitTriplet(const CaloHit *const pHitU, const CaloHit *const pHitV, const CaloHit *const pHitW)
{
    const int nValidHits{(pHitU ? 1 : 0) + (pHitV ? 1 : 0) + (pHitW ? 1 : 0)};
    if (nValidHits < 2)
        return false;

    auto hitTriplet(std::make_unique<HitTriplet>(HitTriplet{pHitU, pHitV, pHitW}));
    HitTriplet *raw = hitTriplet.get();
    m_hitTriplets.emplace_back(std::move(hitTriplet));

    if (pHitU)
    {
        PANDORA_THROW_IF(STATUS_CODE_ALREADY_PRESENT, m_uIndex.find(pHitU) != m_uIndex.end());
        m_uIndex[pHitU] = raw;
    }
    if (pHitV)
    {
        PANDORA_THROW_IF(STATUS_CODE_ALREADY_PRESENT, m_vIndex.find(pHitV) != m_vIndex.end());
        m_vIndex[pHitV] = raw;
    }
    if (pHitW)
    {
        PANDORA_THROW_IF(STATUS_CODE_ALREADY_PRESENT, m_wIndex.find(pHitW) != m_wIndex.end());
        m_wIndex[pHitW] = raw;
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

const LArPlaneContextObject::HitTriplet *LArPlaneContextObject::GetHitTriplet(const CaloHit *const pCaloHit) const
{
    if (m_uIndex.count(pCaloHit))
        return m_uIndex.at(pCaloHit);
    if (m_vIndex.count(pCaloHit))
        return m_vIndex.at(pCaloHit);
    if (m_wIndex.count(pCaloHit))
        return m_wIndex.at(pCaloHit);

    return nullptr;
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPlaneContextObject::HitTripletIterator LArPlaneContextObject::begin() const
{
    return HitTripletIterator(m_hitTriplets.begin());
}

//------------------------------------------------------------------------------------------------------------------------------------------

LArPlaneContextObject::HitTripletIterator LArPlaneContextObject::end() const
{
    return HitTripletIterator(m_hitTriplets.end());
}

//------------------------------------------------------------------------------------------------------------------------------------------

size_t LArPlaneContextObject::Size() const
{
    return m_hitTriplets.size();
}

} // namespace lar_content
