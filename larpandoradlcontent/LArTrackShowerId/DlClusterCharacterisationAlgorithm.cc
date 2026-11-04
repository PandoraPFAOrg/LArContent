/**
 *  @file   larpandoracontent/LArTrackShowerId/DlClusterCharacterisationAlgorithm.cc
 *
 *  @brief  Implementation of the deep learning based cluster characterisation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArTrackShowerId/DlClusterCharacterisationAlgorithm.h"

#include "larpandoracontent/LArObjects/LArCaloHit.h"

#include <numeric>

using namespace pandora;
using namespace lar_content;

namespace lar_dl_content
{

DlClusterCharacterisationAlgorithm::DlClusterCharacterisationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

DlClusterCharacterisationAlgorithm::~DlClusterCharacterisationAlgorithm()
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

bool DlClusterCharacterisationAlgorithm::IsClearTrack(const Cluster *const pCluster) const
{
    const OrderedCaloHitList &orderedCaloHitList{pCluster->GetOrderedCaloHitList()};
    CaloHitList caloHits;
    orderedCaloHitList.FillCaloHitList(caloHits);
    const CaloHitList &isolatedHits{pCluster->GetIsolatedCaloHitList()};
    caloHits.insert(caloHits.end(), isolatedHits.begin(), isolatedHits.end());
    
    FloatVector trackLikelihoods;
    try
    {
        for (const CaloHit *pCaloHit : caloHits)
        {
            const LArCaloHit *pLArCaloHit{dynamic_cast<const LArCaloHit *>(pCaloHit)};
            if (pLArCaloHit)
            {
                const float pTrack{pLArCaloHit->GetTrackScore()};
                const float pShower{pLArCaloHit->GetShowerScore()};
                
                if ((pTrack + pShower) > std::numeric_limits<float>::epsilon())
                    trackLikelihoods.emplace_back(pTrack / (pTrack + pShower));
            }
        }

        const unsigned long N{trackLikelihoods.size()};
        if (N > 0)
        {
            float mean{std::accumulate(std::begin(trackLikelihoods), std::end(trackLikelihoods), 0.f) / N};
            if (mean >= 0.5f)
                return true;
            else
                return false;
        }
    }
    catch (const StatusCodeException &)
    {
    }

    return true;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode DlClusterCharacterisationAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    return ClusterCharacterisationBaseAlgorithm::ReadSettings(xmlHandle);
}

} // namespace lar_dl_content
