/**
 *  @file   larpandoracontent/LArControlFlow/PreProcessingAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "Pandora/AlgorithmHeaders.h"

#include "larpandoracontent/LArControlFlow/PreProcessingAlgorithm.h"

#include "larpandoracontent/LArHelpers/LArClusterHelper.h"

#include "larpandoracontent/LArUtility/KDTreeLinkerAlgoT.h"

using namespace pandora;

namespace lar_content
{

PreProcessingAlgorithm::PreProcessingAlgorithm() :
    m_mipEquivalentCut(std::numeric_limits<float>::epsilon()),
    m_minCellLengthScale(std::numeric_limits<float>::epsilon()),
    m_maxCellLengthScale(3.f),
    m_searchRegion1D(0.1f),
    m_maxEventHits(std::numeric_limits<unsigned int>::max()),
    m_onlyAvailableCaloHits(true),
    m_inputCaloHitListName("Input")
{
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PreProcessingAlgorithm::Reset()
{
    m_processedHits.clear();
    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PreProcessingAlgorithm::Run()
{
    if (!this->GetPandora().GetSettings()->SingleHitTypeClusteringMode())
    {
        std::cout << "PreProcessingAlgorithm: expect Pandora to be configured in SingleHitTypeClusteringMode." << std::endl;
        return STATUS_CODE_FAILURE;
    }

    try
    {
        this->ProcessCaloHits();
    }
    catch (const StatusCodeException &statusCodeException)
    {
        if (STATUS_CODE_OUT_OF_RANGE == statusCodeException.GetStatusCode())
        {
            std::cout << "PreProcessingAlgorithm: Excessive number of hits in event, skipping the reconstruction" << std::endl;
            this->PopulateVoidCaloHitLists();
        }
    }

    if (!m_currentCaloHitListReplacement.empty())
    {
        if (STATUS_CODE_SUCCESS != PandoraContentApi::ReplaceCurrentList<CaloHit>(*this, m_currentCaloHitListReplacement))
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PreProcessingAlgorithm: could not replace current calo hit list with list named: " << m_currentCaloHitListReplacement
                          << std::endl;
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PreProcessingAlgorithm::ProcessCaloHits()
{
    const CaloHitList *pCaloHitList(nullptr);
    PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::GetList(*this, m_inputCaloHitListName, pCaloHitList));

    std::cout << "Input name '" << m_inputCaloHitListName << "' with " << pCaloHitList->size() << " hits" << std::endl;

    if (pCaloHitList->empty())
        return;

    if (pCaloHitList->size() > m_maxEventHits)
        throw StatusCodeException(STATUS_CODE_OUT_OF_RANGE);

    CaloHitList selectedCaloHitListU, selectedCaloHitListV, selectedCaloHitListW, selectedCaloHitListOp;
    for (const CaloHit *const pCaloHit : *pCaloHitList)
    {
        if (m_processedHits.count(pCaloHit))
            continue;
        m_processedHits.insert(pCaloHit);

        switch (pCaloHit->GetHitType())
        {
            case OPTICAL_TPC:
            case OPTICAL_TRAP:
            case OPTICAL_SIPM:
                this->ProcessOpticalHit(pCaloHit, selectedCaloHitListOp);
                break;
            case TPC_VIEW_U:
                this->ProcessChargeHit(pCaloHit, selectedCaloHitListU);
                break;
            case TPC_VIEW_V:
                this->ProcessChargeHit(pCaloHit, selectedCaloHitListV);
                break;
            case TPC_VIEW_W:
                this->ProcessChargeHit(pCaloHit, selectedCaloHitListW);
                break;
            default:
                break;
        }

    }

    CaloHitList filteredCaloHitListU, filteredCaloHitListV, filteredCaloHitListW, filteredCaloHitListOp;
    this->GetFilteredCaloHitList(selectedCaloHitListU, filteredCaloHitListU);
    this->GetFilteredCaloHitList(selectedCaloHitListV, filteredCaloHitListV);
    this->GetFilteredCaloHitList(selectedCaloHitListW, filteredCaloHitListW);
    filteredCaloHitListOp.insert(filteredCaloHitListOp.end(), selectedCaloHitListOp.begin(), selectedCaloHitListOp.end());

    CaloHitList filteredInputList;
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListU.begin(), filteredCaloHitListU.end());
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListV.begin(), filteredCaloHitListV.end());
    filteredInputList.insert(filteredInputList.end(), filteredCaloHitListW.begin(), filteredCaloHitListW.end());

    std::cout << "Saw " << selectedCaloHitListOp.size() << " optical hits." << std::endl;

    if (!filteredInputList.empty() && !m_filteredCaloHitListName.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredInputList, m_filteredCaloHitListName));

    if (!filteredCaloHitListU.empty() && !m_outputCaloHitListNameU.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListU, m_outputCaloHitListNameU));

    if (!filteredCaloHitListV.empty() && !m_outputCaloHitListNameV.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListV, m_outputCaloHitListNameV));

    if (!filteredCaloHitListW.empty() && !m_outputCaloHitListNameW.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListW, m_outputCaloHitListNameW));

     if (!filteredCaloHitListOp.empty() && !m_outputCaloHitListNameOp.empty())
        PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, filteredCaloHitListOp, m_outputCaloHitListNameOp));
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PreProcessingAlgorithm::ProcessChargeHit(const CaloHit *const pCaloHit, CaloHitList &hitList) const
{
    if (m_onlyAvailableCaloHits && !PandoraContentApi::IsAvailable(*this, pCaloHit))
        return;

    if (pCaloHit->GetMipEquivalentEnergy() < m_mipEquivalentCut)
        return;

    if (pCaloHit->GetInputEnergy() < std::numeric_limits<float>::epsilon())
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
            std::cout << "PreProcessingAlgorithm: found a hit with zero energy, will remove it" << std::endl;

        return;
    }

    if ((pCaloHit->GetCellLengthScale() < m_minCellLengthScale) || (pCaloHit->GetCellLengthScale() > m_maxCellLengthScale))
    {
        if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
        {
            std::cout << "PreProcessingAlgorithm: found a hit with extent " << pCaloHit->GetCellLengthScale() << ", require ("
                      << m_minCellLengthScale << " - " << m_maxCellLengthScale << "), will remove it" << std::endl;
        }

        return;
    }

    hitList.push_back(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PreProcessingAlgorithm::ProcessOpticalHit(const CaloHit *const pCaloHit, CaloHitList &hitList) const
{
    hitList.push_back(pCaloHit);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PreProcessingAlgorithm::PopulateVoidCaloHitLists() noexcept
{
    CaloHitList emptyList;

    try
    {
        if (!m_filteredCaloHitListName.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_filteredCaloHitListName));

        if (!m_outputCaloHitListNameU.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_outputCaloHitListNameU));

        if (!m_outputCaloHitListNameV.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_outputCaloHitListNameV));

        if (!m_outputCaloHitListNameW.empty())
            PANDORA_THROW_RESULT_IF(STATUS_CODE_SUCCESS, !=, PandoraContentApi::SaveList(*this, emptyList, m_outputCaloHitListNameW));
    }
    catch (...)
    {
        std::cout << "PreProcessingAlgorithm: Unable to create void calo hit lists" << std::endl;
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

void PreProcessingAlgorithm::GetFilteredCaloHitList(const CaloHitList &inputList, CaloHitList &outputList)
{
    HitKDTree2D kdTree;
    HitKDNode2DList hitKDNode2DList;

    KDTreeBox hitsBoundingRegion2D = fill_and_bound_2d_kd_tree(inputList, hitKDNode2DList);
    kdTree.build(hitKDNode2DList, hitsBoundingRegion2D);

    // Remove hits that are in the same physical location!
    for (const CaloHit *const pCaloHit1 : inputList)
    {
        bool isUnique(true);
        KDTreeBox searchRegionHits(build_2d_kd_search_region(pCaloHit1, m_searchRegion1D, m_searchRegion1D));

        HitKDNode2DList found;
        kdTree.search(searchRegionHits, found);

        for (const auto &hit : found)
        {
            const CaloHit *const pCaloHit2(hit.data);

            if (pCaloHit1 == pCaloHit2)
                continue;

            const float displacementSquared((pCaloHit2->GetPositionVector() - pCaloHit1->GetPositionVector()).GetMagnitudeSquared());

            if (displacementSquared < std::numeric_limits<float>::epsilon())
            {
                const float deltaMip(pCaloHit2->GetMipEquivalentEnergy() > pCaloHit1->GetMipEquivalentEnergy());

                if ((deltaMip > std::numeric_limits<float>::epsilon()) ||
                    ((std::fabs(deltaMip) < std::numeric_limits<float>::epsilon()) &&
                        (outputList.end() != std::find(outputList.begin(), outputList.end(), pCaloHit2))))
                {
                    isUnique = false;
                    break;
                }
            }
        }

        if (isUnique)
        {
            outputList.push_back(pCaloHit1);
        }
        else
        {
            if (PandoraContentApi::GetSettings(*this)->ShouldDisplayAlgorithmInfo())
                std::cout << "PreProcessingAlgorithm: found two hits in same location, will remove lowest pulse height" << std::endl;
        }
    }
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode PreProcessingAlgorithm::ReadSettings(const TiXmlHandle xmlHandle)
{
    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MipEquivalentCut", m_mipEquivalentCut));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MinCellLengthScale", m_minCellLengthScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxCellLengthScale", m_maxCellLengthScale));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "SearchRegion1D", m_searchRegion1D));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "MaxEventHits", m_maxEventHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OnlyAvailableCaloHits", m_onlyAvailableCaloHits));

    PANDORA_RETURN_RESULT_IF_AND_IF(
        STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "InputCaloHitListName", m_inputCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameU", m_outputCaloHitListNameU));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameV", m_outputCaloHitListNameV));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameW", m_outputCaloHitListNameW));

    PANDORA_RETURN_RESULT_IF_AND_IF(STATUS_CODE_SUCCESS, STATUS_CODE_NOT_FOUND, !=, XmlHelper::ReadValue(xmlHandle, "OutputCaloHitListNameOp", m_outputCaloHitListNameOp));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "FilteredCaloHitListName", m_filteredCaloHitListName));

    PANDORA_RETURN_RESULT_IF(STATUS_CODE_SUCCESS, !=, XmlHelper::ReadValue(xmlHandle, "CurrentCaloHitListReplacement", m_currentCaloHitListReplacement));

    return STATUS_CODE_SUCCESS;
}

} // namespace lar_content
