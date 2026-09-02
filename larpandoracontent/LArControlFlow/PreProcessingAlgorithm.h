/**
 *  @file   larpandoracontent/LArControlFlow/PreProcessingAlgorithm.h
 *
 *  @brief  Header file for the pre processing algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PRE_PROCESSING_ALGORITHM_H
#define LAR_PRE_PROCESSING_ALGORITHM_H 1

#include "Pandora/Algorithm.h"

namespace lar_content
{

template <typename, unsigned int>
class KDTreeLinkerAlgo;
template <typename, unsigned int>
class KDTreeNodeInfoT;

//------------------------------------------------------------------------------------------------------------------------------------------

/**
 *  @brief  PreProcessingAlgorithm class
 */
class PreProcessingAlgorithm : public pandora::Algorithm
{
public:
    /**
     *  @brief  Default constructor
     */
    PreProcessingAlgorithm();

private:
    typedef KDTreeLinkerAlgo<const pandora::CaloHit *, 2> HitKDTree2D;
    typedef KDTreeNodeInfoT<const pandora::CaloHit *, 2> HitKDNode2D;
    typedef std::vector<HitKDNode2D> HitKDNode2DList;

    pandora::StatusCode Reset();
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    /**
     *  @brief Build separate CaloHitLists for charge and light, separating charge hits by view
     */
    void ProcessCaloHits();

    /**
     *  @brief Build separate CaloHitLists for charge hits for each view
     */
    void ProcessChargeHit(const pandora::CaloHit *const pCaloHit, pandora::CaloHitList &hitList) const;

    /**
     *  @brief Build CaloHitList for optical hits
     */
    void ProcessOpticalHit(const pandora::CaloHit *const pCaloHit, pandora::CaloHitList &hitList) const;

    /**
     *  @brief Build empty calo hit lists
     */
    void PopulateVoidCaloHitLists() noexcept;

    /**
     *  @brief Clean up the input CaloHitList
     *
     *  @param inputList the input CaloHitList
     *  @param outputList the output CaloHitList
     */
    void GetFilteredCaloHitList(const pandora::CaloHitList &inputList, pandora::CaloHitList &outputList);

    /**
     *  @brief Build separate MCParticleLists for each view
     */
    void ProcessMCParticles();

    pandora::CaloHitSet m_processedHits; ///< The set of all previously processed calo hits

    float m_mipEquivalentCut;    ///< Minimum mip equivalent energy for calo hit
    float m_minCellLengthScale;  ///< The minimum length scale for calo hit
    float m_maxCellLengthScale;  ///< The maximum length scale for calo hit
    float m_searchRegion1D;      ///< Search region, applied to each dimension, for look-up from kd-trees
    unsigned int m_maxEventHits; ///< The maximum number of hits in an event to proceed with the reconstruction

    bool m_onlyAvailableCaloHits;                ///< Whether to only include available calo hits
    std::string m_inputCaloHitListName;          ///< The input calo hit list name
    std::string m_outputCaloHitListNameU;        ///< The output calo hit list name for TPC_VIEW_U hits
    std::string m_outputCaloHitListNameV;        ///< The output calo hit list name for TPC_VIEW_V hits
    std::string m_outputCaloHitListNameW;        ///< The output calo hit list name for TPC_VIEW_W hits
    std::string m_outputCaloHitListNameOp;       ///< The output calo hit list name for optical hits
    std::string m_filteredCaloHitListName;       ///< The output calo hit list name for all U, V and W hits
    std::string m_currentCaloHitListReplacement; ///< The name of the calo hit list to replace the current list (optional)
};

} // namespace lar_content

#endif // #ifndef LAR_PRE_PROCESSING_ALGORITHM_H
