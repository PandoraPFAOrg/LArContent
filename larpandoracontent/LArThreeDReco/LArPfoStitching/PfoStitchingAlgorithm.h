/**
 *  @file   larpandoracontent/LArThreeDReco/LArPfoStitching/PfoStitchingAlgorithm.h
 *
 *  @brief  Header file for the track pfo stitching algorithm class.
 *
 *  $Log: $
 */
#ifndef LAR_PFO_STITCHING_ALGORITHM_H
#define LAR_PFO_STITCHING_ALGORITHM_H 1

#include "Api/PandoraApi.h"
#include "Geometry/LArTPC.h"
#include "Pandora/Algorithm.h"

#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingBaseTool.h"
#include "larpandoracontent/LArThreeDReco/LArPfoStitching/StitchingPfoOperations.h"

#include <memory>
#include <unordered_map>
#include <vector>

namespace lar_content
{

/**
 *  @brief  PfoStitchingAlgorithm class
 */
class PfoStitchingAlgorithm : public pandora::Algorithm, public StitchingPfoOperations
{
public:
    /**
     *  @brief  Lightweight LArTPC used only by worker-side stitching
     */
    class InferredLArTPC : public pandora::LArTPC
    {
    public:
        explicit InferredLArTPC(const PandoraApi::Geometry::LArTPC::Parameters &parameters);
        ~InferredLArTPC() override = default;
    };

    /**
     *  @brief  Default constructor
     */
    PfoStitchingAlgorithm();

    /**
     *  @brief  Shift a Pfo hierarchy by a specified x0 value
     *
     *  @param  pParentPfo the address of the parent pfo
     *  @param  pfoToLArTPCMap the pfo to lar tpc map
     *  @param  x0 the x0 correction relative to the input pfo
     */
    void ShiftPfoHierarchy(const pandora::ParticleFlowObject *const pParentPfo, const PfoToLArTPCMap &pfoToLArTPCMap, const float x0) const override;

    /**
     *  @brief  Stitch together a pair of pfos
     *
     *  @param  pPfoToEnlarge the address of the pfo to enlarge
     *  @param  pPfoToDelete the address of the pfo to delete
     *  @param  pfoToLArTPCMap the mapping between pfos and inferred tpcs
     */
    void StitchPfos(const pandora::ParticleFlowObject *const pPfoToEnlarge, const pandora::ParticleFlowObject *const pPfoToDelete,
        PfoToLArTPCMap &pfoToLArTPCMap) const override;

private:
    typedef std::vector<StitchingBaseTool *> StitchingToolVector;
    typedef std::vector<std::unique_ptr<InferredLArTPC>> InferredLArTPCVector;

    pandora::StatusCode Run() override;

    /**
     *  @brief  Infer physical TPC volumes by removing drift-gap intervals from the worker TPC
     *
     *  @param[out]  inferredTPCs the inferred TPC objects
     */
    pandora::StatusCode InferLArTPCs(InferredLArTPCVector &inferredTPCs) const;

    /**
     *  @brief  Build a mapping of track pfos whose 3D hits originate from one physical TPC volume
     *
     *  @param[in]   pPfoList the input pfo list
     *  @param[in]   inferredTPCs the inferred TPC objects
     *  @param[out]  pfoToLArTPCMap the map of to inferred tpc
     */
    void GetPfoToLArTPCMap(const pandora::PfoList *const pPfoList, const InferredLArTPCVector &inferredLArTPCs, PfoToLArTPCMap &pfoToLArTPCMap) const;

    /**
     *  @brief  Find the named list containing a cluster used by a pfo in the input list
     *
     *  @param[in]  pCluster the cluster whose list is required
     *
     *  @return the list name
     */
    const std::string GetListName(const pandora::Cluster *const pCluster) const;

    /**
     *  @brief  Find the named list containing a vertex used by a pfo in the input list
     *
     *  @param[in]  pVertex the vertex whose list is required
     *
     *  @return the list name
     */
    const std::string GetListName(const pandora::Vertex *const pVertex) const;

    /**
     *  @brief  Find the named list containing the input pfo
     *
     *  @param[in]  pPfo the input pfo whose list is required
     *
     *  @return the list name
     */
    const std::string GetListName(const pandora::ParticleFlowObject *const pPfo) const;

    /**
     *  @brief  Find the named list containing the input object
     *
     *  @param[in]  pObject the object whose list is required
     *  @param[in]  listNames the list names to search
     *
     *  @return the list name
     */
    template <typename TObject, typename TList>
    const std::string GetListName(const TObject *const pObject, const pandora::StringVector &listNames) const;

    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle) override;

    pandora::StringVector m_inputPfoListNames; ///< The input pfo list names
    pandora::StringVector m_daughterListNames; ///< Cluster and vertex lists used by the input pfos, updated when pfos are merged
    StitchingToolVector m_stitchingToolVector; ///< The stitching tools to be applied
    bool m_visualise;                          ///< Visualise merges
};

} // namespace lar_content

#endif // #ifndef LAR_PFO_STITCHING_ALGORITHM_H
