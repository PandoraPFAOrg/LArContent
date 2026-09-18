/**
 *  @file   larpandoracontent/LArObjects/LArPlaneContextObject.h
 *
 *  @brief  Header file for the LArPlaneContextObject class.
 *
 *  $Log: $
 */
#ifndef LAR_PLANE_CONTEXT_OBJECT_H
#define LAR_PLANE_CONTEXT_OBJECT_H 1

#include "Objects/CaloHit.h"
#include "Objects/EventContext.h"

#include <memory>

namespace lar_content
{

/**
 *  @brief  LArPlaneContextObject class
 */
class LArPlaneContextObject : public pandora::EventContextObject
{
public:
    struct HitTriplet
    {
        const pandora::CaloHit *m_uHit;
        const pandora::CaloHit *m_vHit;
        const pandora::CaloHit *m_wHit;
    };

    class HitTripletIterator
    {
    public:
        using InnerIterator = std::vector<std::unique_ptr<HitTriplet>>::const_iterator;

        explicit HitTripletIterator(InnerIterator it) :
            m_it(it)
        {
        }

        const HitTriplet *operator*() const
        {
            return m_it->get();
        }

        HitTripletIterator &operator++()
        {
            ++m_it;
            return *this;
        }

        bool operator!=(const HitTripletIterator &other) const
        {
            return m_it != other.m_it;
        }

    private:
        InnerIterator m_it;
    };

    /**
     *  @brief Adds a triplet of hits to the context object, where any one of the hits in this triplet could be null.
     *         The triplet is indexed by each of the hits it contains (if non-null) for retrieval. If more than one hit is null, the triplet
     *         is not added to the context object, as it cannot be valid.
     *
     *  @param  pUHit the U view hit in the triplet (can be null)
     *  @param  pVHit the V view hit in the triplet (can be null)
     *  @param  pWHit the W view hit in the triplet (can be null)
     *
     *  @return boolean indicating whether the triplet was successfully added to the context object (i.e. at most one hit is null)
     */
    bool AddHitTriplet(const pandora::CaloHit *const pUHit, const pandora::CaloHit *const pVHit, const pandora::CaloHit *const pWHit);

    /**
     *  @brief  Retrieve the triplet of hits associated with the specified hit. Note that any or all of the hits in this triplet could be null.
     *
     *  @param  pCaloHit the hit for which to retrieve the associated triplet
     *
     *  @return The triplet of hits associated with the specified hit.
     */
    const HitTriplet *GetHitTriplet(const pandora::CaloHit *const pCaloHit) const;

    /**
     *  @brief  Provides access to an iterator at the beginning of the hit triplets to allow iteration over the hit triplets in the context object.
     *
     *  @return A const iterator at the beginning of the triplets that can be used to iterate over the hit triplets in the context object.
     */
    HitTripletIterator begin() const;

    /**
     *  @brief  Provides access to an iterator at the end of the hit triplets to allow iteration over the hit triplets in the context object.
     *
     *  @return A const iterator at the end of the triplets that can be used to iterate over the hit triplets in the context object.
     */
    HitTripletIterator end() const;

    /**
     *  @brief  Retrieve the number of hit triplets.
     *
     *  @return The number of hit triplets.
     */
    size_t Size() const;

private:
    std::vector<std::unique_ptr<HitTriplet>> m_hitTriplets;              ///< Storage for all of the triplets of hits
    std::unordered_map<const pandora::CaloHit *, HitTriplet *> m_uIndex; ///< A map from each U hit to the triplet of hits associated with it
    std::unordered_map<const pandora::CaloHit *, HitTriplet *> m_vIndex; ///< A map from each V hit to the triplet of hits associated with it
    std::unordered_map<const pandora::CaloHit *, HitTriplet *> m_wIndex; ///< A map from each W hit to the triplet of hits associated with it
};

} // namespace lar_content

#endif // #ifndef LAR_PLANE_CONTEXT_OBJECT_H
