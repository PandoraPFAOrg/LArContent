/**
 *  @file   larpandoracontent/LArPersistency/LArTPCFactory.h
 *
 *  @brief  Header file for the lar tpc factory class, responsible for persisting the readout volume/unit/channel
 *          geometry associated with a LArTPC.
 *
 *  $Log: $
 */
#ifndef LAR_TPC_FACTORY_H
#define LAR_TPC_FACTORY_H 1

#include "Objects/CartesianVector.h"

#include "Pandora/ObjectCreation.h"
#include "Pandora/PandoraObjectFactories.h"

#include "Persistency/FieldMap.h"

namespace lar_content
{

/**
 *  @brief  LArTPCFactory class. Responsible for persisting and reconstructing the readout volume/unit/channel geometry associated with a
 *          LArTPC.
 *
 *          Only the generative parameters of each readout unit are persisted (view, channel count, a reference coordinate for channel 0,
 *          and the signed pitch between adjacent channels); the per-channel cross-view intervals are recomputed on read, rather than being
 *          persisted to keep geometry files compact.
 */
class LArTPCFactory : public pandora::PandoraObjectFactory<object_creation::Geometry::LArTPC::Parameters, object_creation::Geometry::LArTPC::Object>
{
public:
    /**
     *  @brief  Read the readout volume/unit/channel parameters from the supplied field map, regenerating the per-channel cross-view
     *          intervals from the persisted generative parameters.
     *
     *  @param  parameters the parameters to populate
     *  @param  fields the field map from which to read the parameters
     */
    pandora::StatusCode Read(Parameters &parameters, const pandora::FieldMap &fields) const;

    /**
     *  @brief  Write the readout volume/unit/channel parameters into the supplied field map.
     *
     *  @param  pObject the address of the lar tpc to persist
     *  @param  fields the field map into which to write the parameters
     */
    pandora::StatusCode Write(const Object *const pObject, pandora::FieldMap &fields) const;

    /**
     *  @brief  Compute the cross-view channel interval for a single channel against a single other view, from the generative
     *          parameters of both readout units. Used both when regenerating persisted intervals on read, and by LArPandora when
     *          computing intervals live, so both paths agree exactly.
     *
     *  @param  selfCoordinate the wire coordinate of the channel under consideration
     *  @param  thetaSelf the wire angle of the channel's own view
     *  @param  selfUnitCenter the centre of the channel's own readout unit active-area box (X unused)
     *  @param  selfUnitSize the extent of the channel's own readout unit active-area box (X unused)
     *  @param  thetaOther the wire angle (to the vertical) of the other view
     *  @param  otherReferenceCoordinate the reference coordinate of the other readout unit
     *  @param  otherPitch the signed pitch of the other readout unit
     *  @param  otherNChannels the number of channels in the other readout unit
     *
     *  @return the channel interval in the other view
     */
    static pandora::LArReadoutChannel::ChannelInterval ComputeChannelInterval(const float selfCoordinate, const float thetaSelf,
        const pandora::CartesianVector &selfUnitCenter, const pandora::CartesianVector &selfUnitSize, const float thetaOther,
        const float otherReferenceCoordinate, const float otherPitch, const unsigned int otherNChannels);

private:
    /**
     *  @brief  Collection of per-readout-unit information required while regenerating channel intervals.
     */
    class UnitInfo
    {
    public:
        pandora::HitType m_view;     ///< The view of the readout unit
        unsigned int m_nChannels;    ///< The number of channels in the readout unit
        float m_referenceCoordinate; ///< The wire-coordinate of channel 0's midpoint
        float m_pitch;               ///< The signed wire-coordinate difference between adjacent channels
        pandora::CartesianVector m_unitCenter{0.f, 0.f, 0.f}; ///< The center of the readout unit's own active-area box (X unused)
        pandora::CartesianVector m_unitSize{0.f, 0.f, 0.f};   ///< The size of the readout unit's own active-area box (X unused)
    };

    typedef std::vector<UnitInfo> UnitInfoVector;

    /**
     *  @brief  Project a point in the Y-Z plane onto a view's wire-coordinate axis.
     *
     *  @param  y the y coordinate
     *  @param  z the z coordinate
     *  @param  thetaView the wire angle (to the vertical) of the view being projected onto
     *
     *  @return the wire coordinate
     */
    static float ProjectCoordinate(const float y, const float z, const float thetaView);

    /**
     *  @brief  Find the two points at which the line of constant wire-coordinate (for the given view) crosses the Y-Z bounding box of a
     *          readout volume.
     *
     *  @param  coordinate the wire coordinate defining the line
     *  @param  thetaView the wire angle (to the vertical) of the view the line belongs to
     *  @param  boxCenter the centre of the readout volume's bounding box
     *  @param  boxSize the extent of the readout volume's bounding box
     *  @param  point1 to receive the first crossing point
     *  @param  point2 to receive the second crossing point
     */
    static void ClipLineAgainstBox(const float coordinate, const float thetaView, const pandora::CartesianVector &boxCenter,
        const pandora::CartesianVector &boxSize, pandora::CartesianVector &point1, pandora::CartesianVector &point2);

    /**
     *  @brief  Get the wire angle associated with a given view.
     *
     *  @param  view the view
     *  @param  thetaU the u wire angle
     *  @param  thetaV the v wire angle
     *  @param  thetaW the w wire angle
     *
     *  @return the wire angle associated with the specified view
     */
    static float GetWireAngle(const pandora::HitType view, const float thetaU, const float thetaV, const float thetaW);
};

} // namespace lar_content

#endif // #ifndef LAR_TPC_FACTORY_H
