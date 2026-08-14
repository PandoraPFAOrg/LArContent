/**
 *  @file   larpandoracontent/LArPersistency/LArTPCFactory.cc
 *
 *  @brief  Implementation of the lar tpc factory class.
 *
 *  $Log: $
 */

#include "Geometry/LArReadoutChannel.h"
#include "Geometry/LArReadoutUnit.h"
#include "Geometry/LArReadoutVolume.h"
#include "Geometry/LArTPC.h"

#include "larpandoracontent/LArPersistency/LArTPCFactory.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

using namespace pandora;

namespace lar_content
{

StatusCode LArTPCFactory::Write(const Object *const pObject, FieldMap &fields) const
{
    fields.Set("nReadoutVolumes", static_cast<unsigned int>(pObject->GetReadoutVolumes().size()));

    for (const auto &volumeEntry : pObject->GetReadoutVolumes())
    {
        const LArReadoutVolume &volume(volumeEntry.second);
        const std::string vPrefix("readoutVolume" + std::to_string(volume.GetId()) + "_");

        fields.Set(vPrefix + "center", volume.GetCenter());
        fields.Set(vPrefix + "size", volume.GetSize());
        fields.Set(vPrefix + "nUnits", static_cast<unsigned int>(volume.GetReadoutUnits().size()));

        for (const LArReadoutUnit &unit : volume.GetReadoutUnits())
        {
            const std::string uPrefix(vPrefix + "unit" + std::to_string(unit.GetId()) + "_");

            fields.Set(uPrefix + "view", unit.GetView());
            fields.Set(uPrefix + "nChannels", static_cast<unsigned int>(unit.GetReadoutChannels().size()));
            fields.Set(uPrefix + "referenceCoordinate", unit.GetReferenceCoordinate());
            fields.Set(uPrefix + "pitch", unit.GetPitch());
            fields.Set(uPrefix + "unitCenter", unit.GetUnitCenter());
            fields.Set(uPrefix + "unitSize", unit.GetUnitSize());
        }
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

StatusCode LArTPCFactory::Read(Parameters &parameters, const FieldMap &fields) const
{
    // ATTN Wire angles are read directly from the field map (rather than from parameters), since the caller populates the base LArTPC
    // fields on parameters only after this method returns.
    const float thetaU(fields.GetOrDefault<float>("wireAngleU", 0.f));
    const float thetaV(fields.GetOrDefault<float>("wireAngleV", 0.f));
    const float thetaW(fields.GetOrDefault<float>("wireAngleW", 0.f));

    const unsigned int nReadoutVolumes(fields.GetOrDefault<unsigned int>("nReadoutVolumes", 0u));

    for (unsigned int v = 0; v < nReadoutVolumes; ++v)
    {
        const std::string vPrefix("readoutVolume" + std::to_string(v) + "_");

        object_creation::LArReadoutVolumeParameters volumeParams;
        volumeParams.m_id = v;
        volumeParams.m_center = fields.GetOrDefault<CartesianVector>(vPrefix + "center", CartesianVector(0.f, 0.f, 0.f));
        volumeParams.m_size = fields.GetOrDefault<CartesianVector>(vPrefix + "size", CartesianVector(0.f, 0.f, 0.f));

        const unsigned int nUnits(fields.GetOrDefault<unsigned int>(vPrefix + "nUnits", 0u));

        UnitInfoVector unitInfoVector;
        for (unsigned int u = 0; u < nUnits; ++u)
        {
            const std::string uPrefix(vPrefix + "unit" + std::to_string(u) + "_");

            UnitInfo unitInfo;
            unitInfo.m_view = fields.GetOrDefault<HitType>(uPrefix + "view", TPC_VIEW_U);
            unitInfo.m_nChannels = fields.GetOrDefault<unsigned int>(uPrefix + "nChannels", 0u);
            unitInfo.m_referenceCoordinate = fields.GetOrDefault<float>(uPrefix + "referenceCoordinate", 0.f);
            unitInfo.m_pitch = fields.GetOrDefault<float>(uPrefix + "pitch", 1.f);
            unitInfo.m_unitCenter = fields.GetOrDefault<CartesianVector>(uPrefix + "unitCenter", CartesianVector(0.f, 0.f, 0.f));
            unitInfo.m_unitSize = fields.GetOrDefault<CartesianVector>(uPrefix + "unitSize", CartesianVector(0.f, 0.f, 0.f));
            unitInfoVector.push_back(unitInfo);
        }

        for (unsigned int u = 0; u < nUnits; ++u)
        {
            const UnitInfo &unitInfo(unitInfoVector.at(u));
            const float theta(LArTPCFactory::GetWireAngle(unitInfo.m_view, thetaU, thetaV, thetaW));

            object_creation::LArReadoutUnitParameters unitParams;
            unitParams.m_id = u;
            unitParams.m_view = unitInfo.m_view;
            unitParams.m_referenceCoordinate = unitInfo.m_referenceCoordinate;
            unitParams.m_pitch = unitInfo.m_pitch;
            unitParams.m_unitCenter = unitInfo.m_unitCenter;
            unitParams.m_unitSize = unitInfo.m_unitSize;

            // Regenerate the per-channel cross-view intervals for this readout unit
            for (unsigned int c = 0; c < unitInfo.m_nChannels; ++c)
            {
                const float selfCoordinate(unitInfo.m_referenceCoordinate + static_cast<float>(c) * unitInfo.m_pitch);

                CartesianVector point1(0.f, 0.f, 0.f), point2(0.f, 0.f, 0.f);
                LArTPCFactory::ClipLineAgainstBox(selfCoordinate, theta, unitInfo.m_unitCenter, unitInfo.m_unitSize, point1, point2);

                object_creation::LArReadoutChannelParameters channelParams;
                channelParams.m_id = c;
                std::size_t slot(0);

                for (unsigned int other = 0; other < nUnits; ++other)
                {
                    if (other == u)
                        continue;

                    const UnitInfo &otherInfo(unitInfoVector.at(other));
                    const float otherTheta(LArTPCFactory::GetWireAngle(otherInfo.m_view, thetaU, thetaV, thetaW));

                    const float coordinate1(LArTPCFactory::ProjectCoordinate(point1.GetY(), point1.GetZ(), otherTheta));
                    const float coordinate2(LArTPCFactory::ProjectCoordinate(point2.GetY(), point2.GetZ(), otherTheta));

                    const auto toIndex = [&](const float coordinate) -> unsigned int
                    {
                        const int rawIndex(static_cast<int>(std::round((coordinate - otherInfo.m_referenceCoordinate) / otherInfo.m_pitch)));
                        const int clampedIndex(std::max(0, std::min(rawIndex, static_cast<int>(otherInfo.m_nChannels) - 1)));
                        return static_cast<unsigned int>(clampedIndex);
                    };

                    const unsigned int index1(toIndex(coordinate1)), index2(toIndex(coordinate2));

                    channelParams.m_channelIntervalArray[slot++] = {otherInfo.m_view,
                        LArReadoutChannel::ChannelInterval{std::min(index1, index2), std::max(index1, index2)}};
                }

                unitParams.m_channelParametersVector.push_back(channelParams);
            }

            volumeParams.m_readoutUnitParametersVector.push_back(unitParams);
        }

        parameters.m_readoutVolumeParametersVector.push_back(volumeParams);
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArTPCFactory::ProjectCoordinate(const float y, const float z, const float thetaView)
{
    return z * std::cos(thetaView) - y * std::sin(thetaView);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArTPCFactory::ClipLineAgainstBox(const float coordinate, const float thetaView, const CartesianVector &boxCenter,
    const CartesianVector &boxSize, CartesianVector &point1, CartesianVector &point2)
{
    const float cosTheta(std::cos(thetaView)), sinTheta(std::sin(thetaView));
    const float dirY(cosTheta), dirZ(sinTheta);

    const float y0(boxCenter.GetY());
    const float z0((std::fabs(cosTheta) > 1e-6f) ? (coordinate + y0 * sinTheta) / cosTheta : boxCenter.GetZ());

    const float yMin(boxCenter.GetY() - 0.5f * boxSize.GetY()), yMax(boxCenter.GetY() + 0.5f * boxSize.GetY());
    const float zMin(boxCenter.GetZ() - 0.5f * boxSize.GetZ()), zMax(boxCenter.GetZ() + 0.5f * boxSize.GetZ());

    float tMin(-std::numeric_limits<float>::max()), tMax(std::numeric_limits<float>::max());

    const auto clip = [&](const float p, const float d, const float lo, const float hi)
    {
        if (std::fabs(d) < 1e-9f)
            return;

        const float t0((lo - p) / d), t1((hi - p) / d);
        tMin = std::max(tMin, std::min(t0, t1));
        tMax = std::min(tMax, std::max(t0, t1));
    };

    clip(y0, dirY, yMin, yMax);
    clip(z0, dirZ, zMin, zMax);

    point1 = CartesianVector(0.f, y0 + tMin * dirY, z0 + tMin * dirZ);
    point2 = CartesianVector(0.f, y0 + tMax * dirY, z0 + tMax * dirZ);
}

//------------------------------------------------------------------------------------------------------------------------------------------

float LArTPCFactory::GetWireAngle(const HitType view, const float thetaU, const float thetaV, const float thetaW)
{
    switch (view)
    {
        case TPC_VIEW_U:
            return thetaU;
        case TPC_VIEW_V:
            return thetaV;
        case TPC_VIEW_W:
        default:
            return thetaW;
    }
}

} // namespace lar_content
