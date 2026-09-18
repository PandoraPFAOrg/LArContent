/**
 *  @file   include/BasicVertexingMetrics.h
 *
 *  @brief  Header file for the Basic Vertexing Metrics.
 *
 *  $Log: $
 */
#include <iostream>
#ifndef LAR_DL_BASIC_VERTEXING_METRICS_H
#define LAR_DL_BASIC_VERTEXING_METRICS_H 1

#include "Pandora/Algorithm.h"
#include "Pandora/AlgorithmHeaders.h"

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

using namespace lar_content;

namespace lar_dl_content
{
/**
 *  @brief  BasicVertexingMetrics class
 */
class BasicVertexingMetrics : public pandora::Algorithm
{
public:
    /**
     *  @brief Default constructor
     */
    BasicVertexingMetrics();

    virtual ~BasicVertexingMetrics();

private:
    pandora::StatusCode Run();
    pandora::StatusCode ReadSettings(const pandora::TiXmlHandle xmlHandle);

    std::string m_caloHitListName;      ///< The name of the input calo hit list.
    std::string m_inputVertexListName;  ///< The name of the input vertex list

    std::string m_fileName;             ///< The name of the output ROOT file to write metrics to.
    std::string m_treeName;             ///< The name of the output ROOT tree to write
};

} // namespace lar_dl_content

#endif // LAR_DL_BASIC_VERTEXING_METRICS_H
