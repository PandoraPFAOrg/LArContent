/**
 *  @file   larpandoradlcontent/LArHelpers/LArDLHelper.h
 *
 *  @brief  Header file for the lar deep learning helper helper class.
 *
 *  $Log: $
 */
#ifndef LAR_DL_HELPER_H
#define LAR_DL_HELPER_H 1

#include <torch/script.h>
#include <torch/torch.h>

#include "Pandora/StatusCodes.h"

namespace lar_dl_content
{

/**
 *  @brief  LArDLHelper class
 */
class LArDLHelper
{
public:
    typedef torch::jit::script::Module TorchModel;
    typedef torch::Tensor TorchInput;
    typedef std::vector<torch::jit::IValue> TorchInputVector;
    typedef at::Tensor TorchOutput;
    typedef at::IValue TorchMultiOutput;
    typedef c10::impl::GenericDict TorchDict;

    /**
     *  @brief  Loads a deep learning model
     *
     *  @param  filename the filename of the model to load
     *  @param  model the TorchModel in which to store the loaded model
     *
     *  @return STATUS_CODE_SUCCESS upon successful loading of the model. STATUS_CODE_FAILURE otherwise.
     */
    static pandora::StatusCode LoadModel(const std::string &filename, TorchModel &model);

    /**
     *  @brief  Create a torch input tensor
     *
     *  @param  dimensions the size of each dimension of the tensor: pass as {a, b, c, d} for example
     *  @param  tensor the tensor to be initialised
     */
    static void InitialiseInput(const at::IntArrayRef dimensions, TorchInput &tensor);

    /**
     *  @brief  Creates a output dictionary with string keys and tensor values
     *
     *  @return the created output dictionary
     */
    static TorchDict CreateOutputDict();

    /**
     *  @brief  Run a deep learning model
     *
     *  @param  model the model to run
     *  @param  input the input to run over
     *  @param  output the tensor to store the output in
     */
    static void Forward(TorchModel &model, const TorchInputVector &input, TorchOutput &output);

    /**
     *  @brief  Run a deep learning model
     *
     *  @param  model the model to run
     *  @param  input the input to run over
     *  @param  output the at::IValue to store the output in
     */
    static void Forward(TorchModel &model, const TorchInputVector &input, TorchMultiOutput &output);

    /**
     *  @brief  Run a deep learning model
     *
     *  @param  model the model to run
     *  @param  input the input to run over
     *  @param  output the dictionary to store the output in
     */
    static void Forward(TorchModel &model, const TorchInputVector &input, TorchDict &output);
};

} // namespace lar_dl_content

#endif // #ifndef LAR_DL_HELPER_H
