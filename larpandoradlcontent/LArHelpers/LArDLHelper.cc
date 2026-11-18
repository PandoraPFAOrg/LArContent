/**
 *  @file   larpandoradlcontent/LArHelpers/LArDLHelper.cc
 *
 *  @brief  Implementation of the lar deep learning helper helper class.
 *
 *  $Log: $
 */

#include "larpandoradlcontent/LArHelpers/LArDLHelper.h"

namespace lar_dl_content
{

using namespace pandora;

StatusCode LArDLHelper::LoadModel(const std::string &filename, LArDLHelper::TorchModel &model)
{
    try
    {
        model = torch::jit::load(filename);
        std::cout << "Loaded the TorchScript model \'" << filename << "\'" << std::endl;

        // Set the model to evaluation mode.
        // This should have been done during the model export, but we do it here just in case.
        // This ensures that layers like dropout and batch normalization behave correctly during inference.
        model.eval();
    }
    catch (const std::exception &e)
    {
        std::cout << "Error loading the TorchScript model \'" << filename << "\':\n" << e.what() << std::endl;
        return STATUS_CODE_FAILURE;
    }

    return STATUS_CODE_SUCCESS;
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::InitialiseInput(const at::IntArrayRef dimensions, TorchInput &tensor)
{
    tensor = torch::zeros(dimensions);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::InitialiseInput(const at::IntArrayRef dimensions, TorchInput &tensor, const torch::TensorOptions &options)
{
    tensor = torch::zeros(dimensions, options);
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::Forward(TorchModel &model, const TorchInputVector &input, TorchOutput &output)
{
    // Set torch to no_grad mode to avoid tracking gradients, which are not
    // needed during inference.
    // This uses RAII, so the guard is only active within this scope.
    torch::NoGradGuard guard;

    output = model.forward(input).toTensor();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::Forward(TorchModel &model, const TorchInputVector &input, TorchOutput &output, const std::string &funcName)
{
    // Set torch to no_grad mode to avoid tracking gradients, which are not
    // needed during inference.
    // This uses RAII, so the guard is only active within this scope.
    torch::NoGradGuard guard;

    output = model.get_method(funcName)(input).toTensor();
}

//------------------------------------------------------------------------------------------------------------------------------------------

void LArDLHelper::Forward(TorchModel &model, const TorchInputVector &input, TorchMultiOutput &output)
{
    // Set torch to no_grad mode to avoid tracking gradients, which are not
    // needed during inference.
    // This uses RAII, so the guard is only active within this scope.
    torch::NoGradGuard guard;

    output = model.forward(input);
}

} // namespace lar_dl_content
