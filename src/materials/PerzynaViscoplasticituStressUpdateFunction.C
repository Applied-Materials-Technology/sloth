#include "PerzynaViscoplasticityStressUpdateFunction.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("SolidMechanicsApp", PerzynaViscoplasticityStressUpdateFunction);
registerMooseObject("SolidMechanicsApp", ADPerzynaViscoplasticityStressUpdateFunction);

template <bool is_ad>
InputParameters
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::validParams()
{
  InputParameters params = RadialReturnStressUpdateTempl<is_ad>::validParams();
  params.addClassDescription("This class uses the discrete material for a Perzyna type "
                             "viscoplasticity model in which the effective plastic strain is "
                             "solved for using a creep approach.");

  // Non-linear function strain hardening parameters
  params.addRequiredParam<Real>("yield_stress",
                                "The point at which plastic strain begins accumulating");
  params.addRequiredParam<FunctionName>("hardening_function",
                                        "True Stress as a function of plastic stain");

  // Viscoplasticity constitutive equation parameters
  params.addRequiredParam<Real>("n", "Viscoplasticity coefficient, power law exponent");
  params.addRequiredParam<Real>("eta", "Viscoplasticity coefficient viscosity / drag stress");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

  return params;
}

template <bool is_ad>
PerzynaViscoplasticityStressUpdateFunctionTempl<
    is_ad>::PerzynaViscoplasticityStressUpdateFunctionTempl(const InputParameters & parameters)
  : RadialReturnStressUpdateTempl<is_ad>(parameters),
    _plastic_prepend(this->template getParam<std::string>("plastic_prepend")),
    _yield_stress(this->isParamValid("yield_stress") ? this->template getParam<Real>("yield_stress")
                                                     : 0),
    _hardening_function(this->isParamValid("hardening_function")
                            ? &this->getFunction("hardening_function")
                            : nullptr),
    _hardening_slope(0.0),
    _n(parameters.get<Real>("n")),
    _eta(parameters.get<Real>("eta")),
    _yield_condition(-1.0), // set to a non-physical value to catch uninitalized yield condition
    _hardening_variable(this->template declareGenericProperty<Real, is_ad>("hardening_variable")),
    _hardening_variable_old(this->template getMaterialPropertyOld<Real>("hardening_variable")),

    _plastic_strain(this->template declareGenericProperty<RankTwoTensor, is_ad>(
        _base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(this->template getMaterialPropertyOld<RankTwoTensor>(
        _base_name + _plastic_prepend + "plastic_strain"))
{
}

template <bool is_ad>
void
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::initQpStatefulProperties()
{
  _hardening_variable[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

template <bool is_ad>
void
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  RadialReturnStressUpdateTempl<is_ad>::propagateQpStatefulPropertiesRadialReturn();
}

template <bool is_ad>
void
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::computeStressInitialize(
    const GenericReal<is_ad> & effective_trial_stress,
    const GenericRankFourTensor<is_ad> & elasticity_tensor)
{
  RadialReturnStressUpdateTempl<is_ad>::computeStressInitialize(effective_trial_stress,
                                                                elasticity_tensor);

  _yield_condition = effective_trial_stress - _hardening_variable_old[_qp] - _yield_stress;

  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

template <bool is_ad>
GenericReal<is_ad>
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::computeResidual(
    const GenericReal<is_ad> & effective_trial_stress, const GenericReal<is_ad> & scalar)
{
  GenericReal<is_ad> residual = 0.0;

  mooseAssert(_yield_condition != -1.0,
              "the yield stress was not updated by computeStressInitialize");

  if (_yield_condition > 0.0)
  {
    _hardening_slope = computeHardeningDerivative(scalar);

    const GenericReal<is_ad> xflow = ((effective_trial_stress - (_three_shear_modulus * scalar)) /
                                      (computeHardeningValue(scalar) + _yield_stress)) -
                                     1;
    const GenericReal<is_ad> xphi = _eta * std::pow(xflow, _n);

    _xphidp =
        ((-_three_shear_modulus * _eta * _n) / (computeHardeningValue(scalar) + _yield_stress)) *
        std::pow(xflow, _n - 1);
    _xphir = -_eta * _n *
             ((effective_trial_stress - (_three_shear_modulus * scalar)) /
              std::pow((computeHardeningValue(scalar) + _yield_stress), 2)) *
             std::pow(xflow, _n - 1);
    residual = xphi * _dt - scalar;
  }
  return residual;
}

template <bool is_ad>
GenericReal<is_ad>
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::computeDerivative(
    const GenericReal<is_ad> & /*effective_trial_stress*/, const GenericReal<is_ad> & /*scalar*/)
{
  GenericReal<is_ad> derivative = 1.0;
  if (_yield_condition > 0.0)
    derivative = _xphidp * _dt + _hardening_slope * _xphir * _dt - 1.0;

  return derivative;
}

template <bool is_ad>
void
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::iterationFinalize(
    const GenericReal<is_ad> & scalar)
{
  if (_yield_condition > 0.0)
    _hardening_variable[_qp] = computeHardeningValue(scalar);
  // std::cout << _xphidp;
  // std::cout << ';';
  // std::cout << _xphir;
  // std::cout << '\n';

  // std::cout << scalar;
  // std::cout << ';';
  // std::cout << computeHardeningDerivative(scalar);
  // std::cout << ';';
  // std::cout << computeHardeningValue(scalar);
  // std::cout << '\n';
}

template <bool is_ad>
GenericReal<is_ad>
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::computeHardeningValue(
    const GenericReal<is_ad> & scalar)
{
  // return _hardening_variable_old[_qp] + (_hardening_slope * scalar);
  return _hardening_function->value(this->_effective_inelastic_strain_old[_qp] +
                                    scalar); // - _yield_stress; // Should be the main change gets
                                             // the hardening value using the function
}

template <bool is_ad>
GenericReal<is_ad>
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::computeHardeningDerivative(
    const GenericReal<is_ad> & scalar)
{
  return _hardening_function->timeDerivative(this->_effective_inelastic_strain_old[_qp]);
}

template <bool is_ad>
void
PerzynaViscoplasticityStressUpdateFunctionTempl<is_ad>::computeStressFinalize(
    const GenericRankTwoTensor<is_ad> & plasticStrainIncrement)
{
  _plastic_strain[_qp] += plasticStrainIncrement;
}

template class PerzynaViscoplasticityStressUpdateFunctionTempl<false>;
template class PerzynaViscoplasticityStressUpdateFunctionTempl<true>;