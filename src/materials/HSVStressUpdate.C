#include "HSVStressUpdate.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("SolidMechanicsApp", HSVStressUpdate);

InputParameters
HSVStressUpdate::validParams()
{
  InputParameters params = RadialReturnStressUpdate::validParams();
  params.addClassDescription("This class uses the discrete material for a hyperbolic sine "
                             "viscoplasticity model in which the effective plastic strain is "
                             "solved for using a creep approach. Voce model is hard coded in. Parameters should be material properties that vary spatially with"
                             "temperature.");

  // Non-linear Voce function strain hardening parameters
  params.addRequiredParam<MaterialPropertyName>("yield_stress",
                                "The point at which plastic strain begins accumulating");
  params.addRequiredParam<MaterialPropertyName>("sat_stress",
                                  "Saturation Stress of the Voce Model");
  params.addRequiredParam<MaterialPropertyName>("exp_rate",
                                  "Exponential rate of the saturation of the Voce Model");
  params.addRequiredParam<MaterialPropertyName>("lin_rate",
                                  "Linear stress increase of the Voce Model");
  // Viscoplasticity constitutive equation parameters
  params.addRequiredParam<MaterialPropertyName>("c_alpha",
                                "Viscoplasticity coefficient, scales the hyperbolic function");
  params.addRequiredParam<MaterialPropertyName>("c_beta",
                                "Viscoplasticity coefficient inside the hyperbolic sin function");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

  return params;
}

HSVStressUpdate::HSVStressUpdate(
    const InputParameters & parameters)
  : RadialReturnStressUpdate(parameters),
    _plastic_prepend(getParam<std::string>("plastic_prepend")),
    _yield_stress(getMaterialProperty<Real>("yield_stress")),
    _sat_stress(getMaterialProperty<Real>("sat_stress")),
    _exp_rate(getMaterialProperty<Real>("exp_rate")),
    _lin_rate(getMaterialProperty<Real>("lin_rate")),
    _hardening_slope(0.0),
    _c_alpha(getMaterialProperty<Real>("c_alpha")),
    _c_beta(getMaterialProperty<Real>("c_beta")),
    _yield_condition(-1.0), // set to a non-physical value to catch uninitalized yield condition
    _hardening_variable(declareProperty<Real>("hardening_variable")),
    _hardening_variable_old(getMaterialPropertyOld<Real>("hardening_variable")),

    _plastic_strain(
        declareProperty<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_strain")),
    _plastic_strain_old(
        getMaterialPropertyOld<RankTwoTensor>(_base_name + _plastic_prepend + "plastic_strain"))
{
}

void
HSVStressUpdate::initQpStatefulProperties()
{
  _hardening_variable[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

void
HSVStressUpdate::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  propagateQpStatefulPropertiesRadialReturn();
}

void
HSVStressUpdate::computeStressInitialize(
    const Real & effective_trial_stress, const RankFourTensor & elasticity_tensor)
{
  RadialReturnStressUpdate::computeStressInitialize(effective_trial_stress, elasticity_tensor);

  _yield_condition = effective_trial_stress - _hardening_variable_old[_qp] - _yield_stress[_qp];

  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

Real
HSVStressUpdate::computeResidual(const Real & effective_trial_stress,
                                                       const Real & scalar)
{
  Real residual = 0.0;

  mooseAssert(_yield_condition != -1.0,
              "the yield stress was not updated by computeStressInitialize");

  if (_yield_condition > 0.0)
  {
    _hardening_slope = computeHardeningDerivative(scalar);

    const Real xflow = _c_beta[_qp] * (effective_trial_stress - (_three_shear_modulus * scalar) -
                                  computeHardeningValue(scalar) - _yield_stress[_qp]);
    const Real xphi = _c_alpha[_qp] * std::sinh(xflow);

    _xphidp = -_three_shear_modulus * _c_alpha[_qp] * _c_beta[_qp] * std::cosh(xflow);
    _xphir = -_c_alpha[_qp] * _c_beta[_qp] * std::cosh(xflow);
    residual = xphi * _dt - scalar;
  }
  return residual;
}

Real
HSVStressUpdate::computeDerivative(const Real & /*effective_trial_stress*/,
                                                         const Real & /*scalar*/)
{
  Real derivative = 1.0;
  if (_yield_condition > 0.0)
    derivative = _xphidp * _dt + _hardening_slope * _xphir * _dt - 1.0;

  return derivative;
}

void
HSVStressUpdate::iterationFinalize(const Real & scalar)
{
  if (_yield_condition > 0.0)
    _hardening_variable[_qp] = computeHardeningValue(scalar);
    //std::cout << scalar;
    //std::cout << ';';
    //std::cout << computeHardeningDerivative(scalar);
    //std::cout << ';';
    //std::cout << computeHardeningValue(scalar);
    //std::cout << '\n';

}

Real
HSVStressUpdate::computeHardeningValue(Real scalar)
{
  // Voce model
  const Real current_strain = _effective_inelastic_strain_old[_qp] + scalar;
  return _sat_stress[_qp]*(1-std::exp(-1*_exp_rate[_qp]*current_strain))+_lin_rate[_qp]*current_strain;
}

Real
HSVStressUpdate::computeHardeningDerivative(Real scalar)
{

    return _lin_rate[_qp] + _sat_stress[_qp]*_exp_rate[_qp]*std::exp(-1*_exp_rate[_qp]*_effective_inelastic_strain_old[_qp]);
}

void
HSVStressUpdate::computeStressFinalize(
    const RankTwoTensor & plasticStrainIncrement)
{
  _plastic_strain[_qp] += plasticStrainIncrement;
}