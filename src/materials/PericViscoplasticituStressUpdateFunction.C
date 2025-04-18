#include "PericViscoplasticityStressUpdateFunction.h"

#include "Function.h"
#include "ElasticityTensorTools.h"

registerMooseObject("SolidMechanicsApp", PericViscoplasticityStressUpdateFunction);

InputParameters
PericViscoplasticityStressUpdateFunction::validParams()
{
  InputParameters params = RadialReturnStressUpdate::validParams();
  params.addClassDescription("This class uses the discrete material for a Peric type "
                             "viscoplasticity model in which the effective plastic strain is "
                             "solved for using a creep approach.");

  // Non-linear function strain hardening parameters
  params.addRequiredParam<Real>("yield_stress",
                                "The point at which plastic strain begins accumulating");
  params.addRequiredParam<FunctionName>("hardening_function","True Stress as a function of plastic stain");

  // Viscoplasticity constitutive equation parameters
  params.addRequiredParam<Real>("n",
                                "Viscoplasticity coefficient, power law exponent");
  params.addRequiredParam<Real>("eta",
                                "Viscoplasticity coefficient viscosity / drag stress");
  params.addDeprecatedParam<std::string>(
      "plastic_prepend",
      "",
      "String that is prepended to the plastic_strain Material Property",
      "This has been replaced by the 'base_name' parameter");
  params.set<std::string>("effective_inelastic_strain_name") = "effective_plastic_strain";

  return params;
}

PericViscoplasticityStressUpdateFunction::PericViscoplasticityStressUpdateFunction(
    const InputParameters & parameters)
  : RadialReturnStressUpdate(parameters),
    _plastic_prepend(getParam<std::string>("plastic_prepend")),
    _yield_stress(parameters.get<Real>("yield_stress")),
    _hardening_function(getFunction("hardening_function")), // Added but no clue if this will work as other class is based on templates.
    _hardening_slope(0.0),
    _n(parameters.get<Real>("n")),
    _eta(parameters.get<Real>("eta")),
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
PericViscoplasticityStressUpdateFunction::initQpStatefulProperties()
{
  _hardening_variable[_qp] = 0.0;
  _plastic_strain[_qp].zero();
}

void
PericViscoplasticityStressUpdateFunction::propagateQpStatefulProperties()
{
  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];

  propagateQpStatefulPropertiesRadialReturn();
}

void
PericViscoplasticityStressUpdateFunction::computeStressInitialize(
    const Real & effective_trial_stress, const RankFourTensor & elasticity_tensor)
{
  RadialReturnStressUpdate::computeStressInitialize(effective_trial_stress, elasticity_tensor);

  _yield_condition = effective_trial_stress - _hardening_variable_old[_qp] - _yield_stress;

  _hardening_variable[_qp] = _hardening_variable_old[_qp];
  _plastic_strain[_qp] = _plastic_strain_old[_qp];
}

Real
PericViscoplasticityStressUpdateFunction::computeResidual(const Real & effective_trial_stress,
                                                       const Real & scalar)
{
  Real residual = 0.0;

  mooseAssert(_yield_condition != -1.0,
              "the yield stress was not updated by computeStressInitialize");

  if (_yield_condition > 0.0)
  {
    _hardening_slope = computeHardeningDerivative(scalar);

    const Real xflow =  ((effective_trial_stress - (_three_shear_modulus * scalar))/
                                  (computeHardeningValue(scalar) + _yield_stress));
    const Real xphi = _eta*(std::pow(xflow,_n)-1);

    _xphidp = ((-_three_shear_modulus * _eta * _n)/(computeHardeningValue(scalar) + _yield_stress))*std::pow(xflow,_n-1);
    _xphir = -_eta * _n*((effective_trial_stress - (_three_shear_modulus * scalar))/std::pow((computeHardeningValue(scalar) + _yield_stress),2)) * std::pow(xflow,_n-1);
    residual = xphi * _dt - scalar;
  }
  return residual;
}

Real
PericViscoplasticityStressUpdateFunction::computeDerivative(const Real & /*effective_trial_stress*/,
                                                         const Real & /*scalar*/)
{
  Real derivative = 1.0;
  if (_yield_condition > 0.0)
    derivative = _xphidp * _dt + _hardening_slope * _xphir * _dt - 1.0;

  return derivative;
}

void
PericViscoplasticityStressUpdateFunction::iterationFinalize(const Real & scalar)
{
  if (_yield_condition > 0.0)
    _hardening_variable[_qp] = computeHardeningValue(scalar);
    //std::cout << _xphidp;
    //std::cout << ';';
    //std::cout << _xphir;
    //std::cout << '\n';

    //std::cout << scalar;
    //std::cout << ';';
    //std::cout << computeHardeningDerivative(scalar);
    //std::cout << ';';
    //std::cout << computeHardeningValue(scalar);
    //std::cout << '\n';

}

Real
PericViscoplasticityStressUpdateFunction::computeHardeningValue(Real scalar)
{
  //return _hardening_variable_old[_qp] + (_hardening_slope * scalar);
  return _hardening_function.value(_effective_inelastic_strain_old[_qp] + scalar,_q_point[_qp]);// - _yield_stress; // Should be the main change gets the hardening value using the function 
}

Real
PericViscoplasticityStressUpdateFunction::computeHardeningDerivative(Real scalar)
{
    return _hardening_function.timeDerivative(_effective_inelastic_strain_old[_qp],_q_point[_qp]);
}

void
PericViscoplasticityStressUpdateFunction::computeStressFinalize(
    const RankTwoTensor & plasticStrainIncrement)
{
  _plastic_strain[_qp] += plasticStrainIncrement;
}