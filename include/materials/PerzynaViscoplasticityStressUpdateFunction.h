#pragma once

#include "RadialReturnStressUpdate.h"

/**
 * This class uses the Discrete material in an isotropic radial return hyperbolic
 * sine viscoplasticity model.
 *
 * This class inherits from RadialReturnStressUpdate and must be used
 * in conjunction with ComputeReturnMappingStress. This uniaxial viscoplasticity
 * class computes the plastic strain as a stateful material property.  The
 * constitutive equation for scalar plastic strain rate used in this model is
 * /f$ \dot{p} = \phi (\sigma_e , r) = \alpha sinh \beta (\sigma_e -r - \sigma_y) f/$
 *
 * This class is based on the implicit integration algorithm in F. Dunne and N.
 * Petrinic's Introduction to Computational Plasticity (2004) Oxford University
 * Press, pg. 162 - 163.
 *
 * Expanded to permit function-based strain hardening rather than linear hardening.
 * Basing chages on IsotropicPlasticityStressUpdate approach, which also derives
 * from Dunne & Petrinic.
 */
template <bool is_ad>
class PerzynaViscoplasticityStressUpdateFunctionTempl : public RadialReturnStressUpdateTempl<is_ad>
{
public:
  static InputParameters validParams();

  PerzynaViscoplasticityStressUpdateFunctionTempl(const InputParameters & parameters);

  using Material::_qp;
  using RadialReturnStressUpdateTempl<is_ad>::_base_name;
  using RadialReturnStressUpdateTempl<is_ad>::_three_shear_modulus;
  using RadialReturnStressUpdateTempl<is_ad>::_dt;

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  virtual void
  computeStressInitialize(const GenericReal<is_ad> & effective_trial_stress,
                          const GenericRankFourTensor<is_ad> & elasticity_tensor) override;
  virtual GenericReal<is_ad> computeResidual(const GenericReal<is_ad> & effective_trial_stress,
                                             const GenericReal<is_ad> & scalar) override;
  virtual GenericReal<is_ad> computeDerivative(const GenericReal<is_ad> & effective_trial_stress,
                                               const GenericReal<is_ad> & scalar) override;
  virtual void iterationFinalize(const GenericReal<is_ad> & scalar) override;
  virtual void
  computeStressFinalize(const GenericRankTwoTensor<is_ad> & plasticStrainIncrement) override;
  virtual GenericReal<is_ad> computeHardeningValue(const GenericReal<is_ad> & scalar);
  virtual GenericReal<is_ad> computeHardeningDerivative(const GenericReal<is_ad> & scalar);

  /// a string to prepend to the plastic strain Material Property name
  const std::string _plastic_prepend;

  ///@{ Strain hardening parameters
  const Real _yield_stress; // Constant for now, but in IsotropicPlasticityStressUpdate it can be a
  // function of temperature.
  const Function * const _hardening_function; // Modified to function.
  GenericReal<is_ad> _hardening_slope;        // Added, stores current gradient of function
  ///@}

  ///@{ Viscoplasticity constitutive equation parameters
  const Real _n;
  const Real _eta;
  ///@}

  /// Elastic properties
  GenericReal<is_ad> _yield_condition;

  ///@{ Viscoplasticity terms corresponding to Dunne and Petrinic eqn 5.64
  GenericReal<is_ad> _xphir;
  GenericReal<is_ad> _xphidp;
  ///@}

  GenericMaterialProperty<Real, is_ad> & _hardening_variable;
  const MaterialProperty<Real> & _hardening_variable_old;

  /// plastic strain of this model
  GenericMaterialProperty<RankTwoTensor, is_ad> & _plastic_strain;

  /// old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  std::vector<GenericMaterialProperty<Real, is_ad> *> _properties;
};

typedef PerzynaViscoplasticityStressUpdateFunctionTempl<false>
    PerzynaViscoplasticityStressUpdateFunction;
typedef PerzynaViscoplasticityStressUpdateFunctionTempl<true>
    ADPerzynaViscoplasticityStressUpdateFunction;