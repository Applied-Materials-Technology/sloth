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
class HSVStressUpdate : public RadialReturnStressUpdate
{
public:
  static InputParameters validParams();

  HSVStressUpdate(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void propagateQpStatefulProperties() override;

  virtual void computeStressInitialize(const Real & effective_trial_stress,
                                       const RankFourTensor & elasticity_tensor) override;
  virtual Real computeResidual(const Real & effective_trial_stress, const Real & scalar) override;
  virtual Real computeDerivative(const Real & effective_trial_stress, const Real & scalar) override;
  virtual void iterationFinalize(const Real & scalar) override;
  virtual void computeStressFinalize(const RankTwoTensor & plasticStrainIncrement) override;
  virtual Real computeHardeningValue(Real scalar);
  virtual Real computeHardeningDerivative(Real scalar); // Added. 

  /// a string to prepend to the plastic strain Material Property name
  const std::string _plastic_prepend;

  ///@{ Strain hardening parameters
  const MaterialProperty<Real> & _yield_stress; // Material property now
  const MaterialProperty<Real> & _sat_stress; // Hard coded voce now 
  const MaterialProperty<Real> & _exp_rate;
  const MaterialProperty<Real> & _lin_rate;
  Real _hardening_slope; // Added, stores current gradient of function
  ///@}

  ///@{ Viscoplasticity constitutive equation parameters
  const MaterialProperty<Real> & _c_alpha;
  const MaterialProperty<Real> & _c_beta;
  ///@}

  /// Elastic properties
  Real _yield_condition;

  ///@{ Viscoplasticity terms corresponding to Dunne and Petrinic eqn 5.64
  Real _xphir;
  Real _xphidp;
  ///@}

  MaterialProperty<Real> & _hardening_variable;
  const MaterialProperty<Real> & _hardening_variable_old;

  /// plastic strain of this model
  MaterialProperty<RankTwoTensor> & _plastic_strain;

  /// old value of plastic strain
  const MaterialProperty<RankTwoTensor> & _plastic_strain_old;
};