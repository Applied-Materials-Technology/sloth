#include "KRDamage.h"

registerMooseObject("SolidMechanicsApp", KRDamage);
registerMooseObject("SolidMechanicsApp", ADKRDamage);

template <bool is_ad>
InputParameters
KRDamageTempl<is_ad>::validParams()
{
  InputParameters params = ScalarDamageBaseTempl<is_ad>::validParams();
  params.addClassDescription(
      "Scalar damage model for which the damage is prescribed by another material");
  params.addRequiredParam<Real>("a", "Stress Scaling Parameter");
  params.addRequiredParam<Real>("phi", "Power for previous damage");
  params.addRequiredParam<Real>("zeta", "Stress power");

  return params;
}

template <bool is_ad>
KRDamageTempl<is_ad>::KRDamageTempl(const InputParameters & parameters)
  : ScalarDamageBaseTempl<is_ad>(parameters),
   _a(parameters.get<Real>("a")),
   _phi(parameters.get<Real>("phi")),
   _zeta(parameters.get<Real>("zeta")),
   _stress(this->template getGenericMaterialProperty<RankTwoTensor, is_ad>(_base_name + "stress"))
{
}

template <bool is_ad>
void
KRDamageTempl<is_ad>::updateQpDamageIndex()
{
  //_damage_index[_qp] = _damage_property[_qp];

  GenericRankTwoTensor<is_ad> devstress = _stress[_qp].deviatoric();
  GenericReal<is_ad> vmstress = std::sqrt((3/2)*devstress.doubleContraction(devstress));
  //std::cout << vmstress;
  //std::cout << " ";
  GenericReal<is_ad> damageIncrement = std::pow(vmstress/_a,_zeta)*std::pow(1-_damage_index_old[_qp],-_phi);
  
  _damage_index[_qp] = _damage_index_old[_qp] +_dt * damageIncrement;

  if (MooseUtils::absoluteFuzzyLessThan(_damage_index[_qp], 0.0) ||
      MooseUtils::absoluteFuzzyGreaterThan(_damage_index[_qp], 1.0))
    mooseError(_base_name + "damage_index ",
               "must be between 0 and 1. Current value is: ",
               _damage_index[_qp]);
  
}

template class KRDamageTempl<false>;
template class KRDamageTempl<true>;