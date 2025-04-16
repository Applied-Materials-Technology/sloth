#pragma once

#include "ScalarDamageBase.h"

template <bool is_ad>
class KRDamageTempl : public ScalarDamageBaseTempl<is_ad>
{
public:
  static InputParameters validParams();

  KRDamageTempl(const InputParameters & parameters);

protected:
  virtual void updateQpDamageIndex() override;

  ///@{ Material properties for the damage model.
  const Real _a;
  const Real _phi;
  const Real _zeta;

  ///@}

  using Material::_current_elem;
  using Material::_dt;
  using Material::_q_point;
  using Material::_qp;

  const GenericMaterialProperty<RankTwoTensor, is_ad> & _stress;

  using ScalarDamageBaseTempl<is_ad>::_damage_index;
  using ScalarDamageBaseTempl<is_ad>::_damage_index_old;
  using ScalarDamageBaseTempl<is_ad>::_base_name;
  using ScalarDamageBaseTempl<is_ad>::_damage_index_name;
};

typedef KRDamageTempl<false> KRDamage;
typedef KRDamageTempl<true> ADKRDamage;