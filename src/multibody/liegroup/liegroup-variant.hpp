//
// Copyright (c) 2018 CNRS
//
// This file is part of Pinocchio
// Pinocchio is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// Pinocchio is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// Pinocchio If not, see
// <http://www.gnu.org/licenses/>.

#ifndef __se3_lie_group_variant_hpp__
#define __se3_lie_group_variant_hpp__

#include "pinocchio/multibody/liegroup/vector-space.hpp"
#include "pinocchio/multibody/liegroup/cartesian-product.hpp"
#include "pinocchio/multibody/liegroup/special-orthogonal.hpp"
#include "pinocchio/multibody/liegroup/special-euclidean.hpp"

#include <boost/variant.hpp>

namespace se3
{
  typedef boost::variant< SpecialOrthogonalOperationTpl<2>
                         ,SpecialOrthogonalOperationTpl<3>
                         ,SpecialEuclideanOperationTpl<2>
                         ,SpecialEuclideanOperationTpl<3>
                         ,VectorSpaceOperationTpl<1,double>
                         ,VectorSpaceOperationTpl<2,double>
                         ,VectorSpaceOperationTpl<3,double>
                         ,VectorSpaceOperationTpl<Eigen::Dynamic,double>
                        > LieGroupVariant;
  
}

#endif // ifndef __se3_lie_group_variant_hpp__
