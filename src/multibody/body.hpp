//
// Copyright (c) 2016 CNRS
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

#ifndef __se3_body_hpp__
#define __se3_body_hpp__

#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/spatial/inertia.hpp"
#include "pinocchio/multibody/fwd.hpp"
#include "pinocchio/tools/string-generator.hpp"

#include <Eigen/StdVector>
#include <iostream>
#include <vector>

namespace se3 {
  
  template<typename _Scalar, int _Options>
  struct BodyTpl
  {
    typedef SE3Tpl<_Scalar, _Options> SE3;
    typedef InertiaTpl<_Scalar, _Options> Inertia;
    
    typedef se3::JointIndex JointIndex;
    
    ///
    /// \brief Default constructor
    ///
    BodyTpl()
      : parent_id(0)
      , bodyPlacement(1)
      , Y(Inertia::Zero())
      , name("body_" + random(8))
    {}
    
    BodyTpl(const JointIndex parent_id,
            const SE3 & bodyPlacement,
            const Inertia & Y,
            const std::string & body_name = "")
      : parent_id(parent_id)
      , bodyPlacement(bodyPlacement)
      , Y(Y)
      , visual_geometry_paths(0)
    {
      name = (body_name!="")?body_name:random(8);
    }
    
    /// \brief Index of the parent joint.
    JointIndex parent_id;
    
    /// \brief bodyPlacement regarding to the supporting joint.
    SE3 bodyPlacement;
    
    /// \brief Spatial inertia characterizing the body.
    Inertia Y;
    
    /// \brief Name of the body
    std::string name;
    
    /// \brief List of visuals describing the body
    std::vector <std::string> visual_geometry_paths;
    
  };
  
  typedef BodyTpl<double,0> Body;
  
} // namespace se3

EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(se3::Body)

#endif // ifndef __se3_body_hpp__
