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

#ifndef __se3_spatial_hpp__
#define __se3_spatial_hpp__

#include "pinocchio/spatial/fwd.hpp"
#include <iostream>

namespace se3 {
  
  template <class Derived>
  class SpatialBase
  {
  public:
    
    typedef SpatialBase<Derived> Base;
    
    inline Derived & derived() { return *static_cast<Derived*> (this); }
    inline const Derived & derived() const { return *static_cast<const Derived*> (this); }
    inline Derived & const_cast_derived() const
    {
      return *static_cast<Derived*>(const_cast<SpatialBase*>(this));
    }
    
    template<typename SE3Derived>
    inline typename traits<Derived>::SE3ActionReturnType SE3ActOn(const SE3Base<SE3Derived> & M) const
    { return derived().SE3ActOn(M.derived()); }
    
    template<typename SE3Derived>
    inline typename traits<Derived>::SE3ActionReturnType SE3InvActOn(const SE3Base<SE3Derived> & M) const
    { return derived().SE3InvActOn(M.derived()); }
    
    inline void disp(std::ostream & os) const { derived().disp(os); }
    
    friend std::ostream & operator << (std::ostream & os, const SpatialBase<Derived> & X)
    { X.derived().disp(os); return os; }
    
  };
  
}

#endif // ifndef __se3_spatial_hpp__
