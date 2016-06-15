//
// Copyright (c) 2015-2016 CNRS
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

#include "pinocchio/spatial/fwd.hpp"
#include "pinocchio/spatial/se3.hpp"
#include "pinocchio/multibody/joint.hpp"
#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/multibody/model.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/rnea.hpp"
#include "pinocchio/multibody/parser/urdf.hpp"
#include "pinocchio/multibody/parser/sample-models.hpp"

#include <iostream>

#include "pinocchio/tools/timer.hpp"

#include <Eigen/StdVector>
EIGEN_DEFINE_STL_VECTOR_SPECIALIZATION(Eigen::VectorXd)


int main(int argc, const char ** argv)
{
  using namespace Eigen;
  using namespace se3;
  
  _MM_SET_DENORMALS_ZERO_MODE(_MM_DENORMALS_ZERO_ON);
  
  StackTicToc timer(StackTicToc::US);
#ifdef NDEBUG
  const int NBT = 100000;
#else
  const int NBT = 10;
  std::cout << "(the time score in debug mode is not relevant) " << std::endl;
#endif
  
  se3::Model model;
  
  std::string filename = PINOCCHIO_SOURCE_DIR"/models/simple_humanoid.urdf";
  if(argc>1) filename = argv[1];
  if( filename == "HS")
    se3::buildModels::humanoidSimple(model,true);
  else if( filename == "H2" )
    se3::buildModels::humanoid2d(model);
  else
    model = se3::urdf::buildModel(filename,JointModelFreeFlyer());
  std::cout << "nq = " << model.nq << std::endl;
  
  se3::Data data(model);
//  VectorXd q = VectorXd::Random(model.nq);
//  VectorXd qdot = VectorXd::Random(model.nv);
//  VectorXd qddot = VectorXd::Random(model.nv);
  
  std::vector<VectorXd> qs     (NBT);
  std::vector<VectorXd> qdots  (NBT);
  std::vector<VectorXd> qddots (NBT);
  for(size_t i=0;i<NBT;++i)
  {
    qs[i]     = Eigen::VectorXd::Random(model.nq);
    qs[i].segment<4>(3) /= qs[i].segment<4>(3).norm();
    qdots[i]  = Eigen::VectorXd::Random(model.nv);
    qddots[i] = Eigen::VectorXd::Random(model.nv);
  }
  
  timer.tic();
  SMOOTH(NBT)
  {
    crba(model,data,qs[_smooth]);
  }
  std::cout << "CRBA = \t\t"; timer.toc(std::cout,NBT);
  
//  timer.tic();
//  SMOOTH(NBT)
//  {
//    rnea(model,data,qs[_smooth],qdots[_smooth],qddots[_smooth]);
//  }
//  std::cout << "RNEA = \t\t"; timer.toc(std::cout,NBT);
  return -1;
}
