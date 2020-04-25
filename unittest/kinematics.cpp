//
// Copyright (c) 2018-2019 CNRS INRIA
//

#include "pinocchio/multibody/model.hpp"
#include "pinocchio/multibody/data.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/crba.hpp"
#include "pinocchio/algorithm/joint-configuration.hpp"
#include "pinocchio/parsers/sample-models.hpp"

#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/utility/binary.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

BOOST_AUTO_TEST_CASE(test_kinematics_constant_vector_input)
{
  using namespace Eigen;
  using namespace pinocchio;
  
  Model model;
  buildModels::humanoidRandom(model);
  
  Data data(model);
  
  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  
  forwardKinematics(model,data,Model::ConfigVectorType::Ones(model.nq));
}

BOOST_AUTO_TEST_CASE(test_kinematics_zero)
{
  using namespace Eigen;
  using namespace pinocchio;
  
  Model model;
  buildModels::humanoidRandom(model);
  
  Data data(model), data_ref(model);
  
  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  
  forwardKinematics(model,data_ref,q);
  crba(model,data,q);
  updateGlobalPlacements(model,data);
  
  for(Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    BOOST_CHECK(data.oMi[i] == data_ref.oMi[i]);
    BOOST_CHECK(data.liMi[i] == data_ref.liMi[i]);
  }
}

BOOST_AUTO_TEST_CASE(test_kinematics_first)
{
  using namespace Eigen;
  using namespace pinocchio;
  
  Model model;
  buildModels::humanoidRandom(model);
  
  Data data(model);
  
  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v(VectorXd::Zero(model.nv));
  
  forwardKinematics(model,data,q,v);
  
  for(Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    BOOST_CHECK(data.v[i] == Motion::Zero());
  }
}

BOOST_AUTO_TEST_CASE(test_kinematics_second)
{
  using namespace Eigen;
  using namespace pinocchio;
  
  Model model;
  buildModels::humanoidRandom(model);
  
  Data data(model);
  
  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v(VectorXd::Zero(model.nv));
  VectorXd a(VectorXd::Zero(model.nv));
  
  forwardKinematics(model,data,q,v,a);
  
  for(Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    BOOST_CHECK(data.v[i] == Motion::Zero());
    BOOST_CHECK(data.a[i] == Motion::Zero());
  }
}

BOOST_AUTO_TEST_CASE(test_get_velocity)
{
  using namespace Eigen;
  using namespace pinocchio;

  Model model;
  buildModels::humanoidRandom(model);

  Data data(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v(VectorXd::Random(model.nv));

  forwardKinematics(model,data,q,v);

  for(Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    BOOST_CHECK(data.v[i].isApprox(getVelocity(model,data,i)));
    BOOST_CHECK(data.v[i].isApprox(getVelocity(model,data,i,ReferenceFrame::LOCAL)));
    BOOST_CHECK(data.oMi[i].act(data.v[i]).isApprox(getVelocity(model,data,i,ReferenceFrame::WORLD)));
    BOOST_CHECK(SE3(data.oMi[i].rotation(), Eigen::Vector3d::Zero()).act(data.v[i]).isApprox(getVelocity(model,data,i,ReferenceFrame::LOCAL_WORLD_ALIGNED)));
  }
}

BOOST_AUTO_TEST_CASE(test_get_acceleration)
{
  using namespace Eigen;
  using namespace pinocchio;

  Model model;
  buildModels::humanoidRandom(model);

  Data data(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v(VectorXd::Random(model.nv));
  VectorXd a(VectorXd::Random(model.nv));

  forwardKinematics(model,data,q,v,a);

  for(Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    BOOST_CHECK(data.a[i].isApprox(getAcceleration(model,data,i)));
    BOOST_CHECK(data.a[i].isApprox(getAcceleration(model,data,i,ReferenceFrame::LOCAL)));
    BOOST_CHECK(data.oMi[i].act(data.a[i]).isApprox(getAcceleration(model,data,i,ReferenceFrame::WORLD)));
    BOOST_CHECK(SE3(data.oMi[i].rotation(), Eigen::Vector3d::Zero()).act(data.a[i]).isApprox(getAcceleration(model,data,i,ReferenceFrame::LOCAL_WORLD_ALIGNED)));
  }
}

BOOST_AUTO_TEST_CASE(test_get_classical_acceleration)
{
  using namespace Eigen;
  using namespace pinocchio;

  Model model;
  buildModels::humanoidRandom(model);

  Data data(model);

  model.lowerPositionLimit.head<3>().fill(-1.);
  model.upperPositionLimit.head<3>().fill(1.);
  VectorXd q = randomConfiguration(model);
  VectorXd v(VectorXd::Random(model.nv));
  VectorXd a(VectorXd::Random(model.nv));

  forwardKinematics(model,data,q,v,a);

  for(Model::JointIndex i = 1; i < (Model::JointIndex)model.njoints; ++i)
  {
    SE3 T = data.oMi[i];
    Motion vel = data.v[i];
    Motion acc = data.a[i];
    Vector3d linear;

    Motion acc_classical_local = acc;
    linear = acc.linear() + vel.angular().cross(vel.linear());
    acc_classical_local.linear() = linear;

    BOOST_CHECK(acc_classical_local.isApprox(getClassicalAcceleration(model,data,i)));
    BOOST_CHECK(acc_classical_local.isApprox(getClassicalAcceleration(model,data,i,ReferenceFrame::LOCAL)));

    Motion vel_world = T.act(vel);
    Motion acc_classical_world = T.act(acc);
    linear = acc_classical_world.linear() + vel_world.angular().cross(vel_world.linear());
    acc_classical_world.linear() = linear;

    BOOST_CHECK(acc_classical_world.isApprox(getClassicalAcceleration(model,data,i,ReferenceFrame::WORLD)));

    Motion vel_aligned = SE3(data.oMi[i].rotation(), Eigen::Vector3d::Zero()).act(vel);
    Motion acc_classical_aligned = SE3(data.oMi[i].rotation(), Eigen::Vector3d::Zero()).act(acc);
    linear = acc_classical_aligned.linear() + vel_aligned.angular().cross(vel_aligned.linear());
    acc_classical_aligned.linear() = linear;

    BOOST_CHECK(acc_classical_aligned.isApprox(getClassicalAcceleration(model,data,i,ReferenceFrame::LOCAL_WORLD_ALIGNED)));
  }
}

BOOST_AUTO_TEST_SUITE_END()
