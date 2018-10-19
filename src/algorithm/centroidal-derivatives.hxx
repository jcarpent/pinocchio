//
// Copyright (c) 2018 INRIA
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

#ifndef __se3_algorithm_centroidal_derivatives_hxx__
#define __se3_algorithm_centroidal_derivatives_hxx__

#include "pinocchio/multibody/visitor.hpp"
#include "pinocchio/spatial/act-on-set.hpp"
#include "pinocchio/algorithm/kinematics.hpp"
#include "pinocchio/algorithm/check.hpp"

/// @cond DEV

namespace se3
{
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl,
           typename ConfigVectorType, typename TangentVectorType1, typename TangentVectorType2>
  struct CentroidalDynDerivativesForwardStep
  : public fusion::JointVisitorBase<CentroidalDynDerivativesForwardStep<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1,TangentVectorType2> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef typename Model::JointIndex JointIndex;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &,
                                  const ConfigVectorType &,
                                  const TangentVectorType1 &,
                                  const TangentVectorType2 &
                                  > ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     JointDataBase<typename JointModel::JointDataDerived> & jdata,
                     const Model & model,
                     Data & data,
                     const Eigen::MatrixBase<ConfigVectorType> & q,
                     const Eigen::MatrixBase<TangentVectorType1> & v,
                     const Eigen::MatrixBase<TangentVectorType2> & a)
    {
      typedef typename Model::JointIndex JointIndex;
      typedef typename Data::Motion Motion;
      
      const JointIndex & i = jmodel.id();
      const JointIndex & parent = model.parents[i];
      Motion & ov = data.ov[i];
      Motion & oa = data.oa[i];
      
      jmodel.calc(jdata.derived(),q.derived(),v.derived());
      
      data.liMi[i] = model.jointPlacements[i]*jdata.M();
      
      data.v[i] = jdata.v();
      
      if(parent > 0)
      {
        data.oMi[i] = data.oMi[parent] * data.liMi[i];
        data.v[i] += data.liMi[i].actInv(data.v[parent]);
      }
      else
        data.oMi[i] = data.liMi[i];
      
      data.a[i] = jdata.S() * jmodel.jointVelocitySelector(a) + jdata.c() + (data.v[i] ^ jdata.v());
      if(parent > 0)
      {
        data.a[i] += data.liMi[i].actInv(data.a[parent]);
      }
      
      data.oYcrb[i] = data.oMi[i].act(model.inertias[i]);
      ov = data.oMi[i].act(data.v[i]);
      oa = data.oMi[i].act(data.a[i]) + data.oa[0]; // add gravity contribution
      
      data.oh[i] = data.oYcrb[i] * ov;
      data.of[i] = data.oYcrb[i] * oa + ov.cross(data.oh[i]);
      
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;
      ColsBlock J_cols = jmodel.jointCols(data.J);
      ColsBlock dJ_cols = jmodel.jointCols(data.dJ);
      ColsBlock dVdq_cols = jmodel.jointCols(data.dVdq);
      ColsBlock dAdq_cols = jmodel.jointCols(data.dAdq);
      ColsBlock dAdv_cols = jmodel.jointCols(data.dAdv);
      
      J_cols = data.oMi[i].act(jdata.S());
      motionSet::motionAction(ov,J_cols,dJ_cols);
      motionSet::motionAction(data.oa[parent],J_cols,dAdq_cols);
      
      if(parent > 0)
      {
        motionSet::motionAction<ADDTO>(data.ov[parent],dJ_cols,dAdq_cols);
        dVdq_cols = dJ_cols;
        dAdv_cols.noalias() = (Scalar)2*dJ_cols;
      }
      else
      {
        dVdq_cols.setZero();
        dAdv_cols = dJ_cols;
      }
      
      // computes variation of inertias
      data.doYcrb[i] = data.oYcrb[i].variation(ov);
      
      addForceCrossMatrix(data.oh[i],data.doYcrb[i]);
    }
    
    template<typename ForceDerived, typename M6>
    static void addForceCrossMatrix(const ForceDense<ForceDerived> & f,
                                    const Eigen::MatrixBase<M6> & mout)
    {
      M6 & mout_ = EIGEN_CONST_CAST(M6,mout);
      addSkew(-f.linear(),mout_.template block<3,3>(ForceDerived::LINEAR,ForceDerived::ANGULAR));
      addSkew(-f.linear(),mout_.template block<3,3>(ForceDerived::ANGULAR,ForceDerived::LINEAR));
      addSkew(-f.angular(),mout_.template block<3,3>(ForceDerived::ANGULAR,ForceDerived::ANGULAR));
    }
    
  }; // struct CentroidalDynDerivativesForwardStep
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl>
  struct CentroidalDynDerivativesBackwardStep
  : public fusion::JointVisitorBase<CentroidalDynDerivativesBackwardStep<Scalar,Options,JointCollectionTpl> >
  {
    typedef ModelTpl<Scalar,Options,JointCollectionTpl> Model;
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    
    typedef boost::fusion::vector<const Model &,
                                  Data &
                                  > ArgsType;
    
    template<typename JointModel>
    static void algo(const JointModelBase<JointModel> & jmodel,
                     const Model & model,
                     Data & data)
    {
      typedef typename Model::JointIndex JointIndex;
      
      const JointIndex & i = jmodel.id();
      const JointIndex & parent = model.parents[i];
      
      typedef typename SizeDepType<JointModel::NV>::template ColsReturn<typename Data::Matrix6x>::Type ColsBlock;
      
      ColsBlock J_cols = jmodel.jointCols(data.J);
      ColsBlock dVdq_cols = jmodel.jointCols(data.dVdq);
      ColsBlock dAdq_cols = jmodel.jointCols(data.dAdq);
      ColsBlock dAdv_cols = jmodel.jointCols(data.dAdv);
      ColsBlock dFdq_cols = jmodel.jointCols(data.dFdq);
      ColsBlock dFdv_cols = jmodel.jointCols(data.dFdv);
      ColsBlock dFda_cols = jmodel.jointCols(data.dFda);
      
      // tau
      jmodel.jointVelocitySelector(data.tau).noalias() = J_cols.transpose()*data.of[i].toVector();
      
      // dtau/da similar to data.M
      motionSet::inertiaAction(data.oYcrb[i],J_cols,dFda_cols);

      // dtau/dv
      dFdv_cols.noalias() = data.doYcrb[i] * J_cols;
      motionSet::inertiaAction<ADDTO>(data.oYcrb[i],dAdv_cols,dFdv_cols);

      // dtau/dq
      if(parent>0)
      {
        dFdq_cols.noalias() = data.doYcrb[i] * dVdq_cols;
        motionSet::inertiaAction<ADDTO>(data.oYcrb[i],dAdq_cols,dFdq_cols);
      }
      else
        motionSet::inertiaAction(data.oYcrb[i],dAdq_cols,dFdq_cols);

      motionSet::act<ADDTO>(J_cols,data.of[i],dFdq_cols);
      
      data.oYcrb[parent] += data.oYcrb[i];
      data.doYcrb[parent] += data.doYcrb[i];
      data.oh[parent] += data.oh[i];
      data.of[parent] += data.of[i];
    }
    
    template<typename Min, typename Mout>
    static void lhsInertiaMult(const typename Data::Inertia & Y,
                               const Eigen::MatrixBase<Min> & J,
                               const Eigen::MatrixBase<Mout> & F)
    {
      Mout & F_ = EIGEN_CONST_CAST(Mout,F);
      motionSet::inertiaAction(Y,J.derived().transpose(),F_.transpose());
    }
  }; // struct CentroidalDynDerivativesBackwardStep
  
  namespace
  {
    // TODO: should me moved to ForceSet
    template<typename Matrix6xLikeIn, typename Vector3Like, typename Matrix6xLikeOut>
    inline void translateForceSet(const Eigen::MatrixBase<Matrix6xLikeIn> & Fin,
                                  const Eigen::MatrixBase<Vector3Like> & v3,
                                  const Eigen::MatrixBase<Matrix6xLikeOut> & Fout)
    {
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix6xLikeIn,6,Eigen::Dynamic)
      EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(Vector3Like,3)
      EIGEN_STATIC_ASSERT_MATRIX_SPECIFIC_SIZE(Matrix6xLikeOut,6,Eigen::Dynamic)

      assert(Fin.cols() == Fout.cols() && "Fin and Fout do not have the same number of columns");

      for(Eigen::Index k = 0; k < Fin.cols(); ++k)
      {
        typedef ForceRef<typename Matrix6xLikeIn::ColXpr> ForceTypeIn;
        typedef ForceRef<typename Matrix6xLikeOut::ColXpr> ForceTypeOut;
        ForceTypeOut fout(EIGEN_CONST_CAST(Matrix6xLikeOut,Fout).col(k));
        const ForceTypeIn fin(EIGEN_CONST_CAST(Matrix6xLikeIn,Fin).col(k));
        fout.linear() = fin.linear();
        fout.angular().noalias() = fin.angular() - v3.cross(fin.linear());
      }
    }
  } // internal namespace
  
  template<typename Scalar, int Options, template<typename,int> class JointCollectionTpl,
  typename ConfigVectorType, typename TangentVectorType1, typename TangentVectorType2,
  typename Matrix6xLike1, typename Matrix6xLike2, typename Matrix6xLike3>
  inline void
  computeCentroidalDynamicsDerivatives(const ModelTpl<Scalar,Options,JointCollectionTpl> & model,
                                       DataTpl<Scalar,Options,JointCollectionTpl> & data,
                                       const Eigen::MatrixBase<ConfigVectorType> & q,
                                       const Eigen::MatrixBase<TangentVectorType1> & v,
                                       const Eigen::MatrixBase<TangentVectorType2> & a,
                                       const Eigen::MatrixBase<Matrix6xLike1> & dhdot_dq,
                                       const Eigen::MatrixBase<Matrix6xLike2> & dhdot_dv,
                                       const Eigen::MatrixBase<Matrix6xLike3> & dhdot_da)
  {
    assert(q.size() == model.nq && "The joint configuration vector is not of right size");
    assert(v.size() == model.nv && "The joint velocity vector is not of right size");
    assert(a.size() == model.nv && "The joint acceleration vector is not of right size");
    assert(dhdot_dq.cols() == model.nv);
    assert(dhdot_dq.rows() == 6);
    assert(dhdot_dv.cols() == model.nv);
    assert(dhdot_dv.rows() == 6);
    assert(dhdot_da.cols() == model.nv);
    assert(dhdot_da.rows() == 6);
    assert(model.check(data) && "data is not consistent with model.");
    
    typedef CentroidalDynDerivativesForwardStep<Scalar,Options,JointCollectionTpl,ConfigVectorType,TangentVectorType1,TangentVectorType2> Pass1;
    for(JointIndex i=1; i<(JointIndex) model.njoints; ++i)
    {
      Pass1::run(model.joints[i],data.joints[i],
                 typename Pass1::ArgsType(model,data,q.derived(),v.derived(),a.derived()));
    }
    
    data.oYcrb[0].setZero();
    typedef CentroidalDynDerivativesBackwardStep<Scalar,Options,JointCollectionTpl> Pass2;
    for(JointIndex i=(JointIndex)(model.njoints-1); i>0; --i)
    {
      Pass2::run(model.joints[i],typename Pass2::ArgsType(model,data));
    }
    
    // expressed all the quantities around the center of mass
    typedef DataTpl<Scalar,Options,JointCollectionTpl> Data;
    typedef typename Data::Inertia Inertia;

    const Inertia & Ytot = data.oYcrb[0];
    const typename Inertia::Vector3 & com = Ytot.lever();
    
    // Center of mass of the system
    data.com[0] = Ytot.lever();
    
    // Compute the centroidal quantities
    data.hg = data.oh[0];
    data.hg.angular() += data.hg.linear().cross(data.com[0]);
    
    data.dhg = data.of[0];
    data.dhg.angular() += data.dhg.linear().cross(data.com[0]);
    
    // Compute centroidal inertia
    data.Ig.mass() = Ytot.mass();
    data.Ig.lever().setZero();
    data.Ig.inertia() = Ytot.inertia();
    
    // Compute the partial derivatives
    translateForceSet(data.dFdq,com,EIGEN_CONST_CAST(Matrix6xLike1,dhdot_dq));
    translateForceSet(data.dFdv,com,EIGEN_CONST_CAST(Matrix6xLike2,dhdot_dv));
    translateForceSet(data.dFda,com,EIGEN_CONST_CAST(Matrix6xLike3,dhdot_da));
  }
  
} // namespace se3

/// @endcond

#endif // ifndef __se3_algorithm_centroidal_derivatives_hxx__
