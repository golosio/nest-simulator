/*
 *  stdp_connection.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef STDP_CONNECTION_H
#define STDP_CONNECTION_H

// C++ includes:
#include <cmath>

// Includes from nestkernel:
#include "common_synapse_properties.h"
#include "connection.h"
#include "connector_model.h"
#include "event.h"

// Includes from sli:
#include "dictdatum.h"
#include "dictutils.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Synapses
@ingroup stdp

Name: stdp_synapse - Synapse type for spike-timing dependent
plasticity.

Description:

stdp_synapse is a connector to create synapses with spike time
dependent plasticity (as defined in [1]). Here the weight dependence
exponent can be set separately for potentiation and depression.

Examples:

    multiplicative STDP [2]  mu_plus = mu_minus = 1.0
    additive STDP       [3]  mu_plus = mu_minus = 0.0
    Guetig STDP         [1]  mu_plus = mu_minus = [0.0,1.0]
    van Rossum STDP     [4]  mu_plus = 0.0 mu_minus = 1.0

Parameters:
\verbatim embed:rst
========= =======  ======================================================
 tau_plus  ms      Time constant of STDP window, potentiation
                   (tau_minus defined in post-synaptic neuron)
 lambda    real    Step size
 alpha     real    Asymmetry parameter (scales depressing increments as
                   alpha*lambda)
 mu_plus   real    Weight dependence exponent, potentiation
 mu_minus  real    Weight dependence exponent, depression
 Wmax      real    Maximum allowed weight
========= =======  ======================================================
\endverbatim

Transmits: SpikeEvent

References:

\verbatim embed:rst
.. [1] Guetig et al. (2003). Learning input correlations through nonlinear
       temporally asymmetric hebbian plasticity. Journal of Neuroscience,
       23:3697-3714 DOI: https://doi.org/10.1523/JNEUROSCI.23-09-03697.2003
.. [2] Rubin J, Lee D, Sompolinsky H (2001). Equilibrium
       properties of temporally asymmetric Hebbian plasticity. Physical Review
       Letters, 86:364-367. DOI: https://doi.org/10.1103/PhysRevLett.86.364
.. [3] Song S, Miller KD, Abbott LF (2000). Competitive Hebbian learning
       through spike-timing-dependent synaptic plasticity. Nature Neuroscience
       3(9):919-926.
       DOI: https://doi.org/10.1038/78829
.. [4] van Rossum MCW, Bi G-Q, Turrigiano GG (2000). Stable Hebbian learning
       from spike timing-dependent plasticity. Journal of Neuroscience,
       20(23):8812-8821.
       DOI: https://doi.org/10.1523/JNEUROSCI.20-23-08812.2000
\endverbatim

FirstVersion: March 2006

Author: Moritz Helias, Abigail Morrison

Adapted by: Philipp Weidel

SeeAlso: synapsedict, tsodyks_synapse, static_synapse
*/
// connections are templates of target identifier type (used for pointer /
// target index addressing) derived from generic connection template
template < typename targetidentifierT >
class STDPConnection : public Connection< targetidentifierT >
{

public:
  typedef CommonSynapseProperties CommonPropertiesType;
  typedef Connection< targetidentifierT > ConnectionBase;

  /**
   * Default Constructor.
   * Sets default values for all parameters. Needed by GenericConnectorModel.
   */
  STDPConnection();


  /**
   * Copy constructor.
   * Needs to be defined properly in order for GenericConnector to work.
   */
  STDPConnection( const STDPConnection& );

  // Explicitly declare all methods inherited from the dependent base
  // ConnectionBase. This avoids explicit name prefixes in all places these
  // functions are used. Since ConnectionBase depends on the template parameter,
  // they are not automatically found in the base class.
  using ConnectionBase::get_delay_steps;
  using ConnectionBase::get_delay;
  using ConnectionBase::get_rport;
  using ConnectionBase::get_target;

  /**
   * Get all properties of this connection and put them into a dictionary.
   */
  void get_status( DictionaryDatum& d ) const;

  /**
   * Set properties of this connection from the values given in dictionary.
   */
  void set_status( const DictionaryDatum& d, ConnectorModel& cm );

  /**
   * Send an event to the receiver of this connection.
   * \param e The event to send
   * \param cp common properties of all synapses (empty).
   */
  void send( Event& e, thread t, const CommonSynapseProperties& cp );

  double get_t_lastspike() { return t_lastspike_; }
  
  double get_den_delay() { return den_delay_; }

  void weight_update( double Dt, double tau_minus );

  void post_spike( double t_spike_post, double tau_minus );
  
  class ConnTestDummyNode : public ConnTestDummyNodeBase
  {
  public:
    // Ensure proper overriding of overloaded virtual functions.
    // Return values from functions are ignored.
    using ConnTestDummyNodeBase::handles_test_event;
    port
    handles_test_event( SpikeEvent&, rport )
    {
      return invalid_port_;
    }
  };

  void
  check_connection( Node& s, Node& t, rport receptor_type, const CommonPropertiesType& )
  {
    ConnTestDummyNode dummy_target;

    ConnectionBase::check_connection_( dummy_target, s, t, receptor_type );

    t.register_stdp_connection( t_lastspike_ - get_delay(),
				fabs( get_delay() - den_delay_ ) );
  }

  void
  set_previous_weight()
  {
    weight_ = previous_weight_;
  }

  void
  set_weight( double w )
  {
    weight_ = w;
  }

  double get_weight()
  {
    return weight_;
  }

private:
  double
  facilitate_( double w, double kplus )
  {
    double norm_w = ( w / Wmax_ ) + ( lambda_ * std::pow( 1.0 - ( w / Wmax_ ), mu_plus_ ) * kplus );
    return norm_w < 1.0 ? norm_w * Wmax_ : Wmax_;
  }

  double
  depress_( double w, double kminus )
  {
    double norm_w = ( w / Wmax_ ) - ( alpha_ * lambda_ * std::pow( w / Wmax_, mu_minus_ ) * kminus );
    return norm_w > 0.0 ? norm_w * Wmax_ : 0.0;
  }

  // data members of each connection
  double weight_;
  double tau_plus_;
  double lambda_;
  double alpha_;
  double mu_plus_;
  double mu_minus_;
  double Wmax_;
  //double Kplus_;
  double den_delay_;

  double t_lastspike_;
  double previous_weight_;
};


/**
 * Send an event to the receiver of this connection.
 * \param e The event to send
 * \param t The thread on which this connection is stored.
 * \param cp Common properties object, containing the stdp parameters.
 */
template < typename targetidentifierT >
inline void
STDPConnection< targetidentifierT >::send( Event& e, thread t, const CommonSynapseProperties& )
{
  // synapse STDP depressing/facilitation dynamics
  Node* target = get_target( t );

  double tau_minus = target->get_tau_minus();

  double delay_pre = get_delay();
  double delay_post = den_delay_;

  double t_spike_pre_last = e.get_stamp().get_ms();
  double t_spike_pre_lastm1 = t_lastspike_;

  double t_spike_post_last = 0;
  double t_spike_post_k = 0;
  double t_spike_post_km1 = 0;

  ///// CRITICAL SECTION /////
  std::mutex *target_spike_times_mtx = target->get_spike_times_mtx(); 
  target_spike_times_mtx->lock();
  std::vector<double> target_spike_times = target->get_spike_times();
  const uint n_post_spikes = target_spike_times.size();
  if (n_post_spikes>0) {
    t_spike_post_last = target_spike_times[n_post_spikes - 1];
  }
  uint k = n_post_spikes;
  if (n_post_spikes>0 && delay_post > delay_pre) {
    while (k>0 && (target_spike_times[k - 1] + delay_post)>
	   (t_spike_pre_last + delay_pre) ) {
      k--;
    }
    t_spike_post_k = target_spike_times[k];
    if (k>0) {
      t_spike_post_km1 = target_spike_times[k-1];
    }
  }
  target_spike_times_mtx->unlock();
  ///// END CRITICAL SECTION /////
  
  if (n_post_spikes>0) {
    if (delay_post > delay_pre) {
      // depression
      if (k>0) {
	if ( (t_spike_post_km1 + delay_post) >
	     (t_spike_pre_lastm1 + delay_pre) ) {
	  double Dt = t_spike_post_km1 + delay_post -
	    ( t_spike_pre_last + delay_pre );
	  weight_update(Dt, tau_minus);
	}
      }
      // facilitation
      if (k < n_post_spikes) {      
	if (t_spike_pre_lastm1>0) {
	  if ( k==0 || (t_spike_pre_lastm1 + delay_pre) >
	       (t_spike_post_km1 + delay_post) ) {
	    set_previous_weight();
	  }
	}
	double Dt = t_spike_post_k + delay_post -
	  (t_spike_pre_last + delay_pre);
	weight_update(Dt, tau_minus);
      }
    }
    else {
      if (t_spike_pre_lastm1==0.0 ||
	  (t_spike_post_last + delay_post >=
	   t_spike_pre_lastm1 + delay_pre ) ) {
	double Dt = t_spike_post_last + delay_post -
	  ( t_spike_pre_last + delay_pre );
	weight_update(Dt, tau_minus);
      }
    }
  }
  
  e.set_receiver( *target );
  e.set_weight( weight_ );
  // use accessor functions (inherited from Connection< >) to obtain delay in
  // steps and rport
  e.set_delay_steps( get_delay_steps() );
  e.set_rport( get_rport() );
  e();
  
  t_lastspike_ = t_spike_pre_last;
}


template < typename targetidentifierT >
STDPConnection< targetidentifierT >::STDPConnection()
  : ConnectionBase()
  , weight_( 1.0 )
  , tau_plus_( 20.0 )
  , lambda_( 0.01 )
  , alpha_( 1.0 )
  , mu_plus_( 1.0 )
  , mu_minus_( 1.0 )
  , Wmax_( 100.0 )
  //, Kplus_( 0.0 )
  , den_delay_( 1.0 )
  , t_lastspike_( 0.0 )
  , previous_weight_( 0.0 )
{
}

template < typename targetidentifierT >
STDPConnection< targetidentifierT >::STDPConnection( const STDPConnection< targetidentifierT >& rhs )
  : ConnectionBase( rhs )
  , weight_( rhs.weight_ )
  , tau_plus_( rhs.tau_plus_ )
  , lambda_( rhs.lambda_ )
  , alpha_( rhs.alpha_ )
  , mu_plus_( rhs.mu_plus_ )
  , mu_minus_( rhs.mu_minus_ )
  , Wmax_( rhs.Wmax_ )
  //, Kplus_( rhs.Kplus_ )
  , den_delay_( rhs.den_delay_ )
  , t_lastspike_( rhs.t_lastspike_ )
  , previous_weight_( rhs.previous_weight_ )
{
}

template < typename targetidentifierT >
void
STDPConnection< targetidentifierT >::get_status( DictionaryDatum& d ) const
{
  ConnectionBase::get_status( d );
  def< double >( d, names::weight, weight_ );
  def< double >( d, names::tau_plus, tau_plus_ );
  def< double >( d, names::lambda, lambda_ );
  def< double >( d, names::alpha, alpha_ );
  def< double >( d, names::mu_plus, mu_plus_ );
  def< double >( d, names::mu_minus, mu_minus_ );
  def< double >( d, names::Wmax, Wmax_ );
  def< double >( d, names::den_delay, den_delay_ );
  def< long >( d, names::size_of, sizeof( *this ) );
}

template < typename targetidentifierT >
void
STDPConnection< targetidentifierT >::set_status( const DictionaryDatum& d, ConnectorModel& cm )
{
  ConnectionBase::set_status( d, cm );
  updateValue< double >( d, names::weight, weight_ );
  updateValue< double >( d, names::tau_plus, tau_plus_ );
  updateValue< double >( d, names::lambda, lambda_ );
  updateValue< double >( d, names::alpha, alpha_ );
  updateValue< double >( d, names::mu_plus, mu_plus_ );
  updateValue< double >( d, names::mu_minus, mu_minus_ );
  updateValue< double >( d, names::Wmax, Wmax_ );
  updateValue< double >( d, names::den_delay, den_delay_ );

  // check if weight_ and Wmax_ has the same sign
  if ( not( ( ( weight_ >= 0 ) - ( weight_ < 0 ) ) == ( ( Wmax_ >= 0 ) - ( Wmax_ < 0 ) ) ) )
  {
    throw BadProperty( "Weight and Wmax must have same sign." );
  }
}

template < typename targetidentifierT >
void
STDPConnection< targetidentifierT >::weight_update(double Dt, double tau_minus)
{
  previous_weight_ = weight_;
  if (Dt>=0) {
    double fact = lambda_ * std::exp( -Dt / tau_plus_ );
    weight_ = weight_ + fact * Wmax_ * std::pow( 1.0 - weight_/Wmax_, mu_plus_);
  }
  else {
    double fact = -alpha_ * lambda_ * std::exp( Dt / tau_minus);
    weight_ = weight_ + fact * Wmax_ * std::pow( weight_/Wmax_, mu_minus_);
  }
  
  weight_ = weight_ >0.0 ? weight_ : 0.0;
  weight_ = weight_ < Wmax_ ? weight_ : Wmax_;
}

} // of namespace nest

#endif // of #ifndef STDP_CONNECTION_H
