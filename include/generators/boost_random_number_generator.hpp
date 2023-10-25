// VolEsti (volume computation and sampling library)

// Copyright (c) 2020 Vissarion Fisikopoulos

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP
#define GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP

#include <chrono>
#include <boost/random.hpp>
#include <boost/function.hpp>
#include <vector>

/////////////////// Random numbers generator
///
/// \tparam RNGType
/// \tparam NT
/// \tparam Ts

template <typename RNGType, typename NT, int ... Ts>
struct BoostRandomNumberGenerator;

template <typename RNGType, typename NT>
struct BoostRandomNumberGenerator<RNGType, NT>
{
    BoostRandomNumberGenerator(int d)
            :   _rng(std::chrono::system_clock::now().time_since_epoch().count())
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
            ,   _mndists()
    {
      for(int i=0; i<d; i++)
        _mndists.push_back(boost::random::normal_distribution<NT>(0,1));
    }

  BoostRandomNumberGenerator(int d, std::vector<NT> variances)
            :   _rng(std::chrono::system_clock::now().time_since_epoch().count())
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
            ,   _mndists()
    {
      for(int i=0; i<d; i++)
        _mndists.push_back(boost::random::normal_distribution<NT>(0,variances[i]));
    }

    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

    std::vector<NT> sample_mndists()
    {
      std::vector<NT> sample;
      for(auto& dist : _mndists)
        sample.push_back(dist(_rng));
      return sample;
    }

    void set_variances(int d, std::vector<NT> variances)
    {
      for (int i=0; i<d; i++)
        _mndists[i]=boost::random::normal_distribution<NT>(0,variances[i]);
    }

  
    void set_seed(unsigned rng_seed)
    {
      _rng.seed(rng_seed);
    }

private :
  RNGType _rng;
  boost::random::uniform_real_distribution<NT> _urdist;
  boost::random::uniform_int_distribution<> _uidist;
  boost::random::normal_distribution<NT> _ndist;
  std::vector<boost::random::normal_distribution<NT>> _mndists;
};


template <typename RNGType, typename NT, int Seed>
struct BoostRandomNumberGenerator<RNGType, NT, Seed>
{
    BoostRandomNumberGenerator(int d)
            :   _rng(Seed)
            ,   _urdist(0, 1)
            ,   _uidist(0, d-1)
            ,   _ndist(0, 1)
            ,   _mndists()
  {
    for(int i=0; i<d; i++)
      _mndists.push_back(boost::random::normal_distribution<NT>(0,1));
  }

  BoostRandomNumberGenerator(int d, std::vector<NT> variances)
    :   _rng(Seed)
    ,   _urdist(0, 1)
    ,   _uidist(0, d-1)
    ,   _ndist(0, 1)
    ,   _mndists()
  {
    for(int i=0; i<d; i++)
      _mndists.push_back(boost::random::normal_distribution<NT>(0,variances[i]));
  }
  
    NT sample_urdist()
    {
        return _urdist(_rng);
    }

    NT sample_uidist()
    {
        return _uidist(_rng);
    }

    NT sample_ndist()
    {
        return _ndist(_rng);
    }

    std::vector<NT> sample_mndist()
    {
      std::vector<NT> sample;
      for(auto& dist : _mndists)
        sample.push_back(dist(_rng));
      return sample;
    }
    void set_variances(int d, std::vector<NT> variances)
    {
      for (int i=0; i<d; i++)
        _mndists[i]=boost::random::normal_distribution<NT>(0,variances[i]);
    }

    void set_seed(unsigned rng_seed){
        _rng.seed(rng_seed);
    }

private :
    RNGType _rng;
    boost::random::uniform_real_distribution<NT> _urdist;
    boost::random::uniform_int_distribution<> _uidist;
    boost::random::normal_distribution<NT> _ndist;
    std::vector<boost::random::normal_distribution<NT>> _mndists;
    
};

#endif // GENERATORS_BOOST_RANDOM_NUMBER_GENERATOR_HPP