// VolEsti (volume computation and sampling library)

// Copyright (c) 2012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 program.

// Licensed under GNU LGPL.3, see LICENCE file

#ifndef RANDOM_WALKS_MIRROR_HPP
#define RANDOM_WALKS_MIRROR_HPP

#include <Eigen/Eigen>
#include <tuple>

#include "convex_bodies/ball.h"
#include "convex_bodies/ballintersectconvex.h"
#include "convex_bodies/hpolytope.h"
#include "convex_bodies/vpolytope.h"
#include "convex_bodies/vpolyintersectvpoly.h"
#include "convex_bodies/zpolytope.h"
#include "convex_bodies/zonoIntersecthpoly.h"
#include "generators/boost_random_number_generator.hpp"
#include "sampling/random_point_generators.hpp"
#include "sampling/sphere.hpp"
#include "volume/sampling_policies.hpp"


// Generate Length and direction according to rng


template <typename NT, typename Point>
struct GetLengthDirection_Mirror
{
    //typedef typename Point::FT NT;

    template <typename RandomNumberGenerator>
    inline static std::tuple<NT, Point> apply(unsigned int const& dim,
                              RandomNumberGenerator &rng)
    {
        NT length = NT(0);
        Point p(dim);
        std::vector<NT> multivariate = rng.sample_mndists();
        NT* data = p.pointerToData();

        for (unsigned int i=0; i<dim; ++i)
        {
            *data = multivariate[i];
            length += *data * *data;
            data++;
        }

        length = std::sqrt(length);
        p *= NT(1)/length;

        return {length,p};
    }
};



// Billiard walk for uniform distribution

struct Mirror
{
    Mirror(double L)
            :   param(L, true)
    {}

    Mirror()
            :   param(0, false)
    {}

    struct parameters
    {
        parameters(double L, bool set)
                :   m_L(L), set_L(set)
        {}
        double m_L;
        bool set_L;
    };

    parameters param;


template
<
    typename Polytope,
    typename RandomNumberGenerator
>
struct Walk
{
    typedef typename Polytope::PointType Point;
    typedef typename Point::FT NT;
    typedef HPolytope<Point> Hpolytope;
    typedef Zonotope<Point> zonotope;
    typedef ZonoIntersectHPoly <zonotope, Hpolytope> ZonoHPoly;
    typedef Ball<Point> BallType;
    typedef BallIntersectPolytope<Polytope,BallType> BallPolytope;

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p,  RandomNumberGenerator &rng)
    {
        _Len = compute_diameter<GenericPolytope>
                ::template compute<NT>(P);
        initialize(P, p, rng);
    }

    template <typename GenericPolytope>
    Walk(GenericPolytope const& P, Point const& p, RandomNumberGenerator &rng,
         parameters const& params)
    {
        _Len = params.set_L ? params.m_L
                          : compute_diameter<GenericPolytope>
                            ::template compute<NT>(P);
        initialize(P, p,  rng);
    }

    template
    <
        typename GenericPolytope
    >
    inline void apply(GenericPolytope const& P,
                      Point& p,   // a point to start
                      unsigned int const& walk_length,
                      RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        NT T;
        Point _v;
        const NT dl = 0.995;

        for (auto j=0u; j<walk_length; ++j)
        {
            /* T = rng.sample_urdist()*_Len;
            _v = GetDirection<Point>::apply(n, rng); */
            std::tie(T, _v)=GetLengthDirection_Mirror<NT, Point>::apply(n, rng);
            Point p0 = _p;
            int it = 0;
            while (it < 500*n)
            {
                auto pbpair = P.line_positive_intersect(_p, _v, _lambdas,
                                                        _Av, _lambda_prev);
                if (T <= pbpair.first) {
                    _p += (T * _v);
                    _lambda_prev = T;
                    break;
                }
                _lambda_prev = dl * pbpair.first;
                _p += (_lambda_prev * _v);
                T -= _lambda_prev;
                P.compute_reflection(_v, _p, pbpair.second);
                it++;
            }
            if (it == 50*n){
                _p = p0;
            }
        }
        p = _p;
    }

    inline void update_delta(NT L)
    {
        _Len = L;
    }

private :

    template
    <
        typename GenericPolytope
    >
    inline void initialize(GenericPolytope const& P,
                           Point const& p,
                           RandomNumberGenerator &rng)
    {
        unsigned int n = P.dimension();
        const NT dl = 0.995;
        _lambdas.setZero(P.num_of_hyperplanes());
        _Av.setZero(P.num_of_hyperplanes());
        _p = p;
        NT T;
        Point _v;


/*         _v = GetDirection<Point>::apply(n, rng);

        NT T = rng.sample_urdist()*_Len; */
        std::tie(T, _v)=GetLengthDirection_Mirror<NT, Point>::apply(n, rng);
        Point p0 = _p;
        int it = 0;

        std::pair<NT, int> pbpair
                = P.line_positive_intersect(_p, _v, _lambdas, _Av);
        if (T <= pbpair.first) {
            _p += (T * _v);
            _lambda_prev = T;
            return;
        }
        _lambda_prev = dl * pbpair.first;
        _p += (_lambda_prev * _v);
        T -= _lambda_prev;
        P.compute_reflection(_v, _p, pbpair.second);

        while (it <= 500*n) // change the limit from 50n to 500n
        {
            std::pair<NT, int> pbpair
                    = P.line_positive_intersect(_p, _v, _lambdas, _Av, _lambda_prev);
            if (T <= pbpair.first) {
                _p += (T * _v);
                _lambda_prev = T;
                break;
            }else if (it == 50*n) {
                _lambda_prev = rng.sample_urdist() * pbpair.first;
                _p += (_lambda_prev * _v);
                break;
            }
            _lambda_prev = dl * pbpair.first;
            _p += (_lambda_prev * _v);
            T -= _lambda_prev;
            P.compute_reflection(_v, _p, pbpair.second);
            it++;
        }
        //if (it == 30*n) _p = p0;

    }

    NT _Len;
    Point _p;
    Point _v;
    NT _lambda_prev;
    typename Point::Coeff _lambdas;
    typename Point::Coeff _Av;
};

};








#endif // RANDOM_WALKS_UNIFORM_BILLIARD_WALK_HPP
