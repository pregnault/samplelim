
// Copyright (c) 20012-2018 Vissarion Fisikopoulos
// Copyright (c) 2018 Apostolos Chalkis


#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include "cartesian_geom/cartesian_kernel.h"
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "generators/boost_random_number_generator.hpp"
#include "convex_bodies/hpolytope.h"
//#include "volume/volume_sequence_of_balls.hpp"

// ' Compute an inscribed ball of a convex polytope
// '
// ' For a H-polytope described by a \eqn{m\times d} matrix \eqn{A} and a \eqn{m}-dimensional vector \eqn{b}, s.t.: \eqn{P=\{x\ |\  Ax\leq b\} }, this function computes the largest inscribed ball (Chebychev ball) by solving the corresponding linear program.
// ' For both zonotopes and V-polytopes the function computes the minimum \eqn{r} s.t.: \eqn{ r e_i \in P} for all \eqn{i=1, \dots ,d}. Then the ball centered at the origin with radius \eqn{r/ \sqrt{d}} is an inscribed ball.
// '
// ' @param P A convex polytope. It is an object from class (a) Hpolytope or (b) Vpolytope or (c) Zonotope or (d) VpolytopeIntersection.
// '
// ' @return A \eqn{(d+1)}-dimensional vector that describes the inscribed ball. The first \eqn{d} coordinates corresponds to the center of the ball and the last one to the radius.
// '
// ' @examples
// ' # compute the Chebychev ball of the 2d unit simplex
// ' P = gen_simplex(2,'H')
// ' ball_vec = inner_ball(P)
// '
// ' # compute an inscribed ball of the 3-dimensional unit cube in V-representation
// ' P = gen_cube(3, 'V')
// ' ball_vec = inner_ball(P)

 // [[Rcpp::export]]
Rcpp::NumericVector inner_ball(Rcpp::Reference P) {

    typedef double NT;
    typedef Cartesian<NT>    Kernel;
    typedef typename Kernel::Point    Point;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,1> VT;
    typedef Eigen::Matrix<NT,Eigen::Dynamic,Eigen::Dynamic> MT;

    unsigned int n, type_num;
    n = Rcpp::as<MT>(P.slot("A")).cols();

    

    std::pair <Point, NT> InnerBall;
    // Hpolytope
    Hpolytope HP;
    HP.init(n, Rcpp::as<MT>(P.slot("A")), Rcpp::as<VT>(P.slot("b")));
    InnerBall = HP.ComputeInnerBall();
            
    Rcpp::NumericVector vec(n + 1);
    for (unsigned int k = 0; k < n; ++k){
        vec[k] = InnerBall.first[k];
    }

    vec[n] = InnerBall.second;
    return vec;
}
