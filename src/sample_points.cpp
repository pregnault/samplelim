// [[Rcpp::depends(BH)]]

// VolEsti (volume computation and sampling library)

// Copyright (c) 20012-2020 Vissarion Fisikopoulos
// Copyright (c) 2018-2020 Apostolos Chalkis

// Contributed and/or modified by Apostolos Chalkis, as part of Google Summer of Code 2018 and 2019 program.

#include <Rcpp.h>
#include <RcppEigen.h>
#include <chrono>
#include <boost/random.hpp>
#include <boost/random/uniform_int.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/uniform_real_distribution.hpp>
#include "cartesian_geom/cartesian_kernel.h"
#include "random_walks/random_walks.hpp"
#include "convex_bodies/hpolytope.h"
#include "sampling/sampling.hpp"

template <typename Polytope, typename RNGType, typename PointList, typename NT, typename Point>
void sample_from_polytope(Polytope &P, RNGType &rng, PointList &randPoints, unsigned int const &walkL, unsigned int const &numpoints,
                          NT const &L, Point const &StartingPoint, unsigned int const &nburns,
                          bool const &set_L, bool const &billiard, bool const &mirror)
{
    
    if (billiard)
    {
        if (set_L)
        {
            BilliardWalk WalkType(L);
            uniform_sampling(randPoints, P, rng, WalkType, walkL, numpoints, StartingPoint, nburns);
        }
        else
        {
            uniform_sampling<BilliardWalk>(randPoints, P, rng, walkL, numpoints,
                                           StartingPoint, nburns);
        }
    }
    else if (mirror)
    {
        if (set_L)
        {
            Mirror WalkType(L);
            uniform_sampling(randPoints, P, rng, WalkType, walkL, numpoints, StartingPoint, nburns);
        }
        else
        {
            uniform_sampling<Mirror>(randPoints, P, rng, walkL, numpoints,
                                     StartingPoint, nburns);
        }
    }
    
}

// ' Sample uniformly distributed points from a convex Polytope (H-polytope, V-polytope, zonotope or intersection of two V-polytopes).
// '
// ' Sample n points with uniform distribution as the target distribution.
// '
// ' @param P A convex polytope. It is an object from class \link{Hpolytope}.
// ' @param n The number of points that the function is going to sample into the convex polytope P.
// ' @param random_walk Optional. A list that declares the random walk and some related parameters as follows:
// ' \itemize{
// ' \item{\code{walk} }{ A string to declare the random walk: i) \code{'CDHR'} for Coordinate Directions Hit-and-Run, ii) \code{'RDHR'} for Random Directions Hit-and-Run, iii) \code{'BaW'} for Ball Walk, iv) \code{'BiW'} for Billiard walk, v) \code{'BCDHR'} boundary sampling by keeping the extreme points of CDHR or vi) \code{'BRDHR'} boundary sampling by keeping the extreme points of RDHR. The default walk is \code{'BiW'} for the uniform distribution or \code{'CDHR'} for the Gaussian distribution.}
// ' \item{\code{walk_length} }{ The number of the steps per generated point for the random walk. The default value is 1.}
// ' \item{\code{nburns} }{ The number of points to burn before start sampling.}
// ' \item{\code{starting_point} }{ A \eqn{d}-dimensional numerical vector that declares a starting point in the interior of the polytope for the random walk. The default choice is the center of the ball as that one computed by the function \code{inner_ball()}.}
// ' \item{\code{BaW_rad} }{ The radius for the ball walk.}
// ' \item{\code{L} }{ The maximum length of the billiard trajectory.}
// ' \item{\code{seed} }{ A fixed seed for the number generator.}
// ' }
// ' @param distribution Optional. A list that declares the target density and some related parameters as follows:
// ' \itemize{
// ' \item{\code{density} }{ A string: (a) \code{'uniform'} for the uniform distribution or b) \code{'gaussian'} for the multidimensional spherical distribution. The default target distribution is uniform.}
// ' \item{\code{variance} }{ The variance of the multidimensional spherical gaussian. The default value is 1.}
// '  \item{\code{mode} }{ A \eqn{d}-dimensional numerical vector that declares the mode of the Gaussian distribution. The default choice is the center of the as that one computed by the function \code{inner_ball()}.}
// ' }
// '
// ' @return A \eqn{d\times n} matrix that contains, column-wise, the sampled points from the convex polytope P.
// ' @examples
// ' # uniform distribution from the 3d unit cube in H-representation using ball walk
// ' P = gen_cube(3, 'H')
// ' points = sample_points(P, n = 100, random_walk = list("walk" = "BaW", "walk_length" = 5))
// '
// ' # gaussian distribution from the 2d unit simplex in H-representation with variance = 2
// ' A = matrix(c(-1,0,0,-1,1,1), ncol=2, nrow=3, byrow=TRUE)
// ' b = c(0,0,1)
// ' P = Hpolytope(A = A, b = b)
// ' points = sample_points(P, n = 100, distribution = list("density" = "gaussian", "variance" = 2))
// '
// ' # uniform points from the boundary of a 2-dimensional random H-polytope
// ' P = gen_rand_hpoly(2,20)
// ' points = sample_points(P, n = 100, random_walk = list("walk" = "BRDHR"))
// '
// [[Rcpp::export]]
Rcpp::NumericMatrix sample_points(Rcpp::Reference P,
                                  int n,
                                  Rcpp::Nullable<Rcpp::List> random_walk = R_NilValue)
{

    typedef double NT;
    typedef Cartesian<NT> Kernel;
    typedef BoostRandomNumberGenerator<boost::mt19937, NT> RNGType;
    typedef typename Kernel::Point Point;
    typedef HPolytope<Point> Hpolytope;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, 1> VT;
    typedef Eigen::Matrix<NT, Eigen::Dynamic, Eigen::Dynamic> MT;
    unsigned int dim, walkL = 1;

 

    dim = Rcpp::as<MT>(P.slot("A")).cols();



    RNGType rng(dim);
    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("seed"))
    {
        unsigned seed2 = Rcpp::as<double>(Rcpp::as<Rcpp::List>(random_walk)["seed"]);
        rng.set_seed(seed2);
    }

    Hpolytope HP;

    unsigned int numpoints, nburns = 0;
    NT radius = 1.0, L;
    bool set_mode = false, billiard = false, set_starting_point = false, set_L = false, mirror = false;
    std::list<Point> randPoints;
    std::pair<Point, NT> InnerBall;
    Point mode(dim);

    numpoints = n;
    if (numpoints <= 0)
        throw Rcpp::exception("The number of samples has to be a positice integer!");


// Point de départ de la marche 

    Point StartingPoint;
    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("starting_point"))
    {
        if (Rcpp::as<VT>(Rcpp::as<Rcpp::List>(random_walk)["starting_point"]).size() != dim)
        {
            throw Rcpp::exception("Starting Point has to lie in the same dimension as the polytope P");
        }
        else
        {
            set_starting_point = true;
            VT temp = Rcpp::as<VT>(Rcpp::as<Rcpp::List>(random_walk)["starting_point"]);
            StartingPoint = Point(dim, std::vector<NT>(&temp[0], temp.data() + temp.cols() * temp.rows()));
        }
    }

    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("walk_length"))
    {
        walkL = Rcpp::as<int>(Rcpp::as<Rcpp::List>(random_walk)["walk_length"]);
        if (walkL <= 0)
        {
            throw Rcpp::exception("The walk length has to be a positive integer!");
        }
    }

    if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("nburns"))
    {
        nburns = Rcpp::as<int>(Rcpp::as<Rcpp::List>(random_walk)["nburns"]);
        if (nburns < 0)
        {
            throw Rcpp::exception("The number of points to burn before sampling has to be a positive integer!");
        }
    }

// Polytope definition


    // Hpolytope
        HP.init(dim, Rcpp::as<MT>(P.slot("A")),
                Rcpp::as<VT>(P.slot("b")));

        if (!set_starting_point)
        {
            InnerBall = HP.ComputeInnerBall();
            StartingPoint = InnerBall.first;

        }
        if (HP.is_in(StartingPoint) == 0)
        {
            throw Rcpp::exception("The given point is not in the interior of the polytope!");
        }
        HP.normalize();
   

// Random walk definition

    if (!random_walk.isNotNull() || !Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("walk") || Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("BiW")) == 0)
    {
        
        billiard = true;
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("L"))
        {
            L = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(random_walk)["L"]);
            set_L = true;
            if (L <= 0.0)
                throw Rcpp::exception("L must be a postitive number!");
        }
    
    }
    else if (Rcpp::as<std::string>(Rcpp::as<Rcpp::List>(random_walk)["walk"]).compare(std::string("mirror")) == 0)
    {
        billiard = false;
        mirror = true;
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("L"))
        {
            L = Rcpp::as<NT>(Rcpp::as<Rcpp::List>(random_walk)["L"]);
            set_L = true;
            if (L <= 0.0)
                throw Rcpp::exception("L must be a postitive number!");
        }

        std::vector<NT> jump;
        if (Rcpp::as<Rcpp::List>(random_walk).containsElementNamed("jump")) {
            jump = Rcpp::as<std::vector<NT>>(Rcpp::as<Rcpp::List>(random_walk)["jump"]);
            if (jump.size()==1){
                NT variance = jump.front();
                for(int i=1; i<dim; i++){
                    jump.push_back(variance);
                    }
            }
            else if (jump.size() != dim){
                throw Rcpp::exception("The jump must be either a unique value or a vector of size the dimension of the polytope!");
            }
        } 
        else {

            for(int i=0; i<dim; i++){
                jump.push_back(NT(1));
            }

        }
        rng.set_variances(dim, jump);
    }
    else
    {
        throw Rcpp::exception("Unknown walk type!");
    }

// Sampling function

    sample_from_polytope(HP, rng, randPoints, walkL, numpoints, L, StartingPoint, nburns,
                             set_L, billiard, mirror);

    MT RetMat(dim, numpoints);
    unsigned int jj = 0;

    for (typename std::list<Point>::iterator rpit = randPoints.begin(); rpit != randPoints.end(); rpit++, jj++)
    {
        RetMat.col(jj) = (*rpit).getCoefficients();
    }

    return Rcpp::wrap(RetMat);
}
