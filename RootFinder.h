/*
*	Copyright 2017 Harry Dickson
*	
*	Licensed under the Apache License, Version 2.0 (the "License");
*	you may not use this file except in compliance with the License.
*	You may obtain a copy of the License at
*	
*	http ://www.apache.org/licenses/LICENSE-2.0
*	
*	Unless required by applicable law or agreed to in writing, software
*	distributed under the License is distributed on an "AS IS" BASIS,
*	WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
*	See the License for the specific language governing permissions and
*	limitations under the License.
*/

#include <algorithm>
#include <cassert>

/*! Root finding algorithms
*/

namespace Numeric {
	class Error : public std::exception {
	public:
		const char* what() const throw () {
			return "A Numeric Error occured";
		}
	};

	/*! Hold a function point for a real valued function of a single variable
	*/
	class Eval {
	public:
		double x = 0;
		double fx = 0;

		/*! Create an Eval calling the function
		*/
		template <typename Function>
		Eval(double x, Function&& f) : x(x), fx(f(x)) {}
		
		/*! Create an Eval given a function evaluation
		*	@param x the variable/ordinate
		*	@param fx the function value/co-ordinate
		*/
		Eval(double x, double fx) : x(x), fx(fx) {}

		Eval(Eval const& rhs) = default;
		Eval& operator = (Eval const& rhs) = default;
		Eval() = default;
	};

	//! Compute a tolerance adapted to the magnitude of the problem
	// @param x the value around which tolerance is required (x - tol, x + tol)
	// @param xtol additional absolute tolerance
	double AdjustedTolerance(double x, double xtol=0) {
		double abs_tolerance = 2 * std::numeric_limits<double>::epsilon();
		double rel_tolerance = fabs(x) * abs_tolerance;
		return std::max(abs_tolerance, rel_tolerance) + xtol;
	}


	/*! A function evaluation on a interval end points
	*/
	class Interval {
	public:
		/*! Construct an Interval from end points
		*	
		*/
		Interval(Eval const& a, Eval const& b) : _a(a), _b(b) { assert(a.x < b.x); }

		/*! Construct an interval from range and function
		*/
		template <typename Function>
		Interval(double a, double b, Function&& f) : _a(a, f), _b(b, f) { assert(this->_a.x < this->_b.x); }

		Interval(Interval const& rhs) = default;
		Interval() = default;
		Interval& operator = (Interval const&) = default;

		/*! Test if end points have opposite sign
		*/
		bool valid_zero() const {
			return this->_a.fx < 0 != this->_b.fx < 0;
		}
		
		//! Query for the lower bound
		Eval const& aval() const { return this->_a; }

		//! Query for the upper bound
		Eval const& bval() const { return this->_b; }

		//! The variable value at lower bound
		double a() const { return this->_a.x; }

		//! The variable value at upper bound
		double b() const { return this->_b.x; }

		//! The function value at lower bound f(a)
		double fa() const { return this->_a.fx; }

		//! The function value at upper bound f(b)
		double fb() const { return this->_b.fx; }

		//! Interval width
		double dx() const { return this->_b.x - this->_a.x; }

		//! Difference in function values at end points
		double df() const { return this->_b.fx - this->_a.fx; }

		//! Test if a value lies within interval
		// @param m the value to test
		// @param open if true test if m is in the open interval (a, b) else [a, b]
		bool contains(double m, bool open=true) const {
			if (open) {
				return m > this->_a.x && m < this->_b.x;
			}
			return m >= this->_a.x && m <= this->_b.x;
		}

		//! Compute the mid-point of the interval
		double interp_mid() const {
			return (this->_a.x + this->_b.x) / 2;
		}

		//! Compute the zero of the line from (a, fa) to (b, fb)
		double interp_linear() const {
			return this->_a.x - this->_a.fx * this->dx() / this->df();
		}

		//! Compute the point eqidistant from linear interp zero as end point
		//  of minimum function value
		// f(b) -|                       + 
		//       |                   '   |
		//       |               *       |
		//       |           '           |
		//    0 -|- - - -+ - - - - - - - | 
		//       |   '                   |
		// f(a) -+                       |
		//       |-------|-------|-------|
		//       a       i       x       b
		// 
		// z is L(z) = 0 for line, L, from (a, f(a)) to (b, f(b))
		// x is the inserted point s.t. x - z = z - a
		// x is guaranteed to lie withing interval.
		// @return the value x as described
		double interp_linear_mirrored() const {
			Eval const& u = fabs(this->fa()) < fabs(this->fb()) ? this->_a : this->_b;
			return u.x - 2 * u.fx * this->dx() / this->df();
		}

		//! Given an extra point (d, f(d)) form the quadratic fit to
		//	(a, f(a)), (b, f(b)) and (d, f(d)) and perform a number of
		//  Newton steps to find zero of quadratic
		//	@param k the number of steps to perform
		//	@param d te extra point
		double interp_quadratic(int k, Eval const& d) const {
			double a = this->a();
			double b = this->b();

			double a0 = this->fa();
			double a1 = this->df() / this->dx();
			double a2 = (d.fx - this->fb() / (d.x - b) - a1) / (d.x - a);
			
			// Check for degeneracy
			if (a2 == 0 || !std::isfinite(a2)) {
				return a - a0 / a1; // Linear interp
			}

			// Determine the starting point of newton steps.
			//
			double c = (a2 < 0 == a0 < 0) ? a : b;
			// Perform the steps
			for (int i = 0; i < k; ++i) {
				double pc  = a0 + (a1 + a2) * (c - b) * (c - a);
				double pdc = a1 + a2 * (2 * c) - (a + b);
				if (pdc == 0 || !std::isfinite(pdc)) {
					return a - a0 / a1; // Linear interp
				}
				c -= pc / pdc;
			}
			return c;
		}

		//! Given two extra points compute the zero of the cubic
		//  fitting the points (a, f(a)), (b, f(b)), (d, f(d)), (e, f(e))
		double interp_cubic(Eval const& d, Eval const& e) const {
			Eval const& a = this->_a;
			Eval const& b = this->_b;

			double q11 = (d.x - e.x)*d.fx / (e.fx - d.fx);
			double q21 = (b.x - d.x)*b.fx / (d.fx - b.fx);
			double q31 = (a.x - b.x)*a.fx / (b.fx - a.fx);
			double d21 = (b.x - d.x)*d.fx / (d.fx - b.fx);
			double d31 = (a.x - b.x)*b.fx / (b.fx - a.fx);
			double q22 = (d21 - q11)*b.fx / (e.fx - b.fx);
			double q32 = (d31 - q21)*a.fx / (d.fx - a.fx);
			double d32 = (d31 - q21)*d.fx / (d.fx - a.fx);
			double q33 = (d32 - q22)*a.fx / (e.fx - a.fx);

			return a.x + q31 + q32 + q33;
		}

		//! Force a value to lie with the reduced interval [a + tol, b - tol]
		//  If the tolerance reduced interval disappears the mid-point is taken
		//	@param m the value to constrain
		//	@param tol the amount to reduce the interval
		double constrain(double m, double tol) const {
			double imin = std::max(this->_a.x * (this->_a.x < 0 ? 1 - tol : 1 + tol), this->_a.x + tol);
			double imax = std::min(this->_b.x * (this->_b.x < 0 ? 1 + tol : 1 - tol), this->_b.x - tol);
			if (imin > imax) {
				return this->interp_mid();
			}
			if (m < imin) {
				return imin;
			}
			if (m > imax) {
				return imax;
			}
			return m;
		}

		//! Reduce the interval [a, b] to [m, b]
		//	@return the trimmed point, a
		Eval ltrim(Eval const& m) {
			Eval x = this->_a;
			this->_a = m;
			return x;
		}

		//! Reduce the interval [a, b] to [a, m]
		//	@return the trimmed point, b
		Eval rtrim(Eval const& m) {
			Eval x = this->_b;
			this->_b = m;
			return x;
		}

		//! Reduce the interval [a,b] to contain a 
		//	zero of the function
		Eval trim(Eval const& m) {
			return (m.fx < 0 == this->_a.fx < 0) ? this->ltrim(m) : this->rtrim(m);
		}
	private:
		Eval _a;
		Eval _b;
	};

	/*! Search for a simple zero of a function f:R=>R on an interval [a, b]
	*	Uses the method published as ACM 748 (more or less)
	*
	*   Algorithm 748: Enclosing Zeros of Continuous Functions, G. E. Alefeld, F. A. Potra and Yixun Shi, ACM Transactions on Mathematica1 Software, Vol. 21. No. 3. September 1995. Pages 327-344.
	*
	*	Some routines (interp_quadratic and interp_cubic are translations of the fortran code available from netlibs
	*
	*	This has also been implemented in boost as boost::math::tools::toms748_solve
	*/
	template <typename Function>
	class RootFinder {
	public:
		RootFinder(Function& f) : _function(f) {}

		/*! Query for the solution
		*	The caller should use solved() to check this is valid
		*/
		Eval const& solution() const {
			return this->_root;
		}

		/*! Compute the best interpolation
		*	@param k number of quadratic interpolations to perform if cubic fails
		*	@param d first extra point for cubic
		*	@param e second extra point for cubic
		*/
		double interp(int k, Eval const& d, Eval const& e) {

			// Check values are valid for interp_cubic
			//
			Eval const& a = this->_interval.aval();
			Eval const& b = this->_interval.bval();
			double test = (a.x - e.x) * (a.x - d.x) * (a.x - b.x) *
				          (b.x - e.x) * (b.x - d.x) *
						  (d.x - e.x);

			if (test == 0) {
				return this->_interval.interp_quadratic(k, d);
			}
			double m = this->_interval.interp_cubic(d, e);
			if (!this->_interval.contains(m)) {
				m = this->_interval.interp_quadratic(k, d);
			}
			return m;
		}

		/*! Record a success
		*	@param r the root found
		*/
		void set_solved(Eval const& r) {
			this->_root = r;
			this->_solved = true;
		}

		/*! Query if a root was found
		*/
		bool solved() const {
			return this->_solved;
		}

		/*! Compute a reduced interval and check for termination 
		*	@param m the candidate value
		*	@param c [out] the evaluation at m
		*	@param d [out] the redundant end point
		*	@param e [out] the original value of d
		*	@return true if a root was found
		*/
		bool bracket(double m, Eval& c, Eval& d, Eval& e) {
			double xtol = AdjustedTolerance(m, this->_xtol);
			
			m = this->_interval.constrain(m, xtol);
			c = Eval(m, this->_function(m));
			
			if (c.fx == 0) { //(c.fx) <= this->_ftol) {
				this->set_solved(c);
				return true;
			}
			
			e = d;
			d = this->_interval.trim(c);

			if (this->_interval.dx() < xtol) {
				// One last gasp attempt to get even closer
				c = Eval(this->_interval.interp_linear(), this->_function(m));
				this->set_solved(c);
				return true;
			}
			return false;
		}

		/*! Search for a root over an interval
		*	Root is found if function evaluates to (near) zero at any time
		*	or if interval is reduced to effectively zero width.
		*
		*	@param a lower bound of interval
		*	@param b upper bound of interval
		*	@param ftol halting condition on function if |f(x)| < ftol terminate successfully
		*	@param xtol halting condition on interval width if |b - a| < xtol terminate successfully
		*/
		bool root(double a, double b, double ftol, double xtol=0) {
			this->_interval = Interval(a, b, this->_function);
			this->_xtol = xtol;
			this->_ftol = ftol;
			this->_solved = false;

			if (!this->_interval.valid_zero()) {
				return false;
			}

			Eval c, d, e;

			// Factor to determine if an interval has
			// not been reduced satisfactorily. In which
			// case an extra reduction is applied
			//
			double const mu = 0.5;

			// We have 2 points (a, f(a)) and (b, f(b))
			// Do a linear interp (a bisection would also be valid)
			// This reduces the interval and yields
			// one exterior point, (d, f(d))
			//
			double m = this->_interval.interp_linear();
			if (this->bracket(m, c, d, e)) {
				return true;
			}

			// We have 3 points (a, f(a)), (b, f(b)) and (d, f(d))
			// Make a quadratic interpolation which reduces the
			// interval and yields an extra exterior point
			// (e, f(e)) becomes a copy of (d, f(d)) and (d, f(d))
			// is the new exterior point
			//
			m = this->_interval.interp_quadratic(2, d);
			if (this->bracket(m, c, d, e)) {
				return true;
			}

			// From now on we always have 4 function evaluations
			// (a, f(a)), (b, f(b)), (d, f(d)) and (e, f(e)) so
			// we can always attempt the cubic interpolation
			//
			for (int intnum = 2;; ++intnum) {
				// Remember interval width for later 
				double dx0 = this->_interval.dx();

				// Make the cubic interpolation guess
				// This yields a reduced interval and
				// a new exterior point (d, f(d)) with
				// (e, f(e)) set to prior (d, f(d))
				// 
				m = this->interp(3, d, e);
				if (this->bracket(m, c, d, e)) {
					return true;
				}

				// Now we introduce extra points to enhance 
				// the next cubic interpolation (presumably)

				// Introduce a new point equidistant from linear interpolated
				// zero as end point of least magnitude
				// Interval is now reduced to [a, c] with exterior points 
				// e = prior d, and d = redundant end point
				//
				m = this->_interval.interp_linear_mirrored();
				if (this->bracket(m, c, d, e)) {
					return true;
				}

				// If reduced interval is still large compared 
				// to original interval, introduce another point
				// which is simply the bisection of the interval
				// We reduce the interval again as usual
				//
				if (this->_interval.dx() > mu * dx0) {
					m = this->_interval.interp_mid();
					if (this->bracket(m, c, d, e)) {
						return true;
					}
				}
			}
		}

		/*! Check for a valid solution
		*/
		void assert_solution() {
			if (!this->_solved) {
				throw Error();
			}
			if (fabs(this->_root.fx) < this->_ftol) {
				return;
			}
			if (this->_interval.dx() <= AdjustedTolerance(this->_root.x, this->_xtol)) {
				return;
			}
			throw Error();
		}
	private:
		Interval  _interval;
		Function& _function;
		double    _xtol = 0;
		double    _ftol = 0;
		Eval      _root;
		bool      _solved = false;
	};

	/*! Find a simple zero of f:R=>R over an interval
	*	Uses the trival bisection algorithm.
	*	@param ax lower limit of interval
	*	@param bx upper limit of interval
	*	@param root [out] solution
	*	@param ftol terminate successfully if |f(x)| <= ftol. Can be zero to terminate on interval width
	*	@param xtol additional tolerance on interval width
	*/
	template <typename Function>
	bool FindRootBisect(double ax, double bx, Function&& fun, double& root, double ftol, double xtol=0) {
		Eval a(ax, fun);
		Eval b(bx, fun);

		Interval i(a, b);

		if (! i.valid_zero()) {
			return false;
		}

		if (fabs(a.fx) < ftol) {
			root = a.x;
			return true;
		}

		if (fabs(b.fx) < ftol) {
			root = b.x;
			return true;
		}

		double tol = AdjustedTolerance(i.interp_mid(), xtol);

		while (i.dx() > tol) {
			Eval m(i.interp_mid(), fun);

			// Test for termination
			if (fabs(m.fx) <= ftol) {
				root = m.x;
				return true;
			}

			// Compute new reduced interval
			i.trim(m);

			// Adjust tolerance
			tol = AdjustedTolerance(m.x, xtol);
		}

		root = i.interp_linear();
		return true;
	}
}