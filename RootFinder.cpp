// RootFinder.cpp : Defines the entry point for the console application.
//

#include <cmath>
#include "RootFinder.h"
#include <functional>
#include <vector>
#include <iostream>

template <typename Function>
class Func {
public:
	int neval = 0;
	Func(Function& f) :_f(f) {}
	double operator () (double x) {
		++this->neval;
		return this->_f(x);
	}
private:
	Function& _f;
};

/*! Run a root finding mission
*/
class Test {
public:
	/*! Test identifier: suite corresponds to one (parameterised) function
	*	test to a particular interval or parameter
	*/
	class Naming {
	public:
		int suite = 1;
		int test = 0;

		//! Go to next test id
		Naming& operator ++ () {
			++this->test;
			return *this;
		}

		//! Go to next suite
		Naming& next() {
			++this->suite;
			this->test = 0;
			return *this;
		}
	};

	Naming        name;
	Numeric::Eval root;

	Test(Naming const& naming, double a, double b, std::function<double(double)> f)
		: _f(f)
		, _a(a)
		, _b(b)
		, name(naming)
	{}

	int run(bool bisect=false) {
		using ftype = Func<std::function<double(double)>>;
		auto&& f = ftype(this->_f);
		if (bisect) {
			double x = std::numeric_limits<double>::quiet_NaN();
			if (! Numeric::FindRootBisect(this->_a, this->_b, f, x, 1E-15, 1E-15)) {
				throw Numeric::Error();
			}
			this->root = Numeric::Eval(x, this->_f);
		}
		else {
			auto&& solver = Numeric::RootFinder<ftype>(f);
			solver.root(this->_a, this->_b, 1E-15, 0);
			solver.assert_solution();
			this->root = solver.solution();
		}

		return f.neval;
	}
private:
	double _a, _b;
	std::function<double(double)> _f;
};

namespace Math {
	double const Pi = 2 * acos(0);
	double const e = exp(1);
}

inline double sqr(double x) { return x * x; }
inline double cub(double x) { return x * sqr(x); }
inline double pw4(double x) { return sqr(sqr(x)); }

static double f1(double x) {
	return sin(x) - x / 2;
}

static double f2(double x) {
	double f = 0;
	for (int i = 1; i <= 20; ++i) {
		f += sqr(2 * i - 5) / cub(x - i * i);
	}
	return -2 * f;
}

static double f3(double x, double a, double b) {
	return a * x * std::exp(b * x);
}

static double f4(double x, double a, double n) {
	return std::pow(x, n) - a;
}

static double f5(double x) {
	return sin(x) - 0.5;
}

static double f6(double x, double n) {
	return 2 * x * exp(-n) - 2 * exp(-n*x) + 1;
}

static double f7(double x, double n) {
	return (1 + sqr(1 - n))*x - sqr(1 - n*x);
}

static double f8(double x, double n) {
	return x * x - pow(1 - x, n);
}

static double f9(double x, double n) {
	return (1 + pw4(1 - n)) * x - pw4(1 - n * x);
}

static double f10(double x, double n) {
	return exp(-n*x)*(x - 1) + pow(x, n);
}

static double f11(double x, double n) {
	return (n*x - 1) / ((n - 1)*x);
}

static double f12(double x, double n) {
	return pow(x, 1 / n) - pow(n, 1 / n);
}

static double f13(double x) {
	if (x == 0) {
		return 0;
	}
	return x / exp(1 / sqr(x));
}

static double f14(double x, double n) {
	if (x >= 0) {
		return (n / 20)*((x / 1.5) + sin(x) - 1);
	}
	return -n / 20;
}

static double f15(double x, double n) {
	if (x < 0) {
		return -0.859;
	}
	if (x <= 2e-3 / (1 + n)) {
		return exp(1000 * x * (n + 1) / 2) - 1.859;
	}
	return Math::e - 1.859;
}

std::vector<Test> init() {
	std::vector<Test> tests;

	Test::Naming name;
	// 1
	tests.emplace_back(++name, Math::Pi / 2, Math::Pi, f1);

	// 2
	name.next();
	for (int n = 1; n <= 10; ++n) {
		tests.emplace_back(++name, sqr(n) + 1e-9, sqr(n + 1) - 1e-9, f2);
	}

	// 3
	name.next();
	tests.emplace_back(++name, -9, 31, [] (double x) { return f3(x, -40, -1); });
	tests.emplace_back(++name, -9, 31, [] (double x) { return f3(x, -100, -2); });
	tests.emplace_back(++name, -9, 31, [] (double x) { return f3(x, -200, -3); });

	// 4
	name.next();
	double a4[] = { 0.2, 1 };
	for (double a : a4) {
		for (int n = 4; n <= 12; n += 2) {
			tests.emplace_back(++name, 0, 5, [a, n] (double x) { return f4(x, a, n); });
		}
	}
	for (int n = 8; n <= 14; n += 12) {
		tests.emplace_back(++name, -0.95, 4.05, [n] (double x) { return f4(x, 1, n); });
	}

	// 5
	name.next();
	tests.emplace_back(++name, 0, 1.5, f5);

	// 6
	name.next();
	for (int n = 1; n <= 5; ++n) {
		tests.emplace_back(++name, 0, 1, [n] (double x) { return f6(x, n); });
	}
	for (int n = 20; n <= 100; n += 20) {
		tests.emplace_back(++name, 0, 1, [n] (double x) { return f6(x, n); });
	}

	// 7
	name.next();
	double n7[3] = { 5., 10, 20 };
	for (double n : n7) {
		tests.emplace_back(++name, 0, 1, [n] (double x) { return f7(x, n); });
	}

	// 8
	name.next();
	double n8[] = { 2, 5, 10, 15, 20 };
	for (double n : n8) {
		tests.emplace_back(++name, 0, 1, [n] (double x) { return f8(x, n); });
	}

	// 9
	name.next();
	double n9[] = { 1, 2, 4, 5, 8, 15, 20 };
	for (double n : n9) {
		tests.emplace_back(++name, 0, 1, [n] (double x) { return f9(x, n); });
	}

	// 10
	name.next();
	double n10[] = { 1, 5, 10, 15, 20 };
	for (double n : n10) {
		tests.emplace_back(++name, 0, 1, [n] (double x) { return f10(x, n); });
	}

	// 11
	name.next();
	double n11[] = { 2, 5, 15, 20 };
	for (double n : n11) {
		tests.emplace_back(++name, 0.01, 1, [n] (double x) { return f11(x, n); });
	}

	// 12
	name.next();
	for (int n = 2; n <= 6; ++n) {
		tests.emplace_back(++name, 1, 100, [n] (double x) { return f12(x, n); });
	}
	for (int n = 7; n <= 33; n += 2) {
		tests.emplace_back(++name, 1, 100, [n] (double x) { return f12(x, n); });
	}

	// 13
	name.next();
	tests.emplace_back(++name, -1, 4, f13);

	// 14
	name.next();
	for (int n = 1; n <= 40; ++n) {
		tests.emplace_back(++name, 1e-4, Math::Pi / 2, [n] (double x) { return f14(x, n); });
	}

	// 15
	name.next();
	for (int n = 20; n <= 40; ++n) {
		tests.emplace_back(++name, -1e4, 1e-4, [n] (double x) { return f15(x, n); });
	}
	for (int n = 100; n <= 1000; n += 100) {
		tests.emplace_back(++name, -1e4, 1e-4, [n] (double x) { return f15(x, n); });
	}
	return tests;
}

void run_tests(std::vector<Test> const& tests, bool bisect) {
	int neval = 0;
	for (auto&& t : init()) {
		int teval = t.run(bisect);
		auto& z = t.root;
		std::cout << "Test: " << t.name.suite << '.' << t.name.test << " NEval: " << teval << " Root: " << z.x << " Zero: " << z.fx << '\n';
		neval += teval;
	}
	std::cout << "Tests: " << tests.size() << " NEval: " << neval << std::endl;
}

int main(int argc, char* argv[]) {
	auto&& tests = init();
	run_tests(tests, true);
	run_tests(tests, false);
	return 0;
}

