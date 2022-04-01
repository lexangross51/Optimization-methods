#pragma once
#ifndef METHODS_HPP
#define METHODS_HPP

#include "utilities.hpp"

enum class method
{
	barrier,
	fine
};

//=====================================================================================
// Методы прямого поиска
class direct_search
{
	struct result
	{
		double r0 = 0;
		double eps = 0;
		std::vector<double> x0;
		uint32_t iter_cnt = 0;
		uint32_t call_to_func = 0;
		std::vector<double> x;
		double fx = 0;

		void save(std::string filename);
	};

public:
	direct_search(const uint32_t dimension);

	void test1();
	void test2();
	void test3();
	void test4();

	void test5();
	void test6();
	void test7();
	void test8();

	result rosenbrock(const VectorFunc& f, const VectorFunc& Q, const VectorFunc& G, std::vector<double>& x0, double r0, double r_mult, double eps, method method);

	result gauss(const VectorFunc& f, const VectorFunc& Q, const VectorFunc& G, std::vector<double>& x0, double r0, double r_mult, double eps, method method);

private:
	void save(std::string filename, result& res);

	std::vector<double> x;

	std::vector<double> S1, S2;	// Система ортогональных направлений

	std::vector<double> A1, A2;	// Для Грамма - Жмыха
};


//=====================================================================================
// Методы одномерного поиска
class one_dimensional_search
{
	struct interval
	{
		double a;
		double b;

		interval(double _a = 0.0, double _b = 0.0) :
			a(_a), b(_b) {};
	};

	struct result
	{
		uint32_t iters_cnt = 0;
		uint32_t call_to_func = 0;
		double value = 0;
		interval interval;
	};

public:
	static result golden_ratio(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result parabola(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result fibonacci(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result find_interval(const VectorFunc& f, std::vector<double>& x, 
								std::vector<double>& S, double eps);

private:
	static double fibonacci_number(uint32_t n);
};
#endif