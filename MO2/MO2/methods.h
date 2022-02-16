#pragma once
#ifndef METHODS_H
#define METHODS_H

#include "utilities.h"

// Методы спуска
class descent_methods
{
	struct result
	{
		std::vector<double> x0;
		double eps = 0;
		uint32_t iter_cnt = 0;
		uint32_t call_to_func = 0;
		std::vector<double> x;
		double fx = 0;

		void save(std::string filename);
	};

public:
	descent_methods(const uint32_t dimension);

	result broyden(const VectorFunc& f, std::vector<double>& x0, double eps);

	result CGMFR(const VectorFunc& f, std::vector<double>& x0, double eps);

private:
	void print_iter(std::ofstream& res, uint32_t iter, std::vector<double>& x, const VectorFunc& f, double lambda);

	std::vector<double> x;
	std::vector<double> S;
	std::vector<double> grad_f;
	std::vector<std::vector<double>> etta_k;
};

//=====================================================================================
// Методы одномерного поиска
class one_dimensional_search_methods
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
	static result dichotomy(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result golden_ratio(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result fibonacci(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result parabola(const VectorFunc& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result find_interval(const VectorFunc& f, std::vector<double>& x, 
								std::vector<double>& S, double eps);

private:
	static double fibonacci_number(uint32_t n);
};

#endif