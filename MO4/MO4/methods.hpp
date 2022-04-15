#pragma once
#ifndef METHODS_HPP
#define METHODS_HPP

#include "utilities.hpp"

#pragma region Из старых лаб
//=====================================================================================
// Методы прямого поиска
class direct_search
{
	struct result
	{
		double eps = 0;
		std::vector<double> x0;
		uint32_t iter_cnt = 0;
		uint32_t call_to_func = 0;
		std::vector<double> x;
		double fx = 0;
	};

public:
	direct_search(const uint32_t dimension);

	result rosenbrock(const function2D& f, std::vector<double>& x0, double eps);

private:
	std::vector<double> x;

	std::vector<double> S1, S2;	// Система ортогональных направлений

	std::vector<double> A1, A2;	// Для Грамма - Жмыха
};


//=====================================================================================
// Методы одномерного поиска

struct interval
{
	double a;
	double b;

	interval(double _a = 0.0, double _b = 0.0) :
		a(_a), b(_b) {};
};

class one_dimensional_search
{
	struct result
	{
		uint32_t iters_cnt = 0;
		uint32_t call_to_func = 0;
		double value = 0;
		interval interval;
	};

public:
	static result golden_ratio(const function2D& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result parabola(const function2D& f, std::vector<double>& x,
		interval& interval, std::vector<double>& S, double eps);

	static result find_interval(const function2D& f, std::vector<double>& x,
		std::vector<double>& S, double eps);
};

#pragma endregion

//=====================================================================================
// Методы статистического поиска
class statistical_search
{
public:
	struct result {
		double eps = 0;
		double P = 0;
		uint32_t N = 0;
		std::vector<double> x;
		double fx = 0;
		uint32_t call_to_func = 0;
	};

public:
	statistical_search(const uint32_t dimenesion) {
		ds = new direct_search(dimenesion);
	}

	result simple_random_search(
		const function2D& f, 
		interval& x_int, 
		interval& y_int, 
		double eps, double P
	);

	result global_search_1(
		const function2D& f, 
		double eps, 
		uint32_t m
	);

	result global_search_2(
		const function2D& f,
		double eps,
		uint32_t m
	);

	result global_search_3(
		const function2D& f,
		interval& x_int,
		interval& y_int,
		double eps,
		uint32_t m
	);


private:
	direct_search* ds;
};

#endif