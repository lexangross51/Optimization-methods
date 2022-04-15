#include "methods.hpp"

extern uint32_t calls_to_func;

std::mt19937_64 gen(1); // <- // ВОТ ТУТ МЕНЯТЬ ДЛЯ ИССЛЕДОВАНИЯ СО ЗВЕЗДОЧКОЙ
							 // ЖЕЛАТЕЛЬНО БРАТЬ СИЛЬНО ОТЛИЧАЮЩИЕСЯ ЗНАЧЕНИЯ
							 // МОЖНО ДАЖЕ БРАТЬ ОГРОМНЫЕ ЗНАЧЕНИЯ
std::uniform_real_distribution<> dist(-10, 10);

//=====================================================================================
// МЕТОДЫ СТАТИСТИЧЕСКОГО ПОИСКА
// Простой случайный поиск
statistical_search::result statistical_search::simple_random_search(
	const function2D& f,
	interval& x_int,
	interval& y_int,
	double eps, double P
)
{
	result res;
	res.eps = eps;
	res.P = P;

	// Определяем число испытаний
	// ---------------------------------------------------
	double V = (x_int.b - x_int.a) * (y_int.b - y_int.a);

	double V_eps = eps * eps;

	double P_eps = V_eps / V;

	uint32_t N_tests = log(1 - P) / log(1 - P_eps);

	res.N = N_tests;
	// ---------------------------------------------------

	std::vector<double> min_x = { dist(gen), dist(gen) };

	double min_fx = f(min_x);

	std::vector<double> x1, x2;

	double fx1, fx2;


	for (uint32_t test = 0; test < N_tests; test++)
	{
		x1 = { dist(gen), dist(gen) };
		x2 = { dist(gen), dist(gen) };

		fx1 = f(x1);
		fx2 = f(x2);

		if (fx1 < fx2 && fx1 < min_fx) {
			min_fx = fx1;
			min_x = x1;
		}
		else if (fx2 < fx1 && fx2 < min_fx) {
			min_fx = fx2;
			min_x = x2;
		}
	}

	res.x = min_x;
	res.fx = min_fx;

	return res;
}

// Метод глобального поиска №1
statistical_search::result statistical_search::global_search_1(
	const function2D& f,
	double eps,
	uint32_t m
)
{
	calls_to_func = 0;

	result res;
	res.eps = eps;

	std::vector<double> min_x = { dist(gen), dist(gen) };

	double min_fx = f(min_x);

	std::vector<double> x;

	for (uint32_t test = 0; test < m; )
	{
		x = { dist(gen), dist(gen) };

		auto direct_res = ds->rosenbrock(f, x, eps);
		res.call_to_func += direct_res.call_to_func;

		if (direct_res.fx < min_fx) 
		{
			min_x = direct_res.x;
			min_fx = direct_res.fx;
			test = 0;
		}
		else test++;
	}

	res.x = min_x;
	res.call_to_func += calls_to_func;
	res.fx = f(min_x);

	return res;
}

// Метод глобального поиска №2
statistical_search::result statistical_search::global_search_2(
	const function2D& f,
	double eps,
	uint32_t m
)
{
	calls_to_func = 0;

	result res;
	res.eps = eps;

	std::vector<double> min_x = { dist(gen), dist(gen) };
	double min_fx = f(min_x);

	std::vector<double> x;

	double fx;

	for (uint32_t test = 0; test < m; )
	{
		do {
			x = { dist(gen), dist(gen) };
			fx = f(x);
			test++;
		} while (fx > min_fx && test < m);

		if (test < m)
		{
			auto direct_res = ds->rosenbrock(f, x, eps);
			res.call_to_func += direct_res.call_to_func;

			if (direct_res.fx < min_fx)
			{
				min_x = direct_res.x;
				min_fx = direct_res.fx;
				test = 0;
			}
		}
	}

	res.x = min_x;
	res.call_to_func += calls_to_func;
	res.fx = f(min_x);

	return res;
}

// Метод глобального поиска №3
statistical_search::result statistical_search::global_search_3(
	const function2D& f,
	interval& x_int,
	interval& y_int,
	double eps,
	uint32_t m
)
{
	calls_to_func = 0;

	result res;
	res.eps = eps;

	std::vector<double> min_x = { dist(gen), dist(gen) };
	double min_fx = f(min_x);

	std::vector<double> x = min_x;

	double fx;

	for (uint32_t test = 0; test < m; )
	{
		auto direct_res = ds->rosenbrock(f, x, eps);
		res.call_to_func += direct_res.call_to_func;

		x = { dist(gen), dist(gen) };

		uint32_t step = 1;

		auto tmp = direct_res.x + step * (direct_res.x - x);

		double tmp_fx = f(tmp);

		while (
			tmp_fx >= direct_res.fx && 
			x_int.a < tmp[0] && tmp[0] < x_int.b && 
			y_int.a < tmp[1] && tmp[1] < y_int.b
		)
		{
			tmp = direct_res.x + step * (direct_res.x - x);
			tmp_fx = f(tmp);
			step++;
			test++;
		}

		if (direct_res.fx < min_fx)
		{
			min_x = direct_res.x;
			min_fx = direct_res.fx;
			test = 0;
		}
	}

	res.x = min_x;
	res.call_to_func += calls_to_func;
	res.fx = f(min_x);

	return res;
}


//=====================================================================================
// МЕТОДЫ ПРЯМОГО ПОИСКА 
direct_search::direct_search(const uint32_t dimension)
{
	x.resize(dimension);

	S1.resize(dimension);
	S2.resize(dimension);

	A1.resize(dimension);
	A2.resize(dimension);
}

direct_search::result direct_search::rosenbrock(
	const function2D& f, 
	std::vector<double>& x0, 
	double eps
)
{
	result res;

	res.eps = eps;
	res.x0 = x0;

	S1 = { 1.0, 0.0 };
	S2 = { 0.0, 1.0 };

	do
	{
		x = x0;

		// Минимизируем по каждому из направлений
		auto interval = one_dimensional_search::find_interval(f, x0, S1, eps).interval;
		double lambda1 = one_dimensional_search::golden_ratio(f, x0, interval, S1, eps).value;
		x0 = x0 + lambda1 * S1;

		interval = one_dimensional_search::find_interval(f, x0, S2, eps).interval;
		double lambda2 = one_dimensional_search::golden_ratio(f, x0, interval, S2, eps).value;
		x0 = x0 + lambda2 * S2;

		// Производим ортогонализацию
		A1 = lambda1 * S1 + lambda2 * S2;

		// Лямбды должны располагаться в порядке убывния по абсолютному значению
		if (lambda1 >= lambda2)
			A2 = lambda2 * S2;
		else
			A2 = lambda1 * S1;

		// Грамм - Жмых
		S1 = A1 / norm(A1);

		auto B2 = A2 - (A2 * S2) * S2;

		if (norm(B2) > 0.0)
			S2 = B2 / norm(B2);


		res.iter_cnt++;
		res.x = x0;

		//save("result_in_iter.txt", res);

	} while (abs(f(x) - f(x0)) > eps && abs(norm(x0) - norm(x)) > eps);

	res.call_to_func = calls_to_func;
	calls_to_func = 0;
	res.x = x0;
	res.fx = f(x0);

	return res;
}


//=====================================================================================
// МЕТОДЫ ОДНОМЕРНОГО ПОИСКА 
// Поиска интервала, содержащего минимум функции n переменных по направлению S
one_dimensional_search::result one_dimensional_search::find_interval(
	const function2D& f, 
	std::vector<double>& x, 
	std::vector<double>& S, 
	double eps
)
{
	result res;

	double delta = eps / 2.0;
	double h;

	double lambda = 0.0;
	double lambda_k1, lambda_k_1;

	double f_xk = f(x + lambda * S);
	double f_xk1 = f(x + (lambda + delta) * S);
	res.call_to_func += 2;

	if (f_xk > f_xk1)
	{
		lambda_k1 = lambda + delta;
		h = delta;
	}
	else
	{
		lambda_k1 = lambda - delta;
		h = -delta;
	}

	f_xk = f(x + lambda_k1 * S);
	do
	{
		lambda_k_1 = lambda;
		lambda = lambda_k1;
		f_xk = f_xk1;

		h *= 2;
		lambda_k1 = lambda + h;
		f_xk1 = f(x + lambda_k1 * S);

		res.call_to_func++;

	} while (f_xk > f_xk1);

	double a = lambda_k_1;
	double b = lambda_k1;

	if (b < a)
		std::swap(a, b);

	res.interval = interval(a, b);

	return res;
}

// Метод золотого сечения одномерного поиска
one_dimensional_search::result one_dimensional_search::golden_ratio(
	const function2D& f, 
	std::vector<double>& x, 
	interval& interval, 
	std::vector<double>& S, 
	double eps
)
{
	result res;

	double a = interval.a;
	double b = interval.b;

	double lambda1 = a + (3.0 - sqrt(5.0)) / 2.0 * (b - a);
	double lambda2 = a + (sqrt(5.0) - 1.0) / 2.0 * (b - a);

	double f_x1 = f(x + lambda1 * S);
	double f_x2 = f(x + lambda2 * S);

	for (; abs(b - a) >= eps; )
	{
		if (f_x1 > f_x2)
		{
			a = lambda1;
			lambda1 = lambda2;
			f_x1 = f_x2;
			lambda2 = a + (sqrt(5.0) - 1.0) / 2.0 * (b - a);
			f_x2 = f(x + lambda2 * S);
		}
		else
		{
			b = lambda2;
			lambda2 = lambda1;
			f_x2 = f_x1;
			lambda1 = a + (3.0 - sqrt(5.0)) / 2.0 * (b - a);
			f_x1 = f(x + lambda1 * S);
		}
	}

	res.value = (a + b) / 2.0;
	res.call_to_func = calls_to_func;
	return res;
}

// Метод парабол одномерного поиска
one_dimensional_search::result one_dimensional_search::parabola(
	const function2D& f, 
	std::vector<double>& x, 
	interval& interval, 
	std::vector<double>& S, 
	double eps
)
{
	result res;

	double a = interval.a;
	double b = interval.b;

	double lambda1 = a;
	double lambda2 = (a + b) / 2.0;
	double lambda3 = b;
	double lambda = 0, lambda_prev = 0;

	double f1 = f(x + lambda1 * S);
	double f2 = f(x + lambda2 * S);
	double f3 = f(x + lambda3 * S);

	while (true)
	{
		double c1 = f1;
		double c2 = (f2 - f1) / (lambda2 - lambda1);
		double c3 = ((f3 - f1) / (lambda3 - lambda1) - (f2 - f1) / (lambda2 - lambda1)) / (lambda3 - lambda2);

		lambda = (lambda1 + lambda2 - c2 / c3) / 2.0;

		double fx = f(x + lambda * S);

		if (abs(lambda - lambda_prev) < eps)
			break;

		lambda_prev = lambda;

		if (lambda > lambda2)
		{
			if (fx > f2)
			{
				lambda3 = lambda;
				f3 = fx;
			}
			else
			{
				lambda1 = lambda2;
				f1 = f2;
				lambda2 = lambda;
				f2 = fx;
			}
		}
		else
		{
			if (fx < f2)
			{
				lambda3 = lambda2;
				f3 = f2;
				lambda2 = lambda;
				f2 = fx;
			}
			else
			{
				lambda1 = lambda;
				f1 = fx;
			}
		}
		res.iters_cnt++;

	}

	res.call_to_func = calls_to_func;
	res.value = lambda;
	return res;
}