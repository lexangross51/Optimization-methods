#include "methods.hpp"

extern uint32_t calls_to_func;
double r;

//------------------------------------------------------------------------------------------------------------------
// МЕТОД ШТРАФНЫХ И БАРЬЕРНЫХ ФУНКЦИЙ
// Исходная функция
VectorFunc f = [](const std::vector<double>& x)
{
	calls_to_func++;
	return 4 * pow((x[1] - x[0]), 2) + 3 * pow((x[0] - 1), 2);
};


// Пункт а
VectorFunc barrier_f = [](const std::vector<double>& x) { return x[0] + x[1] + 1; };

// Пункт б
VectorFunc fine_f = [](const std::vector<double>& x) { return fabs(x[1] - x[0] - 1); };

//------------------------------------------------------------------------------------------------------------------
// Барьерные функции
VectorFunc G1 = [](const std::vector<double>& x)
{
	double val = 0.0;

	if ((val = barrier_f(x)) <= 0)
		return -1.0 / val;
	else
		return std::numeric_limits<double>::infinity();
};

VectorFunc G2 = [](const std::vector<double>& x)
{
	double val = 0.0;

	if ((val = barrier_f(x)) <= 0)
		return -log(-val);
	else
		return std::numeric_limits<double>::infinity();
};

// Штрафные функции
VectorFunc G3 = [](const std::vector<double>& x)
{
	double val = 0.0;

	if ((val = fine_f(x)) > 0)
		return (val + fabs(val)) / 2.0;
	else
		return 0.0;
};

VectorFunc G4 = [](const std::vector<double>& x)
{
	double val = 0.0;

	if ((val = fine_f(x)) > 0)
		return pow( (val + fabs(val)) / 2.0, 2);
	else
		return 0.0;
};

VectorFunc G5 = [](const std::vector<double>& x)
{
	double val = 0.0;

	if ((val = fine_f(x)) > 0)
		return pow( (val + fabs(val)) / 2.0, 4);
	else
		return 0.0;
};
//------------------------------------------------------------------------------------------------------------------

// Минимизируемые функции
VectorFunc Q1 = [](const std::vector<double>& x) { return f(x) + r * G1(x); };
VectorFunc Q2 = [](const std::vector<double>& x) { return f(x) + r * G2(x); };
VectorFunc Q3 = [](const std::vector<double>& x) { return f(x) + r * G3(x); };
VectorFunc Q4 = [](const std::vector<double>& x) { return f(x) + r * G4(x); };
VectorFunc Q5 = [](const std::vector<double>& x) { return f(x) + r * G5(x); };

//------------------------------------------------------------------------------------------------------------------
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
	const VectorFunc& f, const VectorFunc& Q, const VectorFunc& G,
	std::vector<double>& _x0, double r0, double r_mult, double eps, method method
)
{
	r = r0;

	result res;

	calls_to_func = 0;

	auto x0 = _x0;

	res.eps = eps;
	res.x0 = x0;

	S1 = { 1.0, 0.0 };
	S2 = { 0.0, 1.0 };

	do
	{
		x = x0;

		// Минимизируем по каждому из направлений
		auto interval = one_dimensional_search::find_interval(Q, x0, S1, eps).interval;
		double lambda1 = one_dimensional_search::fibonacci(Q, x0, interval, S1, eps).value;
		x0 = x0 + lambda1 * S1;

		interval = one_dimensional_search::find_interval(Q, x0, S2, eps).interval;
		double lambda2 = one_dimensional_search::fibonacci(Q, x0, interval, S2, eps).value;
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
		res.r0 = r;

		save("result_in_iter.txt", res);


		if (method == method::fine)
			r *= r_mult;
		else
			r /= r_mult;

	} while (abs(Q(x) - Q(x0)) > eps && abs(norm(x0) - norm(x)) > eps && r * G(x0) > eps);

	res.call_to_func = calls_to_func;
	calls_to_func = 0;
	res.x = x0;
	res.r0 = r0;
	res.fx = f(x0);

	return res;
}

direct_search::result direct_search::gauss(
	const VectorFunc& f, const VectorFunc& Q, const VectorFunc& G,
	std::vector<double>& _x0, double r0, double r_mult, double eps, method method
)
{
	calls_to_func = 0;

	result res;

	auto x0 = _x0;

	res.eps = eps;
	res.x0 = x0;

	// Задаем ортонормированный базис
	S1 = { 1.0, 0.0 };
	S2 = { 0.0, 1.0 };

	do
	{
		x = x0;

		auto interval = one_dimensional_search::find_interval(Q, x0, S1, eps).interval;
		double lambda1 = one_dimensional_search::fibonacci(Q, x0, interval, S1, eps).value;
		x0 = x0 + lambda1 * S1;

		interval = one_dimensional_search::find_interval(Q, x0, S2, eps).interval;
		double lambda2 = one_dimensional_search::fibonacci(Q, x0, interval, S1, eps).value;
		x0 = x0 + lambda2 * S2;


		res.iter_cnt++;
		res.x = x0;
		res.r0 = r;

		save("result_in_iter.txt", res);


		if (method == method::fine)
			r *= r_mult;
		else
			r /= r_mult;

	} while (abs(Q(x) - Q(x0)) > eps && abs(norm(x) - norm(x0)) > eps && r * G(x0) > eps);

	res.call_to_func = calls_to_func;
	calls_to_func = 0;
	res.x = x0;
	res.r0 = r0;
	res.fx = f(x0);

	return res;
}

void direct_search::result::save(std::string filename)
{
	FILE* res;

	if (!fopen_s(&res, filename.c_str(), "a"))
	{
		fprintf(res, "%-12s|%-24s|%-12s|%-12s|%-13s|%-24s|%-12s\n",
			"r0", "x0", "eps", "iters cnt", "call to func", "x", "f(x)");
		fprintf(res, "------------------------------------");
		fprintf(res, "------------------------------------");
		fprintf(res, "------------------------------------\n");

		fprintf(res, "%-12lf|(%-10lf, %-10lf)|%-12e|%-12ld|%-13ld|(%-10lf, %-10lf)|%-12e\n\n\n\n",
			r0, x0[0], x0[1], eps, iter_cnt, call_to_func, x[0], x[1], fx);

		fclose(res);
	}
}

void direct_search::save(std::string filename, result& rslt)
{
	FILE* res;

	if (!fopen_s(&res, filename.c_str(), "a"))
	{
		if (rslt.iter_cnt == 1)
		{
			fprintf(res, "%-12s|%-24s|%-12s|%-12s|%-13s|%-24s|%-12s\n",
				"r0", "x0", "eps", "iters cnt", "call to func", "x", "f(x)");
			fprintf(res, "------------------------------------");
			fprintf(res, "------------------------------------");
			fprintf(res, "------------------------------------\n");
		}

		fprintf(res, "%-12lf|(%-10lf, %-10lf)|%-12e|%-12ld|%-13ld|(%-10lf, %-10lf)|%-12e\n\n",
			rslt.r0, rslt.x0[0], rslt.x0[1], rslt.eps, rslt.iter_cnt, rslt.call_to_func, rslt.x[0], rslt.x[1], rslt.fx);

		fclose(res);
	}
}


// ==========================================================================================================
// ==================================== ТЕСТИРОВАНИЕ И ИССЛЕДОВАНИЯ =========================================
// ==========================================================================================================

// Здесь сразу исследования дял разной точности и разных штрафных функций 
void direct_search::test1()
{
	double epsilons[] = { 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 };

	std::vector<double> x0 = { -2, -1 };
	double r0 = 10;
	double r_mult = 2;

	for (const auto eps : epsilons)
	{
		auto res = rosenbrock(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
		res.save("result1.txt");

		res = rosenbrock(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
		res.save("result2.txt");

		res = rosenbrock(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
		res.save("result3.txt");

		res = rosenbrock(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
		res.save("result4.txt");

		res = rosenbrock(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
		res.save("result5.txt");
	}
}

// Исследование со стратегией изменения коэффициента штрафа
void direct_search::test2()
{
	std::vector<double> x0 = { -2, -1 };
	double r0 = 10;
	double eps = 1e-12;	// Можно поиграться с этим значением

	for (double r_mult = 2; r_mult <= 5; r_mult++)
	{
		auto res = rosenbrock(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
		res.save("result1.txt");

		res = rosenbrock(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
		res.save("result2.txt");

		res = rosenbrock(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
		res.save("result3.txt");

		res = rosenbrock(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
		res.save("result4.txt");

		res = rosenbrock(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
		res.save("result5.txt");
	}
}

// Исследования с изменением начального значения коэффициента штрафа
void direct_search::test3()
{
	std::vector<double> x0 = { -2, -1 };
	double r_mult = 2;
	double eps = 1e-12;	// Можно поиграться с этим значением

	for (double r0 = 2; r0 <= 32; r0 *= 2)
	{
		auto res = rosenbrock(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
		res.save("result1.txt");

		res = rosenbrock(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
		res.save("result2.txt");

		res = rosenbrock(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
		res.save("result3.txt");

		res = rosenbrock(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
		res.save("result4.txt");

		res = rosenbrock(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
		res.save("result5.txt");
	}
}

// Исследования с разными начальными точками
void direct_search::test4()
{
	std::vector<double> x0 = { 0, -2 };
	double r0 = 10;
	double r_mult = 2;
	double eps = 1e-12;	// Можно поиграться с этим значением

	x0 = { 0, -2 };
	auto res = rosenbrock(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
	res.save("result1.txt");

	res = rosenbrock(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
	res.save("result2.txt");

	res = rosenbrock(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
	res.save("result3.txt");

	res = rosenbrock(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
	res.save("result4.txt");

	res = rosenbrock(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
	res.save("result5.txt");

	// -------------------------------------------------------------------------------

	x0 = { 2, 3 };
	res = rosenbrock(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
	res.save("result1.txt");

	res = rosenbrock(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
	res.save("result2.txt");

	res = rosenbrock(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
	res.save("result3.txt");

	res = rosenbrock(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
	res.save("result4.txt");

	res = rosenbrock(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
	res.save("result5.txt");
}


// Все те же самые тесты, только для метода Гаусса
// Здесь сразу исследования дял разной точности и разных штрафных функций 
void direct_search::test5()
{
	double epsilons[] = { 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 };

	std::vector<double> x0 = { -2, -1 };
	double r0 = 10;
	double r_mult = 2;

	for (const auto eps : epsilons)
	{
		auto res = gauss(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
		res.save("result1.txt");

		res = gauss(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
		res.save("result2.txt");

		res = gauss(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
		res.save("result3.txt");

		res = gauss(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
		res.save("result4.txt");

		res = gauss(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
		res.save("result5.txt");
	}
}

// Исследование со стратегией изменения коэффициента штрафа
void direct_search::test6()
{
	std::vector<double> x0 = { -2, -1 };
	double r0 = 10;
	double eps = 1e-12;	// Можно поиграться с этим значением

	for (double r_mult = 2; r_mult <= 5; r_mult++)
	{
		auto res = gauss(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
		res.save("result1.txt");

		res = gauss(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
		res.save("result2.txt");

		res = gauss(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
		res.save("result3.txt");

		res = gauss(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
		res.save("result4.txt");

		res = gauss(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
		res.save("result5.txt");
	}
}

// Исследования с изменением начального значения коэффициента штрафа
void direct_search::test7()
{
	std::vector<double> x0 = { -2, -1 };
	double r_mult = 2;
	double eps = 1e-12;	// Можно поиграться с этим значением

	for (double r0 = 2; r0 <= 32; r0 *= 2)
	{
		auto res = gauss(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
		res.save("result1.txt");

		res = gauss(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
		res.save("result2.txt");

		res = gauss(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
		res.save("result3.txt");
;
		res = gauss(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
		res.save("result4.txt");

		res = gauss(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
		res.save("result5.txt");
	}
}

// Исследования с разными начальными точками
void direct_search::test8()
{
	std::vector<double> x0 = { 0, -2 };
	double r0 = 10;
	double r_mult = 2;
	double eps = 1e-12;	// Можно поиграться с этим значением

	x0 = { 0, -2 };
	auto res = gauss(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
	res.save("result1.txt");

	res = gauss(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
	res.save("result2.txt");

	res = gauss(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
	res.save("result3.txt");

	res = gauss(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
	res.save("result4.txt");

	res = gauss(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
	res.save("result5.txt");

	// -------------------------------------------------------------------------------

	x0 = { 2, 3 };
	res = gauss(f, Q1, G1, x0, r0, r_mult, eps, method::barrier);
	res.save("result1.txt");

	res = gauss(f, Q2, G2, x0, r0, r_mult, eps, method::barrier);
	res.save("result2.txt");

	res = gauss(f, Q3, G3, x0, r0, r_mult, eps, method::fine);
	res.save("result3.txt");

	res = gauss(f, Q4, G4, x0, r0, r_mult, eps, method::fine);
	res.save("result4.txt");

	res = gauss(f, Q5, G5, x0, r0, r_mult, eps, method::fine);
	res.save("result5.txt");
}



//------------------------------------------------------------------------------------------------------------------
// МЕТОДЫ ОДНОМЕРНОГО ПОИСКА 
// Поиска интервала, содержащего минимум функции n переменных по направлению S
one_dimensional_search::result one_dimensional_search::find_interval(
	const VectorFunc& f, std::vector<double>& x, std::vector<double>& S, double eps)
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
	const VectorFunc& f, std::vector<double>& x, interval& interval, std::vector<double>& S, double eps)
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
	const VectorFunc& f, std::vector<double>& x, interval& interval, std::vector<double>& S, double eps)
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

double one_dimensional_search::fibonacci_number(uint32_t n)
{
	double F1 = 1.0;
	double F2 = 1.0;

	for (uint32_t i = 2; i < n; i++)
	{
		double tmp = F2;
		F2 += F1;
		F1 = tmp;
	}

	return F2;
}

// Метод Фибоначчи одномерного поиска
one_dimensional_search::result one_dimensional_search::fibonacci(
	const VectorFunc& f, std::vector<double>& x, interval& interval, std::vector<double>& S, double eps)
{
	result res;

	uint32_t n = 0;

	double a = interval.a;
	double b = interval.b;
	double length = b - a;

	while (fibonacci_number(n) < length / eps) n++;

	double lambda1 = a + (fibonacci_number(n - 2) / fibonacci_number(n)) * (b - a);
	double f1 = f(x + lambda1 * S);

	double lambda2 = a + (fibonacci_number(n - 1) / fibonacci_number(n)) * (b - a);
	double f2 = f(x + lambda2 * S);

	for (uint32_t k = 0; k < n - 3; k++)
	{
		if (f1 <= f2)
		{
			b = lambda2;
			lambda2 = lambda1;
			f2 = f1;
			lambda1 = a + (fibonacci_number(n - k - 3) / fibonacci_number(n - k - 1)) * (b - a);
			f1 = f(x + lambda1 * S);
		}
		else
		{
			a = lambda1;
			lambda1 = lambda2;
			f1 = f2;
			lambda2 = a + (fibonacci_number(n - k - 2) / fibonacci_number(n - k - 1)) * (b - a);
			f2 = f(x + lambda2 * S);
		}
		length = b - a;
	}
	lambda2 = lambda1 + eps;
	f2 = f(x + lambda2 * S);

	if (f1 <= f2)
		b = lambda1;
	else
		a = lambda1;

	res.value = (a + b) / 2.0;
	res.call_to_func = calls_to_func;
	return res;
}