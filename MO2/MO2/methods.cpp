#include "methods.h"

extern uint32_t calls_to_func;

//------------------------------------------------------------------------------------------------------------------
descent_methods::descent_methods(const uint32_t dimension)
{
	x.resize(dimension);

	etta_k =
	{
		{0.0, 0.0},
		{0.0, 0.0}
	};

	grad_f.resize(dimension);
	S.resize(dimension);
}

// Метод Бройдена
descent_methods::result descent_methods::broyden(const VectorFunc& f, std::vector<double>& x0, double eps)
{
	calls_to_func = 0;

	// Сюда будем выводить точки для отрисовки хода метода
	std::ofstream steps("steps_broyden.txt");
	std::ofstream table("broyden.txt", std::ios::app);

	result res;
	res.x0 = x0;
	res.eps = eps;

	// Точки для отрисовки хода
	steps << x0[0] << " " << x0[1] << std::endl;

	// Положительно определенная матрица
	etta_k =
	{
		{1.0, 0.0},
		{0.0, 1.0}
	};

	std::vector<std::vector<double>> delta_etta_k;

	grad(f, x0, grad_f);

	do
	{
		if (res.iter_cnt % 2 == 0)
		{
			etta_k =
			{
				{1.0, 0.0},
				{0.0, 1.0}
			};
		}
		S = etta_k * grad_f;

		// Найдем lambda_k: для этого воспользуемся методом Фибоначии одномерного поиска
		// Но сначала определим интервал, на котором достигает своего минимума функция
		auto res_int = one_dimensional_search_methods::find_interval(f, x0, S, eps);
		auto res_fib = one_dimensional_search_methods::fibonacci(f, x0, res_int.interval, S, eps);
		 
		double lambda_k = res_fib.value;

		auto dx = lambda_k * S;

		x = x0;

		x0 = x0 + dx;

		auto grad_f_prev = grad_f;
	
		grad(f, x0, grad_f);

		auto delta_grad_f = grad_f - grad_f_prev;

		auto tmp_vector = dx - etta_k * delta_grad_f;

		vec_to_vecT(tmp_vector, tmp_vector, delta_etta_k);

		double scalar_prod = (tmp_vector * delta_grad_f);

		if (scalar_prod == 0.0)
		{
			etta_k =
			{
				{1.0, 0.0},
				{0.0, 1.0}
			};
		}
		else
		{
			delta_etta_k / scalar_prod;
			etta_k = etta_k + delta_etta_k;
		}

		res.iter_cnt++;

		steps << x0[0] << " " << x0[1] << std::endl;

		print_iter(table, res.iter_cnt, x0, f, lambda_k);

	} while (norm(grad_f) > eps && abs(f(x0) - f(x)) > eps);

	res.call_to_func = calls_to_func - res.iter_cnt * 3;
	res.x = x0;
	res.fx = f(x0);

	steps.close();
	table.close();

	return res;
}

// Метод сопряженных градиентов с модификацией Флетчера-Ривса
descent_methods::result descent_methods::CGMFR(const VectorFunc& f, std::vector<double>& x0, double eps)
{
	calls_to_func = 0;

	std::ofstream steps("steps_CGMFR.txt");
	std::ofstream table("CGMFR.txt", std::ios::app);

	steps << x0[0] << " " << x0[1] << std::endl;

	result res;
	res.eps = eps;
	res.x0 = x0;

	double omega_k;

	while (true)
	{
		grad(f, x0, grad_f);
		
		S = -grad_f;

		if (norm(S) < eps)
			break;

		for (uint32_t k = 0; k <= 2; k++)
		{
			omega_k = 1.0 / (norm(grad_f) * norm(grad_f));

			// Найдем lambda_k: для этого воспользуемся методом Фибоначии одномерного поиска
			// Но сначала определим интервал, на котором достигает своего минимума функция
			auto res_int = one_dimensional_search_methods::find_interval(f, x0, S, eps);
			auto res_fib = one_dimensional_search_methods::fibonacci(f, x0, res_int.interval, S, eps);
			double lambda_k = res_fib.value;

			x = x0 + lambda_k * S;

			res.iter_cnt++;

			print_iter(table, res.iter_cnt, x0, f, lambda_k);

			grad(f, x, grad_f);

			omega_k *= norm(grad_f) * norm(grad_f);

			S = -grad_f + omega_k * S;
		
			std::swap(x, x0);
			steps << x0[0] << " " << x0[1] << std::endl;

			if (norm(S) < eps)
				break;
		}
	}

	res.call_to_func = calls_to_func - res.iter_cnt * 3;
	res.x = x0;
	res.fx = f(x0);

	steps.close();
	table.close();

	return res;
}

// Сохранить результат в файл
void descent_methods::result::save(std::string filename)
{
	FILE* res;

	if (!fopen_s(&res, filename.c_str(), "a"))
	{
		fprintf(res, "%-24s|%-12s|%-12s|%-13s|%-24s|%-12s\n",
			"x0", "eps", "iters cnt", "call to func", "x", "f(x)");
		fprintf(res, "------------------------------------");
		fprintf(res, "------------------------------------");
		fprintf(res, "------------------------------------\n");

		fprintf(res, "(%-10lf, %-10lf)|%-12e|%-12ld|%-13ld|(%-10lf, %-10lf)|%-12e\n\n\n",
			x0[0], x0[1], eps, iter_cnt, call_to_func, x[0], x[1], fx);

		fclose(res);
	}
}

// Напечатать строку в таблицу для одной итерации
void descent_methods::print_iter(std::ofstream& res, uint32_t iter,
	std::vector<double>& x0, const VectorFunc& f,
	double lambda)
{
	res << "Iteration #" << iter << std::endl;
	res << "(x_i, y_i): "; res << x0 << std::endl;
	res << "f(x_i, y_i): " << f(x0) << std::endl;
	res << "(s1, s2): "; res << S << std::endl;
	res << "lambda: " << lambda << std::endl;
	res << "|x_i - x_i-1|: " << abs(x0[0] - x[0]) << std::endl;
	res << "|y_i - y_i-1|: " << abs(x0[1] - x[1]) << std::endl;
	res << "|f_i - f_i-1|: " << abs(f(x0) - f(x)) << std::endl;
	double num = (x0[0] * S[0] + x0[1] * S[1]);
	double denum = (sqrt(x0[0] * x0[0] + x0[1] * x0[1]) * sqrt(S[0] * S[0] + S[1] * S[1]));

	double angle = acos(num / denum);
	res << "(x_i, y_i)^(s1, s2): " << angle << std::endl;
	res << "grad_f: "; res << grad_f << std::endl;
	res << "etta_k: " << std::endl; res << etta_k;
	res << "---------------------------------------------" << std::endl << std::endl;
}

//------------------------------------------------------------------------------------------------------------------
// Поиска интервала, содержащего минимум функции n переменных по направлению S
one_dimensional_search_methods::result one_dimensional_search_methods::find_interval(
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

// Посчиать n-ое число Фибоначчи
double one_dimensional_search_methods::fibonacci_number(uint32_t n)
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
one_dimensional_search_methods::result one_dimensional_search_methods::fibonacci(
	const VectorFunc& f, std::vector<double>& x, interval& interval, std::vector<double>& S, double eps)
{
	result res;

	uint32_t n = 0;

	double a = interval.a;
	double b = interval.b;
	double length = b - a;

	for (; fibonacci_number(n + 2) < (b - a) / eps; n++);

	double Fn2 = fibonacci_number(n + 2);
	double lambda1 = a + fibonacci_number(n) / Fn2 * length;
	double lambda2 = a + fibonacci_number(n + 1) / Fn2 * length;

	double f1 = f(x + lambda1 * S);
	double f2 = f(x + lambda2 * S);

	for (uint32_t k = 1; k <= n; k++)
	{
		if (f1 > f2)
		{
			a = lambda1;
			lambda1 = lambda2;
			f1 = f2;
			lambda2 = a + fibonacci_number(n - k + 2) / Fn2 * length;
			f2 = f(x + lambda2 * S);
		}
		else
		{
			b = lambda2;
			lambda2 = lambda1;
			f2 = f1;
			lambda1 = a + fibonacci_number(n - k + 1) / Fn2 * length;
			f1 = f(x + lambda1 * S);
		}
	}
	
	res.value = (a + b) / 2.0;
	res.call_to_func = calls_to_func;
	return res;
}

// Метод дихотомии одномерного поиска
one_dimensional_search_methods::result one_dimensional_search_methods::dichotomy(
	const VectorFunc& f, std::vector<double>& x, interval& interval, std::vector<double>& S, double eps)
{
	result res;

	double a = interval.a;
	double b = interval.b;

	double delta = eps / 2.0;

	for (; abs(b - a) >= eps; )
	{
		double lambda1 = (a + b - delta) / 2.0;
		double lambda2 = (a + b + delta) / 2.0;

		double f_x1 = f(x + lambda1 * S);
		double f_x2 = f(x + lambda2 * S);

		res.iters_cnt++;

		if (f_x1 > f_x2)
			a = lambda1;
		else
			b = lambda2;
	}

	res.value = (a + b) / 2.0;
	res.call_to_func = calls_to_func;
	return res;
}

// Метод золотого сечения одномерного поиска
one_dimensional_search_methods::result one_dimensional_search_methods::golden_ratio(
	const VectorFunc& f, std::vector<double>& x, interval& interval, std::vector<double>& S, double eps)
{
	result res;

	double a = interval.a;
	double b = interval.b;

	double delta = eps / 2.0;

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