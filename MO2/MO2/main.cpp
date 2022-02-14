#include "methods.h"

uint32_t calls_to_func;

int main()
{
	VectorFunc functions[] =
	{
		// Функция по варианту
		[](const std::vector<double>& x)
		{
			calls_to_func++;
			return	3 * exp(-(pow((x[0] - 2) / 1.0, 2) - pow((x[1] - 3) / 2.0, 2))) +
					1 * exp(-(pow((x[0] - 1) / 2.0, 2) - pow((x[1] - 1) / 1.0, 2)));
		},

		// Квадратичная функция
		[](const std::vector<double>& x)
		{
			calls_to_func++;
			return 100 * pow(x[1] - x[0], 2) + pow(1 - x[0], 2);
		},

		// Функция Розенброка
		[](const std::vector<double>& x)
		{
			calls_to_func++;
			return 100 * pow(x[1] - x[0] * x[0], 2) + pow(1 - x[0], 2);
		}
	};

	descent_methods methods(2);

	// Начальное приближение
	std::vector<double> x0 = {0, 0};

	double epsilons[] = {1e-3, 1e-4, 1e-5, 1e-6, 1e-7};
	double eps = 1e-7;
	
	//for (const auto& eps : epsilons)
	{
		auto res1 = methods.broyden(functions[2], x0, eps);
		res1.save("broyden.txt");
		x0 = { 0, 0 };

		auto res2 = methods.CGMFR(functions[2], x0, eps);
		res2.save("CGMFR.txt");
		x0 = { 0, 0 };
	}

	return 0;
}