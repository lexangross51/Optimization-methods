#include "methods.hpp"

uint32_t calls_to_func;

const uint8_t DIMENSION = 2;

int main()
{
	function2D f = [](const std::vector<double>& x)
	{
		calls_to_func++;

		return -(
			1.0  / (1 + (x[0] - 0) * (x[0] - 0) + (x[1] + 1)  * (x[1] + 1))  +
			2.0	 / (1 + (x[0] - 0) * (x[0] - 0) + (x[1] + 4)  * (x[1] + 4))  +
			10.0 / (1 + (x[0] - 3) * (x[0] - 3) + (x[1] + 2)  * (x[1] + 2))  +
			5.0  / (1 + (x[0] + 7) * (x[0] + 7) + (x[1] + 6)  * (x[1] + 6))  +
			7.0  / (1 + (x[0] - 6) * (x[0] - 6) + (x[1] + 10) * (x[1] + 10)) +
			9.0  / (1 + (x[0] - 6) * (x[0] - 6) + (x[1] - 1)  * (x[1] - 1))
		);
	};

	// ====================================================================================
	// ================================== ИССЛЕДОВАНИЯ ====================================
	// ====================================================================================
	statistical_search stat_search(DIMENSION);

	std::ofstream table("table.txt", std::ios::app);

	if (table.is_open())
	{
		#pragma region Исследование №1
		{
			table << "SIMPLE RANDOM SEARCH" << std::endl;
			table << std::left
			<< std::setw(7) << "eps" << "|"
			<< std::setw(10) << "P" << "|"
			<< std::setw(10) << "N" << "|"
			<< std::setw(22) << "(x, y)" << "|"
			<< std::setw(10) << "f(x, y)" << "|"
			<< std::endl;
			table << "--------------------------------";
			table << "--------------------------------" << std::endl;

			interval x_int(-10, 10);
			interval y_int(-10, 10);

			double epsilons[] = { 1, 1e-1, 1e-2 };
			double probabilities[] = { 0.1, 0.4, 0.7, 0.9 };


			for (const auto& eps : epsilons)
			{
				for (const auto& P : probabilities)
				{
					auto res = stat_search.simple_random_search(f, x_int, y_int, eps, P);

					table << std::left
						<< std::setw(7) << res.eps << "|"
						<< std::setw(10) << res.P << "|"
						<< std::setw(10) << res.N << "|"
						<< std::setw(10) << res.x[0] << ", "
						<< std::setw(10) << res.x[1] << "|"
						<< std::setw(10) << res.fx << "|"
						<< std::endl;
				}
				table << "--------------------------------";
				table << "--------------------------------" << std::endl;
			}
			table << std::endl << std::endl;
		}
		#pragma endregion


		#pragma region Исследование №2
		{
			interval x_int(-10, 10);
			interval y_int(-10, 10);

			auto print_global_search_result = [&](uint32_t num) -> void
			{
				table << "GLOBAL SEARCH №" << num << std::endl;
				table << std::left
					<< std::setw(7) << "eps" << "|"
					<< std::setw(10) << "calls to f" << "|"
					<< std::setw(10) << "m" << "|"
					<< std::setw(22) << "(x, y)" << "|"
					<< std::setw(10) << "f(x, y)" << "|"
					<< std::endl;
				table << "--------------------------------";
				table << "--------------------------------" << std::endl;

				double epsilons[] = { 1e-5, 1e-6, 1e-7, 1e-8 };

				for (const auto& eps : epsilons)
				{
					for (uint32_t m = 10; m < 50; m += 10)
					{
						statistical_search::result res;

						if (num == 1) res = stat_search.global_search_1(f, eps, m);
						else if (num == 2) res = stat_search.global_search_2(f, eps, m);
						else res = stat_search.global_search_3(f, x_int, y_int, eps, m);

						table << std::left
							<< std::setw(7) << res.eps << "|"
							<< std::setw(10) << res.call_to_func << "|"
							<< std::setw(10) << m << "|"
							<< std::setw(10) << res.x[0] << ", "
							<< std::setw(10) << res.x[1] << "|"
							<< std::setw(10) << res.fx << "|"
							<< std::endl;
					}
					table << "--------------------------------";
					table << "--------------------------------" << std::endl;
				}
				table << std::endl << std::endl;
			};		

			print_global_search_result(1);
			print_global_search_result(2);
			print_global_search_result(3);
		}

		#pragma endregion

		table.close();
	}

	return 0;
}