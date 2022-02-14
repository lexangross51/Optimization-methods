#include <iostream>
#include <functional>
#include <fstream>
#include <iomanip>

typedef std::function<double(double x)> func1D;

// Вывод в файлы для метода дихотомии, золотого сечения и фибоначчи
void print(std::ofstream& out, uint32_t iter, double x1, double x2, double f1, double f2, double a, double b, double a_prev, double b_prev)
{
    out.precision(10);

    if (iter == 1)
    {
        out << std::left << std::setw(6) << "iter"
            << std::setw(15) << "|x1" << std::setw(15) << "|x2"
            << std::setw(20) << "|f1" << std::setw(20) << "|f2"
            << std::setw(15) << "|a_i" << std::setw(15) << "|b_i"
            << std::setw(20) << "|b_i - a_i"
            << std::setw(15)  << "|(b_prev - a_prev)/(b - a)"
            << std::endl;
        out << "--------------------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
    }

    out << std::left << std::setw(6) << iter
        << "|" << std::setw(14) << x1 << "|" << std::setw(14) << x2
        << "|" << std::setw(19) << f1 << "|" << std::setw(19) << f2
        << "|" << std::setw(14) << a_prev << "|" << std::setw(14) << b_prev
        << "|" << std::setw(19) << b_prev - a_prev
        << "|" << std::setw(14) << (b_prev - a_prev) / (b - a)
        << std::endl;
}

//============================================================================
// Метод дихотомии
void dichotomy(const func1D &f, double a, double b, double eps)
{
    double x1 = 0, x2 = 0;
    double f_x1, f_x2;

    double a_prev, b_prev;

    uint32_t call_to_func = 0;
    uint32_t iter = 0;

    double delta = eps / 2.0;

    std::ofstream out("dichotomy.txt", std::ios::app);

    for ( ; abs(b - a) >= eps; iter++)
    {
        a_prev = a; 
        b_prev = b;

        x1 = (a + b - delta) / 2.0;
        x2 = (a + b + delta) / 2.0;

        f_x1 = f(x1);
        f_x2 = f(x2);
        call_to_func += 2;

        if (f_x1 > f_x2)
            a = x1;
        else
            b = x2;

        print(out, iter + 1, x1, x2, f_x1, f_x2, a, b, a_prev, b_prev);
    }

    out << "call to function: " << call_to_func
        << std::endl << std::endl << std::endl;

    out.close();
}

//============================================================================
// Метод золотого сечения
void golden_ratio(const func1D& f, double a, double b, double eps)
{
    uint32_t call_to_func = 0;
    uint32_t iter = 0;

    double delta = eps / 2.0;

    std::ofstream out("golden_ratio.txt", std::ios::app);

    double x1 = a + (3.0 - sqrt(5.0)) / 2.0 * (b - a);
    double x2 = a + (sqrt(5.0) - 1.0) / 2.0 * (b - a);

    double a_prev, b_prev;

    double f_x1 = f(x1);
    double f_x2 = f(x2);
    call_to_func += 2;

    for (; abs(b - a) >= eps; iter++)
    {
        a_prev = a;
        b_prev = b;

        if (f_x1 > f_x2)
        {
            a = x1;
            x1 = x2;
            f_x1 = f_x2;
            x2 = a + (sqrt(5.0) - 1.0) / 2.0 * (b - a);
            f_x2 = f(x2);
            call_to_func++;
        }
        else
        {
            b = x2;
            x2 = x1;
            f_x2 = f_x1;
            x1 = a + (3.0 - sqrt(5.0)) / 2.0 * (b - a);
            f_x1 = f(x1);
            call_to_func++;
        }

        print(out, iter + 1, x1, x2, f_x1, f_x2, a, b, a_prev, b_prev);
    }

    out << "call to function: " << call_to_func
        << std::endl << std::endl << std::endl;

    out.close();
}

//============================================================================
// Поиск интервала, содержащего минимум функции
void find_interval(const func1D& f, double a, double b, double x0, double eps)
{
    std::ofstream out("interval.txt", std::ios::app);
    out.precision(10);

    double delta = eps / 2.0;

    double x, xk_1, xk1;
    double h;

    double f1, f2;

    f1 = f(x0);
    f2 = f(x0 + delta);

    if (f1 > f2)
    {
        xk1 = x0 + delta;
        h = delta;
    }
    else
    {
        xk1 = x0 - delta;
        h = -delta;
    }
    
    x = x0;
    f2 = f(xk1);

    out << std::left << std::setw(15) << "x"
        << std::setw(15) << "f(x)" << std::endl;
    do
    {
        xk_1 = x;
        x = xk1;
        f1 = f2;

        h *= 2;
        xk1 = x + h;
        f2 = f(xk1);

        out << std::left << std::setw(15) << x
            << std::setw(15) << f1 << std::endl;

    } while (f1 > f2);

    out << "Interval with a minimum of a function: [" << xk_1 << ", " << xk1 << "]" << std::endl << std::endl;
    out.close();
}

//============================================================================
// Посчитать число Фибоначчи
double fibonacci_number(int n)
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

// Метод Фибоначчи
void fibonacci(const func1D& f, double a, double b, double eps)
{
    std::ofstream out("fibonacci.txt", std::ios::app);

    uint32_t call_to_func = 0;

    int n = 0; 
    for (; fibonacci_number(n + 2) < (b - a) / eps; n++);

    double Fn2 = fibonacci_number(n + 2);
    double x1 = a + fibonacci_number(n) / Fn2 * (b - a);
    double x2 = a + fibonacci_number(n + 1) / Fn2 * (b - a);

    double f_x1 = f(x1);
    double f_x2 = f(x2);
    call_to_func += 2;

    double a_prev, b_prev;

    double length = b - a;

    for (uint32_t k = 1; k <= n; k++)
    {
        a_prev = a;
        b_prev = b;

        if (f_x1 > f_x2)
        {
            a = x1;
            x1 = x2;
            f_x1 = f_x2;
            x2 = a + fibonacci_number(n - k + 2) / Fn2 * length;
            f_x2 = f(x2);
            call_to_func++;
        }
        else
        {
            b = x2;
            x2 = x1;
            f_x2 = f_x1;
            x1 = a + fibonacci_number(n - k + 1) / Fn2 * length;
            f_x1 = f(x1);
            call_to_func++;
        }

        print(out, k, x1, x2, f_x1, f_x2, a, b, a_prev, b_prev);
    }

    out << "call to function: " << call_to_func
        << std::endl << std::endl << std::endl;

    out.close();
}

int main()
{
    func1D f = [](double x) { return (x - 5) * (x - 5); };

    double a = -2, b = 20, eps = 1e-7;

    double epsilons[] = { 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7 };
    
    for (const auto& i : epsilons)
    {
        dichotomy(f, a, b, i);

        golden_ratio(f, a, b, i);

        find_interval(f, a, b, 0, i);

        fibonacci(f, a, b, i);
    }

    return 0;
}