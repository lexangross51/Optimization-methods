#pragma once
#ifndef UTILITIES_H
#define UTILITIES_H

#include <iostream>
#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>

typedef std::function<double(const std::vector<double>&)> VectorFunc;

// Посчитать градиент функции
inline void grad(const VectorFunc& f, std::vector<double>& x, std::vector<double>& gradient)
{
	gradient.clear();
	gradient.resize(x.size());

	const double h = 1e-14;

	const double fx = f(x);

	for (uint32_t i = 0; i < x.size(); i++)
	{
		x[i] += h;
		gradient[i] = (f(x) - fx) / h;
		x[i] -= h;
	}
}

// Сложение векторов
inline std::vector<double> operator+(const std::vector<double>&a, const std::vector<double>&b)
{
	std::vector<double> res = a;
	for (uint32_t i = 0; i < b.size(); i++)
		res[i] += b[i];
	return res;
}

// Вычитание векторов
inline std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<double> res = a;
	for (uint32_t i = 0; i < b.size(); i++)
		res[i] -= b[i];
	return res;
}

// Умножение векторов (скалярное произведение)
inline double operator*(const std::vector<double>& a, const std::vector<double>& b)
{
	double scalar = 0.0;
	for (uint32_t i = 0; i < b.size(); i++)
		scalar += a[i] * b[i];
	return scalar;
}

// Умножить константу на вектор
inline std::vector<double> operator*(const double& c, const std::vector<double>& a)
{
	std::vector<double> res = a;
	for (auto& it : res)
		it *= c;
	return res;
}

// Умножить вектор на транспонированный вектор
inline void vec_to_vecT(const std::vector<double>& a, const std::vector<double>& b, std::vector<std::vector<double>>& matrix)
{
	matrix.clear();

	matrix.resize(a.size());
	for (auto& it : matrix)
		it.resize(a.size());

	for (uint32_t i = 0; i < a.size(); i++)
		for (uint32_t j = 0; j < a.size(); j++)
			matrix[i][j] = a[i] * b[j];
}

// Норма вектора
inline double norm(const std::vector<double>& a)
{
	return sqrt(a * a);
}

// Умножить матрицу на вектор
inline std::vector<double> operator*(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector)
{
	std::vector<double> result(vector.size());

	for (uint32_t i = 0; i < matrix.size(); i++)
		for (uint32_t j = 0; j < matrix[i].size(); j++)
			result[i] += matrix[i][j] * vector[j];

	return result;
}

// Разделить матрицу на число
inline void operator/(std::vector<std::vector<double>>& matrix, const double c)
{
	for (uint32_t i = 0; i < matrix.size(); i++)
		for (uint32_t j = 0; j < matrix[i].size(); j++)
			matrix[i][j] /= c;
}

// Сложить две матрицы
inline std::vector<std::vector<double>> operator+(const std::vector<std::vector<double>>& A, const std::vector<std::vector<double>>& B)
{
	std::vector<std::vector<double>> res = A;

	for (uint32_t i = 0; i < A.size(); i++)
		for (uint32_t j = 0; j < A.size(); j++)
			res[i][j] += B[i][j];

	return res;
}

// Вывод вектора в поток
inline std::ofstream& operator <<(std::ofstream& out, std::vector<double>& a)
{
	out << "(" << a[0] << ", " << a[1] << ")";
	
	return out;
}

// Вывод матрицы в поток
inline std::ofstream& operator <<(std::ofstream& out, std::vector<std::vector<double>>& A)
{
	out << A[0][0] << "\t" << A[0][1] << std::endl;
	out << A[1][0] << "\t" << A[1][1] << std::endl;

	return out;
}

// Изменить все значения в векторе на противоположные
inline std::vector<double> operator-(const std::vector<double>& a)
{
	std::vector<double> neg_a = a;
	for (auto& it : neg_a)
		it = -it;
	return neg_a;
}

#endif