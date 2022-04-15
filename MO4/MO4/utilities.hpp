#pragma once
#ifndef UTILITIES_HPP
#define UTILITIES_HPP

#include <iostream>
#include <vector>
#include <functional>
#include <fstream>
#include <iomanip>
#include <random>

typedef std::function<double(const std::vector<double>&)> function2D;

// �������� ��������
inline std::vector<double> operator+(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<double> res = a;
	for (uint32_t i = 0; i < b.size(); i++)
		res[i] += b[i];
	return res;
}

// ��������� ��������
inline std::vector<double> operator-(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<double> res = a;
	for (uint32_t i = 0; i < b.size(); i++)
		res[i] -= b[i];
	return res;
}

// ��������� �������� (��������� ������������)
inline double operator*(const std::vector<double>& a, const std::vector<double>& b)
{
	double scalar = 0.0;
	for (uint32_t i = 0; i < b.size(); i++)
		scalar += a[i] * b[i];
	return scalar;
}

// �������� ��������� �� ������
inline std::vector<double> operator*(const double& c, const std::vector<double>& a)
{
	std::vector<double> res = a;
	for (auto& it : res)
		it *= c;
	return res;
}

// �������� ������ �� ����������������� ������
inline void vec_to_vecT(const std::vector<double>& a, const std::vector<double>& b, std::vector<std::vector<double>>& matrix)
{
	matrix.clear();

	matrix.resize(a.size());
	for (auto& it : matrix)
		it.resize(a.size());

	for (uint32_t i = 0; i < a.size(); i++)
		for (uint32_t j = 0; j < a.size(); j++)
			matrix[i][j] += a[i] * b[j];
}

// ����� �������
inline double norm(const std::vector<double>& a)
{
	return sqrt(a * a);
}

// �������� ������� �� ������
inline std::vector<double> operator*(const std::vector<std::vector<double>>& matrix, const std::vector<double>& vector)
{
	std::vector<double> result(vector.size());

	for (uint32_t i = 0; i < matrix.size(); i++)
		for (uint32_t j = 0; j < matrix[i].size(); j++)
			result[i] += matrix[i][j] * vector[j];

	return result;
}

// ��������� ������ �� �����
inline std::vector<double> operator/(std::vector<double>& vector, const double c)
{
	auto res = vector;

	for (uint32_t i = 0; i < vector.size(); i++)
		res[i] /= c;

	return res;
}

// ����� ������� � �����
inline std::ostream& operator <<(std::ostream& out, std::vector<double>& a)
{
	out << "(" << a[0] << ", " << a[1] << ")";

	return out;
}

// �������� ��� �������� � ������� �� ���������������
inline std::vector<double> operator-(const std::vector<double>& a)
{
	std::vector<double> neg_a = a;
	for (auto& it : neg_a)
		it = -it;
	return neg_a;
}

#endif