#pragma once

#include <vector>
#include <string>
#include <algorithm>
using namespace std;
class BigInt
{
public:
	static const int osn = 10000; // или 100, 1000 и т.д.
	static const int len = 4;////////Цифра в 4 разряда

	vector <int> digits;
	BigInt() {};
	BigInt(int i);

	void input(string str);

	void output();

	void MultOsn();

	int& operator [] (const int i);

	BigInt operator + (const int &b) const;

	BigInt operator + (const BigInt &b) const;

	BigInt operator - (const BigInt &b) const;

	BigInt operator * (const int &n) const;

	BigInt  operator * (const BigInt &b) const;

	BigInt operator / (const int &n) const;

	BigInt operator / (const BigInt &b) const;

	int operator % (const int &n) const;

	BigInt operator % (const BigInt &b) const;

	bool operator == (const BigInt &b) const;

	bool operator != (const BigInt &b) const;

	bool operator >(const BigInt &b) const;

	bool operator < (const BigInt &b) const;

	bool operator <= (const BigInt &b) const;
};