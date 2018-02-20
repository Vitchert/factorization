#include "BigInt.h"
#include <string>
#include <iostream>
using namespace std;
BigInt::BigInt(int num)
{
	digits.clear();
	digits.push_back(num % osn);
	num = num / osn;
	if (num > 0)
	{
		while (num > osn)
		{
			digits.push_back(num % osn);
			num = num / osn;
		}
		digits.push_back(num % osn);
	}
	digits.shrink_to_fit();
}

void BigInt::input(string str)
{
	digits.clear();
	for (int i = str.size() - 1; i >= 0; i -= len)
	{
		int start = i - len + 1;
		if (start<0) start = 0;
		string dig = str.substr(start, i - start + 1);
		digits.push_back(atoi(dig.c_str()));
	}
	digits.shrink_to_fit();
}

void BigInt::output()
{
	int i = digits.size() - 1;
	if (i < 0)
	{
		cout << 0 << endl;
		return;
	}
	printf("%d", digits[i]);
	for (i = i - 1; i >= 0; i--)
		printf("%.4d", digits[i]);
	cout << endl;
}

void BigInt::MultOsn()
{
	if (digits.size() == 0) return;
	digits.push_back(digits.back());
	for (int i = digits.size() - 2; i >= 1; i--)
		digits[i] = digits[i - 1];
}

int& BigInt::operator [] (const int i)
{
	if ((i >= (int)digits.size()) || (i < 0)) cout << "Out of borders!\n";
	return digits[i];
}

BigInt BigInt::operator + (const int &b) const
{
	BigInt res = *this;
	int pos = 0, num = b;
	int size = res.digits.size();
	if (size == 0) res.digits.push_back(num%osn);
	else res[0] += num%osn;
	num /= osn;
	while (num)
	{
		if (pos<res.digits.size())  res.digits.push_back(num%osn);
		else res[pos + 1] += num%osn;
		pos++;
		num /= osn;
	}
	pos = 0;
	size = res.digits.size();
	while ((pos<size - 1) && (res[pos] >= osn))
	{
		res[pos + 1]++;
		res[pos++] -= osn;
	}
	if (res[pos] >= osn)
	{
		res[pos] -= osn;
		res.digits.push_back(1);
	}
	return res;
}

BigInt BigInt::operator + (const BigInt &b) const
{
	BigInt res;
	int size = max(this->digits.size(), b.digits.size());
	int thissize = this->digits.size();
	int bsize = b.digits.size();
	int r = 0, na, nb;
	for (int i = 0; i<size || r; i++)
	{
		if (i < thissize) na = this->digits[i];
		else na = 0;
		if (i < bsize) nb = b.digits[i];
		else nb = 0;

		res.digits.push_back(na + nb + r);
		if (res.digits[i] >= osn)
		{
			res.digits[i] -= osn;
			r = 1;
		}
		else
			r = 0;
	}
	return res;
}

BigInt BigInt::operator - (const BigInt &b) const
{
	BigInt res = *this;
	int r = 0, nb;
	int size = res.digits.size();
	int bsize = b.digits.size();
	for (int i = 0; i<size; i++)
	{
		if (i < bsize) nb = b.digits[i];
		else nb = 0;
		res.digits[i] -= nb + r;
		if (res.digits[i]<0)
		{
			res.digits[i] += osn;
			res.digits[i + 1]--;
		}
	}
	int pos = size - 1;
	while ((pos >= 0) && !res.digits[pos])
	{
		res.digits.pop_back();
		pos--;
	}
	return res;
}

BigInt BigInt::operator * (const int &n) const //Доделать переполнение при инте больше 4 знаков
{
	BigInt res;
	int size = this->digits.size();
	int an;
	if (size == 0)return res;

	int r = 0;
	for (int i = 0; i<size || r; i++)
	{
		if (i < size) an = this->digits[i];
		else an = 0;
		res.digits.push_back(an * n + r);
		r = res.digits[i] / osn;
		res.digits[i] -= r*osn;
	}
	return res;
}

BigInt BigInt::operator * (const BigInt &b) const
{
	BigInt res;
	int thissize = this->digits.size();
	int bsize = b.digits.size();
	int nb;
	res.digits.resize(thissize + bsize);
	for (int i = 0; i<thissize; i++)
	{
		int r = 0;
		for (int j = 0; j<bsize || r; j++)
		{
			if (j < bsize) nb = b.digits[j];
			else nb = 0;
			res.digits[i + j] += this->digits[i] * nb + r;
			r = res.digits[i + j] / osn;
			res.digits[i + j] -= r*osn;
		}
	}
	int pos = thissize + bsize;
	while (pos>0 && !res.digits[pos - 1])
	{
		res.digits.pop_back();
		pos--;
	}
	return res;
}

BigInt BigInt::operator / (const int &n) const //Проверит переполнение если больше 4 знаков
{
	BigInt res;
	int thissize = this->digits.size();
	if (thissize == 0)
	{
		return res;
	}
	int ost = 0;
	for (int i = thissize - 1; i >= 0; i--)
	{
		int cur = ost * osn + this->digits[i];
		res.digits.push_back(cur / n);
		ost = cur % n;
	}
	reverse(res.digits.begin(), res.digits.end());
	int pos = res.digits.size();
	while (pos>0 && !res.digits[pos - 1])
	{
		res.digits.pop_back();
		pos--;
	}
	return res;
}

BigInt BigInt::operator / (const BigInt &b) const
{
	BigInt res;
	BigInt curValue;
	int thissize = this->digits.size();
	for (int i = thissize - 1; i >= 0; i--)
	{
		curValue.MultOsn(); // * osn
		if (curValue.digits.size() == 0) curValue.digits.push_back(this->digits[i]);
		else curValue.digits[0] = this->digits[i];
		// подбираем максимальное число x, такое что b * x <= curValue
		int x = 0;
		int l = 0, r = osn;
		while (l <= r)
		{
			int m = (l + r) >> 1;
			BigInt cur = b * m;
			if (cur <= curValue)
			{
				x = m;
				l = m + 1;
			}
			else
				r = m - 1;
		}
		res.digits.push_back(x);
		curValue = curValue - b * x;
	}
	// избавляемся от лидирующих нулей
	reverse(res.digits.begin(), res.digits.end());
	int pos = res.digits.size();
	while (pos>0 && !res.digits[pos - 1])
	{
		res.digits.pop_back();
		pos--;
	}

	return res;
}

int BigInt::operator % (const int &n) const
{
	int thissize = this->digits.size();
	int ost = 0;
	for (int i = thissize - 1; i >= 0; i--)
	{
		int cur = ost * osn + this->digits[i];
		ost = cur % n;
	}
	return ost;
}

BigInt BigInt::operator % (const BigInt &b) const
{
	BigInt curValue;
	int thissize = this->digits.size();
	int bsize = b.digits.size();
	if ((thissize == 1) && (bsize == 1)) return curValue = this->digits[0] % b.digits[0];
	if (bsize > thissize) return curValue = *this;
	for (int i = thissize - 1; i >= 0; i--)
	{
		curValue.MultOsn(); // * osn
		if (curValue.digits.size() == 0) curValue.digits.push_back(this->digits[i]);
		else curValue.digits[0] = this->digits[i];
		// подбираем максимальное число x, такое что b * x <= curValue
		int x = 0;
		int l = 0, r = osn;
		while (l <= r)
		{
			int m = (l + r) >> 1;
			BigInt cur = b * m;
			if (cur <= curValue)
			{
				x = m;
				l = m + 1;
			}
			else
				r = m - 1;
		}
		curValue = curValue - b * x;
	}
	return curValue;
}

bool BigInt::operator == (const BigInt &b) const
{
	if ((this->digits.size() == 0) && (b == 0))
		return true;

	int thissize = this->digits.size();
	int bsize = b.digits.size();

	if (thissize != bsize)
		return false;
	for (int i = 0; i<thissize; i++)
	{
		if (this->digits[i] != b.digits[i])
			return false;
	}
	return true;
}

bool BigInt::operator != (const BigInt &b) const
{
	if ((this->digits.size() == 0) && (b == 0))
		return false;

	int thissize = this->digits.size();
	int bsize = b.digits.size();

	if (thissize != bsize)
		return true;
	for (int i = 0; i<thissize; i++)
	{
		if (this->digits[i] != b.digits[i])
			return true;
	}
	return false;

}

bool  BigInt::operator >(const BigInt &b) const
{

	int thissize = this->digits.size();
	int bsize = b.digits.size();

	if (thissize != bsize)
		return thissize>bsize;
	for (int i = thissize - 1; i >= 0; i--)
	{
		if (this->digits[i] != b.digits[i])
			return this->digits[i]>b.digits[i];
	}
	return false;

}

bool BigInt::operator < (const BigInt &b) const
{

	int thissize = this->digits.size();
	int bsize = b.digits.size();

	if (thissize != bsize)
		return thissize<bsize;
	for (int i = thissize - 1; i >= 0; i--)
	{
		if (this->digits[i] != b.digits[i])
			return this->digits[i]<b.digits[i];
	}
	return false;

}

bool BigInt::operator <= (const BigInt &b) const
{
	if ((this->digits.size() == 0) && (b == 0))
		return true;

	int thissize = this->digits.size();
	int bsize = b.digits.size();

	if (thissize != bsize)
		return thissize<bsize;
	for (int i = thissize - 1; i >= 0; i--)
	{
		if (this->digits[i] != b.digits[i])
			return this->digits[i]<b.digits[i];
	}
	return true;

}

BigInt pow(const BigInt &a, const int &N)
{
	BigInt res = 1, cur = a;
	int n = N;
	if (n == 0) return 1;
	while (n)
	{
		if (n & 1)
			res = res * cur;
		cur = cur * cur;
		n >>= 1;
	}
	return res;
}

BigInt pow(BigInt aa, BigInt kk)
{
	BigInt res = 1, a = aa, k = kk;

	while (k != 0) {
		if (k % 2 == 0) {
			k = k / 2;
			a = a * a;
		}
		else {
			k = k - 1;
			res = res * a;
		}
	}
	return res;
}

BigInt powmod(BigInt aa, BigInt kk, BigInt nn)
{
	BigInt res = 1, a = aa, k = kk, n = nn;

	while (k != 0) {
		if (k % 2 == 0) {
			k = k / 2;
			a = (a*a) % n;
		}
		else {
			k = k - 1;
			res = (res*a) % n;
		}
	}
	return res;
}

BigInt sqrt(BigInt n)
{
	BigInt cur;
	int nsize = n.digits.size();
	int pos = (nsize + 1) / 2;
	pos--;
	for (int i = 0; i <= pos; i++)
		cur.digits.push_back(0);
	while (pos >= 0)
	{
		int l = 0, r = BigInt::osn;
		int curDigit = 0;
		while (l <= r) // подбираем текущую цифру
		{
			int m = (l + r) >> 1;
			cur.digits[pos] = m;
			if ((cur * cur) <= n)
			{
				curDigit = m;
				l = m + 1;
			}
			else
				r = m - 1;
		}
		cur.digits[pos] = curDigit;
		pos--;
	}
	// избавляемся от лидирующих нулей
	if (cur.digits.size() != 1 && !cur.digits[cur.digits.size() - 1])
		cur.digits.pop_back();
	return cur;
}

BigInt Nuton(BigInt n)
{
	BigInt  A = n;
	BigInt x1, x2, razn;
	x1 = n / 2;
	x2 = (x1 + A / x1) / 2;
	razn = x1 - x2;
	while (razn > 1)
	{
		x1 = x2;
		x2 = (x1 + A / x1) / 2;
		razn = x1 - x2;
	}
	return x2;
}

BigInt Nod(BigInt a, BigInt b)
{
	if ((a.digits.size() == 0) || (b.digits.size() == 0)) return 0;

	while (b != 0) {
		a = a % b;
		swap(a, b);
	}
	return a;

}