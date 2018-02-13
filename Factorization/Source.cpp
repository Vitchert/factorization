#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "BigInt.h"

using namespace std;

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

int mil_rab(BigInt n) // 1-prostoe
{

	int k, s = 0, buf = 0;;
	BigInt t, tt, a, x, i;
	k = 10; //к-во итераций

	if (n <= 3)
	{
		return 1;
	}

	t = n - 1;

	while (t % 2 == 0) {  //представление n-1 = (2^s)t, где t нечетно
		t = t / 2;
		s++;
	}

	for (int i = 0; i < k; i++)
	{
		a = rand();
		a = a % (n - 2) + 2; //случайное число из [2,n-1]

		x = powmod(a, t, n);
		if ((x == 1) || (x == n - 1)) continue;
		for (int j = 1; j < s; j++)
		{
			//x = pow(x, 2) % n;
			x = powmod(x, 2, n);
			if (x == 1)
			{
				return 0;
			}
			if (x == n - 1) break;
		}
		if (x != n - 1)
		{
			return 0;
		}
	}
	return 1;
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

int Miller(BigInt n)
{
	BigInt f = sqrt(sqrt(sqrt(n))) * 2;
	BigInt a = 2;
	BigInt buf;
	bool fl = false;
	int v, k = 0;
	for (a; a <= f; a = a + 1)
	{
		if (n%a == 0) return 0;
		if (powmod(a, n - 1, n) != 1) return 0;
		buf = 2;
		k = 0;
		while (buf < n)
		{
			if ((n - 1) % buf == 0) v = k;
			k++;
			buf = buf * 2;
		}
		for (k = 1; k <= v; k++)
		{
			buf = powmod(a, (n - 1) / pow((BigInt)2, k), n);
			if (buf == 0) buf = n - 1;
			else buf = buf - 1;
			//buf = (pow(a, ((n - 1) / pow((BigInt)2, k))) - 1)%n;
			buf = Nod(buf, n);
			if ((buf > 1) && (buf < n))
			{
				fl = true;
				break;
			}
		}
		if (fl == true) return 0;
	}
	return 1;
}

int main(int argc, char* argv[])
{
	BigInt n;
	//string str = argv[2];
	string str = "18014398241046523";
	n.input(str);

	if (!mil_rab(n))
	{
		cout << 0 << endl;
		return 0;
	}
	else
		if (!Miller(n))
		{
			cout << 0 << endl;
			return 0;
		}
		else
		{
			cout << 1 << endl;
			return 1;
		}
}

int QS(int n) {
	//step 1
	int P = exp(sqrt(log(n)*log(log(n))));
	int A = P*P / 2;
	//step 2
	vector<int> listing(A);
	int sqrtN = sqrt(n);
	for (int t = 0; t < A; ++t)
		listing[t] = (sqrtN + t + 1)*(sqrtN + t + 1) - n;

}