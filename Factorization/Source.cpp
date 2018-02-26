#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "AtkinSieve.cpp"
#include <map>
#include <set>
#include <mpirxx.h>
#include <cstdint>
#include <mpfr.h>
using namespace std;

int mil_rab(mpz_class n) // 1-prostoe
{

	mpz_class k, s = 0, buf = 0;;
	mpz_class t, tt, a, x, i;
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

	for (mpz_class i = 0; i < k; i++)
	{
		a = rand();
		a = a % (n - 2) + 2; //случайное число из [2,n-1]

		//x = powmod(a, t, n);
		mpz_powm(x.get_mpz_t(), a.get_mpz_t(), t.get_mpz_t(), n.get_mpz_t());
		if ((x == 1) || (x == n - 1)) continue;
		for (mpz_class j = 1; j < s; j++)
		{
			//x = pow(x, 2) % n;
			//x = powmod(x, 2, n);
			mpz_powm(x.get_mpz_t(), x.get_mpz_t(), mpz_class(2).get_mpz_t(), n.get_mpz_t());
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

int Miller(mpz_class n)
{
	mpz_class f = sqrt(sqrt(sqrt(n))) * 2;
	mpz_class a = 2;
	mpz_class buf;
	bool fl = false;
	mpz_class v, k = 0;
	for (a; a <= f; a = a + 1)
	{
		if (n%a == 0) return 0;
		mpz_class temp;
		mpz_powm(temp.get_mpz_t(), a.get_mpz_t(), mpz_class(n - 1).get_mpz_t(), n.get_mpz_t());
		if (temp != 1) return 0;
		
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
			//buf = powmod(a, (n - 1) / pow((int)2, k), n);
			mpz_class power;
			mpz_pow_ui(power.get_mpz_t(), mpz_class(2).get_mpz_t(), k.get_ui());
			power = (n - 1) / power;
			mpz_powm(buf.get_mpz_t(), a.get_mpz_t(), power.get_mpz_t(), n.get_mpz_t());
			if (buf == 0) buf = n - 1;
			else buf = buf - 1;

			//buf = Nod(buf, n);
			mpz_gcd(buf.get_mpz_t(), buf.get_mpz_t(), n.get_mpz_t());
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

mpz_class LegendreSymbol(mpz_class a, mpz_class p) {
	if (a >= p)
		a %= p;
	mpz_class res;// = powmod(a, (p - 1) / 2, p);
	mpz_powm(res.get_mpz_t(), a.get_mpz_t(), mpz_class((p-1)/2).get_mpz_t(), p.get_mpz_t());
	return res > 1 ? -1 : res;
}

mpz_class SqrtMod(mpz_class a, mpz_class p) {
	if (a >= p)
		a %= p;
	mpz_class n = 1;
	while (LegendreSymbol(n, p) >= 0) {
		++n;
	}
	mpz_class alpha = 0;
	mpz_class s = p - 1;
	while (s % 2 == 0) {
		++alpha;
		s /= 2;
	}
	mpz_class b;// = powmod(n, s, p);
	mpz_powm(b.get_mpz_t(), n.get_mpz_t(), s.get_mpz_t(), p.get_mpz_t());
	mpz_class r;// = powmod(a, (s + 1) / 2, p);
	mpz_powm(r.get_mpz_t(), a.get_mpz_t(), mpz_class((s + 1) / 2).get_mpz_t(), p.get_mpz_t());
	mpz_class rCalc;// = powmod(a, s, p);
	mpz_powm(rCalc.get_mpz_t(), a.get_mpz_t(), s.get_mpz_t(), p.get_mpz_t());
	mpz_class j;	
	mpz_class power;
	mpz_pow_ui(power.get_mpz_t(), mpz_class(2).get_mpz_t(), mpz_class(alpha - 2).get_ui());
	mpz_class check;// = powmod(rCalc, pow(2, alpha - 2), p);
	mpz_powm(check.get_mpz_t(), rCalc.get_mpz_t(), power.get_mpz_t(), p.get_mpz_t());
	if (check > 1)
		check -= p;
	if (check == 1) {
		j = 0;
	}else
		if (check == -1) {
			j = 1;
		}
		else
			cout << "check error" << endl;
	for (mpz_class i = 1; i < alpha - 1; ++i) {
		//check = powmod(powmod(powmod(b, j, p), 2, p) * rCalc, pow(2, alpha - i - 2), p);
		mpz_class base, exp;
		mpz_pow_ui(exp.get_mpz_t(), mpz_class(2).get_mpz_t(), mpz_class(alpha - i - 2).get_ui());
		mpz_powm(base.get_mpz_t(), b.get_mpz_t(), j.get_mpz_t(), p.get_mpz_t());
		mpz_powm(base.get_mpz_t(), base.get_mpz_t(), mpz_class(2).get_mpz_t(), p.get_mpz_t());
		base = base * rCalc;
		mpz_powm(check.get_mpz_t(), base.get_mpz_t(), exp.get_mpz_t(), p.get_mpz_t());
		if (check > 1)
			check -= p;
		if (check == 1) {
			//j = 0;
		}
		else
			if (check == -1) {
				mpz_class addition;
				mpz_pow_ui(addition.get_mpz_t(), mpz_class(2).get_mpz_t(), i.get_ui());
				//j += pow(2, i);
				j += addition;
			}
			else
				std::cout << "check error" << endl;
	}
	mpz_class result;
	mpz_powm(result.get_mpz_t(), b.get_mpz_t(), j.get_mpz_t(), p.get_mpz_t());
	//return (powmod(b, j, p) * r) % p;
	return (result * r) % p;
}

mpz_class gcd(mpz_class a, mpz_class b) {
	while (b != 0) {
		a %= b;
		swap(a, b);
	}
	return a;
}

mpz_class gcdex(mpz_class a, mpz_class b, mpz_class & x, mpz_class & y) {//–асширенный алгоритм евклида emaxx
	if (a == 0) {
		x = 0; y = 1;
		return b;
	}
	mpz_class x1, y1;
	mpz_class d = gcdex(b%a, a, x1, y1);
	x = y1 - (b / a) * x1;
	y = x1;
	return d;
}

mpz_class ReverseMod(mpz_class a, mpz_class m) {//обратное по модулю число
	mpz_class x = 0, y;
	mpz_class g = gcdex(a, m, x, y);
	if (g != 1)
		std::cout << "ReverseMod no solution";
	else {
		x = (x % m + m) % m;
	}
	return x;
}

mpz_class sqrtModPowerLiftingPlusOne(mpz_class root, mpz_class quadResidue, mpz_class mod, mpz_class modBase) {
	mpz_class newmod = mod*modBase;
	mpz_class yk = ReverseMod(2 * root, mod)% newmod;	
	quadResidue %= newmod;
	root = root - (((root*root) % newmod - quadResidue)*yk) % newmod;

	if (root < 0) {
		root = -root;
	}
	root %= newmod;
	//root = min(root, -(root - newmod));
	return root;
}
/*
mpz_class primeCheck() {
	mpz_class n;
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
*/
bool chekNuberInRange(mpz_class lowerBound, mpz_class upperBound, mpz_class t, mpz_class p, mpz_class A) {// check if there is range lower <=T<= upper that T = t1 or t2 ( mod p)
	if (A >= p) {
		return true;
	}
	mpz_class l = lowerBound % p;
	mpz_class u = upperBound % p;
	if (l >= u) {
		if (t >= u || t <= l)
			return true;
	}
	else {
		if (t >= l && t <= u)
			return true;
	}
	return false;
}

mpz_class Step4(mpz_class n, mpz_class p, mpz_class A, mpz_class& res1, mpz_class& res2) {
	mpz_class lowerBound = sqrt(n) + 1;
	mpz_class upperBound = sqrt(n) + A;

	mpz_class residue = n;
	mpz_class modBase = p;
	mpz_class currentMod = p;

	mpz_class root = SqrtMod(residue, currentMod);
	mpz_class root2 = -(root - modBase);
	res1 = root;
	res2 = root2;
	mpz_class b = 1;
	//while (root < lowerBound && root2 < lowerBound) {
	while (!chekNuberInRange(lowerBound, upperBound, root, currentMod,A) 
		&& !chekNuberInRange(lowerBound, upperBound, root2, currentMod, A)) {
		root = sqrtModPowerLiftingPlusOne(root, residue, currentMod, modBase);
		currentMod *= modBase;
		root2 = -(root - currentMod);
		++b;
		res1 = root;
		res2 = root2;
		//cout << "power " << currentMod;
		//cout << endl << root << endl;
	}
	//while ((root >= lowerBound && root <= upperBound)
	//	|| (root2 >= lowerBound && root2 <= upperBound)) {
	do {
		res1 = root;
		res2 = root2;
		root = sqrtModPowerLiftingPlusOne(root, residue, currentMod, modBase);
		currentMod *= modBase;
		root2 = -(root - currentMod);
		++b;
	} while (chekNuberInRange(lowerBound, upperBound, root, currentMod, A)
		|| chekNuberInRange(lowerBound, upperBound, root2, currentMod, A));
	return b - 1;
}

mpz_class liftRootEvenModulo(mpz_class root, mpz_class a, mpz_class p) {
	mpz_class i = ((root * root - a) / p) % 2;
	return root + i * (p / 2);
}

void step7(mpz_class n, mpz_class A, vector<mpz_class>& listingT, vector<mpz_class>& listingTSqr, vector<vector<mpz_class>>& exponentMatrix) {
	mpz_class sizeT = listingT.size();
	if (n % 8 != 1) {
		for (uint64_t j = 0; j < sizeT; ++j) {
			if (listingT[j] % 2 == 1) {
				exponentMatrix[0][j] = 1;
				listingTSqr[j] /= 2;
			}
		}
	}
	else {
		mpz_class lowerBound = sqrt(n) + 1;
		mpz_class upperBound = sqrt(n) + A;

		mpz_class residue = n;
		mpz_class modBase = 2;
		mpz_class currentMod = 2;

		mpz_class b = 3;

		mpz_class root1 = 1;
		mpz_class root2 = 3;
		mpz_class root3 = 5;
		mpz_class root4 = 7;
		mpz_class t1 = 1;
		mpz_class t2 = 3;
		mpz_class t3 = 5;
		mpz_class t4 = 7;
		while (!chekNuberInRange(lowerBound, upperBound, root1, currentMod, A)
			&& !chekNuberInRange(lowerBound, upperBound, root2, currentMod, A)
			&& !chekNuberInRange(lowerBound, upperBound, root3, currentMod, A)
			&& !chekNuberInRange(lowerBound, upperBound, root4, currentMod, A)) {
			t1 = root1;
			t2 = root2;
			t3 = root3;
			t4 = root4;
			root1 = liftRootEvenModulo(root1, n, currentMod);
			root2 = liftRootEvenModulo(root2, n, currentMod);
			root3 = liftRootEvenModulo(root3, n, currentMod);
			root4 = liftRootEvenModulo(root4, n, currentMod);
			currentMod *= modBase;
			++b;
		}
		do {
			t1 = root1;
			t2 = root2;
			t3 = root3;
			t4 = root4;
			root1 = liftRootEvenModulo(root1, n, currentMod);
			root2 = liftRootEvenModulo(root2, n, currentMod);
			root3 = liftRootEvenModulo(root3, n, currentMod);
			root4 = liftRootEvenModulo(root4, n, currentMod);
			currentMod *= modBase;
			++b;
		} while (chekNuberInRange(lowerBound, upperBound, root1, currentMod, A)
			|| chekNuberInRange(lowerBound, upperBound, root2, currentMod, A)
			|| chekNuberInRange(lowerBound, upperBound, root3, currentMod, A)
			|| chekNuberInRange(lowerBound, upperBound, root4, currentMod, A));

		--b;

		for (uint64_t j = 0; j < sizeT; ++j) {
			for (uint64_t k = 0; k < b; ++k) {
				mpz_class differ = abs(listingT[j] - t1);
				if (differ % (mpz_class)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t2);
				if (differ % (mpz_class)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t3);
				if (differ % (mpz_class)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t4);
				if (differ % (mpz_class)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
			}
		}
	}
}

void step45678(mpz_class n, mpz_class A, vector<mpz_class>& listingT, vector<mpz_class>& listingTSqr, vector<mpz_class>& factorBase, vector<vector<mpz_class>>& matrix) {
	//step 4
	uint64_t size = factorBase.size();
	uint64_t sizeT = listingT.size();
	vector<vector<mpz_class>> exponentMatrix(size, vector<mpz_class>(sizeT, 0));
	for (uint64_t i = 1; i < size; ++i) {//for each odd prime
		mpz_class t1, t2;
		mpz_class p = factorBase[i];
		mpz_class b = Step4(n, p, A, t1, t2);
		//step 5
		for (uint64_t j = 0; j < sizeT; ++j) {
			for (mpz_class k = 0; k < b; ++k) {
				mpz_class differ = abs(listingT[j] - t1);
				mpz_class mod;
				mpz_pow_ui(mod.get_mpz_t(), p.get_mpz_t(), mpz_class(k+1).get_ui());
				//if (differ % (mpz_class)pow(p, k + 1) == 0) {
				if (differ % mod == 0) {
					if (exponentMatrix[i][j] < k + 1) {
						exponentMatrix[i][j] = k + 1;
						//step 6
						listingTSqr[j] /= p;
					}
				}
				differ = abs(listingT[j] - t2);
				//if (differ % (mpz_class)pow(p, k + 1) == 0) {
				if (differ % mod == 0) {
					if (exponentMatrix[i][j] < k + 1) {
						exponentMatrix[i][j] = k + 1;
						//step 6
						listingTSqr[j] /= p;
					}
				}
			}
		}
	}
	//step 7
	step7(n, A, listingT, listingTSqr, exponentMatrix);
	//step 8
	uint64_t matrixSize = 0;
	for (uint64_t i = 0; i < sizeT; ++i) {
		if (listingTSqr[i] == 1)
			++matrixSize;
	}
	matrix = vector<vector<mpz_class>>(size + 1, vector<mpz_class>(matrixSize, 0));
	uint64_t pos = 0;
	for (uint64_t i = 0; i < sizeT; ++i) {
		if (listingTSqr[i] == 1) {
			matrix[0][pos] = listingT[i];
			for (uint64_t j = 0; j < size; ++j) {
				matrix[j + 1][pos] = exponentMatrix[j][i];
			}
			++pos;
		}
	}
}

void step9(vector<vector<mpz_class>>& matrix, vector<vector<mpz_class>>& factorizationMembers) {
	//step 9
	uint64_t cols = matrix.size();
	uint64_t rows = matrix[0].size();
	vector<vector<mpz_class>> gaussianEliminationMatrix(rows, vector<mpz_class>(cols, 0));
	for (uint64_t i = 0; i < rows; ++i) {
		gaussianEliminationMatrix[i][0] = matrix[0][i];
		for (uint64_t j = 1; j < cols; ++j) {
			gaussianEliminationMatrix[i][j] = matrix[j][i] % 2;
		}
	}
	//Gaussian elimination	
	map<mpz_class, multiset<mpz_class>> linearCompositionMembers;
	for (uint64_t i = 0; i < rows; ++i) {
		linearCompositionMembers.insert({ gaussianEliminationMatrix[i][0],{ gaussianEliminationMatrix[i][0] } });
	}

	uint64_t curpos = 1;
	for (uint64_t i = 0; (i < rows) && (curpos < cols); ++i, ++curpos) {
		uint64_t onePos = i;
		while (onePos < rows &&
			gaussianEliminationMatrix[onePos][curpos] != 1)
			++onePos;
		if (onePos == rows) {
			if (i)
				--i;
			continue;
		}
		if (onePos != i) {
			swap(gaussianEliminationMatrix[onePos], gaussianEliminationMatrix[i]);
		}

		for (uint64_t j = 0; j < rows; ++j) {
			if (j == i
				|| gaussianEliminationMatrix[j][curpos] == 0)
				continue;
			for (uint64_t k = 1; k < cols; ++k) {
				gaussianEliminationMatrix[j][k] = (gaussianEliminationMatrix[j][k] + gaussianEliminationMatrix[i][k]) % 2;
			}
			if (j>i)
				linearCompositionMembers[gaussianEliminationMatrix[j][0]].insert(
					linearCompositionMembers[gaussianEliminationMatrix[i][0]].begin(),
					linearCompositionMembers[gaussianEliminationMatrix[i][0]].end());

		}
	}
	//get list of posible combinations
	vector<multiset<mpz_class>> posibleFactorizationMembers;

	for (uint64_t i = 0; i < rows; ++i) {
		uint64_t j = 1;
		while (j < cols && gaussianEliminationMatrix[i][j] == 0)
			++j;
		if (j < cols)
			continue;
		posibleFactorizationMembers.push_back(linearCompositionMembers[gaussianEliminationMatrix[i][0]]);
	}
	//remove all even encounters
	mpz_class size = posibleFactorizationMembers.size();
	for (uint64_t i = 0; i < size; ++i) {
		vector<mpz_class> factors;
		while (posibleFactorizationMembers[i].size() > 0) {
			multiset<mpz_class>::iterator it = posibleFactorizationMembers[i].begin();
			if (posibleFactorizationMembers[i].count(*it) % 2 == 1)
				factors.push_back(*it);
			posibleFactorizationMembers[i].erase(*it);
		}
		factorizationMembers.push_back(factors);
	}
}

mpz_class QS(mpz_class n) {//n is odd
	//step 1
	string numberString = n.get_str();
	mpfr_rnd_t rnd = mpfr_get_default_rounding_mode();
	mpfr_t exp;
	mpfr_init_set_str(exp, numberString.c_str(), 10, rnd);
	mpfr_t logN;	
	mpfr_t loglogN;
	mpfr_init(logN);
	mpfr_init(loglogN);
	mpfr_log(logN, exp, rnd);
	mpfr_log(loglogN, logN, rnd);
	mpfr_mul(exp, logN, loglogN, rnd);
	mpfr_sqrt(exp, exp, rnd);
	mpfr_exp(exp, exp, rnd);
	mpfr_rint_ceil(exp, exp, rnd);
	mp_exp_t exponent;
	char* str = mpfr_get_str(NULL, &exponent, 10, 0, exp, rnd);
	numberString = string(str);
	mpfr_free_str(str);
	mpfr_clear(exp);
	uint64_t P = stoi(numberString.substr(0,exponent));

	uint64_t A = P*10;
	//step 2
	vector<mpz_class> listingT(A);
	vector<mpz_class> listingTSqr(A);
	mpz_class sqrtN = sqrt(n);
	for (uint64_t t = 0; t < A; ++t) {
		listingT[t] = sqrtN + t + 1;
		listingTSqr[t] = listingT[t]* listingT[t] - n;
	}
	Atkin atkin(P);
	//step 3
	vector<mpz_class> factorBase;
	mpz_class size = atkin.primes.size();
	factorBase.push_back(2);
	for (uint64_t i = 1; i < size; ++i) {//for each odd prime
		if (LegendreSymbol(n, atkin.primes[i]) == 1)
			factorBase.push_back(atkin.primes[i]);
	}
	//step 4 5 6 7 8
	vector<vector<mpz_class>> matrix;
	step45678(n, A, listingT, listingTSqr, factorBase, matrix);

	//step 9
	vector<vector<mpz_class>> factorizationMembers;
	step9(matrix, factorizationMembers);
	//test for factorization
	size = factorizationMembers.size();
	if (size == 0)
		std::cout << "No factorization members" << endl;
	for (uint64_t i = 0; i < size; ++i) {
		vector<mpz_class> factors = factorizationMembers[i];
		vector<mpz_class> powers(matrix.size()-1,0);
		uint64_t j = 0;
		uint64_t pos = 0;
		mpz_class left = 1;
		while (pos < factors.size()) {
			while (matrix[0][j] != factors[pos])
				++j;
			for (uint64_t k = 1; k < matrix.size(); ++k) {
				powers[k - 1] += matrix[k][j];
			}
			left = (left*factors[pos]) % n;
			++pos;
		}
		mpz_class right = 1;
		for (uint64_t k = 0; k < powers.size(); ++k) {
			mpz_class multiplier;
			mpz_pow_ui(multiplier.get_mpz_t(), factorBase[k].get_mpz_t(), mpz_class(powers[k] / 2).get_ui());
			//right = (right*((mpz_class)pow(factorBase[k], powers[k] / 2)) % n) % n;
			right = (right*multiplier) % n;
		}
		std::cout << "left " << left << " right " << right << endl;
		if (left == right) {			
			std::cout << "trivial\n factors: ";
			for (uint64_t k = 0; k < factorizationMembers[i].size(); ++k)
				std::cout << factorizationMembers[i][k] << " ";
			std::cout << endl<<endl;
		}
		else {
			std::cout << "multiplier " << gcd(abs(left - right), n) << endl;
			std::cout << "factors: ";
			for (uint64_t k = 0; k < factorizationMembers[i].size(); ++k)
				std::cout << factorizationMembers[i][k] << " ";
			std::cout << endl << endl;
		}
	}

	return 1;
}


int main(mpz_class argc, char* argv[])
{
	QS(1042385);
	std::cout << endl;

}

void SQRTMODCHECK() {
	Atkin atkin(1000);
	for (int i = 2; i < 20; ++i) {
		mpz_class prime = atkin.primes[i];
		for (mpz_class j = 2; j < prime; ++j) {
			if (LegendreSymbol(j, prime) == 1) {
				mpz_class res = SqrtMod(j, prime);
				std::cout << prime << " " << j << " " << res << " ";
				if ((res*res) % prime == j)
					std::cout << "ok";
				else
					std::cout << "err";
				std::cout << endl;
			}
		}
	}
}
