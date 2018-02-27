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
#include <numeric>
#include <thread>
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

mpz_class sqrtModPowerLiftingPlusOne(mpz_class root, mpz_class quadResidue, mpz_class mod, uint64_t modBase) {
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

uint64_t Step4(mpz_class n, uint64_t p, mpz_class A, mpz_class& res1, mpz_class& res2) {
	mpz_class lowerBound = sqrt(n) + 1;
	mpz_class upperBound = sqrt(n) + A;

	mpz_class residue = n;
	uint64_t modBase = p;
	mpz_class currentMod = p;

	mpz_class root = SqrtMod(residue, currentMod).get_ui();
	mpz_class root2 = -(root - modBase);
	res1 = root;
	res2 = root2;
	uint64_t b = 1;
	while (!chekNuberInRange(lowerBound, upperBound, root, currentMod,A) 
		&& !chekNuberInRange(lowerBound, upperBound, root2, currentMod, A) 
		&& (res1 <= upperBound || res2 <= upperBound)) {
		root = sqrtModPowerLiftingPlusOne(mpz_class(root), residue, currentMod, modBase).get_ui();
		currentMod *= modBase;
		root2 = -(root - currentMod.get_ui());
		++b;
		res1 = root;
		res2 = root2;
	}
	if (res1 > upperBound && res2 > upperBound)
		return 0;
	do {
		res1 = root;
		res2 = root2;
		root = sqrtModPowerLiftingPlusOne(mpz_class(root), residue, currentMod, modBase).get_ui();
		currentMod *= modBase;
		root2 = -(root - currentMod.get_ui());
		++b;
	} while (chekNuberInRange(lowerBound, upperBound, root, currentMod, A)
		|| chekNuberInRange(lowerBound, upperBound, root2, currentMod, A));
	return b - 1;
}

mpz_class liftRootEvenModulo(mpz_class root, mpz_class a, mpz_class p) {
	mpz_class i = ((root * root - a) / p) % 2;
	return root + i * (p / 2);
}

void step7(mpz_class n, mpz_class A, vector<mpz_class>& listingT, vector<uint64_t>& listingTSqr, vector<vector<uint64_t>>& exponentMatrix) {
	uint64_t sizeT = listingT.size();
	if (n % 8 != 1) {
		for (uint64_t j = 0; j < sizeT; ++j) {
			if (listingT[j] % 2 == 1) {
				exponentMatrix[j][0] = 1;
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
		mpz_class mod;
		mpz_class differ;
		for (uint64_t j = 0; j < sizeT; ++j) {
			for (uint64_t k = 0; k < b; ++k) {
				differ = abs(listingT[j] - t1);
				mpz_pow_ui(mod.get_mpz_t(), mpz_class(2).get_mpz_t(), k + 1);
				if (differ % mod == 0) {
					if (exponentMatrix[j][0] < k + 1) {
						exponentMatrix[j][0] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t2);
				if (differ % mod == 0) {
					if (exponentMatrix[j][0] < k + 1) {
						exponentMatrix[j][0] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t3);
				if (differ % mod == 0) {
					if (exponentMatrix[j][0] < k + 1) {
						exponentMatrix[j][0] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t4);
				if (differ % mod == 0) {
					if (exponentMatrix[j][0] < k + 1) {
						exponentMatrix[j][0] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
			}
		}
	}
}

void step45678Thread(uint64_t start, uint64_t stop, mpz_class n, mpz_class A, vector<mpz_class>& listingT, vector<uint64_t>& listingTSqr, vector<uint64_t>& factorBase, vector<vector<uint64_t>>& exponentMatrix) {
	uint64_t sizeT = listingT.size();
	for (uint64_t i = start; i < stop; ++i) {//for each odd prime
		mpz_class t1, t2;
		uint64_t p = factorBase[i];
		uint64_t b = Step4(n, p, A, t1, t2);
		//step 5
		for (uint64_t j = 0; j < sizeT; ++j) {
			mpz_class mod;
			mpz_class differ;
			for (uint64_t k = 0; k < b; ++k) {
				differ = abs(listingT[j] - t1);
				mpz_pow_ui(mod.get_mpz_t(), mpz_class(p).get_mpz_t(), k + 1);
				//if (differ % (mpz_class)pow(p, k + 1) == 0) {
				if (differ % mod == 0) {
					if (exponentMatrix[j][i] < k + 1) {
						exponentMatrix[j][i] = k + 1;
						//step 6
						listingTSqr[j] /= p;
					}
				}
				differ = abs(listingT[j] - t2);
				//if (differ % (mpz_class)pow(p, k + 1) == 0) {
				if (differ % mod == 0) {
					if (exponentMatrix[j][i] < k + 1) {
						exponentMatrix[j][i] = k + 1;
						//step 6
						listingTSqr[j] /= p;
					}
				}
			}
		}
	}
}

void step45678(mpz_class n, mpz_class A, vector<mpz_class>& listingT, vector<uint64_t>& listingTSqr, vector<uint64_t>& factorBase, vector<vector<uint64_t>>& matrix) {
	//step 4
	cout << "step 4\n";
	uint64_t factorBaseSize = factorBase.size();
	uint64_t sizeT = listingT.size();
	vector<vector<uint64_t>> exponentMatrix(sizeT, vector<uint64_t>(factorBaseSize, 0));
	cout << "size " << factorBaseSize << endl;
	cout << "sizeT " << sizeT << endl;
	//create threads
	vector<std::thread> threads;
	uint64_t pos = factorBaseSize / 4;
	threads.push_back(thread(step45678Thread, 1, 1 + pos, n, A, ref(listingT), ref(listingTSqr), ref(factorBase), ref(exponentMatrix)));
	threads.push_back(thread(step45678Thread, 1 + pos, 1 + 2*pos, n, A, ref(listingT), ref(listingTSqr), ref(factorBase), ref(exponentMatrix)));
	threads.push_back(thread(step45678Thread, 1 + 2 * pos, 1 + 3*pos, n, A, ref(listingT), ref(listingTSqr), ref(factorBase), ref(exponentMatrix)));
	threads.push_back(thread(step45678Thread, 1 + 3 * pos, min(1 + 4*pos, factorBaseSize), n, A, ref(listingT), ref(listingTSqr), ref(factorBase), ref(exponentMatrix)));
	//wait for them to complete
	for (auto& th : threads)
		th.join();
	//step 7
	cout << "step 7\n";
	step7(n, A, listingT, listingTSqr, exponentMatrix);
	//step 8
	vector<mpz_class> listingTInMatrix;
	cout << "step 8\n";
	uint64_t matrixSize = 0;
	for (uint64_t i = 0; i < sizeT; ++i) {
		if (listingTSqr[i] == 1)
			++matrixSize;
		listingTInMatrix.push_back(listingT[i]);
	}
	matrix = vector<vector<uint64_t>>(matrixSize, vector<uint64_t>(factorBaseSize, 0));
	pos = 0;
	for (uint64_t i = 0; i < sizeT; ++i) {
		if (listingTSqr[i] == 1) {
			for (uint64_t j = 0; j < factorBaseSize; ++j) {
				matrix[pos][j] = exponentMatrix[i][j];
			}
			++pos;
		}
	}
	listingT = listingTInMatrix;
}

mpz_class factorsCheck(uint64_t rowNumber, mpz_class& n, mpz_class t, vector<uint64_t>& factorBase, vector<uint8_t>& factorsVector,
	vector<vector<uint64_t>>& matrix, vector<mpz_class>& listingTIdentity) {
	
	uint64_t j = 0;
	uint64_t pos = 0;
	mpz_class left = t;
	uint64_t size = factorsVector.size();
	vector<uint64_t> powers(size, 0);
	for (uint64_t i = 0; i < size; ++i) {
		if (factorsVector[i] == 1) {
			std::transform(powers.begin(), powers.end(), matrix[i].begin(), powers.begin(), std::plus<uint64_t>());
			left = left * listingTIdentity[i];
		}
	}
	std::transform(powers.begin(), powers.end(), matrix[rowNumber].begin(), powers.begin(), std::plus<uint64_t>());

	mpz_class right = 1;
	for (uint64_t k = 0; k < powers.size(); ++k) {
		mpz_class multiplier;
		mpz_pow_ui(multiplier.get_mpz_t(), mpz_class(factorBase[k]).get_mpz_t(), mpz_class(powers[k] / 2).get_ui());
		//right = (right*((mpz_class)pow(factorBase[k], powers[k] / 2)) % n) % n;
		right = (right*multiplier) % n;
	}
	if (left != right) {
		mpz_class multiplier = gcd(abs(left - right), n);
		if (abs(multiplier) == 1 || n % multiplier != 0)
			return 0;
		std::cout << "multiplier " << multiplier << endl;
		std::cout << endl << endl;
		return multiplier;
	}
	return 0;
}

mpz_class step9(mpz_class& n, vector<uint64_t>& factorBase, vector<vector<uint64_t>>& matrix, vector<mpz_class>& listingT) {
	//step 9
	uint64_t cols = matrix[0].size();
	uint64_t rows = matrix.size();
	vector<vector<uint8_t>> gaussianEliminationMatrix(cols, vector<uint8_t>(rows, 0));
	for (uint64_t i = 0; i < rows; ++i) {
		for (uint64_t j = 0; j < cols; ++j) {
			gaussianEliminationMatrix[j][i] = matrix[i][j] % 2;
		}
	}
	//Gaussian elimination	
	vector<mpz_class> listingTIdentity(cols,0);
	uint64_t curpos = 0;
	uint64_t lastRow = 0;
	for (uint64_t i = 0; (i < cols) && (curpos < rows); ++curpos) {//possible that first is all zeros?		
		uint64_t onePos = i;
		lastRow = i;
		while (onePos < cols &&
			gaussianEliminationMatrix[onePos][curpos] != 1)
			++onePos;
		if (onePos == cols) {
			continue;
		}		
		if (onePos != i) {
			gaussianEliminationMatrix[onePos].swap(gaussianEliminationMatrix[i]);			
		}
		cout << i << endl;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				cout << (int)gaussianEliminationMatrix[j][i] << " ";
			}
			cout << endl;
		}
		for (uint64_t j = 0; j < cols; ++j) {
			if (j == i
				|| gaussianEliminationMatrix[j][curpos] == 0)
				continue;
			for (uint64_t k = 0; k < rows; ++k) {
				gaussianEliminationMatrix[j][k] ^= gaussianEliminationMatrix[i][k];
			}
		}	
		cout << "sustract" << endl;
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < cols; ++j) {
				cout << (int)gaussianEliminationMatrix[j][i] << " ";
			}
			cout << endl;
		}
		++i;
	}
	vector<uint8_t> factorIndex(cols);
	//get list of combinations	
	for (uint64_t i = lastRow+2; i < rows; ++i) {		
		for (int j = 0; j < cols; ++j)
			factorIndex[j] = gaussianEliminationMatrix[j][i];
		mpz_class result = factorsCheck(i, n, listingT[i], factorBase, factorIndex, matrix, listingTIdentity);
		if (result>0)
			return result;
	}
	return 0;
}

mpz_class QS(mpz_class n, Atkin& atkin) {//n is odd
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
	P= 50;
	A = 500;
	cout << "A = " << A << endl << "P = " << P << endl;
	//step 2
	cout << "step 2\n";
	vector<mpz_class> listingT(A);
	vector<uint64_t> listingTSqr(A);
	mpz_class sqrtN = sqrt(n);
	for (uint64_t t = 0; t < A; ++t) {
		listingT[t] = sqrtN + t + 1;
		listingTSqr[t] = mpz_class(listingT[t]* listingT[t] - n).get_ui();
	}
	if (atkin.primes.back() < P)
		atkin = Atkin(P);
	//step 3
	cout << "step 3\n";
	vector<uint64_t> factorBase;
	uint64_t size = atkin.primes.size();
	factorBase.push_back(2);
	for (uint64_t i = 1; i < size; ++i) {//for each odd prime
		if (atkin.primes[i] > P)
			break;
		if (LegendreSymbol(n, atkin.primes[i]) == 1)
			factorBase.push_back(atkin.primes[i]);
	}
	//step 4 5 6 7 8
	vector<vector<uint64_t>> matrix;
	step45678(n, A, listingT, listingTSqr, factorBase, matrix);

	//step 9
	cout << "step 9\n";
	return step9(n, factorBase, matrix, listingT);
}

int primecheck(mpz_class n) {//0-composite, 1-prime
	if (!mil_rab(n))
		return 0;
	else
		if (!Miller(n))
		{
			return 0;
		}
		else
			return 1;
}

uint64_t dumbcheck(mpz_class num, Atkin& atkin) {
	int size = atkin.primes.size();
	for (int i=0;i<size;++i)
		if (num % atkin.primes[i] == 0)
			return atkin.primes[i];
	return 0;
}

int main(int argc, char* argv[])
{
	Atkin atkin(1000);
	string n = "1042387";
	mpz_class num;
	mpz_set_str(num.get_mpz_t(), n.c_str(), 10);
	vector<mpz_class> factors;
	vector<mpz_class> numbers = { num };
	while (numbers.size() > 0) {
		num = numbers.back();
		cout << "factorizing " << num << endl;
		numbers.pop_back();
		if (primecheck(num)) {
			factors.push_back(num);
			continue;
		}	
		mpz_class factor;
		//actor = dumbcheck(num, atkin);
		//f (factor > 0) {
		//	factors.push_back(factor);
		//	numbers.push_back(num / factor);
		//	continue;
		//
		factor = QS(num, atkin);
		if (factor == 0) {
			cout << "factorization error\n";
			return 1;
		}
		numbers.push_back(factor);
		numbers.push_back(num/factor);
		
	}	
	cout << "factors\n";
	for (auto it = factors.begin(); it != factors.end(); ++it)
		cout << *it << endl;
	std::cout << endl;
	return 0;
}

