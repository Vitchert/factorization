#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include "BigInt.h"
#include "AtkinSieve.cpp"
#include <map>
#include <set>
using namespace std;



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





int powmod(int base, int power, int mod)
{
	int res = 1, a = base, k = power, n = mod;

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

int LegendreSymbol(int a, int p) {
	if (a >= p)
		a %= p;
	int res =  powmod(a, (p - 1) / 2, p);
	return res > 1 ? -1 : res;
}

int SqrtMod(int a, int p) {
	if (a >= p)
		a %= p;
	int n = 1;
	while (LegendreSymbol(n, p) >= 0) {
		++n;
	}
	int alpha = 0;
	int s = p - 1;
	while (s % 2 == 0) {
		++alpha;
		s /= 2;
	}
	int b = powmod(n, s, p);
	int r = powmod(a, (s + 1) / 2, p);
	int rCalc = powmod(a, s, p);
	int j;
	//int check = powmod((r*r) / a, pow(2, alpha - 2), p);
	int check = powmod(rCalc, pow(2, alpha - 2), p);
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
	for (int i = 1; i < alpha - 1; ++i) {
		//check = powmod(pow((pow(b, j)*r), 2) / a, pow(2, alpha - i - 2), p);
		check = powmod(powmod(powmod(b, j, p), 2, p) * rCalc, pow(2, alpha - i - 2), p);
		if (check > 1)
			check -= p;
		if (check == 1) {
			//j = 0;
		}
		else
			if (check == -1) {
				j += pow(2, i);
			}
			else
				cout << "check error" << endl;
	}
	return (powmod(b, j, p) * r) % p;
}


int gcd(int a, int b) {
	while (b) {
		a %= b;
		swap(a, b);
	}
	return a;
}

int gcdex(int a, int b, int & x, int & y) {//–асширенный алгоритм евклида emaxx
	if (a == 0) {
		x = 0; y = 1;
		return b;
	}
	int x1, y1;
	int d = gcdex(b%a, a, x1, y1);
	x = y1 - (b / a) * x1;
	y = x1;
	return d;
}

int ReverseMod(int a, int m) {//обратное по модулю число
	int x = 0, y;
	int g = gcdex(a, m, x, y);
	if (g != 1)
		cout << "ReverseMod no solution";
	else {
		x = (x % m + m) % m;
	}
	return x;
}

int sqrtModPowerLiftingPlusOne(int root, int quadResidue, int mod, int modBase) {
	int newmod = mod*modBase;
	int yk = ReverseMod(2 * root, mod)% newmod;	
	quadResidue %= newmod;
	root = root - (((root*root) % newmod - quadResidue)*yk) % newmod;

	if (root < 0) {
		root = -root;
	}
	root %= newmod;
	//root = min(root, -(root - newmod));
	return root;
}

int primeCheck() {
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

bool chekNuberInRange(int lowerBound, int upperBound, int t, int p, int A) {// check if there is range lower <=T<= upper that T = t1 or t2 ( mod p)
	if (A >= p) {
		return true;
	}
	int l = lowerBound % p;
	int u = upperBound % p;
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

int Step4(int n, int p, int A, int& res1, int& res2) {
	int lowerBound = sqrt(n) + 1;
	int upperBound = sqrt(n) + A;

	int residue = n;
	int modBase = p;
	int currentMod = p;

	int root = SqrtMod(residue, currentMod);
	int root2 = -(root - modBase);
	res1 = root;
	res2 = root2;
	int b = 1;
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

int liftRootEvenModulo(int root, int a, int p) {
	int i = ((root * root - a) / p) % 2;
	return root + i * (p / 2);
}

vector<double> gauss(vector< vector<double> > A) {
	int n = A.size();

	for (int i = 0; i<n; i++) {
		// Search for maximum in this column
		double maxEl = abs(A[i][i]);
		int maxRow = i;
		for (int k = i + 1; k<n; k++) {
			if (abs(A[k][i]) > maxEl) {
				maxEl = abs(A[k][i]);
				maxRow = k;
			}
		}

		// Swap maximum row with current row (column by column)
		for (int k = i; k<n + 1; k++) {
			double tmp = A[maxRow][k];
			A[maxRow][k] = A[i][k];
			A[i][k] = tmp;
		}

		// Make all rows below this one 0 in current column
		for (int k = i + 1; k<n; k++) {
			double c = -A[k][i] / A[i][i];
			for (int j = i; j<n + 1; j++) {
				if (i == j) {
					A[k][j] = 0;
				}
				else {
					A[k][j] += c * A[i][j];
				}
			}
		}
	}

	// Solve equation Ax=b for an upper triangular matrix A
	vector<double> x(n);
	for (int i = n - 1; i >= 0; i--) {
		x[i] = A[i][n] / A[i][i];
		for (int k = i - 1; k >= 0; k--) {
			A[k][n] -= A[k][i] * x[i];
		}
	}
	return x;
}

void step7(int n, int A, vector<int>& listingT, vector<int>& listingTSqr, vector<vector<int>>& exponentMatrix) {
	int sizeT = listingT.size();
	if (n % 8 != 1) {
		for (int j = 0; j < sizeT; ++j) {
			if (listingT[j] % 2 == 1) {
				exponentMatrix[0][j] = 1;
				listingTSqr[j] /= 2;
			}
		}
	}
	else {
		int lowerBound = sqrt(n) + 1;
		int upperBound = sqrt(n) + A;

		int residue = n;
		int modBase = 2;
		int currentMod = 2;

		int b = 3;

		int root1 = 1;
		int root2 = 3;
		int root3 = 5;
		int root4 = 7;
		int t1 = 1;
		int t2 = 3;
		int t3 = 5;
		int t4 = 7;
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

		for (int j = 0; j < sizeT; ++j) {
			for (int k = 0; k < b; ++k) {
				int differ = abs(listingT[j] - t1);
				if (differ % (int)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t2);
				if (differ % (int)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t3);
				if (differ % (int)pow(2, k + 1) == 0) {
					if (exponentMatrix[0][j] < k + 1) {
						exponentMatrix[0][j] = k + 1;
						//step 6
						listingTSqr[j] /= 2;
					}
				}
				differ = abs(listingT[j] - t4);
				if (differ % (int)pow(2, k + 1) == 0) {
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
int QS(int n) {//n is odd
	//step 1
	int P = exp(sqrt(log(n)*log(log(n))));
	int A = P*10;
	//int P = 50;
	//int A = 500;
	//step 2
	vector<int> listingT(A);
	vector<int> listingTSqr(A);
	int sqrtN = sqrt(n);
	for (int t = 0; t < A; ++t) {
		listingT[t] = sqrtN + t + 1;
		listingTSqr[t] = listingT[t]* listingT[t] - n;
	}
	Atkin atkin(P);
	//step 3
	vector<int> factorBase;
	int size = atkin.primes.size();
	factorBase.push_back(2);
	for (int i = 1; i < size; ++i) {//for each odd prime
		if (LegendreSymbol(n, atkin.primes[i]) == 1)
			factorBase.push_back(atkin.primes[i]);
	}
	//step 4
	size = factorBase.size();
	int sizeT = listingT.size();
	vector<vector<int>> exponentMatrix(size, vector<int>(sizeT,0));
	for (int i = 1; i < size; ++i) {//for each odd prime
		int t1, t2;
		int p = factorBase[i];
		int b = Step4(n, p, A, t1, t2);
		//step 5
		for (int j = 0; j < sizeT; ++j) {
			for (int k = 0; k < b; ++k) {
				int differ = abs(listingT[j] - t1);
				if (differ % (int)pow(p,k+1) == 0) {
					if (exponentMatrix[i][j] < k + 1) {
						exponentMatrix[i][j] = k + 1;
						//step 6
						listingTSqr[j] /= p;
					}
				}
				differ = abs(listingT[j] - t2);
				if (differ % (int)pow(p, k + 1) == 0) {
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
	int matrixSize = 0;
	for (int i = 0; i < sizeT; ++i) {
		if (listingTSqr[i] == 1)
			++matrixSize;
	}
	vector<vector<int>> matrix(size + 1, vector<int>(matrixSize, 0));
	int pos = 0;
	for (int i = 0; i < sizeT; ++i) {
		if (listingTSqr[i] == 1) {
			matrix[0][pos] = listingT[i];
			for (int j = 0; j < size; ++j) {
				matrix[j+1][pos] = exponentMatrix[j][i];
			}
			++pos;
		}
	}
	//step 9
	int cols = matrix.size();
	int rows = matrix[0].size();
	vector<vector<int>> gaussianEliminationMatrix(rows, vector<int>(cols, 0));
	for (int i = 0; i < rows; ++i) {
		gaussianEliminationMatrix[i][0] = matrix[0][i];
		for (int j = 1; j < cols; ++j) {
			gaussianEliminationMatrix[i][j] = matrix[j][i]%2;
		}
	}
	//Gaussian elimination

	map<int,multiset<int>> linearCompositionMembers;
	for (int i = 0; i < rows; ++i) {
		linearCompositionMembers.insert({ gaussianEliminationMatrix[i][0], {gaussianEliminationMatrix[i][0]} });
	}

	int curpos = 1;
	for (int i = 0; (i < rows) && (curpos < cols); ++i,++curpos) {
		int onePos = i;
		while (onePos < rows && 
			gaussianEliminationMatrix[onePos][curpos]!= 1)
			++onePos;
		if (onePos == rows) {
			if (i)
				--i;
			continue;
		}
		if (onePos != i) {
			swap(gaussianEliminationMatrix[onePos], gaussianEliminationMatrix[i]);
		}
		
		for (int j = 0; j < rows; ++j) {
			if (j == i
				|| gaussianEliminationMatrix[j][curpos] ==0)
				continue;
			for (int k = 1; k < cols; ++k) {
				gaussianEliminationMatrix[j][k] = (gaussianEliminationMatrix[j][k] + gaussianEliminationMatrix[i][k]) % 2;
			}
			if(j>i)
				linearCompositionMembers[gaussianEliminationMatrix[j][0]].insert(
					linearCompositionMembers[gaussianEliminationMatrix[i][0]].begin(),
					linearCompositionMembers[gaussianEliminationMatrix[i][0]].end());

		}
	}
	//get list of posible combinations
	vector<multiset<int>> posibleFactorizationMembers;

	for (int i = 0; i < rows; ++i) {
		int j = 1;
		while (j < cols && gaussianEliminationMatrix[i][j] == 0)
			++j;
		if (j < cols)
			continue;
		posibleFactorizationMembers.push_back(linearCompositionMembers[gaussianEliminationMatrix[i][0]]);
	}
	//remove all even encounters
	vector<vector<int>> factorizationMembers;
	
	size = posibleFactorizationMembers.size();
	for (int i = 0; i < size; ++i) {
		vector<int> factors;
		while (posibleFactorizationMembers[i].size() > 0) {
			multiset<int>::iterator it = posibleFactorizationMembers[i].begin();
			if (posibleFactorizationMembers[i].count(*it) % 2 == 1)
				factors.push_back(*it);
			posibleFactorizationMembers[i].erase(*it);
		}
		factorizationMembers.push_back(factors);
	}
	//test for factorization
	size = factorizationMembers.size();
	if (size == 0)
		cout << "No factorization members" << endl;
	for (int i = 0; i < size; ++i) {
		vector<int> factors = factorizationMembers[i];
		vector<int> powers(matrix.size()-1,0);
		int j = 0;
		pos = 0;
		int left = 1;
		while (pos < factors.size()) {
			while (matrix[0][j] != factors[pos])
				++j;
			for (int k = 1; k < matrix.size(); ++k) {
				powers[k - 1] += matrix[k][j];
			}
			left = (left*factors[pos]) % n;
			++pos;
		}
		int right = 1;
		for (int k = 0; k < powers.size(); ++k) {
			right = (right*((int)pow(factorBase[k], powers[k] / 2)) % n) % n;
		}
		cout << "left " << left << " right " << right << endl;
		if (left == right) {			
			cout << "trivial\n factors: ";
			for (int k = 0; k < factorizationMembers[i].size(); ++k)
				cout << factorizationMembers[i][k] << " ";
			cout << endl<<endl;
		}
		else {
			cout << "multiplier " << gcd(left - right, n) << endl;
			cout << "factors: ";
			for (int k = 0; k < factorizationMembers[i].size(); ++k)
				cout << factorizationMembers[i][k] << " ";
			cout << endl << endl;
		}
	}

	return 1;
}


int main(int argc, char* argv[])
{
	int a = -1;
	int b = a % 2;
	//int r = 1021 % 729;
	//int b = 1520 % 729;
	//bool q = chekNuberInRange(1021, 1520, 589, 2187, 500);
	//Atkin atkin(1000);
	//QS(1042387);
	QS(1042387);
	cout << endl;

}

void SQRTMODCHECK() {
	Atkin atkin(1000);
	for (int i = 2; i < 20; ++i) {
		int prime = atkin.primes[i];
		for (int j = 2; j < prime; ++j) {
			if (LegendreSymbol(j, prime) == 1) {
				int res = SqrtMod(j, prime);
				cout << prime << " " << j << " " << res << " ";
				if ((res*res) % prime == j)
					cout << "ok";
				else
					cout << "err";
				cout << endl;
			}
		}
	}
}
