#include <iostream>
#include <algorithm>
#include <iomanip>
using namespace std;

// a struct used by the 3 operations
struct Union
{
	int expRow;
	int coefColumn;
	int value;
};

// a class for the "merge" operation
class list
{
public:
	void merge(Union*, Union*, Union*, int, int, string);
	bool Compare(Union, Union, string );
	void NewTerm(Union* MergeArr, Union* tobeMerged, int pMrg, int ptbM)
	{
		MergeArr[pMrg] = tobeMerged[ptbM];
	}

};

// compare element value for merging
bool list::Compare(Union lefti, Union righti, string mode)
{
	if (mode == "MergeSort")
	{
		if (lefti.value > righti.value)
			return true;
		return false;
	}
	else if (mode == "SPmatrix")
	{
		if (lefti.coefColumn > righti.coefColumn)
			return true;
		return false;
	}
	else  // (mode == "Poly")
	{
		if (lefti.expRow > righti.expRow)
			return true;
		return false;
	}
}


// Merge operation: from small to large
void list::merge(Union* Left, Union* Right, Union* Merge, 
			int leftBoundForMerge, int rightBoundForMerge, string mode)   
{

	int pL = 0, pR = 0, pMrg = 0;
	string check;

	// compare and move
	while (pL < leftBoundForMerge && pR < rightBoundForMerge)
	{
		if (Compare(Left[pL], Right[pR], mode))   // Lefti > Righti
		{
			NewTerm(Merge, Right, pMrg, pR);
			pR++;
		}
		else   // Lefti < or == Righti
		{
			NewTerm(Merge, Left, pMrg, pL);
			pL++;
		}

		pMrg++;
	}

	//move the rest of Array
	for (; pR < rightBoundForMerge; pR++)
	{
		NewTerm(Merge, Right, pMrg, pR);	
		pMrg++;
	}
	for (; pL < leftBoundForMerge; pL++)
	{
		NewTerm(Merge, Left, pMrg, pL);
		pMrg++;
	}

}



class MergeSort
{
public:
	void sort(Union*, Union*, int);
};

// Merge Sort: Bottom-up
void MergeSort::sort(Union* OriginalArr, Union* MergeArr, int length) // OriginalArr: original array 
{																      // MergeArr: array after merging
	Union* Left = new Union[length]; 
	Union* Right = new Union[length];
	int i, count, leftStart, rightStart, leftBound, rightBound, rightBoundForMerge, leftBoundForMerge;
	list var;

	for (count = 1; count < length; count *= 2) //count = numbers taken Every time
	{
		for (leftStart = 0; leftStart < length; leftStart += 2 * count) // exactly decide which elements will be taken
		{
			// divide OriginalArr into many subarrays(called "Left" and "Right"). 
			leftBound = leftStart + count;
			rightStart = leftStart + count;
			rightBound = rightStart + min(count, length - rightStart);
			if (rightBound <= rightStart)
				break;

			for (i = leftStart; i < leftBound; i++)   //decide Left array
				Left[i - leftStart] = OriginalArr[i];

			for (i = rightStart; i < rightBound; i++)   //decide Right array
				Right[i - rightStart] = OriginalArr[i];

			leftBoundForMerge = leftBound - leftStart;
			rightBoundForMerge = rightBound - rightStart;

			//Merge
			var.merge(Left, Right, MergeArr, leftBoundForMerge, rightBoundForMerge, "MergeSort");	

			for (i = leftStart; i < rightBound; i++)  // 
				OriginalArr[i] = MergeArr[i - leftStart];
		}
	}
	delete[] Left;
	delete[] Right;

}



class spmatrix
{
public:
	Union* mul(Union*, Union*, int, int, int&, int);

};

// Sparse Matrix Multiplication 
// mul = spLeft * spRight
Union* spmatrix::mul(Union* spLeft, Union* spRight,int NumSpLeft, int NumSpRight, int &pMnum, int totalNum)
{
	Union *xchgRight = new Union[NumSpRight];
	int i;

	// row and column exchange in Right
	for (i = 0; i < NumSpRight; i++)
	{
		xchgRight[i].expRow = spRight[i].coefColumn;
		xchgRight[i].coefColumn = spRight[i].expRow;
		xchgRight[i].value = spRight[i].value;
	};


	Union *mul = new Union[totalNum], *secMerge = new Union[NumSpLeft + NumSpRight];
	Union *secLeft = new Union[NumSpLeft], *secXchgRight = new Union[NumSpRight];
	int pL = 0, pR = 0, pMul = 0, leftRowIdx, rightClmnIdx, sum = 0;
	int pLStart, pRStart, mgLength, pMrg;
	list var;

	// Left Matrix
	while (pL < NumSpLeft)
	{
		leftRowIdx = spLeft[pL].expRow;
		pLStart = pL;

		// decide where to stop for the pointer of "spLeft"
		while (leftRowIdx == spLeft[pL].expRow && pL < NumSpLeft)
		{
			secLeft[pL - pLStart] = spLeft[pL];
			pL++;
		}

		// Right Matrix
		while (pR < NumSpRight)
		{
			rightClmnIdx = xchgRight[pR].expRow;
			pRStart = pR;

			// assign elements for merging, and decide where to stop for spRight pointer
			while (rightClmnIdx == xchgRight[pR].expRow && pR < NumSpRight) 
			{
				secXchgRight[pR - pRStart] = xchgRight[pR];
				pR++;
			} 

			var.merge(secLeft, secXchgRight, secMerge, (pL - pLStart), (pR - pRStart), "SPmatrix");  // Merge

			mgLength = (pL - pLStart) + (pR - pRStart); 
			pMrg = 0;

			// multiply the value whose index are the same
			while (pMrg < (mgLength-1) )  
				if (secMerge[pMrg].coefColumn == secMerge[pMrg + 1].coefColumn)  
				{
					sum = sum + secMerge[pMrg].value * secMerge[pMrg + 1].value;
					pMrg = pMrg + 2;
				}
				else
					pMrg++;

			if (sum != 0)  // sum = final element value
			{
				mul[pMul].expRow = leftRowIdx;
				mul[pMul].coefColumn = rightClmnIdx;
				mul[pMul].value = sum;
				sum = 0;  // reset sum
				pMul++;
			}
		} // end while (b < m) 
		while (pL < NumSpLeft && spLeft[pL].expRow == leftRowIdx) // move spLeft pointer to the next start
			pL++;

		pR = 0;  // move spRight pointer to the start of array 
	}
	pMnum = pMul;
	delete[] secMerge, secLeft, secXchgRight, xchgRight;

	return mul;
} // end spmatrix::mul()


class poly
{
public:
	Union* add(Union*, Union*, int, int, int&);

};

// Polynomial Addition
// Sum = Left + Right
Union* poly::add(Union* Left, Union* Right, int n, int m, int &pS)
{

	int pSum = 0, i, pMrg = 0;
	int sumCoef, MrgLength = n + m;
	Union *Sum = new Union[n + m], *tmpMrg = new Union[n + m];
	Union *reverseLeft = new Union[n], *reverseRight = new Union[m];
	list var;

	// reverse array -> from small to Large
	for (i = 0; i < n; i++)	
		reverseLeft[i] = Left[n - i - 1];
	for (i = 0; i < m; i++)  
		reverseRight[i] = Right[m - i - 1];

	var.merge(reverseLeft, reverseRight, tmpMrg, n, m, "Poly");    //Merge

	// caculate the sum, if exp are the same
	while (pMrg < MrgLength)
		if (tmpMrg[pMrg].expRow == tmpMrg[pMrg + 1].expRow)
		{
			sumCoef = tmpMrg[pMrg].coefColumn + tmpMrg[pMrg + 1].coefColumn;
			if (sumCoef)
			{
				Sum[pSum].coefColumn = sumCoef;
				Sum[pSum].expRow = tmpMrg[pMrg].expRow;
				pSum++;
			}
			pMrg = pMrg + 2;
		}
		else
		{
			Sum[pSum] = tmpMrg[pMrg];
			pMrg++;
			pSum++;
		}
	pS = pSum;
	delete[] tmpMrg, reverseLeft, reverseRight;

	return Sum;
}



int main()
{
	int i, j;

	//======= Merge Sort =======

	// Merge-- exmaple 1:
	int lenghtLeft1 = 4, lengthRight1 = 4, intSubLeft1[4] = { 2,4,5,7 }, intSubRight1[4] = { 1,2,3,6 };
	Union subLeft1[4], subRight1[4], Merge1[8];
	list merge1;

	cout << "Merge-- exmaple 1:" << endl;
	cout << "A: ";
	for (i = 0; i < lenghtLeft1; i++)    // show initialization
	{
		subLeft1[i].value = intSubLeft1[i];
		cout << setw(3) << subLeft1[i].value;
	}

	cout << "      B: ";
	for (i = 0; i < lengthRight1; i++)    // show initialization
	{
		subRight1[i].value = intSubRight1[i];
		cout << setw(3) << subRight1[i].value;
	}
	cout << endl << endl;

	merge1.merge(subLeft1, subRight1, Merge1, lenghtLeft1, lengthRight1, "MergeSort"); //Merge

	cout << "Merge A and B:  ";    // show results
	for (i = 0; i < lenghtLeft1 + lengthRight1; i++)
		cout << setw(3) << Merge1[i].value;
	cout << endl << endl << endl;


	// Merge-- exmaple 2:
	int lenghtLeft2 = 5, lengthRight2 = 6, intSubLeft2[5] = { -3,-1,0,4,6 }, intSubRight2[6] = { -6,-1,3,6,8,9 };
	Union subLeft2[5], subRight2[6], Merge2[11];
	list merge2;

	cout << "Merge-- exmaple 2:" << endl;
	cout << "C: ";
	for (i = 0; i < lenghtLeft2; i++)    // show initialization
	{
		subLeft2[i].value = intSubLeft2[i];
		cout << setw(3) << subLeft2[i].value;
	}

	cout << "      D: ";
	for (i = 0; i < lengthRight2; i++)    // show initialization
	{
		subRight2[i].value = intSubRight2[i];
		cout << setw(3) << subRight2[i].value;
	}
	cout << endl << endl;

	merge2.merge(subLeft2, subRight2, Merge2, lenghtLeft2, lengthRight2, "MergeSort"); //Merge

	cout << "Merge C and D:  ";    // show results
	for (i = 0; i < lenghtLeft2 + lengthRight2; i++)
		cout << setw(3) << Merge2[i].value;
	cout << endl << endl << endl;


	// MergeSort-- exmaple 1:
	int lenght1 = 8, intOriginalArr1[8] = { 2,4,5,7,1,2,3,6 };
	Union MergeArr1[8], OriginalArr1[8];
	MergeSort MergeSort1;

	cout << "MergeSort-- exmaple 1:" << endl;
	cout << "Before: ";
	for (i = 0; i < lenght1; i++)
		cout << setw(3) << intOriginalArr1[i];
	cout << endl;

	for (i = 0; i < lenght1; i++)
		OriginalArr1[i].value = intOriginalArr1[i];

	MergeSort1.sort(OriginalArr1, MergeArr1, lenght1); //Merge
	
	cout << "After:  ";    //output
	for (i = 0; i < lenght1; i++)
		cout << setw(3) << MergeArr1[i].value;
	cout << endl << endl;

	// MergeSort-- exmaple 2:
	int lenght2 = 12, intOriginalArr2[12] = {-25, 598, 0, 55, -3, 25, 63, 0, -78, 10, 55, 258};
	Union MergeArr2[12], OriginalArr2[12];
	MergeSort MergeSort2;

	cout << "MergeSort-- exmaple 2:" << endl;
	cout << "Before: ";
	for (i = 0; i < lenght2; i++)
		cout << setw(4) << intOriginalArr2[i];
	cout << endl;

	for (i = 0; i < lenght2; i++)
		OriginalArr2[i].value = intOriginalArr2[i];

	MergeSort2.sort(OriginalArr2, MergeArr2, lenght2); //Merge
	
	cout << "After:  ";    //output
	for (i = 0; i < lenght2; i++)
		cout << setw(4) << MergeArr2[i].value;
	cout << endl << endl << endl << endl;



	//======= Polynomial addition =======

	poly PolyEx1;
	int polyL1Num = 4, polyR1Num = 4, pS1num;
	cout << "polynomial addition-- example 1:" << endl << endl;

	// initialization
	Union polyL1[8], polyR1[11], *Sum1;
	polyL1[0].coefColumn =  2;    polyL1[0].expRow = 2000;
	polyL1[1].coefColumn =  5;    polyL1[1].expRow = 3;
	polyL1[2].coefColumn = -3;    polyL1[2].expRow = 2;
	polyL1[3].coefColumn =  1;    polyL1[3].expRow = 0;


	cout << "f1(x) = ";    // show initialization
	for (i = 0; i < polyL1Num; i++)
	{
		if (polyL1[i].expRow != 0)
			cout << polyL1[i].coefColumn << "x^" << polyL1[i].expRow;
		else
			cout << polyL1[i].coefColumn;
		if (i < polyL1Num - 1)
			cout << "+ ";
	}
	cout << endl;

	// initialization
	polyR1[0].coefColumn =  1;    polyR1[0].expRow = 4;
	polyR1[1].coefColumn = 10;    polyR1[1].expRow = 3;
	polyR1[2].coefColumn =  3;    polyR1[2].expRow = 2;
	polyR1[3].coefColumn =  1;    polyR1[3].expRow = 0;


	cout << "f2(x) = ";    // show initialization
	for (i = 0; i < polyR1Num; i++)
	{
		if (polyR1[i].expRow != 0)
			cout << polyR1[i].coefColumn << "x^" << polyR1[i].expRow;
		else
			cout << polyR1[i].coefColumn;
		if (i < polyR1Num - 1)
			cout << "+ ";
	}
	cout << endl;

	Sum1 = PolyEx1.add(polyL1, polyR1, polyL1Num, polyR1Num, pS1num);  // poly::add()

	cout << endl << "f1(x) + f2(x) = ";    // show results
	for (i = 0; i < pS1num; i++)
	{
		if (Sum1[pS1num - i - 1].expRow != 0)
			cout << Sum1[pS1num - i - 1].coefColumn << "x^" << Sum1[pS1num - i - 1].expRow;
		else
			cout << Sum1[pS1num - i - 1].coefColumn;
		if (i < pS1num - 1)
			cout << "+ ";
	}
	cout << endl << endl << endl;


	//////////  Polynomial example 2
	poly PolyEx2;
	int polyL2Num = 8, polyR2Num = 11, pS2num;

	cout << "polynomial addition-- example 2:" << endl << endl;

	// initialization
	Union polyL2[8], polyR2[11], *Sum2;
	polyL2[0].coefColumn = 2;   polyL2[0].expRow = 1000;
	polyL2[1].coefColumn = 55;  polyL2[1].expRow = 100;
	polyL2[2].coefColumn = 5;   polyL2[2].expRow = 10;
	polyL2[3].coefColumn = -3;  polyL2[3].expRow = 7;
	polyL2[4].coefColumn = 5;   polyL2[4].expRow = 6;
	polyL2[5].coefColumn = 9;   polyL2[5].expRow = 5;
	polyL2[6].coefColumn = 7;   polyL2[6].expRow = 4;
	polyL2[7].coefColumn = 1;   polyL2[7].expRow = 2;

	cout << "f1(x) = ";    // show initialization
	for (i = 0; i < polyL2Num; i++)
	{
		if(polyL2[i].expRow !=0)
			cout << polyL2[i].coefColumn << "x^" << polyL2[i].expRow;
		else 
			cout << polyL2[i].coefColumn;
		if (i < polyL2Num - 1)
			cout << "+ ";
	}
	cout << endl;

	// initialization
	polyR2[0].coefColumn = 6;    polyR2[0].expRow = 13;
	polyR2[1].coefColumn = 1;    polyR2[1].expRow = 12;
	polyR2[2].coefColumn = -5;   polyR2[2].expRow = 10;
	polyR2[3].coefColumn = 10;   polyR2[3].expRow = 9;
	polyR2[4].coefColumn = 3;    polyR2[4].expRow = 8;
	polyR2[5].coefColumn = -7;   polyR2[5].expRow = 7;
	polyR2[6].coefColumn = 15;   polyR2[6].expRow = 6;
	polyR2[7].coefColumn = -9;   polyR2[7].expRow = 5;
	polyR2[8].coefColumn = 2;    polyR2[8].expRow = 3;
	polyR2[9].coefColumn = 22;   polyR2[9].expRow = 1;
	polyR2[10].coefColumn = -6;  polyR2[10].expRow = 0;

	cout << "f2(x) = ";    // show initialization
	for (i = 0; i < polyR2Num; i++)
	{
		if (polyR2[i].expRow != 0)
			cout << polyR2[i].coefColumn << "x^" << polyR2[i].expRow;
		else
			cout << polyR2[i].coefColumn;
		if (i < polyR2Num - 1)
			cout << "+ ";
	}
	cout << endl;

	Sum2 = PolyEx2.add(polyL2, polyR2, polyL2Num, polyR2Num, pS2num);  // poly::add()

	cout << endl << "f1(x) + f2(x) = ";    // show results
	for (i = 0; i < pS2num; i++)
	{
		if (Sum2[pS2num -i-1].expRow != 0)
			cout << Sum2[pS2num -i-1].coefColumn << "x^" << Sum2[pS2num -i-1].expRow;
		else
			cout << Sum2[pS2num -i-1].coefColumn;
		if (i < pS2num - 1)
			cout << "+ ";
	}
	cout << endl << endl << endl << endl;
	


	//======= Sparse matrix multipication =======
	//Matrix A * B
	cout << "Sparse matrix multipication-- exmaple 1:" << endl << endl;

	spmatrix spEx1;

	int NumSp1 = 6, totalElemtNumSp1 = NumSp1*NumSp1;
	int sp1LeftMatrix[6][6] =     // Matrix A
	{
		{  15,   0,   0,  22,   0, -15},
		{   0,  11,   3,   0,   0,   0},
		{   0,   0,   0,  -6,   0,   0},
		{   0,   0,   0,   0,   0,   0},
		{  91,   0,   0,   0,   0,   0},
		{   0,   0,  28,   0,   0,   0}
	};

	int sp1RightMatrix[6][6] =     // Matrix B
	{
		{   0,  -8,  -1,   0,   0,   0},
		{   0,   0,   2,   0,   0,   0},
		{   5,   0,   2,   4,   0,   0},
		{   0,   4,   0,   0,   0,   0},
		{   0,   0,   5,   0,   0,   8},
		{   0,   1,   0,   0,  -7,   0}
	};

	cout << "Matrix A = " << "                                  " << "Matrix B = " << endl;

	// show Matrix A & B on the screen
	for (i = 0; i < NumSp1; i++)
	{
		for (j = 0; j < NumSp1; j++)
		{
			cout << setw(5) << sp1LeftMatrix[i][j];
		}
		cout << "            ";
		for (j = 0; j < NumSp1; j++)
		{
			cout << setw(5) << sp1RightMatrix[i][j];
		}
		cout << endl;
	}

	
	Union sp1Left[8], sp1Right[11], *mul;
	int pL1 = 0, pR1 = 0, pL1Num, pR1Num, pMnum, pM = 0;

	// non-zero elements taken
	for (i = 0; i < NumSp1; i++)
		for (j = 0; j < NumSp1; j++)
			if (sp1LeftMatrix[i][j])
			{
				sp1Left[pL1].expRow = i;
				sp1Left[pL1].coefColumn = j;
				sp1Left[pL1].value = sp1LeftMatrix[i][j];
				pL1++;
			}

	for (j = 0; j < NumSp1; j++)
		for (i = 0; i < NumSp1; i++)
			if (sp1RightMatrix[i][j])
			{
				sp1Right[pR1].expRow = i;
				sp1Right[pR1].coefColumn = j;
				sp1Right[pR1].value = sp1RightMatrix[i][j];
				pR1++;
			}

	pL1Num = pL1;
	pR1Num = pR1;

	mul = spEx1.mul(sp1Left, sp1Right, pL1Num, pR1Num, pMnum, totalElemtNumSp1); // Multiply

	cout << endl <<endl;

	// show multiplication results on the screen
	cout << "Matrix (AxB) = " << endl;
	for (i = 0; i < NumSp1; i++)
	{
		for (j = 0; j < NumSp1; j++)
		{
			if (pM < pMnum && mul[pM].expRow == i && mul[pM].coefColumn == j)
			{
				cout << setw(5) << mul[pM].value;
				pM++;
			}
			else
				cout << setw(5) << 0;
		}
		cout << endl;
	}
	cout << endl << endl << endl;


	///// Sparse Matrix example 2
	cout << "Sparse matrix multipication-- exmaple 2:" << endl << endl;
	spmatrix spEx2;

	int NumSp2 = 7, totalElemtNumSp2 = NumSp2*NumSp2;
	int sp2LeftMatrix[7][7] =     // Matrix A
	{
		{   2,   0,   8,   0,   0,   0,   0 },
		{   0,  -3,   0,   0,   0,   0,   0 },
		{   0,   0,   0,   0,   0,   0,   2 },
		{   0,   0,   5,   0,   4,   4,   0 },
		{   0,   0,   0,   7,   0,   0,   0 },
		{  -1,   0,   0,  -5,   0,   0,   3 },
		{   0,   0,   1,   0,   2,   0,   0 }
	};

	int sp2RightMatrix[7][7] =     // Matrix B
	{
		{   0,   0,   0,  -5,   0,   0,   7 },
		{   0,   0,   0,   0,   0,  -9,   0 },
		{   7,   0,  18,   0,   0,   0,   0 },
		{   0,   0,   0,   0,   0,   0,   7 },
		{  11,   0,   0,   0,  -6,   0,   0 },
		{   0,  -8,   0,   0,   0,   6,   5 },
		{   0,   0,   0,  17,   0,   0,   0 }
	};

	cout << "Matrix A = " << "                                  " << "Matrix B = " << endl;

	// show Matrix A & B on the screen
	for (i = 0; i < NumSp2; i++)
	{
		for (j = 0; j < NumSp2; j++)
		{
			cout << setw(5) << sp2LeftMatrix[i][j];
		}
		cout << "            ";
		for (j = 0; j < NumSp2; j++)
		{
			cout << setw(5) << sp2RightMatrix[i][j];
		}
		cout << endl;
	}


	Union sp2Left[13], sp2Right[12], *mul2;
	int pL2 = 0, pR2 = 0, pL2Num, pR2Num, pM2num, pM2 = 0;

	// non-zero elements taken
	for (i = 0; i < NumSp2; i++)
		for (j = 0; j < NumSp2; j++)
			if (sp2LeftMatrix[i][j])
			{
				sp2Left[pL2].expRow = i;
				sp2Left[pL2].coefColumn = j;
				sp2Left[pL2].value = sp2LeftMatrix[i][j];
				pL2++;
			}

	for (j = 0; j < NumSp2; j++)
		for (i = 0; i < NumSp2; i++)
			if (sp2RightMatrix[i][j])
			{
				sp2Right[pR2].expRow = i;
				sp2Right[pR2].coefColumn = j;
				sp2Right[pR2].value = sp2RightMatrix[i][j];
				pR2++;
			}

	pL2Num = pL2;
	pR2Num = pR2;

	mul2 = spEx2.mul(sp2Left, sp2Right, pL2Num, pR2Num, pM2num, totalElemtNumSp2);    // Multiply

	cout << endl << endl;

	// show multiplication results on the screen
	cout << "Matrix (AxB) = " << endl;
	for (i = 0; i < NumSp2; i++)
	{
		for (j = 0; j < NumSp2; j++)
		{
			if (pM2 < pM2num && mul2[pM2].expRow == i && mul2[pM2].coefColumn == j)
			{
				cout << setw(5) << mul2[pM2].value;
				pM2++;
			}
			else
				cout << setw(5) << 0;
		}
		cout << endl;
	}
	cout << endl;


	system("pause");
	return 0;
}


