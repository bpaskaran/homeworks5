
/*
 * Homework 5 – Pointer Operations
 * CSE 220 – Spring 2026
 *
 * Name: Brendan Paskaran
 * SBU ID: 116603319
 */

#include "matrix.h"

// SPARSE MATRIX 
int SparseMatrix(int D[2], int M[D[0]][D[1]],
                 int S[3][(D[1] > D[0] ? D[1] : D[0])]){

    int rows = D[0];
    int cols = D[1];

    // m = largest dimension of matrix
    int m = (rows > cols) ? rows : cols;

    // zero S
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < m; j++){
            *(*(S + i) + j) = 0;
        }
    }

    int count = 0;

    // count non zero
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            if(*(*(M + i) + j) != 0){
                count++;
            }
        }
    }

    if(count > m) return -1;

    int k = 0;

    // store sparse
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            int val = *(*(M + i) + j);
            if(val != 0){
                *(*(S + 0) + k) = i;
                *(*(S + 1) + k) = j;
                *(*(S + 2) + k) = val;
                k++;
            }
        }
    }

    return count;
}


// HADAMARD 
int HadamardProduct(const int D[6],
                    const int M[D[0]][D[1]],
                    const int N[D[2]][D[3]],
                    int A[D[4]][D[5]]){

    int mRows = D[0], mCols = D[1];
    int nRows = D[2], nCols = D[3];
    int aRows = D[4], aCols = D[5];

    // zero A
    for(int i = 0; i < aRows; i++){
        for(int j = 0; j < aCols; j++){
            *(*(A + i) + j) = 0;
        }
    }

    int rLimit = (mRows < nRows) ? mRows : nRows;
    int cLimit = (mCols < nCols) ? mCols : nCols;

    // compute overlap
    for(int i = 0; i < rLimit && i < aRows; i++){
        for(int j = 0; j < cLimit && j < aCols; j++){
            *(*(A + i) + j) =
                (*(*(M + i) + j)) * (*(*(N + i) + j));
        }
    }

    // same dimensions
    if(mRows == nRows && mCols == nCols){
        if(aRows == mRows && aCols == mCols) return 1;
        if(aRows >= mRows && aCols >= mCols) return 2;
        return -3;
    }

    // different dimensions
    if(aRows >= rLimit && aCols >= cLimit) return -1;
    return -2;
}


// MULTIPLICATION 
int Multiplication(const int D[6],
                   const int M[D[0]][D[1]],
                   const int N[D[2]][D[3]],
                   int A[D[4]][D[5]]){

    int mRows = D[0], mCols = D[1];
    int nRows = D[2], nCols = D[3];
    int aRows = D[4], aCols = D[5];

    // zero A
    for(int i = 0; i < aRows; i++){
        for(int j = 0; j < aCols; j++){
            *(*(A + i) + j) = 0;
        }
    }

    int resultRows = mRows;
    int resultCols = nCols;


    int kLimit = mCols;

    if(nRows < mCols){
    kLimit = nRows;
    }


    // compute
    for(int i = 0; i < resultRows && i < aRows; i++){
        for(int j = 0; j < resultCols && j < aCols; j++){

            int sum = 0;

            for(int k = 0; k < kLimit; k++){
                sum += (*(*(M + i) + k)) * (*(*(N + k) + j));
            }

            *(*(A + i) + j) = sum;
        }
    }

    // compatible
    if(mCols == nRows){
        if(aRows == resultRows && aCols == resultCols) return 1;
        if(aRows >= resultRows && aCols >= resultCols) return 2;
        return -3;
    }

    // not compatible
    if(aRows >= resultRows && aCols >= resultCols) return -1;
    return -2;
}


// DIAGONAL SUM 
int DiagonalSum(const int D[4],
                const int A[D[0]][D[1]],
                int DS[D[2]][D[3]]){

    int rows = D[0], cols = D[1];
    int dsRows = D[2], dsCols = D[3];

    // zero DS
    for(int i = 0; i < dsRows; i++){
        for(int j = 0; j < dsCols; j++){
            *(*(DS + i) + j) = 0;
        }
    }

    int mainDiag = 0, antiDiag = 0;

    // compute sums
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){

            int val = *(*(A + i) + j);

            // row sums → DS[i+2][0]
            if(i + 2 < dsRows){
                *(*(DS + (i + 2)) + 0) += val;
            }

            // column sums → DS[1][j]
            if(j < dsCols){
                *(*(DS + 1) + j) += val;
            }

            if(i == j) mainDiag += val;
            if(i + j == cols - 1) antiDiag += val;
        }
    }

    if(dsRows > 0 && dsCols > 0){
        *(*(DS + 0) + 0) = mainDiag;
    }

    if(rows == cols && dsCols > 1){
        *(*(DS + 0) + 1) = antiDiag;
    }

    // minimum required size
    int needRows = rows + 2;

    int needCols = cols;
    if(rows > cols){
    needCols = rows;
    }

    if(dsRows == needRows && dsCols == needCols) return 1;
    if(dsRows >= needRows && dsCols >= needCols) return 2;
    return -1;
    }