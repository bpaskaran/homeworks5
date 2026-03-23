
/*
 * Homework 5 – Pointer Operations
 * CSE 220 – Spring 2026
 *
 * Name:
 * SBU ID:
 */

#include "matrix.h"

int SparseMatrix(int D[2], int M[D[0]][D[1]],
                 int S[3][(D[1] > D[0] ? D[1] : D[0])]){

    int rows = D[0];
    int cols = D[1];

    int m = (rows > cols) ? rows : cols;

    // set S to 0
    for(int i = 0; i < 3; i++){
        for(int j = 0; j < m; j++){
            S[i][j] = 0;
        }
    }

    int count = 0;

    // count non-zero
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){
            if(M[i][j] != 0){
                count++;
            }
        }
    }

    if(count > m){
        return -1;
    }

    int k = 0;

    // store sparse
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){

            if(M[i][j] != 0){
                S[0][k] = i;
                S[1][k] = j;
                S[2][k] = M[i][j];
                k++;
            }
        }
    }

    return count;
}


int HadamardProduct(const int D[6],
                    const int M[D[0]][D[1]],
                    const int N[D[2]][D[3]],
                    int A[D[4]][D[5]]){

    int mRows = D[0];
    int mCols = D[1];
    int nRows = D[2];
    int nCols = D[3];
    int aRows = D[4];
    int aCols = D[5];

    // zero A
    for(int i = 0; i < aRows; i++){
        for(int j = 0; j < aCols; j++){
            A[i][j] = 0;
        }
    }

    int rLimit = (mRows < nRows) ? mRows : nRows;
    int cLimit = (mCols < nCols) ? mCols : nCols;

    for(int i = 0; i < rLimit && i < aRows; i++){
        for(int j = 0; j < cLimit && j < aCols; j++){
            A[i][j] = M[i][j] * N[i][j];
        }
    }

    if(mRows == nRows && mCols == nCols){
        if(aRows == mRows && aCols == mCols) return 1;
        if(aRows >= mRows && aCols >= mCols) return 2;
        return -3;
    }

    if(aRows >= rLimit && aCols >= cLimit) return -1;
    return -2;
}


int Multiplication(const int D[6],
                   const int M[D[0]][D[1]],
                   const int N[D[2]][D[3]],
                   int A[D[4]][D[5]]){

    int mRows = D[0];
    int mCols = D[1];
    int nRows = D[2];
    int nCols = D[3];
    int aRows = D[4];
    int aCols = D[5];

    // zero A
    for(int i = 0; i < aRows; i++){
        for(int j = 0; j < aCols; j++){
            A[i][j] = 0;
        }
    }

    int resultRows = mRows;
    int resultCols = nCols;

    int kLimit = (mCols < nRows) ? mCols : nRows;

    for(int i = 0; i < resultRows && i < aRows; i++){
        for(int j = 0; j < resultCols && j < aCols; j++){

            int sum = 0;

            for(int k = 0; k < kLimit; k++){
                sum += M[i][k] * N[k][j];
            }

            A[i][j] = sum;
        }
    }

    if(mCols == nRows){
        if(aRows == resultRows && aCols == resultCols) return 1;
        if(aRows >= resultRows && aCols >= resultCols) return 2;
        return -3;
    }

    if(aRows >= resultRows && aCols >= resultCols) return -1;
    return -2;
}


int DiagonalSum(const int D[4],
                const int A[D[0]][D[1]],
                int DS[D[2]][D[3]]){

    int rows = D[0];
    int cols = D[1];
    int dsRows = D[2];
    int dsCols = D[3];

    // zero DS
    for(int i = 0; i < dsRows; i++){
        for(int j = 0; j < dsCols; j++){
            DS[i][j] = 0;
        }
    }

    int mainDiag = 0;
    int antiDiag = 0;

    for(int i = 0; i < rows; i++){
        for(int j = 0; j < cols; j++){

            int val = A[i][j];

            if(i + 2 < dsRows){
                DS[i + 2][0] += val;
            }

            if(j < dsCols){
                DS[1][j] += val;
            }

            if(i == j){
                mainDiag += val;
            }

            if(i + j == cols - 1){
                antiDiag += val;
            }
        }
    }

    if(dsRows > 0 && dsCols > 0){
        DS[0][0] = mainDiag;
    }

    if(rows == cols && dsCols > 1){
        DS[0][1] = antiDiag;
    }

    return 1;
}