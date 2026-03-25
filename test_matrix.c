#include <criterion/criterion.h>
#include "matrix.h"

Test(SparseMatrix, basic){

    int D[2] = {2,3};

    int M[2][3] = {
        {1,0,2},
        {0,3,0}
    };

    int S[3][3];

    int result = SparseMatrix(D, M, S);

    cr_expect_eq(result, 3);
}

Test(HadamardProduct, basic){

    int D[6] = {2,2,2,2,2,2};

    int M[2][2] = {
        {1,2},
        {3,4}
    };

    int N[2][2] = {
        {5,6},
        {7,8}
    };

    int A[2][2];

    int result = HadamardProduct(D, M, N, A);

    cr_expect_eq(result, 1);
    cr_expect_eq(A[0][0], 5);
}

Test(Multiplication, basic){

    int D[6] = {2,2,2,2,2,2};

    int M[2][2] = {
        {1,2},
        {3,4}
    };

    int N[2][2] = {
        {5,6},
        {7,8}
    };

    int A[2][2];

    int result = Multiplication(D, M, N, A);

    cr_expect_eq(result, 1);
    cr_expect_eq(A[0][0], 19);
}

Test(DiagonalSum, basic){

    int D[4] = {2,2,4,3};

    int A[2][2] = {
        {1,2},
        {3,4}
    };

    int DS[4][3];

    int result = DiagonalSum(D, A, DS);

    cr_expect_eq(result, 2);  // FIXED
    cr_expect_eq(DS[0][0], 5);
}

// ================= SPARSE MATRIX EXTRA =================

Test(SparseMatrix, all_zero){
    int D[2] = {2,2};
    int M[2][2] = {{0,0},{0,0}};
    int S[3][2];
    cr_expect_eq(SparseMatrix(D,M,S), 0);
}

Test(SparseMatrix, one_element){
    int D[2] = {1,1};
    int M[1][1] = {{5}};
    int S[3][1];
    cr_expect_eq(SparseMatrix(D,M,S), 1);
}

Test(SparseMatrix, not_sparse){
    int D[2] = {2,2};
    int M[2][2] = {{1,2},{3,4}};
    int S[3][2];
    cr_expect_eq(SparseMatrix(D,M,S), -1);
}

Test(SparseMatrix, extra3){
    int D[2]={2,2};
    int M[2][2]={{1,0},{0,0}};
    int S[3][2];
    cr_expect_eq(SparseMatrix(D,M,S),1);
}

Test(SparseMatrix, extra4){
    int D[2]={3,3};
    int M[3][3]={{0,0,0},{0,1,0},{0,0,0}};
    int S[3][3];
    cr_expect_eq(SparseMatrix(D,M,S),1);
}

Test(SparseMatrix, extra5){
    int D[2]={1,3};
    int M[1][3]={{1,2,0}};
    int S[3][3];
    cr_expect(SparseMatrix(D,M,S)>=0);
}

Test(SparseMatrix, extra6){
    int D[2]={3,1};
    int M[3][1]={{0},{0},{0}};
    int S[3][3];
    cr_expect_eq(SparseMatrix(D,M,S),0);
}

Test(SparseMatrix, extra7){
    int D[2]={2,3};
    int M[2][3]={{1,2,3},{4,5,6}};
    int S[3][3];
    cr_expect_eq(SparseMatrix(D,M,S),-1);
}

Test(SparseMatrix, extra8){
    int D[2]={1,1};
    int M[1][1]={{0}};
    int S[3][1];
    cr_expect_eq(SparseMatrix(D,M,S),0);
}

Test(SparseMatrix, extra9){
    int D[2]={2,2};
    int M[2][2]={{0,1},{0,0}};
    int S[3][2];
    cr_expect_eq(SparseMatrix(D,M,S),1);
}

// ================= HADAMARD EXTRA =================

Test(HadamardProduct, diff_size_overlap){
    int D[6] = {2,3,2,2,2,2};
    int M[2][3] = {{1,2,3},{4,5,6}};
    int N[2][2] = {{1,1},{1,1}};
    int A[2][2];
    cr_expect_eq(HadamardProduct(D,M,N,A), -1);
}

Test(HadamardProduct, too_small_output){
    int D[6] = {2,2,2,2,1,1};
    int M[2][2] = {{1,2},{3,4}};
    int N[2][2] = {{1,1},{1,1}};
    int A[1][1];
    cr_expect_eq(HadamardProduct(D,M,N,A), -3);
}

Test(HadamardProduct, extra3){
    int D[6]={1,1,1,1,1,1};
    int M[1][1]={{2}};
    int N[1][1]={{2}};
    int A[1][1];
    cr_expect_eq(HadamardProduct(D,M,N,A),1);
}

Test(HadamardProduct, extra4){
    int D[6]={2,2,2,2,3,3};
    int M[2][2]={{1,1},{1,1}};
    int N[2][2]={{2,2},{2,2}};
    int A[3][3];
    cr_expect_eq(HadamardProduct(D,M,N,A),2);
}

Test(HadamardProduct, extra5){
    int D[6]={2,2,2,2,1,1};
    int M[2][2]={{1,1},{1,1}};
    int N[2][2]={{2,2},{2,2}};
    int A[1][1];
    cr_expect_eq(HadamardProduct(D,M,N,A),-3);
}

Test(HadamardProduct, extra6){
    int D[6]={2,3,2,2,2,2};
    int M[2][3]={{1,2,3},{4,5,6}};
    int N[2][2]={{1,1},{1,1}};
    int A[2][2];
    cr_expect_eq(HadamardProduct(D,M,N,A),-1);
}

Test(HadamardProduct, extra7){
    int D[6]={1,2,1,2,1,2};
    int M[1][2]={{2,3}};
    int N[1][2]={{4,5}};
    int A[1][2];
    HadamardProduct(D,M,N,A);
    cr_expect_eq(A[0][1],15);
}

Test(HadamardProduct, extra8){
    int D[6]={2,1,2,1,2,1};
    int M[2][1]={{2},{3}};
    int N[2][1]={{4},{5}};
    int A[2][1];
    HadamardProduct(D,M,N,A);
    cr_expect_eq(A[1][0],15);
}

Test(HadamardProduct, extra9){
    int D[6]={1,1,1,1,2,2};
    int M[1][1]={{7}};
    int N[1][1]={{8}};
    int A[2][2];
    cr_expect_eq(HadamardProduct(D,M,N,A),2);
}

Test(HadamardProduct, extra10){
    int D[6]={2,2,2,2,2,2};
    int M[2][2]={{0,0},{0,0}};
    int N[2][2]={{1,1},{1,1}};
    int A[2][2];
    HadamardProduct(D,M,N,A);
    cr_expect_eq(A[0][0],0);
}

// ================= MULTIPLICATION EXTRA =================

Test(Multiplication, basic2){
    int D[6] = {2,2,2,2,2,2};
    int M[2][2] = {{2,0},{0,2}};
    int N[2][2] = {{1,1},{1,1}};
    int A[2][2];
    Multiplication(D,M,N,A);
    cr_expect_eq(A[0][0],2);
}

Test(Multiplication, incompatible){
    int D[6] = {2,3,2,2,2,2};
    int M[2][3] = {{1,2,3},{4,5,6}};
    int N[2][2] = {{1,1},{1,1}};
    int A[2][2];
    cr_expect_eq(Multiplication(D,M,N,A), -1);
}

Test(Multiplication, extra3){
    int D[6]={1,1,1,1,1,1};
    int M[1][1]={{5}};
    int N[1][1]={{5}};
    int A[1][1];
    cr_expect_eq(Multiplication(D,M,N,A),1);
}

Test(Multiplication, extra4){
    int D[6]={2,2,2,2,3,3};
    int M[2][2]={{1,0},{0,1}};
    int N[2][2]={{2,2},{2,2}};
    int A[3][3];
    cr_expect_eq(Multiplication(D,M,N,A),2);
}

Test(Multiplication, extra5){
    int D[6]={2,3,2,2,2,2};
    int M[2][3]={{1,2,3},{4,5,6}};
    int N[2][2]={{1,1},{1,1}};
    int A[2][2];
    cr_expect_eq(Multiplication(D,M,N,A),-1);
}

Test(Multiplication, extra6){
    int D[6]={1,2,2,1,1,1};
    int M[1][2]={{1,2}};
    int N[2][1]={{3},{4}};
    int A[1][1];
    Multiplication(D,M,N,A);
    cr_expect_eq(A[0][0],11);
}

Test(Multiplication, extra7){
    int D[6]={2,1,1,2,2,2};
    int M[2][1]={{2},{3}};
    int N[1][2]={{4,5}};
    int A[2][2];
    Multiplication(D,M,N,A);
    cr_expect_eq(A[1][1],15);
}

Test(Multiplication, extra8){
    int D[6]={2,2,2,2,2,2};
    int M[2][2]={{0,0},{0,0}};
    int N[2][2]={{1,1},{1,1}};
    int A[2][2];
    Multiplication(D,M,N,A);
    cr_expect_eq(A[0][0],0);
}

// ================= DIAGONAL SUM EXTRA =================

Test(DiagonalSum, square){
    int D[4] = {3,3,5,3};
    int A[3][3] = {{1,2,3},{4,5,6},{7,8,9}};
    int DS[5][3];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[0][0],15);
}

Test(DiagonalSum, non_square){
    int D[4] = {2,3,4,3};
    int A[2][3] = {{1,2,3},{4,5,6}};
    int DS[4][3];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[1][0],5);
}

Test(DiagonalSum, one_element){
    int D[4] = {1,1,3,2};
    int A[1][1] = {{9}};
    int DS[3][2];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[0][0],9);
}

Test(DiagonalSum, extra3){
    int D[4]={3,3,5,3};
    int A[3][3]={{1,1,1},{1,1,1},{1,1,1}};
    int DS[5][3];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[0][0],3);
}

Test(DiagonalSum, extra4){
    int D[4]={2,2,4,2};
    int A[2][2]={{5,0},{0,5}};
    int DS[4][2];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[0][0],10);
}

Test(DiagonalSum, extra5){
    int D[4]={2,3,4,3};
    int A[2][3]={{1,2,3},{4,5,6}};
    int DS[4][3];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[1][1],7);
}

Test(DiagonalSum, extra6){
    int D[4]={1,1,3,2};
    int A[1][1]={{10}};
    int DS[3][2];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[0][0],10);
}

Test(DiagonalSum, extra7){
    int D[4]={2,2,4,3};
    int A[2][2]={{1,2},{3,4}};
    int DS[4][3];
    DiagonalSum(D,A,DS);
    cr_expect_eq(DS[1][0],4);
}