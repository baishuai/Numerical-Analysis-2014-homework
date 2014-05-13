/*
* @Author: Bai Shuai
* @Date:   2014-05-12 13:57:18
* @Last Modified by:   Bai Shuai
* @Last Modified time: 2014-05-13 13:13:37
 */

//Iterative Methods for Sparse Linear Systems
package main

import (
	"fmt"
	"math"
)

var (
	n = 100
	h = 1.0 / 100
)

func matrixA(epsilon float64) [][]float64 {
	mA := make([][]float64, n)
	eh := epsilon + h
	n2eh := -epsilon*2.0 - h
	mA[0] = make([]float64, n)
	mA[0][0] = n2eh
	mA[0][1] = eh
	for i := 1; i < n-1; i++ {
		mA[i] = make([]float64, n)
		mA[i][i-1] = epsilon
		mA[i][i] = n2eh
		mA[i][i+1] = eh
	}
	mA[n-1] = make([]float64, n)
	mA[n-1][n-2] = epsilon
	mA[n-1][n-1] = n2eh
	return mA
}

func linerB(a float64) []float64 {
	ah2 := a * h * h
	b := make([]float64, n)
	for i := 0; i < n; i++ {
		b[i] = ah2
	}
	return b
}

func accuAns(x []float64, a float64, epsilon float64) []float64 {
	xrt := make([]float64, len(x))
	for i, xi := range x {
		xrt[i] = xi*a + (1-a)*(1-math.Pow(math.E, (-xi/epsilon)))/(1-math.Pow(math.E, (-1.0/epsilon)))
	}
	return xrt
}

func Norm1(x1, x2 []float64) (xo float64) {
	for i := 0; i < len(x1); i++ {
		xo += math.Abs(x1[i] - x2[i])
	}
	return xo
}

func Jacobi(a [][]float64, b []float64, x0 []float64, accu float64) []float64 {
	size := len(b)
	xrt := make([]float64, size)
	x := make([]float64, size)
	copy(x, x0)
	var i, j int
	for {
		for i = 0; i < size; i++ {
			xrt[i] = b[i]
			for j = 0; j < i; j++ {
				xrt[i] -= a[i][j] * x[j]
			}
			for j = i + 1; j < size; j++ {
				xrt[i] -= a[i][j] * x[j]
			}
			xrt[i] /= a[i][i]
		}
		//test accuracy
		if Norm1(xrt, x) < accu {
			break
		}
		copy(x, xrt)
		// fmt.Println(xrt)
	}
	return xrt
}

func GaussSeidel(a [][]float64, b []float64, x0 []float64, accu float64) []float64 {
	size := len(b)
	xrt := make([]float64, size)
	x := make([]float64, size)
	copy(x, x0)
	var i, j int
	for {
		for i = 0; i < size; i++ {
			xrt[i] = b[i]
			for j = 0; j < i; j++ {
				xrt[i] -= a[i][j] * xrt[j]
			}
			for j = i + 1; j < size; j++ {
				xrt[i] -= a[i][j] * x[j]
			}
			xrt[i] /= a[i][i]
		}
		if Norm1(xrt, x) < accu {
			break
		}
		copy(x, xrt)
	}
	return xrt
}

func SOR(a [][]float64, b []float64, omega float64, x0 []float64, accu float64) []float64 {
	size := len(b)
	xrt := make([]float64, size)
	x := make([]float64, size)
	copy(x, x0)
	var i, j int
	for {
		for i = 0; i < size; i++ {
			xrt[i] = b[i]
			for j = 0; j < i; j++ {
				xrt[i] -= a[i][j] * xrt[j]
			}
			for j = i; j < size; j++ {
				xrt[i] -= a[i][j] * x[j]
			}
			xrt[i] = xrt[i] * omega / a[i][i]
			xrt[i] += x[i]
		}
		if Norm1(xrt, x) < accu {
			break
		}
		copy(x, xrt)
	}
	return xrt
}

func main() {

	a := matrixA(1.0)
	b := linerB(0.5)
	x := make([]float64, n)
	for i := 1; i < n; i++ {
		x[i] = float64(i) / 100
	}
	fmt.Println(x)
	fmt.Println(accuAns(x, 0.5, 1.0))

	fmt.Println("---------Jacobi-----------")
	fmt.Println(Jacobi(a, b, x, 1e-8))
	fmt.Println("---------GaussSeidel-----------")
	fmt.Println(GaussSeidel(a, b, x, 1e-8))
	fmt.Println("---------SOR-----------")
	fmt.Println(SOR(a, b, 1.1, x, 1e-8))
}
