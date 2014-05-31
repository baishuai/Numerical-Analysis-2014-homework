/*
* @Author: Bai Shuai
* @Date:   2014-05-12 18:45:55
* @Last Modified by:   Bai Shuai
* @Last Modified time: 2014-05-31 15:35:18
 */

package main

import (
	"fmt"
	"math"
	"time"
)

// Hilbert 返回n阶Hilbert矩阵
func Hilbert(n int) [][]float64 {
	mxh := make([][]float64, n)
	for i := 0; i < n; i++ {
		mxh[i] = make([]float64, n)
		for j := 0; j < n; j++ {
			mxh[i][j] = 1.0 / float64(i+j+1)
		}
	}
	return mxh
}

// Comb 计算组合数C(m,n)
func Comb(m, n int) int {
	//=m!/((m-n)!*n!)
	if n == 0 {
		return 1
	}
	mnf := n + 1
	for i := n + 2; i <= m; i++ {
		mnf *= i
	}
	msnf := 1
	for i := 2; i <= m-n; i++ {
		msnf *= i
	}
	return mnf / msnf
}

// InHilbert 返回n阶Hilbert矩阵的逆矩阵
func InHilbert(n int) [][]float64 {
	imxh := make([][]float64, n)
	for i := 0; i < n; i++ {
		imxh[i] = make([]float64, n)
		for j := i; j < n; j++ {
			imxh[i][j] = float64((i + j + 1) * Comb(n+i, n-j-1) *
				Comb(n+j, n-i-1) * Comb(i+j, i) * Comb(i+j, i))
			if (i+j)%2 == 1 {
				imxh[i][j] = -imxh[i][j]
			}
		}
	}
	for i := 1; i < n; i++ {
		for j := 0; j < i; j++ {
			imxh[i][j] = imxh[j][i]
		}
	}
	return imxh
}

// NormInf 两个向量误差的无穷范数
func NormInf(x1, x2 []float64) (xo float64) {
	for i := 0; i < len(x1); i++ {
		if math.Abs(x1[i]-x2[i]) > xo {
			xo = math.Abs(x1[i] - x2[i])
		}
	}
	return xo
}

// NormMInf 矩阵m的无穷范数
func NormMInf(mx [][]float64) (x float64) {
	for i := 0; i < len(mx); i++ {
		tmp := 0.0
		for j := 0; j < len(mx); j++ {
			tmp += math.Abs(mx[i][j])
		}
		if tmp > x {
			x = tmp
		}
	}
	return x
}

// CondInfWithInverse 计算条件数
func CondInfWithInverse(mx, imx [][]float64) float64 {
	return NormMInf(mx) * NormMInf(imx)
}

// A2LL 对矩阵进行LL分解（改进平方根法）
func A2LL(mx [][]float64) [][]float64 {
	n := len(mx)
	for i := 1; i < n; i++ {
		for j := 0; j < i; j++ {
			// mx[j][i] = mx[i][j]
			// tmp := mx[j][i]
			for k := 0; k < j; k++ {
				mx[j][i] -= mx[k][i] * mx[j][k]
			}
			// mx[j][i] = tmp
		}
		for j := 0; j < i; j++ {
			mx[i][j] = mx[j][i] / mx[j][j]
		}
		//tmp := mx[i][i]
		for k := 0; k < i; k++ {
			mx[i][i] -= mx[k][i] * mx[i][k]
		}
		//mx[i][i] = tmp
	}
	return mx
}

// LLUb2x 使用平方根分解求解线性方程组 mxa*x = b
func LLUb2x(mxa [][]float64, b []float64) []float64 {
	n := len(b)
	ll := A2LL(mxa)
	y := make([]float64, n)

	y[0] = b[0]
	for i := 1; i < n; i++ {
		tmp := b[i]
		for k := 0; k < i; k++ {
			tmp -= ll[i][k] * y[k]
		}
		y[i] = tmp
	}

	x := make([]float64, n)
	x[n-1] = y[n-1] / ll[n-1][n-1]
	for i := n - 2; i >= 0; i-- {
		tmp := y[i] / ll[i][i]
		for k := i + 1; k < n; k++ {
			tmp -= ll[k][i] * x[k]
		}
		x[i] = tmp
	}
	return x
}

// A2LU 对矩阵进行LU分解
func A2LU(mx [][]float64) [][]float64 {
	n := len(mx)
	//计算U第1行，L第1列
	for i := 1; i < n; i++ {
		mx[i][0] = mx[i][0] / mx[0][0]
	}
	for r := 1; r < n; r++ {
		for i := r; i < n; i++ {
			for k := 0; k < r; k++ {
				mx[r][i] -= mx[r][k] * mx[k][i]
			}
		}
		for i := r + 1; i < n; i++ {
			for k := 0; k < r; k++ {
				mx[i][r] -= mx[i][k] * mx[k][r]
			}
			mx[i][r] /= mx[r][r]
		}
	}
	return mx
}

// LUb2x 使用LU分解法求解方程组
func LUb2x(mxa [][]float64, b []float64) []float64 {
	n := len(b)
	lu := A2LU(mxa)
	y := make([]float64, n)
	y[0] = b[0]
	for i := 1; i < n; i++ {
		y[i] = b[i]
		for k := 0; k < i; k++ {
			y[i] -= lu[i][k] * y[k]
		}
	}
	x := make([]float64, n)
	x[n-1] = y[n-1] / lu[n-1][n-1]
	for i := n - 2; i >= 0; i-- {
		x[i] = y[i]
		for k := i + 1; k < n; k++ {
			x[i] -= lu[i][k] * x[k]
		}
		x[i] /= lu[i][i]
	}
	return x
}

func martixMulLiner(hi [][]float64, x []float64) []float64 {
	n := len(x)
	b := make([]float64, n)
	for i := 0; i < n; i++ {
		for j := 0; j < n; j++ {
			b[i] += hi[i][j] * x[i]
		}
	}
	return b
}

func lab62(n int) (r, dx float64) {
	hn := Hilbert(n)
	x := make([]float64, n)
	for i := 0; i < n; i++ {
		x[i] = 1.0
	}
	b := martixMulLiner(hn, x)
	xhat := LUb2x(hn, b)
	fmt.Println(xhat)
	r = NormInf(b, martixMulLiner(Hilbert(n), xhat))
	dx = NormInf(x, xhat)
	fmt.Println(n, r, dx)
	return
}

func main() {

	fmt.Println("Hilbert(3)条件数", CondInfWithInverse(Hilbert(3), InHilbert(3)))
	fmt.Println("Hilbert(4)条件数", CondInfWithInverse(Hilbert(4), InHilbert(4)))

	x := []float64{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
	b := martixMulLiner(Hilbert(10), x)
	//h10 = A2LU(h10)
	xhat := LUb2x(Hilbert(10), b)
	fmt.Println(xhat)
	r := NormInf(b, martixMulLiner(Hilbert(10), xhat))
	dx := NormInf(x, xhat)
	fmt.Println(r, dx)

	b[0] = b[0] + 1e-7
	xhat = LUb2x(Hilbert(10), b)
	fmt.Println(xhat)
	r = NormInf(b, martixMulLiner(Hilbert(10), xhat))
	dx = NormInf(x, xhat)
	fmt.Println("1e-7", r, dx)
	fmt.Println()

	x = []float64{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
	b = martixMulLiner(Hilbert(10), x)
	xhat = LLUb2x(Hilbert(10), b)
	fmt.Println(xhat)
	r = NormInf(b, martixMulLiner(Hilbert(10), xhat))
	dx = NormInf(x, xhat)
	fmt.Println(r, dx)

	fmt.Printf("\nLU分解法与改进平方根分解法各运行100000次的时间\n")

	now := time.Now()
	for i := 0; i < 100000; i++ {
		LUb2x(Hilbert(10), b)
	}
	fmt.Println(time.Now().Sub(now).String())

	now = time.Now()
	for i := 0; i < 100000; i++ {
		LLUb2x(Hilbert(10), b)
	}
	fmt.Println(time.Now().Sub(now).String())

	fmt.Printf("\n增大希尔伯特矩阵的阶，观察误差\n")
	for i := 3; i < 25; i++ {
		r, dx := lab62(i)
		if r > 1 || dx > 1 {
			break
		}
	}
}
