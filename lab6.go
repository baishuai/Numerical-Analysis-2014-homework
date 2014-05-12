package main

import (
	"fmt"
	"math"
)

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

func Comb(m, n int) int {
	//=m!/((m-n)!*n!)
	if n == 0 {
		return 1
	} else {
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
}

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
func NormInf(x1, x2 []float64) (xo float64) {
	for i := 0; i < len(x1); i++ {
		if math.Abs(x1[i]-x2[i]) > xo {
			xo = math.Abs(x1[i] - x2[i])
		}
	}
	return xo
}

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

func CondInfWithInverse(mx, imx [][]float64) float64 {
	return NormMInf(mx) * NormMInf(imx)
}

func A2LU(mx [][]float64) [][]float64 {
	n := len(mx)
	//计算U第1行，L第1列
	for i := 1; i < n; i++ {
		mx[i][0] = mx[i][0] / mx[0][0]
	}
	for r := 1; r < n; r++ {
		//mx[r][i]
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

func LUb2x(lu [][]float64, b []float64) []float64 {
	n := len(b)
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
	hn = A2LU(hn)
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

	h10 := Hilbert(10)
	x := []float64{1, 1, 1, 1, 1, 1, 1, 1, 1, 1}
	b := martixMulLiner(h10, x)
	h10 = A2LU(h10)
	xhat := LUb2x(h10, b)
	fmt.Println(xhat)
	r := NormInf(b, martixMulLiner(Hilbert(10), xhat))
	dx := NormInf(x, xhat)
	fmt.Println(r, dx)

	b[0] = b[0] + 1e-7
	x = LUb2x(h10, b)
	fmt.Println(x)
	fmt.Println()

	for i := 3; i < 25; i++ {
		r, dx := lab62(i)
		if r > 1 || dx > 1 {
			break
		}
	}
}
