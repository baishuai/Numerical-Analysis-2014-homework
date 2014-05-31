/*
* @Author: Bai Shuai
* @Date:   2014-05-12 18:45:55
* @Last Modified by:   bai
* @Last Modified time: 2014-05-31 16:12:03
 */

/* 运行结果
Hilbert(3)条件数 748
Hilbert(4)条件数 28374.999999999996
[0.9999999987548338 1.0000001067849733 0.9999977378614773 1.0000204794185086 0.9999026418473493 1.0002669070133299 0.9995630884702723 1.0004214011575985 0.9997791407962119 1.0000484987218523]
0.00040664622819364116 0.00043691152972769043

引入微小偏差
[1.0000099979202124 0.9995051783319866 1.0079162240487853 0.9399741640249581 1.252089684043749 0.36981472847701546 1.9602324057847962 0.12554140511875153 1.4372125629498806 0.9077026506617476]
1e-7 0.8937161403019374 0.9602324057847962

[0.9999999987517899 1.0000001068204887 0.9999977407626763 1.0000204270666115 0.999902990991373 1.0002657230669296 0.9995653379322409 1.000418976008124 0.9997805271605786 1.0000481722546402]
0.0004045525887202972 0.0004346620677591062

LU分解法与改进平方根分解法各运行100000次的时间
LU 488.518679ms
平方根 403.02228ms

增大希尔伯特矩阵的阶，观察误差
[0.9999999999999981 1.0000000000000122 0.9999999999999875]
3 1.3322676295501878e-14 1.2545520178264269e-14
[0.9999999999999775 1.0000000000002485 0.999999999999402 1.0000000000003886]
4 5.681011217006926e-13 5.979661210631093e-13
[0.9999999999999655 1.0000000000005238 0.9999999999980628 1.0000000000025913 0.9999999999988525]
5 2.292166456641098e-12 2.5912605394751154e-12
[0.999999999999228 1.000000000021937 0.9999999998517923 1.0000000003853693 0.999999999574584 1.00000000016768]
6 3.836870821061211e-10 4.254160357319847e-10
[0.9999999999944524 1.0000000002215987 0.9999999978643687 1.0000000083033338 0.9999999847780691 1.0000000131521882 0.9999999956820236]
7 1.4256008373791929e-08 1.5221930937947548e-08
[0.9999999999662691 1.00000000180906 0.9999999763726758 1.0000001278681026 0.9999996557641158 1.0000004870421637 0.9999996534271247 1.0000000977747472]
8 4.3677961769628126e-07 4.870421637104272e-07
[0.9999999997602121 1.0000000164521592 0.9999997226859642 1.0000019734507875 0.9999927795056808 1.0000147132093888 0.9999831308861128 1.0000101748012682 0.9999974890861828]
9 1.4646253767125472e-05 1.686911388720791e-05
[0.9999999987548338 1.0000001067849733 0.9999977378614773 1.0000204794185086 0.9999026418473493 1.0002669070133299 0.9995630884702723 1.0004214011575985 0.9997791407962119 1.0000484987218523]
10 0.00040664622819364116 0.00043691152972769043
[0.9999999947511966 1.0000005467463535 0.9999858683437002 1.0001575494686312 0.9990635370043287 1.0032863331278048 0.9928557892293705 1.0097264868815559 0.9919301559258117 1.0037298503490202 0.999263885025643]
11 0.008775731886447868 0.009726486881555862
[0.9999999780153819 1.0000027186586886 0.9999161379832394 1.0011249174234058 0.9918589744462714 1.0353847131912473 0.9023148128133737 1.1754194756694725 0.7957677562849788 1.1486625704207052 0.9385254722941213 1.011022484909254]
12 0.1797003800983541 0.20423224371502124
[1.00000010084578 0.9999841078014012 1.000614076745965 0.9897633935558622 1.0919750464302016 0.5010034238206807 2.740851874285045 -3.0360583613106233 7.283884268701919 -5.493529661418356 5.270874315267408 -0.618191404610789 1.2688289098048169]
13 5.828312474054975 6.493529661418356
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
