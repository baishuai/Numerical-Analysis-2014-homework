package main

import (
	"fmt"
)

var (
	mSize, nSize int       //m为点集数，n为拟合曲线次数
	xVal, omega  []float64 //拟合数据以及权重
	fx           map[float64]float64
	aStar        []float64 //拟合曲线系数
	alpha, beta  []float64 //计算Pk需要的中间值
	pFunc        []func(float64) float64
)

func innerProduct(fa, fb func(float64) float64) (ip float64) {
	for i := 0; i < mSize; i++ {
		ip += fa(xVal[i]) * fb(xVal[i]) * omega[i]
	}
	return ip
}

func Solve(x, y, om []float64, n int) (func(float64) float64, string) {
	xVal = x
	mSize = len(xVal)
	nSize = n
	if om == nil {
		omega = make([]float64, mSize)
		for i := 0; i < mSize; i++ {
			omega[i] = 1.0
		}
	} else {
		omega = om
	}
	fx = make(map[float64]float64)
	for i := 0; i < mSize; i++ {
		fx[xVal[i]] = y[i]
	}

	alpha = make([]float64, nSize)
	beta = make([]float64, nSize)
	pFunc = make([]func(float64) float64, nSize)
	aStar = make([]float64, nSize)

	for i := 0; i < nSize; i++ {
		pFunc[i] = calcPunck(i)
		aStar[i] = innerProduct(func(x float64) float64 {
			return fx[x]
		}, pFunc[i]) / innerProduct(pFunc[i], pFunc[i])
	}

	fans := func(x float64) (val float64) {
		for i := 0; i < nSize; i++ {
			val += aStar[i] * pFunc[i](x)
		}
		return val
	}
	return fans, ""
}

func calcAlpha(k int) float64 {
	return innerProduct(func(x float64) float64 {
		return x * pFunc[k-1](x)
	}, pFunc[k-1]) / innerProduct(pFunc[k-1], pFunc[k-1])
}

func calcBeta(k int) float64 {
	return innerProduct(pFunc[k], pFunc[k]) /
		innerProduct(pFunc[k-1], pFunc[k-1])
}

func calcPunck(i int) (p func(float64) float64) {
	if i == 0 {
		p = func(x float64) float64 {
			return 1
		}
	} else if i == 1 {
		alpha[i] = calcAlpha(i)
		p = func(x float64) float64 {
			return (x - alpha[i]) * pFunc[0](x)
		}
	} else {
		alpha[i] = calcAlpha(i)
		beta[i-1] = calcBeta(i - 1)
		p = func(x float64) float64 {
			return (x-alpha[i])*pFunc[i-1](x) - beta[i-1]*pFunc[i-2](x)
		}
	}
	return p
}

func main() {
	x := []float64{20, 25, 30, 35, 40, 45, 50, 55, 60}
	y := []float64{805, 985, 1170, 1365, 1570, 1790, 2030, 2300, 2610}
	fAns, _ := Solve(x, y, nil, 4)
	for i := 0; i < len(x); i++ {
		fmt.Println(y[i], fAns(x[i]))
	}
}
