package main

import (
	"fmt"
)

var (
	xinput = []float64{20, 25, 30, 35, 40, 45, 50, 55, 60}
	yinput = []float64{805, 985, 1170, 1365, 1570, 1790, 2030, 2300, 2610}
)

var (
	mSize, nSize      int       //m为点集数，n为拟合曲线次数
	xVal, yVal, omega []float64 //拟合数据以及权重
	fx                map[float64]float64
	aStar             []float64 //拟合曲线系数
	alpha, beta       []float64 //计算Pk需要的中间值
	pFunc             []func(float64) float64
)

func innerProduct(fa, fb func(float64) float64) (ip float64) {
	for i := 0; i < mSize; i++ {
		ip += fa(xVal[i]) * fb(xVal[i]) * omega[i]
	}
	return ip
}


func prepare(n int) {
	xVal = xinput
	yVal = yinput
	mSize = len(xVal)
	nSize = n
	omega = make([]float64, mSize)
	fx = make(map[float64]float64)
	for i := 0; i < mSize; i++ {
		omega[i] = 1.0
		fx[xVal[i]] = yVal[i]
	}

	alpha = make([]float64, nSize)
	beta = make([]float64, nSize)
	pFunc = make([]func(float64) float64, nSize)
	aStar = make([]float64, nSize)

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

func check() {
	fcheck := func(x float64) (val float64) {
		for i := 0; i < nSize; i++ {
			val += aStar[i] * pFunc[i](x)
		}
		return val
	}
	for i := 0; i < mSize; i++ {
		fmt.Println(yVal[i], fcheck(xVal[i]))
	}

}

func main() {
	prepare(4)
	for i := 0; i < nSize; i++ {
		pFunc[i] = calcPunck(i)
		aStar[i] = innerProduct(func(x float64) float64 {
			return fx[x]
		}, pFunc[i]) / innerProduct(pFunc[i], pFunc[i])
	}

	fmt.Println("satar", aStar)
	fmt.Println("alpha", alpha)
	fmt.Println("beta", beta)
	check()
}
