/*
* @Author: Bai Shuai
* @Date:   2014-05-11 20:42:29
* @Last Modified by:   Bai Shuai
* @Last Modified time: 2014-05-13 13:13:40
 */

//非线性方程求根
//二分法、牛顿法、迭代法
//1. 求解非线性方程
//2. 研究迭代函数、迭代初值对函数收敛性及收敛速度的影响
package main

import (
	"fmt"
	"math"
)

const (
	accuracy = 1e-6
)

func f1(x float64) float64 {
	return math.Pow(x, 3)*2 - math.Pow(x, 2) + 3*x - 1
}

func df1(x float64) float64 {
	return math.Pow(x, 2)*6 - 2*x + 3
}

// BinCalc 二分法解非线性方程
func BinCalc(left, right float64, f func(float64) float64, accu float64) {
	var mid, fa, fmid float64
	var k int
	fa = f(left)
	// fb = f(right)
	for {
		mid = (left + right) / 2
		fmid = f(mid)
		k++
		fmt.Printf("迭代次数 %d, 迭代值 %.10f\n", k, mid)

		if right-left < accu {
			break
		} else {
			if fmid*fa < 0 { //fa的值不需要更新
				right = mid
			} else {
				left = mid
			}
		}
	}
}

func iterf1(x float64) float64 {
	return math.Pow((x+1)/2, 1.0/3)
}

func iterf2(x float64) float64 {
	return math.Pow(x, 3)*2 - 1
}

// Iteration 实现了不动点迭代法解非线性方程
func Iteration(init float64, f func(float64) float64, times int) {
	fmt.Printf("迭代法 初值 %.10f, 次数 %d\n", init, times)
	n := 0
	for n < times {
		n++
		init = f(init)
		fmt.Printf("迭代次数 %d, 迭代值 %.10f\n", n, init)
	}
}

func f2b(x float64) float64 {
	return math.Pow(x, 3) - x - 1
}

func df2b(x float64) float64 {
	return math.Pow(x, 2)*3 - 1
}

// Newton 实现了牛顿法解非线性方程
func Newton(init float64, f, df func(float64) float64, accu float64) {
	fmt.Printf("Newton迭代法, 初值%.10f\n", init)
	k := 0
	xk, xk1 := init, init
	var fxk, dfxk, delta float64
	fxk, dfxk = f(xk), df(xk)
	for {
		k++
		xk = xk1
		xk1 = xk - fxk/dfxk
		fmt.Printf("迭代次数 %d, 迭代值 %.10f\n", k, xk1)
		fxk = f(xk1)
		dfxk = df(xk1)

		//检查精度
		if math.Abs(xk1) < 1 {
			delta = math.Abs(xk1 - xk)
		} else {
			delta = math.Abs((xk1 - xk) / xk1)
		}
		if delta < accu || math.Abs(fxk) < accu {
			break
		}
	}
}

// 打印实验信息
func init() {
	fmt.Println(`实验五-非线性方程求根
二分法、牛顿法、迭代法
1. 求解非线性方程
2. 研究迭代函数、迭代初值对函数收敛性及收敛速度的影响`)
}

func main() {
	BinCalc(-3.0, 3.0, f1, accuracy)
	fmt.Println("-------------------我是分割线----------------")
	Newton(0.0, f1, df1, accuracy)
	fmt.Println("-------------------我是分割线----------------")
	Iteration(0.0, iterf1, 10)
	fmt.Println("-------------------我是分割线----------------")
	Iteration(0.0, iterf2, 10)
	fmt.Println("-------------------我是分割线----------------")
	Newton(1.5, f2b, df2b, accuracy)
	fmt.Println("-------------------我是分割线----------------")
	Newton(0.0, f2b, df2b, accuracy)
}
