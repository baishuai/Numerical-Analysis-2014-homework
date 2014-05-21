/*
* @Author: Bai Shuai
* @Date:   2014-05-20 15:52:35
* @Last Modified by:   bai
* @Last Modified time: 2014-05-21 22:44:17
 */

package main

import (
	"fmt"
	"math"
)

var (
	mA = [][]float64{
		{5, -4, 1},
		{-4, 6, -4},
		{1, -4, 7}}
	mB = [][]float64{
		{25, -41, 10, -6},
		{-41, 68, -17, 10},
		{10, -17, 5, -3},
		{-6, 10, -3, 2}}
)

// matMulLi x = Ab
func matMulLi(mta [][]float64, b, x []float64) {
	length := len(b)
	var tmp float64
	for i := 0; i < length; i++ {
		tmp = 0
		for j := 0; j < length; j++ {
			tmp += mta[i][j] * b[j]
		}
		x[i] = tmp
	}
}

func macv(v []float64) float64 {
	max := 0.0
	for i := 0; i < len(v); i++ {
		if math.Abs(v[i]) > max {
			max = math.Abs(v[i])
		}
	}
	return max
}

func lineDivF(v, u []float64, f float64) {
	for i := 0; i < len(v); i++ {
		u[i] = v[i] / f
	}
}

func Mifa(mx [][]float64, init []float64, accu float64) (lambda float64, x []float64) {

	length := len(init)
	k, k1 := 0, 1
	v := make([][]float64, 2)
	u := make([][]float64, 2)
	v[0] = make([]float64, length)
	v[1] = make([]float64, length)
	u[0] = make([]float64, length)
	u[1] = make([]float64, length)
	copy(v[0], init)
	copy(u[0], init)

	for math.Abs(macv(v[k1])-macv(v[k])) >= accu {
		//v[k1] = mx*u[k]
		matMulLi(mx, u[k], v[k1])
		//uk = vk/(vk)
		lineDivF(v[k1], u[k1], macv(v[k1]))

		k1, k = k, k1
	}
	return macv(v[k]), u[k]
}

func main() {

	lmd, x := Mifa(mA, []float64{1, 1, 1}, 1e-5)
	fmt.Println(lmd, x)
	lmd, x = Mifa(mB, []float64{1, 1, 1, 1}, 1e-5)
	fmt.Println(lmd, x)

}
