/*
* @Author: Bai Shuai
* @Date:   2014-05-14 21:12:12
* @Last Modified by:   Bai Shuai
* @Last Modified time: 2014-05-15 10:53:36
 */

package main

import (
	"fmt"
	"math"
)

func fLn2(x float64) float64 {
	return 1.0 / x
}

func fPi(x float64) float64 {
	return 4.0 / (x*x + 1)
}

// ComSimpson 复合辛普森求积公式
func ComSimpson(f func(float64) float64, a, b float64, n int) float64 {
	h := (b - a) / float64(n)
	ans := f(a) + f(b) + 4*f(a+h/2)
	var xk float64
	for i := 1; i < n; i++ {
		xk = a + h*float64(i)
		ans = ans + 4*f(xk+h/2) + 2*f(xk)
	}
	ans = ans * h / 6
	return ans
}

func main() {
	ln2 := ComSimpson(fLn2, 1, 2, 26)
	fmt.Println(ln2, ln2-math.Ln2)
}
