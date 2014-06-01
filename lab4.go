/*
* @Author: Bai Shuai
* @Date:   2014-05-14 21:12:12
* @Last Modified by:   bai
* @Last Modified time: 2014-06-01 15:07:41
 */

/* 运行结果
Ln(2)
复合辛普森 0.6931471848300214 4.270076070511664e-09
[0.75]
[0.7083333333333333 0.6944444444444443]
[0.6970238095238095 0.6932539682539682 0.6931746031746032]
[0.6941218503718504 0.6931545306545307 0.6931479014812348 0.6931474776448321]
[0.6933912022075268 0.6931476528194189 0.6931471942970782 0.6931471830719329 0.6931471819167451]
[0.693208208269249 0.6931472102898231 0.69314718078785 0.6931471805734177 0.6931471805636197 0.693147180562297]
Romberg外推 0.693147180562297 2.3516744107610066e-12
复合高斯 0.6931471772301613 -3.3297840040802384e-09


Pi
复合辛普森 3.1415926535877863 -2.006839139312433e-12
[3]
[3.1 3.1333333333333337]
[3.131176470588235 3.1415686274509804 3.1421176470588232]
[3.138988494491089 3.141592502458707 3.141594094125889 3.141585783761874]
[3.140941612041389 3.1415926512248227 3.1415926611425635 3.1415926383967965 3.141592665277718]
[3.1414298931749745 3.141592653552836 3.1415926537080368 3.1415926535900285 3.1415926536496097 3.1415926536382432]
[3.141551963485656 3.1415926535892167 3.1415926535916423 3.141592653589795 3.141592653589794 3.1415926535897354 3.1415926535897234]
Romberg外推 3.1415926535897234 -2.4484454763596317
复合高斯 3.1415926535915997 3.141592653589793
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

// Romberg 外推求积分
func Romberg(f func(float64) float64, a, b, accu float64) float64 {
	h := b - a
	T := make([][]float64, 0, 5)
	k := 0
	T = append(T, []float64{(f(a) + f(b)) * h / 2})
	fmt.Println(T[0])
	for {
		tmp := 0.0
		for i := 0; i < 1<<uint(k); i++ {
			tmp += f(a + h*float64(i) + h/2)
		}
		k++
		h = h / 2
		Tk := make([]float64, k+1)
		Tk[0] = h*tmp + T[k-1][0]/2

		for i := 1; i <= k; i++ {
			tmp = math.Pow(4.0, float64(i))
			Tk[i] = Tk[i-1]*tmp/(tmp-1.0) - T[k-1][i-1]/(tmp-1.0)
		}
		T = append(T, Tk)
		fmt.Println(T[k])
		if math.Abs(T[k][k]-T[k-1][k-1]) < accu {
			break
		}
	}
	return T[k][k]
}

// ComGauss 求积分
func ComGauss(f func(float64) float64, a, b float64, n int) float64 {
	h := (b - a) / float64(n)
	h2sq3 := h / (2.0 * math.Sqrt(3.0))

	var xk, ans float64
	for i := 0; i < n; i++ {
		xk = a + h*float64(i) + h/2.0
		ans += f(xk+h2sq3) + f(xk-h2sq3)
	}
	return ans * h / 2.0
}

func main() {
	fmt.Println("Ln(2)")
	ln2 := ComSimpson(fLn2, 1, 2, 26)
	fmt.Println("复合辛普森", ln2, ln2-math.Ln2)

	ln2 = Romberg(fLn2, 1, 2, 0.5e-8)
	fmt.Println("Romberg外推", ln2, ln2-math.Ln2)

	ln2 = ComGauss(fLn2, 1, 2, 25)
	fmt.Println("复合高斯", ln2, ln2-math.Ln2)

	fmt.Println("\n\nPi")
	pi := ComSimpson(fPi, 0, 1, 26)
	fmt.Println("复合辛普森", pi, pi-math.Pi)

	pi = Romberg(fPi, 0, 1, 0.5e-8)
	fmt.Println("Romberg外推", pi, ln2-math.Pi)

	pi = ComGauss(fPi, 0, 1, 25)
	fmt.Println("复合高斯", pi, math.Pi)

}
