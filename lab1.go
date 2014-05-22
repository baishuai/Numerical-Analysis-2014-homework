/*
* @Author: Bai Shuai
* @Date:   2014-05-12 13:53:54
* @Last Modified by:   bai
* @Last Modified time: 2014-05-22 22:49:14
 */

/* 运行结果
截断误差与舍入误差 "精确ln2值:"	     0.693147182465
误差限制:5.000000e-06 n=40871 ln2=0.693152129650 dx=4.9471855e-06
误差限制:5.000000e-07 n=61343 ln2=0.693147659302 dx=4.7683716e-07
*/

package main

import "fmt"

const (
	ln2value float32 = 0.693147180560
	epsilon1 float32 = 0.5e-5 //1/2 x 10^-5
	epsilon2 float32 = 0.5e-6 //1/1 x 10^-6
)

//abs 绝对值函数
func abs(x float32) float32 {
	switch {
	case x < 0:
		return -x
	case x == 0:
		return 0 // return correctly abs(-0)
	}
	return x
}

//ln2 根据级数公式求ln2值，根据级数n-1时的ln2值求级数n时的ln2值
func ln2(n int, value float32) float32 {
	if n&1 == 1 {
		value += 1.0 / float32(n)
	} else {
		value -= 1.0 / float32(n)
	}
	return value
}

//ln2Accu 计算满足给定精度的ln2值，返回n值与计算的ln2值
//para：
//	n		级数公式的n值
//	value	级数公式在n-1时求到的ln2值
//	accu	要求的精度
func ln2Accu(n int, value, accu float32) (int, float32, float32) {
	value = ln2(n, value)
	for abs(value-ln2value) >= accu {
		n++
		value = ln2(n, value)
	}
	return n, value, abs(value - ln2value)
}

func main() {
	fmt.Printf("截断误差与舍入误差 \"精确ln2值:\" %16.12f\n", ln2value)
	n, v, d := ln2Accu(1, 0.0, epsilon1)
	fmt.Printf("误差限制:%e n=%d ln2=%.12f dx=%v\n", epsilon1, n, v, d)
	n, v, d = ln2Accu(n+1, v, epsilon2)
	fmt.Printf("误差限制:%e n=%d ln2=%.12f dx=%v\n", epsilon2, n, v, d)
}
