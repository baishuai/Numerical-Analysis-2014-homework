package main

import (
	"fmt"
)

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
	if n%2 == 1 {
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

//Info 打印此实验相关信息
func Info() {
	fmt.Println("实习一 截断误差与舍入误差")
	fmt.Printf("\"精确ln2值:\"\t\t %19.8f\n", ln2value)
}

func main() {
	n, v, d := ln2Accu(1, 0.0, epsilon1)
	fmt.Printf("误差限制:%e n=%d ln2=%.8f dx=%v\n", epsilon1, n, v, d)
	n, v, d = ln2Accu(n+1, v, epsilon2)
	fmt.Printf("误差限制:%e n=%d ln2=%.8f dx=%v\n", epsilon2, n, v, d)
}
