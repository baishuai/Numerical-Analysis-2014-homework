/*
* @Author: Bai Shuai
* @Date:   2014-05-11 20:56:50
* @Last Modified by:   Bai Shuai
* @Last Modified time: 2014-05-13 13:12:54
 */

package main

import (
	"fmt"
	"math"
	"strconv"
)

var (
	mSize, nSize int       //m为点集数，n为拟合曲线次数
	xVal, omega  []float64 //拟合数据以及权重
	fx           map[float64]float64
	aStar        []float64 //拟合曲线系数
	alpha, beta  []float64 //计算Pk需要的中间值
	pFunc        []*polyFunc
)

type polyFunc struct {
	rank   int
	factor []float64
}

func (p *polyFunc) String() (str string) {
	str = strconv.FormatFloat(p.factor[0], 'f', -1, 64)
	for i := 1; i < p.rank; i++ {
		if p.factor[i] > 0 {
			str += "+"
		}
		str += strconv.FormatFloat(p.factor[i], 'f', -1, 64)
		str += "X^" + strconv.Itoa(i)
	}
	return str
}

func (p *polyFunc) Val() func(float64) float64 {
	return func(x float64) (val float64) {
		for i := 0; i < p.rank; i++ {
			val += p.factor[i] * math.Pow(x, float64(i))
		}
		return val
	}
}

func (p *polyFunc) Mult(q *polyFunc) (s *polyFunc) {
	s = new(polyFunc)
	s.rank = p.rank + q.rank - 1
	s.factor = make([]float64, s.rank)
	for i := 0; i < p.rank; i++ {
		for j := 0; j < q.rank; j++ {
			s.factor[i+j] += p.factor[i] * q.factor[j]
		}
	}
	return s
}

func (p *polyFunc) Add(q *polyFunc) (s *polyFunc) {
	s = new(polyFunc)
	if p.rank > q.rank {
		p, q = q, p
	}
	s.rank = q.rank
	s.factor = make([]float64, s.rank)
	for i := 0; i < p.rank; i++ {
		s.factor[i] = p.factor[i] + q.factor[i]
	}
	for i := p.rank; i < q.rank; i++ {
		s.factor[i] = q.factor[i]
	}
	if s.factor[s.rank-1] == 0 {
		s.rank--
	}
	return s
}

func (p *polyFunc) Neg() *polyFunc {
	for i := 0; i < p.rank; i++ {
		p.factor[i] = -p.factor[i]
	}
	return p
}

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
	pFunc = make([]*polyFunc, nSize)
	aStar = make([]float64, nSize)

	for i := 0; i < nSize; i++ {
		pFunc[i] = calcPunck(i)
		aStar[i] = innerProduct(func(x float64) float64 {
			return fx[x]
		}, pFunc[i].Val()) / innerProduct(pFunc[i].Val(), pFunc[i].Val())
	}

	fans := &polyFunc{1, []float64{0}}
	for i := 0; i < nSize; i++ {
		fans = pFunc[i].Mult(&polyFunc{1, []float64{aStar[i]}}).Add(fans)
	}

	return fans.Val(), fans.String()
}

func calcAlpha(k int) float64 {
	return innerProduct(func(x float64) float64 {
		return x * pFunc[k-1].Val()(x)
	}, pFunc[k-1].Val()) / innerProduct(pFunc[k-1].Val(), pFunc[k-1].Val())
}

func calcBeta(k int) float64 {
	return innerProduct(pFunc[k].Val(), pFunc[k].Val()) /
		innerProduct(pFunc[k-1].Val(), pFunc[k-1].Val())
}

func calcPunck(i int) (p *polyFunc) {
	if i == 0 {
		p = &polyFunc{1, []float64{1.0}}
	} else if i == 1 {
		alpha[i] = calcAlpha(i)
		p = pFunc[0].Mult(&polyFunc{2, []float64{-alpha[i], 1}})
	} else {
		alpha[i] = calcAlpha(i)
		beta[i-1] = calcBeta(i - 1)
		p = pFunc[i-1].Mult(&polyFunc{2, []float64{-alpha[i], 1}}).Add(
			pFunc[i-2].Mult(&polyFunc{1, []float64{-beta[i-1]}}))
	}
	return p
}

func main() {
	x := []float64{20, 25, 30, 35, 40, 45, 50, 55, 60}
	y := []float64{805, 985, 1170, 1365, 1570, 1790, 2030, 2300, 2610}
	fAns, fstr := Solve(x, y, nil, 4)
	fmt.Println(fstr)
	for i := 0; i < len(x); i++ {
		fmt.Println(x[i], y[i], fAns(x[i]), math.Abs(fAns(x[i])-y[i]))
	}
}
