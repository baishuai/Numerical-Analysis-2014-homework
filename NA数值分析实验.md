2014年数值分析（蔡懿慈） 实验题解

运行实验代码需要安装go  
安装方法  
* ubuntu  `sudo apt-get install golang`  
* archlinux `sudo pacman -S go`
* osx  `brew install go`  
* windows  下载二进制文件安装

每个实验独立一个go文件，相互之间没有依赖  
对于无输入参数的可直接运行`go run labx.go`执行程序,`labx.go`为文件名  
对于需要输入运行时参数的程序，先执行`go build labx.go`编译代码，
再输入`./labx -flag=xx`，具体的参数参考对应文件的说明

## 截断误差与舍入误差 [lab1.go](#file-lab1-go)

利用级数逼近在限定误差范围求ln(2)的值

## 多项式插值法 [lab2.go](#file-lab2-go) 未完成

Lagrange 差值, Runge 现象  
三次样条差值

## 曲线拟合的最小二乘法   [lab3.go](#file-lab3-go)

使用正交多项式作最小二乘拟合，这种方法不需要求解线性方程组，只用递推公式，
并且当逼近次数增加一次时，只需要把程序中循环数加1，其余不用改变，
这是目前用多项式做曲线拟合最好的计算方法


## 数值积分 [lab4.go](#file-lab4-go)

* 复合 Simpson 求积公式  
* Romberg 外推方法求积分
* 复合 Gauss 公式作近似积分

## 非线性方程求根 [lab5.go](#file-lab5-go)

二分法、牛顿法（切线法）、迭代法  
迭代函数、迭代初值对收敛性的影响

## 线性方程组的直接解法 [lab6.go](#file-lab6-go)

* 矩阵三角分解(LU 分解)的方法  
* 微小偏差对解的影响（病态矩阵）  
* [x] 平方根法(Cholesky 分解)

## 线性方程组的迭代解法 [lab7.go](#file-lab7-go) 未完成  

用Jacobi法，Gauss-Seidel法和SOR法求解线性方程组的解

## 矩阵特征值问题 [lab8.go](#file-lab8-go)

幂法求出矩阵按模最大的特征值与对应的特征向量

## 常微分方程初值问题数值解法 [lab9.go](#file-lab9-go)

经典的四阶 Runge-Kutta 公式  
四阶 Hamming 公式  
Milne 公式  
预测矫正系统

本实习题目会遇到两个问题:  
a)按给定步长积分时,积分终点与区间端点不重合。  
b)多步法改变步长。
试提出相应解决方法并实现之。
