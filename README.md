# numalg

数值代数 C++ 语言实现

- 作者：Quinn
- 学号：xxx

## 简要项目目录

- include/: 头文件
- src/: 算法实现
- test/: 测试与作业程序
- notes/: 实验笔记
- plot_time_consuming.py: 绘制耗时折线图脚本

## 文件说明

### `include/`

- matrix.hpp: 定义了 `Matrix` 类，包含矩阵的基本属性与方法（行列数、范数、三角矩阵检测等）。
- operations.hpp: 定义了矩阵的基本运算，如加法、减法、乘法、转置等。
- solve.hpp: 上（下）三角矩阵求解函数声明。
- gauss.hpp: 高斯消元法相关函数声明（消元法求解、Cholesky 分解、矩阵求逆等）。
- inverse_norm_estimate.hpp: 矩阵逆的无穷范数估计（Hager 优化法）函数声明。

### `src/`

- solve.cpp: 上（下）三角矩阵求解的实现。
- gauss.cpp: 各种消元法的实现，以及列主元 Gauss 消元法求逆。
- inverse_norm_estimate.cpp: Hager 优化法估计 $\|B\|_1$（当 $B = A^{-T}$ 时即 $\|A^{-1}\|_\infty$）。

## 编译说明

**使用 Makefile 进行编译与运行**

在项目根目录执行：

```bash
make
```

以对所有程序进行编译。

### 测试程序

#### `test/test.cpp`

用于测试矩阵类的基本功能以及算法的正确性。

```bash
make run-test
```

#### `test/homework02.cpp`

作业2相关程序。

```bash
make run-hw2
```

#### `test/homework03.cpp`

作业3相关程序。

```bash
make run-hw3
```

#### `test/homework05.cpp`

作业5：解的精度估计。

- 题目1：估计 5~20 阶 Hilbert 矩阵的无穷范数条件数。
- 题目2：对特殊下三角矩阵 A_n（n=5~30），用列主元 Gauss 消元法求解并估计计算解的精度。

```bash
make run-hw5
```

#### `test/test_inverse_norm.cpp`

逆范数估计算法的单元测试。

```bash
make run-inverseNorm
```

#### `test/timeConsumingTest.cpp`

用于测试不同算法的耗时情况，并将结果输出到 CSV 文件中。

```bash
make run-timeConsumingTest
```

可通过 Make 变量传入起始值、结束值与步长（默认 10, 500, 10）：

```bash
make run-timeConsumingTest START_N=10 END_N=500 STEP=10
```

## Python 绘图

在项目根目录下建立 Python 环境，安装 matplotlib，并运行绘图脚本：

```bash
pip install matplotlib
```

使用以下命令运行程序：

```bash
python plot_time_consuming.py
```

输出文件：

- 绘图结果会写入项目的输出目录
