# numalg

数值代数 C++ 语言实现

- 作者：Quinn
- 学号：xxx

## 简要项目目录

- include/: 头文件
- src/: 算法实现
- test/: 测试与作业程序
- plot_time_consuming.py: 绘制耗时折线图脚本

## 文件说明

### `include/`

- matrix.hpp: 定义了 `Matrix` 类，包含矩阵的基本属性与方法。
- operations.hpp: 定义了矩阵的基本运算，如加法、乘法等。
- solve.hpp: 最基本的上（下）三角矩阵求解函数声明。
- gauss.hpp: 不同消元法的函数声明。

### `src/`

- solve.cpp: 上（下）三角矩阵求解的实现。
- gauss.cpp: 各种消元法的实现。

## 编译说明

**使用 Makefile 进行编译与运行**

在项目根目录执行：

```bash
make
```

以对所有程序进行编译.

### `test/test.cpp`

用于测试矩阵类的基本功能，如元素访问、矩阵属性等。用于验证 `Matrix` 类以及算法的正确性。

```bash
make run-test
```

### `test/timeConsumingTest.cpp`

用于测试不同算法的耗时情况，并将结果输出到 CSV 文件中。

```bash
make run-timeConsumingTest
```

可通过 Make 变量传入起始值、结束值与步长（默认 10, 500, 10）：

```bash
make run-timeConsumingTest START_N=10 END_N=500 STEP=10
```

输出文件：
- 耗时测试结果会写入项目的输出目录

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
