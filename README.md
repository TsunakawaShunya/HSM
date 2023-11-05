# 超球面変調シミュレータ

## random_collection.h

`random_collection.h`はC++で疑似乱数生成をサポートするヘッダーファイルで、異なる確率分布をサポートするテンプレートクラスが含まれている。

### distribution クラス

- `T operator()() const`: このメソッドは、指定された分布に従って乱数を生成し、それを返す。

### uniform_int_distribution クラス

- `void init(T min, T max, unsigned long int seed)`: 最小値と最大値、seed値を設定して初期化

### uniform_real_distribution クラス

- `void init(T min, T max, unsigned long int seed)`: 最小値と最大値、seed値を設定して初期化

### normal_distribution クラス

- `void init(T mean, T sd, unsigned long int seed)`: 平均と標準偏差、seed値を設定して初期化

### cnormal_distribution クラス

- `void init(T mean, T sd, unsigned long int seed)`: 平均と標準偏差、seed値を設定して初期化
- `std::complex<T> operator()() const`: 実数の乱数を実部と虚部に割り当てて、複素数として返す




## Simulator.h

`Simulator.h`はC++で疑似乱数生成をサポートするヘッダーファイルで、通信シミュレーションのためのSimulatorクラスを実装。このクラスは、異なる通信環境でのビット誤り率（BER）の上限を計算するために使用。

### コンストラクタ

- `Simulator()`: コンストラクタは次元数（K）の設定、データの初期化、乱数生成器の初期化を行う。

### 位相設計メソッド

- `setPhaseCube_()`: Cube設計。
- `setPhaseTwistedCube_()`: Twisted Cube設計。
- `setPhaseCrushedTwistedCube_()`: Crushed Twisted Cube設計。

### 伝送路の初期化メソッド

- `initAWGN()`: AWGN伝送路を初期化。
- `initFlat()`: フラットフェージング伝送路を初期化。
- `initSelective()`: 選択性フェージング伝送路を初期化。

### BERの上界計算メソッド

- `getberUpperBoundAWGN()`: AWGN伝送路でのBERの上限を計算。
- `getberUpperBoundFlat()`: フラットフェージング伝送路でのBERの上限を計算。
- `getberUpperBoundSelective()`: 選択性フェージング伝送路でのBERの上限を計算。

### ビット誤り数カウントメソッド

- `getBitErrorCount()`: 指定された通信環境でのビット誤り数を計算。

### コンソール出力メソッド

- `show2(double num1, double num2)`: 2つ
- `show3(double num1, double num2, double num3)`: 3つの数値を表示。



## main.cpp

Simulator ビット誤り率シミュレーション

### 概要

多次元空間上でのシンボル設計の研究をしていて，このプログラムでは、多次元空間上での特定の設計方法（Cube design、Twisted Cube design、Crushed Twisted Cube design）におけるビット誤り率（BER）とその上界を計算するためのシミュレータ。設計する次元と設計方法を選んで通信のシミュレーションを開始する。ヘッダファイルでは選択性フェージング伝送路で通信するようになっているが、AWGN、フラットフェージング伝送路でもシミュレーション可能。

### 依存関係

このプログラムはEigen C++ライブラリを使用しており、Eigenライブラリのインストールが必要。本コードではEigenライブラリのバージョン3.4.0を使用。

### 使用方法

1. プログラムを実行し、以下の情報をコンソールで入力：
   - 次元を設定（K） 
   - 設計方法を選択（0: Cube design、1: Twisted Cube design、2: Crushed Twisted Cube design）
   - シミュレーション試行回数（N_Tri）
2. プログラムは指定した設計方法と試行回数に対するBERとBERの上界を計算し、コンソールとCSVファイルに結果が出力される。
