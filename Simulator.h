/*
 * File:   Simulator.h
 * Author: Tsunakawa
 *
 * Created on 2023/04/29, 15:49
 */

#ifndef SIMULATOR_H
#define SIMULATOR_H

#define _USE_MATH_DEFINES
#include <C:/eigen-3.4.0/Eigen/Dense>
#include <C:/eigen-3.4.0/Eigen/Eigen>
#include <cmath>
#include <vector>
#include "random_collection.h"

class Simulator {
    public:
    int K;      // 次元数

    /** コンストラクタ **/
    Simulator() {
        // 次元を設定
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cout << "DIMENTIONS" << std::endl;
        std::cout << "--------------------------------------------------------------------" << std::endl;
        std::cin >> K;

        numberOfSymbols_ = pow(2, K);

        // リサイズ
        num_.resize(numberOfSymbols_);
        grayNum_.resize(numberOfSymbols_);
        symbol_.resize(numberOfSymbols_, K);
        phi_.resize(numberOfSymbols_, K - 1);
        cir_.resize(K);
        txSignal_.resize(K);
        rxSignal_.resize(K);
        distance_.resize(numberOfSymbols_);

        // データ（0 ~ 2^K - 1）をセット
        setNum_();
        
        // 乱数の設定 
        unitIntUniformRand_.init(0, numberOfSymbols_ - 1, seed1);
        unitCNormalRand_.init(0, M_SQRT1_2, seed1);
        unitNormalRand_.init(0, 1, seed2);
    }

    /** デストラクタ **/
    virtual ~Simulator() {
    }

    /** 位相設計 **/
    // Cube 設計
    void setPhaseCube_() {
        /* QPSKをK次元に拡張した位相を設計 */
        for (int i = 0; i < numberOfSymbols_; i++) {
            for (int j = 0; j <= K - 2; j++){
                if (j == 0) {
                    phi_(i, j) = pow(-1, (grayNum_[i] >> 1) & 1) * (1 + 2 * ((grayNum_[i] >> 0) & 1)) * M_PI / 4;
                }
                else {
                    phi_(i, j) = acos(pow(-1, (grayNum_[i] >> j + 1) & 1) / sqrt(j + 2));
                }
            }
        }
        setSymbol_();       // 得られた位相をもとにシンボルを設計
    }

    // Twisted Cube 設計
    void setPhaseTwistedCube_() {
        setPhaseCube_();        // Cube 設計
        // ひねった位相を設計
        for(int i = 0; i < numberOfSymbols_; i++) {
            int flag = 0;
		    for (int j = 1; j < K - 1; j++) {
			    flag ^= grayNum_[i] >> (j + 1);
		    }
		    phi_(i, 0) += (pow(-1, flag) * M_PI) / 8.0;
        }
        setSymbol_();
    }

    // Crushed Twisted Cube 設計
    void setPhaseCrushedTwistedCube_() {
        setPhaseTwistedCube_();     // Twisted Cube 設計
        // 上下面を近づけた位相を設計
        for (int i = 0; i < numberOfSymbols_; i++) {
		    double Ck = 0.0;
		    double Cl = 1.0;
		    for (int j = 1; j < K - 1; j++) {
			    Ck = 2.0 / (2.0 + (Cl / sqrt(2.0)));
			    phi_(i, j) = (M_PI / 2.0) - pow(-1, (grayNum_[i] >> j + 1)) * 0.5 * acos(2.0 * Ck - 1.0);
			    Cl *= Ck;
		    }
	    }       
        setSymbol_();
    }

    /** 伝送路の初期化 **/
    // AWGN
    void initAWGN() {
        for(int k = 0; k < K; k++) {
            cir_(k) = 1;        // 周波数応答は 1
        }
    }

    // フラットフェージング伝送路の初期化
    void initFlat() {
        cir_(0) = unitCNormalRand_();       // 周波数応答はランダム
        for(int k = 1; k < K; k++) {
            cir_(k) = cir_(0);      // すべて同じ周波数応答の値
        }
    }

    // 選択性フェージング伝送路の初期化
    void initSelective() {
        for(int k = 0; k < K; k++) {
            cir_(k) = unitCNormalRand_();
        }
    }

    // 加法性雑音の標準偏差を設定
    void set_EbN0dB(double EbN0dB) {
        noiseSD_ = sqrt(0.5 * pow(10.0, -0.1 * EbN0dB));        // Eb/N0 [dB] から変換
    }

    /** ビット誤り数をカウント **/
    int getBitErrorCount() {
        /* 伝送路の初期化
        *  AWGN 伝送路：initAWGN()
        *  フラットフェージング伝送路：initFlat()
        *  選択性フェージング伝送路：initSelective()
        */
        initSelective();

        setRxSignal_();     // 受信
        setRxBit_();        //復調

        // グレイ符号化
        rxBit_ = setGrayCode_(rxBit_);
        txBit_ = setGrayCode_(txBit_);

        // ハミング距離計算
        return hammingDistance(rxBit_, txBit_);
    }

    /** 誤り率の上界 **/

    // AWGN
    double getberUpperBoundAWGN() {
        Eigen::VectorXd s_i(K);
        Eigen::VectorXd s_j(K);

        berUpperBound_ = 0.0;

        for (int i = 0; i < numberOfSymbols_; i++) {
            for (int j = 0; j < numberOfSymbols_; j++) {
                for (int k = 0; k < K; k++) {
                    s_i(k) = symbol_(i, k);
                    s_j(k) = symbol_(j, k);
                }
                berUpperBound_ += (double)hammingDistance(setGrayCode_(i), setGrayCode_(j)) * erfc((s_i - s_j).norm() / 2.0 / (sqrt(2) * noiseSD_));
            }
        }
        return berUpperBound_ / pow(2.0, (K + 1)) / (double)K;
    }

    // フラットフェージング
        double getberUpperBoundFlat() {
        Eigen::VectorXd s_i(K);
        Eigen::VectorXd s_j(K);

        berUpperBound_ = 0.0;

        for (int i = 0; i < numberOfSymbols_; i++) {
            for (int j = 0; j < numberOfSymbols_; j++) {
                for (int k = 0; k < K; k++) {
                    s_i(k) = symbol_(i, k);
                    s_j(k) = symbol_(j, k);
                }
                berUpperBound_ += (double)hammingDistance(setGrayCode_(i), setGrayCode_(j)) * (1.0 - 1.0 / sqrt((1.0 + 4.0 * pow(noiseSD_, 2) / (s_i - s_j).squaredNorm())));
            }
        }
        return berUpperBound_ / pow(2.0, (K + 1)) / (double)K;
    }

    // 選択性フェージング
    double getberUpperBoundSelective() {
        double z = 0;
        berUpperBound_ = 0;

        // 閉じた形で表現できなかったため平均化
        for(int range = 0; range < rangeMax; range++) {
            for(int i = 0; i < numberOfSymbols_; i++) {
                for(int j = 0; j < numberOfSymbols_; j++) {
                    z = 0;
                    for(int k = 0; k < K; k++) {
                        z += norm(unitCNormalRand_()) * pow((symbol_(i, k) - symbol_(j, k)), 2);
                    }
                    z = z / (4 * pow(noiseSD_, 2));
                    berUpperBound_ += hammingDistance(setGrayCode_(i), setGrayCode_(j)) * std::erfc(sqrt(z));
                }
            }
        }

        return berUpperBound_ / pow(2, K + 1) / K / rangeMax;
    }

    /** コンソール出力 **/
    // 2カラム
    void show2(double num1, double num2) {
        std::cout << num1 << "," << num2 << std::endl;
    }

    // 3カラム
    void show3(double num1, double num2, double num3) {
        std::cout << num1 << "," << num2 << "," << num3 << std::endl;
    }


    protected:
    int numberOfSymbols_;       // シンボル数

    /** データ **/
    std::vector<int> num_;      // 数値データベクトル
    std::vector<int> grayNum_;      // グレイ符号ベクトル
    Eigen::MatrixXd symbol_;        // シンボルベクトル
    Eigen::MatrixXd phi_;       // 位相ベクトル

    /** 通信 **/
    int txBit_;     // 送信ビット
    int rxBit_;     // 受信ビット
    Eigen::VectorXd txSignal_;      // 送信信号
    Eigen::VectorXd rxSignal_;      // 受信信号
    Eigen::VectorXcd cir_;          // 周波数応答
    double noiseSD_;        // 雑音の標準偏差

    /** 乱数用 **/
    unsigned long int seed1 = 100;      // seed値
    unsigned long int seed2 = 10;
    /*Normal Distribution*/
    normal_distribution<> unitNormalRand_;      // 正規乱数
    uniform_int_distribution<> unitIntUniformRand_;     // int型一様乱数
    cnormal_distribution<> unitCNormalRand_;        // 複素正規乱数
    
    /** 距離 **/
    Eigen::VectorXd distance_;      // 送信シンボルと受信シンボルのユークリッド距離ベクトル

    /** 上界 **/
    double berUpperBound_;      // BER の上界
    static const int rangeMax = 100000000;      // 選択性での上界計算の積分範囲

    // グレイ符号化
    int setGrayCode_(int num) {
        num = num ^ (num >> 1);
        return num;
    }

    // 数値データとグレイ符号のデータのセット
    void setNum_() {
        for (int i = 0; i < numberOfSymbols_; i++) {
            num_[i] = i;
            grayNum_[i] = setGrayCode_(num_[i]);
        }
    }

    // 位相からシンボル設計
    void setSymbol_() {
        double sin_roop = 1.0;
	    // sin(Φi,K-2) × sin(Φi,K-3) × ... × sin(Φi,1)
        for (int i = 0; i < numberOfSymbols_; i++) {
            for (int j = K - 2; j >= 1; j--) {
		        sin_roop *= sin(phi_(i, j));
	        }

            for (int s = 0; s < K; s++) {
		        switch (s) {
		        case 0:
			        symbol_(i, s) = sqrt(K) * sin_roop * sin(phi_(i, 0));
			        break;
		        case 1:
			        symbol_(i, s) = sqrt(K) * sin_roop * cos(phi_(i, 0));
			        break;
		        default:
			        sin_roop = 1.0;
			        for (int t = K - 2; t >= s; t--) {
				        sin_roop *= sin(phi_(i, t));
			        }
			        symbol_(i, s) = sqrt(K) * sin_roop * cos(phi_(i, (s - 1)));
                    break;
                }
		    }
	    }
    }

    // 受信
    void setRxSignal_() {
        txBit_ = unitIntUniformRand_();     // ランダムに送信データを選ぶ
        for(int k = 0; k < K; k++) {
            txSignal_(k) = symbol_(txBit_, k);
            rxSignal_(k) = txSignal_(k) + noiseSD_ / std::abs(cir_(k)) * unitNormalRand_();     // 通信
        }
    }

    // 復調
    void setRxBit_() {
        distance_.setZero();
        // 最尤推定（最短距離のシンボルを受信シンボルとする）
        for (int i = 0; i < numberOfSymbols_; i++) {
            for (int k = 0; k < K; k++) {
                distance_(i) += norm(cir_(k)) * pow(rxSignal_(k) - symbol_(i, k), 2);
            }
        }
        Eigen::VectorXd::Index col;
        distance_.minCoeff(&col);       // 距離が最も小さいインデックスが受信ビット
        rxBit_ = (int)col;
    }

    // ハミング距離計算
    int hammingDistance(int num1, int num2) {
        int ham = 0;
        int xorResult;
        int bitMask = 1;

        xorResult = num1 ^ num2;

        for(int i = 0; i < K; i++) {
            ham += (xorResult & bitMask) >> i;
            bitMask <<= 1;
        }
        return ham;
    }
};

#endif /* SIMULATOR_H */