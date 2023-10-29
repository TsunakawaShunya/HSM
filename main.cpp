/*
 * File:   main.cpp
 * Author: Tsunakawa
 *
 * Created on 2023/05/22, 10:49
 */

#include "Simulator.h"
#include <iostream>
#include <fstream>
#include <string>

int N_Tri;       // 試行回数1000000000
static const double EbN0dBmin = 0.0;    // Eb/N0 の最小値 [dB]
static const double EbN0dBmax = 35.0;       // Eb/N0 の最大値 [dB]
static const double EbN0dBstp = 1.0;        // Eb/N0 の間隔 [dB]

int MODE;       //変調方式のモード

std::string filename_ber;       // ファイル名（BER）
std::string filename_berUpperBound;       // ファイル名（BERの上界）

int berCount;       // ビット誤り数
double ber;     // Bit Error Rate
double berUpperBound;       // BER の上界

/** main部 **/
int main() 
{
    Simulator sim;

    // コンソール入力
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "design" << std::endl;
    std::cout << "(Cube design:0, Twisted Cube design:1, Crushed Twisted Cube design:2)" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> MODE;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cout << "TRIAL" << std::endl;
    std::cout << "--------------------------------------------------------------------" << std::endl;
    std::cin >> N_Tri;
    std::cout << "--------------------------------------------------------------------" << std::endl;

    // モードの設定
    switch (MODE) {
        case 0:
            sim.setPhaseCube_();
            filename_ber = "berCube_K=" + std::to_string(sim.K) + ".csv";
            filename_berUpperBound = "berCubeUpperBound_K=" + std::to_string(sim.K) + ".csv";
            break;
        case 1:
            sim.setPhaseTwistedCube_();
            filename_ber = "berTwistedCube_K=" + std::to_string(sim.K) + ".csv";
            filename_berUpperBound = "berTwistedCubeUpperBound_K=" + std::to_string(sim.K) + ".csv";
            break;
        case 2:
            sim.setPhaseCrushedTwistedCube_();
            filename_ber = "berCrushedTwistedCube_K=" + std::to_string(sim.K) + ".csv";
            filename_berUpperBound = "berCrushedTwistedCubeUpperBound_K=" + std::to_string(sim.K) + ".csv";
            break;
        default:
		    exit(1);
    }

    std::ofstream ofs_ber(filename_ber);        // 書き込み用ファイル
    std::ofstream ofs_berUpperBound(filename_berUpperBound);

    /** 実験　**/
    for(double EbN0dB = EbN0dBmin; EbN0dB <= EbN0dBmax; EbN0dB += EbN0dBstp) {
        sim.set_EbN0dB(EbN0dB);
        berCount = 0;

        // BER の計算
        for(auto tri = 0; tri < N_Tri; tri++) {
            berCount += sim.getBitErrorCount();
        }
        ber = (double)berCount / (double)N_Tri / (double)sim.K;

        berUpperBound = sim.getberUpperBoundSelective();

        ofs_ber << EbN0dB << "," << ber << std::endl;
        ofs_berUpperBound << EbN0dB << "," << berUpperBound << std::endl;

        sim.show3(EbN0dB, ber, berUpperBound);
    }
}