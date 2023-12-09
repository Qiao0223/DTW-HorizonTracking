#include <vector>
#include <math.h>
#include <iostream>
#include <assert.h>
#include <algorithm>

#ifndef DTW_H
#define DTW_H

namespace dtw
{
    // 打印数组指针形式的二维数组
    void Print2DMatrix(float **matrix, int row, int columon);

    // 计算两点之间的欧氏距离
    float CalculateDistance(float a, float b);

    // 计算两个序列之间的距离矩阵
    float **CalculateDistanceMatrix(const std::vector<float> trace_1, const std::vector<float> trace_2);

    // 释放二维指针内存
    void FreeSecondaryPointerMemory(float **distance_matrix, int rows);

    // 计算累计距离矩阵
    void CalculateCumulativeDistanceMatrix(float **distance_matrix, int trace_1_size, int trace_2_size);

    // 搜索扭曲路径
    std::vector<std::vector<int>> SearchWarpingPath(float **cumulativedistance_matrix, int trace_1_size, int trace_2_size);

    // 根据得到的扭曲路径和参考点得到另一序列的对应点
    int GetCorrespondPoint(std::vector<std::vector<int>> path, int reference_point);

    // 选取一条地震道的分层点为参考，计算该地震道所在的inline所有地震道的分层点
    std::vector<int> DtwInline(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                               int reference_inline, int reference_xline, int reference_point);

    // 选取一条地震道的分层点为参考，计算该地震道左右各scope条地震道的分层点
    std::vector<int> DtwLateralLocal(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                     int reference_inline, int reference_xline, int reference_point, int scope);

    // 获取横向局部DTW计算道的第一个和最后一个的索引值（左闭右开）
    std::vector<int> GetDtwLateralLocalScope(int reference_xline, int xline_size, int scope);

    // 选取一条地震道的分层点为参考，计算另一地震道的分层点
    int DtwTwoTrace(std::vector<float> reference_trace, std::vector<float> compute_trace, int reference_point);

    // 以某一inline中所有地震道的真实分层点计算左右范围内地震道的分层点
    std::vector<std::vector<int>> DtwInlineLateralLocal(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                                        std::vector<std::vector<int>> *horizon_data_true, int reference_inline, int scope);

    // 根据地震数据中所有地震道的真实分层点计算左右范围内地震道的分层点
    std::vector<std::vector<std::vector<int>>> DtwAllLateralLocal(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                                                  std::vector<std::vector<int>> *horizon_data_true, int scope);

    // 计算一个trace分层精度
    int CalculateDtwOneTraceLateralAccuracy(std::vector<std::vector<int>> *horizon_data_true, std::vector<int> dtw_local_result,
                                            int reference_inline, int reference_xline, int scope, int difference);

    // 计算inline分层精度
    std::vector<int> CalculateDtwInlineLateralAccuracy(std::vector<std::vector<int>> *horizon_data_true,
                                                       std::vector<std::vector<int>> dtw_inline_local, int reference_inline, int scope, int difference);

    // 计算平均距离矩阵
    float **CalculateAverageDistanceMatrix(float **distance_matrix, int trace_1_size, int trace_2_size, int scope);

}

#endif
