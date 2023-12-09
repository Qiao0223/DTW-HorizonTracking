#include "dtw.h"
#include "seismic_data_processing.h"
#include <iostream>
#include <ctime>
using namespace std;

int main()
{
    int inline_size = 601;
    int xline_size = 951;
    int sample_size = 288;
    int scope = 20;
    int difference = 2;
    int reference_inline = 0;
    int reference_xline = 0;
    clock_t start, end;
    vector<vector<float>> seismic_data_2d = dtw::ReadFloatCsvData("Ori_data\\F3_seismic_crop.csv");
    vector<vector<vector<float>>> seismic_data_3d = dtw::Convert2DSeismicDataTo3D(seismic_data_2d, inline_size, xline_size);
    vector<vector<int>> horizon_data_true = dtw::ReadIntCsvData("Ori_data\\horizon_data_true.csv");
    dtw::SeismicDataNormalization(&seismic_data_3d);
    // 释放seismic_data_2d内存
    seismic_data_2d.clear();
    seismic_data_2d.shrink_to_fit();

    // vector<int> c = dtw::DtwInline(&seismic_data_3d,0,475,horizon_data_true[0][475]);
    // dtw::ExpotrOneDimensionalToCsv(c,"Dtw_data\\result\\inline0_normalization.csv");
    start = clock();
    vector<vector<int>> a = dtw::DtwInlineLateralLocal(&seismic_data_3d, &horizon_data_true, reference_inline, scope);
    vector<int> b = dtw::CalculateDtwInlineLateralAccuracy(&horizon_data_true, a, reference_inline, scope, difference);
    end = clock();
    clock_t run_time = (end - start) / 1000;
    dtw::ExpotrOneDimensionalToCsv(b, "Dtw_data\\result\\inline0_amerced0.3_normalization_accuracy1.csv");
    return 0;
}