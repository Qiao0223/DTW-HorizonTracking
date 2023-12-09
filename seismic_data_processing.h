#include <vector>
#include <string>
#include <fstream>
#include <iostream>
#include <assert.h>
#ifndef SEISMIC_DATA_PROCESSING_H
#define SEISMIC_DATA_PROCESSING_H

namespace dtw
{
    std::vector<std::vector<float>> ReadFloatCsvData(std::string file_location);

    std::vector<std::vector<std::vector<float>>> Convert2DSeismicDataTo3D(std::vector<std::vector<float>> seismic_data_2d, int seismic_inline, int seismic_xline);

    std::vector<std::vector<int>> ReadIntCsvData(std::string file_location);

    void ExpotrOneDimensionalToCsv(std::vector<int> data, std::string file_location);

    void SeismicDataNormalization(std::vector<std::vector<std::vector<float>>>* seismic_data_3d);
}
#endif