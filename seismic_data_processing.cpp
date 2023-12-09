#include "seismic_data_processing.h"

std::vector<std::vector<float>> dtw::ReadFloatCsvData(std::string file_location)
{
    std::ifstream file(file_location);
    std::vector<std::vector<float>> seismic_data_2d;
    std::string line;
    while (std::getline(file, line))
    {
        std::vector<float> row;
        size_t comma_position = 0;
        std::string token_string;
        float token_float;

        while ((comma_position = line.find(",")) != std::string::npos)
        {
            token_string = line.substr(0, comma_position);
            token_float = stof(token_string);
            row.push_back(token_float);
            line.erase(0, comma_position + 1);
        }
        token_float = stof(line);
        row.push_back(token_float);
        seismic_data_2d.push_back(row);
    }
    return seismic_data_2d;
}

std::vector<std::vector<std::vector<float>>> dtw::Convert2DSeismicDataTo3D(std::vector<std::vector<float>> seismic_data_2d, int seismic_inline, int seismic_xline)
{
    int total_trace = seismic_data_2d.size();
    std::vector<std::vector<std::vector<float>>> seismic_data_3d;
    std::vector<std::vector<float>> temp_inline;
    // temp_inline.reserve(951);
    std::vector<float> temp_xline;
    if (seismic_inline * seismic_xline != total_trace)
    {
        std::cout << "inline and xline error!";
        assert(0);
    }
    else
    {
        for (int i = 0; i < seismic_inline; i++)
        {
            for (int j = 0; j < seismic_xline; j++)
            {
                temp_xline = seismic_data_2d[j];
                temp_inline.push_back(temp_xline);
            }
            seismic_data_3d.push_back(temp_inline);
            temp_inline.clear();
        }
    }

    return seismic_data_3d;
}

std::vector<std::vector<int>> dtw::ReadIntCsvData(std::string file_location)
{
    std::ifstream file(file_location);
    std::vector<std::vector<int>> data;
    std::string line;
    while (std::getline(file, line))
    {
        std::vector<int> row;
        size_t comma_position = 0;
        std::string token_string;
        int token_int;

        while ((comma_position = line.find(",")) != std::string::npos)
        {
            token_string = line.substr(0, comma_position);
            token_int = stoi(token_string);
            row.push_back(token_int);
            line.erase(0, comma_position + 1);
        }
        token_int = stoi(line);
        row.push_back(token_int);
        data.push_back(row);
    }
    return data;
}

void dtw::ExpotrOneDimensionalToCsv(std::vector<int> data, std::string file_location)
{
    std::ofstream ofs;
    ofs.open(file_location,std::ios::out);
    int i;
    for( i = 0;i<data.size()-1;i++)
    {
        ofs<<data[i]<<",";
    }
    ofs<<data[i];
    return;

}

void dtw::SeismicDataNormalization(std::vector<std::vector<std::vector<float>>> *seismic_data_3d)
{
    int max = 0;
    for(int i = 0;i<(*seismic_data_3d).size();i++)
    {
        for(int j = 0;j<(*seismic_data_3d)[0].size();j++)
        {
            for(int k = 0;k<(*seismic_data_3d)[0][0].size();k++)
            if((*seismic_data_3d)[i][j][k]>max)
            {
                max = (*seismic_data_3d)[i][j][k];
            }
        }
    }
        for(int i = 0;i<(*seismic_data_3d).size();i++)
    {
        for(int j = 0;j<(*seismic_data_3d)[0].size();j++)
        {
            for(int k = 0;k<(*seismic_data_3d)[0][0].size();k++)
            {
                (*seismic_data_3d)[i][j][k] = (*seismic_data_3d)[i][j][k]/max;
            }
        }
    }
}
