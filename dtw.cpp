#include "dtw.h"

void dtw::Print2DMatrix(float **matrix, int row, int columon)
{
    for (int i = 0; i < row; i++)
    {
        for (int j = 0; j < columon; j++)
        {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }
    return;
}

float dtw::CalculateDistance(float a, float b)
{
    float temp = fabs(a - b);
    float distance = temp * temp;
    return distance;
}

float **dtw::CalculateDistanceMatrix(const std::vector<float> trace_1, const std::vector<float> trace_2)
{
    int trace_1_size = trace_1.size();
    int trace_2_size = trace_2.size();
    // 构造二维数组
    float **distance_matrix = new float *[trace_1_size];
    for (int i = 0; i < trace_1_size; i++)
    {
        distance_matrix[i] = new float[trace_2_size];
    }

    for (int i = 0; i < trace_1_size; i++)
    {
        for (int j = 0; j < trace_2_size; j++)
        {
            distance_matrix[i][j] = dtw::CalculateDistance(trace_1[i], trace_2[j]);
        }
    }
    return distance_matrix;
}

void dtw::FreeSecondaryPointerMemory(float **distance_matrix, int rows)
{
    for (int i = 0; i < rows; i++)
    {
        free(distance_matrix[i]);
    }
    free(distance_matrix);
}

void dtw::CalculateCumulativeDistanceMatrix(float **distance_matrix, int trace_1_size, int trace_2_size)
{
    for (int i = 1; i < trace_1_size; i++)
    {
        for (int j = 1; j < trace_2_size; j++)
        {
            distance_matrix[i][j] = distance_matrix[i][j] +
                                    std::fmin(std::fmin(distance_matrix[i][j - 1], distance_matrix[i - 1][j]), distance_matrix[i - 1][j - 1]);
        }
    }
    return;
}

std::vector<std::vector<int>> dtw::SearchWarpingPath(float **cumulativedistance_matrix, int trace_1_size, int trace_2_size)
{
    int i = trace_1_size - 1;
    int j = trace_2_size - 1;
    std::vector<std::vector<int>> path = {{i, j}};
    while (i > 0 || j > 0)
    {
        if (i == 0)
        {
            j -= 1;
        }
        else if (j == 0)
        {
            i -= 1;
        }
        else
        {
            double temp_step = std::fmin(std::fmin(cumulativedistance_matrix[i - 1][j - 1], cumulativedistance_matrix[i - 1][j]), cumulativedistance_matrix[i][j - 1]);
            if (temp_step == cumulativedistance_matrix[i - 1][j - 1])
            {
                i -= 1;
                j -= 1;
            }
            else if (temp_step == cumulativedistance_matrix[i - 1][j])
            {
                i -= 1;
            }
            else
            {
                j -= 1;
            }
        }
        path.push_back({i, j});
    }
    return path;
}

int dtw::GetCorrespondPoint(std::vector<std::vector<int>> path, int reference_point)
{
    for (int i = 0; i < path.size(); i++)
    {
        if (reference_point == path[i][0])
        {
            return path[i][1];
        }
    }
    std::cout << "reference_point error!" << std::endl;
    assert(0);
    return 0;
}

std::vector<int> dtw::DtwInline(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                int reference_inline, int reference_xline, int reference_point)
{
    std::vector<std::vector<float>> seismic_data_inline = (*seismic_data_3d)[reference_inline];
    std::vector<float> reference_trace = seismic_data_inline[reference_xline];
    std::vector<int> dtw_result;
    std::vector<std::vector<int>> path;
    int xline_size = seismic_data_inline.size();
    int sample_size = reference_trace.size();
    for (int i = 0; i < xline_size; i++)
    {
        std::vector<float> compute_trace = seismic_data_inline[i];
        float **distance_matrix = CalculateDistanceMatrix(reference_trace, compute_trace);
        CalculateCumulativeDistanceMatrix(distance_matrix, sample_size, sample_size);
        path = SearchWarpingPath(distance_matrix, sample_size, sample_size);
        int point_temp = GetCorrespondPoint(path, reference_point);
        dtw_result.push_back(point_temp);
    }
    return dtw_result;
}

std::vector<int> dtw::DtwLateralLocal(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                      int reference_inline, int reference_xline, int reference_point, int scope)
{
    std::vector<std::vector<float>> seismic_data_inline = (*seismic_data_3d)[reference_inline];
    std::vector<float> reference_trace = seismic_data_inline[reference_xline];
    std::vector<int> dtw_result;
    std::vector<std::vector<int>> path;
    int xline_size = seismic_data_inline.size();
    int sample_size = reference_trace.size();
    std::vector<int> dtw_scope = GetDtwLateralLocalScope(reference_xline, xline_size, scope);
    int scope_start = dtw_scope[0];
    int scope_end = dtw_scope[1];
    for (int i = scope_start; i < scope_end; i++)
    {
        std::vector<float> compute_trace = seismic_data_inline[i];
        float **distance_matrix = CalculateDistanceMatrix(reference_trace, compute_trace);
        CalculateCumulativeDistanceMatrix(distance_matrix, sample_size, sample_size);
        path = SearchWarpingPath(distance_matrix, sample_size, sample_size);
        int point_temp = GetCorrespondPoint(path, reference_point);
        dtw_result.push_back(point_temp);
        FreeSecondaryPointerMemory(distance_matrix, sample_size);
    }
    return dtw_result;
}

std::vector<int> dtw::GetDtwLateralLocalScope(int reference_xline, int xline_size, int scope)
{
    std::vector<int> dtw_scope;
    int start = reference_xline - scope;
    int end = reference_xline + scope + 1;
    if (start < 0)
    {
        start = 0;
    }
    if (end > xline_size - 1)
    {
        end = xline_size;
    }
    dtw_scope.push_back(start);
    dtw_scope.push_back(end);
    return dtw_scope;
}

int dtw::DtwTwoTrace(std::vector<float> reference_trace, std::vector<float> compute_trace, int reference_point)
{
    int dtw_result;
    std::vector<std::vector<int>> path;
    int sample_size = reference_trace.size();
    float **distance_matrix = CalculateDistanceMatrix(reference_trace, compute_trace);
    CalculateCumulativeDistanceMatrix(distance_matrix, sample_size, sample_size);
    path = SearchWarpingPath(distance_matrix, sample_size, sample_size);
    dtw_result = GetCorrespondPoint(path, reference_point);
    FreeSecondaryPointerMemory(distance_matrix, sample_size);
    return dtw_result;
}

std::vector<std::vector<int>> dtw::DtwInlineLateralLocal(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                                         std::vector<std::vector<int>> *horizon_data_true, int reference_inline, int scope)
{
    std::vector<std::vector<int>> dtw_result;
    std::vector<int> dtw_result_temp;
    std::vector<float> reference_trace;
    std::vector<float> compute_trace;
    int xline_size = (*horizon_data_true)[0].size();
    int reference_point;
    for (int xline = 0; xline < xline_size; xline++)
    {
        reference_point = (*horizon_data_true)[reference_inline][xline];
        dtw_result_temp = DtwLateralLocal(seismic_data_3d, reference_inline, xline, reference_point, scope);
        dtw_result.push_back(dtw_result_temp);
    }
    return dtw_result;
}

std::vector<std::vector<std::vector<int>>> dtw::DtwAllLateralLocal(std::vector<std::vector<std::vector<float>>> *seismic_data_3d,
                                                                   std::vector<std::vector<int>> *horizon_data_true, int scope)
{
    std::vector<std::vector<std::vector<int>>> dtw_result;
    std::vector<std::vector<int>> dtw_result_temp;
    int inline_size = (*horizon_data_true).size();
    for (int i = 0; i < inline_size; i++)
    {
        dtw_result_temp = DtwInlineLateralLocal(seismic_data_3d, horizon_data_true, i, scope);
        dtw_result.push_back(dtw_result_temp);
        std::cout << "当前Inline:" << i << std::endl;
    }
    return dtw_result;
}

int dtw::CalculateDtwOneTraceLateralAccuracy(std::vector<std::vector<int>> *horizon_data_true, std::vector<int> dtw_local_result,
                                             int reference_inline, int reference_xline, int scope, int difference)
{
    int xline_size = (*horizon_data_true)[0].size();
    std::vector<int> dtw_scope = GetDtwLateralLocalScope(reference_xline, xline_size, scope);
    int scope_start = dtw_scope[0];
    int scope_end = dtw_scope[1];
    std::vector<int> horizon_data_local;
    for (int i = scope_start; i < scope_end; i++)
    {
        int temp = (*horizon_data_true)[reference_inline][i];
        horizon_data_local.push_back(temp);
    }
    std::vector<int> horizon_data_difference;
    for (int i = 0; i < (scope_end - scope_start); i++)
    {
        horizon_data_difference.push_back(abs(horizon_data_local[i] - dtw_local_result[i]));
    }
    int scope_center;
    if (horizon_data_local.size() != 2 * scope + 1)
    {
        if (reference_xline - scope < 0)
        {
            scope_center = horizon_data_local.size() - scope - 1;
        }
        else
        {
            scope_center = scope;
        }
    }
    else
    {
        scope_center = scope;
    }
    int left = 0;
    int right = 0;
    while ((scope_center - left) > 0)
    {
        if (horizon_data_difference[scope_center - left - 1] <= difference)
        {
            left++;
        }
        else
        {
            break;
        }
    }
    while ((scope_center + right) < horizon_data_local.size() - 1)
    {
        if (horizon_data_difference[scope_center + right + 1] <= difference)
        {
            right++;
        }
        else
        {
            break;
        }
    }
    if (left == scope_center)
    {
        left = scope;
    }
    if (right == horizon_data_local.size() - scope_center - 1)
    {
        right = scope;
    }
    return left < right ? left : right;
}

std::vector<int> dtw::CalculateDtwInlineLateralAccuracy(std::vector<std::vector<int>> *horizon_data_true,
                                                        std::vector<std::vector<int>> dtw_inline_local,
                                                        int reference_inline, int scope, int difference)
{
    std::vector<int> inline_lateral_accuracy;
    int xline_size = (*horizon_data_true)[0].size();
    for (int i = 0; i < xline_size; i++)
    {
        int temp = CalculateDtwOneTraceLateralAccuracy(horizon_data_true, dtw_inline_local[i],
                                                       reference_inline, i, scope, difference);
        inline_lateral_accuracy.push_back(temp);
    }
    return inline_lateral_accuracy;
}

float **dtw::CalculateAverageDistanceMatrix(float **distance_matrix, int trace_1_size, int trace_2_size, int scope)
{
    float **average_distance_matrix = new float *[trace_1_size];
    for (int i = 0; i < trace_1_size; i++)
    {
        average_distance_matrix[i] = new float[trace_2_size];
    }
    for (int i = 0; i < trace_1_size; i++)
    {
        for (int j = 0; j < trace_2_size; j++)
        {
            int start = j - scope;
            int end = j + scope + 1;
            if (start < 0)
            {
                start = 0;
            }
            if (end > trace_2_size - 1)
            {
                end = trace_2_size;
            }
            float sum = 0;
            for (int k = start; k < end; k++)
            {
                sum = sum + distance_matrix[i][k];
            }
            average_distance_matrix[i][j] = sum / (end - start);
        }
    }
    return average_distance_matrix;
}