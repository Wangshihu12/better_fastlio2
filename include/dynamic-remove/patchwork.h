#pragma once

#ifndef PATCHWORK_H
#define PATCHWORK_H

#include <Eigen/Dense>
#include <boost/format.hpp>

// 标记点的Z值阈值,低于此值的点被视为地面标记点
#define MARKER_Z_VALUE -2.2

// 判断平面是否足够垂直的阈值,法向量与垂直方向的夹角余弦值需大于此值
#define UPRIGHT_ENOUGH 0.55

// 判断平面是否足够平坦的阈值,平面拟合误差需小于此值
#define FLAT_ENOUGH 0.2

// 判断点的高度是否过高的阈值,高于此值的点被视为非地面点
#define TOO_HIGH_ELEVATION 0.0

// 判断平面是否倾斜过大的阈值,法向量与垂直方向的夹角余弦值小于此值表示倾斜过大
#define TOO_TILTED 1.0

// 每个patch中允许的最大点数,超过此值的点将被随机采样
#define NUM_HEURISTIC_MAX_PTS_IN_PATCH 5000 // 原为3000

using Eigen::JacobiSVD;
using Eigen::MatrixXf;
using Eigen::VectorXf;

using namespace std;

static int64_t gtm()
{
    struct timeval tm;
    gettimeofday(&tm, 0);
    // return ms
    int64_t re = (((int64_t)tm.tv_sec) * 1000 * 1000 + tm.tv_usec);
    return re;
}

/*
    @brief PathWork ROS Node.
*/
template <typename PointT>
bool point_z_cmp(PointT a, PointT b)
{
    return a.z < b.z;
}

template <typename PointT>
class PatchWork
{

public:
    typedef std::vector<pcl::PointCloud<PointT>> Ring;
    typedef std::vector<Ring> Zone;

    PatchWork()
    {

        // cout <<  " ----  pathchwork INITIALIZATION  ----" << endl;

        num_sectors_each_zone_ = {16, 32, 54, 32};
        num_rings_each_zone_ = {2, 4, 4, 4};
        elevation_thr_ = {-1.2, -0.9984, -0.851, -0.605};
        flatness_thr_ = {0.0, 0.000125, 0.000185, 0.000185};

        if (num_zones_ != 4 || num_sectors_each_zone_.size() != num_rings_each_zone_.size())
        {
            throw invalid_argument("Some parameters are wrong! Check the num_zones and num_rings/sectors_each_zone");
        }
        if (elevation_thr_.size() != flatness_thr_.size())
        {
            throw invalid_argument("Some parameters are wrong! Check the elevation/flatness_thresholds");
        }

        // cout << (boost::format("Num. sectors: %d, %d, %d, %d") % num_sectors_each_zone_[0] % num_sectors_each_zone_[1] %
        //          num_sectors_each_zone_[2] %
        //          num_sectors_each_zone_[3]).str() << endl;
        // cout << (boost::format("Num. rings: %01d, %01d, %01d, %01d") % num_rings_each_zone_[0] %
        //          num_rings_each_zone_[1] %
        //          num_rings_each_zone_[2] %
        //          num_rings_each_zone_[3]).str() << endl;
        // cout << (boost::format("elevation_thr_: %0.4f, %0.4f, %0.4f, %0.4f ") % elevation_thr_[0] % elevation_thr_[1] %
        //          elevation_thr_[2] %
        //          elevation_thr_[3]).str() << endl;
        // cout << (boost::format("flatness_thr_: %0.4f, %0.4f, %0.4f, %0.4f ") % flatness_thr_[0] % flatness_thr_[1] %
        //          flatness_thr_[2] %
        //          flatness_thr_[3]).str() << endl;
        num_rings_of_interest_ = elevation_thr_.size();

        revert_pc.reserve(NUM_HEURISTIC_MAX_PTS_IN_PATCH);
        ground_pc_.reserve(NUM_HEURISTIC_MAX_PTS_IN_PATCH);
        non_ground_pc_.reserve(NUM_HEURISTIC_MAX_PTS_IN_PATCH);
        regionwise_ground_.reserve(NUM_HEURISTIC_MAX_PTS_IN_PATCH);
        regionwise_nonground_.reserve(NUM_HEURISTIC_MAX_PTS_IN_PATCH);

        min_range_z2_ = (7 * min_range_ + max_range_) / 8.0;
        min_range_z3_ = (3 * min_range_ + max_range_) / 4.0;
        min_range_z4_ = (min_range_ + max_range_) / 2.0;

        min_ranges = {min_range_, min_range_z2_, min_range_z3_, min_range_z4_};
        ring_sizes = {(min_range_z2_ - min_range_) / num_rings_each_zone_.at(0),
                      (min_range_z3_ - min_range_z2_) / num_rings_each_zone_.at(1),
                      (min_range_z4_ - min_range_z3_) / num_rings_each_zone_.at(2),
                      (max_range_ - min_range_z4_) / num_rings_each_zone_.at(3)};
        sector_sizes = {2 * M_PI / num_sectors_each_zone_.at(0), 2 * M_PI / num_sectors_each_zone_.at(1),
                        2 * M_PI / num_sectors_each_zone_.at(2),
                        2 * M_PI / num_sectors_each_zone_.at(3)};

        for (int iter = 0; iter < num_zones_; ++iter)
        {
            Zone z;
            initialize_zone(z, num_sectors_each_zone_.at(iter), num_rings_each_zone_.at(iter));
            ConcentricZoneModel_.push_back(z);
        }
        std::cout << "\n";
    }

    void estimate_ground(
        const pcl::PointCloud<PointT> &cloudIn,
        pcl::PointCloud<PointT> &cloudOut,
        pcl::PointCloud<PointT> &cloudNonground,
        double &time_taken);

    void set_sensor(const double &height);

private:
    int num_iter_ = 3;
    int num_lpr_ = 20;
    int num_min_pts_ = 10;
    int num_rings_ = 30;
    int num_sectors_ = 108;
    int num_zones_ = 4;
    int num_rings_of_interest_ = 2.7;

    double sensor_height_;  // 1.732
    double th_seeds_ = 0.3; // 0.5
    double th_dist_ = 0.1;  // 0.125
    double max_range_ = 80.0;
    double min_range_ = 0.1;
    double uprightness_thr_ = 0.707; // 0.707
    double adaptive_seed_selection_margin_ = -1.1;
    double min_range_z2_; // 12.3625
    double min_range_z3_; // 22.025
    double min_range_z4_; // 41.35

    bool verbose_ = false;

    float d_;
    MatrixXf normal_;
    VectorXf singular_values_;
    float th_dist_d_;
    Eigen::Matrix3f cov_;
    Eigen::Vector4f pc_mean_;
    double ring_size;
    double sector_size;
    // For visualization
    bool visualize_;

    vector<int> num_sectors_each_zone_;
    vector<int> num_rings_each_zone_;

    vector<double> sector_sizes;
    vector<double> ring_sizes;
    vector<double> min_ranges;
    vector<double> elevation_thr_;
    vector<double> flatness_thr_;

    vector<Zone> ConcentricZoneModel_;

    pcl::PointCloud<PointT> revert_pc, reject_pc;
    pcl::PointCloud<PointT> ground_pc_;
    pcl::PointCloud<PointT> non_ground_pc_;

    pcl::PointCloud<PointT> regionwise_ground_;
    pcl::PointCloud<PointT> regionwise_nonground_;

    void initialize_zone(Zone &z, int num_sectors, int num_rings);

    void flush_patches_in_zone(Zone &patches, int num_sectors, int num_rings);

    double calc_principal_variance(const Eigen::Matrix3f &cov, const Eigen::Vector4f &centroid);

    double xy2theta(const double &x, const double &y);

    double xy2radius(const double &x, const double &y);

    void pc2czm(const pcl::PointCloud<PointT> &src, std::vector<Zone> &czm);

    void estimate_plane_(const pcl::PointCloud<PointT> &ground);

    void extract_piecewiseground(
        const int zone_idx, const pcl::PointCloud<PointT> &src,
        pcl::PointCloud<PointT> &dst,
        pcl::PointCloud<PointT> &non_ground_dst);

    void estimate_plane_(const int zone_idx, const pcl::PointCloud<PointT> &ground);

    void extract_initial_seeds_(
        const int zone_idx, const pcl::PointCloud<PointT> &p_sorted,
        pcl::PointCloud<PointT> &init_seeds);
};

template <typename PointT>
inline void PatchWork<PointT>::initialize_zone(Zone &z, int num_sectors, int num_rings)
{
    z.clear();
    pcl::PointCloud<PointT> cloud;
    cloud.reserve(1000);
    Ring ring;
    for (int i = 0; i < num_sectors; i++)
    {
        ring.emplace_back(cloud);
    }
    for (int j = 0; j < num_rings; j++)
    {
        z.emplace_back(ring);
    }
}

template <typename PointT>
inline void PatchWork<PointT>::flush_patches_in_zone(Zone &patches, int num_sectors, int num_rings)
{
    for (int i = 0; i < num_sectors; i++)
    {
        for (int j = 0; j < num_rings; j++)
        {
            if (!patches[j][i].points.empty())
                patches[j][i].points.clear();
        }
    }
}

template <typename PointT>
inline void PatchWork<PointT>::estimate_plane_(const pcl::PointCloud<PointT> &ground)
{
    pcl::computeMeanAndCovarianceMatrix(ground, cov_, pc_mean_);
    // Singular Value Decomposition: SVD
    Eigen::JacobiSVD<Eigen::MatrixXf> svd(cov_, Eigen::DecompositionOptions::ComputeFullU);
    singular_values_ = svd.singularValues();

    // use the least singular vector as normal
    normal_ = (svd.matrixU().col(2));
    // mean ground seeds value
    Eigen::Vector3f seeds_mean = pc_mean_.head<3>();

    // according to normal.T*[x,y,z] = -d
    d_ = -(normal_.transpose() * seeds_mean)(0, 0);
    // set distance threhold to `th_dist - d`
    th_dist_d_ = th_dist_ - d_;
}

/**
 * @brief 提取初始地面种子点。
 *
 * @tparam PointT 点云类型模板参数。
 * @param[in] zone_idx 区域索引，用于选择特定区域的种子点提取策略。
 * @param[in] p_sorted 已排序的输入点云（通常按高度从低到高排序）。
 * @param[out] init_seeds 提取出的初始地面种子点云。
 *
 * @details
 * 该函数通过以下步骤提取地面初始种子点：
 * 1. **清空初始化**：
 *    - 清空输出种子点云 `init_seeds` 的内容。
 * 2. **计算低点代表（LPR，Low Point Representative）**：
 *    - 如果 `zone_idx == 0`，过滤高度低于 `adaptive_seed_selection_margin_ * sensor_height_` 的点，确定种子点的初始索引 `init_idx`。
 *    - 计算从 `init_idx` 开始的点的高度均值，作为 LPR 高度 `lpr_height`。
 * 3. **筛选初始种子点**：
 *    - 遍历点云，将高度小于 `lpr_height + th_seeds_` 的点加入种子点云 `init_seeds`。
 *
 * @note
 * - `adaptive_seed_selection_margin_` 和 `sensor_height_`：用于动态调整种子点的高度范围。
 * - `num_lpr_`：用于计算 LPR 高度的点数。
 * - `th_seeds_`：种子点的高度阈值。
 */
template <typename PointT>
inline void PatchWork<PointT>::extract_initial_seeds_(
    const int zone_idx, const pcl::PointCloud<PointT> &p_sorted,
    pcl::PointCloud<PointT> &init_seeds)
{
    // 1. 清空初始种子点云
    init_seeds.points.clear();

    // LPR 是低点代表的均值
    double sum = 0;
    int cnt = 0;

    // 确定初始索引
    int init_idx = 0;
    if (zone_idx == 0)
    {
        for (int i = 0; i < p_sorted.points.size(); i++)
        {
            // 如果点的高度小于 adaptive_seed_selection_margin_ * sensor_height_，更新初始索引
            if (p_sorted.points[i].z < adaptive_seed_selection_margin_ * sensor_height_)
            {
                ++init_idx;
            }
            else
            {
                break;
            }
        }
    }

    // 2. 计算低点代表（LPR）的高度均值
    for (int i = init_idx; i < p_sorted.points.size() && cnt < num_lpr_; i++)
    {
        sum += p_sorted.points[i].z;
        cnt++;
    }
    // 防止除以 0 的情况
    double lpr_height = cnt != 0 ? sum / cnt : 0;

    // 3. 筛选初始地面种子点
    for (int i = 0; i < p_sorted.points.size(); i++)
    {
        if (p_sorted.points[i].z < lpr_height + th_seeds_)
        {
            init_seeds.points.push_back(p_sorted.points[i]);
        }
    }
}

/*
    @brief Velodyne pointcloud callback function. The main GPF pipeline is here.
    PointCloud SensorMsg -> Pointcloud -> z-value sorted Pointcloud
    ->error points removal -> extract ground seeds -> ground plane fit mainloop
*/

template <typename PointT>
inline void PatchWork<PointT>::estimate_ground(
    const pcl::PointCloud<PointT> &cloud_in,
    pcl::PointCloud<PointT> &cloud_out,
    pcl::PointCloud<PointT> &cloud_nonground,
    double &time_taken)
{

    static double start, t0, t1, t2, end;

    double t_total_ground = 0.0;
    double t_total_estimate = 0.0;
    // 1.Msg to pointcloud
    pcl::PointCloud<PointT> laserCloudIn;
    laserCloudIn = cloud_in;

    start = gtm();

    // 2.Sort on Z-axis value.
    sort(laserCloudIn.points.begin(), laserCloudIn.end(), point_z_cmp<PointT>);

    t0 = gtm();
    // 3.Error point removal
    // As there are some error mirror reflection under the ground,
    // here regardless point under 1.8* sensor_height
    // Sort point according to height, here uses z-axis in default
    auto it = laserCloudIn.points.begin();
    for (int i = 0; i < laserCloudIn.points.size(); i++)
    {
        if (laserCloudIn.points[i].z < -1.8 * sensor_height_)
        {
            it++;
        }
        else
        {
            break;
        }
    }
    laserCloudIn.points.erase(laserCloudIn.points.begin(), it);

    t1 = gtm();
    // 4. pointcloud -> regionwise setting
    for (int k = 0; k < num_zones_; ++k)
    {
        flush_patches_in_zone(ConcentricZoneModel_[k], num_sectors_each_zone_[k], num_rings_each_zone_[k]);
    }
    pc2czm(laserCloudIn, ConcentricZoneModel_);

    t2 = gtm();

    cloud_out.clear();
    cloud_nonground.clear();
    revert_pc.clear();
    reject_pc.clear();

    int concentric_idx = 0;
    for (int k = 0; k < num_zones_; ++k)
    {
        auto zone = ConcentricZoneModel_[k];
        for (uint16_t ring_idx = 0; ring_idx < num_rings_each_zone_[k]; ++ring_idx)
        {
            for (uint16_t sector_idx = 0; sector_idx < num_sectors_each_zone_[k]; ++sector_idx)
            {
                if (zone[ring_idx][sector_idx].points.size() > num_min_pts_)
                {
                    double t_tmp0 = gtm();
                    extract_piecewiseground(k, zone[ring_idx][sector_idx], regionwise_ground_, regionwise_nonground_);
                    double t_tmp1 = gtm();
                    t_total_ground += t_tmp1 - t_tmp0;

                    // Status of each patch
                    // used in checking uprightness, elevation, and flatness, respectively
                    const double ground_z_vec = abs(normal_(2, 0));
                    const double ground_z_elevation = pc_mean_(2, 0);
                    const double surface_variable =
                        singular_values_.minCoeff() /
                        (singular_values_(0) + singular_values_(1) + singular_values_(2));

                    double t_tmp2 = gtm();
                    if (ground_z_vec < uprightness_thr_)
                    {
                        // All points are rejected
                        cloud_nonground += regionwise_ground_;
                        cloud_nonground += regionwise_nonground_;
                    }
                    else
                    { // satisfy uprightness
                        if (concentric_idx < num_rings_of_interest_)
                        {
                            if (ground_z_elevation > elevation_thr_[ring_idx + 2 * k])
                            {
                                if (flatness_thr_[ring_idx + 2 * k] > surface_variable)
                                {
                                    if (verbose_)
                                    {
                                        std::cout << "\033[1;36m[Flatness] Recovery operated. Check "
                                                  << ring_idx + 2 * k
                                                  << "th param. flatness_thr_: " << flatness_thr_[ring_idx + 2 * k]
                                                  << " > "
                                                  << surface_variable << "\033[0m" << std::endl;
                                        revert_pc += regionwise_ground_;
                                    }
                                    cloud_out += regionwise_ground_;
                                    cloud_nonground += regionwise_nonground_;
                                }
                                else
                                {
                                    if (verbose_)
                                    {
                                        std::cout << "\033[1;34m[Elevation] Rejection operated. Check "
                                                  << ring_idx + 2 * k
                                                  << "th param. of elevation_thr_: " << elevation_thr_[ring_idx + 2 * k]
                                                  << " < "
                                                  << ground_z_elevation << "\033[0m" << std::endl;
                                        reject_pc += regionwise_ground_;
                                    }
                                    cloud_nonground += regionwise_ground_;
                                    cloud_nonground += regionwise_nonground_;
                                }
                            }
                            else
                            {
                                cloud_out += regionwise_ground_;
                                cloud_nonground += regionwise_nonground_;
                            }
                        }
                        else
                        {
                            cloud_out += regionwise_ground_;
                            cloud_nonground += regionwise_nonground_;
                        }
                    }
                    double t_tmp3 = gtm();
                    t_total_estimate += t_tmp3 - t_tmp2;
                }
            }
            ++concentric_idx;
        }
    }
    end = gtm();
    time_taken = end - start;
    //    ofstream time_txt("/home/shapelim/patchwork_time_anal.txt", std::ios::app);
    //    time_txt<<t0 - start<<" "<<t1 - t0 <<" "<<t2-t1<<" "<<t_total_ground<< " "<<t_total_estimate<<"\n";
    //    time_txt.close();
}

template <typename PointT>
inline void PatchWork<PointT>::set_sensor(const double &height)
{
    sensor_height_ = height;
}

template <typename PointT>
inline double PatchWork<PointT>::calc_principal_variance(const Eigen::Matrix3f &cov, const Eigen::Vector4f &centroid)
{
    double angle = atan2(centroid(1, 0), centroid(0, 0)); // y, x
    double c = cos(angle);
    double s = sin(angle);
    double var_x_prime = c * c * cov(0, 0) + s * s * cov(1, 1) + 2 * c * s * cov(0, 1);
    double var_y_prime = s * s * cov(0, 0) + c * c * cov(1, 1) - 2 * c * s * cov(0, 1);
    return max(var_x_prime, var_y_prime);
}

template <typename PointT>
inline double PatchWork<PointT>::xy2theta(const double &x, const double &y)
{ // 0 ~ 2 * PI
    if (y >= 0)
    {
        return atan2(y, x); // 1, 2 quadrant
    }
    else
    {
        return 2 * M_PI + atan2(y, x); // 3, 4 quadrant
    }
}

template <typename PointT>
inline double PatchWork<PointT>::xy2radius(const double &x, const double &y)
{
    return sqrt(pow(x, 2) + pow(y, 2));
}

template <typename PointT>
inline void PatchWork<PointT>::pc2czm(const pcl::PointCloud<PointT> &src, std::vector<Zone> &czm)
{

    for (auto const &pt : src.points)
    {
        int ring_idx, sector_idx;
        double r = xy2radius(pt.x, pt.y);
        if ((r <= max_range_) && (r > min_range_))
        {
            double theta = xy2theta(pt.x, pt.y);

            if (r < min_range_z2_)
            { // In First rings
                ring_idx = min(static_cast<int>(((r - min_range_) / ring_sizes[0])), num_rings_each_zone_[0] - 1);
                sector_idx = min(static_cast<int>((theta / sector_sizes[0])), num_sectors_each_zone_[0] - 1);
                czm[0][ring_idx][sector_idx].points.emplace_back(pt);
            }
            else if (r < min_range_z3_)
            {
                ring_idx = min(static_cast<int>(((r - min_range_z2_) / ring_sizes[1])), num_rings_each_zone_[1] - 1);
                sector_idx = min(static_cast<int>((theta / sector_sizes[1])), num_sectors_each_zone_[1] - 1);
                czm[1][ring_idx][sector_idx].points.emplace_back(pt);
            }
            else if (r < min_range_z4_)
            {
                ring_idx = min(static_cast<int>(((r - min_range_z3_) / ring_sizes[2])), num_rings_each_zone_[2] - 1);
                sector_idx = min(static_cast<int>((theta / sector_sizes[2])), num_sectors_each_zone_[2] - 1);
                czm[2][ring_idx][sector_idx].points.emplace_back(pt);
            }
            else
            { // Far!
                ring_idx = min(static_cast<int>(((r - min_range_z4_) / ring_sizes[3])), num_rings_each_zone_[3] - 1);
                sector_idx = min(static_cast<int>((theta / sector_sizes[3])), num_sectors_each_zone_[3] - 1);
                czm[3][ring_idx][sector_idx].points.emplace_back(pt);
            }
        }
    }
}

// For adaptive
/**
 * @brief 提取分段平面地面点云和非地面点云。
 *
 * @tparam PointT 点云类型模板参数。
 * @param[in] zone_idx 区域索引，用于提取初始地面种子点。
 * @param[in] src 输入点云。
 * @param[out] dst 输出的地面点云。
 * @param[out] non_ground_dst 输出的非地面点云。
 *
 * @details
 * 该函数通过以下步骤实现地面提取：
 * 1. **初始化**：
 *    - 清空 `ground_pc_`、`dst` 和 `non_ground_dst`，确保处理前的容器为空。
 * 2. **提取初始地面种子点**：
 *    - 调用 `extract_initial_seeds_` 提取初始地面种子点，并存储在 `ground_pc_` 中。
 * 3. **地面提取迭代**：
 *    - 对 `ground_pc_` 进行多次迭代：
 *      - 使用 `estimate_plane_` 方法估计地面平面参数。
 *      - 根据地面平面模型的点到平面距离阈值 `th_dist_d_`，区分地面点和非地面点。
 *      - 在迭代的最后一次，将符合地面条件的点加入 `dst`，非地面点加入 `non_ground_dst`。
 *
 * @note
 * - `num_iter_`：迭代次数。
 * - `th_dist_d_`：点到平面距离阈值。
 * - `normal_`：平面法向量。
 */
template <typename PointT>
inline void PatchWork<PointT>::extract_piecewiseground(
    const int zone_idx, const pcl::PointCloud<PointT> &src,
    pcl::PointCloud<PointT> &dst,
    pcl::PointCloud<PointT> &non_ground_dst)
{
    // 0. 初始化
    if (!ground_pc_.empty())
        ground_pc_.clear();
    if (!dst.empty())
        dst.clear();
    if (!non_ground_dst.empty())
        non_ground_dst.clear();

    // 1. 提取初始地面种子点
    extract_initial_seeds_(zone_idx, src, ground_pc_);

    // 2. 多次迭代提取地面
    for (int i = 0; i < num_iter_; i++)
    {
        // 估计平面模型
        estimate_plane_(ground_pc_);
        ground_pc_.clear();

        // 将点云转换为矩阵形式
        Eigen::MatrixXf points(src.points.size(), 3);
        int j = 0;
        for (auto &p : src.points)
        {
            points.row(j++) << p.x, p.y, p.z;
        }

        // 计算点到平面的距离
        Eigen::VectorXf result = points * normal_;

        // 根据阈值区分地面点和非地面点
        for (int r = 0; r < result.rows(); r++)
        {
            if (i < num_iter_ - 1)
            {
                // 非最后一次迭代，仅更新地面点云
                if (result[r] < th_dist_d_)
                {
                    ground_pc_.points.push_back(src[r]);
                }
            }
            else
            {
                // 最后一次迭代，区分地面点和非地面点
                if (result[r] < th_dist_d_)
                {
                    dst.points.push_back(src[r]);
                }
                else
                {
                    non_ground_dst.push_back(src[r]);
                }
            }
        }
    }
}
#endif
