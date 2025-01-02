// paper accepted by TGRS @Yixin Fang, but it's just simple version ...

#pragma once

#include "common_lib.h"
#include "patchwork.h"
#include "tool_color_printf.h"

#define SENSOR_HEIGHT 0.4 // FIXME: move it to yaml  # zyf

#define MIN_DIS 1.0
#define MAX_DIS 50.0
#define MIN_ANGLE 0.0
#define MAX_ANGLE 360.0
#define MIN_AZIMUTH -30.0
#define MAX_AZIMUTH 60.0

#define RANGE_RES 0.25
#define SECTOR_RES 2.0
#define AZIMUTH_RES 3.0

#define RANGE_NUM (int)std::ceil((MAX_DIS - MIN_DIS) / RANGE_RES)
#define SECTOR_NUM (int)std::ceil((MAX_ANGLE - MIN_ANGLE) / SECTOR_RES)
#define AZIMUTH_NUM (int)std::ceil((MAX_AZIMUTH - MIN_AZIMUTH) / AZIMUTH_RES)
#define BIN_NUM RANGE_NUM *SECTOR_NUM *AZIMUTH_NUM

#define PD_HEIGHT (double)(SENSOR_HEIGHT + 0.5)

#define HD_RATIO (float)0.7 // FIXME: check

#define VALID_NUM 5

// 定义点的球坐标系格式结构体
struct PointAPRI
{
    float x, y, z;        // 点的笛卡尔坐标(x,y,z)
    float range;          // 点到原点的距离(径向距离)
    float angle;          // 点在水平面上的角度(方位角)
    float azimuth;        // 点的俯仰角(天顶角)
    int range_idx = -1;   // 点所在的径向距离索引
    int sector_idx = -1;  // 点所在的扇区索引
    int azimuth_idx = -1; // 点所在的俯仰角索引
    int voxel_idx = -1;   // 点在体素点云中的索引
};

// 在哈希点云中的体素结构体定义
struct Voxel
{
    int range_idx;          // 体素在径向距离方向上的索引
    int sector_idx;         // 体素在扇区方向上的索引
    int azimuth_idx;        // 体素在俯仰角方向上的索引
    int label = -1;         // 体素的标签,默认为-1表示未标记
    PointType center;       // 体素中心点,其强度值为该点在体素点云中的索引
    std::vector<int> ptIdx; // 存储体素内点在cloud_use点云中的索引数组
    int ptVoxIdx;           // 该体素中心点在体素点云中的索引
};

// frame_SSC类用于处理单帧点云数据的分割和聚类
class SSC
{
public:
    int frame_id; // 帧ID

    // 体素网格参数
    int range_num;   // 径向距离方向的体素数量
    int sector_num;  // 扇区方向的体素数量
    int azimuth_num; // 俯仰角方向的体素数量
    int bin_num;     // 总体素数量

    std::vector<PointAPRI> apri_vec;                       // 存储点云的球坐标系表示
    std::unordered_map<int, Voxel> hash_cloud;             // 体素哈希表,key为体素索引,value为体素信息
    std::unordered_map<int, std::vector<int>> cluster_vox; // 聚类结果,key为聚类ID,value为该类包含的体素索引
    std::vector<int> PD_cluster;                           // 潜在动态物体聚类ID
    std::vector<int> HD_cluster;                           // 高度差异聚类ID
    std::vector<int> AS_cluster;                           // 辅助聚类ID

    boost::shared_ptr<PatchWork<PointType>> PatchworkGroundSeg; // 地面分割器
    pcl::PointCloud<PointType>::Ptr cloud_g;                    // 地面点云
    pcl::PointCloud<PointType>::Ptr cloud_ng;                   // 非地面点云

    pcl::PointCloud<PointType>::Ptr cloud_use; // 实际使用的点云

    pcl::PointCloud<PointType>::Ptr cloud_d;  // 动态点云
    pcl::PointCloud<PointType>::Ptr cloud_nd; // 非动态点云

    pcl::PointCloud<PointType>::Ptr cloud_vox; // 体素中心点云

    // 分配内存空间
    void allocateMemory()
    {
        PatchworkGroundSeg.reset(new PatchWork<PointType>());
        cloud_g.reset(new pcl::PointCloud<PointType>());
        cloud_ng.reset(new pcl::PointCloud<PointType>());
        cloud_use.reset(new pcl::PointCloud<PointType>());
        cloud_d.reset(new pcl::PointCloud<PointType>());
        cloud_nd.reset(new pcl::PointCloud<PointType>());
        cloud_vox.reset(new pcl::PointCloud<PointType>());
    }

    // 使用PatchWork算法提取地面点云
    void extractGroudByPatchWork(const pcl::PointCloud<PointType>::Ptr &cloud_in)
    {
        double time_pw;
        PatchworkGroundSeg->set_sensor(SENSOR_HEIGHT);
        PatchworkGroundSeg->estimate_ground(*cloud_in, *cloud_g, *cloud_ng, time_pw);
        std::cout << "Ground Extract: " << " all pt num: " << cloud_in->points.size()
                  << " sensor height: " << SENSOR_HEIGHT
                  << " ground pt num: " << cloud_g->points.size()
                  << " non-ground pt num: " << cloud_ng->points.size()
                  << " time cost(ms): " << time_pw << std::endl;
    }

    // 将点云转换为球坐标系表示
    void makeApriVec(const pcl::PointCloud<PointType>::Ptr &cloud_)
    {
        for (size_t i = 0; i < cloud_->points.size(); i++)
        {
            PointType pt = cloud_->points[i];
            float dis = pointDistance2d(pt); // 计算2D距离
            float angle = getPolarAngle(pt); // 计算极角
            float azimuth = getAzimuth(pt);  // 计算方位角

            // 过滤掉范围外的点
            if (dis < MIN_DIS || dis > MAX_DIS ||
                angle < MIN_ANGLE || angle > MAX_ANGLE ||
                azimuth < MIN_AZIMUTH || azimuth > MAX_AZIMUTH)
            {
                continue;
            }

            cloud_use->points.push_back(pt);

            // 构建球坐标系点结构
            PointAPRI apri;
            apri.x = pt.x;
            apri.y = pt.y;
            apri.z = pt.z;
            apri.range = dis;
            apri.angle = angle;
            apri.azimuth = azimuth;
            // 计算体素索引
            apri.range_idx = std::ceil((dis - MIN_DIS) / RANGE_RES) - 1;
            apri.sector_idx = std::ceil((angle - MIN_ANGLE) / SECTOR_RES) - 1;
            apri.azimuth_idx = std::ceil((azimuth - MIN_AZIMUTH) / AZIMUTH_RES) - 1;
            apri.voxel_idx = apri.azimuth_idx * RANGE_NUM * SECTOR_NUM + apri.range_idx * SECTOR_NUM + apri.sector_idx;

            if (apri.voxel_idx > BIN_NUM)
            {
                ROS_WARN("pt %d can't find its bin", (int)i);
                continue;
            }
            apri_vec.emplace_back(apri);
        }
        std::cout << "apri vec size: " << apri_vec.size() << " py use num: " << cloud_use->points.size() << std::endl;
    }

    // 构建体素哈希表
    void makeHashCloud(const std::vector<PointAPRI> &apriIn_)
    {
        std::unordered_map<int, Voxel>::iterator it_find;
        int count = 0;
        for (size_t i = 0; i < apriIn_.size(); i++)
        {
            PointAPRI apri = apriIn_[i];
            it_find = hash_cloud.find(apri.voxel_idx);
            if (it_find != hash_cloud.end())
            {
                // 如果体素已存在,添加点索引
                it_find->second.ptIdx.emplace_back(i);
            }
            else
            {
                // 如果体素不存在,创建新体素
                Voxel voxel;
                voxel.ptIdx.emplace_back(i);
                voxel.range_idx = apri.range_idx;
                voxel.sector_idx = apri.sector_idx;
                voxel.azimuth_idx = apri.azimuth_idx;
                // 计算体素中心点坐标
                float range_center = (apri.range_idx * 2 + 1) / 2 * RANGE_RES + MIN_DIS;
                float sector_center = deg2rad((apri.sector_idx * 2 + 1) / 2 * SECTOR_RES) + MIN_ANGLE;
                float azimuth_center = deg2rad((apri.azimuth_idx * 2 + 1) / 2 * AZIMUTH_RES) + deg2rad(MIN_AZIMUTH);
                voxel.center.x = range_center * std::cos(sector_center);
                voxel.center.y = range_center * std::sin(sector_center);
                voxel.center.z = range_center * std::tan(azimuth_center);
                voxel.center.intensity = apri.voxel_idx;
                cloud_vox->points.emplace_back(voxel.center);
                voxel.ptVoxIdx = count;
                count++;
                hash_cloud.insert(std::make_pair(apri.voxel_idx, voxel));
            }
        }
        std::cout << "hash cloud size: " << hash_cloud.size() << std::endl;
    }

    ~SSC() {}

    // 构造函数
    SSC(const pcl::PointCloud<PointType>::Ptr &cloud_in, const int &id)
    {
        frame_id = id;
        std::cout << ANSI_COLOR_GREEN << "current frame id: " << frame_id << ANSI_COLOR_RESET << std::endl;
        allocateMemory();
        extractGroudByPatchWork(cloud_in); // 提取地面点
        makeApriVec(cloud_ng);             // 对非地面点进行球坐标系转换
        makeHashCloud(apri_vec);           // 构建体素哈希表
    }
};

// TGRS类用于处理点云的分割、分类和跟踪
class TGRS
{
public:
    TGRS() {}  // 默认构造函数
    ~TGRS() {} // 默认析构函数

    // 聚类相关函数
    // 查找给定体素周围的邻居体素
    std::vector<int> findVoxelNeighbors(const int &range_idx_, const int &sector_idx_, const int &azimuth_idx_, int size_);
    // 合并两个聚类
    void mergeClusters(std::vector<int> &clusterIdxs_, const int &idx1_, const int &idx2_);
    // 对点云进行聚类
    void cluster(const std::vector<PointAPRI> &apri_vec_,
                 std::unordered_map<int, Voxel> &hash_cloud_,
                 std::unordered_map<int, std::vector<int>> &cluster_vox);
    // 获取点云的包围盒
    std::pair<PointType, PointType> getBoundingBoxOfCloud(const pcl::PointCloud<PointType>::Ptr &cloud_);

    // 分类相关函数
    // 根据索引从点云中提取点
    pcl::PointCloud<PointType>::Ptr getCloudByVec(const std::vector<int> &vec_, const pcl::PointCloud<PointType>::Ptr &cloud_);
    // 识别潜在的动态物体
    void recognizePD(SSC &ssc);

    // 跟踪相关函数
    // 跟踪两帧之间的动态物体
    void trackPD(SSC &ssc_pre, PointTypePose *pose_pre, SSC &ssc_next, PointTypePose *pose_next);

    // 保存相关函数
    // 保存带颜色的点云
    void saveColorCloud(SSC &ssc, const std::string &path);
};