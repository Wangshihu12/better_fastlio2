#include "tgrs.h"
#include "tictoc.hpp"

/**
 * @brief 合并两个聚类
 * @param clusterIdxs_ 聚类标签数组
 * @param idx1_ 第一个聚类的标签
 * @param idx2_ 第二个聚类的标签
 * @details 将标签为idx1_的聚类合并到标签为idx2_的聚类中
 */
void TGRS::mergeClusters(std::vector<int> &clusterIdxs_, const int &idx1_, const int &idx2_)
{
    for (int i = 0; i < clusterIdxs_.size(); i++)
    {
        if (clusterIdxs_[i] == idx1_)
        {
            clusterIdxs_[i] = idx2_;
        }
    }
}

/**
 * @brief 查找给定体素的邻近体素
 * @param range_idx_ 距离索引
 * @param sector_idx_ 扇区索引 
 * @param azimuth_idx_ 方位角索引
 * @param size_ 搜索范围大小
 * @return 返回邻近体素的索引数组
 * @details 
 * 1. 如果距离索引大于RANGE_NUM的60%,则将搜索范围限制为1
 * 2. 在给定的三维搜索范围内遍历所有体素:
 *    - 检查x,y,z方向的索引是否越界
 *    - 将有效的邻近体素索引添加到结果数组中
 * 3. 体素索引的计算方式为: x * SECTOR_NUM + y + z * RANGE_NUM * SECTOR_NUM
 */
std::vector<int> TGRS::findVoxelNeighbors(const int &range_idx_, const int &sector_idx_, const int &azimuth_idx_, int size_)
{
    std::vector<int> neighborIdxs;
    if (range_idx_ > RANGE_NUM * 0.6)
    {
        size_ = 1;
    }
    for (int x = range_idx_ - size_; x <= range_idx_ + size_; x++)
    {
        if (x > RANGE_NUM - 1 || x < 0)
        {
            continue;
        }
        for (int y = sector_idx_ - size_; y <= sector_idx_ + size_; y++)
        {
            if (y > SECTOR_NUM - 1 || y < 0)
            {
                continue;
            }
            for (int z = azimuth_idx_ - size_; z <= azimuth_idx_ + size_; z++)
            {
                if (z > AZIMUTH_NUM - 1 || z < 0)
                {
                    continue;
                }
                neighborIdxs.emplace_back(x * SECTOR_NUM + y + z * RANGE_NUM * SECTOR_NUM);
            }
        }
    }
    return neighborIdxs;
}

/**
 * @brief 对点云进行聚类
 * @param apri_vec_ 输入的APRI点云向量
 * @param hash_cloud_ 体素哈希表
 * @param cluster_vox 聚类结果的体素映射
 */
void TGRS::cluster(const std::vector<PointAPRI> &apri_vec_,
                   std::unordered_map<int, Voxel> &hash_cloud_,
                   std::unordered_map<int, std::vector<int>> &cluster_vox)
{
    // 初始化聚类ID从4开始
    int cluster_name = 4;
    // 初始化点云聚类标签数组,初始值为-1
    std::vector<int> clusterIdxs = std::vector<int>(apri_vec_.size(), -1);

    TicToc cluster_t;

    // 对每个点进行聚类
    for (int i = 0; i < apri_vec_.size(); i++)
    {
        PointAPRI apri = apri_vec_[i];
        std::unordered_map<int, Voxel>::iterator it_find1;
        std::unordered_map<int, Voxel>::iterator it_find2;
        std::vector<int> neighbors; // 存储邻近点的索引

        // 在哈希表中查找当前点所在体素
        it_find1 = hash_cloud_.find(apri.voxel_idx);
        if (it_find1 != hash_cloud_.end())
        {
            // 获取邻近体素
            std::vector<int> neighbor = findVoxelNeighbors(apri.range_idx, apri.sector_idx, apri.azimuth_idx, 1);
            // 遍历邻近体素
            for (int k = 0; k < neighbor.size(); k++)
            {
                it_find2 = hash_cloud_.find(neighbor[k]);
                if (it_find2 != hash_cloud_.end())
                {
                    addVec(neighbors, it_find2->second.ptIdx);
                }
            }
        }

        neighbors.swap(neighbors);
        // 如果存在邻近点
        if (neighbors.size() > 0)
        {
            for (int n = 0; n < neighbors.size(); n++)
            {
                int oc = clusterIdxs[i];            // 当前点的聚类ID
                int nc = clusterIdxs[neighbors[n]]; // 邻近点的聚类ID

                // 如果当前点和邻近点都已经有聚类ID
                if (oc != -1 && nc != -1)
                {
                    if (oc != nc)
                    {
                        mergeClusters(clusterIdxs, oc, nc); // 合并聚类
                    }
                }
                else
                {
                    // 如果邻近点有聚类ID,当前点继承该ID
                    if (nc != -1)
                    {
                        clusterIdxs[i] = nc;
                    }
                    else
                    {
                        // 如果当前点有聚类ID,邻近点继承该ID
                        if (oc != -1)
                        {
                            clusterIdxs[neighbors[n]] = oc;
                        }
                    }
                }
            }
        }

        // 如果当前点仍未分配聚类ID,创建新的聚类
        if (clusterIdxs[i] == -1)
        {
            cluster_name++;                // 新聚类ID
            clusterIdxs[i] = cluster_name; // 为当前点分配聚类ID
            // 为邻近点分配相同的聚类ID
            for (int m = 0; m < neighbors.size(); m++)
            {
                clusterIdxs[neighbors[m]] = cluster_name;
            }
        }
    }

    // 构建体素聚类映射
    std::unordered_map<int, std::vector<int>>::iterator it_v;
    for (size_t i = 0; i < clusterIdxs.size(); i++)
    {
        it_v = cluster_vox.find(clusterIdxs[i]);
        if (it_v != cluster_vox.end())
        {
            // 如果聚类ID已存在,添加体素索引
            it_v->second.emplace_back(apri_vec_[i].voxel_idx);
            hash_cloud_[apri_vec_[i].voxel_idx].label = it_v->first;
        }
        else
        {
            // 如果是新的聚类ID,创建新的体素向量
            std::vector<int> vox_vec;
            vox_vec.emplace_back(apri_vec_[i].voxel_idx);
            cluster_vox.insert(std::make_pair(clusterIdxs[i], vox_vec));
            hash_cloud_[apri_vec_[i].voxel_idx].label = clusterIdxs[i];
        }
    }

    // 对每个聚类的体素进行采样
    for (auto &it : cluster_vox)
    {
        sampleVec(it.second);
    }
}

/**
 * @brief 获取点云的包围盒
 * @param[in] cloud_ 输入点云
 * @return 包含最小点和最大点的pair
 * @details 使用pcl::getMinMax3D获取点云中所有点的最小和最大坐标值,
 *          返回一个pair,first为最小点,second为最大点
 */
std::pair<PointType, PointType> TGRS::getBoundingBoxOfCloud(const pcl::PointCloud<PointType>::Ptr &cloud_)
{
    PointType point_min, point_max;
    pcl::getMinMax3D(*cloud_, point_min, point_max);
    return std::make_pair(point_min, point_max);
}

/**
 * @brief 根据索引向量获取点云子集
 * @param[in] vec_ 点云索引向量
 * @param[in] cloud_ 输入点云
 * @return 包含指定索引点的点云子集
 * @details 遍历输入的索引向量,从原始点云中提取对应索引的点,
 *          构建一个新的点云子集并返回
 */
pcl::PointCloud<PointType>::Ptr TGRS::getCloudByVec(const std::vector<int> &vec_, const pcl::PointCloud<PointType>::Ptr &cloud_)
{
    pcl::PointCloud<PointType>::Ptr cloudReturn(new pcl::PointCloud<PointType>());
    for (auto &it : vec_)
    {
        cloudReturn->points.emplace_back(cloud_->points[it]);
    }
    return cloudReturn;
}

/**
 * @brief 识别潜在的动态物体(PD)
 * @param[in,out] ssc 场景分割结果
 * @details 
 * 该函数通过以下步骤识别潜在的动态物体:
 * 1. 遍历每个聚类
 * 2. 获取聚类中所有体素的点云索引
 * 3. 根据索引提取体素点云
 * 4. 计算体素点云的高度范围
 * 5. 如果高度范围满足以下条件,则认为是潜在动态物体:
 *    - 最低点高度小于等于-(传感器高度-0.2)
 *    - 最高点高度(加上传感器高度)小于等于预设的动态物体高度阈值
 * 6. 将满足条件的聚类ID添加到PD_cluster中
 * 7. 最后输出识别到的动态物体数量
 */
void TGRS::recognizePD(SSC &ssc)
{
    for (auto &it : ssc.cluster_vox)
    {
        std::vector<int> voxIdx;
        for (auto &id : it.second)
        {
            voxIdx.emplace_back(ssc.hash_cloud[id].ptVoxIdx);
        }
        pcl::PointCloud<PointType>::Ptr voxels(new pcl::PointCloud<PointType>());
        *voxels += *getCloudByVec(voxIdx, ssc.cloud_vox);
        std::pair<PointType, PointType> heightPair = getBoundingBoxOfCloud(voxels);
        if ((heightPair.first.z <= -(SENSOR_HEIGHT - 0.2)) && ((heightPair.second.z + SENSOR_HEIGHT) <= PD_HEIGHT))
        {
            ssc.PD_cluster.emplace_back(it.first);
        }
    }
    std::cout << "There are " << ssc.PD_cluster.size() << " PD objects." << std::endl;
}

/**
 * @brief 跟踪潜在动态物体(PD)
 * @param[in] ssc_pre 前一帧场景分割结果
 * @param[in] pose_pre 前一帧位姿
 * @param[in,out] ssc_next 当前帧场景分割结果
 * @param[in] pose_next 当前帧位姿
 * @details
 * 该函数通过以下步骤跟踪潜在动态物体:
 * 1. 获取前后帧的体素点云
 * 2. 将当前帧点云变换到前一帧坐标系下
 * 3. 对每个PD进行投影跟踪:
 *    - 获取PD所在体素的邻域体素
 *    - 计算与前一帧的重叠率
 *    - 根据重叠率判断是否为高动态物体(HD)
 *    - 重叠率低于阈值的PD标记为HD,否则标记为AS
 * 4. 将AS点云添加到非动态点云中
 * 5. 输出PD和HD的数量
 */
void TGRS::trackPD(SSC &ssc_pre, PointTypePose *pose_pre, SSC &ssc_next, PointTypePose *pose_next)
{
    // Step 1: 获取体素点云
    pcl::PointCloud<PointType>::Ptr voxCloud_pre(new pcl::PointCloud<PointType>());
    *voxCloud_pre += *ssc_pre.cloud_vox;
    pcl::PointCloud<PointType>::Ptr voxCloud_next(new pcl::PointCloud<PointType>());
    *voxCloud_next += *ssc_next.cloud_vox;
    // std::cout << "pre vox cloud size: " << voxCloud_pre->points.size();
    // std::cout << "next vox cloud size: " << voxCloud_next->points.size();

    // Step 2: 点云坐标变换
    Eigen::Affine3f trans_pre = pcl::getTransformation(pose_pre->x, pose_pre->y, pose_pre->z, pose_pre->roll, pose_pre->pitch, pose_pre->yaw);
    Eigen::Affine3f trans_next = pcl::getTransformation(pose_next->x, pose_next->y, pose_next->z, pose_next->roll, pose_next->pitch, pose_next->yaw);
    Eigen::Affine3f trans_n2p = trans_pre.inverse() * trans_next;
    pcl::PointCloud<PointType>::Ptr voxCloud_nextTrans(new pcl::PointCloud<PointType>());
    transformPointCloud(voxCloud_next, trans_n2p, voxCloud_nextTrans);
    // pcl::io::savePCDFile("/home/yixin-f/fast-lio2/src/data_dy/vox_next.pcd", *voxCloud_next);
    // pcl::io::savePCDFile("/home/yixin-f/fast-lio2/src/data_dy/vox_pre.pcd", *voxCloud_pre);

    // Step 3: PD投影跟踪
    for (auto &pd : ssc_next.PD_cluster)
    {
        std::vector<int> projIdx;
        for (auto &voxIdx : ssc_next.cluster_vox[pd])
        {
            PointType voxPt = voxCloud_nextTrans->points[ssc_next.hash_cloud[voxIdx].ptVoxIdx];
            float dis = pointDistance2d(voxPt);
            float angle = getPolarAngle(voxPt);
            float azimuth = getAzimuth(voxPt);
            int range_idx = std::ceil((dis - MIN_DIS) / RANGE_RES) - 1;
            int sector_idx = std::ceil((angle - MIN_ANGLE) / SECTOR_RES) - 1;
            int azimuth_idx = std::ceil((azimuth - MIN_AZIMUTH) / AZIMUTH_RES) - 1;
            int voxel_idx = azimuth_idx * RANGE_NUM * SECTOR_NUM + range_idx * SECTOR_NUM + sector_idx;

            std::vector<int> neighbor = findVoxelNeighbors(range_idx, sector_idx, azimuth_idx, 1);

            // 选择1: 寻找邻域
            addVec(projIdx, neighbor);

            // // 选择2: 直接投影
            // projIdx.emplace_back(voxel_idx);
        }
        sampleVec(projIdx);
        std::cout << "cur pd Idx: " << ssc_next.cluster_vox[pd].size() << std::endl;

        // Step 4: HD检测
        int all = projIdx.size();
        int success = 0;
        std::unordered_map<int, Voxel>::iterator it_find;
        for (auto &proj : projIdx)
        {
            it_find = ssc_pre.hash_cloud.find(proj);
            if (it_find != ssc_pre.hash_cloud.end())
            {
                success++;
            }
        }
        float overlapRatio = (float)success / (float)all;
        std::cout << "name: " << pd << " success: " << success << " all: " << all << " overlap ratio: " << overlapRatio << std::endl;
        if (overlapRatio <= HD_RATIO)
        {
            ssc_next.HD_cluster.emplace_back(pd); // PD转为HD
        }
        else
        {
            ssc_next.AS_cluster.emplace_back(pd); // PD转为AS
        }
    }

    std::vector<int> AS_ptIdx;
    for (auto &as : ssc_next.AS_cluster)
    {
        addVec(AS_ptIdx, ssc_next.hash_cloud[as].ptIdx);
    }
    *ssc_next.cloud_nd += *getCloudByVec(AS_ptIdx, ssc_next.cloud_use);
    *ssc_next.cloud_nd += *ssc_next.cloud_g;

    std::cout << ANSI_COLOR_GREEN << " PD num: " << ssc_next.PD_cluster.size() << "\n"
              << ANSI_COLOR_RED << " HD num: " << ssc_next.HD_cluster.size() << ANSI_COLOR_RESET << std::endl;
}

// 保存带颜色的点云
// 输入:
//   ssc: 场景分割结果
//   path: 保存路径
// 功能:
//   1. 创建一个带颜色的点云对象colorCloud
//   2. 为每个聚类生成随机RGB颜色
//   3. 遍历每个聚类:
//      - 获取聚类中所有点的索引
//      - 根据索引获取原始点云
//      - 为每个点赋予聚类对应的颜色
//   4. 设置点云高度为1(无序点云)
//   5. 设置点云宽度为点数
//   6. 输出分割后的点云大小
//   7. 保存为PCD文件
void TGRS::saveColorCloud(SSC &ssc, const std::string &path)
{
    pcl::PointCloud<pcl::PointXYZRGB>::Ptr colorCloud(new pcl::PointCloud<pcl::PointXYZRGB>());
    cv::RNG rng(12345);
    for (auto &it : ssc.cluster_vox)
    {
        int r, g, b;
        r = rng.uniform(20, 200);
        g = rng.uniform(20, 200);
        b = rng.uniform(20, 200);
        std::vector<int> ptIdx;
        for (auto &vox : it.second)
        {
            addVec(ptIdx, ssc.hash_cloud[vox].ptIdx);
        }
        pcl::PointCloud<PointType>::Ptr cloudGrab(new pcl::PointCloud<PointType>());
        *cloudGrab += *getCloudByVec(ptIdx, ssc.cloud_use);
        for (size_t i = 0; i < cloudGrab->points.size(); i++)
        {
            pcl::PointXYZRGB rgb;
            rgb.x = cloudGrab->points[i].x;
            rgb.y = cloudGrab->points[i].y;
            rgb.z = cloudGrab->points[i].z;
            rgb.r = r;
            rgb.g = g;
            rgb.b = b;
            colorCloud->points.emplace_back(rgb);
        }
    }
    colorCloud->height = 1;
    colorCloud->width = colorCloud->points.size();
    std::cout << "segmented cloud size: " << colorCloud->points.size() << std::endl;
    pcl::io::savePCDFile(path, *colorCloud);
}
