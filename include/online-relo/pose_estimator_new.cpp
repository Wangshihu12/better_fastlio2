#include "pose_estimator_new.h"

pose_estimator::pose_estimator(std::string priorPath){
    nh.param<std::string>("common/rootDir", rootDir, " ");
    nh.param<std::string>("location/pointCloudTopic", pointCloudTopic, "points_raw");
    nh.param<std::string>("location/imuTopic", imuTopic, "/imu");
    nh.param<double>("mapping/gyr_cov", gyr_cov, 0.1);
    nh.param<double>("mapping/acc_cov", acc_cov, 0.1);
    nh.param<double>("mapping/b_gyr_cov", b_gyr_cov, 0.0001);
    nh.param<double>("mapping/b_acc_cov", b_acc_cov, 0.0001);

    // if(p_pre->lidar_type == LIVOX){
    //     subCloud = nh.subscribe(pointCloudTopic, 500, &pose_estimator::livox_pcl_cbk, this);
    // }
    // else{
    //     subCloud = nh.subscribe(pointCloudTopic, 500, &pose_estimator::standard_pcl_cbk, this);
    // }
    
    // subIMU = nh.subscribe(imuTopic, 20000, &pose_estimator::imuCBK, this);
    // subPose = nh.subscribe<geometry_msgs::PoseWithCovarianceStamped>("/initialpose", 1, &pose_estimator::poseCBK, this);  // interaction in rviz

    pubOdomAftMapped = nh.advertise<nav_msgs::Odometry>("/Odometry", 100000);

    std::cout << ANSI_COLOR_GREEN << "rostopic is ok" << ANSI_COLOR_RESET << std::endl;

    MultiSession::Session priorKnown_(1, "priorMap", priorPath, true);  // prior knowledge
    priorKnown = &priorKnown_;  // copy
    downSizeFilterSurf.setLeafSize(0.1, 0.1, 0.1);

    std::cout << ANSI_COLOR_GREEN << "prior knowledge is loaded" << ANSI_COLOR_RESET << std::endl;

    path.header.stamp = ros::Time::now();
    path.header.frame_id = "camera_init";

    Lidar_T_wrt_IMU << VEC_FROM_ARRAY(extrinT);
    Lidar_R_wrt_IMU << MAT_FROM_ARRAY(extrinR);
    p_imu->set_extrinsic(Lidar_T_wrt_IMU, Lidar_R_wrt_IMU);        
    p_imu->set_gyr_cov(V3D(gyr_cov, gyr_cov, gyr_cov));           
    p_imu->set_acc_cov(V3D(acc_cov, acc_cov, acc_cov));            
    p_imu->set_gyr_bias_cov(V3D(b_gyr_cov, b_gyr_cov, b_gyr_cov)); 
    p_imu->set_acc_bias_cov(V3D(b_acc_cov, b_acc_cov, b_acc_cov)); 

    // esekf
    fill(epsi, epsi + 23, 0.001);
    kf.init_dyn_share(get_f, df_dx, df_dw, h_share_model, NUM_MAX_ITERATIONS, epsi);
}

void pose_estimator::run(){
    ros::Rate rate(5000);
    bool status = ros::ok();
    while (status){
        if (flg_exit)
            break;
        ros::spinOnce();

        if (sync_packages(Measures)){
            if (flg_first_scan){
                first_lidar_time = Measures.lidar_beg_time;
                p_imu->first_lidar_time = first_lidar_time;
                flg_first_scan = false;
                continue;
            }
        }

        double t0, t1, t2, t3, t4, t5, match_start, solve_start, svd_time;
        double match_time = 0, kdtree_search_time = 0.0, solve_time = 0, solve_const_H_time = 0;
        t0 = omp_get_wtime();

        p_imu->Process(Measures, kf, feats_undistort);
        state_point = kf.get_x();  
        pos_lid = state_point.pos + state_point.rot * state_point.offset_T_L_I;

        if (feats_undistort->empty() || (feats_undistort == NULL)){
            ROS_WARN("No point, skip this scan!\n");
            continue;
        }

        flg_EKF_inited = (Measures.lidar_beg_time - first_lidar_time) < INIT_TIME ? false : true;

        lasermap_fov_segment();

        downSizeFilterSurf.setInputCloud(feats_undistort);
        downSizeFilterSurf.filter(*feats_down_body); 

        t1 = omp_get_wtime();  
        feats_down_size = feats_down_body->points.size();

        if (ikdtree.Root_Node == nullptr)
        {
            if (feats_down_size > 5)
            {
                ikdtree.set_downsample_param(filter_size_map_min); 
                feats_down_world->resize(feats_down_size);         
                for (int i = 0; i < feats_down_size; i++)
                {
                    pointBodyToWorld(&(feats_down_body->points[i]), &(feats_down_world->points[i]));
                }
                ikdtree.Build(feats_down_world->points); 
            }
            continue;      
        }

        int featsFromMapNum = ikdtree.validnum();
        kdtree_size_st = ikdtree.size();

        /*** ICP and iterated Kalman filter update ***/
        if (feats_down_size < 5)
        {
            ROS_WARN("No point, skip this scan!\n");
            continue;
        }

        normvec->resize(feats_down_size);
        feats_down_world->resize(feats_down_size);

        // if (1) 
        // {
        //     PointVector().swap(ikdtree.PCL_Storage);                     
        //     ikdtree.flatten(ikdtree.Root_Node, ikdtree.PCL_Storage, NOT_RECORD); 
        //     featsFromMap->clear();
        //     featsFromMap->points = ikdtree.PCL_Storage;
        //     publish_map(pubIkdTree);
        // }

        pointSearchInd_surf.resize(feats_down_size);
        Nearest_Points.resize(feats_down_size);

        int rematch_num = 0;
        bool nearest_search_en = true;

        t2 = omp_get_wtime();

        double t_update_start = omp_get_wtime();
        double solve_H_time = 0;
        kf.update_iterated_dyn_share_modified(LASER_POINT_COV, solve_H_time); 
        state_point = kf.get_x();
        euler_cur = SO3ToEuler(state_point.rot);
        pos_lid = state_point.pos + state_point.rot * state_point.offset_T_L_I; // TODO: give initial pose here
        geoQuat.x = state_point.rot.coeffs()[0];                   
        geoQuat.y = state_point.rot.coeffs()[1];
        geoQuat.z = state_point.rot.coeffs()[2];
        geoQuat.w = state_point.rot.coeffs()[3];

        double t_update_end = omp_get_wtime();

        getCurPose(state_point);

        bool useRelo = easyToRelo();

        if(useRelo){
            std::cout << "relo mode" << std::endl;
            recontructIKdTree();  // TODO: get current measurement in prior map
        }
        else{
            std::cout << "lio mode" << std::endl;
            map_incremental(); // TODO: mapping
        }

        publish_odometry(pubOdomAftMapped);

        t3 = omp_get_wtime();

        publish_path(pubPath);

        publish_frame_world(pubCurCloud); 

        status = ros::ok();
        rate.sleep();
    }
}

bool pose_estimator::sync_packages(MeasureGroup &meas){
    if (lidar_buffer.empty() || imu_buffer.empty())
    {
        return false;
    }

    if (!lidar_pushed)
    {
        meas.lidar = lidar_buffer.front();        
        meas.lidar_beg_time = time_buffer.front(); 

        if (meas.lidar->points.size() <= 1)
        {
            lidar_end_time = meas.lidar_beg_time + lidar_mean_scantime;
            ROS_WARN("Too few input point cloud!\n");
        }
        else if (meas.lidar->points.back().curvature / double(1000) < 0.5 * lidar_mean_scantime)
        {
            lidar_end_time = meas.lidar_beg_time + lidar_mean_scantime;
        }
        else
        {
            scan_num++; 
            lidar_end_time = meas.lidar_beg_time + meas.lidar->points.back().curvature / double(1000);
            lidar_mean_scantime += (meas.lidar->points.back().curvature / double(1000) - lidar_mean_scantime) / scan_num;
        }

        meas.lidar_end_time = lidar_end_time; 

        lidar_pushed = true; 
    }
    if (last_timestamp_imu < lidar_end_time)
    {
        return false;
    }

    double imu_time = imu_buffer.front()->header.stamp.toSec(); 
    meas.imu.clear();
    while ((!imu_buffer.empty()) && (imu_time < lidar_end_time))
    {
        imu_time = imu_buffer.front()->header.stamp.toSec(); 
        if (imu_time > lidar_end_time)
            break;
        meas.imu.push_back(imu_buffer.front()); 
        imu_buffer.pop_front();
    }

    lidar_buffer.pop_front(); 
    time_buffer.pop_front();  
    lidar_pushed = false;    
    return true;
}

void pose_estimator::cloudCBK(){
    // TODO: unified
}

void pose_estimator::livox_pcl_cbk(const livox_ros_driver::CustomMsg::ConstPtr &msg){
    mtx_buffer.lock(); 
    double preprocess_start_time = omp_get_wtime();
    scan_count++; 
    if (msg->header.stamp.toSec() < last_timestamp_lidar)
    {
        ROS_ERROR("lidar loop back, clear buffer");
        lidar_buffer.clear();
    }
    last_timestamp_lidar = msg->header.stamp.toSec();
    if (!time_sync_en && abs(last_timestamp_imu - last_timestamp_lidar) > 10.0 && !imu_buffer.empty() && !lidar_buffer.empty())
    {
        printf("IMU and LiDAR not Synced, IMU time: %lf, lidar header time: %lf \n", last_timestamp_imu, last_timestamp_lidar);
    }
    if (time_sync_en && !timediff_set_flg && abs(last_timestamp_lidar - last_timestamp_imu) > 1 && !imu_buffer.empty())
    {
        timediff_set_flg = true;
        timediff_lidar_wrt_imu = last_timestamp_lidar + 0.1 - last_timestamp_imu; 
        printf("Self sync IMU and LiDAR, time diff is %.10lf \n", timediff_lidar_wrt_imu);
    }
    pcl::PointCloud<PointType>::Ptr ptr(new pcl::PointCloud<PointType>());
    p_pre->process(msg, ptr);              
    lidar_buffer.push_back(ptr);                 
    time_buffer.push_back(last_timestamp_lidar); 

    mtx_buffer.unlock();    
    sig_buffer.notify_all();
}

void pose_estimator::standard_pcl_cbk(const sensor_msgs::PointCloud2::ConstPtr &msg){
    mtx_buffer.lock();
    scan_count++;
    double preprocess_start_time = omp_get_wtime();
    if (msg->header.stamp.toSec() < last_timestamp_lidar)
    {
        ROS_ERROR("lidar loop back, clear buffer");
        lidar_buffer.clear();
    }
    pcl::PointCloud<PointType>::Ptr ptr(new pcl::PointCloud<PointType>());
    p_pre->process(msg, ptr);                       
    lidar_buffer.push_back(ptr);                      
    time_buffer.push_back(msg->header.stamp.toSec()); 
    last_timestamp_lidar = msg->header.stamp.toSec(); 

    mtx_buffer.unlock();     
    sig_buffer.notify_all(); 
}

void pose_estimator::imuCBK(const sensor_msgs::Imu::ConstPtr &msg_in){
    publish_count++;
    sensor_msgs::Imu::Ptr msg(new sensor_msgs::Imu(*msg_in));
    if (abs(timediff_lidar_wrt_imu) > 0.1 && time_sync_en)
    {
        msg->header.stamp = ros::Time().fromSec(timediff_lidar_wrt_imu + msg_in->header.stamp.toSec());
    }
    double timestamp = msg->header.stamp.toSec(); // IMU时间戳

    mtx_buffer.lock(); 
    if (timestamp < last_timestamp_imu)
    {
        ROS_WARN("imu loop back, clear buffer");
        imu_buffer.clear();
    }

    last_timestamp_imu = timestamp; 

    imu_buffer.push_back(msg); 
    mtx_buffer.unlock();       
    sig_buffer.notify_all();   
}

void pose_estimator::poseCBK(const geometry_msgs::PoseWithCovarianceStampedConstPtr &msg){
    if(initpose_flag){
        return ;
    }
    PointType p;
    initpose.x = msg->pose.pose.position.x;
    initpose.y = msg->pose.pose.position.y;
    initpose.z = msg->pose.pose.position.z;
    double roll, pitch, yaw;
    tf::Quaternion q;
    tf::quaternionMsgToTF(msg->pose.pose.orientation, q);
    tf::Matrix3x3(q).getRPY(roll, pitch, yaw);
    initpose.roll = roll;
    initpose.pitch = pitch;
    initpose.yaw = yaw;
    std::cout << ANSI_COLOR_RED << "Get initial pose: " << initpose.x << " " << initpose.y << " "
              << initpose.z << " " << roll << " " << pitch << " " << yaw << ANSI_COLOR_RESET << std::endl;
    // initpose_flag = true;
}

void pose_estimator::lasermap_fov_segment(){
    cub_needrm.clear(); 
    kdtree_delete_counter = 0;
    kdtree_delete_time = 0.0;
    V3D pos_LiD = pos_lid; 

    if (!Localmap_Initialized)
    {
        for (int i = 0; i < 3; i++)
        {
            LocalMap_Points.vertex_min[i] = pos_LiD(i) - cube_len / 2.0; 
            LocalMap_Points.vertex_max[i] = pos_LiD(i) + cube_len / 2.0;// 边界距离当前位置100米
        }
        Localmap_Initialized = true;
        return;
    }

    float dist_to_map_edge[3][2];
    bool need_move = false;
    for (int i = 0; i < 3; i++)
    {
        dist_to_map_edge[i][0] = fabs(pos_LiD(i) - LocalMap_Points.vertex_min[i]);
        dist_to_map_edge[i][1] = fabs(pos_LiD(i) - LocalMap_Points.vertex_max[i]);
        if (dist_to_map_edge[i][0] <= MOV_THRESHOLD * DET_RANGE || dist_to_map_edge[i][1] <= MOV_THRESHOLD * DET_RANGE)
            need_move = true;
    }

    if (!need_move)
        return;
    
    BoxPointType New_LocalMap_Points, tmp_boxpoints;
    New_LocalMap_Points = LocalMap_Points;
    float mov_dist = max((cube_len - 2.0 * MOV_THRESHOLD * DET_RANGE) * 0.5 * 0.9, double(DET_RANGE * (MOV_THRESHOLD - 1)));
    for (int i = 0; i < 3; i++)
    {
        tmp_boxpoints = LocalMap_Points;
        if (dist_to_map_edge[i][0] <= MOV_THRESHOLD * DET_RANGE)
        {
            New_LocalMap_Points.vertex_max[i] -= mov_dist;
            New_LocalMap_Points.vertex_min[i] -= mov_dist;
            tmp_boxpoints.vertex_min[i] = LocalMap_Points.vertex_max[i] - mov_dist;
            cub_needrm.push_back(tmp_boxpoints); 
        }
        else if (dist_to_map_edge[i][1] <= MOV_THRESHOLD * DET_RANGE)
        {
            New_LocalMap_Points.vertex_max[i] += mov_dist;
            New_LocalMap_Points.vertex_min[i] += mov_dist;
            tmp_boxpoints.vertex_max[i] = LocalMap_Points.vertex_min[i] + mov_dist;
            cub_needrm.push_back(tmp_boxpoints);
        }
    }
    LocalMap_Points = New_LocalMap_Points;

    double delete_begin = omp_get_wtime();
    if (cub_needrm.size() > 0)
        kdtree_delete_counter = ikdtree.Delete_Point_Boxes(cub_needrm);
    kdtree_delete_time = omp_get_wtime() - delete_begin;

}

void pose_estimator::publish_map(const ros::Publisher &pubLaserCloudMap){
    sensor_msgs::PointCloud2 laserCloudMap;
    pcl::toROSMsg(*featsFromMap, laserCloudMap);
    laserCloudMap.header.stamp = ros::Time().fromSec(lidar_end_time);
    laserCloudMap.header.frame_id = "camera_init";
    pubLaserCloudMap.publish(laserCloudMap);
}

void pose_estimator::getCurPose(state_ikfom cur_state){
    Eigen::Vector3d eulerAngle = cur_state.rot.matrix().eulerAngles(2, 1, 0); 
    transformTobeMapped[0] = eulerAngle(2);    
    transformTobeMapped[1] = eulerAngle(1);    
    transformTobeMapped[2] = eulerAngle(0);    
    transformTobeMapped[3] = cur_state.pos(0); 
    transformTobeMapped[4] = cur_state.pos(1); 
    transformTobeMapped[5] = cur_state.pos(2); 
}

void pose_estimator::recontructIKdTree(){
    pcl::PointCloud<PointType>::Ptr subMapKeyFrames(new pcl::PointCloud<PointType>());
    pcl::PointCloud<PointType>::Ptr subMapKeyFramesDS(new pcl::PointCloud<PointType>());

    for(int i = 0; i < pointSearchIndGlobalMap.size(); i++){
        *subMapKeyFrames += *transformPointCloud(priorKnown->cloudKeyFrames[pointSearchIndGlobalMap[i]].all_cloud, &priorKnown->cloudKeyPoses6D->points[pointSearchIndGlobalMap[i]]);
    }

    downSizeFilterSurf.setInputCloud(subMapKeyFrames);
    downSizeFilterSurf.filter(*subMapKeyFramesDS);
    ikdtree.reconstruct(subMapKeyFramesDS->points);
    int featsFromMapNum = ikdtree.validnum();
    kdtree_size_st = ikdtree.size();

    std::cout << "current measurement in prior map: " << " featsFromMapNum  =  " << featsFromMapNum << "\t"
                  << " kdtree_size_st   =  " << kdtree_size_st << std::endl;

    featsFromMap->clear();  // publish kdtree
    featsFromMap->points = subMapKeyFramesDS->points;
    publish_map(pubIkdTree);
}
 
void pose_estimator::pointBodyToWorld(PointType const *const pi, PointType *const po)
{
    V3D p_body(pi->x, pi->y, pi->z);
    V3D p_global(state_point.rot * (state_point.offset_R_L_I * p_body + state_point.offset_T_L_I) + state_point.pos);

    po->x = p_global(0);
    po->y = p_global(1);
    po->z = p_global(2);
    po->intensity = pi->intensity;
}

void pose_estimator::publish_odometry(const ros::Publisher &pubOdomAftMapped)
{
    odomAftMapped.header.frame_id = "camera_init";
    odomAftMapped.child_frame_id = "body";
    odomAftMapped.header.stamp = ros::Time().fromSec(lidar_end_time);
    set_posestamp(odomAftMapped.pose);
    pubOdomAftMapped.publish(odomAftMapped);
    auto P = kf.get_P();
    for (int i = 0; i < 6; i++)
    {
        int k = i < 3 ? i + 3 : i - 3;
        odomAftMapped.pose.covariance[i * 6 + 0] = P(k, 3);
        odomAftMapped.pose.covariance[i * 6 + 1] = P(k, 4);
        odomAftMapped.pose.covariance[i * 6 + 2] = P(k, 5);
        odomAftMapped.pose.covariance[i * 6 + 3] = P(k, 0);
        odomAftMapped.pose.covariance[i * 6 + 4] = P(k, 1);
        odomAftMapped.pose.covariance[i * 6 + 5] = P(k, 2);
    }

    static tf::TransformBroadcaster br;
    tf::Transform transform;
    tf::Quaternion q;
    transform.setOrigin(tf::Vector3(odomAftMapped.pose.pose.position.x,
                                    odomAftMapped.pose.pose.position.y,
                                    odomAftMapped.pose.pose.position.z));
    q.setW(odomAftMapped.pose.pose.orientation.w);
    q.setX(odomAftMapped.pose.pose.orientation.x);
    q.setY(odomAftMapped.pose.pose.orientation.y);
    q.setZ(odomAftMapped.pose.pose.orientation.z);
    transform.setRotation(q);
    br.sendTransform(tf::StampedTransform(transform, odomAftMapped.header.stamp, "camera_init", "body"));
}

void pose_estimator::map_incremental()
{
    PointVector PointToAdd;           
    PointVector PointNoNeedDownsample; 
    PointToAdd.reserve(feats_down_size);
    PointNoNeedDownsample.reserve(feats_down_size);

    for (int i = 0; i < feats_down_size; i++)
    {
        pointBodyToWorld(&(feats_down_body->points[i]), &(feats_down_world->points[i]));
        if (!Nearest_Points[i].empty() && flg_EKF_inited)
        {
            const PointVector &points_near = Nearest_Points[i];
            bool need_add = true;
            BoxPointType Box_of_Point;
            PointType downsample_result, mid_point;
            mid_point.x = floor(feats_down_world->points[i].x / filter_size_map_min) * filter_size_map_min + 0.5 * filter_size_map_min;
            mid_point.y = floor(feats_down_world->points[i].y / filter_size_map_min) * filter_size_map_min + 0.5 * filter_size_map_min;
            mid_point.z = floor(feats_down_world->points[i].z / filter_size_map_min) * filter_size_map_min + 0.5 * filter_size_map_min;
            float dist = calc_dist(feats_down_world->points[i], mid_point); 

            if (fabs(points_near[0].x - mid_point.x) > 0.5 * filter_size_map_min && fabs(points_near[0].y - mid_point.y) > 0.5 * filter_size_map_min && fabs(points_near[0].z - mid_point.z) > 0.5 * filter_size_map_min)
            {
                PointNoNeedDownsample.push_back(feats_down_world->points[i]); 
                continue;
            }

            for (int readd_i = 0; readd_i < NUM_MATCH_POINTS; readd_i++)
            {
                if (points_near.size() < NUM_MATCH_POINTS)
                    break;
                if (calc_dist(points_near[readd_i], mid_point) < dist) 
                {
                    need_add = false;
                    break;
                }
            }
            if (need_add)
                PointToAdd.push_back(feats_down_world->points[i]);
        }
        else
        {
            PointToAdd.push_back(feats_down_world->points[i]);
        }
    }

    double st_time = omp_get_wtime();
    add_point_size = ikdtree.Add_Points(PointToAdd, true); 
    ikdtree.Add_Points(PointNoNeedDownsample, false);      
    add_point_size = PointToAdd.size() + PointNoNeedDownsample.size();
    kdtree_incremental_time = omp_get_wtime() - st_time;
    publish_map(pubIkdTree);

    std::cout << "map incremental: " << " add_point_size  =  " << add_point_size << "\t"
                  << " kdtree_incremental_time   =  " << kdtree_incremental_time << std::endl;
}

void pose_estimator::publish_path(const ros::Publisher pubPath){
    set_posestamp(msg_body_pose);
    msg_body_pose.header.stamp = ros::Time().fromSec(lidar_end_time);
    msg_body_pose.header.frame_id = "camera_init";

    path.poses.push_back(msg_body_pose);
    pubPath.publish(path);
}

void pose_estimator::publish_frame_world(const ros::Publisher &pubLaserCloudFull){
    pcl::PointCloud<PointType>::Ptr laserCloudFullRes(feats_down_body);
    int size = laserCloudFullRes->points.size();  
    pcl::PointCloud<PointType>::Ptr laserCloudWorld(new pcl::PointCloud<PointType>(size, 1));
    for (int i = 0; i < size; i++)
    {
        pointBodyToWorld(&laserCloudFullRes->points[i], &laserCloudWorld->points[i]);
    }
    sensor_msgs::PointCloud2 laserCloudmsg;
    pcl::toROSMsg(*laserCloudWorld, laserCloudmsg);
    laserCloudmsg.header.stamp = ros::Time().fromSec(lidar_end_time);
    laserCloudmsg.header.frame_id = "camera_init";
    pubLaserCloudFull.publish(laserCloudmsg);
}

bool pose_estimator::easyToRelo(){
    pointSearchIndGlobalMap.clear();
    pointSearchSqDisGlobalMap.clear();

    kdtreeGlobalMapPoses->setInputCloud(priorKnown->cloudKeyPoses3D);  // find nearest poses in prior map
    PointType curPose;
    curPose.x = transformTobeMapped[3];
    curPose.y = transformTobeMapped[4];
    curPose.z = transformTobeMapped[5];
    kdtreeGlobalMapPoses->radiusSearch(curPose, 5.0, pointSearchIndGlobalMap, pointSearchSqDisGlobalMap, 0);  // TODO: add the radius into yaml

    if(pointSearchIndGlobalMap.size() >= 2){  // TODO: add the radius into yaml
        return true;
    }
    else{
        return false;
    }
}

void pose_estimator::getInitPose(){
    TicToc time_count;

    Eigen::MatrixXd curSC = priorKnown->scManager.makeScancontext(*feats_down_body);  // FIXME: just use a single scan !! please stay static before get accuracy pose
    Eigen::MatrixXd ringkey = priorKnown->scManager.makeRingkeyFromScancontext(curSC);
    Eigen::MatrixXd sectorkey = priorKnown->scManager.makeSectorkeyFromScancontext(curSC);
    std::vector<float> polarcontext_invkey_vec = ScanContext::eig2stdvec(ringkey);

    auto detectResult = priorKnown->scManager.detectLoopClosureIDBetweenSession(polarcontext_invkey_vec, curSC);  // idx & yaw_diff
    
    int reloIdx = detectResult.first;
    float yawDiff = detectResult.second;

    PointTypePose diffPose;
    diffPose.x = 0;
    diffPose.y = 0;
    diffPose.z = 0;
    diffPose.roll = 0;
    diffPose.pitch = 0;
    diffPose.yaw = yawDiff;

    pcl::PointCloud<PointType>::Ptr feats_down_body_diff(new pcl::PointCloud<PointType>());
    *feats_down_body_diff += *transformPointCloud(feats_down_body, &diffPose);

    // TODO: teaser ++ to get the precise initial pose in global map and reset it in the ieskf !!
}
