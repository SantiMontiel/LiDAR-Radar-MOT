/********************************************************************
 * 
 *  Multi-Object Tracking system based on LiDAR and Radar for 
 *  Intelligent Vehicles Applications
 * 
 *  @authors: Santiago Montiel Marín and Carlos Gómez Huélamo
 *
 *  RobeSafe Research Group
 *  Department of Electronics, University of Alcalá
 * 
 * ******************************************************************/

// -- General purpose includes
#include <iostream>
#include <vector>

// -- ROS modules includes
#include <ros/ros.h>
#include <std_msgs/ColorRGBA.h>
#include <sensor_msgs/PointCloud2.h>
#include <visualization_msgs/Marker.h>
#include <visualization_msgs/MarkerArray.h>

// -- Submodules includes
#include "../include/lidar.hpp"
#include "../include/object.hpp"
#include "../include/kdtree.hpp"

// -- ROS publishers definitions
ros::Publisher pub_LiDAR_FilteredCloud;
ros::Publisher pub_LiDAR_RoadCloud;
ros::Publisher pub_LiDAR_ObstaclesMarkers;
//ros::Publisher pub_LiDAR_BoundingBoxes;

// -- ROS subscribers definitions
ros::Subscriber sub_LiDAR_RawPointCloud;



// -- LiDAR callback
void LiDAR_CB (const sensor_msgs::PointCloud2::ConstPtr& LidarMsg) {

    // 0. CLOUD PRE-PROCESSING
    // 0.a. Auxiliar definitions for publishing objects with colours
    std_msgs::ColorRGBA red;
    red.a = 1.0, red.r = 1.0, red.g = 0.0, red.b = 0.0;
    std_msgs::ColorRGBA green;
    green.a = 1.0, green.r = 0.0, green.g = 1.0, green.b = 0.0;
    std_msgs::ColorRGBA blue;
    blue.a = 1.0, blue.r = 0.0, blue.g = 0.0, blue.g = 1.0;

    // 0.b. Transforming ROS message into PCL point cloud
    pcl::PointCloud<pcl::PointXYZ>::Ptr nonFilteredCloud (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::fromROSMsg(*LidarMsg, *nonFilteredCloud);

    // 1. CLOUD FILTERING (by XYZ and angle)
    pcl::PointCloud<pcl::PointXYZ> auxFilteredCloud = CloudFiltering(nonFilteredCloud);
    pcl::PointCloud<pcl::PointXYZ>::Ptr FilteredCloud (new pcl::PointCloud<pcl::PointXYZ>);
    *FilteredCloud = auxFilteredCloud;

    // 1.b. Publishing intermediate output 1 -> filtered cloud
    sensor_msgs::PointCloud2 msgFilteredCloud;
    pcl::toROSMsg(*FilteredCloud, msgFilteredCloud);
    msgFilteredCloud.header.frame_id = LidarMsg->header.frame_id;
    msgFilteredCloud.header.stamp = LidarMsg->header.stamp;

    pub_LiDAR_FilteredCloud.publish(msgFilteredCloud);

    // 2. PLANE SEGMENTATION (with Ransac3D algorithm)
    // std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> segCloud (std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr>);
    // pcl::PointCloud<pcl::PointXYZ>::Ptr roadCloud (new pcl::PointCloud<pcl::PointXYZ>);
    // pcl::PointCloud<pcl::PointXYZ>::Ptr obstCloud (new pcl::PointCloud<pcl::PointXYZ>);
    // segCloud = PlaneSegmentation(FilteredCloud, 100, 0.2);
    // obstCloud = *segCloud.first;
    // roadCloud = *segCloud.second;

    // // 2.b. Publishing intermediate output 2 -> segmented plane
    // sensor_msgs::PointCloud2 msgRoadCloud;
    // pcl::toROSMsg(*roadCloud, msgRoadCloud);
    // msgRoadCloud.header.frame_id = LidarMsg->header.frame_id;
    // msgRoadCloud.header.stamp = LidarMsg->header.stamp;

    // pub_LiDAR_RoadCloud.publish(msgRoadCloud);
    
    // 3. CLUSTERING EXTRACTION
    std::vector<Object> Obstacles;
    int numObstacles = 0;

    pcl::PointCloud<pcl::PointXYZ>::Ptr ObstaclesCloud (new pcl::PointCloud<pcl::PointXYZ>);
    *ObstaclesCloud = auxFilteredCloud;
    // ClusteringExtraction(ObstaclesCloud, 1, 10, 400, &Obstacles, &numObstacles);
    cluster_filter(ObstaclesCloud, 1.4, 20, 150, &Obstacles, &numObstacles);

    visualization_msgs::MarkerArray ObstaclesMarkerArray;

    // 3.b. Publishing clusters
    for (int i = 0; i < Obstacles.size(); i++) {

		visualization_msgs::Marker ObstaclesMarker;

        ObstaclesMarker.header.frame_id = "/ego_vehicle/lidar/lidar1";                          // map == global coordinates. Base_link == local coordinates
        // ObstaclesMarker.header.frame_id = "lidar_vlp16";
        ObstaclesMarker.header.stamp = LidarMsg->header.stamp;
        ObstaclesMarker.ns = "liramot_lidar_bounding_boxes";
        ObstaclesMarker.action = visualization_msgs::Marker::ADD;
        ObstaclesMarker.type = visualization_msgs::Marker::CUBE;
        ObstaclesMarker.lifetime = ros::Duration(0.2);
        ObstaclesMarker.id = Obstacles[i].object_id;

        ObstaclesMarker.points.clear();

        ObstaclesMarker.color = red;
        ObstaclesMarker.scale.x = Obstacles[i].d;
        ObstaclesMarker.scale.y = Obstacles[i].w;
        ObstaclesMarker.scale.z = Obstacles[i].h;
        ObstaclesMarker.lifetime = ros::Duration(0.0);

        // std::cout << "XYZ: " << Obstacles[i].centroid_x << " " << Obstacles[i].centroid_y << " " << Obstacles[i].centroid_z << std::endl;
        ObstaclesMarker.pose.position.x = Obstacles[i].centroid_x;
        ObstaclesMarker.pose.position.y = Obstacles[i].centroid_y;
        ObstaclesMarker.pose.position.z = Obstacles[i].centroid_z;

        ObstaclesMarkerArray.markers.push_back(ObstaclesMarker);
	}

    pub_LiDAR_ObstaclesMarkers.publish(ObstaclesMarkerArray);

    // 3.c. Publishing bounding boxes
    
}

// -- Main function
int main (int argc, char** argv) {

    // -- ROS initialization
    ros::init(argc, argv, "lidar_radar_mot_node");
    ros::NodeHandle nh;

    // -- ROS publishers initialization
    pub_LiDAR_FilteredCloud = nh.advertise<sensor_msgs::PointCloud2>("/t4ac/perception/detection/liramot/filtered_cloud", 1, true);
    pub_LiDAR_RoadCloud = nh.advertise<sensor_msgs::PointCloud2>("/t4ac/perception/detection/liramot/road_cloud", 1, true);
    pub_LiDAR_ObstaclesMarkers = nh.advertise<visualization_msgs::MarkerArray>("/t4ac/perception/detection/liramot/obstacle_markers", 1, true);

    // -- ROS subscribers initialization
    sub_LiDAR_RawPointCloud = nh.subscribe("/carla/ego_vehicle/lidar/lidar1/point_cloud", 1, &LiDAR_CB);
    // sub_LiDAR_RawPointCloud = nh.subscribe("/velodyne_points", 1, &LiDAR_CB);

    // -- ROS spin
    ros::spin();

}

