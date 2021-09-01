#include "../include/lidar.hpp"
#include "../include/kdtree.hpp"
#include "../include/object.hpp"

/****************************************************
 *  XYZ and angle cloud filtering
 * 
 *  I- Non filetered cloud
 *  O- Filtered cloud by XYZ and angle
 * **************************************************/

pcl::PointCloud<pcl::PointXYZ> CloudFiltering (pcl::PointCloud<pcl::PointXYZ>::Ptr nonFilteredCloud) {

    pcl::PointCloud<pcl::PointXYZ> auxFilteredCloud;
    
    // -- XYZ filtering
    for (int i = 0; i < nonFilteredCloud->points.size(); i++) {

        pcl::PointXYZ Point = nonFilteredCloud->points[i];

        // If the point is above the sidewalk -> push back
        if (Point.z > -1.90) {
            auxFilteredCloud.points.push_back(Point);
        }
    }

    // -- Angle filtering
    /*
    pcl::PointCloud<pcl::PointXYZ> FilteredCloud;
    float FieldOfView = 80;
    double FovMin = -(FieldOfView/2)*(3.14/180);;
    double FovMax = (FieldOfView/2)*(3.14/180);

    for (int i = 0; i < auxFilteredCloud.points.size(); i++) {

        pcl::PointXYZ Point = auxFilteredCloud.points[i];
        double PointAngle = atan2(Point.y, Point.x);

        // If the point is between Fov/2 and -Fov/2
        if (PointAngle < FovMax && PointAngle > FovMin) {
            FilteredCloud.points.push_back(Point);
        }

    }*/

    return auxFilteredCloud;
    // return FilteredCloud;
}

/****************************************************
 *  Plane segmentation with RANSAC-3D algorithm
 *  (needs subsequent cloud separation)
 * 
 *  I- Cloud, MaxIterations, Threshold
 *  O- Pair of clouds (plane and obstacles)
 * **************************************************/
   std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> segCloud (std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr>);

std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> PlaneSegmentation (pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud, int MaxIterations, float Threshold) {

    auto StartTime = std::chrono::steady_clock::now();

    std::unordered_set<int> InlierPoints;
    while (MaxIterations--) {

        // 1. Choose three random points
        std::unordered_set<int> TargetPoints;
        while(TargetPoints.size() < 3) {
            TargetPoints.insert(rand()%Cloud->size());
        }

        // 2. Parameters required for the plane equation
        //      ax + by + cz + d = 0 
        // 2.a. Points characterization
        float x1, y1, z1, x2, y2, z2, x3, y3, z3;
        auto itr = TargetPoints.begin();
        x1 = Cloud->points[*itr].x;
        y1 = Cloud->points[*itr].y;
        z1 = Cloud->points[*itr].z;
        itr++;
        x2 = Cloud->points[*itr].x;
        y2 = Cloud->points[*itr].y;
        z2 = Cloud->points[*itr].z;
        itr++;
        x3 = Cloud->points[*itr].x;
        y3 = Cloud->points[*itr].y;
        z3 = Cloud->points[*itr].z;
        
        // 2.b. Plane characterization
        float a, b, c, d, den;
        a = ((y2 - y1) * (z3 - z1) - (z2 - z1) * (y3 - y2));
        b = ((z2 - z1) * (x3 - x1) - (x2 - x1) * (z3 - z2));
        c = ((x2 - x1) * (y3 - y1) - (y2 - y1) * (x3 - x2));
        d = - (a * x1 + b * y1 + c * z1);
        den = sqrt(a * a + b * b + c * c);

        // 3. For all the points in the cloud, estimate distance to
        //    the plane
        for (int i = 0; i<Cloud->points.size(); i++) {
            
            // If the point is already an inlier, continue
            if (TargetPoints.count(i) > 0) {    continue;   }

            // Distance from a point to the plane
            pcl::PointXYZ Point = Cloud->points[i];
            float xi = Point.x;
            float yi = Point.y;
            float zi = Point.z;
            float Dist = fabs(a * xi + b * yi + c * zi)/den;

            // If the dist < threshold -> point is an inlier
            if (Dist < Threshold) {
                TargetPoints.insert(i);
            }
            
            // 4. Store the results of the plane that has more points
            if (TargetPoints.size() > InlierPoints.size()) {
                InlierPoints = TargetPoints;
            }
            
        }

    }

    // 5. Creating two new point clouds, one with obstacles and other with plane
    pcl::PointCloud<pcl::PointXYZ>::Ptr obstCloud (new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PointCloud<pcl::PointXYZ>::Ptr planeCloud (new pcl::PointCloud<pcl::PointXYZ>);

    // Divide the floor and the rest in 2 different point clouds
    for (int i = 0; i < (int) Cloud->points.size(); i++) {

        pcl::PointXYZ Point = Cloud->points[i];
        if (InlierPoints.count(i)) {
            planeCloud->points.push_back(Point);
        } else {
            obstCloud->points.push_back(Point);
        }

    }

    auto EndTime = std::chrono::steady_clock::now();
    auto ElapsedTime = std::chrono::duration_cast<std::chrono::milliseconds>(EndTime - StartTime);
    std::cout << "Plane segmentation (Ransac3D) took: " << ElapsedTime.count() << " ms" << std::endl;

    std::pair<pcl::PointCloud<pcl::PointXYZ>::Ptr, pcl::PointCloud<pcl::PointXYZ>::Ptr> segResult(obstCloud, planeCloud);
    return segResult;

}

 /****************************************************
 *  Clustering extraction
 * 
 *  I- Original cloud, Inliers to split
 *  O- New cloud with only inliers
 * **************************************************/

 void ClusteringExtraction (pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud, float Tolerance, int MinSize, int MaxSize, std::vector<Object>* outputObjects, int* numOutputObjects) 
 {
    // -- KD-tree object definition
    pcl::search::KdTree<pcl::PointXYZ>::Ptr Tree (new pcl::search::KdTree<pcl::PointXYZ>);
    Tree->setInputCloud(Cloud);

    // -- Configuration of the search method of extraction
    std::vector<pcl::PointIndices> ClusterIndices;
    pcl::EuclideanClusterExtraction<pcl::PointXYZ> EC;
    EC.setClusterTolerance(Tolerance);
    EC.setMinClusterSize(MinSize);
    EC.setMaxClusterSize(MaxSize);
    EC.setSearchMethod(Tree);
    EC.setInputCloud(Cloud);
    EC.extract(ClusterIndices);

    // -- Clusters storage
    for (std::vector<pcl::PointIndices>::const_iterator it = ClusterIndices.begin (); it != ClusterIndices.end (); ++it) {

        pcl::PointCloud<pcl::PointXYZ>::Ptr CloudCluster (new pcl::PointCloud<pcl::PointXYZ>);
        for (std::vector<int>::const_iterator pit = it->indices.begin (); pit != it->indices.end (); ++pit) {
            CloudCluster->points.push_back (Cloud->points[*pit]);
        }
        
        CloudCluster->width = CloudCluster->points.size ();
        CloudCluster->height = 1;
        CloudCluster->is_dense = true;

        // Initialize point cloud vertices
        float xmin = INFINITY, xmax = -INFINITY, ymin = INFINITY, ymax = -INFINITY, zmin = INFINITY, zmax = -INFINITY;
        float xcen = -INFINITY, ycen = -INFINITY, zcen = -INFINITY;

        // Cloud vertices search
        for (int i = 0; i < CloudCluster->points.size(); i++)
        {
            if (CloudCluster->points[i].x < xmin) {xmin = CloudCluster->points[i].x;}
            if (CloudCluster->points[i].x > xmax) {xmax = CloudCluster->points[i].x;}
            if (CloudCluster->points[i].y < ymin) {ymin = CloudCluster->points[i].y;}
            if (CloudCluster->points[i].y > ymax) {ymax = CloudCluster->points[i].y;}
            if (CloudCluster->points[i].z < zmin) {zmin = CloudCluster->points[i].z;}
            if (CloudCluster->points[i].z > zmax) {zmax = CloudCluster->points[i].z;}

        }

        // Calculus of the centroid
        xcen = (xmax + xmin) / 2.0;
        ycen = (ymax + ymin) / 2.0;
        zcen = (zmax + zmin) / 2.0;

        // Creation of an object type
        Object object;

        object.x_min = xmin;
        object.x_max = xmax;
        object.y_min = ymin;
        object.y_max = ymax;
        object.z_min = zmin;
        object.z_max = zmax;

        object.d = xmax - xmin;
        object.w = ymax - ymin;
        object.h = zmax - zmin;

        object.centroid_x = xcen;
        object.centroid_y = ycen;
        object.centroid_z = zcen;

        object.type = "none";
        object.cloud = CloudCluster;

        outputObjects->push_back(object);
        *numOutputObjects = *numOutputObjects + 1;
    }
}

// Extract clusters from the coloured XYZ filtered LiDAR point cloud according to the input cluster parameters
void cluster_filter(pcl::PointCloud<pcl::PointXYZ>::Ptr filtered_cloud, float tolerance, int min_cluster, int max_cluster, std::vector<Object> *output_objects, int *number_output_objects)
{
	// Parameters:
	// filtered_cloud: XYZ and angle filtered LiDAR point cloud that contains the clusters
	// tolerance: Tolerance of clusters
	// min_cluster: Minimum size of a cluster
	// max_cluster: Maximum size of a cluster
	// only_laser_objects: Pointer that points to the array that contains the clusters
	// only_laser_objects_number: Number of objects

	// This function only takes into account the size of the clusters

    // Extract clusters from point cloud

	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree (new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud (filtered_cloud);
	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
	ec.setClusterTolerance (tolerance); 
	ec.setMinClusterSize (min_cluster); 
	ec.setMaxClusterSize (max_cluster);
	ec.setSearchMethod (tree);
	ec.setInputCloud (filtered_cloud);
	ec.extract (cluster_indices);

    // Store the clusters

	for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
	{
		pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster (new pcl::PointCloud<pcl::PointXYZ>);
 
		for (std::vector<int>::const_iterator pit = it->indices.begin(); pit != it->indices.end(); ++pit)
		{
			cloud_cluster->points.push_back (filtered_cloud->points[*pit]); 
		}

		cloud_cluster->width = cloud_cluster->points.size();
		cloud_cluster->height = 1;
		cloud_cluster->is_dense = true;
 
		// Initialize point cloud vertices. Set to +/- INFINITY to ensure a proper behaviour for the first cluster

		float x_min = INFINITY; 
		float y_min = INFINITY;
		float z_min = INFINITY;
		float x_max = -INFINITY;
		float y_max = -INFINITY;
		float z_max = -INFINITY;
 
		float centroid_x = -INFINITY;
		float centroid_y = -INFINITY;
		float centroid_z = -INFINITY;
 
		for (int i = 0; i < cloud_cluster->points.size(); i++)
		{
			if (cloud_cluster->points[i].x < x_min)		
			{
				x_min = cloud_cluster->points[i].x;
			}
 
			if (cloud_cluster->points[i].y < y_min)		
			{
				y_min = cloud_cluster->points[i].y;
			}
 
			if (cloud_cluster->points[i].z < z_min)	
			{
				z_min = cloud_cluster->points[i].z;		
			}
			if (cloud_cluster->points[i].x > x_max)
			{
				x_max = cloud_cluster->points[i].x;
			}
			if (cloud_cluster->points[i].y > y_max)
			{
				y_max = cloud_cluster->points[i].y;
			}
			if (cloud_cluster->points[i].z > z_max)
			{
				z_max = cloud_cluster->points[i].z;
			}
		}

		// Centroid

		centroid_x = (x_max+x_min)/2.0;
		centroid_y = (y_max+y_min)/2.0;
		centroid_z = (z_max+z_min)/2.0;

		Eigen::Vector4f centroid;
		pcl::compute3DCentroid(*cloud_cluster,centroid);

		geometry_msgs::PointStamped local_centroid;
		geometry_msgs::Point32 global_centroid;

		local_centroid.point.x = centroid_x;
		local_centroid.point.y = centroid_y;
		local_centroid.point.z = centroid_z;

        Object object;

        object.x_max = x_max;
        object.x_min = x_min;
        object.y_max = y_max;
        object.y_min = y_min;
        object.z_max = z_max;
        object.z_min = z_min;

        object.d = x_max - x_min;
        object.w = y_max - y_min;
        object.h = z_max - z_min;

        object.centroid_x = centroid[0];
        object.centroid_y = centroid[1];
        object.centroid_z = centroid[2];

        // Type

        object.type = "none";

        // Cloud

        object.cloud = cloud_cluster;

        output_objects->push_back(object);
        *number_output_objects = *number_output_objects + 1;        
	}
}



// -- Da error: al recorrer con un puntero no analizas el resto de puntos. Carlos lo hace con un vector
void ClusterHelper (int idx, std::vector<std::vector<float>> Cloud, std::vector<int>& Cluster,
                     std::vector<bool>& Processed, KdTree* Tree, float DistanceTol) {

    Processed[idx] = true;
    Cluster.push_back(idx);
    std::vector<int> nearest = Tree->search(Cloud[idx], DistanceTol);

    for (auto id : nearest) {
        if(!Processed[id]) {     ClusterHelper(id, Cloud, Cluster, Processed, Tree, DistanceTol);     }
    }

}

std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> EuclideanClustering (pcl::PointCloud<pcl::PointXYZ>::Ptr Cloud, float Tolerance, int MinSize, int MaxSize) {

    // 1. Create and fill KdTree based on cloud points
    KdTree* Tree = new KdTree;
    std::vector<std::vector<float>> treePoints;

    for (int idx = 0; idx < Cloud->points.size(); idx++) {

        std::vector<float> Point = {Cloud->points[idx].x, Cloud->points[idx].y, Cloud->points[idx].z};
        treePoints.push_back(Point);
        Tree->insert(Point, idx);

    }

    // 2. Perform Euclidean Clustering based on hyperparameters
    //      Obtaining a vector of indices for each cluster
    std::vector<std::vector<int>> clustersIndicesVector;
    std::vector<bool> processedPoints(treePoints.size(), false);

    int idy = 0;
    while (idy < treePoints.size()) {

        // If the point was previously processed, continue
        if (processedPoints[idy]) {
            idy++;
            continue;
        }

        // Else, create a new cluster
        std::vector<int> Cluster;
        ClusterHelper(idy, treePoints, Cluster,  processedPoints, Tree, Tolerance);
        clustersIndicesVector.push_back(Cluster);
        idy++;

    }

    // 3. Obtain individual clusters based on clusters indices vector
    std::vector<pcl::PointCloud<pcl::PointXYZ>::Ptr> clustersVector;
    for (std::vector<int> clusterIndex : clustersIndicesVector) {

        // Conversion from vector of indices to pcl::PointCloud
        pcl::PointCloud<pcl::PointXYZ>::Ptr cloudCluster (new pcl::PointCloud<pcl::PointXYZ>);
        for (int idz : clusterIndex) {  cloudCluster->points.push_back(Cloud->points[idz]);  }

        // If the cluster meets the size requirements, push it to pcl::PointCloud output vector. Else, discard.
        if ((cloudCluster->points.size() > MinSize) && (cloudCluster->points.size() < MaxSize)) {

            // Adding extra information to the cluster
            cloudCluster->width = cloudCluster->points.size();
            cloudCluster->height = 1;
            cloudCluster->is_dense = true;

            // Pushing to output vector
            clustersVector.push_back(cloudCluster);

        }
        
    }

    return clustersVector;

}
