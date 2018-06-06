"""
Created on Sat Dec 16 11:37:05 2017

@author: Yusuf Tatlier
"""

import math 
import alg_cluster

#%%
# Help functions

def scl_collisions(cluster_list,factor=1):
    """
    Generalized version of the sort_cluster_list algorithm that can deal with collisions
    """
    cluster_list_copy=[x*factor for x in list(cluster_list)]
    order_list=[]
    
    while(len(cluster_list_copy)>0):
        min_value=min(cluster_list_copy)
        min_value_indices=[i for i, x in enumerate(cluster_list) if x == min_value*factor]
        for element_i in min_value_indices:
            order_list.append(element_i)
            cluster_list_copy.remove(min_value)
            
    return order_list
    

def min_cluster_tuple(closest_pair_tuple_1,closest_pair_tuple_2):
    """
    Chooses the tuple with the lowest distance between the pairs (elements 1 and 2)
    """
    if(closest_pair_tuple_1[0]<closest_pair_tuple_2[0]):
        return closest_pair_tuple_1
    else:
        return closest_pair_tuple_2

def get_coordinates(cluster_list,dimension):
    """
    Gets horizontal or vertical coordinates for a cluster list
    """
    if(dimension=="horizontal"):
        return [x.horiz_center() for x in cluster_list]
    else:
        return [x.vert_center() for x in cluster_list]
    
def sort_cluster_list(cluster_list):
    """
    Returns the indices of cluster list values sorted in increasing order
    """
    #Check for collisions 
    return [cluster_list.index(x) for x in sorted(cluster_list)]

def sort_cluster(cluster_list,dimension,factor=1):    
    """
    Returns a sorted cluster list, based on given sort dimension (horizontal, vertical coordinates)
    """
    index_list=scl_collisions(get_coordinates(cluster_list,dimension),factor)
    #index_list=sort_cluster_list(get_coordinates(cluster_list,dimension))
    return [cluster_list[x] for x in index_list]
    
#%%
#Algorithms built for Cluster class

def slow_closest_pair(cluster_list):
    """
    Implementation of the (brute-force) slow closest pairs algorithm.
    Takes a list of Cluster objects and returns a closest pair where the pair is represented by the tuple (dist, idx1, idx2) 
    with idx1 < idx2 where dist is the distance between the closest pair cluster_list[idx1] and cluster_list[idx2].
    """
    closest_pair_tuple=(float('inf'),-1,-1)
    cluster_list_length=len(cluster_list)
    for index_i in range(cluster_list_length):
        for index_j in [x for x in range(cluster_list_length) if x != index_i]:
            dist_pairs=cluster_list[index_i].distance(cluster_list[index_j])
            if(dist_pairs<closest_pair_tuple[0]):
                closest_pair_tuple=(dist_pairs,min(index_i,index_j),max(index_i,index_j))
    return closest_pair_tuple

#%%

def fast_closest_pair(cluster_list):
    """
    Implementation of the fast closest pairs algorithm.
    Takes a list of Cluster objects and returns a closest pair where the pair is represented by the tuple (dist, idx1, idx2) 
    with idx1 < idx2 where dist is the distance between the closest pair cluster_list[idx1] and cluster_list[idx2].
    """
    cluster_list_indices_hor=scl_collisions(get_coordinates(cluster_list,"horizontal"))
    cluster_list=sort_cluster(cluster_list,"horizontal")
    cluster_length=len(cluster_list)
    
    #Base case
    closest_pair_tuple=(float('inf'),-1,-1)
    if(cluster_length<4):
        closest_pair_tuple=slow_closest_pair(cluster_list)
        return closest_pair_tuple
    else:
    #Divide the problem (recursion)
        half_cluster_length=cluster_length//2
        cluster_list_left=cluster_list[:half_cluster_length]
        cluster_list_right=cluster_list[half_cluster_length:]
        pair_tuple_left=fast_closest_pair(cluster_list_left)
        pair_tuple_right=fast_closest_pair(cluster_list_right)
    #Merge 
        closest_pair_tuple=min_cluster_tuple(pair_tuple_left,(pair_tuple_right[0],pair_tuple_right[1]+half_cluster_length,pair_tuple_right[2]+half_cluster_length))
        mid=0.5*(cluster_list[half_cluster_length-1].horiz_center()+cluster_list[half_cluster_length].horiz_center())
        closest_pair_tuple=min_cluster_tuple(closest_pair_tuple,closest_pair_strip(cluster_list,mid,closest_pair_tuple[0]))
        index_i, index_j = closest_pair_tuple[1], closest_pair_tuple[2]
        tr_index_i_h, tr_index_j_h = cluster_list_indices_hor[index_i], cluster_list_indices_hor[index_j]
        return closest_pair_tuple[0],tr_index_i_h,tr_index_j_h


def prepare_strip(cluster_list,horiz_center,half_width):
    """
    Help function for closest_pair_strip, as function was getting too big. In this way the total was split into two smaller parts
    """
    #List of points in the strip
    strip_cluster_list=[x for x in cluster_list if abs(horiz_center-x.horiz_center())<half_width]
    strip_length=len(strip_cluster_list)
    
    #Determine startindex, for first point in the strip
    left_number_points_total=sum([y.horiz_center()<horiz_center for y in cluster_list])
    left_number_points_strip=sum([y.horiz_center()<horiz_center for y in strip_cluster_list])
    start_index=left_number_points_total-left_number_points_strip
    
    #sort the strip vertically
    cluster_list_indices_ver=scl_collisions(get_coordinates(strip_cluster_list,"vertical"))
    vert_sorted_strip_cluster_list=sort_cluster(strip_cluster_list,"vertical")
    return start_index,strip_length,cluster_list_indices_ver,vert_sorted_strip_cluster_list
    
def closest_pair_strip(cluster_list,horiz_center,half_width):
    """
    Takes a list of Cluster objects and two floats horiz_center and half_width. horiz_center specifies the horizontal position of the center 
    line for a vertical strip. half_width specifies the maximal distance of any point in the strip from the center line.
    """

    start_index,strip_length,cluster_list_indices_ver,vert_sorted_strip_cluster_list=prepare_strip(cluster_list,horiz_center,half_width)
    
    closest_pair_tuple=(float('inf'),-1,-1)
    for index_i in range(strip_length-1):
        for index_j in [x for x in range(index_i+1,min(index_i+4,strip_length)) if x !=index_i]:
            distance=vert_sorted_strip_cluster_list[index_i].distance(vert_sorted_strip_cluster_list[index_j])
            tr_index_i_v, tr_index_j_v = cluster_list_indices_ver[index_i]+start_index, cluster_list_indices_ver[index_j]+start_index
            closest_pair_tuple=min_cluster_tuple(closest_pair_tuple,(distance,min(tr_index_i_v,tr_index_j_v),max(tr_index_i_v,tr_index_j_v)))
           
    return closest_pair_tuple
#%%
#Hierarchical clustering

def hierarchical_clustering(cluster_list, num_clusters):
    """
    Implementation of hierarchical clustering algorithm. Clusters based on the first closest pair.
    """
    
    cluster_len=len(cluster_list)
    cluster_list_copy=list(cluster_list)
    
    while(cluster_len>num_clusters):
        closest_pair_tuple=fast_closest_pair(cluster_list_copy)
        cluster_list_copy[closest_pair_tuple[1]].merge_clusters(cluster_list_copy[closest_pair_tuple[2]])
        del cluster_list_copy[closest_pair_tuple[2]]
        cluster_len=len(cluster_list_copy)
    return cluster_list_copy


#%%
#K-means-clustering

def largest_population_cluster(cluster_list,num_clusters):
    """
    Returns the indices of the num_clusters most populated clusters
    """
    population = [x.total_population() for x in cluster_list]
    population_sort = scl_collisions(population,-1)
    return population_sort[0:num_clusters]

def closest_cluster(point,centers_list):
    """
    Returns the cluster with the center that is closest to point
    We go for a linear search, as it is O(n)
    """
    closest_cluster_index=(float('inf'),-1)
    cluster_list_length=len(centers_list)
    for index_i in range(cluster_list_length):
        dist_to_cluster=point.distance(centers_list[index_i])
        if(dist_to_cluster<closest_cluster_index[0]):
            closest_cluster_index=(dist_to_cluster,index_i)
    return closest_cluster_index

def kmeans_clustering(cluster_list, num_clusters, num_iterations):
    """
    Implementation of the k-means clustering algorithm
    """
    cluster_list_copy=list(cluster_list)
    cluster_list_len=len(cluster_list)
    key_fun=lambda x: x.total_population()
    main_clusters_list=sorted(cluster_list_copy,key=key_fun, reverse=True)[:num_clusters]
    #population_sort=largest_population_cluster(cluster_list_copy,num_clusters)
    #main_clusters_list=[cluster_list_copy[index_c_i] for index_c_i in population_sort]
    #centers_list=[(x.horiz_center(),x.vert_center()) for x in main_clusters_list]
    #centers_list=[(cluster_list[index_c_i]._horiz_center,cluster_list[index_c_i]._vert_center) for index_c_i in population_sort]
    
    for _ in range(1,num_iterations+1):
        c_sets=[alg_cluster.Cluster(set([]),0,0,0,0) for _ in range(num_clusters)]
        for index_j in range(cluster_list_len):
            closest_cluster_index=closest_cluster(cluster_list_copy[index_j],main_clusters_list)
            c_sets[closest_cluster_index[1]].merge_clusters(cluster_list_copy[index_j])
            
        main_clusters_list=list(c_sets)
    
    return main_clusters_list

#%%
# Application part
#Q1
import random
import time
import matplotlib.pyplot as plt

#Create random clusters in the square x,y in [-1,1]
def gen_random_clusters(num_clusters):
    xy_coordinates=[(random.uniform(0,1),random.uniform(-1,1)) for _ in range(num_clusters)]
    random_cluster=[alg_cluster.Cluster(set([]),xy_coordinates[index_i][0],xy_coordinates[index_i][1],1,0) for index_i in range(num_clusters)]
    return random_cluster

def running_times(ref_fun,value):

    start_time = time.clock()
    ref_fun(value)
    return(time.clock() - start_time)

calculational_time_slow, calculational_time_fast=[],[]
random_cluster_list=gen_random_clusters(200)

for index_i in range(200):
    use_list=random_cluster_list[:index_i]
    calculational_time_slow.append(running_times(slow_closest_pair,use_list))
    if(index_i>3):
        calculational_time_fast.append(running_times(fast_closest_pair,use_list))
    
plt.xlabel('Number of clusters')
plt.ylabel('Running times')
plt.title('Running times fast vs slow cluster pair algorithms')
    
plt.subplot(111)
plt.plot(calculational_time_slow,'b',
           calculational_time_fast,'r')
plt.legend(('slow_closest_pair','fast_closest_pair'),loc='upper right')
plt.show()

#%%

def compute_distortion(cluster_list,data_table):
        """
        Takes a list of clusters and uses cluster_error to compute its distortion.
        
        input: list of clusters, original data table
        output: cluster distortion int
        """
        return sum([x.cluster_error(data_table) for x in cluster_list])

from Cancer_Risk_Data import DATA_111_URL, DATA_290_URL, DATA_896_URL, load_data_table, sequential_clustering, run_example

def load_data(data_to_load,number):
    data_table = load_data_table(data_to_load)
    
    singleton_list = []
    for line in data_table:
        singleton_list.append(alg_cluster.Cluster(set([line[0]]), line[1], line[2], line[3], line[4]))
        
    cluster_list = sequential_clustering(singleton_list, number)	
    
    return data_table, cluster_list
#%%
#Import cluster_list and data_table

def make_distortion_plots(data_to_load,title,number):

    distortion_list_kmeans=[]
    distortion_list_hierarchical=[]

    for index_i in range(6,21):
        data_table, cluster_list = load_data(data_to_load,number)
        new_clusters_km=kmeans_clustering(cluster_list,index_i,5)
        distortion_list_kmeans.append(compute_distortion(new_clusters_km,data_table))
    
    for index_i in range(6,21):
        data_table, cluster_list = load_data(data_to_load,number)
        new_clusters_h=hierarchical_clustering(cluster_list,index_i)
        distortion_list_hierarchical.append(compute_distortion(new_clusters_h,data_table))
    
    x_axis=[x for x in range(6,21)]
    plt.xlabel('Number of clusters')
    plt.ylabel('Distortion')
    plt.title('Distortion vs number of clusters for ' + str(title))
    
    plt.subplot(111)
    plt.plot(x_axis,distortion_list_kmeans,'b',x_axis,distortion_list_hierarchical,'r')
    plt.legend(('kmeans','hierarchical'),loc='upper right')
    plt.show()

data_to_load=DATA_896_URL
make_distortion_plots(DATA_896_URL,"DATA_896_URL",896)
