import os
from collections import defaultdict

def load_id_to_index(file_path):
    id_to_index = {}
    with open(file_path, 'r') as f:
        for line in f:
            if '->' in line:
                id_str, index_str = line.strip().split('->')
                id_str = id_str.strip()
                index = int(index_str.strip())
                id_to_index[index] = id_str
    return id_to_index

def load_index_to_cluster(file_path):
    index_to_cluster = {}
    with open(file_path, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) == 2:
                index = int(parts[0])
                cluster = int(parts[1])
                index_to_cluster[index] = cluster
    return index_to_cluster

def split_ids_by_cluster(id_to_index_path, index_to_cluster_path, output_dir):
    id_map = load_id_to_index(id_to_index_path)
    cluster_map = load_index_to_cluster(index_to_cluster_path)

    cluster_to_ids = defaultdict(list)

    for index, cluster in cluster_map.items():
        if index in id_map:
            cluster_to_ids[cluster].append(id_map[index])

    os.makedirs(output_dir, exist_ok=True)
    for cluster, ids in cluster_to_ids.items():
        output_file = os.path.join(output_dir, f'cluster_{cluster}.txt')
        with open(output_file, 'w') as f:
            for id_str in ids:
                f.write(f"{id_str}\n")

# 示例调用
split_ids_by_cluster('node_mapping.txt', '/home/lijy/repo/BipartiteGraphClustering/Community_latest/communities.txt', 'clusters')
