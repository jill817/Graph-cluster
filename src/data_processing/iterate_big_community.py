# 文件名设置
vertex_file = "/home/lijy/repo/BipartiteGraphClustering/Community_latest/community3.txt"
edge_file = "graph.txt"
output_file = "filtered_graph.txt"

# 1. 读取顶点集合
with open(vertex_file, 'r') as f:
    vertex_ids = set()
    for line in f:
        if line.strip():
            vertex_id = int(line.strip().split()[0])
            vertex_ids.add(vertex_id)

# 2. 读取边文件，筛选相关边
with open(edge_file, 'r') as f_in, open(output_file, 'w') as f_out:
    for line in f_in:
        if line.strip():
            parts = line.strip().split()
            u = int(parts[0])
            v = int(parts[1])
            if u in vertex_ids or v in vertex_ids:
                f_out.write(line)
