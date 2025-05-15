def convert_hypergraph_to_edges(input_file, output_file):
    """
    将超图格式转换为边列表格式，并对节点进行重新编号
    
    参数：
    input_file: 输入文件路径，每行格式为：id1 id2 weight node1 node2 node3 ...
    output_file: 输出文件路径，每行格式为：src dest weight
    """
    # 用于存储节点映射关系
    node_mapping = {}
    next_id = 0
    
    # 第一遍读取：收集所有节点并建立映射
    with open(input_file, 'r') as f_in:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
                
            id1 = parts[0]
            id2 = parts[1]
            nodes = parts[3:]
            
            # 处理源节点
            src = f"{id1}_{id2}"
            if src not in node_mapping:
                node_mapping[src] = next_id
                next_id += 1
            
            # 处理目标节点
            for node in nodes:
                if node not in node_mapping:
                    node_mapping[node] = next_id
                    next_id += 1
    
    # 第二遍读取：使用映射后的节点编号写入文件
    with open(input_file, 'r') as f_in, open(output_file, 'w') as f_out:
        for line in f_in:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
                
            id1 = parts[0]
            id2 = parts[1]
            weight = float(parts[2])
            nodes = parts[3:]
            
            src = f"{id1}_{id2}"
            mapped_src = node_mapping[src]
            
            for dest in nodes:
                mapped_dest = node_mapping[dest]
                f_out.write(f"{mapped_src} {mapped_dest} {weight}\n")
    
    # 可选：保存映射关系到文件
    with open('node_mapping.txt', 'w') as f_map:
        for original, mapped in node_mapping.items():
            f_map.write(f"{original} -> {mapped}\n")
    
    print(f"节点总数: {len(node_mapping)}")
    print(f"映射关系已保存到 node_mapping.txt")

# 使用示例
input_file = "supply_0416_connect.txt"  # 替换为您的输入文件路径
output_file = "graph.txt"  # 输出文件路径

convert_hypergraph_to_edges(input_file, output_file)