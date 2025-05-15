def generate_hypergraph(supply_file, demand_file, output_file):
    """
    生成超图文件
    
    参数：
    supply_file: supply文件路径，每行格式为：id1 id2 weight node1 node2 node3 ...
    demand_file: demand文件路径，每行格式为：node_id`weight
    output_file: 输出文件路径，每行格式为：demand_node_id weight supply_node_id1 supply_node_id2 ...
    """
    # 用于存储demand节点到supply节点的映射
    demand_to_supply = {}
    
    # 读取supply文件，建立映射关系
    with open(supply_file, 'r') as f:
        for line in f:
            parts = line.strip().split()
            if len(parts) < 4:
                continue
                
            id1 = parts[0]
            id2 = parts[1]
            supply_id = f"{id1}_{id2}"
            demand_nodes = parts[3:]
            
            # 为每个demand节点添加对应的supply节点
            for demand_node in demand_nodes:
                if demand_node not in demand_to_supply:
                    demand_to_supply[demand_node] = set()
                demand_to_supply[demand_node].add(supply_id)
    
    # 读取demand文件，获取权重
    demand_weights = {}
    with open(demand_file, 'r') as f:
        for line in f:
            parts = line.strip().split('`')
            if len(parts) != 2:
                continue
            demand_node = parts[0]
            weight = float(parts[1])
            demand_weights[demand_node] = weight
    
    # 生成超图文件
    with open(output_file, 'w') as f:
        for demand_node, supply_nodes in demand_to_supply.items():
            if demand_node in demand_weights:
                # 将supply节点集合转换为排序后的列表
                supply_list = sorted(list(supply_nodes))
                # 写入格式：demand_node_id weight supply_node_id1 supply_node_id2 ...
                f.write(f"{demand_node} {demand_weights[demand_node]} {' '.join(supply_list)}\n")

if __name__ == "__main__":
    supply_file = "supply_0416_connect.txt"
    demand_file = "demand_330_0416.txt"
    output_file = "hypergraph.txt"
    
    generate_hypergraph(supply_file, demand_file, output_file)
    print("超图文件生成完成！") 