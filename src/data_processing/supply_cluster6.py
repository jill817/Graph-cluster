def load_ids(id_file):
    supply_ids = set()
    demand_ids = set()
    with open(id_file, 'r') as f:
        for line in f:
            id_str = line.strip()
            if '_' in id_str:
                supply_ids.add(id_str)
            else:
                demand_ids.add(id_str)
    return supply_ids, demand_ids

def filter_supply_demand(supply_demand_file, supply_ids, demand_ids, output_file):
    with open(supply_demand_file, 'r') as fin, open(output_file, 'w') as fout:
        for line in fin:
            parts = line.strip().split()
            if len(parts) < 4:
                continue  # 跳过无效行
            supply_node_id, supply_location_id, weight = parts[:3]
            full_supply_id = f"{supply_node_id}_{supply_location_id}"
            if full_supply_id in supply_ids:
                filtered_demands = [d for d in parts[3:] if d in demand_ids]
                if filtered_demands:  # 至少有一个 demand id 被保留
                    new_line = ' '.join([supply_node_id, supply_location_id, weight] + filtered_demands)
                    fout.write(new_line + '\n')

# 示例调用
supply_ids, demand_ids = load_ids('/pub/netdisk1/lijy/alimama/clusters/cluster_6.txt')
filter_supply_demand('supply_0416_connect.txt', supply_ids, demand_ids, 'supply_0416_cluster6.txt')
