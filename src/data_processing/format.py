def convert_line(line):
    items = line.strip().split()
    if len(items) < 4:
        return '`'.join(items)
    else:
        first_part = '`'.join(items[:3])
        second_part = ','.join(items[3:])
        return f"{first_part}`{second_part}"

def process_file(input_filename, output_filename):
    with open(input_filename, 'r') as infile, open(output_filename, 'w') as outfile:
        for line in infile:
            converted_line = convert_line(line)
            outfile.write(converted_line + '\n')

# 使用方法：
input_filename = '/pub/netdisk1/lijy/alimama/supply_0416_cluster6.txt'  # 输入文件名
output_filename = '/pub/netdisk1/lijy/alimama/supply_0416_cluster6_format.txt'  # 输出文件名
process_file(input_filename, output_filename)
