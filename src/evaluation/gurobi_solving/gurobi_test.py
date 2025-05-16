from gurobipy import *

# 定义输入文件路径和输出文件路径
lp_file = "/pub/netdisk1/lijy/alimama/lp_learning/test_0416_cleaned.lp"  # 输入的 LP 文件路径
log_file = "/pub/netdisk1/lijy/alimama/lp_learning/test_0416_cleaned_gurobi.log"     # 日志输出文件路径
solution_file = "/pub/netdisk1/lijy/alimama/lp_learning/test_0416_cleaned_gurobi.sol"  # 变量值输出文件路径

# 创建模型并读取 LP 文件
model = read(lp_file)

# 设置 Gurobi 参数
model.setParam("LogFile", log_file)   # 设置日志文件
# model.setParam("TimeLimit", 300)      # 设置求解时间限制为 300 秒
model.setParam("OutputFlag", 1)       # 开启输出日志

# 求解模型
model.optimize()

# 如果模型成功求解，则输出每个变量的值到文件
if model.status == GRB.OPTIMAL:
    with open(solution_file, "w") as f:
        for var in model.getVars():
            f.write(f"{var.VarName}: {var.X}\n")
    print(f"最优解已保存到 {solution_file}")
else:
    print("模型未找到最优解或求解未完成。")

print(f"日志已保存到 {log_file}")
