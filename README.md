# Graph-cluster
## Louvain
## Evaluation
## Temp
* 原始数据：supply 文件
* 转化为一般的图文件，保存节点的映射关系
* 聚类：用到 Louvain 代码，输入一般的图文件，输出分类结果
* 分类结果可以通过节点的映射关系重新建模成多个 supply 文件
* 评估比较
    * 建模成 lp 文件（聚类之前和之后）
    * Gurobi 求解
    * 比较两者的结果