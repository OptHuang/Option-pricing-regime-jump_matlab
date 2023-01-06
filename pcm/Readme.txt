代码说明

主函数 (直接运行部分)
1. Pic3DPCM.m: 绘制三维图像
2. Err.m: 计算给定分划下末端数值解的离散L2误差
3. ErrOrder.m: 绘制空间方向的收敛阶
4. ComputeVtrue.m: 设置细密分划计算数值解, 视其为真解

子函数
1. A_K.m/A_M.m: 分别计算Kou和Merton下对应的积分余量函数A(j,n)
2. AssembleMatrixKou.m/AssembleMatrixMerton.m: 分别形成Kou和Merton下对应的LCP问题对应的系数矩阵和向量 (PCM求解时需要)
3. AssembleRKou.m/AssembleRMerton.m: 分别形成Kou和Merton下对应的离散矩阵R (Toeplitz矩阵)
4. f.m/h.m: 形成积分矩阵R时需要, 与文章中符号对应
5. DataMatrixKou.m/DataMatrixMerton.m: 求解数值解并存入矩阵
6. FreeBoundry.m: 求解最佳实施边界
7. Koupara.m/Merton_func.m/Mertonpara.m/ParaImput.m: 都是用于参数设定以及一些文章中的中间变量的计算
8. MeshGeneration.m: 时空网格划分 (绘三维图)
9. Newton_perMerton.m: 牛顿法求解Merton永久美式的自由边界
10. PCMPro.m: PCM算法 (改进版的)
11. RevValueKou.m/RevValueMerton.m: 对 PCMPro.m 求解出的数据进行还原(PCM求解的是做了变换之后的模型)
12. TruncationTech.m: 截断技巧函数, 输出阶段后矩形区间空间方向对应的L
13. V_star.m: 边界函数
14. Vkoutrue20002560.mat/Vmertrue20002560.mat: 分别为Kou和Merton下的(时间分划2000,空间2560)数值解,即“真解”