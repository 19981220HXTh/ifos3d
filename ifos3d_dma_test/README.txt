编译：
先到IFOS3D_DMA_TEST目录下面
执行指令 scons

运行：
cd par/toy-1/
可以在json_crt.py修改进程分配数目
./run.py fw
./run.py inv

画图：
cd model
可以在plot_vp.sh修改画图的切面
./plot_vp.sh

清除编译：
scons -c
