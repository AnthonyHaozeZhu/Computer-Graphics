# Conditional GAN
本项目利用使用清华大学开源代码框架Jittor利用主流的手写数字数据集MNIST实现了Conditional GAN网络，目标是用计算机生成手写数字。

### 项目框架如下

```.
├── CGAN.py
├── example
│   ├── 0.png
│   ├── 10.png
│   ├── 20.png
│   ├── 30.png
│   ├── 40.png
│   ├── 50.png
│   ├── 60.png
│   ├── 70.png
│   ├── 80.png
│   └── 90.png
├── log
├── result.png
└── saved_models
	  ├── discriminator_last.pkl
	  └── generator_last.pkl
```

其中 CGAN.py 包含了模型和训练函数

example 中存储了训练时生成的样例图片

log 记录了训练时的数据

saved_models 中保存了

result.png 为您想生成的数字串的图片



### 项目部署

**环境依赖**

```
python                 3.8.10
numpy                  1.22.3
jittor                 1.3.4.3
tqdm                   4.64.0
```

推荐运行环境 Ubuntu Server 20.04

关于jitter安装参考 https://nbviewer.jupyter.org/github/Jittor/LearnJittorBasicIn60Min/tree/master/



### 项目运行

```bash
python CGAN.py --epochs 训练轮次 --batch_size your batch size --lr 学习率 --latent_dim 隐向量维度 --n_class 分类个数 --img_size 图片的大小(正方形图片边长) --channels 图片的通道 --device 训练设备(CPU/GPU)
```



### 生成示例

![result](result.png)
