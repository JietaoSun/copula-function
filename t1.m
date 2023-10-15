clear;
clc;
% 1、提取各个监测点的残差数据
Monitor=readmatrix('监测点残差数据-正式2.xlsx');%注意将残差数据读取为数值矩阵形式
X1=Monitor(:,1);%提取第1列残差数据
X2=Monitor(:,2);
X3=Monitor(:,3);
X4=Monitor(:,4);%提取第4列残差数据


%2、初步判断是否服从正态分布，包括偏度（skewness函数，为负则左偏，即密度函数左尾巴长，顶点偏右；为0，即对称；为正则右偏）（反映是否对称）、峰度（kurtosis函数，正态分布峰度为3，若大于3说明密度函数曲线在峰值附近陡于正态分布；否则平缓于正态分布）（反映峰值陡峭程度，峰（尖、平）尾（轻、薄）））

x1s=skewness(X1); %计算X1偏度
x2s=skewness(X2); %计算X2偏度
x3s=skewness(X3); %计算X3偏度
x4s=skewness(X4); %计算X4偏度


x1k=kurtosis(X1); %计算X1峰度
x2k=skewness(X2); %计算X2峰度
x3k=skewness(X3); %计算X3峰度
x4k=skewness(X4); %计算X4峰度

% 3、数据正态性检验进一步判断是否服从正态分布（包括jbtest（Jarque-Bera检验）、ktest（Kolmogorov-Smirnov检验）、lilletest（Lilliefors检验)，
% 若检验均通过，则h=0，否则h=1(默认p=0.01)则使用参数估计中的正态估计(normfit函数或Matlab统计工具箱等进行未知参数估计)，若检验不通过则使用非参数估计中的核密度估计、核分布估计）

%针对Xi进行jbtest（Jarque-Bera检验），hi=0正态性检验通过，hi=1检验不通过，并给出对应的显著性水平pi
[h1_jbtest,p1_jbtest]=jbtest(X1); 
[h2_jbtest,p2_jbtest]=jbtest(X2); 
[h3_jbtest,p3_jbtest]=jbtest(X3); 
[h4_jbtest,p4_jbtest]=jbtest(X4); 


% %针对Xi进行ktest（Kolmogorov-Smirnov检验），hi=0正态性检验通过，hi=1检验不通过，并给出对应的显著性水平pi
% [h1_ktest,p1_ktest]=kstest(X1,[X1,normcdf(X1,mean(X1),std(X1))]);
% [h2_ktest,p2_ktest]=kstest(X2,[X2,normcdf(X2,mean(X2),std(X2))]);
% [h3_ktest,p3_ktest]=kstest(X3,[X3,normcdf(X3,mean(X3),std(X3))]);
% [h4_ktest,p4_ktest]=kstest(X4,[X4,normcdf(X4,mean(X4),std(X4))]);


%针对Xi进行lilletest（Lilliefors检验)，hi=0正态性检验通过，hi=1检验不通过，并给出对应的显著性水平pi
[h1_lilletest,p1_lilletest]=lillietest(X1);
[h2_lilletest,p2_lilletest]=lillietest(X2);
[h3_lilletest,p3_lilletest]=lillietest(X3);
[h4_lilletest,p4_lilletest]=lillietest(X4);

 % 4、核密度估计函数与核分布估计函数
% 核密度估计函数计算
[f_kx1,~,u1]=ksdensity(X1);%x1c是在第一列残差数据X1的取值范围内等间隔选取的的100个点构成的向量,f_kx1是与x1c相应的核密度估计值向量,ui是Gaussian核函数由样本数量、标准差计算的的最佳窗宽
[f_kx2,~,u2]=ksdensity(X2);
[f_kx3,~,u3]=ksdensity(X3);
[f_kx4,~,u4]=ksdensity(X4);
u1;   %查看最佳窗宽值
u2;
u3;
u4;
% 核分布估计函数
[F_kx1,x1c]=ksdensity(X1,'function','cdf');%x1c是在第一列残差数据X1的取值范围内等间隔选取的的100个点构成的向量,F_kx1是与x1c相应的核分布估计值向量
[F_kx2,x2c]=ksdensity(X2,'function','cdf');%核分布估计值向量F_kx2
[F_kx3,x3c]=ksdensity(X3,'function','cdf');%核分布估计值向量F_kx3
[F_kx4,x4c]=ksdensity(X4,'function','cdf');%核分布估计值向量F_kx4、


% 5、核密度估计函数与频率直方图的拟合
% 第一列残差数据的频率分布直方图与核密度估计函数的拟合
figure;
histogram(X1,30,'Normalization','pdf', 'FaceColor',[142/255 160/255 199/255]);  % 绘制直方图
hold on;
plot(x1c,f_kx1,'--','Color',[1/255 1/255 1/255],'LineWidth',6); % 绘制核密度估计函数图
xlabel('Residual(mm)', 'FontName', 'Times New Roman','FontSize', 56);
ylabel('Frequency/group spacing','FontName', 'Times New Roman','FontSize', 56);
legend('Frequency distribution histogram','Kernel density function {\itf}_1({\itx}_1)','Location','northeast','FontSize',38, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50,'LineWidth', 2,'FontName', 'Times New Roman','LooseInset', [0,0,0,0]);


% 第二列残差数据的频率分布直方图与核密度估计函数的拟合
figure;
histogram(X2,30,'Normalization','pdf', 'FaceColor',[142/255 160/255 199/255]);  % 绘制直方图
hold on;
plot(x2c,f_kx2,'--','Color',[0 0 0],'LineWidth',6); % 绘制核密度估计函数图
xlabel('Residual(mm)', 'FontName', 'Times New Roman','FontSize', 56);
ylabel('Frequency/group spacing','FontName', 'Times New Roman','FontSize', 56);
legend('Frequency distribution histogram','Kernel density function {\itf}_2({\itx}_2)','Location','northeast','FontSize',40, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50,'LineWidth', 2,'Fontname', 'Times new Roman','LooseInset', [0,0,0,0]);

% 第三列残差数据的频率分布直方图与核密度估计函数的拟合
figure;
histogram(X3,30,'Normalization','pdf', 'FaceColor',[142/255 160/255 199/255]);  % 绘制直方图
hold on;
plot(x3c,f_kx3,'--','Color',[0 0 0],'LineWidth',6); % 绘制核密度估计函数图
xlabel('Residual(mm)', 'FontName', 'Times New Roman','FontSize', 56);
ylabel('Frequency/group spacing','FontName', 'Times New Roman','FontSize', 56);
legend('Frequency distribution histogram','Kernel density function {\itf}_3({\itx}_3)','Location','northeast','FontSize',40, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50,'LineWidth', 2,'Fontname', 'Times new Roman','LooseInset', [0,0,0,0]);

% 第四列残差数据的频率分布直方图与核密度估计函数的拟合
figure;
histogram(X4,30,'Normalization','pdf', 'FaceColor',[142/255 160/255 199/255]);  % 绘制直方图
hold on;
plot(x4c,f_kx4,'--','Color',[0 0 0],'LineWidth',6); % 绘制核密度估计函数图
xlabel('Residual(mm)', 'FontName', 'Times New Roman','FontSize', 56);
ylabel('Frequency/group spacing','FontName', 'Times New Roman','FontSize', 56);
legend('Frequency distribution histogram','Kernel density function {\itf}_4({\itx}_4)','Location','northeast','FontSize',40, 'FontName', 'Times New Roman');
set(gca, 'FontSize', 50,'LineWidth', 2,'Fontname', 'Times new Roman','LooseInset', [0,0,0,0]);


%  6、根据已求的核密度函数、核分布函数，求残差数据X1取任一点xiq处的核密度函数值、核分布函数值，将各个测点的核分布函数值联立为随机向量，方便后续投喂到已求解出未知参数的多维Copula函数中

%求残差数据Xi取任一点xiq处的核密度函数值
f_kx1q=interp1(x1c,f_kx1,X1);%interp1函数用于求出变量X1取任一点处的核密度函数值f_kx1q
f_kx2q=interp1(x2c,f_kx2,X2);
f_kx3q=interp1(x3c,f_kx3,X3);
f_kx4q=interp1(x4c,f_kx4,X4);



%求残差数据Xi取任一点xiq处的核分布函数值

F_kx1q=interp1(x1c,F_kx1,X1,'previous');%interp1函数用于求出变量X1取任一点处的核分布函数值F_kx1q，其中的'previous'参数表示取离x1q最近的一个数。
F_kx2q=interp1(x2c,F_kx2,X2,'previous');
F_kx3q=interp1(x3c,F_kx3,X3,'previous');
F_kx4q=interp1(x4c,F_kx4,X4,'previous');


%将各个测点的核分布函数值联立成随机向量，方便后续投喂到已求解出未知参数的多维Copula函数中
T=[F_kx1q ,F_kx2q, F_kx3q, F_kx4q]; %T为各个测点的核分布函数值联立成的随机向量

