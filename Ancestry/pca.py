import pandas as pd
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
#读取文件
pca = pd.read_csv(sys.argv[1]," ",header=None)
pca = pca.rename(columns=dict([(0,"Superpopulation"),(1,"Individual ID")]+[(x,"PC"+str(x-1)) for x in range(2,22)]))
pca_val = pd.read_csv(sys.argv[2]," ",header=None)
pca_ratio = round(pca_val*100/pca_val.iloc[:,0].sum(),2)

plt.figure(figsize=(20,15))
#sns_plot=sns.scatterplot(data=pca,x="PC1("+str(pca_ratio.iloc[0,0])+"%)",y="PC2("+str(pca_ratio.iloc[1,0])+"%)",hue="Superpopulation",s=50)
sns.set(font_scale=2)
sns.set_style("whitegrid")
sns_plot=sns.scatterplot(data=pca,x="PC1",y="PC2",hue="Superpopulation",s=250)
sns_plot.set_ylabel("Principal Component 2("+str(pca_ratio.iloc[1,0])+"%)", fontsize=30) #设置Y坐标轴标签字体
sns_plot.set_xlabel("Principal Component 1("+str(pca_ratio.iloc[0,0])+"%)", fontsize=30) #设置X坐标轴标签字体
sns_plot.legend(title="Superpopulation",fontsize = 20, title_fontsize = 30)
sns_plot.figure.savefig(sys.argv[3]+".pdf")
