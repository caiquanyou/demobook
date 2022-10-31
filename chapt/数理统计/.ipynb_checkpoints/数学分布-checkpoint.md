# 数学分布(汇总)
by hjfeng@<a href="https://github.com/fenghuijian" target="_blank"><img src="https://img.shields.io/github/stars/fenghuijian?style=social" style="padding-bottom:0.2em;" /></a> 

## 目的：

由于现在对于数学中随机变量的分布知道的太少，如果将来要自己写代码可能会受到影响。所以计划将之前遇到过的数据分布汇总，并作出详细的解释。**所以本文会出现大量的公式，希望可以做形象的解释，并联系到能解什么问题**。同时，希望本文可以持续的更新，学习更多有趣的东西。

**有机会在把分布对应的检验方法讲一下**

那么就直接进入主题。。。



**参考资料：**

概率论与数理统计（陈希孺）

Hand-book on STATISTICAL DISTRIBUTIONS for experimentalists

PROBABILITY AND MATHEMATICAL STATISTICS

Introduction to Probability



## 0. 基础概念：

### 0.1.离散型随机变量：

[参考网址](https://www.jianshu.com/p/b570b1ba92bb)

> 研究一个随机变量，不只是要看它能取哪些值，更重要的是它取各种值的概率如何！

* **概率函数**：就是用**函数的形式来表达概率。** 

  也有称为**概率质量函数**， 离散随机变量在各特定取值上的概率：

$$
p_i 
= P(X = a_i) \qquad (i = 1,2,3,4,5,6)  \tag{0.1.1}
$$

  这个简单的概率函数中，自变量 $X$ 是**随机变量**的**取值**，因变量$p_i$取值的**概率**。

* **概率分布**：顾名思义就是**概率的分布**

  以下的表可以叫作**离散型随机变量的概率分布**，严格来说叫做**离散型随机变量的值分布和值的概率分布列表**。列表中，**是把所有的可能出现的情况都列出来，包括取值与取值对应的概率**，也就是所谓的**概率分布。**

  | $X$   | $x_1$ | $x_2$ | $x_3$ | $\ldots$ | $x_n$ | $\ldots$ |
  | ----- | ----- | ----- | ----- | -------- | ----- | -------- |
  | $p_i$ | $p_1$ | $p_2$ | $p_3$ | $\ldots$ | $p_n$ | $\ldots$ |

* **概率分布函数**：就是把概率函数的取值的累加结果，也叫做**累积概率函数**

  $P\{X = X_i\} = p_i, i = 1,2,3,\ldots$ ，则 ，$F(x) = P(X \leq x) = \sum_{x_i \leq x}\ P_i$。因为$F(x)$是$X\leq x$在的诸值$x_i$的概率之和，又$F(x)$称为**累积概率函数**。

### 0.2.连续型随机变量：

连续型随机变量也有**概率函数**和**概率分布函数**，

但是连续型随机变量的“概率函数“，换了名字，叫做**概率密度函数**。

**概率分布函数**也叫**累积概率函数**

![概率密度函数和累积概率函数](B:/Mathematical_Model/Mathematical_Distribution/pictures/p1.PNG)

### 0.3.统计分布的均值和方差的计算公式：

1. 均值（means）：用来描述样本集合的中间点

   离散型随机变量：设$P(x)$是**离散型概率质量函数**，其均值被定义为：
   
   $$
   {\Bbb E}(x)
   = \sum_{k = 1}^\infty x_kp_k \tag{0.3.1}
   $$
   
   连续型随机变量：设$P(x)$是**连续型概率密度函数**，其均值被定义为：
   
   $$
   {\Bbb E}(x)
   = \int_{-\infty}^{+\infty} xp(x) \, {\rm d}x \tag{0.3.2}
   $$
   
   性质：1，线性运算；2，函数期望；3，乘积期望


2. 方差（variance）：用来描述样本的离散程度。

   方差是一种特殊的期望：
   
   $$
   {\Bbb D}(x)
   = {\Bbb E} ((x - {\Bbb E}(x))^2) \tag{0.3.3}
   $$
   
   离散型方差：
   
   $$
  \begin{align}
   {\Bbb D}(x) &=& {\Bbb E}((x - {\Bbb E}(x))^2) \\
   &=& \sum_{i = 1}^n (x_i - {\Bbb E}(x))^2 P(x_i) \\
   &=& \sum_{i = 1}^n (P_i \cdot x_i^2) - {\Bbb E}^2(x) \\
   &=& {\Bbb E}(x^2) - {\Bbb E}^2(x) \tag{0.3.4}
   \end{align} 
   $$
   
   连续型方差：
   
   $$
  \begin{align}
   {\Bbb D}(x) &=& \sigma \\
   &=& \int_{-\infty}^{+\infty} ((x - {\Bbb E}(x))^2p(x) \, {\rm d}x) \\
   &=& \int_{-\infty}^{+\infty} x^2p(x) \, {\rm d}x - {\Bbb E}^2 (x) \\
   &=& {\Bbb E}(x^2) - {\Bbb E}^2(x) \tag{0.3.5}
   \end{align} 
   $$
   
   两种不同随机变量的计算方法基本上是一样。接下来我们推到方差的计算公式：
    
   $$
   \begin{align}
   {\Bbb D}(X) &=& {\Bbb E}\{(X - {\Bbb E}(X))^2\} \\
   &=& {\Bbb E}\{X^2 - 2X{\Bbb E}(X) + [{\Bbb E}(X)]^2\} \\
   &=& {\Bbb E}(X^2) - {\Bbb E}[2X{\Bbb E}(X)] + 
   {\Bbb E}[{\Bbb E}^2(X)] \\
   &=& {\Bbb E}(X^2) - 2{\Bbb E}(X){\Bbb E}(X) + {\Bbb E}^2(X) \\
   &=& {\Bbb E}(X^2) - {\Bbb E}^2(X) \tag{0.3.6}
   \end{align} 
   $$

3. 待续

### 0.4.多元高斯分布的小结：

#### 0.4.1.标准高斯分布（一元）：Standardized Gaussian Distribution

**概率密度函数：**

$$
f(x) 
= \frac{1}{\sqrt{2\pi}} e^{-\frac{x^2}{2}} \tag{0.4.1}
$$

#### 0.4.2.一般高斯分布（一元）：Gaussian Distribution

**概率密度函数：**

$$
f(x) 
= \frac{1}{\sqrt{2\pi} \cdot \sigma} e^{-\frac{(x-\mu)^2}{2\sigma^2}} \tag{0.4.2}
$$

我们可以利用z-score的方法，一般的高斯分布进行标准化，$z = \frac{x - \mu}{\sigma}$ 这种标准化的方法。可以将原始的分布转化为标准正态分布。



#### 0.4.3.独立分多元高斯分布：Multivariate Gaussian Distribution

**推导公式的过程：理解（在机器学习中，这个多元的高斯分布是很重要的，以这个多元高斯分布理解，有助于我们理解其他类型的多元分布。**

1. 先假$n$个变量 $\vec{x} = [x_1, x_2, x_3,\cdots,x_n ]^T$他们的分量各自**互不相关**， 而且每个分量服从正态分布的（维度不相关的多元正态分布），各个维度的均值${\Bbb E}(x) = [\mu_1,\mu_2,\mu_3,\cdots,\mu_n]^T$。标准差为$\sigma(x) = [\sigma_1,\sigma_2,\sigma_3,\cdots,\sigma_n]^T$。

   

2. 我们再根据**联合概率密度公式**

   $$
   \begin{eqnarray}
   f(x) 
   &=& p(x_1, x_2, x_3, \cdots, x_n) \\
   &=& p(x_1)p(x_2)p(x_3) \cdots p(x_n) \\
   &=& 
   \frac{1}{(\sqrt{2\pi})^n\sigma_1\sigma_2\sigma_3\cdots\sigma_n}
   e^{{-\frac{(x_1-\mu_1)^2}{2\sigma_1^2}}{-\frac{(x_2-\mu_2)^2}{2\sigma_2^2}}{-\frac{(x_3-\mu_3)^2}{2\sigma_3^2}}\cdots {-\frac{(x_n-\mu_n)^2}{2\sigma_n^2}}} \tag{0.4.3}
   \end{eqnarray} 
   $$
   
   当我们令：
   
   $$
   \begin{eqnarray}
   z^2 &=& {\frac{(x_1-\mu_1)^2}{2\sigma_1^2}} + {\frac{(x_2-\mu_2)^2}{2\sigma_2^2}} + {\frac{(x_3-\mu_3)^2}{2\sigma_3^2}}\cdots + {\frac{(x_n-\mu_n)^2}{2\sigma_n^2}} \\
   \sigma_z &=& \sigma_1\sigma_2\sigma_3\cdots\sigma_n \tag{0.4.4}
   \end{eqnarray}
   $$
   
   这样的多元正态分布又可以写成一元那种漂亮的形式
   
   $$
   f(z) = \frac{1}{(\sqrt{2\pi})^n\sigma_z}e^{-\frac{z^2}{2}} \tag{0.4.5}
   $$
   
   因为多元正态分布有着很强的**几何思想**，代数的角度我们很难看出$z$的概率分布规律，但是转化为矩阵的形式就可能比较好一点：
   
   $$
   \begin{eqnarray}
   z^2 &=& z^Tz  \tag{0.4.6} \\
   &=& 
   [x_1-\mu_1, x_2-\mu_2, x_3-\mu_3, \cdots, x_n-\mu_n] \begin{bmatrix} 
   \frac{1}{\sigma_1^2} & 0 & \cdots & 0 \\
   0 & \frac{1}{\sigma_2^2} & 0 & 0 \\
   0 & 0 & \frac{1}{\sigma_3^2} & 0 \\
   \vdots & \cdots & \cdots &\vdots \\
   0 & 0 & 0 & \frac{1}{\sigma_n^2}
   \end{bmatrix}
   [x_1-\mu_1, x_2-\mu_2, x_3-\mu_3, \cdots, x_n-\mu_n] ^T
   \end{eqnarray}
   $$
   
   这样的公式比较长，我们继续做一个变量的替换，
   
   $$
   \vec{x} - \vec{\mu_x} = [x_1-\mu_1, x_2-\mu_2,x_3-\mu_3,\cdots,x_n-\mu_n]^T \tag{0.4.7}
   $$
   
   定义一个符号$\Sigma$ , 它为多维变量 $\vec{x}$ 的**协方差矩阵**。因为$\vec{x}$各个分量是相互独立的所以，除了对角元素之外，其他的值都是为0；对角元素为$x_i$的方差。矩阵的i行j列的元素值表示$x_i$ 和 $x_j$ 的**协方差**。协方差是用来描述两个变量之间的相关程度，这个概念有这样的解释。
   
   $$
   \Sigma =
      \begin{bmatrix} 
      \sigma_1^2 & 0 & \cdots & 0 \\
      0 & \sigma_2^2 & 0 & 0 \\
      0 & 0 & \sigma_3^2 & 0 \\
      \vdots & \cdots & \cdots &\vdots \\
      0 & 0 & 0 & \sigma_n^2 
      \end{bmatrix}
   $$     
   
   那么它对应的**逆矩阵**为：
   
   $$
   \Sigma^{-1} = 
      \begin{bmatrix} 
      \frac{1}{\sigma_1^2} & 0 & \cdots & 0 \\
      0 & \frac{1}{\sigma_2^2} & 0 & 0 \\
      0 & 0 & \frac{1}{\sigma_3^2} & 0 \\
      \vdots & \cdots & \cdots &\vdots \\
      0 & 0 & 0 & \frac{1}{\sigma_n^2} 
      \end{bmatrix}
   $$
   

   对角矩阵的**行列式** = **对角元素**的乘积：
   
   $$
   \mid \Sigma \mid 
      = \sigma_1^2 \sigma_2^2 \sigma_3^2 \cdots \sigma_n^2 \tag{0.4.10}
   $$
   
   所以有：
   
   $$
   \sigma_z 
      = \mid \Sigma \mid^{\frac{1}{2}} 
      = \sqrt{\mid \Sigma \mid}
      = \sqrt{\sigma_1^2 \sigma_2^2 \sigma_3 ^2 \cdots \sigma_n^2}
      = \sigma_1 \sigma_2 \sigma_3 \cdots \sigma_n \tag{0.4.11}
   $$
   
   利用上述公式$(0.4.5), (0.4.6),(0.4.11)$我们可以得到
   
   $$
   \begin{eqnarray}
      z^Tz &=& (\vec{x} - \vec{\mu_x})^T \Sigma^{-1} (\vec{x} - \vec{\mu_x}) \tag{0.4.6} \\
      \sigma_z &=& \mid \Sigma \mid^{\frac{1}{2}} \tag{0.4.11} \\
      \\
      f(z) &=& \frac{1}{(\sqrt{2\pi})^n\sigma_z}e^{-\frac{z^2}{2}} \tag{0.4.5} \\
      \downarrow &\downarrow& \downarrow \\
      f(\vec{x}) &=& \frac{1}{(\sqrt{2\pi})^n\mid\Sigma\mid^{\frac{1}{2}}}
      e^{-\frac{(\vec{x} - \vec{\mu_x})^T \Sigma^{-1} (\vec{x} - \vec{\mu_x})}{2}}  \tag{0.4.12}
      \end{eqnarray}
   $$
   
   所以，现在的$f(\vec{x})$ 为独立多元高斯分布的公式。

   注意一下前面的系数变化：

   1、从非标准正态分布 $\rightarrow$标准正态分布，需要**概率密度函数**的高度压缩$\mid \Sigma \mid^{\frac{1}{2}} $倍。

    2、从一维 $\rightarrow$$n$维，**每增加一维**，高度将被压缩$\sqrt{2\pi}$倍。

   **前提条件是$\vec{x}$的各个分量是相互独立的。**

   

3. **有相关性的多元高斯分布**

   前面的多元正态分布的前提是多元变量之间是相互独立的，事实上，在很多的情况下，**变量于变量之间是有关联的**。

   当多元变量之间是由相关性的，那么可我们可以利用线性变换，将原来的变量转化为相互独立的新变量，一般来说变量是可以进行线性变换，从而将变量转化为相互独立的新变量。

   在**去相关性的推导**那部分已经说明了，相关性矩阵是可以进行decorrelation的操作的，而且是有解的。这个过程可以理解为，坐标轴的转换。

   

   以二维向量为例子。假设新坐标系为$x'_{1} =[u_{x1}^0, u_{x1}^1,]^T, \; x'_{2}=[u_{x2}^0, u_{x2}^1,]^T$那么原坐标系上的任意一点$[x_1, x_2]^T$投影到新的坐标系上的结果为：
   
   $$
   \begin{bmatrix}
   x'_{1} \\ 
   x'_{2} \\ 
   \end{bmatrix} 
   =
   \begin{bmatrix}
   u_{x1}^0 & u_{x1}^1 \\ 
   u_{x2}^0 & u_{x2}^1 \\ 
   \end{bmatrix}
   \begin{bmatrix}
   x_1 \\ 
   x_2 \\ 
   \end{bmatrix}
   $$
   
   所以我们可以定义一个变换矩阵$U$:
   
   $$
   U = 
   \begin{bmatrix}
   u_{x1}^0 & u_{x2}^0 \\ 
   u_{x1}^1 & u_{x2}^1 \\ 
   \end{bmatrix}
   $$
   
   $U$的列空间有新坐标向量组成，坐标映射之后：$X' = U^{T}X$

   那么我们现在的$X'$就是相互独立的了，满足了维度不相关的高斯分布模型，我们就可以套用公式了：
   
   $$
   \begin{eqnarray}
   f(\vec{x}) &=& \frac{1}{(\sqrt{2\pi})^n\mid\Sigma\mid^{\frac{1}{2}}}
      e^{-\frac{(\vec{x} - \vec{\mu_x})^T \Sigma^{-1} (\vec{x} - \vec{\mu_x})}{2}} \tag{0.4.12}
   \end{eqnarray}
   $$
   
   $\vec{x} \to \vec{x'}$和$\vec{\mu_{x} }\to \vec{\mu_{x'}}$，这两个参数的转变是比较简单的，对$\Sigma'$是未知，它是$X'$的协方差矩阵，我们从定义出发，求解出$\Sigma'$:
   
   $$
   \begin{eqnarray}
   \vec{x} \to \vec{x'}, &\quad& \vec{x'} = U^{T}\vec{x} \tag{0.4.13} \\
   \vec{\mu_{x}} \to \vec{\mu_{x'}}, &\quad& \vec{\mu_{x'}} = {\Bbb E}(U^{T}\vec{x}) = U^T{\Bbb E}(\vec{x}) = U^{T}\vec{\mu_{x}} \tag{0.4.14}
   \end{eqnarray}
   $$
   
   $$
   \begin{eqnarray}
   {\Bbb D}(X') 
   &=& {\Bbb E}[(X' - \mu_{X'})(X' - \mu_{X'})^{T}] \\
   &=& {\Bbb E}[(X' - \mu_{X'})(X'^{T} - \mu_{X'}^{T})] \\
   &=& {\Bbb E}[X'X'^{T} - \mu_{X'}X'^{T} - X'\mu_{X'}^{T} + \mu_{X'}\mu_{X'}^{T}] \\
   &=& {\Bbb E}[U^{T}XX^{T}U - {\Bbb E}[U^{T}X]X^{T}U - 
   U^{T}X{\Bbb E}[U^{T}X]^{T} + {\Bbb E}[U^{T}X]{\Bbb E}[U^{T}X]^{T}] \\
   &=& {\Bbb E}[U^{T}XX^{T}U - U^{T}{\Bbb E}[X]X^{T}U - 
   U^{T}X{\Bbb E}[X]^{T}U + U^{T}{\Bbb E}[X]{\Bbb E}[X]^{T}U] \\
   &=& U^{T}{\Bbb E}[XX^{T} - {\Bbb E}[X]X^{T} - X{\Bbb E}[X]^{T} + {\Bbb E}[X]{\Bbb E}[X]^{T}]U \\
   &=& U^{T}{\Bbb D}(X) U  \tag{0.4.15}
   \end{eqnarray}
   $$
   
   所以我们可以由上述公式得到：经过坐标映射后的协方差矩阵是满足关系：

   我们先整理一下$U$这个变换矩阵的特性，首先他是一个**正交矩阵**，**列向量**为**单位向量**。所以会有$U^{T} = U^{-1}$这样的特性。
   
   $$
   \begin{eqnarray}
   (\Sigma)_{\vec{x'}} &=& U^{T}(\Sigma)_{\vec{x}}U \tag{0.4.16} \\
   (\Sigma)_{\vec{x'}} &=& U^{-1}(\Sigma)_{\vec{x}}U \tag{0.4.17} 
   \end{eqnarray}
   $$
   
   由$(0.4.17)$我们可以得到$(\Sigma)_{\vec{x'}}$和$(\Sigma)_{\vec{x}}$是**相似矩阵**，而相似矩阵的**行列式**是相等的。即：
   
   $$
   |(\Sigma)_{\vec{x'}}| = |(\Sigma)_{\vec{x}}| \tag{0.4.18}
   $$
   
   另外，我们还可以推导出：
   
   $$
   \begin{eqnarray}
   (\Sigma)_{\vec{x'}}^{-1}
   &=& (U^{T}(\Sigma)_{\vec{x}}U)^{-1} \\
   &=& (U^{-1}(\Sigma)_{\vec{x}}U)^{-1} \\
   &=& U^{-1}(\Sigma)_{\vec{x}}U \\
   &=& U^{T}(\Sigma)_{\vec{x}}U \tag{0.4.19}
   \end{eqnarray}
   $$
   
   有了上述的$(0.4.13, \ 0.4.14, \ 0.4.16,\ 0.4.17, \ , 0.4.18, \ , 0.4.19)$我们就可以套入上面的独立的高斯分布的公式$(0.4.12)$中了。
   
   $$
   \begin{eqnarray}
   f(\vec{x}) &=& 
   \frac{1}{(\sqrt{2\pi})^n\mid(\Sigma)_{\vec{x}}\mid^{\frac{1}{2}}}
   e^{-\frac{(\vec{x} - \vec{\mu_{x}})^{T}(\Sigma)_{x}^{-1} (\vec{x} - \vec{\mu_{x}})}{2}} \\
   \downarrow &\downarrow& \downarrow \\
   f(\vec{x'}) &=& 
   \frac{1}{(\sqrt{2\pi})^n\mid(\Sigma)_{\vec{x'}}\mid^{\frac{1}{2}}}
   e^{-\frac{(\vec{x'} - \vec{\mu_{x'}})^{T} (\Sigma)_{x'}^{-1} (\vec{x'} - \vec{\mu_{x'}})}{2}} \\
   f(\vec{x}) &=& 
   \frac{1}{(\sqrt{2\pi})^n\mid(\Sigma)_{\vec{x}}\mid^{\frac{1}{2}}}
   e^{-\frac{(U^{T}\vec{x} - U^{T}\vec{\mu_{x}})^{T} 
   U^{T}(\Sigma)_{x}^{-1}U (U^{T}\vec{x} - U^{T}\vec{\mu_{x}})}{2}} \\
   f(\vec{x}) &=& 
   \frac{1}{(\sqrt{2\pi})^n\mid(\Sigma)_{\vec{x}}\mid^{\frac{1}{2}}}
   e^{-\frac{(\vec{x} - \vec{\mu_{x}})^{T}
   UU^{T}(\Sigma)_{x}^{-1}U U^{T}(\vec{x} - \vec{\mu_{x}})}{2}} \\
   f(\vec{x}) &=& 
   \frac{1}{(\sqrt{2\pi})^n\mid(\Sigma)_{\vec{x}}\mid^{\frac{1}{2}}}
   e^{-\frac{(\vec{x} - \vec{\mu_{x}})^{T}(\Sigma)_{x}^{-1} (\vec{x} - \vec{\mu_{x}})}{2}} \tag{0.4.12}
   \end{eqnarray}
   $$
   
   显然，虽然我们通过了公式的推导，得到的结果依然是公式$(0.4.12)$，**所以多元高斯分布的一般形式是不区分当前的矩阵$X$是否独立**，它适用于所有的多元高斯分布。但是推导的原理**还是基于多元独立变量**推导的，这是**原则条件**。

   

   总结一下：

   1、我们先是定义一个新的坐标系，通过矩阵$U^{T}$的映射，将原来的的元素映射得到新的坐标系，**目的是去相关性**

   2、在新的坐标下，我们可以相应定义出新的期望，新的协方差，新的协方差的逆，它们都是可以通过$U$和$U^{T}$计算得到，当然我们是可以不用计算的。其实这个$U = E^{T}$

   3、套用公式，我们得到的结果，其实与$U，U^{T}$无关的。为什么会这样？

   其实可以很好理解，首先我们的条件：概率模型在矩阵$X$的条件下，已经存在的，只是我们不知道长什么样子。但是它是客观存在的。

   所以在这样一个分布下，对应的点，不管更换怎么样的坐标系，其概率都是固定的，而且只会受到原来数据的影响。

4. 待续

5. 待续

#### 0.5.去相关性的推导 decorrelation

在数据分析中，我们遇到的很多数据都基本上是有相关性的，但是这样的相关性有好有坏，好处是，这样的相关性是反映了数据的真实情况。而坏处是，有相关性的数据，不利于我们进行分析，个人观点，相关性的数据是会将数据分析变复杂的，而且不好利用模型去分析。

那么在这里简单学习一下去相关性 decorrelation的证明过程。**由于网上找到的资料十分有限，（网上的资料那是真的少），所以这个推导过程还不是十分的完善：**

对于一个$K$个维度，$n$个观测变量的矩阵$X$，它的协方差矩阵$\Sigma$如下：

$$
\begin{eqnarray}
\Sigma 
&=& Cov(X) \\ 
&=& {\Bbb E}(XX^T) \\
&\approx& \frac{XX^T}{n} \tag{0.5.1}
\end{eqnarray}
$$

然后对协方差矩阵做特征值分解：Eigenvalue Decomposition

$$
\begin{eqnarray}
\Sigma &=& EDE^{-1} \tag{0.5.2} \\
E^{-1}\Sigma E &=& D \tag{0.5.3} \\
\end{eqnarray}
$$

其中我们知道$E$为$[K \times K]$的正交矩阵，$D$是对角矩$[K \times K]$阵。

那么我们现在的目的是让$X$这个有相关性的矩阵，经过线性变换，转化为一个没有相关性的矩阵$Y$。**没有相关性矩阵$Y$的协方差矩阵是一个对角矩阵**，那么怎么找到这样一个线性变换$W_D$呢，这就要利用我们上面的三个公式的特点了$(0.5.1), \ (0.5.2), \ (0.5.3)$。

显然$D$它是一个对角矩阵，所以我们可以想办法构建一个$W_D$：

$$
\begin{eqnarray}
Y &=& W_D X \tag{0.5.4} \\
\\
D &=& Cov(Y) = {\Bbb E}(YY^T) \tag{0.5.5}) \\
\\
D &=& \frac{W_DX(W_DX)T}{n} \\
&=& W_D\Sigma W_D^T \tag{0.5.6}
\end{eqnarray}
$$

根据公式$(0.5.3)$, 同时我们知道

$$
\begin{eqnarray}
E^{-1} \Sigma E &= &W_D \Sigma W_D^T \tag{0.5.7} \\
\end{eqnarray}
$$

然后不知道怎么推导最后就得到了：

$$
W_D = E^T
$$

也就是说去相关性的那个变换矩阵就是$X$协方差矩阵的特征值分解得到的那个特征矩阵的转置。



#### 0.6. standard deviation and standard error of mean

这个内容是用来说一下standard deviation 和 standard error的区别，顺便吐槽一下，中文的翻译真的操蛋，还有一旦查询一些深一点的知识，中文资料是真的不可靠。

1. SD是一个比较清晰的概念，**它是用来描述总体（可以用来描述样本，其实在大数据上，个人认为总体和样本的界限分的不是很清晰）的离散程度**

   所以但是值得注意一下，自由度的$n-1$的问题，其实在大样本的情况下，这个问题就显得没那么重要了。

2. SE是一个比较蛋疼的概念，首先**它是用来描述样本均值的离散程度**，注意**是样本均值**。所以这样就与SD有了本质的区别了，**SD是针对总体的**，**SE是针对样本均值的总体**。

   计算的来源：

   首先我们要知道**抽样是怎么来的**。假设有总体mean为$\mu$，SD为$\sigma$。然后$X_1, X_2, X_3,\cdots, X_n$为我们从这个总体中随机抽样得到的$n$个独立的观测值，**而这些观测值又是随机变量，只是这样的随机变量是我们已经知道的。**

   而这些随机变量的SD为$\sigma$，另外我们可以计算随机变量之和的方差：

   $T = X_1+X_2+X_3+\cdots+X_n$的方差为$n\sigma^2$。

   而对于**抽样样本的均值**：$\frac{T}{n}$， 的方差为：$\frac{1}{n_2}n\sigma^2 = \frac{\sigma^2}{n}$ ； 对应的SD为$\frac{\sigma}{\sqrt{n}}$

   所以我们要理解SE，可以把他描述为对于**样本均值的标准差**，在做假设检验和区间估计我们处理的对象是**样本均值**，所以其对应的SD就是$\frac{\sigma}{\sqrt{n}}$，现在就对于SE有了一个比较清晰的理解。

   当然，我们是基本不知道总体mean为$\mu$，SD为$\sigma$。所以就会用样本均值的$\overline{x}$和样本SD的$s$来做估计，这个过程是不冲突的。

3. 首先我们，



---



## 1. 伯努利分布 Bernoulli Distribution

1. 简介：

   伯努利分布 是一种**离散分布**,有两种可能的结果。

   **1**表示成功，出现的概率为$p,\ 0<p<1$

   **0**表示失败，出现的概率为$q=1-p$。

2. **概率质量函数：**

   $$
   p(r;p) = \begin{cases}
   p&if&r=1&(success)\\
   1-p=q&if&r=0&(failure)
   \end{cases}
   $$

3. 应用：

   这种分布在人工智能里很有用，比如你问机器今天某飞机是否起飞了，它的回复就是Yes或No，非常明确，这个分布在分类算法里使用比较多，因此在这里先学习 一下。

4. 现在遇到有关**伯努利分布**的**生物学数据分析的文章**，在这里提一下**scVI**和**ZIFA**。

5. 待续

___



## 2. 二项分布 Binomial Distribution

1. 简介：

   二项分布的基本描述：**（离散型随机变量）**

   在每次独立试验中只有取两个值，表示成功的值的概率为$p$，那么表示试验不成功的概率为$1-p$。这样一种判断成功和失败的二值试验又叫做**伯努利试验**。

   在概率论和统计学里面，带有参数**$n$和$p$的二项分布**表示的是$n$次独立试验的成功次数的概率分布。

   &nbsp;

   **二项分布**频繁地用于对以下描述的一种实验进行建模：

   从**总数量大小为N**的两个事物中进行$n$次放回抽样，以**某一事物A**为基准，计算成功抽取**这个事物A**的**次数**的**概率**。要注意的是必须进行的是**放回抽样**，对于**不放回抽样我们一般用超几何分布来对这样的实验进行建模**。

2. 公式：

   一般来说，如果一个**随机变量X**满足二项分布的话，那么它一定有一个参数$n\inℕ$且还有一个参数$p\in [0,1]$。这样的话，我们可以把关于X的二项分布写成$X \thicksim B(n, p)$。

   &nbsp;

   **概率质量函数：**
   
   $$
   p(r;n,p)=p(X=r)= \begin{pmatrix}n\\ r\end{pmatrix}
   p^r(1-p)^{n-r} 
   $$
   
   其中，$r = 1,2,3,\ldots,n$，圆括号是组合的表示形式，组合的计算公式：
   
   $$
   \begin{pmatrix}n\\ r\end{pmatrix}=\frac{n!}{r!(n-r)!}
   $$
   
   **这个公式表示的从$n$个数中取$r$个数构成一个组合能有多少种不同的取法**。整个二项分布我们可以描述为求$n$次独立的**伯努利试验**，**成功r次的概率是多少**。

3. 二项分布的均值和方差公式：

   均值和方差：
   
   $$
   \begin{align}
   E(r) &=& np \tag{2.3}\\
   V(r) &=& np(1-p) \tag{2.4}
   \end{align}
   $$

4. 由于一次成功的概率是已知的，因此我们必须求出$n$次试验中，成功$r$次可能发生在哪次试验中，一共有多少种可能都要求出来，因此求的就是$n$取$k$组合数目。

5. 待续

___



## 3. 负二项分布 Negative Binomial Distribution

1. 简介：

   负二项分布也是有**伯努利试验**得到的。**（离散型随机变量）**

   伯努利试验哪（0-1，失败与成功）。所以对于描述**成功**和描述**失败**的时候可以对**负二项分布**提供两个方向的解释。

2. 概率质量函数：对于**成功**的描述，**其变量是$r$, $r$为试验的总数**。

   $$
   \begin{eqnarray}
   P(X = r | k,p) &=& p(r;k,p)\\
   &=& \begin{pmatrix}r-1\\ k-1\end{pmatrix} p^k(1-p)^{r-k} \tag{3.1}
   \end{eqnarray}
   $$
   
   这个概率质量函数的意思是：**进行了$r$伯努利试验的尝试，终于有$k$次实验是成功的, $p$是成功事件的概率.**

   **根据这样的解释，我们不难发现，第$k$次是必定发生成功事件的，在$r-1$次做$k-1$次的组合**

   所以可以把上面的公式还原最初的来源的形式：
   
   $$
   \begin{eqnarray}
   p(r;k,p)&=&\begin{pmatrix}r-1\\ k-1\end{pmatrix}p^k(1-p)^{r-k}\\
   &=&\begin{pmatrix}r-1\\ k-1\end{pmatrix}p^{k-1}(1-p)^{(r-1)-(k-1)}p \tag{3.2}
   \end{eqnarray}
   $$

3. 概率质量函数：对于**失败**的描述， **其变量为$n$， $n$为失败的次数**

   这个概率质量函数的意思是：**进行了$r$伯努利试验的尝试，终于有$k$次实验是成功的, $p$是成功事件的概率.**

   **根据这样的解释，我们不难发现，第$k$次是必定发生成功事件的，令在这之前经历过的失的次数$n = r-k$**。
   
   $$
   \begin{eqnarray}
   P(X = n | k,p) &=& p(n;k,p) \\
   &=& \begin{pmatrix}n+k-1\\ n\end{pmatrix}p^k(1-p)^{n} \tag{3.3}
   \end{eqnarray}
   $$

4. **对于负二项分布，在生物数据中是时常用到的，但是，为什么生物学数据，如RNA-seq就是符合负二项分布呢？这是我现在不太清楚的。**

5. 负二项分布的均值和方差公式：

   对于公式$(3.1)$最后只与成功的次数$k$和成功的概率$p$有关系。
   
   $$
   \begin{eqnarray}
   E(r) &=& \frac kp \tag{3.4} \\
   D(r) &=& \frac{k(1-p)}{p^2} \tag{3.5}
   \end{eqnarray}
   $$
   

   对于公式$(3.3)$均值公式不一样，方差和$(3.1)$一样。
   
   $$
   \begin{eqnarray}
   E(n) &=& E(r) - k = \frac{k(1-p)}{p} \tag{3.6} \\
   D(n) &=& \frac{k(1-p)}{p^2} \tag{3.7}
   \end{eqnarray}
   $$

6. 待续

___



## 4. 泊松分布 Poisson Distribution

1. 简介：

   泊松分布可以用来描述某段时间内，事件具体发生的概率。它是一种**离散随机变量**。

   泊松分布描述的是一个**离散随机变量**在**单位时间内**发生的**次数**的**概率**。

2. 概率质量函数：

   $$
   p(x = r; \mu) = p(r;\mu)=\frac{\mu^re^{-\mu}}{r!} \tag{4.1}
   $$
   
   这个公式符号代表的含义是，

   $\mu$的意思是我们统计得到的**已知单位时间内发生某事件的平均次数**；

   $r$的意思是我们在**一次单位时间内，某事件发生$r$次**

   $p(r;\mu)$是 **一次单位时间内，某事件发生$r$次，的概率有多大。**

3. 扩展：

   这里对泊松分布做一定的扩展。**泊松分布可以用来描述某段时间内，**

   **事件具体发生的概率**

   具体公式：$\mu$上加上事件可以用来描述某段时间，事件的发生次数的概率。
   
   $$
   p(r;t,\mu)=\frac{(\mu t)^re^{-\mu t}}{r!} \tag{4.2}
   $$

4. 与二项分布的联系：

   [与二项分布的联系:](#二项分布 Binomial Distribution)

   泊松分布可以看作二项分布的极限得到的。这里要满足一些条件。

   * 二项分布成功事件的概率$p$要很小；
   * 二项分布试验次数$N$要非常大。

   ![poisson and binomial](B:/Mathematical_Model/Mathematical_Distribution/pictures/p2.PNG)

   &nbsp;

   说到泊松分布，最好是明白：泊松分布是二项分布$n$很大而$p$很小时的一种极限形式。

   [二项分布](#二项分布 Binomial Distribution)：已知某件事情发生的概率是$p$，那么做$n$次试验，事情发生的次数就服从二项分布。

   **泊松分布式某段连续的时间内事情发生的次数。事情发生的时间是可以忽略的。关注的是事件的发生。泊松分布是离散的变量.**
   这段时间是确定大小的，不是说某两件事件（不知何时发生）的间隔。

   **把连续的时间分割层无数小份，那么每个小份之间都是相互独立的。在每个很小的时间区间内，事情可能发生也可能不发生，因此这就是一个$p$很小的二项分布。连续的时间分成无数小份，也就意味着$n$很大，即：泊松分布是二项分布的一种极限形式。**

   **此外，二项分布是最简单的发生于不发生的分布，那么与此关系密切的泊松分布自然在生活中很常见也可以理解了。**

   泊松分布中的$\mu$意义就是：一个时间段内时间平均发生的次数。

   [指数分布](#指数分布 Exponential Distribution)是两件事情发生的平均间隔时间，时间是连续变量。

   &nbsp;

5. **泊松分布均值与方差的的公式：**

   [参考的证明](https://blog.csdn.net/saltriver/article/details/52969014)
   
   $$
   \begin{eqnarray}
   E(r) &= \mu \tag{4.3} \\
   D(r) &= \mu \tag{4.4}
   \end{eqnarray}
   $$

6. **对与这个分布是十分重要的。因为它与时间挂钩了，而且它可以应用到RNA-seq的数据分布上，至于为什么可以用到RNA-seq上，这不太清楚**

7. 待续

___



## 5. 指数分布 Exponential Distribution

指数分布是一种**连续型概率分布**， 它主要应用在**随机事件之间发生的时间间隔**的**概率**问题。

[泊松分布](#泊松分布 Poisson Distribution)是描述某一个时间内发生随机事件次数的概率，而**指数分布**是描述**两次随机事件发生时间间隔的概率分布**。如旅客进机场的时间间隔，中文维基百科新条目出现的时间间隔。

**概率密度函数：**

$$
f(x; \alpha) = \frac{1}{\alpha}e^{-\frac{x}{\alpha}} \tag{5.1}
$$



___



## 6. 伽马分布 Gamma Distribution

从意义上看，

指数分布解决的问题是 ”要等到一个随机事件的发生，要经历多久时间。“

伽马分布解决的问题是 “要等到第n个随机事件的发生，要经历多久时间。”

伽马分布是一种**连续型概率分布**；

**概率密度函数：**

$$
f(x;a,b) = \frac{a(ax)^{b-1}e{-ax}}{\Gamma(b)} \tag{6.1}
$$


**gamma function：**，$\Gamma(x)$正实数的阶乘。

当如果$n \in \{ 1,2,3,...\}$，那么有：

$$
\Gamma(n) = (n - 1)! \tag{6.2}
$$

更一般的形式：对于一般的正实数：

$$
\Gamma(\alpha) = \int_0^{\infty} x^{\alpha - 1}e^{-x} {\rm d}x, \qquad for \quad \alpha > 0 \tag{6.3}
$$

可以把$\Gamma(\alpha)$理解为一个数的阶乘，但是在连续变量中，它会以一种积分的形式展现出来。这样我们可以比较好理解这样一个函数。

**理解Gamma distribution：**$\Gamma$分布可以理解为Possion分布在正实数集上的连续版本

回顾一下Possion 分布的扩展模式，符号可能有点变化：

$$
P(r) = \frac{(\lambda t)^re^{-\lambda t}}{r!} \tag{6.4}
$$

具体的意思是：在给定的$t$时间端内；$\lambda$为我们统计的在单位时间内平均发生的次数；$P(r)$是一个概率值，指得是在$t$时间端内，发生了$r$事件得概率。

**公式的推导：**

利用Possion 分布我们会得到$k: $ ~  事件发生的概率：

$$
\sum_{r = 0}^{k-1} P(r) = \sum_{r = 0}^{k-1} \frac{(\lambda t)^re^{-\lambda t}}{r!} \tag{6.5}
$$

那么在$t$时间段内，至少有$k$次事件发生的概率：

$$
\begin{eqnarray}
F(t) &=& \sum_{r = 0}^{\infty} P(r) \\
&=& 1 -\sum_{r = 0}^{k-1} \frac{(\lambda t)^re^{-\lambda t}}{r!} \\
&=& \int_0^{\lambda t} \frac{z^{k-1}e^{-z}}{(k - 1)!} {\rm d}z \\
&=& \int_0^t \frac{\lambda^k z^{k-1}e^{-\lambda z}}{(k - 1)!} {\rm d}z \tag{6.6}
\end{eqnarray}
$$

当$a = \lambda, \; b = k$，这个函数可以看作是Gamma分布的**累计概率函数**。所以我们可以看到出来Gamma分布是用来描述要发生第“k”次事件的时候，所需要的时间。

**Gamma 分布的均值和方差：**

$$
\begin{eqnarray}
E(x) &=& \frac{b}{a}, \tag{6.7} \\
V(x) &=& \frac{b}{a^2} \tag{6.8}
\end{eqnarray}
$$



___





## 7.待续

