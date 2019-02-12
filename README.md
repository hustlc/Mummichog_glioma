# Mummichog_glioma
Mummichog is used for metabolomics data processing
为了研究胶质瘤瘤内代谢的异质性程度，我们把一块胶质瘤分为6块，加上癌旁AP，每个病人相当于做7份样本的非靶向质谱，15个病人一共做了105份样本；质谱和metDNA分析由中科院IRCBC（没错跟工商银行只差一个字）朱正江课题组完成，经过metDNA(v1.01)分析过的数据如下，NEG代表阴离子峰（39998个，图1），POS代表阳离子峰（23753个，图2），去掉同位素峰之后，分别剩余4678个（实际4675个）和2319个特征峰。X01.01-14为同一病人6块癌组织的特征峰强度，AP为癌旁特征峰强度，qc是对照样本（抽取各个标本的pooling）重复检测得到的一系列强度值

 
图1
 
图2

为了精简特征峰，我们将阳离子峰和阴离子峰进行合并，对于其中70个重合的特征峰，多数选择了阴离子，舍弃了阳离子；除非阳离子强度为阴离子的5倍以上，才取阳离子峰，否则取阴离子峰。最终特征峰合并为6924个，用于肿瘤整体代谢通路富集，即combine分析/全部阴离子/阳离子峰分析。分析单个癌组织的时候，因为脚本原因要暂时分开阴离子/阳离子峰，此时无法回溯哪个特征峰是阴/阳离子峰，所以去掉了这70个重合的特征峰，最终用于单个癌组织mummichog分析的特征峰只有6854个。

肿瘤整体代谢通路富集：
6854个特征峰，4605个是阴离子峰，2249个是阳离子峰，6个部分的癌组织和1个癌旁做单样本t检验，R脚本如下：

args  <- commandArgs(TRUE)
 input <- args[1]
 out   <- args[2]

dat = read.table(input, header=T)
dat$pvalue = 0
dat$ttest = 0
for (i in 1:nrow(dat))
{
   a=t.test(dat[i, 2:7], mu=dat[i, 8]);
   dat[i, ]$pvalue = a$p.value;
   dat[i, ]$ttest = a$statistic;
}
dat$qvalue = p.adjust(dat$pvalue,method="fdr")
write.table(x=dat, file=out, quote=F, row.names=F, sep="\t")

以上得到的特征峰p值和t值，连同mz/rt一起放入mummichog计算

单个肿瘤组织代谢通路富集（瘤内代谢异质性）：
6854个特征峰，4605个是阴离子峰，2249个是阳离子峰，6个部分的癌组织，每个都要同1个癌旁做t检验，然而两个数字是无法做t检验的啊（质谱应当做3个技术重复的，我们没做），于是引入了qc来计算qc的sd值，再用这sd值使每个峰的强度值变成三个，最后再跟癌旁进行单样本t检验，处理较复杂，请吴新贵博士帮我写了个perl脚本，吴博士是个老派码农，所以喜欢perl：

use warnings;
use strict;

my $first=$ARGV[0];
my $outfile=$ARGV[1];
`Rscript /public2/config/licong//SD-15.R $first tmp`;

open RF,"tmp";
open WF,">tmp2";
while(my $line=<RF>)
{
   chomp($line);
   my @all=split(/\t/,$line);
     if($.==1)
  {
    print WF "$all[0]\t$all[1]\t$all[2]-SD\t$all[2]\t$all[2]+SD\t$all[3]-SD\t$all[3]\t$all[3]+SD\t$all[4]-SD\t$all[4]\t$all[4]+SD\t$all[5]-SD\t$all[5]\t$all[5]+SD\t$all[6]-SD\t$all[6]\t$all[6]+SD\t$all[7]-SD\t$all[7]\t$all[7]+SD\t$all[8]\t$all[9]\t$all[10]\t$all[11]\t$all[12]\t$all[13]\t$all[14]\t$all[15]\t$all[16]\t$all[17]\t$all[18]\t$all[19]\t$all[20]\t$all[21]\t$all[22]\t$all[23]\tSD\n";
  }
  else
  {
    print WF "$all[0]\t$all[1]\t";
   for(my $a=2;$a<8;$a++) 
  {
    my $c=$all[$a]-$all[$#all];
    my $b=$all[$a]+$all[$#all];
    print WF "$c\t$all[$a]\t$b\t";
  }
   print WF "$all[8]\t$all[9]\t$all[10]\t$all[11]\t$all[12]\t$all[13]\t$all[14]\t$all[15]\t$all[16]\t$all[17]\t$all[18]\t$all[19]\t$all[20]\t$all[21]\t$all[22]\t$all[23]\t$all[$#all]\n";
  }
}
#`Rscript /public2/config/licong/ttest.R tmp2 >$outfile`;
#`rm tmp tmp2`;
close RF;
close WF;
`Rscript /public2/config/licong/ttest.R tmp2 $outfile`;
`rm tmp tmp2`;

单个肿瘤组织代谢通路富集（瘤内代谢异质性V2）：
这个跟前面有什么区别呢？我们把15个癌旁组织每个峰的强度取中位值（或者平均值），形成2个超级癌旁，每个峰的强度值（三个），跟超级癌旁进行单样本t检验

得到了所有的P值和t值之后，就可以用Mummichog来富集通路啦~ Mummichog的使用可以参考大神Prof Shuzhao Li在mummichog主页上的说明：
python /home/guests/libo/licong/python-packages/mummichog/main.py -f mummi_01_neg.txt -o myoutput

Example of input file (a test data file, right click to save):
            mz      rtime   p-value t-score
            186.0185697     463     0.000149751400132       3.82
            279.1773473     90      0.000399613326314       3.56
            344.1330624     124     0.000998323061251       -3.31
            215.9641894     132     0.00105418285794        -3.29
            177.0323244     77     0.00121065359218        3.256
            296.0973768     135     0.00171645907855        -3.15
            527.3784209     593     0.00176815004959        -3.14
在t-score的右边可以加上features，也可以不加（我自己加了）：
 

mummichog 2.0版本运行如果出现如下错误：
_tkinter.TclError: no display name and no $DISPLAY environment variable
需要自己在get_user_data.py这个python脚本里面加入下列内容，可以解决报错，无法作图的情况：

# Force matplotlib to not use any Xwindows backend.
matplotlib.use('Agg')
import matplotlib.pyplot as plt

最后提一下，mummichog无论是pathway analysis还是modular analysis，都要模拟100次（默认），因此，一模一样的输入，重复进行的话，结果也会略有不同，但差别不会太大，富集到的通路P值会不一样，不用惊慌…….一开始我还以为不显著的P值也会参与计算，其实不然，以下是大神Prof Shuzhao Li拨冗给我回的email：
Hi Chong, mummichog is independent from your statistical selection of features. You choose how to calculate a p-value per feature prior to mummichog. The "t score" column can be anything.

总结：不会代码的话，最多玩玩一个癌，加一个癌旁，千万别玩本鸟这样上百个样本的大规模组学分析了………
在此感谢中科院生物与化学交叉研0究中心朱正江研究员的博士生申小涛（群主大大）全程指导，边靖文同学完成质谱分析，感谢中山大学医学院李隽教授课题组的吴新贵博士全程编写代码，感谢中山大学医学院李淼新教授课题组代码小牛薛超同学帮助修正mummichog的源代码。
