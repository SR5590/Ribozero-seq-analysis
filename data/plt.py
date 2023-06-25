import seaborn as sns; sns.set()
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from glob import glob
import sys
import io

def fastqcDF(file):
    with open(file) as s:
        string=s.read()
        start=string.find('#Base	Mean	Median	Lower Quartile	Upper Quartile	10th Percentile	90th Percentile')
        string=string[start:]
        end=string.find('>>END_MODULE')
        print(start, end,sep='\t')
        df=string[:end]
        print(df)
        f = io.StringIO(df)
        df=pd.read_csv(f,sep='\t')
    return df


def TPM_FPK_Calculator(fcounttsv):
    df =pd.read_csv(fcounttsv,sep='\t',skiprows=1)
    df2=df[['Geneid','Length']+list(df)[6:]]
    df2.columns=['Geneid','Length',]+[x.split('/')[-1] for x in list(df2)[2:]]
    dfFPKM=df2.copy()
    dfTPM=df2.copy()
    for _ in list(dfFPKM)[2:]:
        dfFPKM[_]=(dfFPKM[_]/((dfFPKM[_].sum()/1_000_000))/(dfFPKM['Length']/1000))

    for _ in list(dfTPM)[2:]:
        millFac=(dfTPM[_]/(dfTPM['Length']/1000)).sum()/1_000_000
        dfTPM[_]=(dfTPM[_]/(dfTPM['Length']/1000))/millFac
    return dfFPKM,dfTPM


ID=sys.argv[1]
#ID='AN14016369'

fig = plt.figure(constrained_layout=False, figsize=(10,15))
fig.suptitle(ID+' - mRNA-seq Quality')
plt.subplots_adjust(top=0.9,hspace = 1.75, wspace=0.5)

gs = fig.add_gridspec(8, 2)
ax1=fig.add_subplot(gs[0:3,0:1])
ax2 = fig.add_subplot(gs[0:3,1:2])
#ax3 = fig.add_subplot(gs[3:6,0:1])
ax5 = fig.add_subplot(gs[3:6,0:2])
ax4 = fig.add_subplot(gs[7:8,:])


#dfR=fastqcDF('RNA_comparison/02_fast_qc/AN127467_2P_fastqc/fastqc_data.txt')
#dfF=fastqcDF('RNA_comparison/02_fast_qc/AN127467_1P_fastqc/fastqc_data.txt')

forward=glob(f'../01_raw/{ID}/fastqc/{ID}_1P_fastqc/fastqc_data.txt')[0]
reverse=glob(f'../01_raw/{ID}/fastqc/{ID}_2P_fastqc/fastqc_data.txt')[0]

dfR=fastqcDF(reverse)
dfF=fastqcDF(forward)

medians = dfF['Mean']
quartiles1 = dfF['Lower Quartile']
quartiles3 = dfF['Upper Quartile']
ten_per = dfF['10th Percentile']
ninty_per = dfF['90th Percentile']
x=dfF['#Base'].astype(str).str.split('-',expand=True)
x=x[0].astype(int)
sns.lineplot(x=x, y=medians, legend='full', ax=ax1)
ax1.set(ylim=(0, 40), title='phred scrore distribution - Forward', xlabel='position', ylabel='phred')
ax1.fill_between(x, ten_per, ninty_per, alpha=0.3, facecolor='red')
ax1.fill_between(x, quartiles1, quartiles3, alpha=0.3, facecolor='green')
ax1.legend(labels=['mean', '10/90 percentile', 'quartile'])


medians = dfR['Mean']
quartiles1 = dfR['Lower Quartile']
quartiles3 = dfR['Upper Quartile']
ten_per = dfR['10th Percentile']
ninty_per = dfR['90th Percentile']
x=dfR['#Base'].astype(str).str.split('-',expand=True)
x=x[0].astype(int)
sns.lineplot(x=x, y=medians, legend='full', ax=ax2)
ax2.set(ylim=(0, 40), title='phred scrore distribution - Reverse', xlabel='position', ylabel='phred')
ax2.fill_between(x, ten_per, ninty_per, alpha=0.3, facecolor='red')
ax2.fill_between(x, quartiles1, quartiles3, alpha=0.3, facecolor='green')
ax2.legend(labels=['mean', '10/90 percentile', 'quartile'])




x=f'../05_FeatureCount/{ID}_FeatureCountTable'

df=pd.read_csv('BioType_HomoSapiensGeneSymbols.tsv',sep='\t')
#df_quant=pd.read_csv(x, sep='\t',skiprows=1)
#df_quant=df_quant.iloc[:,[0,-1]]

df_fpkm,df_tpm=TPM_FPK_Calculator(x)
df_quant=df_tpm
df_quant.columns=['Geneid','Length',ID]
print(df_quant)
c=pd.merge(left=df,right=df_quant,left_on='id',right_on='Geneid')
print(c)
df_quant=c.groupby('type').sum()
df_quant=df_quant.sort_values(ID,ascending=False)
print(df_quant)
df_quant[ID].plot(kind='bar',ax=ax5)
ax5.set(ylabel='kallisto - TPM',ylim=(0,1000000),title='TPM distribution')



tableList=[]
tableListList=[]
df=pd.read_csv(f'../03_star/{ID}/UmiReads/{ID}Log.final.out',sep='\t',names=[1,2],index_col=0)

for xx in ['                          Number of input reads |',
          '                   Uniquely mapped reads number |',
          '                        Uniquely mapped reads % |',
          '                          Average mapped length |',
          '        Number of reads mapped to multiple loci |',
          '             % of reads mapped to multiple loci |',]:
    tableList.append(df.at[xx,2])
tableListList.append(tableList)

c=ax4.table(cellText=tableListList,loc='center',colLabels=['Number of \n input reads',
          'Uniquely mapped \n reads number',
          'Uniquely mapped \n reads %',
          'Average mapped \n length',
          'Number of reads \n mapped to multiple loci',
          '% of reads mapped \n to multiple loci',])
c.auto_set_font_size(False)
c.set_fontsize(10)
c.scale(1.2,2)
ax4.axis("off")
ax4.set(title='Overview Table')


plt.savefig(f'../06_quality/{ID}.pdf',bbox_inches='tight')
