#!/usr/bin/env python
import argparse


def readbed(bedfile):
    with open(bedfile,'r') as f:
        for line in f:
            line=line.strip().split('\t')
            yield line   #读取bed12文件，遍历每行，使用yield返回生成器对象
def getTranscript(bedfile,gtffile):
    with open(gtffile,'w') as w:
            for line in readbed(bedfile):
                chr_,start,end,name,a_,strand,b_,c_,rgb,d_,block,blocksit=line
                gene=name.split(";")[0]  #基因名称
                transcript=name.split(";")[1]  #转录本名称
                block=block.split(",") #返回一个包含每个exon长度的列表
                blocksit=blocksit.split(",")  #返回一个包含每个exon相对于第一个exon的起始位置的列表
                exonlist=[]
                #生成一个包含exon起始位点信息的列表
                for i in range(len(block)):
                    if i==0:
                        exonstart=int(start)
                    else:
                        exonstart=int(start)+int(blocksit[i])+1
                    exonend=int(start)+int(blocksit[i])+int(block[i])
                    exonlist.append((exonstart,exonend))
                #这里做一个if判断，判断来源，本例中一些转录本来自Iso-Seq，由红色RGB代表，一些来自ssRNA-Seq由橙色RGB代表，我们将具有多个来源的转录本默认为Iso-Seq来源
                if rgb=="255,100,0":
                    source="ssRNA-Seq"
                else:
                    source="Iso-Seq"
                w.write(f'{chr_}\t{source}\ttranscript\t{start}\t{end}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{transcript}\"\n')
                #写入GTF文件
                for exon in exonlist:
                    exec(r"w.write(f'{chr_}\t{source}\texon\t{exon[0]}\t{exon[1]}\t.\t{strand}\t.\tgene_id \"{gene}\"; transcript_id \"{transcript}\"\n')")

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('bedfile',type=str,help='输入bed12格式的文件')
    parser.add_argument('gtffile',type=str,help='输出GTF文件')
    args=parser.parse_args()
    getTranscript(args.bedfile,args.gtffile)


if __name__=="__main__":
    main()