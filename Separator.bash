#!/bin/bash
#Seperating Amplicons
echo ===============================================
echo ++++++++++++ Welcome to PopS 0.0.1 ++++++++++++
echo ++++++++++++ Designed by Ruoyao-Ni ++++++++++++
echo ===============================================
echo Project start...
echo Reading primer info...
PRIMER=$(echo Merge_DNA_primer\.txt)   # 这里输入你的引物文件路径
echo Done
echo Reading Sequences info...
FASTA=$(echo Merge_DNA_assembled\.fasta)  #这里输入你overlap拼接后的fasta文件
echo Done
echo ===============================================
echo +++++++++++++ Separating Sequences ++++++++++++
echo ===============================================
for line in `cat $PRIMER | sed 's/\t/_/g' | sed 's/\n/ /g'`;do
	 name=$(echo "$line" | awk -v FS='_' -v OFS='_' '{print $1}')
	 a=$(echo "$line" | awk -v FS='_' '{print $2}')
	 b=$(echo "$line" | awk -v FS='_' '{print $3}' | rev | tr ATGC TACG)
	 grep -B 1 ^"$a".*"$b"$ $FASTA | sed '/--/d' > $name.md
	 echo Picking up $name with head=$a tail=$b
	 seqtk seq -Ar $FASTA | grep -B 1 ^"$a".*"$b"$ $FASTA | sed '/--/d' >> $name.md
	 seqtk trimfq -b 8 -e 8 "$name.md" > $name.trim.md
	 fastx_collapser -i "$name.trim.md" -o $name.fas
	 rm *.md
done


echo ========================================
echo +++++++++ Filtering sequences ++++++++++
echo ========================================
mkdir SeqS
mv *.fas ./SeqS
File_list=$(ls ./SeqS/*fas)                            #读取文件夹中.fasta文件列表
Threshold=0.005                                   #过滤频率设定为0.005（100个个体，即1/200）
mkdir ./SeqS/Threshold_filter                          #新建一个文件夹Threshold用于存放过滤后的文件
for FASTA in $File_list; do                       #对每一个.fasta文件进行迭代处理
	echo Doing $(echo $FASTA | sed 's/\.\/SeqS\///g' | sed 's/\.fas//g')                             #开始处理
	touch $FASTA                                  #在Threshold文件夹下创建一个同名空文本
	ready=$(cat $FASTA | sed '/^$/d')             #去除每个.fasta中的空行并存入变量ready
	Total=0                                       #初始化总序列数
	ratio=0                                       #初始化ratio
	for line in $ready;do                         #逐一读取ready中的没一行
		Found=$(echo $line | grep \>)             #检测字符(>)是否存于在该行
		if [[ "$Found" != "" ]]                   #如果存在
		then                                      #则
			name=$line                            #此行为序列名存于name变量
		else                                      #否则
			count=0                               #初始化count
			seq=$line                             #此行为序列存于seq变量
			count=$(echo $name | awk -v FS='-' '{print $2}')  #将该序列频数传入count
			Total=$[$Total + $count]                          #序列总数增加count
			ratio=$(echo "scale=5;$count/$Total" | bc)        #计算该序列频率ratio,浮点计算调用 bc, scale=5表示保留浮点数小数位5个有效数字
			if [ $(echo "$ratio < $Threshold" | bc) -eq 1 ]   #如果该序列ratio小于Threshold,浮点比较调用 bc,返回值为1即True
			then                                              #则
				echo done                                     #结束
				break                                         #跳出for循环
			else                                              #否则
				echo $name >> ./SeqS/Threshold_filter/$(echo $FASTA | sed 's/\.\/SeqS\///g')            #将该序列名存入文本文件
				echo $seq >> ./SeqS/Threshold_filter/$(echo $FASTA | sed 's/\.\/SeqS\///g')             #将该序列存入文本文件
			fi                                                #结束该条件结构
		fi                                        #结束该条件结构
	done                     #结束该循环
done        #结束该循环

#Merge Same Sequences
echo ========================================
echo +++++++++ Identifying Haplotypes +++++++
echo ========================================
cd ./SeqS/Threshold_filter                             #进入Threshold
Files=$(ls *fas)                                #创建Threshold目录下文件列表
PoPgenes=$(echo $Files | sed 's/ /\n/g' | sed 's/\.fas//g' | awk -v FS='-' -v OFS='-' '{$NF="";print $0}' | sed 's/-$//g' | sort | uniq)  #提取文件列表中的种群基因列表
echo The following genes need to be identified..
echo
echo $PoPgenes
mkdir ./Haplotypes                                #创建一个新的文件夹
for gene in $PoPgenes;do                          #逐一提取基因名
	echo
	echo Identifying $gene
	echo 
	mkdir ./Haplotypes/$gene                      #创建对应基因的文件夹
	
	#构建该基因所有单倍型列表
	for FASTA in $Files; do
		Found=$(echo $FASTA | grep $gene)
		if [[ "$Found" != "" ]]
		then
			cat $FASTA >> ./Haplotypes/$gene/All_seq.fas
			cp $FASTA ./Haplotypes/$gene/$FASTA.md
		fi
	done
	fastx_collapser -i ./Haplotypes/$gene/All_seq.fas -o ./Haplotypes/$gene/All.txt
	rm ./Haplotypes/$gene/*fas
	Alls=$(cat ./Haplotypes/$gene/All.txt)
	
	#不同种群该基因相同单倍型识别
	Files_list=$(ls ./Haplotypes/$gene/*md)
	num=0
	for hap in $Alls;do
		Found=$(echo $hap | grep \>)
		if [[ "$Found" != "" ]]
		then
			num=$[$num+1]
			name=$(echo $gene-$num)
			echo $name done ####
		else
			for FASTA in $Files_list;do
				ready=$(cat $FASTA | sed '/^$/d')
				for seq in $ready;do
					F=$(echo $seq | grep \>)
					if [[ "$F" != "" ]]
					then
						count=$(echo $seq | awk -v FS='-' '{print $2}')
						Haplotype_name=$(echo \>$name Frequency:$count)
					else
						if [ $seq = $hap ]
						then
							echo $Haplotype_name >> $(echo $FASTA | sed 's/\.fas\.md/Pop\.fasta/g')
							echo $seq >> $(echo $FASTA | sed 's/\.fas\.md/Pop\.fasta/g')
						fi
					fi
				done
			done
		fi
	done
	rm ./Haplotypes/$gene/*txt
	rm ./Haplotypes/$gene/*md
done
echo ===============================================
echo +++++++++++++ Project is finished! ++++++++++++
echo ===============================================
echo Copyright by Insect Toxicology and Biotechnology Research Group
echo Institution of Zoology, CAS