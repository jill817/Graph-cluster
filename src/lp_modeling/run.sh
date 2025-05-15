# for i in {0..0}
# {
# 	# z= "{$i*100+1}","${$i*100+100}"
# 	# a=$(i)*100+1;
# 	bound=`expr $i \* 100`
# 	low=`expr $bound + 1`
# 	up=`expr $bound + 100`
# 	filebound=$low","$up""p
# 	sed -n $filebound all.txt | while read line2
# 	do
# 	{
# 		result=$(nohup ./main exp ${line2} 2>&1 &);
# 	    # result="${line2}     ${result}";
# 	    # echo ${result} >>result.txt
#     }&
# 	done
# 	# sleep 305s
# }

# sed -n "1,100p" List.txt | while read line2
# do
    line2=1000
    line3="/home/hexiang/pnia_ls/work/experiment/${line2}.res"
    # line3="/home/hexiang/scip/ali_work/experiment/${line2}.res"
	(nohup ./nia_ls exp ${line2} >${line3} 2>&1 &)
# done