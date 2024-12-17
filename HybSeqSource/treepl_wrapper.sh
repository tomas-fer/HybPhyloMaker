#!/bin/bash
#version: 0.9
#Tongjian Liu, molevoliutongjian@gmail.com
#https://github.com/tongjial/treepl_wrapper

# check the arguments

if [ $# != 3 ]
then
	echo -e "Error: incorrect usage"
	echo -e "usage:\ntreepl_wrapper.sh configuration treefile label"
	exit 1
fi

# check the configure and tree file

if [ -f $1 ] && [ -f $2 ]
then
	echo -e "================================================================"
	date
	echo -e "start running, good luck!"
	echo -e "================================================================"
else
	echo -e "Error: can't find the files"
	echo -e "usage:\ntreepl_wrapper.sh configuration treefile label"
	exit 2
fi

echo -e "treepl_wrapper.sh $1 $2 $3"

# generate the prime configure file

cat $1 |\
awk 'BEGIN{print "treefile = '$2'"}{print}END{print "thorough\nprime"}' \
> configure\_prime\_$3

# run primes

for num in $(seq 100)
do treePL configure\_prime\_$3 |\
sed -n '/^opt/,$p' |\
sed 's/.*\(.\)/\1/' |\
sed ':a;N;$!ba;s/\n/ /g' \
>> prime\_$3
done

# generate the cv configure file

cat configure\_prime\_$3 |\
sed 's/^prime/#&/' |\
awk '{print}END{print "cv\ncvoutfile = cv_'$3'\ncvstart = 0.0001\ncvstop = 10000"}' \
> configure\_cv\_$3

# chose the most frequent cv optimal parameters

sort prime\_$3 | uniq -c | sort -nr | sed q | sed 's/^[ \t]*//' |\

awk '{if($3 != "l") $2=$2" o"; print}' |\
awk '{if($5 != "d") $4=$4" o"; print}' |\
awk '{if($7 != "d") $6=$6" o"; print}' |\

sed 's/ /\n/g' | sed '1d' |\
sed 's/o/#/g' |\
sed '1s/^/opt = &/' |\
sed '2s/l/moredetail/' |\
sed '3s/^/optad = &/' |\
sed '4s/d/moredetailad/' |\
sed '5s/^/optcvad = &/' |\
sed '6s/d/moredetailcvad/' >> configure\_cv\_$3

# perform cross validation

treePL configure\_cv\_$3

# generate the smooth configure file

cat configure\_cv\_$3 |\
sed 's/^cv/#&/' |\
awk '{print}END{print "outfile = treepl_'$3'.tre"}' \
> configure\_smooth\_$3

# find the smallest cv score

cat cv\_$3 |\
awk '{printf "%f\t%s\n",$3,$2}' |\
sort -n | sed q | awk '{print $2}' |\
sed 's/[()]//g;s/^/smoothing = /' \
>> configure\_smooth\_$3

# the last step

treePL configure\_smooth\_$3

# clean temporary files

# rm configure\_cv\_$3 configure\_prime\_$3 configure\_smooth\_$3 cv\_$3 out_dates.tre out_dates.tre.r8s prime\_$3

echo -e "================================================================"
echo -e "running finished, please check the result and report the bugs"
date
echo -e "================================================================"

exit 0