#! /bin/bash/

chunkLen=700000

for i in `seq 1 4667`  ; do
	thisScaf=`head -${i} scaf.lengths | tail -1 | cut -f 1`
	thisLen=`head -${i} scaf.lengths | tail -1 | cut -f 2`
#	echo ${thisScaf} is ${thisLen}bp long

	if [ ${thisLen} -gt ${chunkLen} ] ; then
		theseBreaks=( `seq 1 ${chunkLen} ${thisLen}` )
		nBreaks=`expr ${#theseBreaks[@]} - 1`
		for j in $(seq 0 $(expr ${nBreaks} - 1) ) ; do
			echo "-L ${thisScaf}:${theseBreaks[${j}]}-`expr ${theseBreaks[${j}+1]} - 1`"
		done
		echo "-L ${thisScaf}:${theseBreaks[@]:(-1)}-${thisLen}"
	fi
	if [ ${thisLen} -lt `expr ${chunkLen} / 8` ] ; then
        echo "-L ${thisScaf} 8th"
    elif [ ${thisLen} -lt `expr ${chunkLen} / 2` ] ; then
        echo "-L ${thisScaf} half"
	elif [ ${thisLen} -lt ${chunkLen} ] ; then
        echo "-L ${thisScaf}"
    fi
done > indices.list

sed -i s/\-$//g indices.list

sed -i '/half/{N;s/\n//;}' indices.list
sed -i s/half//g indices.list

for i in `seq 1 3` ; do
	sed -i '/8th/{N;s/\n//;}' indices.list
done
sed -i s/8th//g indices.list

#rm indices.long.list
#touch indices.long.list
#
#for i in `cat indices.list`
#do for j in `seq 1 63`
#  do ind=`head -${j} samples.list | tail -1`
#  echo ${ind} ${i} >> indices.long.list
# done
#done
