#!/bin/bash
set -e

#to be added: 
usage()
{
cat << EOF
MRC is a tool for compressing multiple FASTQ files.

Usage: 
Compression - compresses FASTQ datasets. Output written to '*.MRC' file
./MRC.sh -a m -r test.txt (compress with minicom)
./MRC.sh -a p -r test.txt (compress with PgRC)
./MRC.sh -d file.MRC  
Options:
	-h 		print help message
	-a      compression algortihm, p->PgRC1.2, m->minicom, s->SPRING, f->FaStore
	-r      list of fastq files with same read length, each line contain one file
	-t 		number of threads, default: 12
	-k 		length of k-mer, k <= 10, default: 8
	-e 		threshold percentage, default: 2.0
Decompression - decompresses reads. Output written to 'dec' folder
./minicom -d file.MRC
	-d      compressed file .MRC (for decompression) 
 	-t 		number of threads, default: 8
#See README and more supplementary information at:
EOF
# exit 0
}

compress()
{
	output=${filename%.*}${alg}_MRC
	rm -rf $output ${filename%.*}${alg}.MRC
	mkdir $output
	echo ${alg} > ${output}/info

	echo "./MRC -r ${filename} -t ${num_thr} -k ${k} -e ${threshold_per} -o ${output}"
	./MRC -r ${filename} -t ${num_thr} -k ${k} -e ${threshold_per} -o ${output}

	listVar=( )
	while read a b
	do
	        listVar+=($a)
	done < "${filename}"
	while IFS= read -r line
	do
	        fp=""
	        for a in $line
	        do
	                fp="${fp}${a}_"
	        done
  
			rm -rf ${output}/buff.fastq
        	for a in $line
        	do
            	cat ${listVar[a]} >> ${output}/buff.fastq
        	done
        	mv ${output}/buff.fastq ${output}/${fp}.fastq            	

            if [[ $alg = "p" ]]; then
            	./PgRC -o -i ${output}/${fp}.fastq ${output}/${fp}.pgrc
            elif [[ $alg = "m" ]]; then
            	cd minicom 
            	./minicom -r ../${output}/${fp}.fastq -p
            	cd ../
            elif [[ $alg = "s" ]]; then
	       		./spring -c -i ${output}/${fp}.fastq -o ${output}/${fp}.spring	       		       		   	
	       	elif [[ $alg = "f" ]]; then
	       		echo ${fp}.fastq
	       		_cwd="$PWD"
		        cd FaStore
		        sh ./fastore_compress.sh --lossless --in ../${output}/${fp}.fastq --out ../${output}/${fp} --threads 8
				cd ../ 
		    fi
			rm -rf ${output}/${fp}.fastq
	done < "${output}/cluster"
	tar -cf ${filename%.*}${alg}.MRC ${output}
	#rm ${filename%.*}${alg}.MRC
	echo "compressed file write to: ${filename%.*}${alg}.MRC" 
	rm -rf ${output}
}

check_exist()
{
   if [[ $alg = "p" ]]; then
    	if [ ! -f "PgRC" ]; then
			echo "Required tool PgRC does not exist"
			exit
		fi
   elif [[ $alg = "m" ]]; then
    	if [ ! -d "minicom" ]; then
			echo "Required tool minicom does not exist"
			exit
		fi
    elif [[ $alg = "s" ]]; then
    	if [ ! -f "spring" ]; then
			echo "Required tool spring does not exist"
			exit
		fi   		   	
   	elif [[ $alg = "f" ]]; then
   		if [ ! -d "FaStore" ]; then
			echo "Required tool FaStore does not exist"
			exit
		fi
	fi
}

decompress()
{	
	echo "input file: ${filename}"
	dir=${filename%.MRC*}_MRC
	dir="${dir##*/}"
	output=${filename%.MRC*}_decompress 
	echo "dir is ${dir}"
	rm -rf ${dir}
	echo "decompress, file: ${filename}"
	#echo "overwrite -xvkf ${filename}"
	tar --overwrite -xvkf  ${filename} -C ./

    if [ ! -d "${dir}" ]; then
		echo "Error, unzipped file not in current directory"
		exit
	fi

	rm -rf $output
	mkdir $output
	
	file_list=( )
	length_list=( )

	alg=$(head -n 1 ${dir}/info)
	check_exist
	echo "alg:${alg}"

	{
	read
    while read a b
	do
	        file_list+=($a)
	        length_list+=($b)
	done
	} < "${dir}/info"

	while IFS= read -r line
	do
		fp=""
		IFS=', ' read -r -a array <<< $line
		for element in "${array[@]}"
		do
    		fp="${fp}${element}_"
		done
        start=1
        end=0
        if [[ $alg = "p" ]]; then
        	./PgRC -d ${dir}/${fp}.pgrc
        	for element in "${array[@]}"
			do
				end=$(($end + ${length_list[$element]}))
				sed -n -e "${start},${end} p" -e "${end} q" ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				#sed -n \'\' ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				start=$(($end+1))
			done
			rm  ${dir}/${fp}.pgrc_out
        elif [[ $alg = "m" ]]; then
        	cd minicom 
        	./minicom -d ../${dir}/${fp}_comp_order.minicom
        	for element in "${array[@]}"
			do
				end=$(($end + ${length_list[$element]}))
				sed -n -e "${start},${end} p" -e "${end} q" ${fp}_comp_order_dec.reads > ../${output}/${file_list[$element]}
				#sed -n \'\' ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				start=$(($end+1))
			done
			rm ${fp}_comp_order_dec.reads
        	cd ../
        elif [[ $alg = "s" ]]; then
        	./spring -d -i ${dir}/${fp}.spring -o ${output}/${fp}.fastq
        	for element in "${array[@]}"
			do
				line_num=$((4 * length_list[$element]))
				end=$(($end + ${line_num}))
				sed -n -e "${start},${end} p" -e "${end} q" ${output}/${fp}.fastq > ${output}/${file_list[$element]}
				#sed -n \'\' ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				start=$(($end+1))
			done
			rm ${output}/${fp}.fastq 
        elif [[ $alg = "f" ]]; then
        	cd FaStore
        	sh fastore_decompress.sh --in ../${dir}/${fp} --out ../${output}/${fp}.fastq
        	for element in "${array[@]}"
			do
				end=$(($end + ${length_list[$element]}))
				sed -n -e "${start},${end} p" -e "${end} q" ../${output}/${fp}.fastq > ../${output}/${file_list[$element]}
				#sed -n \'\' ${dir}/${fp}.pgrc_out > ${output}/${file_list[$element]}
				start=$(($end+1))
			done
			cd ../
			rm ${output}/${fp}.fastq 
		fi					
	done < "${dir}/cluster"
	#rm -rf ${filename}

	echo "finished, decompressed file write to ${output}"

}
#Initialize variables to default values.
num_thr=12
threshold_per=2
m_dict=0
k=8
alg="PgRc"
filename=""
#Check the number of arguments. If none are passed, print help and exit.
argnum=$#
if [[ $argnum -eq 0 || $1 == "-h" ]]; then
 usage
 exit 1
fi

mode="c"
if [[ $1 == "-d" ]]; then
	mode="d"
fi

while getopts ":a:r:t:k:e" opt; do
	case "$opt" in
		a) alg=$OPTARG;;
		r) filename=$OPTARG;;
		t) num_thr=$OPTARG;;
		k) k=$OPTARG;; #length of k
		e) threshold_per=$OPTARG;; #different threshold
		#\?) usage; echo -e "\033[31m Error parameters. \033[0m"; exit 0;;
		#*) usage; echo -e "\033[31m Error parameters. \033[0m"; exit 0;;
	esac
done

if [[ $mode == "c" ]]; then
	check_exist
	compress
elif [[ $mode == "d" ]]; then
	decompress
fi