#! /bin/bash

while getopts  ":t:l:c:m:g:h" o; do

    case "${o}" in

        t) t=${OPTARG}
            ;;
	l) l=${OPTARG}
            ;;
	c) c=${OPTARG}
            ;;
	m) m=${OPTARG}
            ;;
	g) g=${OPTARG}
            ;;
	h) echo "
			This script is used to download and format sequences from BOLD.
			
			List of non-optional arguments:

                        -t	input taxa
			-l	minimum aminoacids length cutoff
			-c	gen code (e.g. 1 for plastid and nuclear, 5 for mitochondrial invertebrate)
			-m	marker (e.g. COI-5P)
                        -h	help page.
	
			List of optional arguments

			-g	geographic location                            

"
               exit
          ;;
       \?) echo "WARNING! -$OPTARG isn't a valid option"
           exit
          ;;
       :) echo "WARNING! missing -$OPTARG value"
          exit
          ;;
       esac
 done

	if [ -z "$t" ] || [ -z "$l" ] || [ -z "$c" ] || [ -z "$m" ]
		then
	
		echo " WARNING! non-optional argument/s is missing "
	
		exit
	
		fi
  
################################################################################################################ downolad sequences

if [ -z "$g" ]; 

	then echo -e "\n downloading BOLD data for taxa $t and marker $m";

	else echo -e "\n downloading BOLD data for taxa $t and marker $m in $g"; 

	fi

if [ -z "$g" ];

	then ./bold-cli -taxon $t -marker $m -output tmp1.fna;

	else ./bold-cli -taxon $t -marker $m -geo $g -output tmp1.fna;

	fi

if [ ! -s tmp1.fna ]; then echo -e "\n your research on BOLD didn't produce any result"; exit; fi 

################################################################################################################ filter taxonomy

echo -e "\n filtering taxonomy"

old_IFS=$IFS; 

IFS=$'\n'; 

for line in $(cat tmp1.fna); 

	do 

	genus=$(echo $line | awk -F "\t" '{print $20}'); 

	if [[ ! -z "$genus" ]]; 

		then 

		echo $line >> tmp2.fna; 

	fi; 

	unset genus

done

################################################################################################################ parse csv

echo -e "\n finding orfs $t"

old_IFS=$IFS;
 
IFS=$'\n'; 

	for line in $(cat tmp2.fna); 

		do 

		seq=$(echo $line | awk -F "\t" '{print $72}' | sed 's/^-*//;s/-$*//'); 

		header=$(echo $line | awk -F "\t" '{print $1}'); 

		echo -e ">$header\n$seq"; 

	done > tmp3.fna

sed -i -e '1,2d' tmp3.fna

################################################################################################################ collapse identical sequences

cd-hit -i tmp3.fna -o tmp4.fna -c 1.00 &>/dev/null;

############################################################################################################### find ORFs

getorf -sequence tmp4.fna -outseq tmp1.faa -table $c &>/dev/null

awk '/^>/ {printf("\n%s\n",$0);next; } { printf("%s",$0);}  END {printf("\n");}' < tmp1.faa > tmp2.faa

################################################################################################################ pick longest CDSs

echo -e "\n chosing longest ORF"

for i in $(grep ">" tmp2.faa | awk -F "_" '{print $1}' | sort -u); 

	do 

	for j in $(grep $i tmp2.faa | awk '{print $1}'); 

		do 

		size=$(grep -A 1 $j tmp2.faa | tail -1 | wc -c); 

		echo $j $size >> tmp.length; 

		done; 

	length=$(sort -rk 2 -n tmp.length | head -1 | awk '{print $2}'); equal_length=$(grep -cw $length tmp.length); if [ $equal_length -ge 2 ]; then continue; fi		### if two orfs of identical length -> exclude 

	longest=$(sort -rk 2 -n tmp.length | head -1 | awk '{print $1}'); 

	grep -A 1 $longest tmp2.faa >> tmp3.faa; 

	awk -F "_" '{print $1}' tmp3.faa > tmp4.faa

	rm tmp.length; 

done

################################################################################################################ filter CDS length

echo -e "\n finding ORFs longer than $c amminoacids"

for j in $(grep ">" tmp4.faa); 

	do 

	a=$(grep -A 1 $j tmp4.faa | tail -1 | wc -c); 

	l1=$(grep $j tmp4.faa | tr -d ">"); 

	if [ $a -gt $l ]; 

		then 

		grep -A 1 $l1 tmp4.faa >> tmp5.faa
		grep -A 1 $l1 tmp4.fna >> tmp5.fna	
	
 	fi;

	done

################################################################################################################ align and retrotranslate

	if [ -s tmp5.faa ]; 

		then

		mafft --auto --adjustdirection tmp5.faa > tmp6.faa #&>/dev/null

		awk -F "_" '{print $1}' tmp6.faa > tmp7.faa

		pal2nal.pl tmp7.faa tmp5.fna -codontable $c -output fasta > def.aln #&>/dev/null;

		else

		echo -e "\n no sequence passed the length cutoff \n"; exit;

	fi

rm tmp[0-9]*.*
