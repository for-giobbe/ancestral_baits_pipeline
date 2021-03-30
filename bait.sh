#! /bin/bash

while getopts  ":g:t:c:h" o; do

    case "${o}" in

        t) t=${OPTARG}
            ;;
	c) c=${OPTARG}
            ;;
	g) g=${OPTARG}
            ;;
	h) echo "
			This script is used to download and format sequences from BOLD.
			
			List of non-optional arguments:
                        -t	input taxa
			-c	minimum aminoacids length cutoff
			-g	gen code
                        -h	help page.
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

	if [ -z "$t" ] 
		then
		echo " WARNING! non-optional argument/s is missing "
		exit
	fi
  
################################################################################################################ downolad sequences

echo -e "\n downloading BOLD data for taxa: $t"

./bold-cli -taxon $t -output tmp1.fna

################################################################################################################ find ORFs

echo -e "\n finding orfs $t"

old_IFS=$IFS;
 
IFS=$'\n'; 

for line in $(cat tmp1.fna); do 

seq=$(echo $line | awk -F "\t" '{print $72}' | sed 's/^-*//;s/-$*//'); 

header=$(echo $line | awk -F "\t" '{print $1}'); 

echo -e ">$header\n$seq"; done > tmp2.fna

sed -i -e '1,2d' tmp2.fna

getorf -sequence tmp2.fna -outseq tmp1.faa -table $g &>/dev/null

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

	if [ $a -gt $c ]; 

		then 

		grep -A 1 $l1 tmp4.faa >> tmp5.faa
		grep -A 1 $l1 tmp2.fna >> tmp3.fna	
	
 	fi;

	done

################################################################################################################ align and retrotranslate

	if [ -s tmp5.faa ]; 

		then

		linsi tmp5.faa > tmp6.faa #&>/dev/null

		awk -F "_" '{print $1}' tmp6.faa > tmp7.faa

		pal2nal.pl tmp7.faa tmp3.fna -codontable $g -output fasta > def.aln #&>/dev/null;

		else

		echo -e "\n no sequence passed the length cutoff \n"; exit;

	fi

rm tmp[0-9]*.*
