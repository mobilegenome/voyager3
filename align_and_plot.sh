#!/bin/bash

for f in $(ls *.fa);
	do
    	mkdir -p ${f/.fa/} &&  cd ${f/.fa/};
        mv -u ../$f .;
        RepeatMasker -species mammal $f > RM.out 2>&1;
        rm_len=$(wc -l $f.out | cut -f 1 -d " ")
        echo $rm_len
        if [ "$rm_len" -gt "1" ];
        	then
        		repeats=$(sed  '1,3d' $f.out |awk '{print $10"#"$11}' | sort | uniq)
        		echo $repeats
        		samtools faidx /home/fritjof/fin_whale_genome_assembly/data/2017-09-28_TAREAN_data/RepeatMasker.lib  $repeats > matching_repeats.fa
				lastdb ${f/.fa/} $f && lastal ${f/.fa/} matching_repeats.fa >  ${f/.fa/}.repeats.maf && last-dotplot -y 1000 -x 1000 --lengths1 --lengths2  --border-pixels=50 -s 14 -f /usr/share/fonts/truetype/msttcorefonts/arial.ttf ${f/.fa/}.repeats.maf ${f/.fa/}.repeats.last.png
			fi
        cd ..;
	done
