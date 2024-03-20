#!/bin/bash

# Script to launch polyDFE with four models on the bootstrap datasets

if ! test -d .\Results
then
	mkdir Results
fi

for input in $(ls *_polyDFE*.txt)
	do	
		output=$(basename ${input} .txt)
		for ((k=1;k<=4;k++)) do
			echo "#!/bin/bash" > run_${input}_$k.sh
			echo "../polyDFE -d ${input} -v 100 -i ../init_model_BandC.txt ${k} -r ../range_model_BandC.txt 1 -m C > ./Results/${output}_modelC${k}_res.txt" >> run_${input}_$k.sh
			echo "exit 0" >> run_${input}_$k.sh
			chmod u+x run_${input}_$k.sh
			sbatch -A snic2017-7-190 -p core -n 1 -t 20:00 -J polydfe_${input}_$k run_${input}_$k.sh
		done
		echo ${input}
	done
	
exit 0