#!/bin/bash

alreadyWrittenLine=`cat PolytopeSVRScript.R | grep "Rscript" | wc -l`

tmpFileName=temporaryFile.tmp

if (( alreadyWrittenLine == 1 )); then
	echo "Warning: removing and rewriting header in PolytopeSVRScript.R"
	sed '1d' PolytopeSVRScript.R > $tmpFileName
	mv $tmpFileName PolytopeSVRScript.R
	chmod a+x PolytopeSVRScript.R
fi



if [ $(hostname) == "khachiyan.rutcor.rutgers.edu" ]; then
	firstLine=\#\!/home/ggazzola/R-3.2.1/bin/Rscript
else
	firstLine=\#\!/usr/bin/env\ Rscript
fi


for fileToChange in PolytopeSVRScript.R; do
	
	echo $firstLine >> $tmpFileName
	cat $fileToChange >> $tmpFileName
	cat $tmpFileName > $fileToChange
	rm -f $tmpFileName
done


./PolytopeSVRScript.R
