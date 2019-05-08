#！ /bin/bash

echo "'file ID','Numreds'"
for dir in `ls|grep -v 'pipeline_log.txt'`
do  
  Num=`cat "$dir"/*.id*| cut -f 3 |awk '{sum += $1};END {print sum}'`
  echo "$dir,$Num"

done
