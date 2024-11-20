#!/bin/bash


DIRECTORY='/home/kshone/workspace'
MAX_JOBS=60  


for file in $DIRECTORY/*.txt
do
  
  filename=$(basename "$file" .txt)

  
  mkdir -p "$DIRECTORY/$filename"

  
  LOGFILE="$DIRECTORY/$filename/run-prefetch.log"

  
  while read -r line
  do
    
    first_column=$(echo $line | awk '{print $1}')

    
     nohup prefetch "$first_column" --output-directory "$DIRECTORY/$filename" > "$DIRECTORY/$filename/run-prefetch.log" 2>&1 &

    
    taskNum=`ps -aux | grep prefetch | grep -v grep | wc -l`
    while [ "$taskNum" -gt "$MAX_JOBS" ]; do
      echo "The number of tasks remaining: $taskNum"
      sleep 30
      echo `date`
      taskNum=`ps -aux | grep prefetch | grep -v grep | wc -l`
    done

  done < "$file"
done

echo "All tasks have been initiated and are running in the background."
