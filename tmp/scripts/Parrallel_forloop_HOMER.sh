#!bin/bash
# Multiple HOMER motif analysis instances run in parrallel

for i in {1..5}
do
( echo $'\n Starting primary UP motif analysis in parrallel' $i >> $WorkDir/HOMER_TIMER.txt
  date >> $WorkDir/HOMER_TIMER.txt
sleep $[ ( $RANDOM % 30 )  + 10 ]s
findMotifsGenome.pl $UpARPeaks hg19 $MotifDir/Up/$i -size 200 -p 30
echo $'\n Finished primary UP motif analysis ' $i >> $WorkDir/HOMER_TIMER.txt
date >> $WorkDir/HOMER_TIMER.txt ) &
done
wait
echo "All primary UP processes complete" >> $WorkDir/HOMER_TIMER.txt
date >> $WorkDir/HOMER_TIMER.txt
echo "-------------------" >> $WorkDir/HOMER_TIMER.txt
# repeat each motif analysis several times, because output can vary
# between runs; save each output to a different folder (already set up in R)
# run each instance of motif analysis in parrallel, 
# because there is a lot of single-core 'down time' during each that 
# should be used more efficiently
# but randomly delay each process to stagger them so
# they don't all hit their multi-threaded processes at the same time
# wait for the entire loop to finish before processing the next loop
# record timestamps for each process
# I think this will speed things up considerably
# I forgot to include the sleep command when running this and
# the server didn't explode and it didn't actually seem to slow it down
# I guess the server can handle 150 simultaneous processes
