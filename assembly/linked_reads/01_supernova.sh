

mreads=(636000000*60)/150~270000000

supernova run --fastqs $FastqDir/$sp --id $sp   --maxreads=$mreads --localcores=24 --localmem=185
                                              

 ####### Create pseudohap assembly ########
        supernova mkoutput --style=pseudohap \
        --asmdir=$sp/outs/assembly \
        --outprefix=$sp/outs/pseudohap_merged_output


