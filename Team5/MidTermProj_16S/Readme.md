# QIIME2: DADA2 denoising for PacBio HiFi 

qiime dada2 denoise-ccs \
  --i-demultiplexed-seqs demux.qza \
  --p-front AGRGTTYGATYMTGGCTCAG \
  --p-adapter RGYTACCTTGTTACGACTT \
  --p-min-len 1000 --p-max-len 1600 \
  --p-n-threads 32 \
  --o-table table.qza \
  --o-representative-sequences rep-seqs.qza \
  --o-denoising-stats denoising-stats.qza \
  --verbose

  #Look at the DADA2 results

qiime metadata tabulate --m-input-file denoising-stats.qza --o-visualization denoising-stats.qzv
qiime feature-table summarize --i-table table.qza \
  --m-sample-metadata-file ../metadata/sample-metadata.tsv --o-visualization table.qzv
qiime feature-table tabulate-seqs --i-data rep-seqs.qza --o-visualization rep-seqs.qzv
