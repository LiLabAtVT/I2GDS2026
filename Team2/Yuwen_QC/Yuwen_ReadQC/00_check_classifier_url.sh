#!/bin/bash
# Find a working download URL for the SILVA 138 full-length classifier (QIIME2 2024.10, sklearn 1.4.2)
for u in \
  https://data.qiime2.org/2024.10/common/silva-138-99-nb-classifier.qza \
  https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-classifier.qza \
  https://data.qiime2.org/classifiers/sklearn-1.4.2/silva/silva-138-99-nb-weighted-classifier.qza; do
  code=$(curl -s -o /dev/null -L -r 0-0 -w '%{http_code} %{size_download} %{url_effective}' "$u")
  echo "$code  <- $u"
done
