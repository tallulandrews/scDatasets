# https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE77288 
# Tung et al. (2017) Batch effects and the effective design of single-cell gene expression studies. Sci Rep. (7):39921
# PMID: 28045081
# GEO: GSE77288

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_molecules-raw-single-per-sample.txt.gz
#wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_reads-raw-bulk-per-sample.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/suppl/GSE77288_reads-raw-single-per-sample.txt.gz

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE77nnn/GSE77288/matrix/GSE77288_series_matrix.txt.gz

gunzip GSE77288_series_matrix.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE77288_series_matrix.txt > GSE77288.txt
