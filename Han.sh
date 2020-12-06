# Han X, Chen H, Huang D, Chen H et al. Mapping human pluripotent stem cell differentiation pathways using high throughput single-cell RNA-sequencing. Genome Biol 2018 Apr 5;19(1):47.
# PMID : 29622030
# GEO: GSE107552

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE107nnn/GSE107552/matrix/GSE107552_series_matrix.txt.gz
gunzip GSE107552_series_matrix.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE107552_series_matrix.txt > GSE107552_anno.txt # gives lots of errors but does down load the data.
