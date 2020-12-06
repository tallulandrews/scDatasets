# Liver Organoids - published
# Multilineage Communivation regulates human liver bud development from pluripotency - Camp et al. 2017
#GSE81252
#GSE40823

# scRNASeq
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81252/suppl/GSE81252_data.cast.log2.lineage.csv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81252/suppl/GSE81252_data.cast.log2.liverbud.csv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE81nnn/GSE81252/matrix/GSE81252_series_matrix.txt.gz

gunzip GSE81252_data.cast.log2.lineage.csv.gz
gunzip GSE81252_data.cast.log2.liverbud.csv.gz
gunzip GSE81252_series_matrix.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE81252_series_matrix.txt > GSE81252_Ann.txt

# bulk
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE40nnn/GSE40823/suppl/GSE40823_SpagnoliFM_RNASeq.FPKM.txt.gz


# scRNASeq part 2
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96981/suppl/GSE96981_data.mouse.liver.csv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96981/matrix/GSE96981-GPL17021_series_matrix.txt.gz
#gunzip GSE96981-GPL17021_series_matrix.txt.gz
#perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE96981-GPL17021_series_matrix.txt > GSE96981_mouse_Ann.txt


# Human 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96981/suppl/GSE96981_data.lb.late.csv.gz 
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96981/suppl/GSE96981_data.lb.transplant.csv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96981/suppl/GSE96981_data.human.liver.csv.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE96nnn/GSE96981/matrix/GSE96981-GPL16791_series_matrix.txt.gz

gunzip GSE96981_data.lb.late.csv.gz 
gunzip GSE96981_data.lb.transplant.csv.gz
gunzip GSE96981_data.human.liver.csv.gz
gunzip GSE96981-GPL16791_series_matrix.txt.gz
perl ~/Data_Processing_Scripts/parse_series_matrix.pl  GSE96981-GPL16791_series_matrix.txt > GSE96981_human_Ann.txt
