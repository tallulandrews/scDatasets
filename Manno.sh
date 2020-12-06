# Mouse (adult & embryo), Human embry, & human ipsc-derived dopaminergic neurons
# La Manno G, Gyllborg D, Codeluppi S, Nishimura K et al. Molecular Diversity of Midbrain Development in Mouse, Human, and Stem Cells. Cell 2016 Oct 6;167(2):566-580.e19. PMID: 27716510
# GSE76381

# THIS is good since they did a simplified version of CellTypeProfiles in their paper.

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/suppl/GSE76381_ESMoleculeCounts.cef.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/suppl/GSE76381_EmbryoMoleculeCounts.cef.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/suppl/GSE76381_MouseAdultDAMoleculeCounts.cef.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/suppl/GSE76381_MouseEmbryoMoleculeCounts.cef.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/suppl/GSE76381_iPSMoleculeCounts.cef.txt.gz

wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/matrix/GSE76381-GPL11154_series_matrix.txt.gz
wget ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/matrix/GSE76381-GPL13112_series_matrix.txt.gz

gunzip GSE76381-GPL13112_series_matrix.txt.gz
gunzip GSE76381-GPL11154_series_matrix.txt.gz
gunzip GSE76381_ESMoleculeCounts.cef.txt.gz
gunzip GSE76381_EmbryoMoleculeCounts.cef.txt.gz
gunzip GSE76381_MouseAdultDAMoleculeCounts.cef.txt.gz
gunzip GSE76381_MouseEmbryoMoleculeCounts.cef.txt.gz
gunzip GSE76381_iPSMoleculeCounts.cef.txt.gz

perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE76381-GPL11154_series_matrix.txt > Ann_Hsap.txt
perl ~/Data_Processing_Scripts/parse_series_matrix.pl GSE76381-GPL13112_series_matrix.txt > Ann_Mmus.txt

rm GSE76381-GPL11154_series_matrix.txt
rm GSE76381-GPL13112_series_matrix.txt
