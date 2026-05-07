# Ajmix

## Figures
Pre-made figures are available in the directory `figures/`

## Requirements
`bcftools`: if you choose to do the optional setup step below
`Python 3` with the following libraries:
  - `pandas`
  - `numpy`
  - `matplotlib`
  - `scipy`

## Optional Setup [if you wish to reproduce from scratch - not necessary]
Note that the filtered gnomAD and subject files are both provided in the `tsv/` directory (pre-computed with same methods), so this setup is unnecessary. You may skip this section.

To begin, you will need to download the benchmark files and the gnomAD allele frequencies. The recombination map is already included here (as it's small).

### gnomAD Allele Frequencies
To download the gnomAD allele frequencies, run the following:
```sh
for i in {18..22}; do
	wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/genomes/gnomad.genomes.v4.1.1.sites.chr$i.vcf.bgz
	wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/genomes/gnomad.genomes.v4.1.1.sites.chr$i.vcf.bgz.tbi
done
```
Note that this will require >56 GB storage. If you have just 15 GB available (minimum required for this project), you may download chr18-22 individually, run the extraction (as seen below), and then delete the file.

If you would like to download chr18 individually, for example, run 
```
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/genomes/gnomad.genomes.v4.1.1.sites.chr18.vcf.bgz
wget -c https://storage.googleapis.com/gcp-public-data--gnomad/release/4.1.1/vcf/genomes/gnomad.genomes.v4.1.1.sites.chr18.vcf.bgz.tbi
```

After you have downloaded the gnomAD files, you may proceed with extracting the necessary information. Here is an example for chr18 (replace if analyzing another chromosome):
```sh
bcftools query -r chr18 -f '%POS\t%REF\t%ALT\t%AF_nfe\t%AN_nfe\t%AF_mid\t%AN_mid\n' gnomad.genomes.v4.1.1.sites.chr18.vcf.bgz | \
awk '$4 != $6' > tsv/gnomad_chr18_informative.tsv
awk -v OFS='\t' '{
    # Drop the row if NFE/MID frequency missing
    if ($4 == "." || $6 == ".") next;

    diff = $4 - $6;
    abs = (diff < 0 ? -diff : diff);
    if (abs > 0.1) print $1, $2, $3, $4, $5, $6, $7
}' tsv/gnomad_chr18_informative.tsv > tsv/gnomad_chr18_filtered.tsv```

### Subject HG002 Data
To download the AJ subject (HG002) data, simply run the following:
```sh

wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
wget -c https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/release/AshkenazimTrio/HG002_NA24385_son/latest/GRCh38/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi

mv HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz benchmark.vcf.gz  # just for ease of use
mv HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz.tbi benchmark.vcf.gz.tbi  # "

for i in {18..22}; do
	bcftools query -r chr$i -f '%POS\t%REF\t%ALT\t[%GT]\n' benchmark.vcf.gz > tsv/hg002_chr${i}_genotypes.tsv
done
```

You are now done and ready to proceed.

## To Run
(note that some of these, specifically `sensitivity_combined.py` and `batch.py`, may take a few minutes to run)
### Reproduce table from paper - run each chromosome on different window sizes
`python3 batch.py`

### Run Baum-Welch on individual chromosome, outputting the admixture time in generations
`python3 windowed.py [chr]` 

### Run Baum-Welch on individual chromosome with specific window size (enter as integer)
`python3 windowed.py [chr] local [window size]`

### Run Baum-Welch on the combined segment chr18-22
`python3 windowed_combined.py`

### Run sensitivity analysis on individual chromosome, producing a plot of admixture time by window size
`python3 sensitivity.py [chr]`

### Run sensitivity analysis on combined segment chr18-22
`python3 sensitivity_combined.py`


## Explanation of each file
`windowed.py, window_combined.py` - run Baum-Welch to determine admixture time (g) for individual, combined chromosome segments, respectively.

`sensitivity.py, sensitivity_combined.py` - Run a sensitivity analysis by iterating over window sizes 10kb-100kb in 1kb increments, plotting results fit to an exponential

`batch.py` - runs Baum-Welch on all chromosomes individually and for window sizes 50kb, range from 45-70kb; internally runs windowed.py; outputs table of results

`filter.py` - helper file; reads the filtered gnomAD .tsv file for a chromosome and returns a df after further filtering and adding helper columns

`sex-averaged_noncarrier.rmap.txt` - a recombination rate map across the entire genome

## Expected Output
Running `python3 windowed_combined.py` should give you the following output (along with intermediate calculations/logs):
```
RESULT: 33.23 generations
Estimated Admixture Date: ~1096 AD
```
Running `python3 batch.py` should give you the following table:
```
Chr    | 50kb (Anch)  | Range (45-70kb)      | Max Delta  | Year (50k)
--------------------------------------------------------------------------------
18     | 27.25        | 24.34 - 30.85        | 6.51       | ~1263 AD
19     | 41.28        | 30.86 - 41.28        | 10.42      | ~871 AD
20     | 33.98        | 23.28 - 37.82        | 14.54      | ~1075 AD
21     | 34.74        | 28.32 - 35.74        | 7.42       | ~1054 AD
22     | 30.69        | 24.82 - 33.48        | 8.66       | ~1167 AD
```

*Note: Admixture dates are assuming a 28-year generation*
