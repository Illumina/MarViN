# MarViN
Rapid Genotype Refinement for Whole-Genome Sequencing Data using Multi-Variate Normal Distribution. This software is not commercially supported.

Copyright (c) 2015, Illumina, Inc. All rights reserved.

See our [pre-print](http://biorxiv.org/content/biorxiv/early/2015/11/12/031484.full.pdf) for a full description of the method.

##License
MarViN source code is provided under the [GPLv3](https://git.illumina.com/rarthur/MarViN/blob/master/LICENSE.txt) license. MarViN includes several third-party packages provided under other open source licenses, please see [COPYRIGHT.txt](https://git.illumina.com/rarthur/MarViN/blob/master/COPYRIGHT.txt) for additional details.

MarViN relies on HTSlib and Eigen. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a header-only library for matrix algebra released under the MPL2 license - see the link (https://www.mozilla.org/en-US/MPL/2.0/) and [COPYRIGHT.txt](https://git.illumina.com/rarthur/MarViN/blob/master/COPYRIGHT.txt). [HTSlib](http://www.htslib.org/) is a library for efficently parsing vcf files released under the MIT/Expat License - see the link (https://opensource.org/licenses/MIT) and [COPYRIGHT.txt](https://git.illumina.com/rarthur/MarViN/blob/master/COPYRIGHT.txt).
Both Eigen and HTSlib are included with MarViN.

##Installation
You can install MarViN via the following commands:
```
git clone https://github.com/Illumina/MarViN
cd MarViN/
make
```
Which creates the two executables: `marvin` and `marvin_prep`.

##marvin_prep
Set up means and covariances for marvin. Use this if you have a reference panel. If you want to genotype a large cohort from likelihoods you do not need to run this.

Usage:

```./marvin_prep -f panel.bcf -r 20 -O b -o sites.20.bcf -b 200000 -ov 5000 ```

Expects panel.bcf to contain hard genotypes. 

Options:

* -f : The input vcf file. Can be valid vcf or bcf, compressed or uncompressed. Should contain at least two samples and must have a GT field with no missing entries for any sample/variant combination.
* -o : marvin\_prep outputs a site only vcf containing the sites kept from the input vcf and the correlations above the threshold max_ratio. Specify the site file name with this argument.
* -O : output format for site only vcf. Compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Default is v.
* -r : regions chr or chr:start-end. Section of the genome to operate on. 
* -sigma_reg : use sigmoid function of allele frequency in regularization. Can perform better at low allele frequencies. Default regularization (without this option) is `Sigma(i,i) += lambda` with this option it is `Sigma(i,i) += lambda / (1.0 + exp( lambda2 * ( Sigma(i,i) - pct ) ) )`
* -lambda : regularization parameter, changing this has some effect on the results but 0.06 was optimal for most of our tests. Default is 0.06.
* -lambda2 : regularization parameter, controls steepness of sigma regularization. Default is 4.
* -pct : regularization parameter, controls the midpoint of the sigma regularization. Default is 0.2.
* -b : block size, compute correlation matrix for variants in this block. 
* -ov : overlap between neighbouring blocks. 
* -max_ratio : only keep elements of correlation matrix s.t. |Cij/max(|Cij|)| > r. Default r = 0.01.

##marvin
Imputation from genotype likelihoods

Usage:

With panel

```./marvin -f input.20.bcf -O b -o out.20.bcf -site sites.20.bcf -r 20 -b 200000 -ov 5000```

From likelihoods

```./marvin -f input.bcf -O b -o out.bcf -b 20000 -ov 5000```

Expects input_filename.vcf.gz to contain GL or PL field
* I/O
  * -f : the input vcf file. Can be valid vcf or bcf, compressed or uncompressed. 
  * -o : output file name containing the imputed genotypes and dosages.
  * -O : output format. Compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Default is v.
  * -Ogl : Add GP and GQ field to output.
  * -v : when present output more information at runtime.
* Window params
  * -r : regions chr:start-end. Section of the genome to operate on.
  * -b : block size, compute correlation matrix for variants in this block. 
  * -ov : overlap between neighbouring blocks. 
* Panel params
  * -site : site only vcf output by marvin_prep.
  * -c : collapse snps|indels|both|all|some|none. Controls how intersection of sample and panel is performed. Similar to bcftools isec. Default none.
* Run pararmeters
  * -max_its : Number of ‘outer’ iterations of marvin (re-estimations of the covariance matrix). Only has an effect when not using panel. Default, without panel 5, with panel only 1 allowed.
  * -inner_its : Number of ‘inner’ iterations of marvin (re-estimations of data with fixed covariance). Default 1 without panel, 5 with panel.
  * -EMits : When not using a panel marvin does EM on allele frequencies to compute an initial guess. Specifies how many iterations to perform (fewer is generally better). Default 1.
  * -maxlr : Maximum allowed likelihood ratio. Likelihood ratios larger than specified threshold will set the smaller likelihood to zero.
  * -bias : MarViN works with expected values, when calling hard genotypes if the expected value is less than this parameter it reports hom-ref. If the expected value is between bias and 2-bias it reports het if above 2-bias reports hom-alt. Default 0.5.
* Regularization parameters
  * -sigma_reg : use sigmoid function of allele frequency in regularization. Can perform better at low allele frequencies. Default regularization (without this option) is `Sigma(i,i) += lambda` with this option it is `Sigma(i,i) += lambda / (1.0 + exp( lambda2 * ( Sigma(i,i) - pct ) ) )`.
  * -lambda : regularization parameter, changing has some effect on the results but 0.06 was optimal for most of our tests. Default is 0.06.
  * -lambda2 : regularization parameter, controls steepness of sigma regularization. Default is 4.
  * -pct : regularization parameter, controls the midpoint of the sigma regularization. Default is 0.2.

#Examples

MarViN should be run on a small window (recommend 200Kb). If there are M variants in a window MarViN scales like M<sup>2</sup> for memory consumption and M<sup>3</sup> for speed. Linkage-Disequilibrium, which is responsible for the correlation patterns used by MarViN, typically decays rapidly with distance. In our experiments we found window sizes between 50Kb and 200Kb to be adequate.

##Call genotypes from a population given likelihoods
Assuming that `input.vcf.gz` contains genotype likelihoods in the specified region and is indexed (.csi or .tbi)
```
./marvin -f input.vcf.gz -O z -o out.vcf.gz -r chr20:1000000-2000000 -b 200000 -ov 5000
```
`out.vcf.gz` will contain genotypes imputed under MarViN's LD model between 1Mb and 2Mb on chromosome chr20, in blocks of 200kb, using overlap of 5kb on either side of the window.

##Improve calls in sample given a reference panel
Preprocessing the panel, which must be done once to precalculate the necessary matrix inverses.
```
./marvin_prep -f panel.vcf.gz -O z -o sites.20.vcf.gz -r 20 -b 100000 -ov 5000
```
If panel.vcf.gz contains reference panel data (multisample vcf/bcf of hard genotypes, with no missing sites) for chromosome 20 then the code above computes correlation matrices in blocks of 100kb using overlap of 5kb between neighbouring blocks.
After indexing sites.20.vcf.gz, to impute any new sample run
```
./marvin -f input.vcf.gz -O z -o out.vcf.gz -site sites.20.vcf.gz -r 20 -b 100000 -ov 5000
```
`out.vcf.gz` has the GT field added or overwritten with the imputed genotypes using correlation information from sites.20.vcf.gz.


