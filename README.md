# MarViN
Genotype imputation and refinement

##Installation
MarViN relies on HTSlib and Eigen. [Eigen](http://eigen.tuxfamily.org/index.php?title=Main_Page) is a header only library for matrix algebra, included. [HTSlib](http://www.htslib.org/) is a library for efficently parsing vcf files. To install first install HTSlib and then type make to create the two executables: marvin and marvin prep.

##marvin_prep
Set up means and covariances for marvin

Usage:
```./marvin_prep -f input.vcf.gz -O z -o sites.vcf.gz -om mu.dat -ov sig.dat -ow omega.dat```
Expects input.vcf.gz to contain hard genotypes only. Options:

* -f : The input vcf file. Can be valid vcf or bcf compressed or uncompressed. Should contain at least two samples and must have a GT field with no missing entries for any sample/variant combination.
* -o : marvin prep outputs a site only vcf containing the sites kept from the input vcf. Specify the file name with this argument.
* -O : output format for site only vcf. Compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Default is v.
* -om : file name for the vector of allele frequencies.
* -os : file name for the covariance matrix. This is not required by marvin so can be left out.
* -ow : file name for the matrix of inverses.
* -ov : file name for the vector of variances.
* -sigma_reg : use sigmoid function of allele frequency in regularization. Performs better at low allele frequencies.
* -lambda : regulatization parameter, changing has some effect on the results but 0.06 was optimal for most of our tests.
* -lambda2 : regulatization parameter, controls steepness of sigma regularization.
* -pct : regulatization parameter, controls the midpoint of the sigma regularization.
* -r : regions chr:start-end. Section of the genome to operate on.

##marvin
Imputation from genotype likelihoods

Usage:
```./marvin -f input.vcf.gz -O z -o out.vcf.gz```
Expects input_filename.vcf to contain GL or PL field
* I/O
  * -f : the input vcf file. Can be valid vcf or bcf compressed or uncompressed. 
  * -o : output file name containing the imputed genotypes and dosages.
  * -O : output format. Compressed BCF (b), uncompressed BCF (u), compressed VCF (z), uncompressed VCF (v). Default is v.
  * -Ogl : Add GP and GQ field in input with posterior probabilities computed by marvin.
  * -r : regions chr:start-end. Section of the genome to operate on.
  * -b : padding to the left and right.
* Panel input
  * -fm : vector of allele frequencies output by marvin_prep.
  * -fw : matrix of inverses output by marvin_prep.
  * -fv : vector of variances output by marvin_prep.
  * -site : site only vcf output by marvin_prep.
  * -zm : zero missing rows (assumes 0/0)
  * -c : collapse snps|indels|both|all|some|none. Controls how intersection of sample and panel is performed. Similar to bcftools isec.
* Run pararmeters
  * -max_its : Number of ‘outer’ iterations of marvin (re-estimations of the covariance matrix). Only has an effect when not using panel.
  * -inner_its : Number of ‘inner’ iterations of marvin (re-estimations of data with fixed covariance).
  * -maxlr : Maximum allowed likelihood ratio. Likelihood ratios larger than specified threshold will set the smaller likelihood to zero.
  * -bias : marvin works with expected values, when calling hard genotypes if the expected value is less than this parameter it reports hom-ref. If the expected value is between bias and 2-bias it reports het if above 2-bias reports hom-alt.
  * -EMits : When not using a panel marvin does EM on allele frequencies to compute an initial guess. Specifies how many iterations to perform (fewer is generally better, default is 1).
* Regularization parameters
  * -sigma_reg : use sigmoid function of allele frequency in regularization. Performs better at low allele frequencies.
  * -lambda : regulatization parameter, changing has some effect on the results but 0.06 was optimal for most of our tests.
  * -lambda2 : regulatization parameter, controls steepness of sigma regularization.
  * -pct : regulatization parameter, controls the midpoint of the sigma regularization.

#Examples
##Call genotypes from a population given likelihoods
Assuming that `input.vcf.gz` contains genotype likelihoods and is indexed (.csi or .tbi)
```
./marvin -f input.vcf.gz -O z -o out.vcf.gz -r chr20:200000-200000 -b 5000
```
`out.vcf.gz` will contain genotypes imputed under MarViNs LD model, using overlap of 5000 on either side of the specified window.

##Improve calls in sample given a reference panel
```
./marvin_prep -f panel.vcf.gz -O z -o sites.vcf.gz -om mu.dat -ov sig.dat -ow omega.dat -r chr20:195000-205000
```
Preprocessing the panel, which must be done once to precalculate the necessary matrix inverses.
```
./marvin -f input.vcf.gz -O z -o out.vcf.gz -r chr20:200000-200000 -b 5000 -fm mu.dat -fv sig.dat -fw omega.dat -max_its 1 -inner_its 5
```
`out.vcf.gz` has the GT field added or overwritten with the imputed genotypes.


