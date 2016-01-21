/*  glcall.c fills in GT field using best guess from GLs.  EM is used to estimate the allele frequency.

    Copyright (C) 2015 Illumina

    Author: Jared O'Connell <joconnell@illumina.com>

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
    THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.  */

#include <stdio.h>
#include <stdlib.h>
#include <htslib/hts.h>
#include <htslib/vcf.h>
#include <htslib/vcfutils.h>
#include <math.h>

bcf_hdr_t *in_hdr,*out_hdr;
int *gt = NULL, ngt = 0,n;
float *gl = NULL;
int ngl=0;
float prior[3],tmp_gl[3];
float *dosage=NULL;
float af_start; //starting value for af
const char *about(void)
{
  return "call GTs from GLs.\n";
}


float estimate_gt() {
  int i,j,count=0,g;
  float den,af=af_start,old_af=100;
  dosage=(float *)realloc(dosage,n*sizeof(float));
  prior[0] = (1-af)*(1-af);
  prior[1] = 2 * af * (1-af);
  prior[2] = af * af;

  for( i=0;i<n;i++) 
    if(bcf_float_is_missing(gl[i*3]) || bcf_float_is_missing(gl[i*3+1]) || bcf_float_is_missing(gl[i*3+2]) )
      for( j=0;j<3;j++) gl[i*3+j]=1.;
  
  //Expectation-Maximisation to estimate AF 
  while(1)  {
    af = 0.;
    for(i=0;i<n;i++) {
      dosage[i]=0.;
      den = 0.;
      for(j=0;j<3;j++) {
	tmp_gl[j] = prior[j] *pow(10.,gl[i*3 + j]);
	den += tmp_gl[j];
      }
      for(j=0;j<3;j++)	tmp_gl[j]/=den;
      for(j=1;j<3;j++)	dosage[i] += j * tmp_gl[j];
      af+=dosage[i];
    }
    af /= (2.*n);
    prior[2] = af * af;
    prior[1] = 2 * af * (1-af);
    prior[0] = (1. - af)*(1. - af);
    //    assert(af>=0 && af<=1);
    //    fprintf(stderr,"count=%d af=%.5f old_af=%.5f %.5f\n",count,af,old_af,fabs(af-old_af));
    if( (count>0 && fabs(af-old_af)<.000001) || count>100)
      break;
    else
      old_af=af;
    count++;
    //    fprintf(stderr,"%d\n",count);
  }

  //dosage -> GT
  for(i=0;i<n;i++) {
    g = roundf(dosage[i]);
    if(g==0) {
      gt[i*2] = bcf_gt_unphased(0);
      gt[i*2+1] = bcf_gt_unphased(0);
    }
    else if(g==1) {
      gt[i*2] = bcf_gt_unphased(0);
      gt[i*2+1] = bcf_gt_unphased(1);      
    }
    else {
      gt[i*2] = bcf_gt_unphased(1);
      gt[i*2+1] = bcf_gt_unphased(1);
    }      
  }
  return(af);
}

int init(int argc, char **argv, bcf_hdr_t *in, bcf_hdr_t *out)
{

  in_hdr  = in;
  out_hdr = out;
  n =  bcf_hdr_nsamples(in_hdr);
  ngt=2*n;
  gt = (int *)malloc(ngt*sizeof(int));
  bcf_hdr_append(out_hdr, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"allele frequency estimated with EM\">");
  bcf_hdr_append(out_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
  bcf_hdr_append(out_hdr, "##FORMAT=<ID=DS,Number=1,Type=Float,Description=\"Estimated allelic dosage\">");
  float  Nh = 2. * n;
  af_start = 1/(log(Nh) + 0.5772156649 + 1 / (2*Nh) - 1 /(12*Nh*Nh));
  return 0;
}


bcf1_t *process(bcf1_t *rec)
{
  float af;
  if(rec->n_allele==2) {
    bcf_get_format_float(in_hdr, rec, "GL", &gl, &ngl);
    assert(ngl==(3*n));
    //    bcf_get_genotypes(in_hdr, rec, &gt, &ngt);
    //    assert(ngt==(2*n));
    af =  estimate_gt();
//    fprintf(stderr,"%d AF=%f\n",rec->pos+1,af);
    bcf_update_format_float(out_hdr,rec,"DS",dosage,n);
    bcf_update_info_float(out_hdr,rec,"AF",&af,1); 
    bcf_update_genotypes(out_hdr, rec, gt, ngt);
  }
  return rec;
}

void destroy(void)
{
  free(gl);
  free(gt);
}


