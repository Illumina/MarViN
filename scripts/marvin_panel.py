#!/usr/bin/env python 

import sys,subprocess,os,argparse,time,multiprocessing,tempfile,shutil,itertools,gzip

def parse_line(line):
    return tuple(line.split())

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='performs reference panel based genotype refinement (imputation) using marvin.')
    parser.add_argument('input', metavar='input', type=str, help='input bcf/vcf')
    parser.add_argument('-output', metavar='output', type=str, help='output bcf',required=True)
    parser.add_argument('-panel', metavar='panel', type=str, help='marvin reference panel (bed file)',required=True)
    parser.add_argument('-marvin', metavar='marvin', type=str, help='marvin binary',required=True)
    parser.add_argument('-bcftools', metavar='bcftools', type=str, help='bcftools binary',required=True)
    parser.add_argument('-tmp', metavar='tmp', type=str, help='tmp directory',default="/tmp/")
    parser.add_argument('-nprocess', metavar='nprocess', type=int, help='number of processes to use (defaults to number of CPUs found)',default=None)
    parser.add_argument('--phred', action='store_true')
    args = parser.parse_args()

    ##binaries
    if not os.path.isfile(args.marvin):
        sys.exit(args.marvin + " binary not found!")

    if not os.path.isfile(args.bcftools):
        sys.exit(args.bcftools + " binary not found!")

    if os.path.isfile(args.output):
        sys.exit(args.output+" already exists! will not overwrite")

    if not os.path.isfile(args.input):
        sys.exit(args.input+" does not exist!")

    if args.nprocess==None:
        args.nprocess=multiprocessing.cpu_count()
    print "Using",args.nprocess,"processes"

    tmp_dir = tempfile.mkdtemp(prefix=args.tmp)
    print "tmp dir is",tmp_dir
    assert os.path.isdir(tmp_dir)
    region_list=[]

#~/kimura/marvin/marvin -f NA12878.fb.1000G.v4.bcf  -fm /illumina/scratch/denovo/rarthur/Marvin/bayes_panel2/data/marvin/mu.1.11995000-12205000.dat -fw /illumina/scratch/denovo/rarthur/Marvin/bayes_panel2/data/marvin/omega.1.11995000-12205000.dat -fv /illumina/scratch/denovo/rarthur/Marvin/bayes_panel2/data/marvin/var.1.11995000-12205000.dat -inner_its 5 -O b -o test.marvin.bcf  -site  /illumina/scratch/denovo/rarthur/Marvin/bayes_panel2/data/marvin/sites.1.11995000-12205000.bcf -r 1:11995000-12205000

    for r  in itertools.imap(parse_line,gzip.open(args.panel)):
        assert(len(r)==3)
        region_list.append(r)


    ##collate chunks with multiprocessing pool
    def process_chunk(r):
        region = "%s:%s-%s"%r
        tmp_out = "%s/%s.bcf"%(tmp_dir,region)
        suffix=".%s.%s-%s.dat"%r
        panel_dir=os.path.dirname(args.panel) + "/chr%s/"%r[0]
        sites_file="%s/sites.%s.%s-%s.bcf"%(panel_dir,r[0],r[1],r[2])
        cmd = " OMP_NUM_THREADS=1 "+args.marvin + " -Ogp -O b -f " + args.input + " -inner_its 5 -r "+region+ " -o " + tmp_out + " -fm %s/mu%s "%(panel_dir,suffix) +  "-fw %s/omega%s "%(panel_dir,suffix) + "-fv %s/var%s "%(panel_dir,suffix) + "-site %s"%sites_file

        if args.phred:            cmd += " -phred"

        try:
            print cmd
            subprocess.check_output(cmd,shell=True)
            subprocess.check_output(args.bcftools+" index "+tmp_out,shell=True)
            return(tmp_out)
        except:
            raise Exception("problem running marvin on %s\n"%region)
        

    time0=time.time()
    print "Running marvin..."
    pool = multiprocessing.Pool(processes=args.nprocess)
    output_files = pool.map(process_chunk, region_list)
    print "Marvin took ",time.time()-time0,"seconds"
    output_files.sort()


    ##concatenate the final chunks into our database
    time0=time.time()
    print "Concatenating..."
    open(tmp_dir+"/concat.txt","wb").write("\n".join(output_files))

    cmd = args.bcftools + " concat -a -D -f" + tmp_dir + "/concat.txt -Ou | bcftools view -O b -o "+args.output
    print cmd
    subprocess.call(cmd,shell=True)
    print "Concatenation took ",time.time()-time0,"seconds"

    shutil.rmtree(tmp_dir)



