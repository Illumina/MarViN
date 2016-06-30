#!/usr/bin/env python 

import sys,subprocess,os,argparse,time,multiprocessing,tempfile,shutil


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='performs low coverage imputation with marvin and concatenates teh chunks')
    parser.add_argument('input', metavar='input', type=str, help='input bcf/vcf')
    parser.add_argument('-w', metavar='w', type=int, help='window size in bp. (analyses w + 2b but only outputs w)',default=200000)
    parser.add_argument('-b', metavar='b', type=int, help='buffer size in bp.  ',default=100000)
    #parser.add_argument('-alpha', metavar='alpha', type=float, help='regularization paramter',default=None)
    parser.add_argument('-lamda', metavar='lamda', type=float, help='regularization paramter',default=None)
    parser.add_argument('-output', metavar='output', type=str, help='output bcf',required=True)
    parser.add_argument('-r', metavar='r', type=str, help='region in chromosome:start-stop format',required=True)
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
    for r in args.r.split(","):
        chromosome,tmp=r.split(":")
        start,stop=tmp.split("-")
        start=int(start)
        stop=int(stop)
        for i in range(start,stop,args.w):
            a=i
            b = a+args.w
#            print a,b
            if b>stop:
                b = stop
                a = stop-args.w

            region = "%s:%010d-%010d"%(chromosome,a,b)
            print i,region
            region_list.append(region)


    ##collate chunks with multiprocessing pool
    def process_chunk(region):

        tmp_out = "%s/%s.bcf"%(tmp_dir,region)
        cmd = " OMP_NUM_THREADS=1 "+args.marvin + " -Ogp -O b -f " + args.input + " -r "+region+ " -o " + tmp_out + " -b %d"%args.b
        # if args.phred:            cmd += " -phred"
        # if args.alpha!=None:            cmd += " -alpha %f"%args.alpha
        # if args.lamda!=None:            cmd += " -lambda %f"%args.lamda
        cmd += " -lambda .06 -lambda2 .06 -regc "
        print cmd
        try:
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
    cmd = args.bcftools + " concat -a -D " + " ".join(output_files) + " -O b -o "+args.output
    print cmd
    subprocess.call(cmd,shell=True)
    print "Concatenation took ",time.time()-time0,"seconds"

    shutil.rmtree(tmp_dir)


