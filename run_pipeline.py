import os
import stat
import argparse
import sys
import subprocess
from subprocess import Popen, PIPE

def check(path):
    with open(path,'r') as f:
        for line in f:
            attrs = line.split()
            if float(attrs[4]) >= 1:
                return False
            else:
                return True


ng50 = []
def NG50(lengths,sz = 0):
    genome_size = sz
    if genome_size == 0:
        genome_size = sum(lengths.values())
    contig_lengths = sorted(lengths.values(),reverse=True)
    lensum = 0
    for i in xrange(len(contig_lengths)):
        lensum += contig_lengths[i]
        if lensum >= genome_size/2:
            return contig_lengths[i]

def main():
    bin=os.path.dirname(os.path.abspath(__file__))
    #print bin
    parser = argparse.ArgumentParser(description="SALSA Iterative Pipeline")
    parser.add_argument('-a','--assembly',help='Path to initial assembly',required=True)
    parser.add_argument('-l','--length',help='Length of contigs at start',required=True)
    parser.add_argument('-b','--bed',help='Bed file of alignments sorted by read names',required=True)
    parser.add_argument('-o','--output',help='Output directory to put results',required=False,default='SALSA_output')
    parser.add_argument('-c','--cutoff',help='Minimum contig length to scaffold, default=1000',required=False,default=1000)
    parser.add_argument('-g','--gfa',help='GFA file for assembly',required=False,default='abc')
    #parser.add_argument('-u','--unitigs',help='The tiling of unitigs to contigs in bed format',required=False,default='abc')
    #parser.add_argument('-t','--tenx',help='10x links tab separated file, sorted by last columnls',required=False,default='abc')
    parser.add_argument('-e','--enzyme',help='Restriction Enzyme used for experiment',required=False,default='AAGCTT')
    parser.add_argument('-i','--iter',help='Number of iterations to run, default = 3',required=False,default=3)
    parser.add_argument("-x",'--dup',help='File containing duplicated contig information',required=False,default='abc')
    parser.add_argument("-s",'--exp',help="Expected Genome size of the assembled genome",required=False,default=0)
    parser.add_argument("-m","--clean",help="Set this option to \"yes\" if you want to find misassemblies in input assembly",required=False,default="yes")
    #parser.add_argument("-d","--dist",help="Maximum distance between pairs to consider for misassembly detection",required=False,default=2000000)
    parser.add_argument("-f","--filter",help="Filter bed file for contigs present in the assembly",required=False,default="no")
    parser.add_argument("-p","--prnt",help="Set this option to \"yes\" if you want to output the scaffolds sequence and agp file for each iteration", required=False,default="no")

    args = parser.parse_args()


    #iteration counter
    iter_num = 1
    genome_size = int(args.exp)
    if not os.path.exists(args.output):
        os.mkdir(args.output)


    '''
    Check if misassembly detection needs to be done
    '''
    log = open(args.output+'/commands.log','w',1)


    if not os.path.exists(args.output+'/scaffold_length_iteration_1'):
        cmd = 'cp '+args.length+' '+args.output+'/scaffold_length_iteration_1'

        try:
            p = subprocess.check_output(cmd, shell=True)
            os.chmod(args.output + "/scaffold_length_iteration_1", stat.S_IRUSR | stat.S_IWUSR | stat.S_IRGRP | stat.S_IWGRP | stat.S_IROTH)
        except subprocess.CalledProcessError as err:
            print >> sys.stderr, str(err.output)
            sys.exit(1)


    if not os.path.isfile(args.output+'/alignment_iteration_1.bed'):
        #os.system('cp '+args.bed+' '+args.output+'/alignment_iteration_1.bed')
        if args.filter == "yes":
            os.system("cut -f 1 "+args.length+" | grep -v '>'  > "+args.output+"/contig_names.txt")
            os.system("grep -f "+args.output+"/contig_names.txt -w  "+ args.bed+ " > "+args.output+"/alignment_iteration_1.bed")
        else:
            os.symlink(os.path.abspath(args.bed),args.output+'/alignment_iteration_1.bed')
        #cmd  = 'ln -s '+os.path.abspath(args.bed)+' '+args.output+'/alignment_iteration_1.bed'
        #try:
          #  p = subprocess.check_output(cmd, shell=True)
        #except subprocess.CalledProcessError as err:
         #   print >> sys.stderr, str(err.output)
          #  sys.exit(1)

    os.system('ln -s ' + os.path.abspath(args.assembly) + ' '+args.output+'/assembly.cleaned.fasta')

    if args.clean == 'yes':
       cmd = bin+'/break_contigs_start -a ' + args.output+'/alignment_iteration_1.bed -l ' + args.output+'/scaffold_length_iteration_1 > ' + args.output+'/input_breaks -s 100'
       log.write(cmd)
       try:
           p = subprocess.check_output(cmd,shell=True)
       except subprocess.CalledProcessError as err:
           print >> sys.stderr,str(err.output)

       cmd = 'python ' +bin+'/correct.py  ' + args.assembly + ' '+args.output+'/input_breaks '+args.output+'/alignment_iteration_1.bed '+args.output
       log.write(cmd)
       try:
           p = subprocess.check_output(cmd,shell=True)
       except subprocess.CalledProcessError as err:
           print >> sys.stderr, str(err.output)

       os.system("mv " + args.output+"//alignment_iteration_1.tmp.bed " + args.output+"//alignment_iteration_1.bed")
       os.system("mv " + args.output+"//asm.cleaned.fasta " + args.output+"//assembly.cleaned.fasta")






    #First get RE sites
    if not os.path.isfile(args.output+'/re_counts_iteration_'+str(iter_num)):
        try:
            cmd = 'python '+bin+'/RE_sites.py -a '+args.output + '/assembly.cleaned.fasta -e '+ args.enzyme + ' > '+ args.output+'/re_counts_iteration_'+str(iter_num)
            print cmd
            log.write(cmd+'\n')
            p = subprocess.check_output(cmd,shell=True)

        except subprocess.CalledProcessError as err:
            # if os.path.isfile(args.output+'/re_counts_iteration_'+str(iter_num)):
            #   os.system(args.output+'/RE_sites_iteration_'+str(iter_num))
            print >> sys.stderr, str(err.output)
            sys.exit(1)

    #Now compute normal links with old new_links code
    print >> sys.stderr, "Starting Iteration "+ str(iter_num)
    if not os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
        try:
            cmd = 'python '+bin+'/make_links.py -b '+ args.output+'/alignment_iteration_1.bed' + ' -d '+ args.output +' -i '+str(iter_num) + ' -x ' + args.dup
            print cmd
            log.write(cmd+'\n')
            p = subprocess.check_output(cmd,shell=True)

        except subprocess.CalledProcessError as err:
            print >> sys.stderr, str(err.output)
            if os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
                os.system('rm '+args.output+'/contig_links_iteration_'+str(iter_num))
            sys.exit(1)

    #now use Serge's code to calculate
    if not os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
        try:
            cmd =  'python '+bin+'/fast_scaled_scores.py -d '+args.output+' -i '+str(iter_num)
            log.write(cmd+'\n')
            print cmd
            p = subprocess.check_output(cmd,shell=True)
        except subprocess.CalledProcessError as err:
            if os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
                os.system('rm '+args.output+'/contig_links_scaled_iteration_'+str(iter_num))
            print >> sys.stderr, str(err.output)
            sys.exit(1)

    #Sort the links by column 5
    if not os.path.isfile(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
        try:
            cmd = 'sort -k 5 -gr '+args.output+'/contig_links_scaled_iteration_'+str(iter_num)+ ' > '+ args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)
            log.write(cmd+'\n')
            print cmd
            p = subprocess.check_output(cmd,shell=True)
        except subprocess.CalledProcessError as err:
            if os.path.isfile(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
                os.system('rm '+args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num))
            print >> sys.stderr, str(err.output)
            sys.exit(1)

    if args.gfa != 'abc' and not os.path.isfile(args.output+'/tmp.links'):
        try:
            cmd = bin+'/correct_links -g ' + args.gfa + ' -l ' + args.output+'/contig_links_scaled_sorted_iteration_1 > ' + args.output+'/tmp.links'
            log.write(cmd+'\n')
            p = subprocess.check_output(cmd,shell=True)
        except  subprocess.CalledProcessError as err:
            print >> sys.stderr, str(err.output)
            sys.exit(1)

        os.system('mv ' + args.output+'/tmp.links '+args.output+'/contig_links_scaled_sorted_iteration_1')


    if not os.path.isfile(args.output+'/scaffolds_iteration_1.p'):
        try:
            cmd = 'python '+bin+'/layout_unitigs.py -x '+args.gfa + ' -l '+args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num) +' -c ' +str(args.cutoff)+' -i 1 -d '+args.output
            log.write(cmd+'\n')
            print cmd
            p = subprocess.check_output(cmd,shell=True)

        except subprocess.CalledProcessError as err:
            if os.path.isfile(args.output+'scaffolds_iteration_1.p'):
                os.system('rm '+args.output+'scaffolds_iteration_1.p')
            print >> sys.stderr, str(err.output)
            sys.exit(1)
    if not os.path.isfile(args.output+'/misasm_iteration_'+str(iter_num+1)+'.report'):
        try:
            cmd = bin+'/break_contigs -a ' + args.output+'/alignment_iteration_'+str(iter_num+1)+'.bed -b ' + args.output+'/breakpoints_iteration_'+str(iter_num+1)+'.txt -l '+ args.output+'/scaffold_length_iteration_'+str(iter_num+1) + ' -i '+str(iter_num+1)+' -s 100   > ' + args.output+'/misasm_iteration_'+str(iter_num+1)+'.report'
            p = subprocess.check_output(cmd,shell=True)
            print cmd
            log.write(cmd+'\n')
        except subprocess.CalledProcessError as err:
            print  >> sys.stderr, str(err.output)
            sys.exit(1)

    if not os.path.isfile(args.output+'/misasm_'+str(iter_num+1)+'.DONE'):
        try:
            cmd = 'python '+bin+'/refactor_breaks.py -d ' + args.output + ' -i ' + str(iter_num+1)
            p = subprocess.check_output(cmd,shell=True)
            print cmd
            log.write(cmd+'\n')
            #os.system('mv scaffold_length_iteration_'+str(iter_num+1)+'_tmp scaffold_length_iteration'+str(iter_num+1))
            #os.system('mv re_counts_iteration_'+str(iter_num+1)+'_tmp re_counts_iteration_'+str(iter_num+1))
            #os.system('mv alignment_iteration_'+str(iter_num+1)+'_tmp.bed alignment_iteration_'+str(iter_num+1)+'.bed')
        except subprocess.CalledProcessError as err:
            print  >> sys.stderr, str(err.output)
            sys.exit(1)

    if args.prnt == 'yes':
        cmd = 'python ' + bin+'/get_seq.py -a '+ args.output +'/assembly.cleaned.fasta -f ' + args.output+'/scaffolds_ITERATION_'+str(iter_num)+'.fasta -g ' + args.output+'/scaffolds_ITERATION_'+str(iter_num)+'.agp -p '+ args.output+'/scaffolds_iteration_'+str(iter_num)+'.p'
        log.write(cmd+'\n')
        try:
            p = subprocess.check_output(cmd,shell=True)
        except subprocess.CalledProcessError as err:
            print >> sys.stderr, str(err.output)


    iter_num += 1
    scaffold_length = {}
    with open(args.output+'/scaffold_length_iteration_'+str(iter_num),'r') as f:
        for line in f:
            attrs = line.split()
            scaffold_length[attrs[0]] = int(attrs[1])
        ng50.append(NG50(scaffold_length,genome_size))
    if iter_num - 1 == int(args.iter):
        cmd ='python '+bin+'/get_seq.py -a '+ args.output + '/assembly.cleaned.fasta -f ' + args.output+'/scaffolds_FINAL.fasta -g ' + args.output+'/scaffolds_FINAL.agp -p '+args.output+'/scaffolds_iteration_'+str(args.iter)+'.p'
        log.write(cmd+'\n')
        os.system(cmd)
        sys.exit(0)
    #now do iterative
    while True:
        print >> sys.stderr, "Starting Iteration "+ str(iter_num)
        if not os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
            try:
                cmd = 'python '+bin+'/make_links.py -b '+ args.output+'/alignment_iteration_'+str(iter_num)+'.bed' + ' -d '+ args.output +' -i '+str(iter_num)
                print cmd
                p = subprocess.check_output(cmd,shell=True)
                log.write(cmd+'\n')
                #os.system("rm "+args.output+'/alignment_iteration_'+str(iter_num)+'.bed')

            except subprocess.CalledProcessError as err:
                print >> sys.stderr, str(err.output)
                if os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
                    os.system('rm '+args.output+'/contig_links_iteration_'+str(iter_num))
                sys.exit(1)

        print >> sys.stderr, "Starting Iteration "+ str(iter_num)
        if not os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
            try:
                cmd =  'python '+bin+'/fast_scaled_scores.py -d '+args.output+' -i '+str(iter_num)
                log.write(cmd+'\n')
                p = subprocess.check_output(cmd,shell=True)
            except subprocess.CalledProcessError as err:
                if os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
                    os.system('rm '+args.output+'/contig_links_scaled_iteration_'+str(iter_num))
                print >> sys.stderr, str(err.output)
                sys.exit(1)

        if not os.path.isfile(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
            try:
                cmd = 'sort -k 5 -gr '+args.output+'/contig_links_scaled_iteration_'+str(iter_num)+ ' > '+ args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)
                log.write(cmd+'\n')
                p = subprocess.check_output(cmd,shell=True)
            except subprocess.CalledProcessError as err:
                if os.path.isfile(args.output+'/new_links_scaled_sorted_iteration_'+str(iter_num)):
                    os.system('rm '+args.output+'/new_links_scaled_sorted_iteration_'+str(iter_num))
                print >> sys.stderr, str(err.output)
                sys.exit(1)
        #NOW check if any useful link here
        if check(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
            break

        if not os.path.isfile(args.output+'/scaffolds_iteration_'+str(iter_num)+'.p'):
            try:
                cmd = 'python '+bin+'/layout_unitigs.py -x abc -l '+args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num) +' -c ' +str(args.cutoff)+' -i '+str(iter_num)+' -d '+args.output
                print cmd
                log.write(cmd+'\n')
                p = subprocess.check_output(cmd,shell=True)

            except subprocess.CalledProcessError as err:
                if os.path.isfile(args.output+'scaffolds_iteration_'+str(iter_num)+'.p'):
                    os.system('rm '+args.output+'scaffolds_iteration_'+str(iter_num)+'.p')
                print >> sys.stderr, str(err.output)
                sys.exit(1)


        if not os.path.isfile(args.output+'/misasm_iteration_'+str(iter_num+1)+'.report'):
            try:
                cmd = bin+'/break_contigs -a ' + args.output+'/alignment_iteration_'+str(iter_num+1)+'.bed -b ' + args.output+'/breakpoints_iteration_'+str(iter_num+1)+'.txt -l '+ args.output+'/scaffold_length_iteration_'+str(iter_num+1) + ' -i '+str(iter_num+1)+' -s 100  > ' + args.output+'/misasm_iteration_'+str(iter_num+1)+'.report'
                print cmd
                log.write(cmd+'\n')
                p = subprocess.check_output(cmd,shell=True)
            except subprocess.CalledProcessError as err:
                print  >> sys.stderr, str(err.output)
                sys.exit(1)

        if not os.path.isfile(args.output+'/misasm_'+str(iter_num+1)+'.DONE'):
            try:
                cmd = 'python '+bin+'/refactor_breaks.py -d ' + args.output + ' -i ' + str(iter_num+1) + ' > '+args.output+'/misasm_'+str(iter_num+1)+'.log'
                print cmd
                log.write(cmd+'\n')
                p = subprocess.check_output(cmd,shell=True)
                #os.system('mv scaffold_length_iteration_'+str(iter_num+1)+'_tmp scaffold_length_iteration'+str(iter_num+1))
                #os.system('mv re_counts_iteration_'+str(iter_num+1)+'_tmp re_counts_iteration_'+str(iter_num+1))
                #os.system('mv alignment_iteration_'+str(iter_num+1)+'_tmp.bed alignment_iteration_'+str(iter_num+1)+'.bed')
            except subprocess.CalledProcessError as err:
                print  >> sys.stderr, str(err.output)
                sys.exit(1)


        if args.prnt == 'yes':
            cmd = 'python ' + bin+'/get_seq.py -a '+ args.output +'/assembly.cleaned.fasta -f ' + args.output+'/scaffolds_ITERATION_'+str(iter_num)+'.fasta -g ' + args.output+'/scaffolds_ITERATION_'+str(iter_num)+'.agp -p '+ args.output+'/scaffolds_iteration_'+str(iter_num)+'.p'
            log.write(cmd+'\n')
            try:
                p = subprocess.check_output(cmd,shell=True)
            except subprocess.CalledProcessError as err:
                print >> sys.stderr, str(err.output)

        scaffold_length = {}
        with open(args.output+'/scaffold_length_iteration_'+str(iter_num+1),'r') as f:
            for line in f:
                attrs = line.split()
                scaffold_length[attrs[0]] = int(attrs[1])
            ng50.append(NG50(scaffold_length,genome_size))
            curr_sz = len(ng50)
            if ng50[curr_sz - 1] == ng50[curr_sz - 2]:
                cmd ='python '+bin+'/get_seq.py -a '+ args.output + '/assembly.cleaned.fasta -f ' + args.output+'/scaffolds_FINAL.fasta -g ' + args.output+'/scaffolds_FINAL.agp -p '+args.output+'/scaffolds_iteration_'+str(iter_num-1)+'.p'
                log.write(cmd+'\n')
                os.system(cmd)
                sys.exit(0)


        if iter_num - 1 == int(args.iter):
            cmd ='python '+bin+'/get_seq.py -a '+ args.output + '/assembly.cleaned.fasta -f ' + args.output+'/scaffolds_FINAL.fasta -g ' + args.output+'/scaffolds_FINAL.agp -p '+args.output+'/scaffolds_iteration_'+str(args.iter)+'.p'
            log.write(cmd+'\n')
            os.system(cmd)
            sys.exit(0)

        iter_num += 1

    log.close()

if __name__ == '__main__':
    main()

