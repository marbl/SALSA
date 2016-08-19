import os
import argparse
import sys

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--assembly",help="assembled contigs",required=True)
    parser.add_argument("-m","--mapping", help="mapping of read to contigs in bam format",required=True)
    parser.add_argument("-d","--dir",help="output directory for results",default='out',required=True)
    parser.add_argument("-b",'--missassembly',help="add flag to find and break misassemblies from the contigs",default=False)

    args = parser.parse_args()

    if not os.path.exists(args.dir):
        os.makedirs(args.dir)

    # if os.path.exists(args.dir+'/alignment_unique.bed') == False:
    #     print >> sys.stderr, "convertng bam file to bed file"
    #     os.system('bamToBed -i ' + args.mapping + " > " + args.dir+'/alignment_unique.bed')
    #     print >> sys.stderr, 'finished conversion'
    # os.system('samtools faidx '+args.assembly)
    # os.system('cut -f 1,2 '+ args.assembly+'.fai > '+args.dir+'/contig_length_new')

    # print >> sys.stderr, 'finished conversion'

    final_assembly = args.assembly
    final_mapping = args.mapping


    # if args.missassembly:
    #     print >> sys.stderr, 'started finding misassemblies'
    #     final_assembly = args.dir+'/cleaned.fa'
    #     os.system('./break_contigs -a '+args.dir+'/alignment_unique.bed' +' -d ' + args.dir + ' -l ' + args.dir+'/contig_length_new')
    #     os.system('python break_contigs.py -b '+args.dir+'/breakpoints' + ' -a ' + args.assembly + ' -l ' + args.dir+'/contig_length_new' + ' -o ' + final_assembly)
    #     print >> sys.stderr, 'finished detecting misassemblies'

    print >> sys.stderr, 'started generating links between contigs'
    os.system('python RE_sites.py -a '+ final_assembly + ' > '+ args.dir + '/RE_counts')
    os.system('python create_links1.py -m '+args.dir+'/alignment_unique.bed' + ' -c '+ args.dir+'/RE_counts' + ' -l ' + args.dir+'/contig_length_new > '+args.dir+'/new_links')
    os.system('sort -k 3 -g -r '+args.dir+'/new_links > '+args.dir+'/new_links_sorted') 
    print >> sys.stderr, 'done generating links between contigs'

    print >> sys.stderr, 'started scaffolding'

    os.system('python links_to_graph3.py -a '+final_assembly+' -f '+args.dir+'/scaffolds.fasta'+' -l '+args.dir+'/new_links_sorted'+ ' -n '+args.dir+'/contig_length_new > '+args.dir+'/paths')

    print >> sys.stderr, 'done scaffolding, sequences written to file'

if __name__ == '__main__':
    main()
