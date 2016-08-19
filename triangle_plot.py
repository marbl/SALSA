import argparse, os, sys
import matplotlib.pyplot as plt
import random

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-a","--alignment",help="bed file for alignment")
    #parser.add_argument("-o","--output",help="file to save the plot")
    args = parser.parse_args()
    # ofile = 'coords_'+str(random.randint(0,100))
    # cmd = './triangle_plot -a {alignment} -o {coords}'.format(alignment = args.alignment, coords = ofile)

    # os.system(cmd)

    x = []
    y = []
    curr_contig = ''
    first = True
    with open('coords_32','r') as f:
        lines = f.readlines()
        for line in lines:
            #print line
            attrs = line.split()
            if len(attrs) == 1:
                if first:
                    curr_contig = attrs[0]
                    first = False
                    continue
                plt.scatter(x,y,alpha=0.005)
                plt.xlabel('Midpoints between two mates')
                plt.ylabel('Distance between two mates')
                plt.title('Triangle Plot ' + curr_contig)
                plt.xlim(xmin=-0.1)
                plt.ylim(ymin=-0.1)
                plt.savefig(curr_contig+'_triangle.png')
                curr_contig = attrs[0]
                x = []
                y = []
            else:
	            x.append(float(attrs[0]))
	            y.append(float(attrs[1]))

        


if __name__ == '__main__':
    main()