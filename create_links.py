import operator
import numpy as np
import math

contig_lengths = {}
with open('contig_lengths_new','r') as cfile:
	lines = cfile.readlines()
	for line in lines:
		attrs = line.split()
		contig_lengths[attrs[0]] = float(attrs[1])


sorted_contigs = sorted(contig_lengths.items(), key = operator.itemgetter(1), reverse = True)

sorted_contigs1 = sorted_contigs[0:70]
sorted_contigs = [i[0] for i in sorted_contigs1]

coordinates_map = {}
first_in_pair = {}
second_in_pair = {}

with open('alignment_unique_clean.bed','r') as f:
	lines = f.readlines()
	for line in lines:
		attrs = line.split()
		if attrs[0] in sorted_contigs:
			if attrs[0] not in coordinates_map:
				coordinates_map[attrs[0]] = {}
			read = attrs[3][:-2]
			pos = (long(attrs[1]) + long(attrs[2]))/2.0
			if read in coordinates_map[attrs[0]]:
				coordinates_map[attrs[0]][read].append(pos)
			else:
				coordinates_map[attrs[0]][read] = [pos]
		pos = (long(attrs[1]) + long(attrs[2]))/2.0
		read = attrs[3]
		if read[-1] == '1':
			read = read[:-2]
			rec = (attrs[0],pos)
			first_in_pair[read] = rec
		else:
			read = read[:-2]
			rec = (attrs[0],pos)
			second_in_pair[read] = rec

dist = []
mates = []

lengths = [1000,10000,20000,30000,50000,70000,100000,200000,30000,500000,600000,800000,900000,1000000,1250000,1500000,1750000,2000000,2500000]

count_map = {}
for contig in sorted_contigs:
	alignments = coordinates_map[contig]
	midpoint = contig_lengths[contig]/2.0
	for each in range(100000,6000000,100000):
		if each not in count_map:
			count_map[each] = []
		count = 0
		left = midpoint - each
		right = midpoint + each
		for read in alignments:
			coords = alignments[read]
			if len(coords) == 2:
				start = min(coords[0],coords[1])
				end = max(coords[0],coords[1])
				if left <= start <= midpoint and midpoint <= end <= right:
					count += 1
		count_map[each].append(count)
	

x = []
y = []
for each in count_map:
	average = reduce(lambda x,y : x + y, count_map[each])/len(count_map[each])*1.0
	x.append(math.log(each))
	y.append(average)

x = np.array(x)
y = np.array(y)

z = np.polyfit(x,y,1)

m = z[0]
c = z[1]
		

contig_links = {}
exp_links = {}

for read in first_in_pair:
	if read in second_in_pair:
		rec1 = first_in_pair[read]
		rec2 = second_in_pair[read]
		if rec1[0] != rec2[0]:
			len1 = contig_lengths[rec1[0]]
			len2 = contig_lengths[rec2[0]]
			pos1 = rec1[1]
			pos2 = rec2[1]
			key = ''
			key1 = {}
			tot_len = len1/2 + len2/2
			exp_num = m*math.log(tot_len) + c
			if pos1 <= len1/2 and pos2 <= len2/2:
				key = rec1[0] + ":B$" + rec2[0] + ':B'
				exp_links[key] = exp_num
			if pos1 <= len1/2 and pos2 > len2/2:
				key = rec1[0] + ":B$" + rec2[0] + ':E'
				exp_links[key] = exp_num
			if pos1 > len1/2 and pos2 <= len2/2:
				key = rec1[0] + ":E$" + rec2[0] + ':B'
				exp_links[key] = exp_num
			if pos1 > len1/2 and pos2 > len2/2:
				key = rec1[0] + ":E$" + rec2[0] + ':E'
				exp_links[key] = exp_num
			if key not in contig_links:
				contig_links[key] = 1
			contig_links[key] += 1

for key in contig_links:
	exp_num = exp_links[key]
	score = 1.0/abs(exp_num - contig_links[key])**0.5	
	edge = key.split('$')
	print edge[0],edge[1],score, exp_num, contig_links[key]
