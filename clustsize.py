import re
import sys

clust_fn = sys.argv[1]
clust_file = open(clust_fn, "r")

# "5101-1140-MS515F-926R_clustered.clstr","r")

# >M02696:83:000000000-B82WP:1:1101:16593:1819 1:N:0:82
# 0       251aa, >M02696:83:000000000-B82WP:1:1101:15652:1240... *
for aline in clust_file.readlines():
	striped_line = aline.rstrip()
	if(re.match("^>Cluster (.+)", aline)):
		q = re.search("^>Cluster (.+)", aline)
		clustnum = q.group(1)
	else:
		bits = re.search("^(.*)\t.*>(.*)\.\.\.", striped_line)
		temp = clustnum + "\t" + bits.group(1) + "\t" + bits.group(2) 
		print(temp)
clust_file.close()

