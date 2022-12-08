### chlorobiontas (good genomes)

egrep "bacteria\s+|size" good_genomes_multicellularidad__allDivisions.tsv | awk '{ print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6g"}' | awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobiaceae"}'| cut -f4 | grep -v Lineage | sort
11 => Chlorobiaceae: 11 especies


egrep "bacteria\s+|size" good_genomes_multicellularidad__allDivisions.tsv | awk '{ print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6g"}' | awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobiaceae"}' |  awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobaculum"}' | awk '{ print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$2"\t"$3}' | cut -f1,2,3,4,5
count	mean	var	size_Mbp	gc_percent
2255	843.443458980044	424579.723398855	2.15495	56.5
2452	952.910685154976	439754.103402713	2.79728	56.4
2043	973.828193832599	417599.146276217	2.28925	55.8

egrep "bacteria\s+|size" good_genomes_multicellularidad__allDivisions.tsv | awk '{ print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6g"}' | awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobiaceae"}' |  awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Prosthecochloris"}'  | awk '{ print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$2"\t"$3}' | cut -f1,2,3,4,5
count	mean	var	size_Mbp	gc_percent
2269	953.446452181578	391688.681985864	2.46848	56
2327	969.510958315428	589952.885412968	2.5797	50.1104
2200	999.074545454545	1758164.88848319	2.39985	52.1

emuro@arcturus:~/tmp/goingOn$ egrep "bacteria\s+|size" good_genomes_multicellularidad__allDivisions.tsv | awk '{ print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6g"}' | awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobiaceae"}' |  awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobium"}'  | awk '{ print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$2"\t"$3}' | cut -f1,2,3,4,5
count	mean	var	size_Mbp	gc_percent
2707	953.294421869228	744401.343808489	3.01824	48.1
2434	977.984798685292	745549.556283416	2.76318	51.3
2650	993.618867924528	772769.015876408	3.1339	48.4
2083	1011.17090734518	874560.69027713	2.36484	57.3

egrep "bacteria\s+|size" good_genomes_multicellularidad__allDivisions.tsv | awk '{ print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6g"}' | awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chlorobiaceae"}' |  awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Chloroherpeton"}'  | awk '{ print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$2"\t"$3}' | cut -f1,2,3,4,5
count	mean	var	size_Mbp	gc_percent
2710	1052.11549815498	612754.720136241	3.29346	45








### stramenopiles (good genomes)

egrep "protists\s+|size" good_genomes_multicellularidad__allDivisions.tsv | awk '{ print $1"\t"$2"\t"$3"\t"$9"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8}' | awk 'NR == 1; NR > 1 {print $0 | "sort -k6,6g"}' | awk 'NR == 1; NR > 1 {print $0 | "grep --color=always Stramenopiles"}' | awk '{ print $5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$4"\t"$2"\t"$3}' | cut -f1,2,3,4,5

count	mean	var	size_Mbp	gc_percent
12178	1616.46140581376	1643675.27078597	27.4507	48.8432
11672	1749.45724811515	2232545.87119071	32.4374	46.905
10496	1922.09317835366	1884052.96892535	27.5893	54.5409
16269	6836.10590693958	34733793.3844678	195.811	53.4936







All genomes
###########

Run on:  all_genomes_multicellularidad__allDivisions.tsv 



