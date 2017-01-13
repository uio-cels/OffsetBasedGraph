
sizes = {}
f = open("grch38.chrom.sizes")
for line in f.readlines():
    l = line.split()
    sizes[l[0]] = int(l[1])

def create_alt_loci_file():
    f = open("grch38_alternative_loci.txt")
    lines_out = []
    for line in f.readlines():
        if line.startswith("#"):
            continue

        l = line.split()
        if l[4] != "alt-scaffold":
            continue

        chr = l[1]
        start = l[2]
        stop = l[3]
        id = "chr" + chr + "_" + l[5].replace(".", "v") + "_alt"
        size = sizes[id]

        lines_out.append("%s    chr%s  %s  %s %d\n" % (id, chr, start, stop, size))

    f.close()
    f2 = open("grch38_alt_loci.txt", "w")
    f2.writelines(lines_out)
    f2.close()


def create_alt_loci_file_from_db():
    f = open("hgTables.txt")
    lines_out = []
    for line in f.readlines():
        if line.startswith("#"):
            continue

        l = line.split()
        if "alt" in l[1]:
            continue

        chr = l[1]
        start = int(l[2]) + 1
        stop = l[3]
        id = l[4]
        size = sizes[id]

        lines_out.append("%s    %s  %s  %s %d\n" % (id, chr, start, stop, size))

    f.close()
    f2 = open("grch38_alt_loci_from_db.txt", "w")
    f2.writelines(lines_out)
    f2.close()

def get_genes(alt_loci, chr, start, end):
    # writes all genes on the alt loci and between start and end on the main chromosome to file
    fn = "genes_" + alt_loci + ".txt"
    out = []
    f = open("genes.txt")
    for line in f.readlines():
        d = line.split()
        if d[1] == chr and int(d[3]) >= start and int(d[4]) <= end:
            out.append(line)
        elif d[1] == alt_loci:
            out.append(line)

    o = open(fn, "w")
    o.writelines(out)
    o.close()
    f.close()

get_genes("chr13_KI270838v1_alt", "chr13", 112182875, 112505943)




# create_alt_loci_file_from_db()