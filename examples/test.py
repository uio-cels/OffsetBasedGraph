
sizes = {}
f = open("grch38.chrom.sizes")
for line in f.readlines():
    l = line.split()
    sizes[l[0]] = int(l[1])



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


