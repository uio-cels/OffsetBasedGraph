

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
    id = "chr" + chr + "_" + l[5].replace(".", "v")

    lines_out.append("%s    chr%s  %s  %s\n" % (id, chr, start, stop))

f.close()
f2 = open("grch38_alt_loci.txt", "w")
f2.writelines(lines_out)
f2.close()


