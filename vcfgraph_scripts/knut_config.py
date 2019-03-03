chrom_sizes = {
    1:	249250621,
    2:	243199373,
    3:	198022430,
    4:	191154276,
    5:	180915260,
    6:  4970557,
    7:	159138663,
    8:	146364022,
    9:	141213431,
    10:	135534747,
    11:	135006516,
    12:	133851895,
    13:	115169878,
    14:	107349540,
    15:	102531392,
    16:	90354753,
    17:	81195210,
    18:	78077248,
    20:	63025520,
    19:	59128983,
    22:	51304566,
    21:	48129895}


data_path = "/home/knut/Documents/phd/two_step_graph_mapper/benchmarking/mhc_graph_data/"
out_path = "./"
vcf_path = data_path + "1000genomes_variants%s.vcf"
all_vcf_path = data_path + "1000genomes_variants.vcf"
fasta_path = data_path + "linear_ref.fa"
gpc_path = "/home/knut/Documents/phd/graph_peak_caller/tests/mhc_test_data/chips_%s.json"
obg_base_name = data_path + "%s"
vcf_base_name = out_path + "%s_test"
chroms = [6]
