"""
Simple methods for accessing the Togows API (http://togows.org/)
"""
import urllib.request as urllib2
import urllib
import xml.etree.ElementTree as ET


def parse_ucsc_xml(text):
    root = ET.fromstring(text)
    elem = root.find("SEQUENCE").find("DNA")
    seq = "".join(elem.text.split())
    return seq


def get_sequence(loci_id, start=1, end=0, caching=True):
    return get_sequence_ucsc(loci_id, start, end, caching)
    """
    Gets the sequence (as fasta format)
    of an alternative loci (e.g. chr1 or chr1_GL383519v1_alt)
    """
    url = "http://togows.org/api/ucsc/hg38/%s:%d-%d.fasta" % (loci_id, start, end)
    print("Fetching sequence from " + url)

    try:
        sequence = urllib2.urlopen(url).read()
    except urllib.error.HTTPError:
        print(url)
        raise

    # Remove line breaks
    sequence = sequence.decode('utf8').replace('\n', '')

    return sequence


def get_sequence_ucsc(loci_id, start=1, end=0, caching=True):

    file_name = "data/tmp/sequence_%s_%s_%s.fasta" % (loci_id, start, end)

    if caching:
        import os
        assert os.path.isdir("data/tmp"), "data/tmp folder must exist"
        if os.path.isfile(file_name):
            f = open(file_name)
            seq = f.read()
            f.close()
            return seq

    url = "http://genome.ucsc.edu/cgi-bin/das/hg38/dna?segment=%s:%d,%d" % (loci_id, start, end)
    data = urllib2.urlopen(url).read()
    try:
        sequence = parse_ucsc_xml(data)
    except:
        print(url)
        raise

    if caching:
        f = open(file_name, "w")
        f.write(sequence)
        f.close()

    return sequence


def save_sequence_to_fasta(loci_id, start, end, file_name=""):
    """
    Collects a sequnce from the togows API, and saves it to file (fasta format)
    :param loci_id: The id of the loci (e.g. chr1 or chr11_KI270827v1_alt)
    :param start: Start position >= 1
    :param end: End position (inclusive (?))
    :param file_name: Alternative file name. If empty, one will be generated
    :return: Returns the file name
    """

    start += 1
    end += 1

    if file_name == "":
        file_name = DATA_PATH + "%s:%d-%d.fasta" % (loci_id, start, end)


    import os.path
    if os.path.isfile(file_name):
        return file_name

    curpath = os.path.abspath(os.curdir)
    f = open(file_name, "w")
    f.write(get_sequence(loci_id, start, end))
    return file_name

if __name__ == "__main__":
    xml = """<DASDNA>
    <SEQUENCE id="chr1" start="1" stop="100" version="1.00">
    <DNA length="100">
    nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
    </DNA>
    </SEQUENCE>
    </DASDNA>"""
    parse_ucsc_xml(xml)
    
