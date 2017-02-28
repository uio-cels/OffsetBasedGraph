from .interval import Interval
from .translation import Translation


def align_cigar(cigar, main_interval, alt_interval, graph):
    """Merge main interval and alt interval based on cigar string

    :param cigar: cigar string [(code, n), (code, n),,,(code, n)]
    :param main_interval: Interval on main path
    :param alt_interval: Interval on alt
    :param graph: current graph
    :returns: Translation object for the cigar string
    :rtype: Translation

    """
    main_offset = main_interval.start_position.offset
    alt_offset = alt_interval.start_position.offset

    main_id = main_interval.region_paths[0]
    alt_id = alt_interval.region_paths[0]

    def MainInterval(start, end):
        return Interval(start, end, [main_id], graph=graph)

    def AltInterval(start, end):
        return Interval(start, end, [alt_id], graph=graph)

    alt_len = alt_interval.length()

    main_rps = []
    alt_rps = []
    b_to_a = {}
    _id = graph._next_id()
    #print("#", _id)

    main_n = 0
    alt_n = 0
    if main_offset > 0:
        main_rps.append(_id)
        b_to_a[_id] = [MainInterval(0, main_offset)]
        _id += 1
    if alt_offset > 0:
        alt_rps.append(_id)
        b_to_a[_id] = [AltInterval(0, alt_offset)]
        _id += 1

    for var_type, n in cigar:
        if var_type == "D":
            main_rps.append(_id)
            b_to_a[_id] = [MainInterval(main_offset, main_offset+n)]
            main_n = n
            main_offset += n

        elif var_type == "I":
            alt_rps.append(_id)
            b_to_a[_id] = [AltInterval(alt_offset, alt_offset+n)]
            alt_n = n
            alt_offset += n

        elif var_type == "M":
            alt_rps.append(_id)
            main_rps.append(_id)
            b_to_a[_id] = [MainInterval(main_offset, main_offset+n),
                           AltInterval(alt_offset, alt_offset+n)]
            main_n = n
            alt_n = n
            alt_offset += n
            main_offset += n

        elif var_type == "V":
            id_a = _id
            _id += 1
            alt_rps.append(id_a)
            main_rps.append(_id)
            b_to_a[_id] = [MainInterval(main_offset, main_offset+n)]
            b_to_a[id_a] = [AltInterval(alt_offset, alt_offset+n)]
            main_n = n
            alt_n = n
            alt_offset += n
            main_offset += n
        assert alt_offset <= alt_interval.end_position.offset
        _id += 1

    if main_offset < graph.blocks[main_id].length():
        main_rps.append(_id)
        b_to_a[_id] = [MainInterval(main_offset,
                                    graph.blocks[main_id].length())]

        main_n = graph.blocks[main_id].length()-main_offset
        _id += 1

    alt_len = graph.blocks[alt_id].length()
    if alt_offset < alt_len:
        alt_rps.append(_id)
        b_to_a[_id] = [AltInterval(alt_offset,
                                   alt_len)]  # graph.blocks[alt_id].length())]
        alt_n = alt_len-alt_offset
        _id += 1

    a_to_b = {alt_id: [Interval(0, alt_n, alt_rps)],
              main_id: [Interval(0, main_n, main_rps)]}

    return Translation(a_to_b, b_to_a, graph=graph)


def get_match_cigar(seq1, seq2):
    """Create cigar string for matched sequence
    V for variation
    M for identity

    :param seq1: str
    :param seq2: str
    :returns: cigar string
    :rtype: list(tuple(char, int))

    """
    prev = None
    changes = []
    n = 0
    for a, b in zip(seq1, seq2):
        cur = (a == b)
        if cur != prev:
            changes.append(n)
        n += 1
        prev = cur

    cigar = []
    changes.append(len(seq1))
    cur_symbol = "M" if seq1[0] == seq2[0] else "V"
    for start, end in zip(changes[:-1], changes[1:]):
        cigar.append((cur_symbol, end-start))
        cur_symbol = "V" if cur_symbol == "M" else "M"

    return cigar


def cigar_lens(cigar):
    main_len = 0
    alt_len = 0
    for c, n in cigar:
        if c in ("M", "V", "D"):
            main_len += n
        if c in ("M", "V", "I"):
            alt_len += n
    return main_len, alt_len


def clean_cigar(cigar, alt_seq, main_seq):
    """Clean up and create full cigar string (convetrt M to M or V)

    :param cigar: list of cigar codes
    :param alt_seq: sequence on alt
    :param main_seq: sequence on main
    :param alt_start: start offset on alt
    :param main_start: start offset on main
    :returns: cleaned up cigar string
    :rtype: list(tuple(char, int))

    """
    cigar = cigar.split()
    cleaned_cigar = []
    alt_offset = 0
    main_offset = 0
    for c in cigar:
        first = c[0]
        if first == "M" or first == "D" or first == "I":
            var_type = first
            n = int(c[1:])
        else:
            raise ValueError("Invalid cigar string. Should only contain letters I, D or M.")

        if var_type == "I":
            alt_offset += n
            cleaned_cigar.append((var_type, n))
            continue

        if var_type == "D":
            main_offset += n
            cleaned_cigar.append((var_type, n))
            continue
        seq_on_alt = alt_seq[alt_offset:alt_offset+n]
        seq_on_main = main_seq[main_offset:main_offset+n]
        match_cigar = get_match_cigar(seq_on_main, seq_on_alt)
        cleaned_cigar.extend(match_cigar)
        main_offset += n
        alt_offset += n
    lens = cigar_lens(cleaned_cigar)
    seq_lens = (len(main_seq), len(alt_seq))
    assert lens == seq_lens, "%s %s" % (lens, seq_lens)
    return cleaned_cigar
