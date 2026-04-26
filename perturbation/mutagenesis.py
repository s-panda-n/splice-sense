NUCLEOTIDES = ['A', 'T', 'G', 'C']

def mutate_sequence(sequence, position, alt):
    assert 0 <= position < len(sequence)
    assert alt in NUCLEOTIDES
    return sequence[:position] + alt + sequence[position+1:]

def get_alternatives(sequence, position):
    ref = sequence[position].upper()
    return [n for n in NUCLEOTIDES if n != ref]

def type1_mutations(sequence, site_type):
    """Consensus GT/AG position mutations — model SHOULD flip here."""
    mutations = []
    positions = [200, 201] if site_type == "donor" else [198, 199]
    for pos in positions:
        for alt in get_alternatives(sequence, pos):
            mutations.append({
                "type": "Type1", "position": pos,
                "ref": sequence[pos], "alt": alt,
                "mutated_seq": mutate_sequence(sequence, pos, alt)
            })
    return mutations

def type2_mutations(sequence, n_positions=5, seed=42):
    """ESE region mutations — exonic context near splice site."""
    import random
    random.seed(seed)
    positions = random.sample(range(170, 200), n_positions)
    mutations = []
    for pos in positions:
        for alt in get_alternatives(sequence, pos):
            mutations.append({
                "type": "Type2", "position": pos,
                "ref": sequence[pos], "alt": alt,
                "mutated_seq": mutate_sequence(sequence, pos, alt)
            })
    return mutations

def type3_mutations(sequence, n_positions=5, seed=42):
    """Deep intronic mutations — model should NOT flip here."""
    import random
    random.seed(seed)
    positions = random.sample(range(210, 380), n_positions)
    mutations = []
    for pos in positions:
        for alt in get_alternatives(sequence, pos):
            mutations.append({
                "type": "Type3", "position": pos,
                "ref": sequence[pos], "alt": alt,
                "mutated_seq": mutate_sequence(sequence, pos, alt)
            })
    return mutations

def get_all_mutations(sequence, site_type):
    return (
        type1_mutations(sequence, site_type) +
        type2_mutations(sequence) +
        type3_mutations(sequence)
    )
