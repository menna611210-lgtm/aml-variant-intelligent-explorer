def analyze_chromosomes(variants):
    counts = {}
    for variant in variants:
        chrom = variant["CHROM"]
        counts[chrom] = counts.get(chrom, 0) + 1
    return counts


# __________________________
# 3. QUALITY ANALYSIS
# __________________________


def analyze_quality(variants):
    values = []
    missing = 0

    for variant in variants:
        quality = variant["QUAL"]
        try:
            values.append(float(quality))
        except:
            missing += 1

    return values, missing


# __________________________
# 4. FILTER ANALYSIS
# __________________________


def analyze_filters(variants):
    counts = {}
    for variant in variants:
        filter = variant["FILTER"]
        counts[filter] = counts.get(filter, 0) + 1
    return counts


# __________________________
# 5. MUTATION CLASSIFICATION
# __________________________

def analyze_mutations(variants):
    result = {
        "Transition": 0,
        "Transversion": 0,
        "Insertion": 0,
        "Deletion": 0,
        "Other": 0
    }

    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for variant in variants:
        ref = variant["REF"]
        alt = variant["ALT"].split(",")[0]  # handle multi-allelic safely

        if len(ref) == 1 and len(alt) == 1:
            if (ref, alt) in transitions:
                result["Transition"] += 1
            else:
                result["Transversion"] += 1

        elif len(ref) > len(alt):
            result["Deletion"] += 1

        elif len(ref) < len(alt):
            result["Insertion"] += 1

        else:
            result["Other"] += 1

    return result
