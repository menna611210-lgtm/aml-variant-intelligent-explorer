def parse_vcf(uploaded_file):
    parsed_data = []
    for line in uploaded_file:
        line = line.decode("utf-8").strip()
        if line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            continue
        columns = line.split("\t")
        if len(columns) < 10:
            continue
        parsed_data.append(
            {"CHROM": columns[0],
             "POS": int(columns[1]),
             "ID": columns[2],
             "REF": columns[3],
             "ALT": columns[4],
             "QUAL": columns[5],
             "FILTER": columns[6],
             "INFO": columns[7],
             "FORMAT": columns[8],
             "SAMPLES": columns[9]})
    return parsed_data
