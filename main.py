import pandas as pd


def parse_vcf(file_path):
    parsed_data = []

    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith("##") or line.startswith("#CHROM"):
                continue

            columns = line.split()

            # حماية من نقص الأعمدة
            if len(columns) < 8:
                continue

            parsed_data.append({
                "CHROM": columns[0],
                "POS": columns[1],
                "ID": columns[2],
                "REF": columns[3],
                "ALT": columns[4],
                "QUAL": columns[5],
                "FILTER": columns[6],
                "INFO": columns[7],
                "FORMAT": columns[8] if len(columns) > 8 else ".",
                "SAMPLES": columns[9] if len(columns) > 9 else "."
            })

    return parsed_data


def analyze_chromosomes(variants_list):
    chromosome_counts = {}
    for v in variants_list:
        chrom = v["CHROM"]
        chromosome_counts[chrom] = chromosome_counts.get(chrom, 0) + 1
    return chromosome_counts


def analyze_quality(variants_list):
    qual_values = []
    missing = 0

    for v in variants_list:
        q = v["QUAL"]
        if q != ".":
            try:
                qual_values.append(float(q))
            except:
                missing += 1
        else:
            missing += 1

    return qual_values, missing


def analyze_filters(variants_list):
    filter_counts = {}
    for v in variants_list:
        f = v["FILTER"]
        filter_counts[f] = filter_counts.get(f, 0) + 1
    return filter_counts


def analyze_mutations(variants_list):
    mutation_counts = {"Transition": 0, "Transversion": 0, "Other": 0}

    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for v in variants_list:
        ref = v["REF"]
        alt = v["ALT"]

        if len(ref) == 1 and len(alt) == 1:
            if (ref, alt) in transitions:
                mutation_counts["Transition"] += 1
            else:
                mutation_counts["Transversion"] += 1
        else:
            mutation_counts["Other"] += 1

    return mutation_counts


def main(file_path):

    variants_list = parse_vcf(file_path)
    df = pd.DataFrame(variants_list)

    print("\nPreview:")
    print(df.head())

    chrom = analyze_chromosomes(variants_list)
    qual_values, missing = analyze_quality(variants_list)
    filters = analyze_filters(variants_list)
    mutations = analyze_mutations(variants_list)

    print("\nChromosomes:", chrom)
    print("Mutations:", mutations)

    if qual_values:
        print("Max QUAL:", max(qual_values))
        print("Min QUAL:", min(qual_values))
        print("Avg QUAL:", sum(qual_values) / len(qual_values))

    print("Missing QUAL:", missing)
    print("Filters:", filters)

    with open("vcf_summary_report.txt", "w") as report:
        report.write("VCF Report\n\n")
        report.write(f"Chromosomes: {chrom}\n")
        report.write(f"Mutations: {mutations}\n")


if __name__ == "__main__":
    main("sample_vcf_file.vcf")