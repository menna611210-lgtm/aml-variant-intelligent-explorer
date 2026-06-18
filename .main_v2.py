import matplotlib.pyplot as plt
import pandas as pd
import os

print("AVIE v1.0 \n AML Variant Intelligent Explorer")

# ******************
# File Analysis
# ******************
while True:
    file_path = input("\nEnter VCF File Path : ")

    if not file_path.endswith(".vcf"):
        print("This is not a VCF file")
        continue

    if not os.path.exists(file_path):
        print("The File Does Not Exist")
        continue

    elif file_path.endswith(".vcf"):
        file_name = file_path.replace(".vcf" , "")
        print("Parsing VCF File ....")
        def parse_vcf(file_path):
            parsed_data = []
            with open(file_path, "r") as file :
                for line in file:
                    line = line.strip()
                    if line.startswith("##"):
                      continue
                    if line.startswith("#CHROM"):
                      continue
                    columns = line.split("\t")
                    if len(columns) < 8 :
                        continue
                    parsed_data.append({
                       "CHROM" : columns[0],
                        "POS" : int(columns[1]),
                        "ID" : columns[2],
                        "REF" : columns[3],
                        "ALT" : columns[4],
                        "QUAL" : columns[5],
                        "FILTER" : columns[6],
                        "INFO" : columns[7],
                        "FORMAT" : columns[8],
                        "SAMPLES" : columns[9]
                    })
            return parsed_data

# *************************
# 2. CHROMOSOME ANALYSIS
# *************************
print("Analyzing Chromosomes .....")
def analyze_chromosomes(variants):
    counts = {}
    for v in variants:
        chrom = v["CHROM"]
        counts[chrom] = counts.get(chrom, 0) + 1
    return counts


# *************************
# 3. QUALITY ANALYSIS
# *************************
print("Analyzing Quality SCore ....")
def analyze_quality(variants):
    values = []
    missing = 0

    for v in variants:
        q = v["QUAL"]
        try:
            values.append(float(q))
        except:
            missing += 1

    return values, missing


# *************************
# 4. FILTER ANALYSIS
# *************************
print("Analyzing Mutation Types ....")
def analyze_filters(variants):
    counts = {}
    for v in variants:
        f = v["FILTER"]
        counts[f] = counts.get(f, 0) + 1
    return counts


# **************************
# 5. MUTATION CLASSIFICATION
# **************************
print("Analyzing Mutation Types ....")
def analyze_mutations(variants):
    result = {
        "Transition": 0,
        "Transversion": 0,
        "Insertion": 0,
        "Deletion": 0,
        "Other": 0
    }

    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for v in variants:
        ref = v["REF"]
        alt = v["ALT"].split(",")[0]  # handle multi-allelic safely

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


# *************************
# 6. PLOTTING
# *************************
print("Generating Plots ....")
def plot_chromosome_distribution(counts):
    plt.figure(figsize=(5, 2))
    plt.barh(counts.keys(), counts.values(), color="skyblue", edgecolor="black")
    plt.title("Chromosome Variant Distribution")
    plt.xlabel("Count")
    plt.ylabel("Chromosome")
    plt.tight_layout()
    plt.savefig("results/chromosome_distribution.png" , dpi=300)
    plt.show()


def plot_quality_distribution(qual_values):
    plt.figure(figsize=(5, 2))
    plt.hist(qual_values, bins=30, color="violet", edgecolor="black")
    plt.title("Quality Distribution")
    plt.xlabel("Frequency")
    plt.ylabel("Quality Score")
    plt.tight_layout()
    plt.savefig("results/quality_distribution.png", dpi=300)
    plt.show()


def plot_filter_distribution(filter_values):
    plt.figure(figsize=(5, 2))
    plt.barh(filter_values.keys(), filter_values.values(), color="green", edgecolor="black")
    plt.title("Filter's Quality Distribution")
    plt.xlabel("Values")
    plt.ylabel("Filter's Quality")
    plt.tight_layout()
    plt.savefig("results/filter_distribution.png", dpi=300)
    plt.show()


def plot_mutation_distribution(mutations):
    plt.figure(figsize=(5, 2))
    plt.barh(mutations.keys(), mutations.values(), color="blue", edgecolor="black")
    plt.title("Mutation Distribution")
    plt.xlabel("Values")
    plt.ylabel("Mutation Type")
    plt.tight_layout()
    plt.savefig("results/mutation_distribution.png", dpi=300)
    plt.show()


# *************************
# 7. REPORT GENERATION
# *************************
print("Generating Report ....")
def write_report(path, total, chrom, mutations, qual, missing, filters):
    with open(path, "w") as report:
        report.write("results/" + file_name + "REPORT\n\n")

        report.write(f"Total variants: {total}\n\n")

        report.write("Chromosome distribution:\n")
        for k, v in chrom.items():
            report.write(f"{k}: {v}\n")

        report.write("\nMutation types:\n")
        for k, v in mutations.items():
            report.write(f"{k}: {v}\n")

        report.write("\nQuality metrics:\n")
        if qual:
            report.write(f"Max: {max(qual)}\n")
            report.write(f"Min: {min(qual)}\n")
            report.write(f"Avg: {sum(qual) / len(qual)}\n")
        else:
            report.write("No quality data\n")

        report.write(f"\nMissing QUAL: {missing}\n")

        report.write("\nFilters:\n")
        for k, v in filters.items():
            report.write(f"{k}: {v}\n")


# *************************
# 8. MAIN PIPELINE
# *************************
print("Analysis Completed Successfully :)")
def main(file_path):
    variants = parse_vcf(file_path)

    df = pd.DataFrame(variants)
    print("\nPreview:")
    print(df.head())

    chrom = analyze_chromosomes(variants)
    qual, missing = analyze_quality(variants)
    filters = analyze_filters(variants)
    mutations = analyze_mutations(variants)

    os.makedirs("results" , exist_ok = True)
    total = len(variants)

    print("\nChromosomes:", chrom)
    print("Mutations:", mutations)

    if qual:
        print("Max QUAL:", max(qual))
        print("Min QUAL:", min(qual))
        print("Avg QUAL:", sum(qual) / len(qual))

    print("Missing QUAL:", missing)
    print("Filters:", filters)

    plot_chromosome_distribution(chrom)

    write_report(
        "results/" + file_name + "report.txt",
        total,
        chrom,
        mutations,
        qual,
        missing,
        filters
    )


if __name__ == "__main__":
    main("sample_vcf_file.vcf")
