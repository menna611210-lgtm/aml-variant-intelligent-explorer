import streamlit as st
import matplotlib.pyplot as plt
import pandas as pd
import os

st.title("AVIE v1.0 \n🧬 AML Variant Intelligent Explorer\n")
st.write("Welcome to the VCF analysis tool")
uploaded_file = st.file_uploader("Upload your VCF file here :", type=["VCF", "txt"])

if uploaded_file is not None:
    st.success("VCF File Uploaded Successfully ")
    st.write("File Name : ", uploaded_file.name)
st.subheader(" VCF Preview ")


# __________________________
# File Analysis
# __________________________
def parse_vcf(uploaded_file):
    parsed_data = []
    for line in uploaded_file :
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


if uploaded_file is not None:
    st.success("VCF File Uploaded Successfully ")
    st.write("File Name : ", uploaded_file.name)
    variants = parse_vcf(uploaded_file)
    df = pd.DataFrame(variants)
    st.write(df)
    st.subheader(" VCF Preview ")

# __________________________
# 2. CHROMOSOME ANALYSIS
# __________________________

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


# _________________________
# 6. plotting
# _________________________

os.makedirs("results", exist_ok=True)


def plot_chromosome_distribution(counts):
    plt.figure(figsize=(5, 2))
    plt.barh(counts.keys(), counts.values(), color="skyblue", edgecolor="black")
    plt.title("Chromosome Variant Distribution")
    plt.xlabel("Count")
    plt.ylabel("Chromosome")
    plt.tight_layout()
    plt.savefig("results/chromosome_distribution.png", dpi=300)
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


# __________________________
# 7. REPORT GENERATION
# __________________________

def write_report(file_path, total, chrom, mutations, qual, missing, filters):
    with open(file_path, "w") as report:
        report.write("AVIE Analysis Report\n\n")

        report.write(f"Total variants: {total}\n\n")

        report.write("Chromosome distribution:\n")
        for key, value in chrom.items():
            report.write(f"{key}: {value}\n")

        report.write("\nMutation types:\n")
        for key, value in mutations.items():
            report.write(f"{key}: {value}\n")

        report.write("\nQuality metrics:\n")
        if qual:
            report.write(f"Max: {max(qual)}\n")
            report.write(f"Min: {min(qual)}\n")
            report.write(f"Avg: {sum(qual) / len(qual)}\n")
        else:
            report.write("No quality data\n")

        report.write(f"\nMissing QUAL: {missing}\n")

        report.write("\nFilters:\n")
        for key, value in filters.items():
            report.write(f"{key}: {value}\n")


# __________________________
# 8. MAIN PIPELINE
# __________________________

def main(file_path):
    variants = parse_vcf(file_path)

    df = pd.DataFrame(variants)
    print("\nPreview:")
    print(df.head())
    print("Analyzing Chromosomes .....")
    chrom = analyze_chromosomes(variants)
    print("Analyzing Quality SCore ....")
    qual, missing = analyze_quality(variants)
    print("Analyzing Mutation Types ....")
    filters = analyze_filters(variants)
    print("Analyzing Mutation Types ....")
    mutations = analyze_mutations(variants)
    print("Generating Report ....")
    print("Analysis Completed Successfully :)")

    plot_chromosome_distribution(chrom)
    plot_quality_distribution(qual)
    plot_filter_distribution(filters)
    plot_mutation_distribution(mutations)

    total = len(variants)

    print("\nChromosomes:", chrom)
    print("Mutations:", mutations)

    if qual:
        print("Max QUAL:", max(qual))
        print("Min QUAL:", min(qual))
        print("Avg QUAL:", sum(qual) / len(qual))

    print("Missing QUAL:", missing)
    print("Filters:", filters)

    write_report(
        "results/" + uploaded_file.name + "report.txt",
        total,
        chrom,
        mutations,
        qual,
        missing,
        filters
    )


if __name__ == "__main__":
    main(uploaded_file)
