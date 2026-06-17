import pandas as pd
import matplotlib.pyplot as plt

#parse vcf file and extract variant information
def parse_vcf(file_path):
    parsed_data = []

#Open the file and start analyzing
    with open(file_path, "r") as file:
        for line in file:
            line = line.strip()

            if not line or line.startswith("##") or line.startswith("#CHROM"):
                continue

            columns = line.split()

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


#counts variants per chromosome
def analyze_chromosomes(variants_list):
    chromosome_counts = {}
    for v in variants_list:
        chrom = v["CHROM"]
        chromosome_counts[chrom] = chromosome_counts.get(chrom, 0) + 1
    return chromosome_counts

#determines the value of the quality and counts it
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


#determines the result of passing on filters
def analyze_filters(variants_list):
    filter_counts = {}
    for v in variants_list:
        f = v["FILTER"]
        filter_counts[f] = filter_counts.get(f, 0) + 1
    return filter_counts


#the type of the mutations
def analyze_mutations(variants_list):
    mutation_counts = {"Transition": 0, "Transversion": 0, "Insertion": 0, "Deletion": 0, "Other": 0}

    transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

    for v in variants_list:
        ref = v["REF"]
        alt = v["ALT"]

        if len(ref) == 1 and len(alt) == 1:
            if (ref, alt) in transitions:
                mutation_counts["Transition"] += 1
            else:
                mutation_counts["Transversion"] += 1
        elif len(ref) > 1 and len(alt) == 1:
            mutation_counts["Deletion"] +=1
        elif len(ref)== 1 and len (alt) > 1:
            mutation_counts["Insertion"] += 1
        else:
            mutation_counts["Other"] += 1

    return mutation_counts

#the output of the parsing
def main(file_path):

    variants_list = parse_vcf(file_path)
    df = pd.DataFrame(variants_list)

    print("\nPreview:")
    print(df.head())

    chrom = analyze_chromosomes(variants_list)
    qual_values, missing = analyze_quality(variants_list)
    filters = analyze_filters(variants_list)
    mutations = analyze_mutations(variants_list)
    total_variants = len(variants_list)

    print("\nChromosomes:", chrom)
    print("Mutations:", mutations)

    if qual_values:
        print("Max QUAL:", max(qual_values))
        print("Min QUAL:", min(qual_values))
        print("Avg QUAL:", sum(qual_values) / len(qual_values))

    print("Missing QUAL:", missing)
    print("Filters:", filters)
#the external report
    with open("vcf_summary_report.txt", "w") as report:

        report.write("VCF Analysis Report\n\n")

        report.write(f"Total variants : {total_variants}\n")
        report.write(f"\nChromosomes distribution :\n")
        for key , value in chrom.items():
            report.write(f"chr{key}: {value}\n")
            x = [{key}]
            y = [value]
            plt.bar({key} , value , colour = 'skyblue' , edgecolor = 'black')
            plt.show()
        report.write("\n")

        report.write(f"\nMutations type : \n")
        for key , value in mutations.items():
            report.write(f"{key}: {value}\n")
        report.write("\n")

        report.write(f"\nQuality Matrics :\n")
        if qual_values :
            report.write(f"\nMax Qual : {max(qual_values)}\n")
            report.write(f"\nMin Qual : {min(qual_values)}\n")
            report.write(f"\nAVG Qual : {sum(qual_values) / len(qual_values)}")
        else:
            report.write("\nNo Quality data available\n")

if __name__ == "__main__":
    main("sample_vcf_file.vcf")

