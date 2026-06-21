import streamlit as st
import base64
import matplotlib.pyplot as plt
import pandas as pd
import os
from fpdf import FPDF


def get_base64_image(image_path):
    with open(image_path, "rb") as image_file:
        return base64.b64encode(image_file.read()).decode()
bg_base64 = get_base64_image("img.png")
st.markdown(
    f"""
    <style>
    .stApp {{
        background-image: url("data:image/jpg;base64,{bg_base64}");
        background-size: cover;
        background-position: center;
        background-repeat: no-repeat;
        background-attachment: fixed;
        color : white;
    }}
    h1, h2, h3, h4, h5, h6, p, .stMarkdown, [data-testid="stWidgetLabel"] p {{
        color: #ffffff !important;
    }}
    [data-testid="stFileUploader"] section {{
        background-color: rgba(255, 255, 255, 0.95) !important; 
    }}
    [data-testid="stFileUploader"] * {{
        color:#111111 !important; 
    }}
    [data-testid="stFileUploader"] button {{
        background-color: #f0f2f6 !important;
        border: 1px solid #ccced2 !important;
        color:#111111 !important;
    }}
    </style>
    """,
    unsafe_allow_html=True
)


# ___________________________
# Generating PDF File Report
# ___________________________

def generate_pdf(dataframe, chrom, qual, mutations, filters):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_font("Arial","B" , size=16)
    pdf.cell(0, 10, "AVIE Analysis Report" , ln=True , align="C")
    pdf.ln(5)
    pdf.set_font("Arial" , "B" , 14)
    pdf.cell(0, 10, "1. Chromosome Variant Distribution:", ln=True)
    pdf.ln()
    for key, value in chrom.items():
        pdf.cell(0, 8, f" {key}: {value}", ln=True)
        pdf.ln()
    pdf.image(
      "results/chromosome_distribution.png", x=15 , w=150)


    pdf.add_page()
    pdf.set_font("Arial", "B", size=16)
    pdf.cell(0, 10, "Quality Score Metrics:" , ln=True)
    pdf.ln()
    if qual :
      pdf.cell(0, 10, f"-Max Quality: {max(qual)}" , ln=True)
      pdf.cell(0, 10, f"-Min Quality: {min(qual)}" , ln=True)
      pdf.cell(0, 10, f"-Average Quality: {sum(qual) / len(qual)}" , ln=True)
      pdf.ln()
      pdf.image(
          "results/quality_distribution.png", x=15 , w=150)
    else :
      pdf.cell(0, 10, "No Quality Data Available" , ln=True)

    pdf.add_page()
    pdf.set_font("Arial", "B", size=16)
    pdf.cell(0, 10, "Mutations Type:" , ln=True)
    pdf.ln()
    for key, value in mutations.items():
       pdf.cell(0, 10, f"{key}: {value}\n" , ln=True)
    pdf.ln()
    pdf.image(
      "results/mutation_distribution.png", x=15 , w=150)

    pdf.add_page()
    pdf.set_font("Arial", "B", size=16)
    pdf.cell(0, 10, "Filter Quality:" , ln=True)
    pdf.ln()
    for key, value in filters.items():
       pdf.cell(0, 10, f"{key}: {value}\n" , ln=True)
    pdf.ln()
    pdf.image(
       "results/filter_distribution.png", x=15 , w=150 )
    return bytes(pdf.output(dest = "S"))

# __________________________
# File Analysis
# __________________________
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



# ____________________
# Streamlit Edition
# ____________________

st.title("AVIE v1.0 \n🧬 AML Variant Intelligent Explorer\n")
st.write("Welcome to the VCF analysis tool")
uploaded_file = st.file_uploader("Upload your VCF file here :", type=["VCF", "txt"])

if uploaded_file is not None:
    st.success("VCF File Uploaded Successfully !")
    st.write("File Name : ", uploaded_file.name)
    variants = parse_vcf(uploaded_file)
    df = pd.DataFrame(variants)
    st.write("VCF Preview Data")
    st.dataframe(df)

    chrom = analyze_chromosomes(variants)
    qual , missing = analyze_quality(variants)
    filters = analyze_filters(variants)
    mutations = analyze_mutations(variants)


    plot_chromosome_distribution(chrom)
    plot_quality_distribution(qual)
    plot_mutation_distribution(mutations)
    plot_filter_distribution(filters)

    pdf_data = generate_pdf(df, chrom, qual, mutations, filters)
    st.write("### Visual Dashboard ###")
    st.image("results/chromosome_distribution.png")
    if qual:
        st.image("results/quality_distribution.png")

    st.download_button(
        label="📥 Download Full Analysis Report (PDF)",
        data=pdf_data,
        file_name="AVIE_Analysis_Report.pdf",
        mime="application/pdf"
    )
    print("M")