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
st.markdown("""
<style>
[data-testid="stSidebar"] {
    background-color: #1E293B;
}

[data-testid="stSidebar"] * {
    color: white;
}
</style>
""", unsafe_allow_html=True)


#|Version 1 : VCF Analysis|

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
    plt.tight_layout()
    plt.ylabel("Chromosome")
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

st.title("Genova\n")
st.write("Welcome to the VCF analysis tool \n Version 1 : (AVIE) AML Variant Intelligent Explorer ")
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
#___***___***___***___***___***___***___***___***___***___***___***___***___***___***___***___***___***___***___***___***

# |Version 2 : Clinical Dataset Analysis|
#__________________________________
#Validating & Analyzing the dataset
#__________________________________

def validate_dataset(df):
        required_columns = ["AGE",
                            "SEX",
                            "RACE",
                            "FAB", "WBC",
                            "BM_BLAST_PERCENTAGE",
                            "PB_BLAST_PERCENTAGE",
                            "CYTOGENETICS",
                            "RISK_CYTO",
                            "RISK_MOLECULAR",
                            "MUTATION_COUNT",
                            "FRACTION_GENOME_ALTERED",
                            "TMB_NONSYNONYMOUS",
                            "SUB_CLONES",
                            "OS_MONTHS",
                            "OS_STATUS",
                            "DFS_MONTHS",
                            "DFS_STATUS",
                            "TRANSPLANT_TYPE",
                            "INDUCTION"]
        uploaded_columns = df.columns.tolist()
        missing_columns = []

        for column in required_columns:
          if column not in uploaded_columns:
            missing_columns.append(column)
        return missing_columns




def dataset_overview(df):
    patients = df.shape[0]
    features = df.shape[1]
    missing_values = df.isnull().sum()
    return patients, features, missing_values


def patient_characteristics(df):
   age = df["AGE"]
   min_age = age.min()
   max_age = age.max()
   avg_age = age.mean()

   sex = df["SEX"]
   sex_distribution = sex.value_counts()

   race = df["RACE"]
   race_distribution = race.value_counts()
   return  min_age, max_age, avg_age, sex_distribution, race_distribution



def disease_characteristics(df):
    fab = df["FAB"]
    fab_distribution = fab.value_counts()

    wbc = df["WBC"]
    min_wbc = wbc.min()
    max_wbc = wbc.max()
    avg_wbc = wbc.mean()

    bm_blast_percentage = df["BM_BLAST_PERCENTAGE"]
    min_bm_per = bm_blast_percentage.min()
    max_bm_per = bm_blast_percentage.max()
    avg_bm_per = bm_blast_percentage.mean()

    pb_blast_percentage = df["PB_BLAST_PERCENTAGE"]
    min_pb_per = pb_blast_percentage.min()
    max_pb_per = pb_blast_percentage.max()
    avg_pb_per = pb_blast_percentage.mean()
    return  (fab_distribution, min_wbc, max_wbc, avg_wbc, min_bm_per, max_bm_per, avg_bm_per ,
    min_pb_per,  max_pb_per,  avg_pb_per)



def risk_characteristics(df):
    cytogenetics = df["CYTOGENETICS"]
    cyto_distribution = cytogenetics.value_counts()

    risk_cyto = df["RISK_CYTO"]
    risk_cyto_distribution = risk_cyto.value_counts()

    risk_molecular = df["RISK_MOLECULAR"]
    risk_molecular_distribution = risk_molecular.value_counts()
    return cyto_distribution,  risk_cyto_distribution,  risk_molecular_distribution



def genomic_characteristics(df):
    mutation_count = df["MUTATION_COUNT"]
    min_count = mutation_count.min()
    max_count = mutation_count.max()
    avg_count = mutation_count.mean()

    frc_genome_alt = df["FRACTION_GENOME_ALTERED"]
    min_frc = frc_genome_alt.min()
    max_frc = frc_genome_alt.max()
    avg_frc = frc_genome_alt.mean()

    tmb = df["TMB_NONSYNONYMOUS"]
    min_tmb = tmb.min()
    max_tmb = tmb.max()
    avg_tmb = tmb.mean()


    sub_clones = df["SUB_CLONES"]
    sub_clones_distribution = sub_clones.value_counts()

    return (min_count, max_count, avg_count, min_frc, max_frc, avg_frc, min_tmb, max_tmb, avg_tmb, sub_clones_distribution)



def clinical_outcome(df):
    os_months = df["OS_MONTHS"]
    min_os_months = os_months.min()
    max_os_months = os_months.max()
    avg_os_months = os_months.mean()

    os_status = df["OS_STATUS"]
    os_status_distribution = os_status.value_counts()

    dfs_months = df["DFS_MONTHS"]
    min_dfs_months = dfs_months.min()
    max_dfs_months = dfs_months.max()
    avg_dfs_months = dfs_months.mean()

    dfs_status = df["DFS_STATUS"]
    dfs_status_distribution = dfs_status.value_counts()

    transplant_type = df["TRANSPLANT_TYPE"]
    transplant_distribution = transplant_type.value_counts()

    induction = df["INDUCTION"]
    induction_distribution = induction.value_counts()

    return (min_os_months, max_os_months, avg_os_months, os_status_distribution, min_dfs_months,
             max_dfs_months,avg_dfs_months, transplant_distribution,induction_distribution, dfs_status_distribution)


#________________________________________________________________________________________________________________________

uploaded_dataset = st.file_uploader("Upload your dataset here :",  type=["csv", "txt" , "tsv"])

if uploaded_dataset is not None:
    df = pd.read_csv(uploaded_dataset , sep="\t")
    missing_columns = validate_dataset(df)
    if len(missing_columns) > 0 :
        st.warning(f"Missing Columns : {missing_columns}")
    else :
        st.success("✅ Dataset Validation Passed All required colums are available.")

        section = st.sidebar.radio(
            "Navigate to",
            [
                "Dataset Overview",
                "Patient Characteristics",
                "Disease Characteristics",
                "Risk Characteristics",
                "Genomic Characteristics",
                "Clinical Outcome",
            ]
        )
        if section == "Dataset Overview":
            patients, features, missing_values = dataset_overview(df)
            st.header("📊 Dataset Overview")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("patients", patients)
            with col2:
                st.metric("features", features)
            with col3:
                st.metric("missing_values", missing_values.sum())

        elif section == "Patient Characteristics":
            min_age, max_age, avg_age, sex_distribution, race_distribution = patient_characteristics(df)
            st.header("Patient Characteristics")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum age", min_age)
            with col2:
                st.metric("maximum age", max_age)
            with col3:
                st.metric("average age", avg_age)
            st.subheader("Sex distribution")
            st.bar_chart(sex_distribution)
            st.subheader("Race distribution")
            st.bar_chart(race_distribution)

        elif section == "Disease Characteristics":
            fab_distribution, min_wbc, max_wbc, avg_wbc, min_bm_per, max_bm_per, avg_bm_per, min_pb_per, max_pb_per, avg_pb_per = disease_characteristics(
                df)
            st.header("📊 Disease Characteristics")
            st.subheader("FAB distribution")
            st.bar_chart(fab_distribution)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum wbc", min_wbc)
            with col2:
                st.metric("maximum wbc", max_wbc)
            with col3:
                st.metric("average wbc", avg_wbc)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum bm percentage", min_bm_per)
            with col2:
                st.metric("maximum bm percentage", max_bm_per)
            with col3:
                st.metric("average bm percentage", avg_bm_per)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum pb percentage", min_pb_per)
            with col2:
                st.metric("maximum pb percentage", max_pb_per)
            with col3:
                st.metric("average pb percentage", avg_pb_per)

        elif section == "Risk Characteristics":
            cyto_distribution, risk_cyto_distribution, risk_molecular_distribution = risk_characteristics(df)
            st.header("Risk characteristics")
            st.subheader("Cytogenetics distribution")
            st.bar_chart(cyto_distribution)
            st.subheader("Risk cyto distribution")
            st.bar_chart(risk_cyto_distribution)
            st.subheader("Risk molecular distribution")
            st.bar_chart(risk_molecular_distribution)

        elif section == "Genomic Characteristics":
            min_count, max_count, avg_count, min_frc, max_frc, avg_frc, min_tmb, max_tmb, avg_tmb, sub_clones_distribution = genomic_characteristics(
                df)
            st.header("Genomic characteristics")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum count", min_count)
            with col2:
                st.metric("maximum count", max_count)
            with col3:
                st.metric("average count", avg_count)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum fraction altered count", min_frc)
            with col2:
                st.metric("maximum fraction", max_frc)
            with col3:
                st.metric("average fraction", avg_frc)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum tmb", min_tmb)
            with col2:
                st.metric("maximum tmb", max_tmb)
            with col3:
                st.metric("average tmb", avg_tmb)
            col1, col2, col3 = st.columns(3)
            st.subheader("sub clones distribution")
            st.bar_chart(sub_clones_distribution)

        elif section == "Clinical Outcome":

            min_os_months, max_os_months, avg_os_months, os_status_distribution, min_dfs_months, max_dfs_months, avg_dfs_months, transplant_distribution, induction_distribution, dfs_status_distribution = clinical_outcome(
                df)
            st.header("Clinical Outcome")
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum OS months", min_os_months)
            with col2:
                st.metric("maximum OS months", max_os_months)
            with col3:
                st.metric("average OS months", avg_os_months)
            st.subheader("OS status distribution")
            st.bar_chart(os_status_distribution)
            col1, col2, col3 = st.columns(3)
            with col1:
                st.metric("minimum DFS months", min_dfs_months)
            with col2:
                st.metric("maximum DFS months", max_dfs_months)
            with col3:
                st.metric("average DFS months", avg_dfs_months)
            st.subheader("DFS status distribution")
            st.bar_chart(dfs_status_distribution)
            st.subheader("Transplant distribution")
            st.bar_chart(transplant_distribution)
            st.subheader("Induction distribution")
            st.bar_chart(induction_distribution)







