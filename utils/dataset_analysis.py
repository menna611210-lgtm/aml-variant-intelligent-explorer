

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



