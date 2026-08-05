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
