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
