#استدعاء مكتبات pandas
import pandas as pd
# =========================================================
# 1. مرحلة تعريف الدالة (المصنع اللي بيقرأ وينظف البيانات)
# =========================================================
def parse_vcf(file_path):
    parsed_data = []  # السلة المحلية اللي بتجمع الطفرات جوه الدالة

    with open(file_path, "r") as file:
        for line in file:
            if line.startswith("##"):
                continue  # تخطي الميتا داتا
            elif line.startswith("#CHROM"):
                continue  # تخطي الهيدر الأساسي
            else:
                columns = line.split()
                parsed_data.append({
                    "CHROM": columns[0],
                    "POS": columns[1],
                    "ID": columns[2],
                    "REF": columns[3],
                    "ALT": columns[4],
                    "QUAL": columns[5],
                    "FILTER": columns[6],
                    "INFO": columns[7],
                    "FORMAT": columns[8],
                    "SAMPLES": columns[9]
                })

    return parsed_data, #تسليم السلة المليانة بره ال function
# =========================================================
# 3. تشغيل ال function و تحويل النتيجه لجدول pandas
# =========================================================
variants_list = parse_vcf("sample_vcf_file.vcf")
df = pd.DataFrame(variants_list)

# الإعدادات السحرية لعرض كل الأعمدة
pd.set_option('display.max_columns', None)
pd.set_option('display.max_rows' , None)
pd.set_option('display.width', 1000)
# أمر طباعة الجدول فى ال terminal
print("\n=== Pandas DataFrame Preview ===\n")
print(df)


# =========================================================
# 2. مرحلة تشغيل ال function (تنظيف وقراءة ملف الجينات الأصلي)
# =========================================================
# بننادي على ال function ونديها ملف الجينات، والنتيجة بتتخزن في قائمة بره اسمها variants_list
variants_list = parse_vcf("sample_vcf_file.vcf")

# =========================================================
# 3. مرحلة الحسابات والإحصائيات (الـ Loops والقواميس)
# =========================================================
# حساب أعداد الكروموسومات
chromosome_counts = {}
for variant in variants_list:
    chrom = variant["CHROM"]
    if chrom in chromosome_counts:
        chromosome_counts[chrom] += 1
    else:
        chromosome_counts[chrom] = 1
print("\n Chromosome Counts : \n " , chromosome_counts)

# حساب ال quality
qual_values = []
dots_value = 0
for variant in variants_list:
    if variant["QUAL"] != ".":
        qual_values.append(float(variant["QUAL"]))
    else:
        dots_value += 1
print("\n Quality Values : \n" , qual_values)
print(" Dots value :" , dots_value)
print(" The highest value is : " , max(qual_values))
print(" The lowest value is : " , min(qual_values))
print(" The average value is : " , sum(qual_values) / len(qual_values))

# حساب قيم ال filter
filter_counts = {}
for variant in variants_list:
    filter_val = variant["FILTER"]
    if filter_val in filter_counts:
        filter_counts[filter_val] += 1
    else:
        filter_counts[filter_val] = 1
print("\n Filter counts : \n " , filter_counts)

 # ج) حساب أنواع الطفرات باستخدام REF و ALT مع بعض 🧬
mutation_counts = {"Transition": 0, "Transversion": 0, "Other": 0}
transitions = {("A", "G"), ("G", "A"), ("C", "T"), ("T", "C")}

for variant in variants_list:
    ref = variant["REF"]  # واخد مسافة واحدة لجوه (Tab)
    alt = variant["ALT"]  # واخد مسافة واحدة لجوه (Tab)

    if len(ref) == 1 and len(alt) == 1:  # واخد مسافة واحدة لجوه (Tab)
        if (ref, alt) in transitions:  # واخد مسافتين لجوه (2 Tabs)
            mutation_counts["Transition"] += 1  # واخد 3 مسافات لجوه (3 Tabs)
        else:  # واخد مسافتين لجوه (2 Tabs)
            mutation_counts["Transversion"] += 1  # واخد 3 مسافات لجوه (3 Tabs)
    else:  # واخد مسافة واحدة لجوه (Tab)
        mutation_counts["Other"] += 1  # واخد مسافتين لجوه (2 Tabs)
print("\n Mutation counts : \n " , mutation_counts)

# =========================================================
# 4. مرحلة كتابة التقرير النهائي وحفظه في الملف النصي
# =========================================================
with open("vcf_summary_report.txt", "w") as report:
    report.write("=== VCF Analysis Summary Report ===\n")
    report.write(f"The Highest Quality is: {max(qual_values)}\n")
    report.write(f"The Lowest Quality is: {min(qual_values)}\n")
    report.write(f"The Average Quality is: {sum(qual_values) / len(qual_values)}\n")

    report.write("\n=== Chromosome Counts ===\n")
    for chrom, count in chromosome_counts.items():
        report.write(f"Chromosome {chrom}: {count} variants\n")

print("تم تشغيل الكود بنجاح وكتابة التقرير!")
