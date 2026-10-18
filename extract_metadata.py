import json
import csv

report_file = r"ncbi_dataset\data\assembly_data_report.jsonl"

with open("taxonomy_mapping.csv", "w", newline="") as csvfile:
    writer = csv.writer(csvfile)
    writer.writerow(["Accession", "Organism"])

    with open(report_file, "r") as f:
        for line in f:
            data = json.loads(line)
            organism = data.get("organism", {}).get("organismName", "Unknown")
            accession = data.get("accession", "Unknown")
            writer.writerow([accession, organism])

print("Successfully created taxonomy_mapping.csv")