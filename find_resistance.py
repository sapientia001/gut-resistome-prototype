import os
import csv
import pyhmmer
import pyrodigal

# 1. Load taxonomy
print("Loading taxonomy mapping...")
taxonomy = {}
with open("taxonomy_mapping.csv", "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        taxonomy[row["Accession"]] = row["Organism"]

# 2. Initialize gene predictor in metagenomic mode (meta=True)
print("Initializing gene predictor...")
orf_finder = pyrodigal.GeneFinder(meta=True)

# 3. Find genomes
data_dir = "ncbi_dataset/data"
genome_folders = [d for d in os.listdir(data_dir) if d.startswith("GCA_") or d.startswith("GCF_")]

print(f"Found {len(genome_folders)} genomes to scan.")
hmm_file_path = "Resfams.hmm"
output_file = "resistance_results.csv"

# List to hold results temporarily for the current batch
batch_results = []

# Create the CSV file and write the header at the start
with open(output_file, "w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=["Organism", "Resistance_Gene"])
    writer.writeheader()

# 4. Scan genomes
for i, folder in enumerate(genome_folders):
    accession = folder
    folder_path = os.path.join(data_dir, folder)
    fna_files = [f for f in os.listdir(folder_path) if f.endswith(".fna")]
    
    if not fna_files:
        continue
        
    fna_path = os.path.join(folder_path, fna_files[0])
    organism = taxonomy.get(accession, "Unknown")
    
    if i % 10 == 0:
        print(f"Scanning {i}/{len(genome_folders)}: {organism}...")
    
    # Read DNA sequence from the .fna file
    with open(fna_path, "r") as f:
        seq_lines = [line.strip() for line in f if not line.startswith(">")]
        dna_seq = "".join(seq_lines)
    
    if not dna_seq:
        continue

    try:
        # Predict genes (ORFs) and get protein sequences
        genes = orf_finder.find_genes(dna_seq)
        protein_seqs = [gene.translate() for gene in genes]
        
        # Convert protein strings to pyhmmer DigitalSequence objects
        digital_proteins = []
        for idx, prot in enumerate(protein_seqs):
            if len(prot) > 20: # Filter very short peptides to save time
                seq_obj = pyhmmer.easel.TextSequence(name=f"protein_{idx}", sequence=prot)
                digital_proteins.append(seq_obj.digitize(pyhmmer.easel.Alphabet.amino()))
        
        if not digital_proteins:
            continue

        # Create a block of digital protein sequences
        protein_block = pyhmmer.easel.DigitalSequenceBlock(pyhmmer.easel.Alphabet.amino(), digital_proteins)
        
        # Open HMM file and search
        with pyhmmer.plan7.HMMFile(hmm_file_path) as hmm_file:
            for hits in pyhmmer.hmmsearch(hmm_file, protein_block):
                for hit in hits:
                    if hit.included:
                        gene_name = hit.name.decode() if isinstance(hit.name, bytes) else hit.name
                        batch_results.append({"Organism": organism, "Resistance_Gene": gene_name})
                        
    except Exception as e:
        print(f"Error processing {organism}: {e}")
        continue

    # Save to CSV every 10 genomes
    if (i + 1) % 10 == 0:
        with open(output_file, "a", newline="") as f:
            writer = csv.DictWriter(f, fieldnames=["Organism", "Resistance_Gene"])
            writer.writerows(batch_results)
        batch_results = [] # Clear the batch after saving
        print(f"  --> Saved progress at genome {i + 1}")

# Save any remaining results at the very end
if batch_results:
    with open(output_file, "a", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=["Organism", "Resistance_Gene"])
        writer.writerows(batch_results)

print("Done! Successfully created resistance_results.csv")