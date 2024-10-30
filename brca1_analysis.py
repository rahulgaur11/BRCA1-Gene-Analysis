# brca1_analysis.py

# Step 1: Import Libraries
from Bio import SeqIO
from collections import Counter
import matplotlib.pyplot as plt
import streamlit as st

# Step 2: Load and Parse the BRCA1 FASTA Sequence
# Load sequence from FASTA file
record = SeqIO.read("BRCA1 FASTA", "fasta")
sequence = record.seq

# Sequence Information
sequence_length = len(sequence)
sequence_preview = sequence[:100]

# Base Composition and GC Content
base_counts = Counter(sequence)
gc_content = (base_counts["G"] + base_counts["C"]) / sequence_length * 100
at_content = (base_counts["A"] + base_counts["T"]) / sequence_length * 100
at_gc_ratio = at_content / gc_content if gc_content > 0 else 0

# Motif Search (Example: ATG start codon)
motif = "ATG"
positions = [i for i in range(len(sequence) - len(motif) + 1) if sequence[i:i + len(motif)] == motif]

# Protein Translation
protein_sequence = sequence.translate(to_stop=True)
protein_length = len(protein_sequence)
protein_preview = protein_sequence[:100]

# Codon Usage Analysis
coding_sequence = sequence[:(sequence_length // 3) * 3]
codon_counts = Counter(str(coding_sequence[i:i+3]) for i in range(0, len(coding_sequence), 3))

# Open Reading Frames (ORFs)
def find_orfs(seq, min_length=100):
    start_codon = "ATG"
    stop_codons = {"TAA", "TAG", "TGA"}
    orfs = []
    for i in range(len(seq) - 2):
        if str(seq[i:i+3]) == start_codon:
            for j in range(i + 3, len(seq) - 2, 3):
                if str(seq[j:j+3]) in stop_codons:
                    if (j - i) >= min_length:
                        orfs.append(seq[i:j+3])
                    break
    return orfs

orfs = find_orfs(sequence)

# CpG Island Analysis
def find_cpg_islands(seq, window_size=200, threshold=0.6):
    cpg_islands = []
    for i in range(0, len(seq) - window_size + 1):
        window = seq[i:i + window_size]
        cpg_count = window.count("CG")
        cpg_ratio = cpg_count / (window_size / 2)
        if cpg_ratio >= threshold:
            cpg_islands.append((i, i + window_size))
    return cpg_islands

cpg_islands = find_cpg_islands(sequence)

# GC and AT Skew Analysis
def calculate_skew(seq, window_size=100):
    gc_skew = []
    at_skew = []
    for i in range(0, len(seq) - window_size + 1, window_size):
        window = seq[i:i + window_size]
        g = window.count("G")
        c = window.count("C")
        a = window.count("A")
        t = window.count("T")
        gc_skew.append((g - c) / (g + c) if (g + c) > 0 else 0)
        at_skew.append((a - t) / (a + t) if (a + t) > 0 else 0)
    return gc_skew, at_skew

gc_skew, at_skew = calculate_skew(sequence)

# Step 3: Streamlit Dashboard
st.title("BRCA1 Gene Analysis Dashboard")
st.write("This dashboard provides an interactive overview of the BRCA1 gene, a critical gene involved in DNA repair. Each section includes explanations for better understanding.")

# Sequence Information
st.header("Sequence Information")
st.write("The BRCA1 gene sequence consists of DNA letters (A, T, C, G) that encode instructions for protein production.")
st.write(f"**Sequence Length**: {sequence_length} nucleotides")
st.write(f"**Sequence Preview**: {sequence_preview}")

# Base Composition
st.header("Base Composition")
st.write("Base composition shows the number of each DNA letter (A, T, C, G) in the gene.")
st.write(f"**Base Composition**: {base_counts}")
st.write(f"**GC Content**: {gc_content:.2f}% (Higher GC content indicates greater stability)")
st.write(f"**AT Content**: {at_content:.2f}%")
st.write(f"**AT/GC Ratio**: {at_gc_ratio:.2f}")

# Motif Search
st.header("Motif Search")
st.write("Motifs are short DNA sequences with important functions. Here, we search for the 'ATG' motif, which often marks the start of a gene.")
st.write(f"**Motif '{motif}' found at positions**: {positions[:10]}... (showing first 10 positions)")

# Protein Translation
st.header("Protein Translation")
st.write("Translation converts the DNA sequence into a protein, which helps repair DNA. Below are details of the protein encoded by BRCA1.")
st.write(f"**Protein Sequence Length**: {protein_length}")
st.write(f"**Protein Sequence Preview**: {protein_preview}")

# Codon Usage
st.header("Codon Usage")
st.write("DNA letters are grouped into triplets called codons. Each codon codes for a specific amino acid. This analysis shows how often each codon appears in the BRCA1 gene.")
st.write(f"**Codon Usage**: {codon_counts}")

# Open Reading Frames (ORFs)
st.header("Open Reading Frames (ORFs)")
st.write("Open Reading Frames (ORFs) are sections of DNA that can be translated into protein.")
st.write(f"**Number of ORFs Found**: {len(orfs)}")
st.write(f"**Example ORF**: {orfs[0] if orfs else 'No ORFs found'}")

# CpG Island Analysis
st.header("CpG Island Analysis")
st.write("CpG islands are regions with a high frequency of CG pairs, often found at the start of genes.")
st.write(f"**CpG Islands Found**: {cpg_islands[:5]}... (showing first 5)")

# GC and AT Skew
st.header("GC and AT Skew Analysis")
st.write("GC and AT skew measure the balance between G vs C and A vs T across the sequence. This can highlight important regions.")
fig, ax = plt.subplots()
ax.plot(gc_skew, label="GC Skew")
ax.plot(at_skew, label="AT Skew")
ax.set_xlabel("Position (windowed)")
ax.set_ylabel("Skew")
ax.legend()
st.pyplot(fig)

# Step 6: GC Content Across the Gene Sequence (Sliding Window)
window_size = 100  # Adjust window size as desired
gc_content_values = [
    (sequence[i:i + window_size].count("G") + sequence[i:i + window_size].count("C")) / window_size * 100
    for i in range(0, sequence_length, window_size)
]

# GC Content Visualization in Streamlit
st.header("GC Content Visualization")
st.write("The plot below shows GC content across the BRCA1 gene, which affects gene stability.")
fig, ax = plt.subplots()
ax.plot(gc_content_values)
ax.set_title("GC Content Across BRCA1 Gene Sequence")
ax.set_xlabel("Position (windowed)")
ax.set_ylabel("GC Content (%)")
st.pyplot(fig)


# Glossary Section
st.sidebar.title("Glossary")
st.sidebar.write("**Base Composition**: Counts of A, T, C, and G in the sequence.")
st.sidebar.write("**GC Content**: Percentage of G and C in the gene; high GC content implies stability.")
st.sidebar.write("**Codon**: A group of three DNA letters that codes for a protein part.")
st.sidebar.write("**ORF (Open Reading Frame)**: Section of DNA that can be translated into protein.")
st.sidebar.write("**CpG Island**: DNA regions rich in CG pairs, often found near gene starts.")
st.sidebar.write("**Skew**: Difference between G vs C or A vs T in specific regions.")

st.write("Explore the glossary on the sidebar for definitions of technical terms.")
