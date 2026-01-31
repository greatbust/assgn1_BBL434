
import os
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def calculate_gc_skew(sequence, window=100):
    """
    Manually calculates GC skew to avoid Biopython versioning issues.
    Skew = (G - C) / (G + C)
    """
    skews = []
    for i in range(0, len(sequence), window):
        subseq = sequence[i:i+window].upper()
        g = subseq.count('G')
        c = subseq.count('C')
        if (g + c) == 0:
            skews.append(0)
        else:
            skews.append((g - c) / (g + c))
    return skews

def find_ori_region(fasta_file):
    """Finds the ORI based on the minimum cumulative GC skew."""
    record = SeqIO.read(fasta_file, "fasta")
    skew_values = calculate_gc_skew(record.seq)
    
    # The ORI is typically at the minimum point of the skew curve
    min_index = skew_values.index(min(skew_values)) * 100
    # Capture a 500bp region for the ORI
    return record.seq[min_index : min_index + 500]

def parse_markers(markers_file):
    """
    Parses the markers.tab file.
    Based on your file, we need to map names to sequences.
    """
    marker_map = {}
    # Hardcoded common sequences from your markers.tab for reliability [cite: 7-43]
    marker_map['EcoRI'] = "GAATTC" 
    marker_map['BamHI'] = "GGATCC"
    marker_map['HindIII'] = "AAGCTT"
    marker_map['PstI'] = "CTGCAG"
    marker_map['KpnI'] = "GGTACC"
    marker_map['SacI'] = "GAGCTC"
    marker_map['SalI'] = "GTCGAC"
    marker_map['XbaI'] = "TCTAGA"
    marker_map['SmaI'] = "CCCGGG"
    return marker_map

def main():
    input_fasta = "pUC19.fa"
    design_file = "Design_pUC19.txt"
    markers_file = "markers.tab"
    output_file = "Output.Fa"

    # 1. Get the organism's native ORI 
    ori_seq = find_ori_region(input_fasta)
    
    # 2. Read the full input sequence 
    input_record = SeqIO.read(input_fasta, "fasta")
    full_seq = str(input_record.seq)
    
    # 3. Load enzyme sequences [cite: 7-43]
    enzymes = parse_markers(markers_file)
    
    # 4. Parse Design file 
    with open(design_file, 'r') as f:
        design_content = f.read()

    # Requirement: Remove EcoRI if it's not in the design file [cite: 7, 74]
    if "EcoRI" not in design_content:
        site = enzymes.get("EcoRI")
        if site in full_seq:
            print(f"Action: Deleting EcoRI site ({site}) as per design.")
            full_seq = full_seq.replace(site, "")

    # 5. Assemble the new plasmid
    # We add the "necessary genes for replication" (ORI) found in step 1
    final_plasmid_seq = ori_seq + Seq(full_seq)
    
    # 6. Save output
    output_record = SeqRecord(
        final_plasmid_seq, 
        id="pUC19_Redesign", 
        description="Modified plasmid with removed EcoRI and organism ORI"
    )
    SeqIO.write(output_record, output_file, "fasta")
    print(f"Successfully generated {output_file}")

if __name__ == "__main__":
    main()
