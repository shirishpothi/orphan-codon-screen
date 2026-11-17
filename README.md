# Orphan-codon screen

Looks for bacteria and archaea that use a codon their own tRNAs cannot read, then checks whether their ribosomes skip a base at that codon.

## Why

Uncultured TANB77 gut Clostridia lost their only CGG-reading tRNA and now skip a base at CGG, genome-wide. That case was found by aligning one genome per GTDB family to known proteins with frameshifts allowed. That method is expensive, so it was run on few genomes.

This screen inverts the search. Losing the decoding tRNA is a necessary condition for this kind of skipping, and it's cheap to detect: about 3 seconds per genome to call genes and tRNAs. So every GTDB genus is screened for "orphan" codons first, and only the genomes that have one get the expensive frameshift alignment.

## Pipeline

1. **Sample.** `data/genus_sample.tsv` holds one GTDB (release 232) species representative per genus: 28,344 genomes with CheckM2 completeness ≥ 90% and contamination ≤ 5%. Within a genus, the representative with the highest completeness − 5 × contamination wins.
2. **Scan** (`scripts/pipeline.py`). Download from NCBI, call genes with pyrodigal and tRNAs with pyaragorn, then count codons and anticodons.
3. **Orphan call** (`scripts/orphans.py`). Wobble rules are generous, so fewer codons are called orphan:
   - C34 reads G
   - G34 reads C and U
   - inosine (A34) reads U, C and A
   - U34 reads all four codons in a four-fold box, or A and G in a two-fold box
   - CAU anticodons are also allowed to read AUA

   Ile and Ala codons are ignored, because the Ile-GAU and Ala-UGC tRNAs sit inside rRNA operons and are often missing from MAGs.
4. **Frameshift check.** Genomes with 1–4 orphan codons (≥ 20 uses, ≤ 5 per 1,000 codons, ≥ 25 tRNAs) go through `diamond blastx -F 15 --fast` against Swiss-Prot. `scripts/fscount.py` locates each +1 and −1 frameshift from the BTOP string. `scripts/codonaa.py` tabulates what amino acid the homologs have at each in-frame codon, in the style of Codetta.
5. **Merge** (`scripts/merge.py`) writes:

| File | Contents |
|---|---|
| `hits.tsv` | Orphan codon that is also a +1 frameshift hotspot |
| `hotspots.tsv` | Any codon holding ≥ 30% of a genome's frameshifts |
| `lineage.tsv` | Orphan codons shared by most genera of a GTDB order |
| `summary.md` | Overview of all of the above |

## Running

Actions → *orphan-codon screen* → Run workflow. It runs 60 shards on GitHub-hosted runners, at most 20 at once, then a merge job. Results are in the `merged` artifact and the run summary.

Locally:

```sh
pip install -r requirements.txt   # plus DIAMOND 2.2.9 on PATH
diamond makedb --in uniprot_sprot.fasta.gz -d sprot
python scripts/pipeline.py data/genus_sample.tsv 0 1 out sprot
python scripts/merge.py data/genus_sample.tsv out merged
```

## Validation so far

- **Positive control.** The blind screen recovers TANB77 CGG skipping without being pointed at it: CGG is orphaned in every sampled TANB77 genus, and in the validated genomes 23–84 of their +1 frameshifts fall at CGG, with glycine in the homolog at most of those sites.
- **Negative controls.** Normal Clostridia show 13–45 scattered +1 frameshifts with no codon preference.
- **Known caveat.** Most orphan codons are read as their usual amino acid, which means a tRNA was missed by the assembly or by ARAGORN. ARAGORN does not detect intron-containing archaeal tRNAs, so archaeal Cys orphans are suspect.

## Limitations

- tRNA calls come from ARAGORN only, not tRNAscan-SE.
- Swiss-Prot is sparse for uncultured lineages, so frameshift counts per genome are modest.
- A hit is a candidate, not a result. Confirming one needs more genomes from the same lineage, phage genomes, and raw reads to rule out assembly error, as was done for TANB77.
