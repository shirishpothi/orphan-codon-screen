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

Actions → *orphan-codon screen* → Run workflow. It runs 60 shards on GitHub-hosted runners, at most 20 at once, then a merge job. Results are in the `merged` artifact and the run summary. Inputs `sample` (default `data/genus_sample.tsv`) and `validate_all` (`1` aligns every genome) run other samples; the TANB77 species run used `shards=20 sample=data/tanb77_sample.tsv validate_all=1`, then `python scripts/tanb77.py data/tanb77_sample.tsv merged`.

Locally:

```sh
pip install -r requirements.txt   # plus DIAMOND 2.2.9 on PATH
diamond makedb --in uniprot_sprot.fasta.gz -d sprot
python scripts/pipeline.py data/genus_sample.tsv 0 1 out sprot
python scripts/merge.py data/genus_sample.tsv out merged
```

## Results (GTDB r232, full run, 17 Nov 2026)

28,338 of 28,344 genera were scanned (6 downloads failed). 7,364 genomes had an orphan codon and went through frameshift alignment. Tables are in `findings/`.

**1. TANB77 is the only lineage that skips at an orphan codon.** All 32 orphan-codon frameshift hotspots are TANB77 CGG. The blind screen reproduces the published numbers closely:

| | This screen | Published |
|---|---|---|
| TANB77 genera with no CGG-reading tRNA | 50 of 52 | 54 of 64 |
| +1 frameshifts at CGG | 1,020 of 1,953 (52%), in 40 genomes | 1,026 of 1,744 (59%) |
| Homolog has glycine at the CGG skip site | 897 of 1,020 (88%) | 1,237 of 1,458 |
| Glycine at other +1 frameshifts | 49 of 933 (5%) | 6% |

**2. Within TANB77, keeping the tRNA means not skipping.** This was tested on every TANB77 species representative: 788 genomes (782 aligned), each aligned regardless of orphan status (`data/tanb77_sample.tsv`, run with `validate_all=1`; per-genome table in `findings/tanb77_cgg_by_genome.tsv`, made by `scripts/tanb77.py`).

| Genomes (≥ 25 tRNAs unless noted) | +1 frameshifts at CGG | Homolog Gly at CGG shift | Gly at other +1 shifts |
|---|---|---|---|
| No CGG-reading tRNA: 611 genomes, 54 genera | 20,572 of 37,801 (54%) | 18,160 of 20,572 (88%) | 863 of 17,229 (5%) |
| tRNA-Arg-CCG present: 13 genomes, any tRNA count | 2 of 137 (1.5%) | — | — |

The difference is Fisher p = 3 × 10⁻⁴³. 551 of the 611 tRNA-less genomes individually pass the skipper threshold (≥ 3 CGG shifts and ≥ 20% of +1 shifts). The rest form the low tail of a continuous distribution, scattered across genera, not a separate lineage.

On the GTDB r232 bac120 tree, family RGIG9151 is sister to all other TANB77 families. All 9 RGIG9151 genomes carry tRNA-Arg-CCG and together show 1 CGG shift in 114. The other four families (CAG-508, CAG-465, UBA1234, SIG350) have essentially all lost it. The simplest history is one loss of tRNA-CCG on the stem of the other four families, with skipping arising alongside it. The 4 CCG carriers outside RGIG9151 (2 in *CAG-465*, 2 in *CAG-793*) sit deep inside skipping clades and show 1 CGG shift in 23. They could be secondary regains or contamination; that is untested.

Two CAG-508 genomes (*CAG-269* GCA_947573965.1, *CAG-452* GCA_944392485.1) carry a tRNA-Arg-UCG besides the usual ACG, yet skip strongly (53 of 87 and 133 of 193 +1 shifts at CGG). The generous U34 rule assumes UCG reads CGG, but here it evidently does not. Either the tRNA is inactive or modified, or it isn't the organism's own.

**3. Losing a tRNA almost never leads to skipping.** In the other 7,321 genomes with an orphan codon, 3,000 of 284,125 +1 frameshifts (1.1%) fall at an orphan codon. That is no more than chance for 1–4 codons out of 61. Their homologs show the codon's usual amino acid (e.g., CCA→Pro 82% in Oscillospirales, TTA→Leu 74% in Burkholderiales). So these orphans are mostly tRNAs that the assembly or ARAGORN missed, not new decoding.

**4. CGG has repeatedly lost its tRNA in Mollicutes-related Bacilli_A, without skipping:**

| Order | Genera with no CGG-reading tRNA | +1 frameshifts at CGG |
|---|---|---|
| RF39 | 82 of 142 | 3 of 383 |
| Mycoplasmatales | 32 of 44 | 0 of 184 |
| MAG-NZ | 4 of 4 | — |
| JAJDGJ01 | 3 of 3 | — |
| UBA3375 | 3 of 3 | — |

CGG is very rare in these genomes (~0.1–0.3 per 1,000 codons). This matches the known unassigned, avoided CGG of *Mycoplasma capricolum* (Oba et al. 1991, PNAS 88:921) and the absence of tRNA-Arg-CCG across Mollicutes (Yokobori et al. 2013, NAR 41:6531). This is an extension of known biology to RF39 and three uncultured orders, not a new phenomenon.

**5. The only hotspot not at an orphan codon is an artifact.** In *Virgibacillus alimentarius* (RS_GCF_000709085.1, a 146-contig assembly), 22 of 55 +1 frameshifts sit at AAA, and 21 of those 22 are inside runs of five or more A's. That is homopolymer sequencing error.

**Controls.** Normal Clostridia show 13–45 scattered +1 frameshifts with no codon preference. ARAGORN does not detect intron-containing archaeal tRNAs, so the archaeal Cys orphans (Sulfolobales and others) in `lineage.tsv` are probably artifacts.

## Limitations

- tRNA calls come from ARAGORN only, not tRNAscan-SE.
- Swiss-Prot is sparse for uncultured lineages, so frameshift counts per genome are modest.
- A hit is a candidate, not a result. Confirming one needs more genomes from the same lineage, phage genomes, and raw reads to rule out assembly error, as was done for TANB77.
