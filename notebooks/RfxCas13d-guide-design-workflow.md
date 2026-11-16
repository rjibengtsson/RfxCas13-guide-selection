# RfxCas13 guide design workflow

This notebook describes the workflow for designing RfxCas13d guides, including the steps for filtering and scoring candidate guides. It is intended to be a reference for the methods used in our guide design pipeline.

## 2026-11-13: Mapping PspCas13b pan-sarbecovirus spacers to SARS-CoV-2 MT007544.1

**Aim.** Locate the target site of each previously designed PspCas13b guide on the SARS-CoV-2 genome, and record which gene each site falls in.

**Inputs.**
- `data/pspcas13b_guides/panSarbeCov_guides.tsv` — 8 PspCas13b guides (`Name`, `Spacer`), all spacers 30 nt.
- `data/SARS-CoV-2/MT007544.1/MT007544.1.gbk` and `.fasta` — downloaded from NCBI nuccore on 2026-11-13 via E-utilities `efetch`. Record: SARS-CoV-2 isolate SARS-CoV-2/human/AUS/VIC01/2020, 29,893 bp, 23 annotated features, no ambiguity codes.

**Method.** A Cas13 spacer is the reverse complement of the RNA it targets, so each spacer was reverse-complemented to target-sense DNA before searching. The genome was searched with `blastn -task blastn-short` (word size 11, e-value 1e-3, dust off), and hits annotated against the GenBank feature table.

```
python scripts/blast_guides.py \
    --input data/pspcas13b_guides/panSarbeCov_guides.tsv \
    --seq-column Spacer --id-column Name --revcomp-query \
    --reference data/SARS-CoV-2/MT007544.1/MT007544.1.gbk \
    --output results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_annotated.csv \
    --outfmt6 results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_fmt6.tsv \
    --all-hits results/SARS-CoV-2/panSarbeCov_guides_MT007544.1_all_hits.csv \
    --threads 4
```

**Results.** All 8 guides hit at 100% identity over the full 30 nt, one HSP each, no mismatches or gaps, all on the plus strand. All 8 target sites fall within the N gene (28274–29533; `nucleocapsid phosphoprotein`), each overlapping both the `gene` and `CDS` features.

| Guide | Target coords |
|---|---|
| A8_NCP | 28341–28370 |
| NCP1 | 28433–28462 |
| Cand_2_N_Sarbecov | 28517–28546 |
| Cand_5_N_Sarbecov | 28669–28698 |
| Cand_1_N_Sarbecov | 28747–28776 |
| Cand_4_N_Sarbecov | 29149–29178 |
| Cand_3_N_Sarbecov | 29347–29376 |
| N_cand_93 | 29474–29503 |


## 2026-11-16: RfxCas13d guide design for ACE2 and Cldn5

**Aim.** Design RfxCas13d guides targeting the mouse ACE2 and Cldn5 genes. There are no pre-existing RfxCas13d guides for these genes, so we will design new guides using the RfxCas13d guide design pipeline.

**Method.** First, we will predict candidate guides using the RfxCas13d guide design tool TIGER, which takes place on another server. 

```bash
rsync -avP data/mouse/ACE2.fa rbengtsson@172.26.144.219:/home/unimelb.edu.au/rbengtsson/software/tiger/my_fastas/
rsync -avP data/mouse/Cldn5.fa rbengtsson@172.26.144.219:/home/unimelb.edu.au/rbengtsson/software/tiger/my_fastas/
```

Output from TIGER is stored at `results/mouse/ACE2_on_target.csv` and `results/mouse/Cldn5_on_target.csv`. The next step is to filter and select the candidate guides based on their predicted efficiency and spacing apart along the target gene.

Since for these two genes we have no pre-existing guides, we will select the top 3 predicted guides for each gene based on their predicted efficiency scores. 
