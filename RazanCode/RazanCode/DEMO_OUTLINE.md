# Code review demo: 15 minutes (Razan + Manisha)

| Time | Who | Content |
|---|---|---|
| 0-2 | Razan | Goal: reproduce a published canker 16S analysis from public data (PRJNA1348523). Headline: canker results reproduce, healthy results do not. |
| 2-7 | Razan | Walk the repo: README table, numbered scripts, one step end to end (DADA2 CCS -> taxonomy -> filter). Show `results/pdf/`. |
| 7-11 | Razan | One decision in depth: paper-style vs filtered vs relaxed tables. Evidence: Fig 4A has no chloroplast; PCoA 24.6/7.6% vs paper 23.7/6.9%. |
| 11-13 | Manisha | Student pipeline: one settings file, `run_pipeline.sh`, rerun from a step (`4-11`). |
| 13-15 | Both | Questions. |

## Likely questions to prepare
1. How were the DADA2 settings chosen when the paper gives none?
2. Why are healthy samples ~99.5% chloroplast, and how do you know the ~1,820 bp reads are host SSU?
3. Why FastSpar instead of MENAP, and why are the networks sparser?
4. What does rarefying vs not rarefying change?
5. How would you make every step run from one command?

## Things for reviewers to look at (code review prompts)
- Hardcoded paths and Slurm account in scripts
- Steps that were re-run by hand after fixes (5, 12)
- Unweighted vs weighted distance bug in `12b_network_analysis.R`
- Error handling and logging in the Slurm scripts
