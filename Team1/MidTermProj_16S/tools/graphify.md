```bash
module load python312
python3 -m venv ~/.venvs/graphify
source ~/.venvs/graphify/bin/activate
pip install graphifyy
cd /projects/intro2gds/I2GDS2026/Team1/Kamal/Claude_Code
graphify install
graphify claude install
graphify update .
graphify god-nodes
graphify query "How are reads denoised in QIIME2?"
graphify path "QIIME2 v2024.10" "ASV table (1,817,968 reads, 2,322 ASVs)"
graphify explain "DADA2 denoising / chimera removal (Callahan 2016)"
graphify benchmark
```