# Collaborator setup — HiggsDNA `v7_noIso_bbgg` + bbgg Run3 production

This is the short path from zero to submitting Condor for your assigned era / samples.  
Full inventory and status tables: see analysis repo `scripts/v7_campaign/PRODUCTION_TEAM.md`.

## What you need from us

| Piece | Where |
|---|---|
| DNA code (this repo) | https://github.com/ZhenyuDongLL/HiggsDNA-bbggrun3 — branch **`v7_noIso_bbgg`** |
| Submit scripts + ledgers + team MD | https://github.com/ZhenyuDongLL/bbgg-pnn — path `scripts/v7_campaign/` |
| Group dump (production) | `/eos/cms/store/group/phys_higgs/cmshgg/zhenyudong/bbgg_sample/v7_DNA_noIso_v1/` |
| Optional personal btag dump | `/eos/home-z/zhdong/Xbbgg/v7_DNA_noIso_v1/btag_eff/` (Zhenyu); use your own EOS for new btag trials if preferred |

Do **not** commit or share local `runner_*.json` / `samples_*.json|txt` — scripts recreate them.

---

## 1. Clone DNA (this branch)

```bash
git clone -b v7_noIso_bbgg https://github.com/ZhenyuDongLL/HiggsDNA-bbggrun3.git HiggsDNA-bbggrun3
cd HiggsDNA-bbggrun3
```

Base is official `HHbbgg_v7_parquet` plus our patches (no abs-iso in `store_flag`, 13p6TeV trees, NMSSM/GGJets process map, FNUF/Material in `variations_mc.json`).

---

## 2. Conda + install HiggsDNA

Official pattern (rename the env if you like; Zhenyu uses `higgs-dna-v7`):

```bash
# Miniforge / mamba recommended on lxplus
mamba create -n higgs-dna-v7 python=3.12 xrootd -y
mamba activate higgs-dna-v7

cd /path/to/HiggsDNA-bbggrun3
pip install -e ".[dev,test]"
```

Check:

```bash
python -c "import higgs_dna; print(higgs_dna.__file__)"
which run_analysis.py   # should be in the env bin after editable install
```

Working directory for submit scripts should be the package root that contains `scripts/` and `submission/`, typically:

```bash
cd /path/to/HiggsDNA-bbggrun3/higgs_dna
```

(Our campaign scripts `chdir` there themselves when paths are set correctly.)

---

## 3. Pull correction / SF JSON files (`pull_files`)

JSON under `higgs_dna/systematics/JSONs/` are **gitignored**. After install you must download them.

From the DNA package (where `scripts/pull_files.py` lives), with a CMS proxy if EOS/xrootd requires it:

```bash
mamba activate higgs-dna-v7
export X509_USER_PROXY=$HOME/private/x509up_cms   # or your path
voms-proxy-info -e || voms-proxy-init --rfc --voms cms -valid 192:00 -out $HOME/private/x509up_cms

cd /path/to/HiggsDNA-bbggrun3/higgs_dna

# Full set (slow but safest once):
python scripts/pull_files.py --all --use-xrdcp

# Or the minimum set this campaign actually needs for Run3 HHbbgg:
for t in \
  GoldenJSON PU SS JetMET JEC JER Material FNUF \
  TriggerSF PreselSF eVetoSF LooseMva PhotonID \
  bTag HHbbgg_bTag_WPs HHbbgg_bTag_eff Flows \
  HHbbgg_bpairing HHbbgg_vbfpairing HHbbgg_mbb_reg_model
do
  echo "===== $t ====="
  python scripts/pull_files.py --target "$t" --use-xrdcp
done
```

Confirm at least:

```bash
ls systematics/JSONs/FNUF/2022/
ls systematics/JSONs/Material/2018/
ls systematics/JSONs/bTagEff/2024_Summer24/   # names vary by era
ls systematics/JSONs/scaleAndSmearing/
```

**FNUF** and **Material** are required so production writes those object-syst directories (nominal weight 1 on Run3 fallbacks). Without them, runner patches that list FNUF/Material will fail at runtime.

---

## 4. Clone / sync campaign scripts

```bash
git clone https://github.com/ZhenyuDongLL/bbgg-pnn.git
# or: git pull if you already have it
ls bbgg-pnn/scripts/v7_campaign/
```

Read:

- `PRODUCTION_TEAM.md` — sample counts, eras, phases, who owns what  
- `SETUP_FOR_COLLABORATORS.md` — this file may also live under `scripts/v7_campaign/` as a copy  

Point scripts at **your** DNA checkout if paths differ: edit `HIGGS_DNA = ...` near the top of `produce_batchA_2024.py` / `produce_s2_btagEff.py` (default is Zhenyu’s EOS path).

---

## 5. Proxy + Condor (every submit session)

```bash
mamba activate higgs-dna-v7
export X509_USER_PROXY=$HOME/private/x509up_cms
export PYTHONUNBUFFERED=1
voms-proxy-info -e
```

Put the proxy on **AFS** (`~/private/...`), not node-local `/tmp`, so Condor workers see the same file. HiggsDNA’s lxplus wrapper also copies proxy under `$HOME/.x509up_*`.

---

## 6. Run the process for your assigned era / samples

### 6.1 Current recommended entry: Batch A (2024)

Finishes remaining 2024 data production and submits **btag-eff** for NMSSM MX ∈ {1000,2000,3000,4000} (77 points). Skips anything already in the TSV ledgers.

```bash
python /path/to/bbgg-pnn/scripts/v7_campaign/produce_batchA_2024.py
```

| Output | Path (default in script) |
|---|---|
| Production parquet | `.../cmshgg/.../v7_DNA_noIso_v1/dna_dump/2024` |
| Signal btag pkls | `/eos/home-z/zhdong/Xbbgg/v7_DNA_noIso_v1/btag_eff/2024` (change if you use your EOS) |

Before you run: **claim** your slice in the team chat and check:

```bash
column -t scripts/v7_campaign/batchA_2024_submitted.tsv
column -t scripts/v7_campaign/btag_eff_submitted.tsv
```

If `condor_submit` succeeded but the Python process died before writing the TSV, append the row by hand so nobody double-submits.

### 6.2 btag-eff only (filtered) — do not relaunch full 949

`produce_s2_btagEff.py` was written for all eras × all signals. For teamwork, either:

- extend Batch-A style filters (year / MX), or  
- temporarily comment / filter `todo` in that script  

Then merge pkls into the era `HHbbgg.json.gz` **before** signal (or 22/23 MC) production. See `PRODUCTION_TEAM.md` §5.

### 6.3 Production after keys exist

Use HiggsDNA `produce_one_mc.py` / `produce_one_data.py` patterns wrapped by our scripts: FNUF/Material via `append_fnuf_material.patch_runner_file` for MC; data has no FNUF. Memory: **32 GB** for GGJets production, 16 GB otherwise; btag usually 10 GB.

---

## 7. Monitor

```bash
condor_q -totals
# Held ExitCode 107 + Errno 5 on conda imports = transient EOS read; wait for periodic_release
# Successful jobs leave the queue — “0 completed” is normal
```

Job logs: `higgs_dna/.higgs_dna_vanilla_lxplus/runner_*_<timestamp>/jobs/*.err`

---

## 8. After your batch finishes

1. Update `PRODUCTION_TEAM.md` §9 (or the shared TSV) and push to `bbgg-pnn`.  
2. For btag: run `scripts/btagging_eff.py` then merge into `systematics/JSONs/bTagEff/<era>/HHbbgg.json.gz` (backup `.orig`). Coordinate so only one person merges a given era.  
3. Only then submit signal production for that era.

---

## 9. Troubleshooting

| Symptom | Action |
|---|---|
| `OSError: [Errno 5] Input/output error` during `import` on lxplus | Retry; Batch A already retries submit-node failures. Do not redesign conda path unless the whole team agrees. |
| `exit()` in bTagMultiFixedWP | Dataset name missing from that era’s `HHbbgg.json.gz` — finish btag-eff + merge first. |
| Ambiguous `HHbbgg_v7_parquet` | Branch and tag share the name; use `refs/heads/...` or our `v7_noIso_bbgg` branch only. |
| Empty DAS file list | Check proxy / DAS name; 2024 signal overrides for (1200,70)/(1200,400) are in `signal_mass_points.py`. |
