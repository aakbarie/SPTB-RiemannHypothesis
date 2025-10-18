# Data Directory

## Odlyzko Zero Data

### Source

The zero data used in this analysis comes from Andrew Odlyzko's public tables:

**Website**: https://www-odlyzko.dtc.umn.edu/zeta_tables/  
**File used**: `zeros1` (first 300,000 nontrivial zeros of ζ(s))

### Download Instructions
```bash
# 1. Visit the Odlyzko website
# 2. Download zeros1 file (or use wget)
wget https://www-odlyzko.dtc.umn.edu/zeta_tables/zeros1

# 3. Place in data/ directory
mv zeros1 data/zeros_1.txt
```

### Format

Plain text file, one zero per line (imaginary parts γ only):
```
14.134725141734693790457251983562470270784257115699243175685567460149963429809...
21.022039638771554992628479593896902777334340524902781754629520403587598586068...
25.010857580145688763213790992562821818659549672557996672496542006745092098942...
...
```

All zeros satisfy β = 1/2 (verified computationally).

### File Size

- **zeros_1.txt**: ~50 MB (300,000 zeros, ~170 bytes each)

### Citation

If you use this data, please cite:
```
A.M. Odlyzko, Tables of zeros of the Riemann zeta function,
Available online: https://www-odlyzko.dtc.umn.edu/zeta_tables/
```

---

## Processed Data

After running `code/sptb_analysis.R`, the following files are generated:

### bias_summary.csv

Summary statistics for each (η, T) combination in the bias regime.

**Columns**:
- `Tmax`: Maximum time T
- `eta`: Off-line displacement η = β - σ
- `Delta`: Block width Δ
- `kappa_bias`: Adaptive κ parameter
- `n_t`: Number of time samples
- `blocks`: Number of blocks
- `gamma0`: Synthetic zero frequency
- `lambda`: Derivative penalty λ
- `emp_tail_rate`: Empirical exponential slope
- `theory_2eta`: Theoretical slope 2η

**Rows**: 9 (3 T × 3 η)

### bias_blocks.csv

Per-block cumulative energy differences for bias regime.

**Columns**:
- `blk`: Block index
- `block_start`: Block start time
- `inc_total_diff`: Incremental energy difference
- `cum_total_diff`: Cumulative energy difference
- `eta`: Off-line displacement
- `Tmax`: Maximum time

**Rows**: ~133,980 (many blocks across 9 scenarios)

### variance_table.csv

Variance regime test results (all zeros on β = 1/2).

**Columns**:
- `Tmax`: Maximum time T
- `Delta`: Block width
- `lambda`: Derivative penalty
- `blocks`: Number of blocks
- `sum_eps2`: ∑ ∫|ε|² (amplitude term)
- `sum_d2`: ∑ ∫|∂ε|² (derivative term)
- `F_lambda`: Total SPTB functional
- `T_logT_loglogT`: T log T log log T
- `T_log2T`: T (log T)²

**Rows**: 6 (different T values)

---

## License

Data files are derived from public-domain sources (Odlyzko) and processed 
outputs are released under CC-BY-SA 4.0.