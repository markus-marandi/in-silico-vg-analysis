# Diagnostic Checks for TSS/eQTL Analysis

## Questions to Answer:

### 1. What spatial region does `vg_eqtl` measure?
```python
# Check: Is vg_eqtl gene-wide or promoter-specific?
# Look at the column definition or how it was calculated
# Expected answer: Gene-wide eQTL variance OR TSS±200bp-specific
```

### 2. Sample sizes in Panel B (TSS bins)
```python
# For each TSS bin, how many genes have vg_eqtl data?
for b in ['1–2', '3–5', '6–15', '16+']:
    hi = tss_df[(tss_df['is_hi']) & (tss_df['tss_bin'] == b)]['vg_eqtl'].dropna()
    bg = tss_df[(~tss_df['is_hi']) & (tss_df['tss_bin'] == b)]['vg_eqtl'].dropna()
    print(f"TSS {b}: HI n={len(hi)}, BG n={len(bg)}")
```

### 3. What does the normalized metric actually measure?
```python
# Check if HI has higher OR lower variant density
hi_density = tss_df[tss_df['is_hi']]['variants_per_kb'].median()
bg_density = tss_df[~tss_df['is_hi']]['variants_per_kb'].median()
print(f"HI variant density: {hi_density:.1f} vars/kb")
print(f"BG variant density: {bg_density:.1f} vars/kb")

# If HI has LOWER density but HIGHER normalized effect, that's consistent
# If HI has SIMILAR density with HIGHER normalized effect, that's different biology
```

### 4. Are you comparing apples to apples?
- Panel A: Enformer predictions on **all variants** in TSS±200bp
- Panel B: Empirical **eQTL** Vg (detected effects only) - where?
- Panel C: Normalized by **gene-wide** variant density

These may be measuring different things:
- Enformer: predictive constraint on all variants
- eQTL: empirical variance from detected cis-eQTLs (could be distal!)

## Possible Interpretations:

### Scenario 1: vg_eqtl is gene-wide
- Panel A measures promoter-core constraint ✓
- Panel B measures total eQTL variance (not just proximal) ✗ mismatch
- Panel C normalizes gene-wide eQTL by gene-wide density
- **Interpretation**: Different spatial regions → not "or

thogonal validation"

### Scenario 2: vg_eqtl is promoter-core-specific
- Panel B shows HI has MORE raw promoter eQTL Vg (contradicts Panel A!)
- Panel C shows HI has HIGHER per-variant effect
- This would mean: HI has stronger eQTL effects in promoter despite lower Enformer predictions
- **Interpretation**: eQTL discovery bias OR different biological signal

## Action Items:

1. **Verify what `vg_eqtl` measures spatially** (gene-wide vs promoter)
2. **Check sample sizes** in Panel B bins (low N → high p-values)
3. **Calculate promoter-specific eQTL Vg** if it doesn't exist
4. **Reconsider the "orthogonal validation" claim** if measuring different regions
